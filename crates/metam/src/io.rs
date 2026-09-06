use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
};

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;

use crate::sequence::normalize_base;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record {
    pub id: String,
    pub seq: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct RefRecord {
    pub id: String,
    pub seq: Vec<u8>,
}

pub fn open_reader(path: &str) -> Result<Box<dyn BufRead + Send>> {
    let file = File::open(path).with_context(|| format!("failed to open {path}"))?;
    let reader: Box<dyn Read + Send> = if path.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader)))
}

pub fn read_references(path: &str) -> Result<Vec<RefRecord>> {
    let mut reader = FastaReader::new(open_reader(path)?);
    let mut refs = Vec::new();
    while let Some(record) = reader.next_record()? {
        refs.push(RefRecord {
            id: record.id,
            seq: record.seq,
        });
    }
    if refs.is_empty() {
        bail!("reference FASTA contains no records");
    }
    let mut ids = refs.iter().map(|r| r.id.as_str()).collect::<Vec<_>>();
    ids.sort_unstable();
    if ids.windows(2).any(|w| w[0] == w[1]) {
        bail!("reference FASTA contains duplicate sequence IDs");
    }
    Ok(refs)
}

pub struct SequenceReader {
    inner: ReaderKind,
}

enum ReaderKind {
    Fasta(FastaReader),
    Fastq(FastqReader),
}

impl SequenceReader {
    pub fn from_path(path: &str) -> Result<Self> {
        let mut reader = open_reader(path)?;
        let first = reader
            .fill_buf()
            .with_context(|| format!("failed reading {path}"))?
            .first()
            .copied();
        match first {
            Some(b'>') => Ok(Self {
                inner: ReaderKind::Fasta(FastaReader::new(reader)),
            }),
            Some(b'@') => Ok(Self {
                inner: ReaderKind::Fastq(FastqReader::new(reader)),
            }),
            Some(other) => bail!(
                "could not detect FASTA/FASTQ for {path}: first byte is {:?}",
                other as char
            ),
            None => bail!("{path} is empty"),
        }
    }

    pub fn next_record(&mut self) -> Result<Option<Record>> {
        match &mut self.inner {
            ReaderKind::Fasta(reader) => reader.next_record(),
            ReaderKind::Fastq(reader) => reader.next_record(),
        }
    }

    pub fn read_batch(&mut self, batch_size: usize) -> Result<Vec<Record>> {
        let mut batch = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            match self.next_record()? {
                Some(record) => batch.push(record),
                None => break,
            }
        }
        Ok(batch)
    }
}

pub struct FastaReader {
    reader: Box<dyn BufRead + Send>,
    pending_header: Option<String>,
    finished: bool,
}

impl FastaReader {
    pub fn new(reader: Box<dyn BufRead + Send>) -> Self {
        Self {
            reader,
            pending_header: None,
            finished: false,
        }
    }

    pub fn next_record(&mut self) -> Result<Option<Record>> {
        if self.finished {
            return Ok(None);
        }

        let header = match self.pending_header.take() {
            Some(header) => header,
            None => loop {
                let mut line = String::new();
                let n = self.reader.read_line(&mut line)?;
                if n == 0 {
                    self.finished = true;
                    return Ok(None);
                }
                if line.starts_with('>') {
                    break line;
                }
                if !line.trim().is_empty() {
                    bail!("FASTA input must start records with '>'");
                }
            },
        };

        let mut seq = Vec::new();
        loop {
            let mut line = String::new();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                self.finished = true;
                break;
            }
            if line.starts_with('>') {
                self.pending_header = Some(line);
                break;
            }
            seq.extend(
                line.trim_ascii_end()
                    .bytes()
                    .filter(|b| !b.is_ascii_whitespace())
                    .map(normalize_base),
            );
        }

        let id = parse_id(&header[1..]);
        if id.is_empty() {
            bail!("FASTA record has an empty ID");
        }
        Ok(Some(Record { id, seq }))
    }
}

pub struct FastqReader {
    reader: Box<dyn BufRead + Send>,
}

impl FastqReader {
    pub fn new(reader: Box<dyn BufRead + Send>) -> Self {
        Self { reader }
    }

    pub fn next_record(&mut self) -> Result<Option<Record>> {
        let mut name = String::new();
        let n = self.reader.read_line(&mut name)?;
        if n == 0 {
            return Ok(None);
        }
        if !name.starts_with('@') {
            bail!("FASTQ record must start with '@'");
        }
        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();
        if self.reader.read_line(&mut seq)? == 0
            || self.reader.read_line(&mut plus)? == 0
            || self.reader.read_line(&mut qual)? == 0
        {
            bail!("truncated FASTQ record {}", parse_id(&name[1..]));
        }
        if !plus.starts_with('+') {
            bail!(
                "FASTQ record {} missing '+' separator",
                parse_id(&name[1..])
            );
        }
        let seq = seq
            .trim_ascii_end()
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .map(normalize_base)
            .collect::<Vec<_>>();
        let qual_len = qual.trim_ascii_end().len();
        if seq.len() != qual_len {
            bail!(
                "FASTQ record {} has {} bases but {} quality values",
                parse_id(&name[1..]),
                seq.len(),
                qual_len
            );
        }
        Ok(Some(Record {
            id: parse_id(&name[1..]),
            seq,
        }))
    }
}

/// Estimate a sample's mean predicted per-base accuracy from the head of a
/// FASTQ file.
///
/// Scans up to `max_reads` records (or all of them if the file holds fewer),
/// stopping early once `max_bases` quality bytes have been tallied so a few
/// pathologically long reads cannot stall the scan. Returns `Ok(None)` for
/// FASTA input (no qualities), empty input, or a file with no quality bytes;
/// the caller then keeps its fixed cutoff.
///
/// For each quality `Q`, the predicted error probability is `10^(-Q/10)`.
/// Expected base accuracy is therefore `1 - mean(error_probability)`, not the
/// accuracy obtained by transforming the median quality. A 256-bin histogram
/// keeps the scan O(bases) time and O(1) memory, while the comparatively
/// expensive exponentiation is performed only once per quality value.
pub fn sample_mean_base_accuracy(
    path: &str,
    max_reads: usize,
    max_bases: u64,
) -> Result<Option<f64>> {
    let mut reader = open_reader(path)?;
    mean_base_accuracy_from_reader(&mut *reader, max_reads, max_bases)
}

fn mean_base_accuracy_from_reader(
    reader: &mut dyn BufRead,
    max_reads: usize,
    max_bases: u64,
) -> Result<Option<f64>> {
    // Only FASTQ carries qualities; peek the first byte without consuming it.
    let is_fastq = reader.fill_buf()?.first().is_some_and(|&byte| byte == b'@');
    if !is_fastq {
        return Ok(None);
    }

    // Valid Phred+33 quality bytes span '!' (Q0) through '~' (Q93); anything
    // else (newlines, stray whitespace) is ignored by this best-effort scan.
    let mut histogram = [0u64; 256];
    let mut total: u64 = 0;
    let mut line = String::new();
    let mut reads_seen = 0usize;
    while reads_seen < max_reads && total < max_bases {
        // A FASTQ record is exactly four lines: name / seq / '+' / quality.
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break; // clean end of input
        }
        if !line.starts_with('@') {
            break; // desynchronized or malformed: stop the best-effort scan
        }
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        line.clear();
        if reader.read_line(&mut line)? == 0 || !line.starts_with('+') {
            break;
        }
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        for &byte in line.as_bytes() {
            if total >= max_bases {
                break;
            }
            if (b'!'..=b'~').contains(&byte) {
                histogram[byte as usize] += 1;
                total += 1;
            }
        }
        reads_seen += 1;
    }

    if total == 0 {
        return Ok(None);
    }
    Ok(Some(histogram_mean_base_accuracy(&histogram, total)))
}

/// Convert a Phred+33 histogram to expected base accuracy. `total` must equal
/// the histogram sum over valid quality bytes and be greater than zero.
fn histogram_mean_base_accuracy(histogram: &[u64; 256], total: u64) -> f64 {
    let mut error_probability_sum = 0.0f64;
    for byte in b'!'..=b'~' {
        let count = histogram[byte as usize];
        if count != 0 {
            let phred = f64::from(byte - b'!');
            error_probability_sum += count as f64 * 10f64.powf(-phred / 10.0);
        }
    }
    (1.0 - error_probability_sum / total as f64).clamp(0.0, 1.0)
}

fn parse_id(header: &str) -> String {
    let mut id = header
        .split_ascii_whitespace()
        .next()
        .unwrap_or("")
        .to_string();
    if id.ends_with("/1") || id.ends_with("/2") {
        id.truncate(id.len() - 2);
    }
    id
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn accuracy(fastq: &str, max_reads: usize) -> Option<f64> {
        mean_base_accuracy_from_reader(&mut Cursor::new(fastq.as_bytes()), max_reads, u64::MAX)
            .expect("scan of well-formed input should not error")
    }

    #[test]
    fn mean_accuracy_reads_phred33_qualities() {
        // 'I' = ASCII 73 => Q40 => accuracy 1 - 1e-4.
        let acc = accuracy("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n", 10_000).unwrap();
        assert!((acc - (1.0 - 1e-4)).abs() < 1e-9, "acc={acc}");
    }

    #[test]
    fn mean_accuracy_averages_error_probabilities_not_phred_scores() {
        // Five Q10 bases have error probability 0.1 and five Q40 bases have
        // error probability 0.0001. The expected accuracy is therefore
        // 1 - mean(0.1, 0.0001) = 0.94995. Transforming the median quality
        // would incorrectly return either 0.90 or 0.9999.
        let acc = accuracy("@a\nAAAAA\n+\n+++++\n@b\nCCCCC\n+\nIIIII\n", 10_000).unwrap();
        assert!((acc - 0.94995).abs() < 1e-9, "acc={acc}");
    }

    #[test]
    fn mean_accuracy_is_none_for_fasta_and_empty() {
        assert!(accuracy(">r1\nACGT\n", 10).is_none());
        assert!(
            mean_base_accuracy_from_reader(&mut Cursor::new(&b""[..]), 10, u64::MAX)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn mean_accuracy_honors_the_read_cap() {
        // Record 1 is Q40; records 2-3 are '#' (ASCII 35 => Q2). Capping at one
        // read sees only the high-quality head, ignoring the low-quality tail.
        let acc = accuracy("@a\nAA\n+\nII\n@b\nAA\n+\n##\n@c\nAA\n+\n##\n", 1).unwrap();
        assert!((acc - (1.0 - 1e-4)).abs() < 1e-9, "acc={acc}");
    }

    #[test]
    fn mean_accuracy_honors_the_base_cap() {
        // Two bases per record; a two-base cap stops after record 1 (Q40),
        // before the low-quality record 2 can shift the median.
        let acc = mean_base_accuracy_from_reader(
            &mut Cursor::new(&b"@a\nAA\n+\nII\n@b\nAA\n+\n##\n"[..]),
            10_000,
            2,
        )
        .unwrap()
        .unwrap();
        assert!((acc - (1.0 - 1e-4)).abs() < 1e-9, "acc={acc}");
    }

    #[test]
    fn histogram_mean_accuracy_matches_weighted_error_probability() {
        let mut hist = [0u64; 256];
        hist[b'+' as usize] = 1; // Q10, error 0.1
        hist[b'I' as usize] = 3; // Q40, error 0.0001
        let expected = 1.0 - (0.1 + 3.0 * 0.0001) / 4.0;
        let observed = histogram_mean_base_accuracy(&hist, 4);
        assert!((observed - expected).abs() < 1e-12, "observed={observed}");
    }

    #[test]
    fn sample_mean_base_accuracy_reads_from_a_path() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fastq");
        // '+' => Q10 => 0.90.
        std::fs::write(&path, "@a\nACGT\n+\n++++\n@b\nACGT\n+\n++++\n").unwrap();
        let acc = sample_mean_base_accuracy(path.to_str().unwrap(), 10_000, u64::MAX)
            .unwrap()
            .unwrap();
        assert!((acc - 0.90).abs() < 1e-9, "acc={acc}");
    }
}
