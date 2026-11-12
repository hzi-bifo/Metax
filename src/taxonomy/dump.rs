//! Discovery and managed installation of the NCBI taxonomy dump.

use super::{dump_file_path, Taxonomy};
use anyhow::{anyhow, bail, Context, Result};
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use std::collections::HashSet;
use std::ffi::OsStr;
use std::fs::{self, File, OpenOptions, TryLockError};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant, UNIX_EPOCH};
use tempfile::Builder;

const NCBI_TAXDUMP_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
pub(super) const REQUIRED_DUMP_FILES: [&str; 4] =
    ["nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"];
const MAX_COMPRESSED_BYTES: u64 = 512 * 1024 * 1024;
const MAX_DUMP_FILE_BYTES: u64 = 2 * 1024 * 1024 * 1024;
const MAX_TOTAL_DUMP_BYTES: u64 = 6 * 1024 * 1024 * 1024;
const IO_BUFFER_BYTES: usize = 1024 * 1024;
const PROGRESS_INTERVAL: Duration = Duration::from_secs(10);
const DOWNLOAD_RESOLVE_TIMEOUT: Duration = Duration::from_secs(30);
const DOWNLOAD_CONNECT_TIMEOUT: Duration = Duration::from_secs(30);
const DOWNLOAD_RESPONSE_TIMEOUT: Duration = Duration::from_secs(60);
const DOWNLOAD_BODY_TIMEOUT: Duration = Duration::from_secs(30 * 60);
const BACKUP_DIRECTORY_NAME: &str = ".taxdmp-backup";
const STAGING_DIRECTORY_PREFIX: &str = ".taxdmp-download-";
/// Folder of a database directory that holds the database's taxonomy dump.
pub(crate) const DATABASE_TAXDMP_DIR: &str = "taxdmp";
/// Archive of the dump files in a database directory, read without
/// unpacking; used when the database has no `taxdmp` folder.
pub(crate) const DATABASE_TAXDMP_ARCHIVE: &str = "taxdmp.tar.gz";

/// Where a taxonomy dump was found.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum DumpLocation {
    /// Given with --dmp-dir.
    Explicit,
    /// Stored in the database directory.
    Database,
    /// The managed NCBI dump in the user's home directory.
    Managed,
}

/// A resolved taxonomy location: a folder of dump files or, for a database,
/// an archive of them. Managed locations may be repaired on load; explicit
/// and database locations are read exactly as found and never modified.
#[derive(Clone, Debug)]
pub struct TaxonomySource {
    path: PathBuf,
    location: DumpLocation,
    archive: bool,
}

impl TaxonomySource {
    /// Use `dmp_dir` when given, else the dump stored in `database`, else
    /// the managed NCBI dump, which is downloaded when missing.
    pub(crate) fn prepare(dmp_dir: Option<&str>, database: Option<&Path>) -> Result<Self> {
        if let Some(source) = Self::locate(dmp_dir, database)? {
            return Ok(source);
        }
        let home = dirs::home_dir()
            .context("failed to determine the user home directory; specify --dmp-dir explicitly")?;
        let directory = home.join(".metax").join("taxdmp");
        Self::prepare_managed_with(directory, |target| install_managed_taxdump(target, false))
    }

    /// The explicit or database dump, or `None` when the managed dump is
    /// needed. A database's `taxdmp` folder comes before its archive. A
    /// database dump that is present but incomplete or unreadable is an error,
    /// so a damaged database never falls back to another taxonomy.
    fn locate(dmp_dir: Option<&str>, database: Option<&Path>) -> Result<Option<Self>> {
        if let Some(directory) = dmp_dir {
            return Ok(Some(Self {
                path: PathBuf::from(directory),
                location: DumpLocation::Explicit,
                archive: false,
            }));
        }
        let Some(database) = database else {
            return Ok(None);
        };
        let directory = database.join(DATABASE_TAXDMP_DIR);
        if fs::symlink_metadata(&directory).is_err() {
            let archive = database.join(DATABASE_TAXDMP_ARCHIVE);
            if fs::symlink_metadata(&archive).is_err() {
                return Ok(None);
            }
            log::info!(
                target: "TAXONOMY",
                "Using the taxonomy dump stored with the database in {}.",
                archive.display()
            );
            return Ok(Some(Self {
                path: archive,
                location: DumpLocation::Database,
                archive: true,
            }));
        }
        let missing = REQUIRED_DUMP_FILES
            .iter()
            .filter(|name| !dump_file_path(&directory, name).is_file())
            .copied()
            .collect::<Vec<_>>();
        if !missing.is_empty() {
            bail!(
                "the taxonomy dump stored with the database in {} lacks {}; restore it or pass --dmp-dir",
                directory.display(),
                missing.join(", ")
            );
        }
        log::info!(
            target: "TAXONOMY",
            "Using the taxonomy dump stored with the database in {}.",
            directory.display()
        );
        Ok(Some(Self {
            path: directory,
            location: DumpLocation::Database,
            archive: false,
        }))
    }

    fn prepare_managed_with<F>(directory: PathBuf, installer: F) -> Result<Self>
    where
        F: FnOnce(&Path) -> Result<()>,
    {
        if !dump_files_look_complete(&directory) {
            installer(&directory).with_context(|| {
                format!(
                    "failed to install the NCBI taxonomy dump in {}",
                    directory.display()
                )
            })?;
        }
        Ok(Self {
            path: directory,
            location: DumpLocation::Managed,
            archive: false,
        })
    }

    /// The dump folder, or the database's dump archive.
    pub(crate) fn path(&self) -> &Path {
        &self.path
    }

    /// Whether the dump is the one stored with the database.
    pub(crate) fn is_database_dump(&self) -> bool {
        self.location == DumpLocation::Database
    }

    pub(crate) fn load(&self) -> Result<Taxonomy> {
        self.load_with_repair(|target| install_managed_taxdump(target, true))
    }

    fn load_with_repair<F>(&self, repair: F) -> Result<Taxonomy>
    where
        F: FnOnce(&Path) -> Result<()>,
    {
        let loaded = if self.archive {
            Taxonomy::from_archive(&self.path)
        } else {
            Taxonomy::from_dir(&self.path)
        };
        match loaded {
            Ok(taxonomy) => Ok(taxonomy),
            Err(error) if self.location != DumpLocation::Managed => Err(error).with_context(|| {
                format!(
                    "failed to load the taxonomy dump from {}",
                    self.path.display()
                )
            }),
            Err(error) => {
                log::warn!(
                    target: "TAXONOMY",
                    "Cached taxonomy dump at {} is invalid ({}); downloading a clean copy.",
                    self.path.display(),
                    error
                );
                repair(&self.path).with_context(|| {
                    format!(
                        "failed to replace the invalid taxonomy dump in {}",
                        self.path.display()
                    )
                })?;
                Taxonomy::from_dir(&self.path).with_context(|| {
                    format!(
                        "downloaded taxonomy dump in {} failed validation",
                        self.path.display()
                    )
                })
            }
        }
    }
}

/// Pack the dump files in `source` into the new `.tar.gz` archive `target`,
/// each under its own name and in the form it has there (plain or gzipped).
/// Returns the name, size and SHA-256 of each packed file.
pub(crate) fn pack_dump(source: &Path, target: &Path) -> Result<Vec<(String, u64, String)>> {
    let output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(target)
        .with_context(|| format!("failed to create {}", target.display()))?;
    let encoder = GzEncoder::new(
        BufWriter::with_capacity(IO_BUFFER_BYTES, output),
        Compression::default(),
    );
    let mut archive = tar::Builder::new(encoder);
    let mut packed = Vec::with_capacity(REQUIRED_DUMP_FILES.len());
    for name in REQUIRED_DUMP_FILES {
        let from = dump_file_path(source, name);
        let file =
            File::open(&from).with_context(|| format!("failed to open {}", from.display()))?;
        let metadata = file.metadata()?;
        let mut header = tar::Header::new_gnu();
        header.set_size(metadata.len());
        header.set_mode(0o644);
        header.set_mtime(
            metadata
                .modified()
                .ok()
                .and_then(|time| time.duration_since(UNIX_EPOCH).ok())
                .map_or(0, |age| age.as_secs()),
        );
        let entry = from.file_name().expect("a dump file path ends in its name");
        let mut reader = Sha256Reader {
            inner: BufReader::with_capacity(IO_BUFFER_BYTES, file),
            digest: ring::digest::Context::new(&ring::digest::SHA256),
        };
        archive
            .append_data(&mut header, entry, &mut reader)
            .with_context(|| {
                format!(
                    "failed to pack {} into {}",
                    from.display(),
                    target.display()
                )
            })?;
        let sha256 = reader
            .digest
            .finish()
            .as_ref()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        packed.push((entry.to_string_lossy().into_owned(), metadata.len(), sha256));
    }
    archive
        .into_inner()
        .and_then(|encoder| encoder.finish())
        .and_then(|mut writer| writer.flush())
        .with_context(|| format!("failed to finish {}", target.display()))?;
    Ok(packed)
}

/// A reader that hashes what it reads.
struct Sha256Reader<R> {
    inner: R,
    digest: ring::digest::Context,
}

impl<R: Read> Read for Sha256Reader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.digest.update(&buffer[..read]);
        Ok(read)
    }
}

fn dump_files_look_complete(directory: &Path) -> bool {
    REQUIRED_DUMP_FILES.iter().all(|name| {
        directory
            .join(name)
            .metadata()
            .map(|metadata| metadata.is_file() && metadata.len() > 0)
            .unwrap_or(false)
    })
}

fn install_managed_taxdump(target: &Path, force: bool) -> Result<()> {
    let parent = target
        .parent()
        .context("taxonomy cache path has no parent")?;
    fs::create_dir_all(parent).with_context(|| format!("failed to create {}", parent.display()))?;

    let lock_path = parent.join(".taxdmp.lock");
    let lock_file = OpenOptions::new()
        .create(true)
        // The file is only an inode on which to hold the advisory lock.
        // Its contents are irrelevant and must not be truncated by a waiter.
        .truncate(false)
        .read(true)
        .write(true)
        .open(&lock_path)
        .with_context(|| format!("failed to open install lock {}", lock_path.display()))?;
    match lock_file.try_lock() {
        Ok(()) => {}
        Err(TryLockError::WouldBlock) => {
            log::info!(
                target: "TAXONOMY",
                "Waiting for another Metax process to finish installing the taxonomy dump."
            );
            lock_file
                .lock()
                .with_context(|| format!("failed to lock {}", lock_path.display()))?;
        }
        Err(TryLockError::Error(error)) => {
            return Err(error).with_context(|| format!("failed to lock {}", lock_path.display()));
        }
    }
    cleanup_interrupted_install(parent, target)?;

    // Another process may have completed the installation while this process
    // waited for the lock.
    if !force && dump_files_look_complete(target) {
        return Ok(());
    }

    log::info!(
        target: "TAXONOMY",
        "Downloading the NCBI taxonomy dump to {}.",
        target.display()
    );
    let started = Instant::now();
    let download_config = ureq::Agent::config_builder()
        .timeout_resolve(Some(DOWNLOAD_RESOLVE_TIMEOUT))
        .timeout_connect(Some(DOWNLOAD_CONNECT_TIMEOUT))
        .timeout_recv_response(Some(DOWNLOAD_RESPONSE_TIMEOUT))
        .timeout_recv_body(Some(DOWNLOAD_BODY_TIMEOUT))
        .build();
    let download_agent: ureq::Agent = download_config.into();
    let mut response = download_agent
        .get(NCBI_TAXDUMP_URL)
        .call()
        .map_err(|error| anyhow!("failed to download {NCBI_TAXDUMP_URL}: {error}"))?;
    let content_length = response.body().content_length();
    if content_length.is_some_and(|length| length > MAX_COMPRESSED_BYTES) {
        bail!(
            "refusing taxonomy archive larger than {} MiB",
            MAX_COMPRESSED_BYTES / (1024 * 1024)
        );
    }

    let reader = response.body_mut().as_reader();
    let reader = DownloadReader::new(reader, content_length, MAX_COMPRESSED_BYTES);
    install_archive(target, reader)?;
    log::info!(
        target: "TAXONOMY",
        "NCBI taxonomy dump installed in {:.2}s.",
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Extract, validate, and atomically install one streamed `taxdump.tar.gz`.
fn install_archive<R: Read>(target: &Path, compressed: R) -> Result<()> {
    let parent = target
        .parent()
        .context("taxonomy cache path has no parent")?;
    fs::create_dir_all(parent).with_context(|| format!("failed to create {}", parent.display()))?;
    let staging = Builder::new()
        .prefix(STAGING_DIRECTORY_PREFIX)
        .tempdir_in(parent)
        .with_context(|| format!("failed to create staging directory in {}", parent.display()))?;

    extract_archive(compressed, staging.path())?;
    drop(
        Taxonomy::from_dir(staging.path())
            .context("downloaded NCBI taxonomy dump failed validation")?,
    );
    replace_directory(staging.path(), target)?;
    Ok(())
}

fn extract_archive<R: Read>(compressed: R, destination: &Path) -> Result<()> {
    let buffered = BufReader::with_capacity(IO_BUFFER_BYTES, compressed);
    let decoder = MultiGzDecoder::new(buffered);
    let mut archive = tar::Archive::new(decoder);
    let mut found = HashSet::with_capacity(REQUIRED_DUMP_FILES.len());
    let mut total_size = 0u64;

    for entry in archive
        .entries()
        .context("failed to read taxonomy tar archive")?
    {
        let mut entry = entry.context("failed to read an entry from the taxonomy archive")?;
        if !entry.header().entry_type().is_file() {
            continue;
        }
        let path = entry
            .path()
            .context("taxonomy archive contains an invalid path")?;
        let Some(file_name) = path.file_name() else {
            continue;
        };
        let Some(&expected_name) = REQUIRED_DUMP_FILES
            .iter()
            .find(|name| file_name == OsStr::new(name))
        else {
            continue;
        };
        if !found.insert(expected_name) {
            bail!("taxonomy archive contains duplicate {expected_name}");
        }

        let size = entry.size();
        if size == 0 || size > MAX_DUMP_FILE_BYTES {
            bail!("invalid size for {expected_name} in taxonomy archive: {size} bytes");
        }
        total_size = total_size
            .checked_add(size)
            .context("taxonomy archive size overflow")?;
        if total_size > MAX_TOTAL_DUMP_BYTES {
            bail!(
                "uncompressed taxonomy dump exceeds {} GiB",
                MAX_TOTAL_DUMP_BYTES / (1024 * 1024 * 1024)
            );
        }

        let output_path = destination.join(expected_name);
        let file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&output_path)
            .with_context(|| format!("failed to create {}", output_path.display()))?;
        let mut writer = BufWriter::with_capacity(IO_BUFFER_BYTES, file);
        let copied = io::copy(&mut entry, &mut writer)
            .with_context(|| format!("failed to extract {expected_name}"))?;
        if copied != size {
            bail!(
                "truncated {expected_name} in taxonomy archive: expected {size} bytes, extracted {copied}"
            );
        }
        writer.flush()?;
        writer.get_ref().sync_all()?;
    }

    let missing = REQUIRED_DUMP_FILES
        .iter()
        .filter(|name| !found.contains(**name))
        .copied()
        .collect::<Vec<_>>();
    if !missing.is_empty() {
        bail!(
            "taxonomy archive is missing required file(s): {}",
            missing.join(", ")
        );
    }
    Ok(())
}

fn replace_directory(staging: &Path, target: &Path) -> Result<()> {
    let parent = target
        .parent()
        .context("taxonomy cache path has no parent")?;
    let backup = parent.join(BACKUP_DIRECTORY_NAME);
    remove_path_if_present(&backup)?;

    let had_target = target.exists();
    if had_target {
        fs::rename(target, &backup).with_context(|| {
            format!(
                "failed to move existing taxonomy dump {} aside",
                target.display()
            )
        })?;
    }

    if let Err(error) = fs::rename(staging, target) {
        if had_target {
            if let Err(restore_error) = fs::rename(&backup, target) {
                return Err(anyhow!(
                    "failed to install taxonomy dump ({error}) and failed to restore the previous dump ({restore_error})"
                ));
            }
        }
        return Err(error)
            .with_context(|| format!("failed to install taxonomy dump at {}", target.display()));
    }

    if had_target {
        if let Err(error) = remove_path_if_present(&backup) {
            log::warn!(
                target: "TAXONOMY",
                "Installed the taxonomy dump but could not remove backup {}: {}",
                backup.display(),
                error
            );
        }
    }
    Ok(())
}

/// Recover the narrow rename window from an interrupted installation and
/// remove staging directories left by a terminated process. The install lock
/// guarantees that none of these paths belongs to a live installer.
fn cleanup_interrupted_install(parent: &Path, target: &Path) -> Result<()> {
    let backup = parent.join(BACKUP_DIRECTORY_NAME);
    if backup.exists() {
        if target.exists() {
            remove_path_if_present(&backup)?;
        } else {
            fs::rename(&backup, target).with_context(|| {
                format!(
                    "failed to restore interrupted taxonomy installation from {}",
                    backup.display()
                )
            })?;
        }
    }

    for entry in
        fs::read_dir(parent).with_context(|| format!("failed to inspect {}", parent.display()))?
    {
        let entry = entry?;
        let name = entry.file_name();
        if name
            .to_str()
            .is_some_and(|name| name.starts_with(STAGING_DIRECTORY_PREFIX))
        {
            remove_path_if_present(&entry.path())?;
        }
    }
    Ok(())
}

fn remove_path_if_present(path: &Path) -> Result<()> {
    let metadata = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error.into()),
    };
    if metadata.file_type().is_dir() && !metadata.file_type().is_symlink() {
        fs::remove_dir_all(path)?;
    } else {
        fs::remove_file(path)?;
    }
    Ok(())
}

struct DownloadReader<R> {
    inner: R,
    total: Option<u64>,
    limit: u64,
    bytes_read: u64,
    started: Instant,
    last_report: Instant,
}

impl<R> DownloadReader<R> {
    fn new(inner: R, total: Option<u64>, limit: u64) -> Self {
        let now = Instant::now();
        Self {
            inner,
            total,
            limit,
            bytes_read: 0,
            started: now,
            last_report: now,
        }
    }
}

impl<R: Read> Read for DownloadReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        if self.bytes_read >= self.limit {
            let mut probe = [0u8; 1];
            return match self.inner.read(&mut probe)? {
                0 => Ok(0),
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "taxonomy archive exceeds the compressed-size limit",
                )),
            };
        }
        let remaining = (self.limit - self.bytes_read).min(buffer.len() as u64) as usize;
        let read = self.inner.read(&mut buffer[..remaining])?;
        self.bytes_read += read as u64;

        if self.last_report.elapsed() >= PROGRESS_INTERVAL {
            let elapsed = self.started.elapsed().as_secs_f64().max(f64::EPSILON);
            let mib = self.bytes_read as f64 / (1024.0 * 1024.0);
            let speed = mib / elapsed;
            if let Some(total) = self.total.filter(|total| *total > 0) {
                let percent = self.bytes_read as f64 * 100.0 / total as f64;
                log::info!(
                    target: "TAXONOMY",
                    "Downloading taxonomy dump: {:.1}% ({:.1} MiB, {:.1} MiB/s).",
                    percent.min(100.0),
                    mib,
                    speed
                );
            } else {
                log::info!(
                    target: "TAXONOMY",
                    "Downloading taxonomy dump: {:.1} MiB ({:.1} MiB/s).",
                    mib,
                    speed
                );
            }
            self.last_report = Instant::now();
        }
        Ok(read)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::cell::Cell;
    use std::io::Cursor;

    const NODES: &str = "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tspecies\t|\n";
    const NAMES: &str =
        "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tSpecies two\t|\t\t|\tscientific name\t|\n";
    const MERGED: &str = "3\t|\t2\t|\n";
    const DELNODES: &str = "4\t|\n";

    fn write_dump(directory: &Path) {
        fs::create_dir_all(directory).unwrap();
        fs::write(directory.join("nodes.dmp"), NODES).unwrap();
        fs::write(directory.join("names.dmp"), NAMES).unwrap();
        fs::write(directory.join("merged.dmp"), MERGED).unwrap();
        fs::write(directory.join("delnodes.dmp"), DELNODES).unwrap();
    }

    fn archive(include_all_required: bool) -> Vec<u8> {
        let encoder = GzEncoder::new(Vec::new(), Compression::fast());
        let mut builder = tar::Builder::new(encoder);
        let files = [
            ("nodes.dmp", NODES),
            ("names.dmp", NAMES),
            ("merged.dmp", MERGED),
            ("delnodes.dmp", DELNODES),
            ("readme.txt", "not installed\n"),
        ];
        for (name, contents) in files {
            if !include_all_required && name == "delnodes.dmp" {
                continue;
            }
            let mut header = tar::Header::new_gnu();
            header.set_size(contents.len() as u64);
            header.set_mode(0o644);
            header.set_cksum();
            builder
                .append_data(&mut header, name, contents.as_bytes())
                .unwrap();
        }
        builder.into_inner().unwrap().finish().unwrap()
    }

    #[test]
    fn complete_managed_dump_skips_installation() {
        let temp = tempfile::tempdir().unwrap();
        let directory = temp.path().join("taxdmp");
        write_dump(&directory);
        let called = Cell::new(false);
        let source = TaxonomySource::prepare_managed_with(directory, |_| {
            called.set(true);
            Ok(())
        })
        .unwrap();
        assert!(!called.get());
        assert_eq!(source.load().unwrap().get_name(2), Some("Species two"));
    }

    #[test]
    fn missing_managed_dump_is_installed() {
        let temp = tempfile::tempdir().unwrap();
        let directory = temp.path().join("taxdmp");
        let called = Cell::new(false);
        let source = TaxonomySource::prepare_managed_with(directory, |target| {
            called.set(true);
            write_dump(target);
            Ok(())
        })
        .unwrap();
        assert!(called.get());
        assert_eq!(source.load().unwrap().get_rank(2), Some("species"));
    }

    #[test]
    fn corrupt_managed_dump_is_repaired_but_explicit_dump_is_not() {
        let temp = tempfile::tempdir().unwrap();
        let directory = temp.path().join("taxdmp");
        fs::create_dir_all(&directory).unwrap();
        for name in REQUIRED_DUMP_FILES {
            fs::write(directory.join(name), "corrupt\n").unwrap();
        }
        let managed = TaxonomySource {
            path: directory.clone(),
            location: DumpLocation::Managed,
            archive: false,
        };
        let repaired = managed
            .load_with_repair(|target| {
                write_dump(target);
                Ok(())
            })
            .unwrap();
        assert_eq!(repaired.get_name(2), Some("Species two"));

        for name in REQUIRED_DUMP_FILES {
            fs::write(directory.join(name), "corrupt\n").unwrap();
        }
        let explicit = TaxonomySource {
            path: directory,
            location: DumpLocation::Explicit,
            archive: false,
        };
        assert!(explicit
            .load_with_repair(|_| panic!("explicit dumps must not be repaired"))
            .is_err());
    }

    #[test]
    fn archive_install_keeps_only_required_validated_files() {
        let temp = tempfile::tempdir().unwrap();
        let target = temp.path().join("taxdmp");
        install_archive(&target, Cursor::new(archive(true))).unwrap();
        assert!(dump_files_look_complete(&target));
        assert!(!target.join("readme.txt").exists());
        assert_eq!(fs::read_dir(&target).unwrap().count(), 4);
        assert_eq!(
            Taxonomy::from_dir(&target).unwrap().get_name(2),
            Some("Species two")
        );
    }

    #[test]
    fn incomplete_archive_does_not_replace_existing_dump() {
        let temp = tempfile::tempdir().unwrap();
        let target = temp.path().join("taxdmp");
        write_dump(&target);
        assert!(install_archive(&target, Cursor::new(archive(false))).is_err());
        assert_eq!(
            Taxonomy::from_dir(&target).unwrap().get_name(2),
            Some("Species two")
        );
    }

    #[test]
    fn interrupted_install_restores_backup_and_removes_staging() {
        let temp = tempfile::tempdir().unwrap();
        let target = temp.path().join("taxdmp");
        let backup = temp.path().join(BACKUP_DIRECTORY_NAME);
        write_dump(&backup);
        let stale_staging = temp.path().join(format!("{STAGING_DIRECTORY_PREFIX}stale"));
        fs::create_dir(&stale_staging).unwrap();

        cleanup_interrupted_install(temp.path(), &target).unwrap();

        assert_eq!(
            Taxonomy::from_dir(&target).unwrap().get_name(2),
            Some("Species two")
        );
        assert!(!backup.exists());
        assert!(!stale_staging.exists());
    }

    #[test]
    fn database_dump_is_used_unless_a_directory_is_given() {
        let temp = tempfile::tempdir().unwrap();
        let database = temp.path().join("db");
        let stored = database.join(DATABASE_TAXDMP_DIR);
        write_dump(&stored);
        let explicit = temp.path().join("explicit");
        write_dump(&explicit);

        let source = TaxonomySource::locate(None, Some(&database))
            .unwrap()
            .unwrap();
        assert_eq!(source.location, DumpLocation::Database);
        assert_eq!(source.path, stored);
        assert_eq!(source.load().unwrap().get_name(2), Some("Species two"));

        let source = TaxonomySource::locate(explicit.to_str(), Some(&database))
            .unwrap()
            .unwrap();
        assert_eq!(source.location, DumpLocation::Explicit);
        assert_eq!(source.path, explicit);
    }

    #[test]
    fn database_without_a_dump_needs_the_managed_dump() {
        let temp = tempfile::tempdir().unwrap();
        assert!(TaxonomySource::locate(None, Some(temp.path()))
            .unwrap()
            .is_none());
        assert!(TaxonomySource::locate(None, None).unwrap().is_none());
    }

    #[test]
    fn incomplete_database_dump_is_an_error_not_a_fallback() {
        let temp = tempfile::tempdir().unwrap();
        let stored = temp.path().join(DATABASE_TAXDMP_DIR);
        write_dump(&stored);
        fs::remove_file(stored.join("merged.dmp")).unwrap();
        let error = TaxonomySource::locate(None, Some(temp.path())).unwrap_err();
        assert!(format!("{error:#}").contains("merged.dmp"));
    }

    #[test]
    fn damaged_database_dump_is_not_repaired() {
        let temp = tempfile::tempdir().unwrap();
        let stored = temp.path().join(DATABASE_TAXDMP_DIR);
        write_dump(&stored);
        fs::write(stored.join("nodes.dmp"), "corrupt\n").unwrap();
        let source = TaxonomySource::locate(None, Some(temp.path()))
            .unwrap()
            .unwrap();
        assert!(source
            .load_with_repair(|_| panic!("database dumps must not be repaired"))
            .is_err());
    }

    #[test]
    fn database_archive_is_read_without_unpacking() {
        let temp = tempfile::tempdir().unwrap();
        let archive_path = temp.path().join(DATABASE_TAXDMP_ARCHIVE);
        fs::write(&archive_path, archive(true)).unwrap();

        let source = TaxonomySource::locate(None, Some(temp.path()))
            .unwrap()
            .unwrap();
        assert!(source.is_database_dump());
        assert_eq!(source.path, archive_path);
        assert_eq!(source.load().unwrap().get_name(3), Some("Species two"));
        let names = fs::read_dir(temp.path()).unwrap().count();
        assert_eq!(names, 1, "reading the archive must not write anything");
    }

    #[test]
    fn database_folder_is_used_before_its_archive() {
        let temp = tempfile::tempdir().unwrap();
        let stored = temp.path().join(DATABASE_TAXDMP_DIR);
        write_dump(&stored);
        fs::write(temp.path().join(DATABASE_TAXDMP_ARCHIVE), b"not read").unwrap();
        let source = TaxonomySource::locate(None, Some(temp.path()))
            .unwrap()
            .unwrap();
        assert_eq!(source.path, stored);
    }

    #[test]
    fn damaged_database_archive_is_an_error_not_a_fallback() {
        let temp = tempfile::tempdir().unwrap();
        fs::write(temp.path().join(DATABASE_TAXDMP_ARCHIVE), archive(false)).unwrap();
        let source = TaxonomySource::locate(None, Some(temp.path()))
            .unwrap()
            .unwrap();
        let error = source
            .load_with_repair(|_| panic!("database dumps must not be repaired"))
            .unwrap_err();
        assert!(format!("{error:#}").contains("delnodes.dmp"));
    }
}
