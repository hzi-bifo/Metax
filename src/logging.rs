use std::fs::File;
use std::io::{self, stdout, BufWriter, Write};

use anyhow::Result;
use chrono::Utc;
use env_logger::{fmt::Formatter, Builder, Target};
use log::Level;

struct MultiWriter<StdoutWriter, FileWriter> {
    stdout: StdoutWriter,
    file: FileWriter,
    stdout_available: bool,
}

impl MultiWriter<io::Stdout, BufWriter<File>> {
    fn new(file: File) -> Self {
        Self {
            stdout: stdout(),
            file: BufWriter::new(file),
            stdout_available: true,
        }
    }
}

impl<StdoutWriter: Write, FileWriter: Write> Write for MultiWriter<StdoutWriter, FileWriter> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        // The persistent log is the durable sink. Write it first so a closed
        // stdout pipe cannot suppress the record users need for diagnosis.
        self.file.write_all(buf)?;
        self.file.flush()?;
        if !self.stdout_available {
            return Ok(buf.len());
        }
        if let Err(error) = self.stdout.write_all(buf).and_then(|_| self.stdout.flush()) {
            self.stdout_available = false;
            if error.kind() == io::ErrorKind::BrokenPipe {
                return Ok(buf.len());
            }
            return Err(error);
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()?;
        if !self.stdout_available {
            return Ok(());
        }
        if let Err(error) = self.stdout.flush() {
            self.stdout_available = false;
            if error.kind() == io::ErrorKind::BrokenPipe {
                return Ok(());
            }
            return Err(error);
        }
        Ok(())
    }
}

pub fn init_logging(default_target: &str, logfile: Option<String>) -> Result<()> {
    let mut builder = Builder::new();
    builder.filter_level(log::LevelFilter::Info);
    let default_target_owned = default_target.to_string();
    builder.format(move |buf: &mut Formatter, record| {
        let target = if record.target().is_empty() {
            default_target_owned.as_str()
        } else {
            record.target()
        };
        let level = match record.level() {
            Level::Error => "ERROR",
            Level::Warn => "WARN",
            Level::Info => "INFO",
            Level::Debug => "DEBUG",
            Level::Trace => "TRACE",
        };
        writeln!(
            buf,
            "[{}] [{}] [{}] {}",
            utc_timestamp_millis(),
            target.to_uppercase(),
            level,
            record.args()
        )
    });

    match logfile {
        Some(path) => {
            let file = File::create(path)?;
            builder.target(Target::Pipe(Box::new(MultiWriter::new(file))));
        }
        None => {
            builder.target(Target::Stdout);
        }
    }

    let _ = builder.try_init();
    Ok(())
}

pub fn utc_timestamp_millis() -> String {
    Utc::now().format("%Y-%m-%d %H:%M:%S%.3fZ").to_string()
}

#[cfg(test)]
mod tests {
    use std::io::{self, Write};

    use super::MultiWriter;

    struct BrokenPipeWriter {
        writes: usize,
    }

    impl Write for BrokenPipeWriter {
        fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
            self.writes += 1;
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed pipe"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn broken_stdout_does_not_interrupt_persistent_logging() {
        let mut writer = MultiWriter {
            stdout: BrokenPipeWriter { writes: 0 },
            file: Vec::new(),
            stdout_available: true,
        };

        writer.write_all(b"first\n").unwrap();
        writer.write_all(b"second\n").unwrap();

        assert_eq!(writer.file, b"first\nsecond\n");
        assert_eq!(writer.stdout.writes, 1, "failed stdout must be disabled");
        assert!(!writer.stdout_available);
    }
}
