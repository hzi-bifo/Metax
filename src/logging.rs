use std::fs::File;
use std::io::{self, stdout, BufWriter, Write};

use anyhow::Result;
use chrono::Local;
use env_logger::{fmt::Formatter, Builder, Target};
use log::Level;

struct MultiWriter {
    stdout: io::Stdout,
    file: BufWriter<File>,
}

impl MultiWriter {
    fn new(file: File) -> Self {
        Self {
            stdout: stdout(),
            file: BufWriter::new(file),
        }
    }
}

impl Write for MultiWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.stdout.write_all(buf)?;
        self.stdout.flush()?;
        self.file.write_all(buf)?;
        self.file.flush()?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.stdout.flush()?;
        self.file.flush()
    }
}

pub fn init_logging(default_target: &str, logfile: Option<String>) -> Result<()> {
    let mut builder = Builder::new();
    builder.filter_level(log::LevelFilter::Info);
    let default_target_owned = default_target.to_string();
    builder.format(move |buf: &mut Formatter, record| {
        let ts = Local::now().format("%Y-%m-%d %H:%M:%S,%3f");
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
            "{} {} {} {}",
            ts,
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
