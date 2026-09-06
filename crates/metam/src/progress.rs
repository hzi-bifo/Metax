use std::{
    fmt,
    io::{self, IsTerminal},
    time::{Instant, SystemTime, UNIX_EPOCH},
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressMode {
    Auto,
    On,
    Off,
}

#[derive(Clone, Copy, Debug)]
pub struct ProgressConfig {
    pub mode: ProgressMode,
    pub quiet: bool,
    pub interval: u64,
}

pub struct ProgressReporter {
    config: ProgressConfig,
    label: &'static str,
    phase: &'static str,
    last_count: u64,
    start: Instant,
}

impl ProgressReporter {
    pub fn new(config: ProgressConfig, label: &'static str) -> Self {
        Self {
            config,
            label,
            phase: "starting",
            last_count: 0,
            start: Instant::now(),
        }
    }

    pub fn enabled(&self) -> bool {
        if self.config.quiet {
            return false;
        }
        match self.config.mode {
            ProgressMode::Off => false,
            ProgressMode::On => true,
            ProgressMode::Auto => io::stderr().is_terminal(),
        }
    }

    pub fn phase(&mut self, phase: &'static str) {
        self.phase = phase;
        if self.enabled() {
            log_message(self.label, format_args!("{phase}"));
        }
    }

    pub fn tick(&mut self, count: u64) {
        if !self.enabled() {
            return;
        }
        let interval = self.config.interval.max(1);
        if count == 0 || count.saturating_sub(self.last_count) >= interval {
            self.last_count = count;
            let elapsed = self.start.elapsed().as_secs_f64();
            let rate = if elapsed > 0.0 {
                count as f64 / elapsed
            } else {
                0.0
            };
            log_message(
                self.label,
                format_args!("{}: {} processed ({:.1}/s)", self.phase, count, rate),
            );
        }
    }

    pub fn finish(&self, count: u64) {
        if self.enabled() {
            let elapsed = self.start.elapsed().as_secs_f64();
            let rate = if elapsed > 0.0 {
                count as f64 / elapsed
            } else {
                0.0
            };
            log_message(
                self.label,
                format_args!(
                    "done: {} processed in {:.2}s ({:.1}/s)",
                    count, elapsed, rate
                ),
            );
        }
    }
}

pub fn log_message(label: &str, message: fmt::Arguments<'_>) {
    eprintln!("[{}] [{}] {}", utc_timestamp(), label, message);
}

fn utc_timestamp() -> String {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();
    format_unix_timestamp_utc(now.as_secs(), now.subsec_millis())
}

fn format_unix_timestamp_utc(secs: u64, millis: u32) -> String {
    let days = (secs / 86_400) as i64;
    let seconds_of_day = secs % 86_400;
    let hour = seconds_of_day / 3_600;
    let minute = (seconds_of_day % 3_600) / 60;
    let second = seconds_of_day % 60;
    let (year, month, day) = civil_from_days(days);
    format!("{year:04}-{month:02}-{day:02} {hour:02}:{minute:02}:{second:02}.{millis:03}Z")
}

fn civil_from_days(days_since_unix_epoch: i64) -> (i32, u32, u32) {
    let z = days_since_unix_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    let year = y + i64::from(month <= 2);
    (year as i32, month as u32, day as u32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unix_timestamp_format_is_utc_with_milliseconds() {
        assert_eq!(format_unix_timestamp_utc(0, 0), "1970-01-01 00:00:00.000Z");
        assert_eq!(
            format_unix_timestamp_utc(1_704_067_199, 123),
            "2023-12-31 23:59:59.123Z"
        );
        assert_eq!(
            format_unix_timestamp_utc(1_704_067_200, 7),
            "2024-01-01 00:00:00.007Z"
        );
    }
}
