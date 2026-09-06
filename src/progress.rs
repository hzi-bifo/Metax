use std::io::{self, IsTerminal};
use std::time::Duration;

use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

use crate::logging::utc_timestamp_millis;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProgressOptions {
    pub force: bool,
    pub disable: bool,
    pub interval: u64,
}

impl Default for ProgressOptions {
    fn default() -> Self {
        Self {
            force: false,
            disable: false,
            interval: 100_000,
        }
    }
}

impl ProgressOptions {
    pub fn new(force: bool, disable: bool, interval: u64) -> Self {
        Self {
            force,
            disable,
            interval: interval.max(1),
        }
    }

    pub fn disabled() -> Self {
        Self {
            force: false,
            disable: true,
            interval: 1,
        }
    }

    pub fn enabled(self) -> bool {
        if self.disable {
            return false;
        }
        self.force || io::stderr().is_terminal()
    }

    pub fn interval(self) -> u64 {
        self.interval.max(1)
    }
}

pub fn spinner(options: ProgressOptions, message: impl Into<String>) -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    configure_draw_target(&pb, options);
    refresh_timestamp(&pb);
    if options.enabled() {
        pb.enable_steady_tick(Duration::from_millis(120));
    }
    if let Ok(style) =
        ProgressStyle::with_template("{prefix} {spinner:.green} {msg} [{elapsed_precise}]")
    {
        pb.set_style(style);
    }
    pb.set_message(message.into());
    pb
}

pub fn counter(options: ProgressOptions, message: impl Into<String>) -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    configure_draw_target(&pb, options);
    refresh_timestamp(&pb);
    if options.enabled() {
        pb.enable_steady_tick(Duration::from_millis(120));
    }
    if let Ok(style) = ProgressStyle::with_template(
        "{prefix} {spinner:.green} {msg}: {pos} [{elapsed_precise}] {per_sec}",
    ) {
        pb.set_style(style);
    }
    pb.set_message(message.into());
    pb
}

pub fn bar(options: ProgressOptions, len: u64, message: impl Into<String>) -> ProgressBar {
    let pb = ProgressBar::new(len);
    configure_draw_target(&pb, options);
    refresh_timestamp(&pb);
    if let Ok(style) = ProgressStyle::with_template(
        "{prefix} {spinner:.green} {msg} {bar:40.cyan/blue} {pos}/{len} [{elapsed_precise}<{eta_precise}] {per_sec}",
    ) {
        pb.set_style(style.progress_chars("#>-"));
    }
    pb.set_message(message.into());
    pb
}

pub fn bytes_bar(options: ProgressOptions, len: u64, message: impl Into<String>) -> ProgressBar {
    let pb = ProgressBar::new(len);
    configure_draw_target(&pb, options);
    refresh_timestamp(&pb);
    if let Ok(style) = ProgressStyle::with_template(
        "{prefix} {msg} {bar:32.cyan/blue} {bytes}/{total_bytes} ({percent}%) [{elapsed_precise}<{eta_precise}] {bytes_per_sec}",
    ) {
        pb.set_style(style.progress_chars("#>-"));
    }
    pb.set_message(message.into());
    pb
}

pub fn tick_count(pb: &ProgressBar, label: &str, count: u64, interval: u64) {
    let interval = interval.max(1);
    let current = pb.position();
    if count == 0 || count.saturating_sub(current) >= interval {
        refresh_timestamp(pb);
        pb.set_message(label.to_string());
        pb.set_position(count);
        pb.tick();
    }
}

/// Update a live bar whose whole line is a caller-formatted `message`. Unlike
/// [`tick_count`] this refreshes on every call — callers already tick at a
/// coarse (per-batch) cadence, so even a slow run keeps a moving, current line.
pub fn tick_message(pb: &ProgressBar, message: impl Into<String>, count: u64) {
    refresh_timestamp(pb);
    pb.set_message(message.into());
    pb.set_position(count);
    pb.tick();
}

pub fn refresh_timestamp(pb: &ProgressBar) {
    pb.set_prefix(format!("[{}]", utc_timestamp_millis()));
}

fn configure_draw_target(pb: &ProgressBar, options: ProgressOptions) {
    if options.enabled() {
        pb.set_draw_target(ProgressDrawTarget::stderr());
    } else {
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tick_count_updates_on_elapsed_delta_not_exact_multiple() {
        let pb = counter(ProgressOptions::disabled(), "records");

        tick_count(&pb, "records", 65_536, 100_000);
        assert_eq!(pb.position(), 0);

        tick_count(&pb, "records", 131_072, 100_000);
        assert_eq!(pb.position(), 131_072);

        tick_count(&pb, "records", 196_608, 100_000);
        assert_eq!(pb.position(), 131_072);

        tick_count(&pb, "records", 262_144, 100_000);
        assert_eq!(pb.position(), 262_144);
    }
}
