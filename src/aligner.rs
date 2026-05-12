use std::io::{BufRead, BufReader};
use std::process::{Command, ExitStatus, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;

use anyhow::{Context, Result};

use crate::profiler::ProfileMode;

pub struct Aligner {
    db: String,
    sequences: Vec<String>,
    outprefix: String,
}

fn format_exit_status(status: &ExitStatus) -> String {
    // Standardize exit status reporting across platforms.
    match status.code() {
        Some(code) => format!("exit code {}", code),
        None => {
            #[cfg(unix)]
            {
                use std::os::unix::process::ExitStatusExt;
                if let Some(signal) = status.signal() {
                    return format!("signal {}", signal);
                }
            }
            "terminated by signal".to_string()
        }
    }
}

impl Aligner {
    pub fn new(db: String, sequences: Vec<String>, outprefix: String) -> Self {
        Self {
            db,
            sequences,
            outprefix,
        }
    }

    pub fn run(
        &self,
        threads: usize,
        sequencer: &str,
        is_paired: bool,
        mode: ProfileMode,
        log_command: bool,
        extra_args: &[String],
    ) -> Result<()> {
        // Build maCMD CLI arguments and output path.
        let out_file = format!("{}.sam", self.outprefix);
        let mut args = Vec::new();
        args.push("-t".to_string());
        args.push(threads.to_string());
        args.push("-x".to_string());
        args.push(self.db.clone());

        let (preset, input_params) = self.build_input_params(sequencer, is_paired)?;
        args.push("-p".to_string());
        args.push(preset);
        args.extend(input_params);

        args.push("-o".to_string());
        args.push(out_file.clone());

        let mut mode_params = self.mode_parameters(sequencer, mode);
        self.apply_extra_args(&mut mode_params, extra_args);

        for (k, v) in mode_params {
            args.push(k);
            if !v.is_empty() {
                args.push(v);
            }
        }

        let command_str = format!("maCMD {}", args.join(" "));
        if log_command {
            log::info!(target: "ALIGN", "Running command: {}", command_str);
        }

        // Stream stderr for real-time feedback while keeping a copy for errors.
        let mut cmd = Command::new("maCMD");
        cmd.args(&args);
        cmd.stdout(Stdio::inherit());
        cmd.stderr(Stdio::piped());
        let mut child = cmd.spawn().context("failed to spawn maCMD")?;
        let stderr_buf = Arc::new(Mutex::new(String::new()));
        let stderr_handle = if let Some(mut stderr) = child.stderr.take() {
            let buffer = Arc::clone(&stderr_buf);
            Some(thread::spawn(move || {
                let mut reader = BufReader::new(&mut stderr);
                loop {
                    match reader.fill_buf() {
                        Ok(data) if data.is_empty() => break,
                        Ok(data) => {
                            let len = data.len();
                            let text = String::from_utf8_lossy(data);
                            eprint!("{}", text);
                            if let Ok(mut buf) = buffer.lock() {
                                buf.push_str(&text);
                            }
                            reader.consume(len);
                        }
                        Err(_) => break,
                    }
                }
            }))
        } else {
            None
        };
        let status = child.wait().context("failed to wait for maCMD to finish")?;
        if let Some(handle) = stderr_handle {
            let _ = handle.join();
        }
        let stderr_buf = stderr_buf.lock().map(|buf| buf.clone()).unwrap_or_default();
        if !status.success() {
            let status_desc = format_exit_status(&status);
            let stderr = stderr_buf.trim();
            let mut message = format!("maCMD failed ({status_desc})");
            if !stderr.is_empty() {
                message.push_str(&format!("\nstderr:\n{}", stderr));
            }
            anyhow::bail!("{message}\ncommand: {command_str}");
        }
        log::info!(target: "ALIGN", "Alignment finished.");
        Ok(())
    }

    fn build_input_params(
        &self,
        sequencer: &str,
        is_paired: bool,
    ) -> Result<(String, Vec<String>)> {
        // Map sequencer + pairing to maCMD preset and input args.
        if sequencer.eq_ignore_ascii_case("Illumina") {
            if is_paired {
                if self.sequences.len() != 2 {
                    anyhow::bail!("The number of input sequence file should be 2 for Illumina paired-end reads.");
                }
                let params = vec![
                    "-i".to_string(),
                    self.sequences[0].clone(),
                    "-m".to_string(),
                    self.sequences[1].clone(),
                ];
                Ok(("Illumina_Paired".to_string(), params))
            } else {
                let joined = self.sequences.join(",");
                Ok(("Illumina".to_string(), vec!["-i".to_string(), joined]))
            }
        } else {
            let joined = self.sequences.join(",");
            let preset = if sequencer.eq_ignore_ascii_case("Nanopore")
                || sequencer.eq_ignore_ascii_case("PacBio")
            {
                sequencer.to_string()
            } else {
                "Default".to_string()
            };
            Ok((preset, vec!["-i".to_string(), joined]))
        }
    }

    fn mode_parameters(&self, sequencer: &str, mode: ProfileMode) -> Vec<(String, String)> {
        use ProfileMode::*;
        let mut params = Vec::new();
        // Keep aligner scoring aligned with profile sensitivity mode.
        if sequencer.eq_ignore_ascii_case("Illumina") {
            params.push(("-s".into(), "maxSpan".into()));
            match mode {
                Recall => {
                    params.push(("--Minimal_Alignment_Score".into(), "50".into()));
                    params.push(("-M".into(), "5".into()));
                }
                Precision => {
                    params.push(("--Minimal_Alignment_Score".into(), "60".into()));
                    params.push(("-l".into(), "18".into()));
                    params.push(("-M".into(), "2".into()));
                    params.push(("-N".into(), "16".into()));
                }
                Default => {
                    params.push(("--Minimal_Alignment_Score".into(), "55".into()));
                    params.push(("-l".into(), "17".into()));
                    params.push(("-M".into(), "3".into()));
                    params.push(("-N".into(), "12".into()));
                }
            }
        } else {
            params.push(("-s".into(), "maxSpan".into()));
            match mode {
                Recall => {
                    params.push(("-M".into(), "5".into()));
                }
                Precision => {
                    params.push(("-l".into(), "18".into()));
                    params.push(("-M".into(), "2".into()));
                    params.push(("-N".into(), "12".into()));
                }
                Default => {
                    params.push(("-l".into(), "17".into()));
                    params.push(("-M".into(), "3".into()));
                    params.push(("-N".into(), "15".into()));
                }
            }

            if !sequencer.eq_ignore_ascii_case("Nanopore")
                && !sequencer.eq_ignore_ascii_case("PacBio")
                && matches!(mode, ProfileMode::Precision)
            {
                Self::update_param(&mut params, "-l", "19");
                Self::update_param(&mut params, "-M", "2");
                Self::update_param(&mut params, "-N", "18");
            }
        }

        // Request extended CIGAR when available.
        params.push(("--Use_M_in_CIGAR".into(), "false".into()));
        params
    }

    fn apply_extra_args(&self, params: &mut Vec<(String, String)>, extra_args: &[String]) {
        for chunk in extra_args.chunks(2) {
            if chunk.len() != 2 {
                continue;
            }
            let key = chunk[0].clone();
            let value = chunk[1].clone();
            if let Some(existing) = params.iter_mut().find(|(k, _)| k == &key) {
                existing.1 = value;
            } else {
                params.push((key, value));
            }
        }
    }

    fn update_param(params: &mut Vec<(String, String)>, key: &str, value: &str) {
        if let Some(existing) = params.iter_mut().find(|(k, _)| k == key) {
            existing.1 = value.to_string();
        } else {
            params.push((key.to_string(), value.to_string()));
        }
    }
}
