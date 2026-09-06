use std::env;
use std::ffi::CString;
use std::io;
use std::mem::MaybeUninit;
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use tempfile::Builder;

const DEFAULT_TEMP_SPACE_FLOOR_BYTES: u64 = 4 * 1024 * 1024 * 1024;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StagingPlan {
    pub directory: PathBuf,
    pub fallback_directory: Option<PathBuf>,
    pub used_output_fallback: bool,
    pub available_bytes: Option<u64>,
}

pub fn plan_staging_directory(explicit: Option<&str>, outprefix: &str) -> Result<StagingPlan> {
    if let Some(directory) = explicit {
        return plan_explicit_staging_directory(directory, outprefix);
    }
    plan_automatic_staging_directory(outprefix)
}

fn plan_explicit_staging_directory(directory: &str, _outprefix: &str) -> Result<StagingPlan> {
    let directory = PathBuf::from(directory);
    let available = probe_directory(&directory).with_context(|| {
        format!(
            "explicit temporary directory '{}' is unavailable",
            directory.display()
        )
    })?;
    if available.is_some_and(|bytes| bytes < DEFAULT_TEMP_SPACE_FLOOR_BYTES) {
        log::warn!(
            target: "MAPPER",
            "Explicit temporary directory '{}' has only {} available; honoring the explicit choice, but the run may exhaust temporary storage.",
            directory.display(),
            human_bytes(available.unwrap_or(0))
        );
    }
    Ok(StagingPlan {
        directory,
        fallback_directory: None,
        used_output_fallback: false,
        available_bytes: available,
    })
}

fn plan_automatic_staging_directory(outprefix: &str) -> Result<StagingPlan> {
    plan_automatic_staging_directory_with(
        outprefix,
        automatic_temp_candidates(),
        DEFAULT_TEMP_SPACE_FLOOR_BYTES,
    )
}

fn plan_automatic_staging_directory_with(
    outprefix: &str,
    candidates: impl IntoIterator<Item = PathBuf>,
    minimum_free_bytes: u64,
) -> Result<StagingPlan> {
    let output_directory = output_directory(outprefix);
    let mut usable_network_directory = None;
    for directory in candidates {
        let available = match probe_directory(&directory) {
            Ok(available) if available.is_none_or(|bytes| bytes >= minimum_free_bytes) => available,
            Ok(available) => {
                log::debug!(
                    target: "MAPPER",
                    "Skipping temporary directory '{}' because only {} is available.",
                    directory.display(),
                    human_bytes(available.unwrap_or(0))
                );
                continue;
            }
            Err(error) => {
                log::debug!(
                    target: "MAPPER",
                    "Skipping unavailable temporary directory '{}': {}",
                    directory.display(),
                    error
                );
                continue;
            }
        };
        let network_backed = is_network_filesystem(&directory).unwrap_or(false);
        if network_backed {
            usable_network_directory.get_or_insert((directory, available));
            continue;
        }
        return Ok(automatic_plan(directory, available, &output_directory));
    }

    if let Some((directory, available)) = usable_network_directory {
        log::warn!(
            target: "MAPPER",
            "No suitable local temporary directory was found; using network-backed temporary directory '{}'.",
            directory.display()
        );
        return Ok(automatic_plan(directory, available, &output_directory));
    }

    let output_available = probe_directory(&output_directory).with_context(|| {
        format!(
            "no suitable temporary directory was found and output directory '{}' cannot be used",
            output_directory.display()
        )
    })?;
    log::warn!(
        target: "MAPPER",
        "No writable temporary directory with at least {} available was found; staging mappings in output directory '{}'.",
        human_bytes(minimum_free_bytes),
        output_directory.display()
    );
    Ok(StagingPlan {
        directory: output_directory,
        fallback_directory: None,
        used_output_fallback: true,
        available_bytes: output_available,
    })
}

fn automatic_plan(
    directory: PathBuf,
    available_bytes: Option<u64>,
    output_directory: &Path,
) -> StagingPlan {
    let fallback_directory = if same_directory(&directory, output_directory) {
        None
    } else {
        Some(output_directory.to_path_buf())
    };
    StagingPlan {
        directory,
        fallback_directory,
        used_output_fallback: false,
        available_bytes,
    }
}

fn automatic_temp_candidates() -> Vec<PathBuf> {
    let mut candidates = Vec::new();
    if let Some(directory) = env::var_os("TMPDIR") {
        push_distinct_directory(&mut candidates, PathBuf::from(directory));
    }
    push_distinct_directory(&mut candidates, PathBuf::from("/tmp"));
    push_distinct_directory(&mut candidates, env::temp_dir());
    candidates
}

fn push_distinct_directory(directories: &mut Vec<PathBuf>, candidate: PathBuf) {
    if !directories
        .iter()
        .any(|existing| same_directory(existing, &candidate))
    {
        directories.push(candidate);
    }
}

pub fn is_storage_exhausted(error: &anyhow::Error) -> bool {
    error.chain().any(|cause| {
        cause
            .downcast_ref::<io::Error>()
            .and_then(io::Error::raw_os_error)
            .is_some_and(|code| code == libc::ENOSPC || code == libc::EDQUOT)
    })
}

pub fn output_directory(outprefix: &str) -> PathBuf {
    Path::new(outprefix)
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf()
}

pub fn human_bytes(bytes: u64) -> String {
    const GIB: f64 = (1024_u64 * 1024 * 1024) as f64;
    format!("{:.2} GiB", bytes as f64 / GIB)
}

fn probe_directory(directory: &Path) -> Result<Option<u64>> {
    if !directory.is_dir() {
        bail!("'{}' is not an existing directory", directory.display());
    }
    let probe = Builder::new()
        .prefix(".metax-write-probe-")
        .tempfile_in(directory)
        .with_context(|| format!("directory '{}' is not writable", directory.display()))?;
    drop(probe);
    available_space(directory).map(Some).or_else(|error| {
        log::debug!(
            target: "MAPPER",
            "Could not determine free space for '{}': {}",
            directory.display(),
            error
        );
        Ok(None)
    })
}

// The statvfs field widths differ across libc targets, so the casts are
// redundant only on the ones where they already match.
#[allow(clippy::unnecessary_cast)]
fn available_space(directory: &Path) -> io::Result<u64> {
    let path = CString::new(directory.as_os_str().as_bytes())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "path contains a NUL byte"))?;
    let mut stats = MaybeUninit::<libc::statvfs>::uninit();
    let result = unsafe { libc::statvfs(path.as_ptr(), stats.as_mut_ptr()) };
    if result != 0 {
        return Err(io::Error::last_os_error());
    }
    let stats = unsafe { stats.assume_init() };
    Ok((stats.f_bavail as u64).saturating_mul(stats.f_frsize as u64))
}

// f_type is __fsword_t, whose width likewise varies by target.
#[allow(clippy::unnecessary_cast)]
fn is_network_filesystem(directory: &Path) -> io::Result<bool> {
    // Linux filesystem magic values from statfs(2).
    const NFS_SUPER_MAGIC: i64 = 0x6969;
    const CIFS_MAGIC_NUMBER: i64 = 0xFF53_4D42;

    let path = CString::new(directory.as_os_str().as_bytes())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "path contains a NUL byte"))?;
    let mut stats = MaybeUninit::<libc::statfs>::uninit();
    let result = unsafe { libc::statfs(path.as_ptr(), stats.as_mut_ptr()) };
    if result != 0 {
        return Err(io::Error::last_os_error());
    }
    let filesystem_type = unsafe { stats.assume_init() }.f_type as i64;
    Ok(matches!(
        filesystem_type,
        NFS_SUPER_MAGIC | CIFS_MAGIC_NUMBER
    ))
}

fn same_directory(left: &Path, right: &Path) -> bool {
    match (left.canonicalize(), right.canonicalize()) {
        (Ok(left), Ok(right)) => left == right,
        _ => left == right,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn explicit_directory_is_honored() {
        let directory = tempdir().unwrap();
        let plan = plan_explicit_staging_directory(
            directory.path().to_str().unwrap(),
            directory.path().join("sample").to_str().unwrap(),
        )
        .unwrap();
        assert_eq!(plan.directory, directory.path());
        assert!(plan.fallback_directory.is_none());
        assert!(!plan.used_output_fallback);
    }

    #[test]
    fn unavailable_automatic_directory_falls_back_to_output() {
        let output = tempdir().unwrap();
        let missing = output.path().join("missing");
        let outprefix = output.path().join("sample");
        let plan = plan_automatic_staging_directory_with(
            outprefix.to_str().unwrap(),
            [missing],
            DEFAULT_TEMP_SPACE_FLOOR_BYTES,
        )
        .unwrap();
        assert_eq!(plan.directory, output.path());
        assert!(plan.used_output_fallback);
    }

    #[test]
    fn insufficient_automatic_space_falls_back_to_output() {
        let automatic = tempdir().unwrap();
        let output = tempdir().unwrap();
        let outprefix = output.path().join("sample");
        let plan = plan_automatic_staging_directory_with(
            outprefix.to_str().unwrap(),
            [automatic.path().to_path_buf()],
            u64::MAX,
        )
        .unwrap();
        assert_eq!(plan.directory, output.path());
        assert!(plan.used_output_fallback);
    }

    #[test]
    fn invalid_explicit_directory_is_an_error() {
        let output = tempdir().unwrap();
        let missing = output.path().join("missing");
        let error = plan_explicit_staging_directory(
            missing.to_str().unwrap(),
            output.path().join("sample").to_str().unwrap(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("explicit temporary directory"));
    }
}
