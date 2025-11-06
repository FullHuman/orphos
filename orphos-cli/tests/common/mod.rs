#![allow(dead_code)]

use std::process;

use assert_cmd::Command;
use std::io::ErrorKind;

/// Runs the original prodigal with given arguments
pub fn run_original_prodigal(
    input_file: &str,
    output_file: &str,
    format: &str,
    mode: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = process::Command::new("prodigal");
    cmd.arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-f")
        .arg(format)
        .arg("-p")
        .arg(mode);

    let output = cmd.output()?;
    if !output.status.success() {
        return Err(format!(
            "Original prodigal failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }
    Ok(())
}

/// Runs the Orphos CLI with given arguments
pub fn run_orphos(
    input_file: &str,
    output_file: &str,
    format: &str,
    mode: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("orphos")?;
    cmd.arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-f")
        .arg(format)
        .arg("-p")
        .arg(mode);

    cmd.assert().success();
    Ok(())
}

/// Normalize dynamic fields in outputs to reduce snapshot churn.
/// Currently strips Rust-specific version suffix and collapses whitespace in headers.
pub fn normalize_output(s: &str) -> String {
    s.replace("Orphos.v1.0.0", "Prodigal.v2.6.3")
        .replace("Orphos_v1.0.0", "Prodigal_v2.6.3")
}

/// Compute a quick line-based similarity (ratio 0..1) for two strings after normalization.
pub fn similarity(a: &str, b: &str) -> f32 {
    let a_n = normalize_output(a);
    let b_n = normalize_output(b);
    similar::TextDiff::from_lines(&a_n, &b_n).ratio() as f32
}

/// Detect if the original prodigal binary is available on PATH.
/// More robust than just `--help` because the C prodigal may differ in flags/exit codes.
pub fn prodigal_available() -> bool {
    // Attempt a sequence of lightweight invocations; treat any successful spawn as availability.
    const CANDIDATE_ARGS: &[&[&str]] = &[&["-v"], &["-h"], &["--help"], &[]];
    for args in CANDIDATE_ARGS {
        match std::process::Command::new("prodigal").args(*args).output() {
            Ok(output) => {
                if output.status.success()
                    || String::from_utf8_lossy(&output.stdout).contains("Prodigal")
                    || String::from_utf8_lossy(&output.stderr).contains("Prodigal")
                {
                    return true;
                }
            }
            Err(e) => {
                if e.kind() == ErrorKind::NotFound {
                    return false; // Binary not in PATH
                }
                // Other IO errors: continue trying other arg variants
            }
        }
    }
    false
}
