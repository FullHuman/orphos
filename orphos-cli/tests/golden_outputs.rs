mod common;

use assert_cmd::Command;
use insta::assert_snapshot;
use std::fs;
use tempfile::NamedTempFile;

use crate::common::run_orphos;

// Simple golden snapshot test for --help output (stable CLI surface)
#[test]
fn cli_help_snapshot() {
    let mut cmd = Command::cargo_bin("orphos").unwrap();
    cmd.arg("--help");
    let output = cmd.assert().success().get_output().stdout.clone();
    let text = String::from_utf8(output).unwrap();

    // Normalize binary name for cross-platform compatibility (remove .exe on Windows)
    let normalized = text.replace("orphos-cli.exe", "orphos-cli");

    assert_snapshot!("cli_help", normalized);
}

// Golden snapshot for GFF output (captures only first N lines for stability)
#[test]
fn small_gff_output_snapshot() {
    let input_path = "tests/data/ecoli.fasta"; // existing fixture
    if !std::path::Path::new(input_path).exists() {
        eprintln!("Skipping: fixture missing");
        return;
    }
    let out_tmp = NamedTempFile::new().unwrap();
    run_orphos(
        input_path,
        out_tmp.path().to_str().unwrap(),
        "gff",
        "single",
    )
    .unwrap();
    let raw = fs::read_to_string(out_tmp.path()).unwrap();
    // Take only first 25 lines to keep snapshot concise
    let head: String = raw.lines().take(25).collect::<Vec<_>>().join("\n");
    assert_snapshot!("ecoli_gff_head", head);
}
