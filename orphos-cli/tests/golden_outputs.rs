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
    let normalized = text.replace("orphos.exe", "orphos");

    assert_snapshot!("cli_help", normalized);
}

// Golden snapshot for GFF output (captures only first N lines for stability)
#[test]
fn small_gff_output_snapshot() {
    assert_format_head_snapshot("gff", "ecoli_gff_head", 25);
}

#[test]
fn small_gbk_output_snapshot() {
    assert_format_head_snapshot("gbk", "ecoli_gbk_head", 25);
}

#[test]
fn small_sco_output_snapshot() {
    assert_format_head_snapshot("sco", "ecoli_sco_head", 25);
}

#[test]
fn small_gca_output_snapshot() {
    assert_format_head_snapshot("gca", "ecoli_gca_head", 25);
}

#[test]
fn small_bed_output_snapshot() {
    assert_format_head_snapshot("bed", "ecoli_bed_head", 25);
}

fn assert_format_head_snapshot(format: &str, snapshot_name: &str, max_lines: usize) {
    let input_path = "tests/data/ecoli.fasta";
    if !std::path::Path::new(input_path).exists() {
        eprintln!("Skipping: fixture missing");
        return;
    }

    let out_tmp = NamedTempFile::new().unwrap();
    run_orphos(
        input_path,
        out_tmp.path().to_str().unwrap(),
        format,
        "single",
    )
    .unwrap();
    let raw = fs::read_to_string(out_tmp.path()).unwrap();
    let head: String = raw.lines().take(max_lines).collect::<Vec<_>>().join("\n");
    assert_snapshot!(snapshot_name, head);
}
