mod common;

use assert_cmd::Command;
use std::fs;
use tempfile::NamedTempFile;

use crate::common::{normalize_output, prodigal_available, prodigal_binary};

const MYCOPLASMA: &str = "tests/data/single/NC_000908.2.fasta";
const FRAGMENT: &str = "tests/data/single/fragment2.fasta";

fn run_orphos_with(input: &str, output: &str, args: &[&str]) -> String {
    let mut cmd = Command::cargo_bin("orphos").unwrap();
    cmd.arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output)
        .arg("-f")
        .arg("gff");
    cmd.args(args);
    let assert = cmd.assert().success();
    String::from_utf8_lossy(&assert.get_output().stderr).to_string()
}

fn run_prodigal_with(input: &str, output: &str, args: &[&str]) {
    let status = std::process::Command::new(prodigal_binary())
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output)
        .arg("-f")
        .arg("gff")
        .arg("-q")
        .args(args)
        .status()
        .unwrap();
    assert!(status.success());
}

fn model_line(gff: &str) -> String {
    gff.lines()
        .find(|line| line.starts_with("# Model Data"))
        .unwrap()
        .to_string()
}

#[test]
fn translation_table_option_is_used_for_training() {
    if !std::path::Path::new(MYCOPLASMA).exists() {
        eprintln!("Skipping: fixture missing");
        return;
    }
    let default_out = NamedTempFile::new().unwrap();
    let table4_out = NamedTempFile::new().unwrap();
    run_orphos_with(MYCOPLASMA, default_out.path().to_str().unwrap(), &[]);
    run_orphos_with(
        MYCOPLASMA,
        table4_out.path().to_str().unwrap(),
        &["-g", "4"],
    );

    let default_gff = fs::read_to_string(default_out.path()).unwrap();
    let table4_gff = fs::read_to_string(table4_out.path()).unwrap();
    assert!(model_line(&default_gff).contains("transl_table=11"));
    assert!(model_line(&table4_gff).contains("transl_table=4"));
    // Mycoplasma genitalium reads TGA as tryptophan; table 4 must change the calls.
    assert_ne!(default_gff, table4_gff);

    if prodigal_available() {
        let prodigal_out = NamedTempFile::new().unwrap();
        run_prodigal_with(
            MYCOPLASMA,
            prodigal_out.path().to_str().unwrap(),
            &["-g", "4"],
        );
        let prodigal_gff = fs::read_to_string(prodigal_out.path()).unwrap();
        assert_eq!(
            normalize_output(&prodigal_gff),
            normalize_output(&table4_gff)
        );
    }
}

#[test]
fn boolean_flags_take_no_value() {
    let closed_out = NamedTempFile::new().unwrap();
    let stderr = run_orphos_with(FRAGMENT, closed_out.path().to_str().unwrap(), &["-c", "-q"]);
    assert!(
        stderr.is_empty(),
        "quiet mode should print nothing, got: {stderr}"
    );

    let open_out = NamedTempFile::new().unwrap();
    run_orphos_with(FRAGMENT, open_out.path().to_str().unwrap(), &[]);
    let closed_gff = fs::read_to_string(closed_out.path()).unwrap();
    let open_gff = fs::read_to_string(open_out.path()).unwrap();
    // The fragment starts and ends inside genes, so closed ends drop edge genes.
    assert_ne!(closed_gff, open_gff);

    if prodigal_available() {
        let prodigal_out = NamedTempFile::new().unwrap();
        run_prodigal_with(FRAGMENT, prodigal_out.path().to_str().unwrap(), &["-c"]);
        let prodigal_gff = fs::read_to_string(prodigal_out.path()).unwrap();
        assert_eq!(
            normalize_output(&prodigal_gff),
            normalize_output(&closed_gff)
        );
    }

    Command::cargo_bin("orphos")
        .unwrap()
        .args(["-i", FRAGMENT, "-o", "/dev/null", "-c", "-r"])
        .assert()
        .failure();
}
