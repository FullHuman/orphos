//! Orphos against GFF files written by Prodigal, line by line.
//!
//! The genome fixtures of the other equivalence tests are whole genomes run
//! with default options. These cases cover what those leave out: inputs of
//! many sequences, sequences too short to hold a gene, runs of `N`, closed
//! ends, masking and another translation table. See
//! `tests/data/prodigal/README.md` for the inputs and the Prodigal build.

mod common;

use std::fs;
use std::path::Path;

use assert_cmd::Command;
use tempfile::NamedTempFile;

use crate::common::normalize_output;

struct Case {
    /// Name of the expected file, without its extension.
    name: &'static str,
    input: &'static str,
    options: &'static [&'static str],
}

const CASES: &[Case] = &[
    Case {
        name: "contigs.meta",
        input: "tests/data/prodigal/contigs.fasta",
        options: &["-p", "meta"],
    },
    Case {
        name: "contigs.meta.closed",
        input: "tests/data/prodigal/contigs.fasta",
        options: &["-p", "meta", "-c"],
    },
    Case {
        name: "gaps.meta",
        input: "tests/data/prodigal/gaps.fasta",
        options: &["-p", "meta"],
    },
    Case {
        name: "gaps.meta.masked",
        input: "tests/data/prodigal/gaps.fasta",
        options: &["-p", "meta", "-m"],
    },
    Case {
        name: "gaps.meta.closed.masked",
        input: "tests/data/prodigal/gaps.fasta",
        options: &["-p", "meta", "-c", "-m"],
    },
    Case {
        name: "draft.single",
        input: "tests/data/prodigal/draft.fasta",
        options: &[],
    },
    Case {
        name: "draft.single.masked",
        input: "tests/data/prodigal/draft.fasta",
        options: &["-m"],
    },
    Case {
        name: "draft.single.closed",
        input: "tests/data/prodigal/draft.fasta",
        options: &["-c"],
    },
    Case {
        name: "draft.single.table4",
        input: "tests/data/prodigal/draft.fasta",
        options: &["-g", "4"],
    },
    Case {
        name: "buchnera.single.closed",
        input: "tests/data/single/NC_008512.1.fasta",
        options: &["-c"],
    },
    Case {
        name: "fragment2.meta.closed",
        input: "tests/data/single/fragment2.fasta",
        options: &["-p", "meta", "-c"],
    },
];

fn expected_output(case: &Case) -> String {
    let path = Path::new("tests/data/prodigal/expected").join(format!("{}.gff", case.name));
    fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The GFF output of Orphos for `case` on `threads` threads (all of them
/// for `None`), with the tool name Prodigal writes.
fn orphos_output(case: &Case, threads: Option<&str>) -> String {
    let output = NamedTempFile::new().unwrap();
    let mut command = Command::cargo_bin("orphos").unwrap();
    command
        .args(["-q", "-f", "gff", "-i", case.input, "-o"])
        .arg(output.path())
        .args(case.options);
    match threads {
        Some(threads) => command.env("RAYON_NUM_THREADS", threads),
        None => command.env_remove("RAYON_NUM_THREADS"),
    };
    command.assert().success();
    normalize_output(&fs::read_to_string(output.path()).unwrap())
}

/// Panics, naming the first line that differs, unless `actual` is `expected`.
fn assert_same_lines(expected: &str, actual: &str, what: &str) {
    let expected_lines: Vec<&str> = expected.lines().collect();
    let actual_lines: Vec<&str> = actual.lines().collect();
    let differing = expected_lines
        .iter()
        .zip(&actual_lines)
        .position(|(expected, actual)| expected != actual);
    if let Some(line) = differing {
        panic!(
            "{what}: line {} differs\n  Prodigal: {}\n  Orphos:   {}",
            line + 1,
            expected_lines[line],
            actual_lines[line]
        );
    }
    assert_eq!(
        expected_lines.len(),
        actual_lines.len(),
        "{what}: the outputs agree up to the end of the shorter one"
    );
}

#[test]
fn output_matches_prodigal_line_by_line() {
    for case in CASES {
        let expected = expected_output(case);
        assert!(
            expected.lines().any(|line| !line.starts_with('#')),
            "{}: Prodigal predicted no gene",
            case.name
        );
        for threads in [Some("1"), None] {
            assert_same_lines(
                &expected,
                &orphos_output(case, threads),
                &format!("{} on {} thread(s)", case.name, threads.unwrap_or("all")),
            );
        }
    }
}

/// The sequences of an input are numbered as Prodigal numbers them, in the
/// sequence headers and in the gene identifiers, whether or not a sequence
/// has genes.
#[test]
fn sequences_are_numbered_like_prodigal() {
    let case = &CASES[1];
    let expected = expected_output(case);
    let actual = orphos_output(case, None);

    let sequence_headers = |gff: &str| -> Vec<String> {
        gff.lines()
            .filter(|line| line.starts_with("# Sequence Data"))
            .map(String::from)
            .collect()
    };
    let headers = sequence_headers(&actual);
    assert_eq!(headers.len(), 90);
    assert_eq!(headers, sequence_headers(&expected));
    assert!(headers[89].starts_with("# Sequence Data: seqnum=90;"));

    // The first gene of every sequence that has one is `<sequence>_1`.
    let mut sequence = 0;
    let mut sequences_with_genes = 0;
    let mut first_gene_expected = false;
    for line in actual.lines() {
        if line.starts_with("# Sequence Data") {
            sequence += 1;
            first_gene_expected = true;
        } else if !line.starts_with('#') && first_gene_expected {
            assert!(
                line.contains(&format!("\tID={sequence}_1;")),
                "sequence {sequence}: {line}"
            );
            first_gene_expected = false;
            sequences_with_genes += 1;
        }
    }
    assert!(sequences_with_genes > 40 && sequences_with_genes < 90);
}

/// The expected files are what the Prodigal binary named by
/// `ORPHOS_PRODIGAL_BIN` writes. Without the variable there is nothing to
/// check them against: the 2.6.3 release differs from them (see the README
/// of the fixtures).
#[test]
fn expected_files_are_what_prodigal_writes() {
    let Ok(prodigal) = std::env::var("ORPHOS_PRODIGAL_BIN") else {
        eprintln!("Skipping: ORPHOS_PRODIGAL_BIN is not set");
        return;
    };

    for case in CASES {
        let output = NamedTempFile::new().unwrap();
        let status = std::process::Command::new(&prodigal)
            .args(["-q", "-f", "gff", "-i", case.input, "-o"])
            .arg(output.path())
            .args(case.options)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .unwrap();
        assert!(status.success(), "{}: {prodigal} failed", case.name);
        assert_same_lines(
            &fs::read_to_string(output.path()).unwrap(),
            &expected_output(case),
            &format!(
                "{} (first: {prodigal}, second: the expected file)",
                case.name
            ),
        );
    }
}
