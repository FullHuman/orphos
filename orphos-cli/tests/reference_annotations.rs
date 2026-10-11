mod common;

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use insta::assert_snapshot;
use tempfile::NamedTempFile;

use assert_cmd::Command;

use crate::common::parse_gff_cds_records;

/// Genomes that use a translation table other than 11.
const TRANSLATION_TABLES: &[(&str, &str)] = &[("NC_000908.2", "4")];

/// A CDS identified by its 3' end (the stop codon), which is what gene
/// finders get right far more often than the start.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct GeneEnd {
    stop: u32,
    strand: char,
}

#[derive(Debug, Clone, Copy)]
struct Cds {
    start: u32,
    end: u32,
}

impl Cds {
    fn key(&self, strand: char) -> GeneEnd {
        let stop = if strand == '+' { self.end } else { self.start };
        GeneEnd { stop, strand }
    }
}

/// Every complete-genome FASTA fixture (named after its accession) with a
/// RefSeq CDS table next to it.
fn fixtures() -> Vec<(String, String, String)> {
    let mut found = Vec::new();
    for dir in ["tests/data", "tests/data/single"] {
        let mut entries: Vec<_> = fs::read_dir(dir).unwrap().flatten().collect();
        entries.sort_by_key(|entry| entry.path());
        for entry in entries {
            let path = entry.path();
            if path.extension().is_none_or(|ext| ext != "fasta") {
                continue;
            }
            let accession = first_word_of_header(&path);
            let reference = format!("tests/data/refseq/{accession}.cds.tsv");
            let stem = path.file_stem().unwrap().to_string_lossy();
            if stem == accession && Path::new(&reference).exists() {
                found.push((accession, path.to_string_lossy().to_string(), reference));
            }
        }
    }
    found
}

fn first_word_of_header(path: &Path) -> String {
    let content = fs::read_to_string(path).unwrap();
    content
        .lines()
        .next()
        .unwrap()
        .trim_start_matches('>')
        .split_whitespace()
        .next()
        .unwrap()
        .to_string()
}

#[test]
fn refseq_cds_accuracy_summary() {
    unsafe {
        std::env::set_var("RAYON_NUM_THREADS", "1");
    }

    let fixtures = fixtures();
    assert!(!fixtures.is_empty(), "no fixtures with RefSeq CDS tables");

    let mut summary = vec![
        "fixture\trefseq_cds\torphos_cds\tfound(%)\texact(%)\tstart_ok_of_found(%)\tfalse_positive\tmissed\twrong_start"
            .to_string(),
    ];
    let mut totals = [0usize; 5]; // refseq, found, exact, false positive, predicted

    for (accession, input_path, reference_path) in &fixtures {
        let reference = load_refseq_cds(reference_path);
        let rust_out = NamedTempFile::new().unwrap();
        run_orphos_single(accession, input_path, rust_out.path().to_str().unwrap());
        let predicted = cds_by_end(&fs::read_to_string(rust_out.path()).unwrap());

        let mut found = 0;
        let mut exact = 0;
        for (key, cds) in &reference {
            if let Some(prediction) = predicted.get(key) {
                found += 1;
                if prediction.start == cds.start && prediction.end == cds.end {
                    exact += 1;
                }
            }
        }
        let false_positive = predicted
            .keys()
            .filter(|key| !reference.contains_key(key))
            .count();

        summary.push(format!(
            "{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\t{}",
            accession,
            reference.len(),
            predicted.len(),
            percentage(found, reference.len()),
            percentage(exact, reference.len()),
            percentage(exact, found),
            false_positive,
            reference.len() - found,
            found - exact,
        ));
        totals[0] += reference.len();
        totals[1] += found;
        totals[2] += exact;
        totals[3] += false_positive;
        totals[4] += predicted.len();
    }
    summary.push(format!(
        "TOTAL\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\t{}",
        totals[0],
        totals[4],
        percentage(totals[1], totals[0]),
        percentage(totals[2], totals[0]),
        percentage(totals[2], totals[1]),
        totals[3],
        totals[0] - totals[1],
        totals[1] - totals[2],
    ));

    assert_snapshot!("refseq_cds_accuracy_summary", summary.join("\n"));
}

fn run_orphos_single(accession: &str, input_path: &str, output_path: &str) {
    let mut cmd = Command::cargo_bin("orphos").unwrap();
    cmd.args(["-i", input_path, "-o", output_path, "-f", "gff", "-q"]);
    if let Some((_, table)) = TRANSLATION_TABLES.iter().find(|(acc, _)| acc == &accession) {
        cmd.args(["-g", table]);
    }
    cmd.assert().success();
}

fn cds_by_end(gff: &str) -> HashMap<GeneEnd, Cds> {
    parse_gff_cds_records(gff)
        .into_iter()
        .map(|record| {
            let cds = Cds {
                start: record.start,
                end: record.end,
            };
            (cds.key(record.strand), cds)
        })
        .collect()
}

fn load_refseq_cds(path: &str) -> HashMap<GeneEnd, Cds> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .filter_map(|line| {
            let mut cols = line.split('\t');
            let start = cols.next()?.parse().ok()?;
            let end = cols.next()?.parse().ok()?;
            let strand = cols.next()?.chars().next()?;
            let cds = Cds { start, end };
            Some((cds.key(strand), cds))
        })
        .collect()
}

fn percentage(part: usize, whole: usize) -> f64 {
    if whole == 0 {
        0.0
    } else {
        part as f64 * 100.0 / whole as f64
    }
}
