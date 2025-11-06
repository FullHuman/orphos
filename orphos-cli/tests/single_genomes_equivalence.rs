mod common;

use crate::common::{normalize_output, prodigal_available, run_original_prodigal, run_orphos};
use insta::assert_snapshot;
use sha2::{Digest, Sha256};
use std::fs;
use std::path::PathBuf;
use tempfile::NamedTempFile;

const SIMILARITY_THRESHOLD: f32 = 0.9993; // 99.93%

use std::collections::HashMap;

fn similarity(a: &str, b: &str) -> f32 {
    let a_n = normalize_output(a);
    let b_n = normalize_output(b);
    similar::TextDiff::from_lines(&a_n, &b_n).ratio() as f32
}

fn sha256_normalized(s: &str) -> String {
    let mut hasher = Sha256::new();
    hasher.update(normalize_output(s));
    format!("{:x}", hasher.finalize())
}

fn collect_fixtures(dir: &str) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map(|e| e == "fasta").unwrap_or(false))
        .collect();
    files.sort();
    files
}

// Parse GFF content and collect a multiset of gene coordinates from CDS entries
// Coordinates are compared as (start, end, strand_char). We use a multiset to be robust to duplicates.
fn gff_coord_multiset(s: &str) -> HashMap<(u32, u32, char), usize> {
    let mut map: HashMap<(u32, u32, char), usize> = HashMap::new();
    for line in s.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 7 {
            continue;
        }
        if cols[2] != "CDS" {
            continue;
        }
        // GFF columns: seqid, source, type, start, end, score, strand, phase, attributes
        let start: u32 = match cols[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: u32 = match cols[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let strand_char = cols[6].chars().next().unwrap_or('.');
        *map.entry((start, end, strand_char)).or_insert(0) += 1;
    }
    map
}

// Compute percentage of original genes whose coordinates exactly match a gene in the Rust output.
fn coord_match_percent(orig_gff: &str, rust_gff: &str) -> f32 {
    let orig = gff_coord_multiset(orig_gff);
    let rust = gff_coord_multiset(rust_gff);
    let orig_total: usize = orig.values().sum();
    if orig_total == 0 {
        // If both are empty, count as perfect match; otherwise 0
        let rust_total: usize = rust.values().sum();
        return if rust_total == 0 { 100.0 } else { 0.0 };
    }
    let mut matches: usize = 0;
    for (k, &cnt) in &orig {
        let r = rust.get(k).copied().unwrap_or(0);
        matches += cnt.min(r);
    }
    (matches as f32) * 100.0 / (orig_total as f32)
}

#[test]
fn single_genome_fixtures_equivalence_summary() {
    if !prodigal_available() {
        eprintln!("Skipping: original prodigal not available in PATH");
        return;
    }

    let fixtures = collect_fixtures("tests/data/single");
    assert!(
        !fixtures.is_empty(),
        "No fixtures found in tests/data/single"
    );

    let mut summary_lines = Vec::new();
    summary_lines
        .push("file\tsimilarity(%)\tcoord_match(%)\texact\torig_hash\trust_hash".to_string());

    let mut below_threshold = Vec::new();

    for fixture in fixtures {
        let file_name = fixture.file_name().unwrap().to_string_lossy().to_string();
        println!("Processing fixture: {}", file_name);

        let orig_out = NamedTempFile::new().unwrap();
        run_original_prodigal(
            fixture.to_str().unwrap(),
            orig_out.path().to_str().unwrap(),
            "gff",
            "single",
        )
        .unwrap();

        let rust_out = NamedTempFile::new().unwrap();
        run_orphos(
            fixture.to_str().unwrap(),
            rust_out.path().to_str().unwrap(),
            "gff",
            "single",
        )
        .unwrap();

        let orig_content = fs::read_to_string(orig_out.path()).unwrap();
        let rust_content = fs::read_to_string(rust_out.path()).unwrap();

        let sim = similarity(&orig_content, &rust_content);
        let coord_pct = coord_match_percent(&orig_content, &rust_content);
        let exact = normalize_output(&orig_content) == normalize_output(&rust_content);

        let orig_hash = sha256_normalized(&orig_content)[..12].to_string();
        let rust_hash = sha256_normalized(&rust_content)[..12].to_string();

        if coord_pct < SIMILARITY_THRESHOLD {
            below_threshold.push((file_name.clone(), sim));
        }

        summary_lines.push(format!(
            "{}\t{:.3}\t{:.3}\t{}\t{}\t{}",
            file_name,
            sim * 100.0,
            coord_pct,
            if exact { "yes" } else { "no" },
            orig_hash,
            rust_hash
        ));
    }

    let summary = summary_lines.join("\n");
    assert_snapshot!("single_genomes_gff_summary", summary);

    if !below_threshold.is_empty() {
        let details = below_threshold
            .into_iter()
            .map(|(f, s)| format!("{}: {:.2}%", f, s * 100.0))
            .collect::<Vec<_>>()
            .join(", ");
        panic!(
            "Coordinate percentage below threshold ({:.2}%) for: {}",
            SIMILARITY_THRESHOLD * 100.0,
            details
        );
    }
}
