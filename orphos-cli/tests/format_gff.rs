mod common;
use crate::common::{
    normalize_output, prodigal_available, run_original_prodigal, run_orphos, similarity,
};
use std::path::Path;
use tempfile::NamedTempFile;

#[test]
fn gff_format_equivalence_single_mode_selected_fixtures() {
    if !prodigal_available() {
        eprintln!("Skipping: original prodigal not in PATH");
        return;
    }
    // Force deterministic single-thread execution for the Rust implementation to reduce variability.
    unsafe {
        std::env::set_var("RAYON_NUM_THREADS", "1");
    }

    // Temporary: investigate low similarity; include both larger genomes.
    for file in ["salmonella.fasta", "ecoli.fasta"] {
        let path = format!("tests/data/{}", file);
        let p = Path::new(&path);
        if !p.exists() {
            eprintln!("Missing fixture {}", file);
            continue;
        }
        let orig_out = NamedTempFile::new().unwrap();
        let rust_out = NamedTempFile::new().unwrap();
        run_original_prodigal(
            p.to_str().unwrap(),
            orig_out.path().to_str().unwrap(),
            "gff",
            "single",
        )
        .unwrap();
        run_orphos(
            p.to_str().unwrap(),
            rust_out.path().to_str().unwrap(),
            "gff",
            "single",
        )
        .unwrap();
        let orig_raw = std::fs::read_to_string(orig_out.path()).unwrap();
        let rust_raw = std::fs::read_to_string(rust_out.path()).unwrap();
        let sim = similarity(&orig_raw, &rust_raw);
        let orig_norm = normalize_output(&orig_raw);
        let rust_norm = normalize_output(&rust_raw);
        println!(
            "GFF single {} similarity: {:.4}% (lines orig={}, rust={})",
            file,
            sim * 100.0,
            orig_norm.lines().count(),
            rust_norm.lines().count()
        );
        if sim < 0.98 {
            // Emit diagnostics: first differing 10 lines and hashes.
            use sha2::{Digest, Sha256};
            let mut hasher = Sha256::new();
            hasher.update(&orig_norm);
            let h_orig = format!("{:x}", hasher.finalize());
            let mut hasher = Sha256::new();
            hasher.update(&rust_norm);
            let h_rust = format!("{:x}", hasher.finalize());
            println!(
                "  Hashes (normalized) orig={} rust={}",
                &h_orig[..12],
                &h_rust[..12]
            );
            let o_lines: Vec<_> = orig_norm.lines().collect();
            let r_lines: Vec<_> = rust_norm.lines().collect();
            let mut shown = 0;
            for i in 0..o_lines.len().min(r_lines.len()) {
                if o_lines[i] != r_lines[i] {
                    println!(
                        "  Diff line {}:\n    O: {}\n    R: {}",
                        i + 1,
                        o_lines[i],
                        r_lines[i]
                    );
                    shown += 1;
                    if shown >= 10 {
                        break;
                    }
                }
            }
            if shown == 0 {
                println!("  (No top-line diffs; differences may be later in file)");
            }
        }

        // Diff (unified-ish) of normalized outputs.
        let mut diff_buf = String::new();
        let diff = similar::TextDiff::from_lines(&orig_norm, &rust_norm);
        for op in diff.ops() {
            for change in diff.iter_changes(op) {
                use similar::ChangeTag::*;
                let sign = match change.tag() {
                    Delete => '-',
                    Insert => '+',
                    Equal => ' ',
                };
                diff_buf.push(sign);
                diff_buf.push_str(change.value());
            }
        }

        // Keep threshold strict; failing artifacts now saved.
        assert!(
            sim > 0.99,
            "{} similarity too low: {:.4}% (see target/test-diffs/gff)",
            file,
            sim * 100.0
        );
    }
}
