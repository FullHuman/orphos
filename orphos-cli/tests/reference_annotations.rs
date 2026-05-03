mod common;

use std::collections::HashSet;
use std::fs;

use insta::assert_snapshot;
use tempfile::NamedTempFile;

use crate::common::{parse_gff_cds_records, run_orphos};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct CdsCoord {
    start: u32,
    end: u32,
    strand: char,
}

#[derive(Debug, Clone, Copy)]
struct ReferenceFixture {
    label: &'static str,
    accession: &'static str,
    input_path: &'static str,
    reference_path: &'static str,
}

const SELECTED_FIXTURES: &[ReferenceFixture] = &[
    ReferenceFixture {
        label: "salmonella",
        accession: "NC_003198.1",
        input_path: "tests/data/salmonella.fasta",
        reference_path: "tests/data/refseq/NC_003198.1.cds.tsv",
    },
    ReferenceFixture {
        label: "ecoli",
        accession: "NC_000913.3",
        input_path: "tests/data/ecoli.fasta",
        reference_path: "tests/data/refseq/NC_000913.3.cds.tsv",
    },
    ReferenceFixture {
        label: "helicobacter",
        accession: "NC_000915.1",
        input_path: "tests/data/single/NC_000915.1.fasta",
        reference_path: "tests/data/refseq/NC_000915.1.cds.tsv",
    },
    ReferenceFixture {
        label: "bacillus",
        accession: "NC_000964.3",
        input_path: "tests/data/single/NC_000964.3.fasta",
        reference_path: "tests/data/refseq/NC_000964.3.cds.tsv",
    },
];

#[test]
fn selected_fixtures_full_refseq_cds_summary() {
    unsafe {
        std::env::set_var("RAYON_NUM_THREADS", "1");
    }

    let mut summary = vec![
        "fixture\taccession\trefseq_cds\torphos_cds\torphos_hits(%)\torphos_only_vs_ref\tref_missing_orphos".to_string(),
    ];

    for fixture in SELECTED_FIXTURES {
        let refseq_coords = load_refseq_cds(fixture.reference_path);
        assert!(
            !refseq_coords.is_empty(),
            "No RefSeq CDS coordinates loaded for {}",
            fixture.accession
        );

        let rust_out = NamedTempFile::new().unwrap();

        run_orphos(
            fixture.input_path,
            rust_out.path().to_str().unwrap(),
            "gff",
            "single",
        )
        .unwrap();

        let orphos_coords = cds_coord_set(&fs::read_to_string(rust_out.path()).unwrap());
        let orphos_hits = intersection_size(&refseq_coords, &orphos_coords);

        summary.push(format!(
            "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
            fixture.label,
            fixture.accession,
            refseq_coords.len(),
            orphos_coords.len(),
            percentage(orphos_hits, refseq_coords.len()),
            orphos_coords.difference(&refseq_coords).count(),
            refseq_coords.difference(&orphos_coords).count(),
        ));
    }

    assert_snapshot!("selected_refseq_full_cds_summary", summary.join("\n"));
}

fn cds_coord_set(gff: &str) -> HashSet<CdsCoord> {
    parse_gff_cds_records(gff)
        .into_iter()
        .map(|record| CdsCoord {
            start: record.start,
            end: record.end,
            strand: record.strand,
        })
        .collect()
}

fn load_refseq_cds(path: &str) -> HashSet<CdsCoord> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .filter_map(|line| {
            let mut cols = line.split('\t');
            let start = cols.next()?.parse().ok()?;
            let end = cols.next()?.parse().ok()?;
            let strand = cols.next()?.chars().next()?;
            Some(CdsCoord { start, end, strand })
        })
        .collect()
}

fn intersection_size(left: &HashSet<CdsCoord>, right: &HashSet<CdsCoord>) -> usize {
    left.intersection(right).count()
}

fn percentage(matches: usize, total: usize) -> f32 {
    (matches as f32) * 100.0 / (total as f32)
}
