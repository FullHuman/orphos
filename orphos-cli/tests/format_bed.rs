mod common;

use crate::common::run_orphos;
use std::path::Path;
use tempfile::NamedTempFile;

#[test]
fn bed_format_structure_selected_fixtures() {
    for file in ["salmonella.fasta", "ecoli.fasta"] {
        let path = format!("tests/data/{}", file);
        let p = Path::new(&path);
        if !p.exists() {
            eprintln!("Missing fixture {}", file);
            continue;
        }

        let bed_out = NamedTempFile::new().unwrap();
        run_orphos(
            p.to_str().unwrap(),
            bed_out.path().to_str().unwrap(),
            "bed",
            "single",
        )
        .unwrap();

        let content = std::fs::read_to_string(bed_out.path()).unwrap();
        for (line_no, line) in content.lines().enumerate() {
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                cols.len(),
                6,
                "{} line {} should have 6 BED columns",
                file,
                line_no + 1
            );

            assert!(
                !cols[0].is_empty(),
                "{} line {} expected non-empty chrom",
                file,
                line_no + 1
            );

            let start: usize = cols[1]
                .parse()
                .unwrap_or_else(|_| panic!("{} line {} invalid start", file, line_no + 1));
            let end: usize = cols[2]
                .parse()
                .unwrap_or_else(|_| panic!("{} line {} invalid end", file, line_no + 1));
            assert!(
                start < end,
                "{} line {} expected start < end, got {} >= {}",
                file,
                line_no + 1,
                start,
                end
            );

            let score: i32 = cols[4]
                .parse()
                .unwrap_or_else(|_| panic!("{} line {} invalid score", file, line_no + 1));
            assert!(
                (0..=1000).contains(&score),
                "{} line {} score out of BED range: {}",
                file,
                line_no + 1,
                score
            );

            assert!(
                matches!(cols[5], "+" | "-" | "."),
                "{} line {} invalid strand symbol: {}",
                file,
                line_no + 1,
                cols[5]
            );
        }

        assert!(
            !content.trim().is_empty(),
            "{} produced empty BED output",
            file
        );
    }
}
