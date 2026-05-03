use std::io::Write;

use bio::io::bed;

use crate::{OrphosError, results::OrphosResults};

fn confidence_to_bed_score(confidence: f64) -> i32 {
    if !confidence.is_finite() {
        return 0;
    }

    (confidence * 10.0).round().clamp(0.0, 1000.0) as i32
}

/// Write results in BED6 format.
///
/// Columns: chrom, chromStart, chromEnd, name, score, strand
pub fn write_bed_format<W: Write>(
    writer: &mut W,
    results: &OrphosResults,
) -> Result<(), OrphosError> {
    let mut bed_writer = bed::Writer::new(&mut *writer as &mut dyn Write);
    for (i, gene) in results.genes.iter().enumerate() {
        if gene.coordinates.begin > gene.coordinates.end {
            // BED6 has no native circular wrap representation.
            continue;
        }

        let strand = gene.coordinates.strand;

        let name = if gene.annotation.identifier.is_empty() {
            format!("{}_{}", results.sequence_info.header, i + 1)
        } else {
            gene.annotation.identifier.clone()
        };

        let mut record = bed::Record::new();
        record.set_chrom(&results.sequence_info.header);
        record.set_start(gene.coordinates.begin as u64);
        record.set_end((gene.coordinates.end + 1) as u64);
        record.set_name(&name);
        record.set_score(&confidence_to_bed_score(gene.score.confidence).to_string());
        record.set_strand(strand);
        bed_writer.write(&record).map_err(std::io::Error::other)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        results::{OrphosResults, SequenceInfo},
        types::{Gene, GeneAnnotation, GeneCoordinates, GeneScore, StartType, Training},
    };
    use bio::bio_types::strand::Strand;
    use std::io::Cursor;

    fn create_results(genes: Vec<Gene>) -> OrphosResults {
        OrphosResults {
            sequence_info: SequenceInfo {
                header: "test_seq".to_string(),
                length: 1000,
                description: Some("Test sequence".to_string()),
                gc_content: 0.5,
                num_genes: genes.len(),
            },
            genes,
            training_used: Training::default(),
            metagenomic_model: None,
        }
    }

    #[test]
    fn test_write_bed_format_single_gene() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 99,
                end: 299,
                strand: Strand::Forward,
                ..Default::default()
            },
            score: GeneScore {
                confidence: 95.5,
                ..Default::default()
            },
            annotation: GeneAnnotation::new(
                "gene_1".to_string(),
                false,
                false,
                StartType::Atg,
                0.45,
            ),
            ..Default::default()
        };

        let results = create_results(vec![gene]);
        write_bed_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "test_seq\t99\t300\tgene_1\t955\t+\n");
    }

    #[test]
    fn test_write_bed_format_fallback_name() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 10,
                end: 50,
                strand: Strand::Reverse,
                ..Default::default()
            },
            score: GeneScore {
                confidence: 50.0,
                ..Default::default()
            },
            ..Default::default()
        };

        let results = create_results(vec![gene]);
        write_bed_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "test_seq\t10\t51\ttest_seq_1\t500\t-\n");
    }

    #[test]
    fn test_write_bed_format_unknown_strand() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 0,
                end: 0,
                strand: Strand::Unknown,
                ..Default::default()
            },
            score: GeneScore {
                confidence: 0.0,
                ..Default::default()
            },
            ..Default::default()
        };

        let results = create_results(vec![gene]);
        write_bed_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "test_seq\t0\t1\ttest_seq_1\t0\t.\n");
    }

    #[test]
    fn test_write_bed_format_score_rounding_and_clamp() {
        assert_eq!(confidence_to_bed_score(99.94), 999);
        assert_eq!(confidence_to_bed_score(99.95), 1000);
        assert_eq!(confidence_to_bed_score(120.0), 1000);
        assert_eq!(confidence_to_bed_score(-1.0), 0);
        assert_eq!(confidence_to_bed_score(f64::NAN), 0);
    }

    #[test]
    fn test_write_bed_format_no_genes() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let results = create_results(vec![]);
        write_bed_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "");
    }

    #[test]
    fn test_write_bed_format_skips_wrapped_gene() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 900,
                end: 100,
                strand: Strand::Forward,
                ..Default::default()
            },
            ..Default::default()
        };

        let results = create_results(vec![gene]);
        write_bed_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "");
    }
}
