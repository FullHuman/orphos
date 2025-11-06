use std::io::Write;

use crate::{results::OrphosResults, types::OrphosError};

/// Write results in GenBank format
pub fn write_gbk_format<W: Write>(
    writer: &mut W,
    results: &OrphosResults,
) -> Result<(), OrphosError> {
    for (i, gene) in results.genes.iter().enumerate() {
        writeln!(
            writer,
            "     CDS             {}..{}",
            gene.coordinates.begin + 1, // GenBank uses 1-based coordinates
            gene.coordinates.end + 1
        )?;
        writeln!(
            writer,
            "                     /locus_tag=\"{}_{}\"",
            results.sequence_info.header,
            i + 1
        )?;
        writeln!(
            writer,
            "                     /confidence={:.2}",
            gene.score.confidence
        )?;
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use bio::bio_types::strand::Strand;

    use crate::{
        results::{OrphosResults, SequenceInfo},
        types::{Gene, GeneCoordinates, GeneScore, Training},
    };

    use super::*;

    #[test]
    fn test_write_gbk_format_empty_genes() {
        let mut buffer = Vec::new();
        let results = OrphosResults {
            genes: vec![],
            sequence_info: SequenceInfo {
                header: "test_seq".to_string(),
                length: 1000,
                description: Some("Test sequence".to_string()),
                gc_content: 12.0,
                num_genes: 25,
            },
            training_used: Training::default(),
            metagenomic_model: None,
        };

        let result = write_gbk_format(&mut buffer, &results);
        assert!(result.is_ok());
        assert_eq!(String::from_utf8(buffer).unwrap(), "");
    }

    #[test]
    fn test_write_gbk_format_single_gene() {
        let mut buffer = Vec::new();
        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 99, // 0-based
                end: 299,
                ..Default::default()
            },
            score: GeneScore {
                confidence: 85.75,
                ..Default::default()
            },
            ..Default::default()
        };
        let results = OrphosResults {
            genes: vec![gene],
            sequence_info: SequenceInfo {
                header: "sequence1".to_string(),
                length: 1000,
                description: Some("Test sequence".to_string()),
                gc_content: 12.0,
                num_genes: 25,
            },
            training_used: Training::default(),
            metagenomic_model: None,
        };

        let result = write_gbk_format(&mut buffer, &results);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("     CDS             100..300"));
        assert!(output.contains("/locus_tag=\"sequence1_1\""));
        assert!(output.contains("/confidence=85.75"));
    }

    #[test]
    fn test_write_gbk_format_multiple_genes() {
        let mut buffer = Vec::new();
        let genes = vec![
            Gene {
                coordinates: GeneCoordinates {
                    begin: 0,
                    end: 100,
                    ..Default::default()
                },
                score: GeneScore {
                    confidence: 90.50,
                    ..Default::default()
                },
                ..Default::default()
            },
            Gene {
                coordinates: GeneCoordinates {
                    begin: 200,
                    end: 400,
                    strand: Strand::Forward,
                    start_index: 0,
                    stop_index: 0,
                },
                score: GeneScore {
                    confidence: 75.25,
                    ..Default::default()
                },
                ..Default::default()
            },
        ];
        let results = OrphosResults {
            genes,
            sequence_info: SequenceInfo {
                header: "multi_gene_seq".to_string(),
                length: 500,
                description: Some("Test sequence".to_string()),
                gc_content: 12.0,
                num_genes: 25,
            },
            training_used: Training::default(),
            metagenomic_model: None,
        };

        let result = write_gbk_format(&mut buffer, &results);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("     CDS             1..101"));
        assert!(output.contains("/locus_tag=\"multi_gene_seq_1\""));
        assert!(output.contains("/confidence=90.50"));
        assert!(output.contains("     CDS             201..401"));
        assert!(output.contains("/locus_tag=\"multi_gene_seq_2\""));
        assert!(output.contains("/confidence=75.25"));
    }

    #[test]
    fn test_write_gbk_format_special_characters_in_header() {
        let mut buffer = Vec::new();
        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 10,
                end: 50,
                ..Default::default()
            },
            score: GeneScore {
                confidence: 95.0,
                ..Default::default()
            },
            ..Default::default()
        };
        let results = OrphosResults {
            genes: vec![gene],
            sequence_info: SequenceInfo {
                header: "seq_with-special.chars".to_string(),
                length: 100,
                description: Some("Test sequence".to_string()),
                gc_content: 12.0,
                num_genes: 25,
            },
            training_used: Training::default(),
            metagenomic_model: None,
        };

        let result = write_gbk_format(&mut buffer, &results);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("/locus_tag=\"seq_with-special.chars_1\""));
    }

    #[test]
    fn test_write_gbk_format_confidence_formatting() {
        let mut buffer = Vec::new();
        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 0,
                end: 10,
                strand: Strand::Forward,
                start_index: 2,
                stop_index: 10,
            },
            score: GeneScore {
                confidence: 99.999,
                total_score: 45.0,
                coding_score: 30.0,
                start_score: 15.0,
                ribosome_binding_score: 10.0,
                upstream_score: 5.0,
                type_score: 20.0,
            },
            ..Default::default()
        };
        let results = OrphosResults {
            genes: vec![gene],
            sequence_info: SequenceInfo {
                header: "test".to_string(),
                length: 20,
                description: Some("Test sequence".to_string()),
                gc_content: 12.0,
                num_genes: 25,
            },
            training_used: Training::default(),
            metagenomic_model: None,
        };

        let result = write_gbk_format(&mut buffer, &results);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("/confidence=100.00"));
    }
}
