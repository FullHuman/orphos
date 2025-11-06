use std::io::Write;

use bio::bio_types::strand::Strand;

use crate::{OrphosError, results::OrphosResults};

/// Write results in SCO format (simple coordinate output)
pub fn write_sco_format<W: Write>(
    writer: &mut W,
    results: &OrphosResults,
) -> Result<(), OrphosError> {
    for gene in &results.genes {
        let strand_num = match gene.coordinates.strand {
            Strand::Forward => 1,
            Strand::Reverse => -1,
            Strand::Unknown => 0,
        };

        writeln!(
            writer,
            "{}\t{}\t{}\t{:.2}",
            gene.coordinates.begin + 1,
            gene.coordinates.end + 1,
            strand_num,
            gene.score.confidence
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        results::{OrphosResults, SequenceInfo},
        types::{Gene, GeneCoordinates, GeneScore, Training},
    };
    use bio::bio_types::strand::Strand;
    use std::io::Cursor;

    fn create_test_results_with_genes(genes: Vec<Gene>) -> OrphosResults {
        OrphosResults {
            sequence_info: SequenceInfo {
                header: "test_seq".to_string(),
                length: 1000,
                description: Some("Test sequence".to_string()),
                gc_content: 50.0,
                num_genes: genes.len(),
            },
            genes,
            training_used: Training::default(),
            metagenomic_model: None,
        }
    }

    #[test]
    fn test_write_sco_format_single_gene_forward() {
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
            ..Default::default()
        };

        let results = create_test_results_with_genes(vec![gene]);
        write_sco_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "100\t300\t1\t95.50\n");
    }

    #[test]
    fn test_write_sco_format_single_gene_reverse() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 199,
                end: 399,
                strand: Strand::Reverse,
                ..Default::default()
            },
            score: GeneScore {
                confidence: 87.25,
                ..Default::default()
            },
            ..Default::default()
        };

        let results = create_test_results_with_genes(vec![gene]);
        write_sco_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "200\t400\t-1\t87.25\n");
    }

    #[test]
    fn test_write_sco_format_unknown_strand() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 299,
                end: 599,
                strand: Strand::Unknown,
                ..Default::default()
            },
            score: GeneScore {
                confidence: 42.0,
                ..Default::default()
            },
            ..Default::default()
        };

        let results = create_test_results_with_genes(vec![gene]);
        write_sco_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "300\t600\t0\t42.00\n");
    }

    #[test]
    fn test_write_sco_format_multiple_genes() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let genes = vec![
            Gene {
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
                ..Default::default()
            },
            Gene {
                coordinates: GeneCoordinates {
                    begin: 399,
                    end: 599,
                    strand: Strand::Reverse,
                    ..Default::default()
                },
                score: GeneScore {
                    confidence: 88.75,
                    ..Default::default()
                },
                ..Default::default()
            },
            Gene {
                coordinates: GeneCoordinates {
                    begin: 699,
                    end: 899,
                    strand: Strand::Unknown,
                    ..Default::default()
                },
                score: GeneScore {
                    confidence: 12.34,
                    ..Default::default()
                },
                ..Default::default()
            },
        ];

        let results = create_test_results_with_genes(genes);
        write_sco_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        let expected = "100\t300\t1\t95.50\n400\t600\t-1\t88.75\n700\t900\t0\t12.34\n";
        assert_eq!(output, expected);
    }

    #[test]
    fn test_write_sco_format_no_genes() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let results = create_test_results_with_genes(vec![]);
        write_sco_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "");
    }

    #[test]
    fn test_write_sco_format_confidence_precision() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let gene = Gene {
            coordinates: GeneCoordinates {
                begin: 0,
                end: 100,
                strand: Strand::Forward,
                ..Default::default()
            },
            score: GeneScore {
                confidence: 99.999,
                ..Default::default()
            },
            ..Default::default()
        };

        let results = create_test_results_with_genes(vec![gene]);
        write_sco_format(&mut cursor, &results).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "1\t101\t1\t100.00\n");
    }
}
