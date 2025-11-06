use std::io::Write;

use bio::bio_types::strand::Strand;

use crate::{results::OrphosResults, types::OrphosError};

/// Write results in GCA format
pub fn write_gca_format<W: Write>(
    writer: &mut W,
    results: &OrphosResults,
) -> Result<(), OrphosError> {
    for (i, gene) in results.genes.iter().enumerate() {
        writeln!(writer, ">{}_{}", results.sequence_info.header, i + 1)?;
        writeln!(
            writer,
            "{}\t{}\t{}",
            gene.coordinates.begin + 1,
            gene.coordinates.end + 1,
            match gene.coordinates.strand {
                Strand::Forward => "+",
                Strand::Reverse => "-",
                Strand::Unknown => "?",
            }
        )?;
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use crate::{
        results::SequenceInfo,
        types::{Gene, GeneCoordinates, Training},
    };

    use super::*;

    #[test]
    fn test_write_gca_format_single_gene_forward() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let results = OrphosResults {
            sequence_info: SequenceInfo {
                header: "test_seq".to_string(),
                length: 1000,
                description: Some("Test sequence".to_string()),
                gc_content: 50.0,
                num_genes: 1,
            },
            genes: vec![Gene {
                coordinates: GeneCoordinates {
                    begin: 99,
                    end: 299,
                    strand: Strand::Forward,
                    ..Default::default()
                },
                ..Default::default()
            }],
            training_used: Training::default(),
            metagenomic_model: None,
        };

        write_gca_format(&mut cursor, &results).unwrap();
        let output = String::from_utf8(buffer).unwrap();

        assert_eq!(output, ">test_seq_1\n100\t300\t+\n");
    }

    #[test]
    fn test_write_gca_format_single_gene_reverse() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let results = OrphosResults {
            sequence_info: SequenceInfo {
                header: "test_seq".to_string(),
                length: 1000,
                description: Some("Test sequence".to_string()),
                gc_content: 50.0,
                num_genes: 1,
            },
            genes: vec![Gene {
                coordinates: GeneCoordinates {
                    begin: 199,
                    end: 399,
                    strand: Strand::Reverse,
                    ..Default::default()
                },
                ..Default::default()
            }],
            training_used: Training::default(),
            metagenomic_model: None,
        };

        write_gca_format(&mut cursor, &results).unwrap();
        let output = String::from_utf8(buffer).unwrap();

        assert_eq!(output, ">test_seq_1\n200\t400\t-\n");
    }

    #[test]
    fn test_write_gca_format_multiple_genes() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let results = OrphosResults {
            sequence_info: SequenceInfo {
                header: "multi_gene_seq".to_string(),
                length: 2000,
                description: Some("Test sequence with multiple genes".to_string()),
                gc_content: 60.0,
                num_genes: 3,
            },
            genes: vec![
                Gene {
                    coordinates: GeneCoordinates {
                        begin: 0,
                        end: 299,
                        strand: Strand::Forward,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                Gene {
                    coordinates: GeneCoordinates {
                        begin: 499,
                        end: 799,
                        strand: Strand::Reverse,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                Gene {
                    coordinates: GeneCoordinates {
                        begin: 999,
                        end: 1199,
                        strand: Strand::Unknown,
                        ..Default::default()
                    },
                    ..Default::default()
                },
            ],
            training_used: Training::default(),
            metagenomic_model: None,
        };

        write_gca_format(&mut cursor, &results).unwrap();
        let output = String::from_utf8(buffer).unwrap();

        assert_eq!(
            output,
            ">multi_gene_seq_1\n1\t300\t+\n>multi_gene_seq_2\n500\t800\t-\n>multi_gene_seq_3\n1000\t1200\t?\n"
        );
    }

    #[test]
    fn test_write_gca_format_no_genes() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let results = OrphosResults {
            sequence_info: SequenceInfo {
                header: "empty_seq".to_string(),
                length: 100,
                description: Some("Empty sequence".to_string()),
                gc_content: 0.0,
                num_genes: 0,
            },
            genes: vec![],
            training_used: Training::default(),
            metagenomic_model: None,
        };

        write_gca_format(&mut cursor, &results).unwrap();
        let output = String::from_utf8(buffer).unwrap();

        assert_eq!(output, "");
    }

    #[test]
    fn test_write_gca_format_unknown_strand() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);

        let results = OrphosResults {
            sequence_info: SequenceInfo {
                header: "unknown_strand".to_string(),
                length: 500,
                description: Some("Test sequence with unknown strand".to_string()),
                gc_content: 0.0,
                num_genes: 1,
            },
            genes: vec![Gene {
                coordinates: GeneCoordinates {
                    begin: 49,
                    end: 149,
                    strand: Strand::Unknown,
                    ..Default::default()
                },
                ..Default::default()
            }],
            training_used: Training::default(),
            metagenomic_model: None,
        };

        write_gca_format(&mut cursor, &results).unwrap();
        let output = String::from_utf8(buffer).unwrap();

        assert_eq!(output, ">unknown_strand_1\n50\t150\t?\n");
    }
}
