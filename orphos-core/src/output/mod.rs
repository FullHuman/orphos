//! Output formatting for gene prediction results.
//!
//! This module provides writers for converting [`OrphosResults`] into various
//! standard bioinformatics file formats.
//!
//! ## Supported Formats
//!
//! - **GenBank (GBK)**: Feature-rich annotation format with sequences
//! - **GFF3**: General Feature Format version 3
//! - **GCA**: Gene coordinate annotation (tabular)
//! - **SCO**: Simple coordinate output
//!
//! ## Examples
//!
//! ### Write results to a file
//!
//! ```rust,no_run
//! use orphos_core::{OrphosAnalyzer, config::{OrphosConfig, OutputFormat}};
//! use orphos_core::output::write_results;
//! use std::fs::File;
//!
//! let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
//! let results = analyzer.analyze_sequence("ATGCGATCG...", None)?;
//!
//! // Write as GenBank
//! let mut gbk_file = File::create("output.gbk")?;
//! write_results(&mut gbk_file, &results, OutputFormat::Genbank)?;
//!
//! // Write as GFF3
//! let mut gff_file = File::create("output.gff")?;
//! write_results(&mut gff_file, &results, OutputFormat::Gff)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ### Write to stdout
//!
//! ```rust,no_run
//! use orphos_core::{OrphosAnalyzer, config::{OrphosConfig, OutputFormat}};
//! use orphos_core::output::write_results;
//! use std::io::stdout;
//!
//! let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
//! let results = analyzer.analyze_sequence("ATGCGATCG...", None)?;
//!
//! write_results(&mut stdout(), &results, OutputFormat::Gff)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::{OrphosError, config::OutputFormat, results::OrphosResults};
use std::io::Write;

mod formats {
    pub mod gbk;
    pub mod gca;
    pub mod gff;
    pub mod sco;
}

use formats::{
    gbk::write_gbk_format, gca::write_gca_format, gff::write_gff_format, sco::write_sco_format,
};

/// Writes gene prediction results in the specified format.
///
/// This is the main entry point for output formatting. It delegates to
/// format-specific writers based on the requested output format.
///
/// # Arguments
///
/// * `writer` - Output writer (file, stdout, buffer, etc.)
/// * `results` - Gene prediction results to write
/// * `format` - Desired output format
///
/// # Errors
///
/// Returns [`OrphosError`] if writing fails (I/O errors, invalid data).
///
/// # Examples
///
/// ```rust,no_run
/// use orphos_core::{OrphosAnalyzer, config::{OrphosConfig, OutputFormat}};
/// use orphos_core::output::write_results;
/// use std::fs::File;
///
/// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
/// let results = analyzer.analyze_fasta_file("genome.fasta")?;
///
/// let mut output = File::create("genes.gff")?;
/// for result in &results {
///     write_results(&mut output, result, OutputFormat::Gff)?;
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn write_results<W: Write>(
    writer: &mut W,
    results: &OrphosResults,
    format: OutputFormat,
) -> Result<(), OrphosError> {
    match format {
        OutputFormat::Genbank => write_gbk_format(writer, results),
        OutputFormat::Gff => write_gff_format(writer, results),
        OutputFormat::Sco => write_sco_format(writer, results),
        OutputFormat::Gca => write_gca_format(writer, results),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        config::OutputFormat,
        results::{OrphosResults, SequenceInfo},
        types::{Gene, GeneCoordinates, GeneScore, Training},
    };
    use bio::bio_types::strand::Strand;
    use std::io::Cursor;

    fn create_test_results() -> OrphosResults {
        OrphosResults {
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
                score: GeneScore {
                    confidence: 95.5,
                    ..Default::default()
                },
                ..Default::default()
            }],
            training_used: Training::default(),
            metagenomic_model: None,
        }
    }

    #[test]
    fn test_write_results_genbank_format() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);
        let results = create_test_results();

        let result = write_results(&mut cursor, &results, OutputFormat::Genbank);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("CDS"));
        assert!(output.contains("100..300"));
        assert!(output.contains("test_seq_1"));
        assert!(output.contains("confidence=95.50"));
    }

    #[test]
    fn test_write_results_gff_format() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);
        let results = create_test_results();

        let result = write_results(&mut cursor, &results, OutputFormat::Gff);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("##gff-version"));
        assert!(output.contains("test_seq"));
        assert!(output.contains("CDS"));
        assert!(output.contains("conf=95.50"));
    }

    #[test]
    fn test_write_results_sco_format() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);
        let results = create_test_results();

        let result = write_results(&mut cursor, &results, OutputFormat::Sco);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output, "100\t300\t1\t95.50\n");
    }

    #[test]
    fn test_write_results_gca_format() {
        let mut buffer = Vec::new();
        let mut cursor = Cursor::new(&mut buffer);
        let results = create_test_results();

        let result = write_results(&mut cursor, &results, OutputFormat::Gca);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains(">test_seq_1"));
        assert!(output.contains("100\t300\t+"));
    }

    #[test]
    fn test_write_results_format_consistency() {
        let results = create_test_results();
        let formats = vec![
            OutputFormat::Genbank,
            OutputFormat::Gff,
            OutputFormat::Sco,
            OutputFormat::Gca,
        ];

        for format in formats {
            let mut buffer = Vec::new();
            let mut cursor = Cursor::new(&mut buffer);

            let result = write_results(&mut cursor, &results, format);
            assert!(result.is_ok(), "Failed to write format: {:?}", format);

            let output = String::from_utf8(buffer).unwrap();
            assert!(!output.is_empty(), "Empty output for format: {:?}", format);
        }
    }

    #[test]
    fn test_write_results_empty_genes() {
        let results = OrphosResults {
            sequence_info: SequenceInfo {
                header: "empty_seq".to_string(),
                length: 100,
                description: None,
                gc_content: 40.0,
                num_genes: 0,
            },
            genes: vec![],
            training_used: Training::default(),
            metagenomic_model: None,
        };

        // Test all formats with empty gene list
        let formats = vec![
            OutputFormat::Genbank,
            OutputFormat::Gff,
            OutputFormat::Sco,
            OutputFormat::Gca,
        ];

        for format in formats {
            let mut buffer = Vec::new();
            let mut cursor = Cursor::new(&mut buffer);

            let result = write_results(&mut cursor, &results, format);
            assert!(
                result.is_ok(),
                "Failed to write empty results for format: {:?}",
                format
            );
        }
    }
}
