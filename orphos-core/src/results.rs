use crate::types::{Gene, Training};

/// Gene finding results from Orphos analysis.
///
/// Contains all information from a gene prediction run including
/// predicted genes, training parameters, and sequence statistics.
///
/// # Fields
///
/// - `genes`: Vector of predicted genes sorted by position
/// - `training_used`: Statistical model parameters from training
/// - `sequence_info`: Metadata about the analyzed sequence
/// - `metagenomic_model`: Model name if metagenomic mode was used
///
/// # Examples
///
/// ```rust,no_run
/// use orphos_core::{OrphosAnalyzer, config::OrphosConfig, config::OutputFormat};
/// use orphos_core::output::write_results;
///
/// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
/// let results = analyzer.analyze_sequence("ATGCGATCG...", None)?;
///
/// println!("Sequence: {}", results.sequence_info.header);
/// println!("Length: {} bp", results.sequence_info.length);
/// println!("GC%: {:.2}", results.sequence_info.gc_content * 100.0);
/// println!("Genes: {}", results.genes.len());
///
/// // Write results to file
/// let mut output = std::fs::File::create("output.gbk")?;
/// write_results(&mut output, &results, OutputFormat::Genbank)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct OrphosResults {
    /// Vector of predicted genes sorted by genomic position.
    ///
    /// Each gene contains coordinates, strand, scores, and sequence information.
    pub genes: Vec<Gene>,

    /// Training parameters used for gene prediction.
    ///
    /// Includes statistical models for start codons, RBS motifs, and codon usage.
    pub training_used: Training,

    /// Information about the analyzed sequence.
    ///
    /// Contains length, GC content, gene count, and sequence identifiers.
    pub sequence_info: SequenceInfo,

    /// Name of metagenomic model used, if applicable.
    ///
    /// Set to `Some("Best")` in metagenomic mode, `None` in single genome mode.
    pub metagenomic_model: Option<String>,
}

/// Information about a processed sequence.
///
/// Contains metadata and statistics for a sequence that was analyzed.
///
/// # Examples
///
/// ```rust,no_run
/// # use orphos_core::results::SequenceInfo;
/// let info = SequenceInfo {
///     length: 4_641_652,
///     gc_content: 0.5079,
///     num_genes: 4321,
///     header: "E. coli K12".to_string(),
///     description: Some("Complete genome".to_string()),
/// };
///
/// println!("{}: {} bp, {:.2}% GC, {} genes",
///          info.header,
///          info.length,
///          info.gc_content * 100.0,
///          info.num_genes);
/// ```
#[derive(Debug, Clone)]
pub struct SequenceInfo {
    /// Length of the sequence in base pairs.
    pub length: usize,

    /// GC content as a fraction (0.0 to 1.0).
    ///
    /// Multiply by 100 to get percentage.
    pub gc_content: f64,

    /// Number of genes predicted in the sequence.
    pub num_genes: usize,

    /// Sequence identifier from FASTA header.
    ///
    /// The first word of the FASTA header line (after '>').
    pub header: String,

    /// Full sequence description from FASTA header.
    ///
    /// Everything after the first word in the FASTA header line.
    pub description: Option<String>,
}
