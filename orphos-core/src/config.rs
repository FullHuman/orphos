/// Output format options for gene prediction results.
///
/// Orphos supports multiple output formats for compatibility with
/// different downstream analysis tools.
///
/// # Formats
///
/// - **GenBank**: Feature-rich annotation format with gene sequences
/// - **GFF**: General Feature Format version 3 (widely supported)
/// - **GCA**: Gene coordinate annotation (simple tabular format)
/// - **SCO**: Simple coordinate output (minimal format)
///
/// # Examples
///
/// ```rust
/// use orphos_core::config::{OutputFormat, OrphosConfig};
///
/// let config = OrphosConfig {
///     output_format: OutputFormat::Gff,
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// GenBank format output with full feature annotations and sequences.
    ///
    /// Includes gene coordinates, translation, product names, and sequence data.
    /// Compatible with NCBI submission tools.
    Genbank,

    /// Gene coordinate annotation format.
    ///
    /// Tab-delimited format with gene coordinates and basic metadata.
    /// Lightweight and easy to parse.
    Gca,

    /// Simple coordinate output format.
    ///
    /// Minimal format with just start/stop positions and strand.
    /// Useful for quick gene counting or position extraction.
    Sco,

    /// General Feature Format version 3.
    ///
    /// Standard genome annotation format supported by most bioinformatics tools.
    /// Includes gene coordinates, scores, and attributes.
    Gff,
}

/// Configuration settings for Orphos gene prediction analysis.
///
/// This struct controls all aspects of gene prediction including analysis mode,
/// sequence handling, and output formatting.
///
/// # Examples
///
/// ## Default configuration
///
/// ```rust
/// use orphos_core::config::OrphosConfig;
///
/// let config = OrphosConfig::default();
/// ```
///
/// ## Custom configuration for closed-ended genomes
///
/// ```rust
/// use orphos_core::config::{OrphosConfig, OutputFormat};
///
/// let config = OrphosConfig {
///     closed_ends: true,
///     mask_n_runs: true,
///     output_format: OutputFormat::Gff,
///     ..Default::default()
/// };
/// ```
///
/// ## Metagenomic mode with multiple threads
///
/// ```rust
/// use orphos_core::config::OrphosConfig;
///
/// let config = OrphosConfig {
///     metagenomic: true,
///     num_threads: Some(8),
///     quiet: true,
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct OrphosConfig {
    /// Enable metagenomic mode for fragmented sequences.
    ///
    /// When `true`, uses pre-computed models instead of training on each sequence.
    /// Recommended for:
    /// - Short contigs (< 100 kb)
    /// - Mixed community samples
    /// - Fragmented assemblies
    ///
    /// **Default**: `false` (single genome mode)
    pub metagenomic: bool,

    /// Treat sequences as having closed ends (complete genomes).
    ///
    /// When `true`, prevents genes from extending off sequence edges.
    /// Use for complete, circularized genomes.
    ///
    /// **Default**: `false` (allow edge genes)
    pub closed_ends: bool,

    /// Mask runs of N characters during analysis.
    ///
    /// When `true`, treats stretches of N's as gaps and prevents
    /// genes from spanning them. Useful for draft genomes with gaps.
    ///
    /// **Default**: `false`
    pub mask_n_runs: bool,

    /// Force use of non-Shine-Dalgarno models for start recognition.
    ///
    /// When `true`, disables detection of ribosome binding sites.
    /// Rarely needed except for organisms without canonical RBS.
    ///
    /// **Default**: `false` (auto-detect)
    pub force_non_sd: bool,

    /// Suppress informational output during processing.
    ///
    /// When `true`, prevents progress messages and statistics from
    /// being printed to stderr.
    ///
    /// **Default**: `false`
    pub quiet: bool,

    /// Output format for gene prediction results.
    ///
    /// Controls the format of generated output files. See [`OutputFormat`]
    /// for available options.
    ///
    /// **Default**: [`OutputFormat::Genbank`]
    pub output_format: OutputFormat,

    /// Genetic code translation table number (1-25).
    ///
    /// Specifies which genetic code to use for translation:
    /// - `11`: Bacterial/Archaeal (most common, default)
    /// - `4`: Mycoplasma/Spiroplasma
    /// - Others: See NCBI genetic code tables
    ///
    /// **Default**: `None` (auto-detect, usually table 11)
    pub translation_table: Option<u8>,

    /// Number of threads to use for parallel processing.
    ///
    /// When set, configures Rayon thread pool for parallel analysis
    /// of multiple sequences. Set to `None` for automatic detection.
    ///
    /// **Default**: `None` (use all available cores)
    pub num_threads: Option<usize>,
}

impl Default for OrphosConfig {
    fn default() -> Self {
        Self {
            metagenomic: false,
            closed_ends: false,
            mask_n_runs: false,
            force_non_sd: false,
            quiet: false,
            output_format: OutputFormat::Genbank,
            translation_table: None,
            num_threads: None,
        }
    }
}
