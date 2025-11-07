use bio::io::fasta;
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::io::Cursor;

use orphos_core::config::{OrphosConfig, OutputFormat};
use orphos_core::engine::OrphosAnalyzer;
use orphos_core::output::write_results;

/// Options for configuring Orphos gene prediction
#[pyclass]
#[derive(Clone)]
pub struct OrphosOptions {
    #[pyo3(get, set)]
    /// Operating mode: "single" for single genome mode, "meta" for metagenomic mode
    pub mode: String,

    #[pyo3(get, set)]
    /// Output format: "gbk" (GenBank), "gff" (GFF3), "sco" (simple coordinates), "gca" (gene calls)
    pub format: String,

    #[pyo3(get, set)]
    /// Closed ends - don't allow genes to run off edges
    pub closed_ends: bool,

    #[pyo3(get, set)]
    /// Mask runs of N's in the sequence
    pub mask_n_runs: bool,

    #[pyo3(get, set)]
    /// Force the use of non-Shine-Dalgarno gene model
    pub force_non_sd: bool,

    #[pyo3(get, set)]
    /// Translation table (1-25, excluding 7, 8, and 17-20)
    pub translation_table: Option<u8>,

    #[pyo3(get, set)]
    /// Number of threads to use (None for default)
    pub num_threads: Option<usize>,

    #[pyo3(get, set)]
    /// Suppress informational output
    pub quiet: bool,
}

#[pymethods]
impl OrphosOptions {
    #[new]
    #[pyo3(signature = (
        mode="single",
        format="gbk",
        closed_ends=false,
        mask_n_runs=false,
        force_non_sd=false,
        translation_table=None,
        num_threads=None,
        quiet=true
    ))]
    fn new(
        mode: &str,
        format: &str,
        closed_ends: bool,
        mask_n_runs: bool,
        force_non_sd: bool,
        translation_table: Option<u8>,
        num_threads: Option<usize>,
        quiet: bool,
    ) -> PyResult<Self> {
        Ok(OrphosOptions {
            mode: mode.to_string(),
            format: format.to_string(),
            closed_ends,
            mask_n_runs,
            force_non_sd,
            translation_table,
            num_threads,
            quiet,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "OrphosOptions(mode='{}', format='{}', closed_ends={}, mask_n_runs={}, force_non_sd={}, translation_table={:?}, num_threads={:?}, quiet={})",
            self.mode, self.format, self.closed_ends, self.mask_n_runs,
            self.force_non_sd, self.translation_table, self.num_threads, self.quiet
        )
    }
}

/// Result from Orphos gene prediction
#[pyclass]
pub struct OrphosResult {
    #[pyo3(get)]
    /// The formatted output (GenBank, GFF, etc.)
    pub output: String,

    #[pyo3(get)]
    /// Total number of genes predicted
    pub gene_count: usize,

    #[pyo3(get)]
    /// Number of sequences analyzed
    pub sequence_count: usize,
}

#[pymethods]
impl OrphosResult {
    fn __repr__(&self) -> String {
        format!(
            "OrphosResult(gene_count={}, sequence_count={}, output_length={})",
            self.gene_count,
            self.sequence_count,
            self.output.len()
        )
    }
}

/// Parse FASTA content from a string
fn parse_fasta_string(content: &str) -> PyResult<Vec<(String, Option<String>, Vec<u8>)>> {
    let cursor = Cursor::new(content.as_bytes());
    let reader = fasta::Reader::new(cursor);
    let mut sequences = Vec::new();

    for result in reader.records() {
        let record =
            result.map_err(|e| PyIOError::new_err(format!("FASTA parsing error: {}", e)))?;

        let id = record.id().to_string();
        let description = record.desc().map(String::from);
        let seq = record.seq().to_vec();
        sequences.push((id, description, seq));
    }

    if sequences.is_empty() {
        return Err(PyValueError::new_err("No sequences found in FASTA input"));
    }

    Ok(sequences)
}

/// Convert OrphosOptions to OrphosConfig
fn options_to_config(options: &OrphosOptions) -> PyResult<OrphosConfig> {
    // Validate translation table
    if let Some(tt) = options.translation_table {
        if !(1..=25).contains(&tt) || tt == 7 || tt == 8 || (17..=20).contains(&tt) {
            return Err(PyValueError::new_err(
                "Invalid translation table. Must be 1-25, excluding 7, 8, and 17-20",
            ));
        }
    }

    // Parse output format
    let output_format = match options.format.as_str() {
        "gbk" | "genbank" => OutputFormat::Genbank,
        "gff" => OutputFormat::Gff,
        "sco" => OutputFormat::Sco,
        "gca" => OutputFormat::Gca,
        _ => {
            return Err(PyValueError::new_err(
                "Invalid output format. Must be one of: gbk, genbank, gff, sco, gca",
            ))
        }
    };

    // Validate mode
    if options.mode != "single" && options.mode != "meta" {
        return Err(PyValueError::new_err(
            "Invalid mode. Must be 'single' or 'meta'",
        ));
    }

    Ok(OrphosConfig {
        metagenomic: options.mode == "meta",
        closed_ends: options.closed_ends,
        mask_n_runs: options.mask_n_runs,
        force_non_sd: options.force_non_sd,
        quiet: options.quiet,
        output_format,
        translation_table: options.translation_table,
        num_threads: options.num_threads,
    })
}

/// Analyze DNA sequences for gene prediction
///
/// Args:
///     fasta_content (str): FASTA-formatted sequence(s) as a string
///     options (OrphosOptions, optional): Configuration options for the analysis
///
/// Returns:
///     OrphosResult: Object containing the formatted output and statistics
///
/// Example:
///     >>> import prodigal
///     >>> fasta = ">seq1\\nATGCGATCGATCGATCG..."
///     >>> result = prodigal.analyze_sequence(fasta)
///     >>> print(result.gene_count)
#[pyfunction]
#[pyo3(signature = (fasta_content, options=None))]
fn analyze_sequence(fasta_content: &str, options: Option<OrphosOptions>) -> PyResult<OrphosResult> {
    // Use default options if none provided
    let options = options.unwrap_or_else(|| OrphosOptions {
        mode: "single".to_string(),
        format: "gbk".to_string(),
        closed_ends: false,
        mask_n_runs: false,
        force_non_sd: false,
        translation_table: None,
        num_threads: None,
        quiet: true,
    });

    // Convert options to config
    let config = options_to_config(&options)?;

    // Parse FASTA content
    let sequences = parse_fasta_string(fasta_content)?;

    // Run Orphos analysis
    let mut analyzer = OrphosAnalyzer::new(config);
    let mut all_results = Vec::new();

    for (header, description, seq_bytes) in sequences {
        let result = analyzer
            .analyze_sequence_bytes(&seq_bytes, header, description)
            .map_err(|e| PyValueError::new_err(format!("Analysis error: {}", e)))?;
        all_results.push(result);
    }

    // Generate output
    let mut output = Vec::new();
    for result in &all_results {
        write_results(&mut output, result, analyzer.config.output_format)
            .map_err(|e| PyIOError::new_err(format!("Output error: {}", e)))?;
    }

    let output_str = String::from_utf8(output)
        .map_err(|e| PyValueError::new_err(format!("UTF-8 conversion error: {}", e)))?;

    let gene_count = all_results.iter().map(|r| r.genes.len()).sum();
    let sequence_count = all_results.len();

    Ok(OrphosResult {
        output: output_str,
        gene_count,
        sequence_count,
    })
}

/// Analyze DNA sequences from a FASTA file
///
/// Args:
///     file_path (str): Path to the FASTA file
///     options (OrphosOptions, optional): Configuration options for the analysis
///
/// Returns:
///     OrphosResult: Object containing the formatted output and statistics
///
/// Example:
///     >>> import prodigal
///     >>> result = prodigal.analyze_file("genome.fasta")
///     >>> print(result.gene_count)
#[pyfunction]
#[pyo3(signature = (file_path, options=None))]
fn analyze_file(file_path: &str, options: Option<OrphosOptions>) -> PyResult<OrphosResult> {
    // Read file content
    let fasta_content = std::fs::read_to_string(file_path)
        .map_err(|e| PyIOError::new_err(format!("Failed to read file '{}': {}", file_path, e)))?;

    // Use the analyze_sequence function
    analyze_sequence(&fasta_content, options)
}

/// Orphos - Gene prediction for microbial genomes
///
/// This module provides Python bindings for Orphos, a tool for finding
/// protein-coding genes in bacterial and archaeal genomes.
#[pymodule]
fn orphos(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<OrphosOptions>()?;
    m.add_class::<OrphosResult>()?;
    m.add_function(wrap_pyfunction!(analyze_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_file, m)?)?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", "Orphos gene prediction for microbial genomes")?;

    Ok(())
}
