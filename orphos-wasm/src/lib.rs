use orphos_core::config::{OutputFormat, OrphosConfig};
use orphos_core::engine::OrphosAnalyzer;
use orphos_core::output::write_results;
use serde::{Deserialize, Serialize};
use std::io::Cursor;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[derive(Serialize, Deserialize)]
pub struct WasmOrphosOptions {
    pub mode: String,                  // "single" or "meta"
    pub format: String,                // "gbk", "gff", "sco", "gca"
    pub closed_ends: bool,             // Closed ends (no genes off edges)
    pub mask_n_runs: bool,             // Mask runs of N's
    pub force_non_sd: bool,            // Force non-Shine-Dalgarno
    pub translation_table: Option<u8>, // Translation table (1-25)
}

impl Default for WasmOrphosOptions {
    fn default() -> Self {
        Self {
            mode: "single".to_string(),
            format: "gbk".to_string(),
            closed_ends: false,
            mask_n_runs: false,
            force_non_sd: false,
            translation_table: None,
        }
    }
}

#[wasm_bindgen]
pub struct OrphosResult {
    output: String,
    gene_count: usize,
    sequence_count: usize,
}

#[wasm_bindgen]
impl OrphosResult {
    #[wasm_bindgen(getter)]
    pub fn output(&self) -> String {
        self.output.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn gene_count(&self) -> usize {
        self.gene_count
    }

    #[wasm_bindgen(getter)]
    pub fn sequence_count(&self) -> usize {
        self.sequence_count
    }
}

// Type alias for FASTA records
type FastaRecord = (String, Option<String>, Vec<u8>);

// Helper function to parse FASTA from string
fn parse_fasta_string(content: &str) -> Result<Vec<FastaRecord>, String> {
    use bio::io::fasta;

    let cursor = Cursor::new(content.as_bytes());
    let reader = fasta::Reader::new(cursor);
    let mut sequences = Vec::new();

    for result in reader.records() {
        let record = result.map_err(|e| format!("FASTA parsing error: {}", e))?;
        let id = record.id().to_string();
        let description = record.desc().map(String::from);
        let seq = record.seq().to_vec();
        sequences.push((id, description, seq));
    }

    if sequences.is_empty() {
        return Err("No sequences found in FASTA input".to_string());
    }

    Ok(sequences)
}

#[wasm_bindgen]
pub fn analyze_sequence(
    fasta_content: &str,
    options_js: JsValue,
) -> Result<OrphosResult, JsValue> {
    // Parse options from JavaScript
    let wasm_options: WasmOrphosOptions =
        serde_wasm_bindgen::from_value(options_js).unwrap_or_default();

    // Validate translation table
    if let Some(tt) = wasm_options.translation_table {
        if !(1..=25).contains(&tt) || tt == 7 || tt == 8 || (17..=20).contains(&tt) {
            return Err(JsValue::from_str("Invalid translation table specified"));
        }
    }

    // Convert to OrphosConfig
    let output_format = match wasm_options.format.as_str() {
        "gbk" | "genbank" => OutputFormat::Genbank,
        "gff" => OutputFormat::Gff,
        "sco" => OutputFormat::Sco,
        "gca" => OutputFormat::Gca,
        _ => return Err(JsValue::from_str("Invalid output format")),
    };

    let config = OrphosConfig {
        metagenomic: wasm_options.mode == "meta",
        closed_ends: wasm_options.closed_ends,
        mask_n_runs: wasm_options.mask_n_runs,
        force_non_sd: wasm_options.force_non_sd,
        quiet: true,
        output_format,
        translation_table: wasm_options.translation_table,
        num_threads: None,
    };

    // Parse FASTA content
    let sequences = parse_fasta_string(fasta_content).map_err(|e| JsValue::from_str(&e))?;

    // Run Orphos analysis
    let mut analyzer = OrphosAnalyzer::new(config);
    let mut all_results = Vec::new();

    for (header, description, seq_bytes) in sequences {
        let result = analyzer
            .analyze_sequence_bytes(&seq_bytes, header, description)
            .map_err(|e| JsValue::from_str(&format!("Analysis error: {}", e)))?;
        all_results.push(result);
    }

    // Generate output
    let mut output = Vec::new();
    for result in &all_results {
        write_results(&mut output, result, analyzer.config.output_format)
            .map_err(|e| JsValue::from_str(&format!("Output error: {}", e)))?;
    }

    let output_str =
        String::from_utf8(output).map_err(|e| JsValue::from_str(&format!("UTF-8 error: {}", e)))?;

    let gene_count = all_results.iter().map(|r| r.genes.len()).sum();
    let sequence_count = all_results.len();

    Ok(OrphosResult {
        output: output_str,
        gene_count,
        sequence_count,
    })
}
