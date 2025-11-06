//! # Orphos CLI - Command-Line Gene Finder
//!
//! A command-line interface for the Orphos prokaryotic gene finding algorithm.
//!
//! ## Usage
//!
//! ```bash
//! # Basic gene prediction
//! orphos-cli -i genome.fasta -o genes.gbk
//!
//! # Output in GFF format
//! orphos-cli -i genome.fasta -f gff -o genes.gff
//!
//! # Metagenomic mode for short contigs
//! orphos-cli -i contigs.fasta -p meta -o genes.gff
//!
//! # Closed ends for complete genomes
//! orphos-cli -i genome.fasta -c -o genes.gbk
//! ```
//!
//! ## Options
//!
//! - `-i, --input <FILE>`: Input FASTA file (default: stdin)
//! - `-o, --output <FILE>`: Output file (default: stdout)
//! - `-f, --format <FORMAT>`: Output format: gbk, gff, sco, gca (default: gbk)
//! - `-p, --mode <MODE>`: Analysis mode: single or meta (default: single)
//! - `-c, --closed`: Closed ends (no genes off edges)
//! - `-m, --mask`: Mask runs of N's
//! - `-q, --quiet`: Suppress progress messages
//! - `-g, --translation-table <TABLE>`: Translation table 1-25 (default: auto)
//!
//! ## Examples
//!
//! ### Single Genome Analysis
//!
//! ```bash
//! orphos-cli -i ecoli.fasta -f gff -o ecoli.gff
//! ```
//!
//! ### Metagenomic Assembly
//!
//! ```bash
//! orphos-cli -i metagenome.fasta -p meta -f gff -o predictions.gff
//! ```
//!
//! ### Complete Circular Genome
//!
//! ```bash
//! orphos-cli -i plasmid.fasta -c -o plasmid.gbk
//! ```

use clap::{Arg, Command};
use orphos_core::config::{OutputFormat, OrphosConfig};
use orphos_core::output::write_results;
use orphos_core::*;
use std::fs::File;
use std::io::{self, BufWriter, Write};

/// Main entry point for the Orphos CLI application.
///
/// Parses command-line arguments, configures Orphos, analyzes input sequences,
/// and writes results in the requested format.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("orphos")
        .version(env!("CARGO_PKG_VERSION"))
        .about("Prokaryotic gene finding algorithm")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("FILE")
                .help("Input FASTA file (default: stdin)"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output file (default: stdout)"),
        )
        .arg(
            Arg::new("format")
                .short('f')
                .long("format")
                .value_name("FORMAT")
                .help("Output format: gbk, gff, sco, gca")
                .default_value("gbk"),
        )
        .arg(
            Arg::new("mode")
                .short('p')
                .long("mode")
                .value_name("MODE")
                .help("Analysis mode: single or meta")
                .default_value("single"),
        )
        .arg(
            Arg::new("closed")
                .short('c')
                .long("closed")
                .help("Closed ends (no genes off edges)"),
        )
        .arg(
            Arg::new("mask")
                .short('m')
                .long("mask")
                .help("Mask runs of N's"),
        )
        .arg(
            Arg::new("quiet")
                .short('q')
                .long("quiet")
                .help("Quiet mode"),
        )
        .arg(
            Arg::new("training")
                .short('t')
                .long("training")
                .value_name("FILE")
                .help("Training file"),
        )
        .arg(
            Arg::new("translation-table")
                .short('g')
                .long("translation-table")
                .value_name("TABLE")
                .help("Translation table (1-25)"),
        )
        .get_matches();

    // Parse options
    let mut options = OrphosConfig {
        metagenomic: matches.get_one::<String>("mode").map(|s| s.as_str()) == Some("meta"),
        closed_ends: matches.contains_id("closed"),
        mask_n_runs: matches.contains_id("mask"),
        quiet: matches.contains_id("quiet"),
        // training_file: matches.get_one::<String>("training").cloned(),
        ..Default::default()
    };

    if let Some(tt_str) = matches.get_one::<String>("translation-table") {
        let tt: u8 = tt_str
            .parse()
            .map_err(|_| "Invalid translation table number")?;
        if !(1..=25).contains(&tt) || tt == 7 || tt == 8 || (17..=20).contains(&tt) {
            return Err("Invalid translation table specified".into());
        }
        options.translation_table = Some(tt);
    }

    options.output_format = match matches.get_one::<String>("format").unwrap().as_str() {
        "gbk" | "genbank" => OutputFormat::Genbank,
        "gff" => OutputFormat::Gff,
        "sco" => OutputFormat::Sco,
        "gca" => OutputFormat::Gca,
        _ => return Err("Invalid output format".into()),
    };

    let mut orphos = OrphosAnalyzer::new(options);
    let results = if let Some(input_file) = matches.get_one::<String>("input") {
        orphos.analyze_fasta_file(input_file)?
    } else {
        unimplemented!("Reading from stdin is not yet implemented");
    };

    // Write output
    let mut writer: Box<dyn Write> = if let Some(output_file) = matches.get_one::<String>("output")
    {
        Box::new(BufWriter::new(File::create(output_file)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    for result in &results {
        write_results(&mut writer, result, orphos.config.output_format)?;
    }

    if !matches.contains_id("quiet") {
        eprintln!(
            "Analysis complete! Found {} genes in {} sequences.",
            results.iter().map(|r| r.genes.len()).sum::<usize>(),
            results.len()
        );
        // Debug: show total dicodons if available
        let total_dicodons: u32 = results.iter().map(|r| r.training_used.total_dicodons).sum();
        if total_dicodons > 0 {
            eprintln!("Total dicodons counted: {}", total_dicodons);
        }
    }

    Ok(())
}
