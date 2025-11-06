# Orphos Core

[![Crates.io](https://img.shields.io/crates/v/orphos-core.svg)](https://crates.io/crates/orphos-core)
[![Documentation](https://docs.rs/orphos-core/badge.svg)](https://docs.rs/orphos-core)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Core library for **Orphos**, a high-performance Rust implementation of Prodigal (prokaryotic gene prediction algorithms). This crate provides the foundational gene-finding capabilities for identifying protein-coding genes in microbial genomes.

## Overview

`orphos-core` implements an unsupervised machine learning approach for finding genes in prokaryotic genomes. It uses dynamic programming and statistical models trained on genomic features to predict gene locations with high accuracy.

### Key Features

- 🚀 **High Performance**: Multi-threaded processing using Rayon for parallel analysis
- 🔒 **Type Safety**: Compile-time guarantees using type-state pattern for training states
- 🧬 **Dual Modes**: Single genome mode for complete genomes and metagenomic mode for fragments
- 📊 **Multiple Output Formats**: GenBank, GFF3, GCA, and SCO formats
- 🎯 **Accurate**: Advanced start codon recognition with Shine-Dalgarno detection
- 💾 **Memory Efficient**: Optimized data structures for large-scale genomic analysis

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
orphos-core = "0.1.0"
```

## Quick Start

### Basic Usage

```rust
use orphos_core::{OrphosAnalyzer, config::OrphosConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create analyzer with default configuration
    let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    
    // Analyze a FASTA file
    let results = analyzer.analyze_file("genome.fasta")?;
    
    println!("Found {} genes", results.genes.len());
    println!("{}", results.output);
    
    Ok(())
}
```

### Analyzing Sequences Directly

```rust
use orphos_core::{OrphosAnalyzer, config::OrphosConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    
    let sequence = "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA...";
    let results = analyzer.analyze_sequence(sequence, Some("MyGenome".to_string()))?;
    
    for gene in &results.genes {
        println!("Gene at {}..{} on {} strand", 
            gene.start, gene.end, 
            if gene.strand == bio::bio_types::strand::Strand::Forward { "+" } else { "-" }
        );
    }
    
    Ok(())
}
```

### Custom Configuration

```rust
use orphos_core::config::{OrphosConfig, OutputFormat};

let config = OrphosConfig {
    closed_ends: true,           // Complete genome (circular)
    mask_n_runs: true,           // Mask stretches of N's
    output_format: OutputFormat::Gff,
    num_threads: Some(4),        // Use 4 threads
    ..Default::default()
};

let mut analyzer = OrphosAnalyzer::new(config);
```

### Metagenomic Mode

For analyzing short contigs or mixed community samples:

```rust
use orphos_core::config::OrphosConfig;

let config = OrphosConfig {
    metagenomic: true,
    ..Default::default()
};

let mut analyzer = OrphosAnalyzer::new(config);
let results = analyzer.analyze_file("metagenome.fasta")?;
```

### Module Organization

- **`config`**: Configuration options and output format settings
- **`engine`**: Main analysis engine with training and prediction logic  
- **`types`**: Core data structures (Gene, Training, error types)
- **`results`**: Gene prediction results and sequence information
- **`sequence`**: Sequence encoding, I/O, and processing utilities
- **`algorithms`**: Gene-finding algorithms including:
  - Dynamic programming for gene prediction
  - Gene optimization and overlap resolution
  - Scoring functions for connections
- **`node`**: Gene node management, creation, and scoring
- **`training`**: Training algorithms for Shine-Dalgarno and non-SD models
- **`output`**: Output formatters for GenBank, GFF, GCA, and SCO
- **`bitmap`**: Efficient sequence encoding utilities
- **`metagenomic`**: Metagenomic mode presets and models

## Output Formats

### GenBank (`.gbk`)

Rich annotation format with full feature information:

```
LOCUS       MyGenome               4641652 bp    DNA     linear   BCT
FEATURES             Location/Qualifiers
     CDS             190..255
                     /gene="1"
                     /protein_id="MyGenome_1"
                     /translation="MTKRSAAAAAAVAAGMTSA"
```

### GFF3 (`.gff`)

Standard genome annotation format:

```
##gff-version 3
MyGenome    Orphos  CDS  190  255  .  +  0  ID=MyGenome_1;
```

### GCA (`.gca`)

Tab-delimited gene coordinate annotation.

### SCO (`.sco`)

Simple coordinate output with minimal information.

## Configuration Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `metagenomic` | `bool` | `false` | Enable metagenomic mode for fragments |
| `closed_ends` | `bool` | `false` | Treat sequences as complete genomes |
| `mask_n_runs` | `bool` | `false` | Mask runs of N characters |
| `force_non_sd` | `bool` | `false` | Disable Shine-Dalgarno detection |
| `quiet` | `bool` | `false` | Suppress informational output |
| `output_format` | `OutputFormat` | `Genbank` | Output format selection |
| `translation_table` | `Option<u8>` | `None` | NCBI genetic code table (1-25) |
| `num_threads` | `Option<usize>` | `None` | Number of parallel threads |

## Error Handling

All operations return `Result<T, OrphosError>` with detailed error types:

```rust
use orphos_core::types::OrphosError;

match analyzer.analyze_file("genome.fasta") {
    Ok(results) => println!("Success: {} genes", results.genes.len()),
    Err(OrphosError::SequenceTooShort { length, min }) => {
        eprintln!("Sequence too short: {} bp (minimum: {} bp)", length, min);
    }
    Err(OrphosError::IoError(e)) => {
        eprintln!("I/O error: {}", e);
    }
    Err(e) => eprintln!("Error: {}", e),
}
```

## Contributing

Contributions are welcome! Please see the main [Orphos repository](https://github.com/FullHuman/orphos) for contribution guidelines.

## License

This project is licensed under the GNU General Public License v3.0 or later - see the [LICENSE](../LICENSE) file for details.

## Citation

If you use Orphos in your research, please cite:

```bibtex
@software{orphos,
  title = {Orphos: High-Performance Prokaryotic Gene Prediction},
  author = {Floriel Fedry},
  year = {2025},
  url = {https://github.com/FullHuman/orphos}
}
```

## Related Projects

- **[orphos-cli](../orphos-cli)**: Command-line interface for gene prediction
- **[orphos-python](../orphos-python)**: Python bindings via PyO3
- **[orphos-wasm](../orphos-wasm)**: WebAssembly module for browser/Node.js

## Acknowledgments

This implementation is based on Prodigal, originally developed by Doug Hyatt. Orphos provides a modern, type-safe Rust implementation while maintaining compatibility with the original algorithms.

## Support

- 📖 [Documentation](https://docs.rs/orphos-core)
- 🐛 [Issue Tracker](https://github.com/FullHuman/orphos/issues)
- 💬 [Discussions](https://github.com/FullHuman/orphos/discussions)
