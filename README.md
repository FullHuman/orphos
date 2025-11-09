# Orphos

[![CI](https://github.com/FullHuman/orphos/workflows/CI/badge.svg)](https://github.com/FullHuman/orphos/actions)
[![Coverage](https://github.com/FullHuman/orphos/workflows/Code%20Coverage/badge.svg)](https://github.com/FullHuman/orphos/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A fast, parallel Rust implementation of Prodigal, a tool for finding protein-coding genes in microbial genomes.

## What is Orphos?

Orphos is a high-performance reimplementation of [Prodigal](https://github.com/hyattpd/Prodigal), the widely-used prokaryotic gene prediction tool. Written in Rust, Orphos delivers the same accurate gene finding algorithm with improved performance and modern language features.

### Why Choose Orphos Over Prodigal?

- **Faster Performance**: Multi-threaded processing using Rayon for parallel genome analysis
- **Memory Efficient**: Optimized memory usage for handling large genomes and metagenomes
- **Browser Support**: Unique WebAssembly build - runs in web browsers
- **100% Compatible**: Output formats fully compatible with original Prodigal (GFF3, GenBank, etc.)
- **Modern Codebase**: Written in safe Rust with excellent error handling
- **Multiple Interfaces**: CLI, Rust library, Python bindings, and WebAssembly
- **Easy Installation**: Available via Cargo, pip, Homebrew, and Conda

## Components

Orphos is available in multiple forms:

- **`orphos-cli`**: Command-line interface for gene prediction
- **`orphos-core`**: Rust library for integrating into your own projects
- **`orphos-python`**: Python bindings (via PyO3)
- **`orphos-wasm`**: WebAssembly module for browser/Node.js usage

## Features

- **High Performance**: Multi-threaded processing using Rayon
- **Memory Efficient**: Optimized memory usage for large genomes
- **Compatible**: Output format compatible with original Prodigal
- **Cross-Platform**: Works on Linux, macOS, and Windows

## 📦 Installation

### Homebrew (macOS/Linux)
```bash
brew install FullHuman/tap/orphos
```

### Using Cargo
```bash
cargo install orphos-cli
```

### From Source
```bash
git clone https://github.com/FullHuman/orphos.git
cd orphos
cargo install --path orphos-cli
```

### Python Bindings
```bash
pip install orphos
```

### Rust Library
Add to your `Cargo.toml`:
```toml
[dependencies]
orphos-core = "0.1.0"
```

## 🏃 Quick Start

### Command Line (orphos-cli)

```bash
# Basic gene prediction
orphos -i input.fasta -o output.gbk

# Output in GFF format
orphos -i input.fasta -f gff -o output.gff

# Use custom training file
orphos -i input.fasta -t training.trn -o output.gbk

# Metagenomic mode
orphos -i metagenome.fasta -p meta -o output.gff
```

### Python

```python
import orphos

# Analyze a FASTA file
result = orphos.analyze_file("genome.fasta")
print(f"Found {result.gene_count} genes")
print(result.output)  # GenBank formatted output

# Analyze a sequence string
fasta_string = """>seq1
ATGCGATCGATCGATCGATCG...
"""
result = orphos.analyze_sequence(fasta_string)

# Customize options
options = orphos.OrphosOptions(
    mode="meta",           # Use metagenomic mode
    format="gff",          # Output in GFF format
    closed_ends=True,      # Don't allow genes off edges
    translation_table=11   # Use translation table 11
)
result = orphos.analyze_file("genome.fasta", options)
```

### Rust Library (orphos-core)

```rust
use orphos_core::{OrphosAnalyzer, config::OrphosConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create analyzer with default configuration
    let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    
    // Analyze a genome sequence
    let results = analyzer.analyze_sequence(
        "ATGCGATCGATCG...",
        Some("MyGenome".to_string())
    )?;
    
    println!("Found {} genes", results.genes.len());
    Ok(())
}
```

For more advanced usage with type-safe training:

```rust
use orphos_core::engine::{UntrainedOrphos, Orphos, Untrained};
use orphos_core::config::OrphosConfig;
use orphos_core::sequence::encoded::EncodedSequence;

// Create an untrained analyzer
let mut untrained = UntrainedOrphos::with_config(OrphosConfig::default())?;

// Encode the sequence
let encoded = EncodedSequence::without_masking(b"ATGCGATCGATCG...");

// Train on the genome (type changes to TrainedOrphos)
let trained = untrained.train_single_genome(&encoded)?;
```

## 📖 Documentation

// TODO: Add documentation links

## 🔧 Advanced Features

### Output Formats

Orphos supports multiple output formats:
- **GenBank (GBK)**: Rich feature annotation format (default)
- **GFF3**: General Feature Format version 3
- **GCA**: Gene coordinate annotation  
- **SCO**: Simple coordinate output

### Modes

- **Single Genome Mode**: Train on a complete genome for optimal gene prediction (default)
- **Metagenomic Mode**: Predict genes in fragmented or mixed sequences

### Performance

- **Parallel Processing**: Multi-threaded execution using Rayon
- **Memory Efficient**: Optimized memory usage for large genomes
- **High Performance**: Significantly faster than the original C implementation

## 🧪 Testing

```bash
# Run all tests
cargo test

# Run with coverage
cargo install cargo-tarpaulin
cargo cov-fast

# Run benchmarks
cargo bench
```

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## 📄 License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

This project is inspired by the original [Prodigal](https://github.com/hyattpd/Prodigal) by Doug Hyatt.