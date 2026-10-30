# Orphos

[![CI](https://github.com/FullHuman/orphos/workflows/CI/badge.svg)](https://github.com/FullHuman/orphos/actions)
[![Coverage](https://github.com/FullHuman/orphos/workflows/Code%20Coverage/badge.svg)](https://github.com/FullHuman/orphos/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A fast, parallel Rust implementation of Prodigal, a tool for finding protein-coding genes in microbial genomes.

## What is Orphos?

Orphos is a high-performance reimplementation of [Prodigal](https://github.com/hyattpd/Prodigal), the widely-used prokaryotic gene prediction tool. Written in Rust, Orphos delivers the same accurate gene finding algorithm with improved performance and modern language features.

### Why Choose Orphos Over Prodigal?

- **Faster Performance**: 6x to 18x faster than Prodigal on complete bacterial genomes, with the same predictions (see [Benchmarks](#-benchmarks))
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

## 📊 Benchmarks

Wall-clock time to predict the genes of six complete bacterial genomes, written as GFF. Each number is the mean of five runs after one warm-up run. Orphos was timed on one thread (`RAYON_NUM_THREADS=1`) and with its default thread pool, which uses every core. The speed-up compares Prodigal 2.6.3 with Orphos on the default thread pool.

### Single-genome mode

| Genome | Size | Genes | Prodigal 2.6.3 | Orphos, 1 thread | Orphos, 10 threads | Speed-up | Peak memory, Prodigal / Orphos |
|---|---:|---:|---:|---:|---:|---:|---:|
| *Mycoplasmoides genitalium* G37 (NC_000908.2) | 0.58 Mb | 995 | 0.64 s | 0.11 s | 0.11 s | 6.0x | 82 / 21 MB |
| *Helicobacter pylori* 26695 (NC_000915.1) | 1.67 Mb | 1,579 | 1.44 s | 0.10 s | 0.08 s | 18.3x | 93 / 57 MB |
| *Escherichia coli* K-12 MG1655 (NC_000913.3) | 4.64 Mb | 4,319 | 4.55 s | 0.40 s | 0.30 s | 15.3x | 170 / 198 MB |
| *Mycobacterium tuberculosis* H37Rv (NC_000962.3) | 4.41 Mb | 4,085 | 7.13 s | 0.53 s | 0.44 s | 16.3x | 177 / 236 MB |
| *Pseudomonas aeruginosa* PAO1 (NC_002516.2) | 6.26 Mb | 5,681 | 7.22 s | 0.55 s | 0.43 s | 16.7x | 215 / 253 MB |
| *Streptomyces coelicolor* A3(2) (NC_003888.3) | 8.67 Mb | 7,738 | 9.66 s | 0.78 s | 0.59 s | 16.5x | 276 / 347 MB |

Most of the gain is single-threaded: on one thread Orphos is 6x to 14x faster than Prodigal, and the default thread pool adds up to another 1.3x on the larger genomes.

### Metagenomic mode (`-p meta`)

| Genome | Size | Genes | Prodigal 2.6.3 | Orphos, 1 thread | Orphos, 10 threads | Speed-up |
|---|---:|---:|---:|---:|---:|---:|
| *Mycoplasmoides genitalium* G37 (NC_000908.2) | 0.58 Mb | 494 | 1.03 s | 0.08 s | 0.08 s | 13.6x |
| *Escherichia coli* K-12 MG1655 (NC_000913.3) | 4.64 Mb | 4,319 | 32.86 s | 2.54 s | 2.15 s | 15.3x |

### Same predictions

In single-genome mode the GFF output of Orphos is identical to Prodigal's on all six genomes, down to every score, once the tool name in the source column is ignored. In metagenomic mode it is identical to a build of Prodigal's GitHub master. Prodigal 2.6.3 predicts three genes fewer than master on each of the two genomes because of a short-gene penalty bug fixed upstream after the 2.6.3 release ([hyattpd/Prodigal@b34b63b](https://github.com/hyattpd/Prodigal/commit/b34b63b)).

### Setup

- Apple M1 Max (8 performance and 2 efficiency cores), 32 GB of memory, macOS 27.0.1.
- Orphos 0.4.0 (commit 357cb49) built with `cargo build --release` on Rust 1.94.0.
- Prodigal 2.6.3 from Homebrew. An `-O3` build of Prodigal's GitHub master (commit c1e2d36) ran within 5% of 2.6.3 on every genome.
- Each tool was run as `<tool> -q -p <mode> -i genome.fasta -f gff -o out.gff` and timed with `/usr/bin/time -l`, which also reports the peak memory. The standard deviation over the five runs was below 6% in every case.
- The genomes are the RefSeq sequences in `orphos-cli/tests/data/single/`, so the table can be reproduced from a checkout.

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
orphos-core = "0.4.0"
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
    circular=False,        # Set True for circular topology
    translation_table=11   # Use translation table 11
)
result = orphos.analyze_file("genome.fasta", options)
```

For circular chromosomes/plasmids, enable wraparound detection with `--circular` in CLI or `circular=True` in Python options.

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

Orphos predicts the genes of a complete bacterial genome 6x to 18x faster than Prodigal and produces the same output. Most of the gain is single-threaded. On top of that, Rayon runs the per-strand passes, the k-mer counting and the start-site scoring on the thread pool, and analyses the records of a multi-record FASTA file in parallel. Set `RAYON_NUM_THREADS` to limit the pool. See [Benchmarks](#-benchmarks) for timings on six genomes.

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

This project is based on the original [Prodigal](https://github.com/hyattpd/Prodigal) by Doug Hyatt ([license](https://github.com/hyattpd/Prodigal/blob/GoogleImport/LICENSE)).
