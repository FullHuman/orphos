# Orphos CLI

[![CI](https://github.com/FullHuman/orphos/workflows/CI/badge.svg)](https://github.com/FullHuman/orphos/actions)
[![Coverage](https://github.com/FullHuman/orphos/workflows/Code%20Coverage/badge.svg)](https://github.com/FullHuman/orphos/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Crates.io](https://img.shields.io/crates/v/orphos-cli.svg)](https://crates.io/crates/orphos-cli)

Command-line interface for Orphos, a fast, parallel Rust implementation of Prodigal for finding protein-coding genes in microbial genomes.

## Features

- 🚀 **High Performance**: Multi-threaded processing using Rayon
- 💾 **Memory Efficient**: Optimized for large genomes and metagenomic assemblies
- 🔄 **Compatible**: Output format compatible with original Prodigal
- 🌍 **Cross-Platform**: Works on Linux, macOS, and Windows
- 📊 **Multiple Output Formats**: GenBank, GFF3, SCO, and GCA formats
- 🧬 **Flexible Modes**: Single genome and metagenomic analysis modes

## Installation

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

### Homebrew (macOS/Linux)

```bash
brew tap FullHuman/orphos
brew install orphos
```

### Conda

```bash
conda install -c bioconda orphos
```

## Quick Start

### Basic Usage

```bash
# Analyze a genome and output GenBank format
orphos -i genome.fasta -o genes.gbk

# Analyze with GFF3 output
orphos -i genome.fasta -f gff -o genes.gff

# Metagenomic mode for short contigs
orphos -i metagenome.fasta -p meta -o genes.gff

# Complete circular genome (closed ends)
orphos -i plasmid.fasta -c -o plasmid.gbk
```

### Reading from stdin/stdout

```bash
# Input from stdin
cat genome.fasta | orphos -o genes.gbk

# Output to stdout
orphos -i genome.fasta > genes.gbk

# Pipe both
cat genome.fasta | orphos > genes.gbk
```

## Command-Line Options

### Required/Input

| Option | Short | Long | Description |
|--------|-------|------|-------------|
| Input file | `-i` | `--input` | Input FASTA file (default: stdin) |
| Output file | `-o` | `--output` | Output file (default: stdout) |

### Output Options

| Option | Short | Long | Default | Description |
|--------|-------|------|---------|-------------|
| Format | `-f` | `--format` | `gbk` | Output format: `gbk`, `gff`, `sco`, `gca` |

### Analysis Options

| Option | Short | Long | Default | Description |
|--------|-------|------|---------|-------------|
| Mode | `-p` | `--mode` | `single` | Analysis mode: `single` or `meta` |
| Closed ends | `-c` | `--closed` | false | No genes off edges (for complete genomes) |
| Mask N's | `-m` | `--mask` | false | Mask runs of N's |
| Translation table | `-g` | `--translation-table` | auto | Translation table (1-25) |
| Training file | `-t` | `--training` | - | Use pre-trained parameters |

### Other Options

| Option | Short | Long | Description |
|--------|-------|------|-------------|
| Quiet | `-q` | `--quiet` | Suppress progress messages |
| Help | `-h` | `--help` | Display help information |
| Version | `-V` | `--version` | Display version information |

## Output Formats

### GenBank (gbk)

Rich annotation format with gene features, translations, and metadata.

```bash
orphos -i genome.fasta -f gbk -o genes.gbk
```

### GFF3 (gff)

General Feature Format version 3, widely used in genomics pipelines.

```bash
orphos -i genome.fasta -f gff -o genes.gff
```

### Simple Coordinate Output (sco)

Tab-delimited gene coordinates for easy parsing.

```bash
orphos -i genome.fasta -f sco -o genes.sco
```

### Gene Coordinate Annotation (gca)

Compact coordinate format.

```bash
orphos -i genome.fasta -f gca -o genes.gca
```

## Analysis Modes

### Single Genome Mode (default)

Use for complete or near-complete genomes (>100kb). Orphos will train on the genome to optimize gene prediction accuracy.

```bash
orphos -i complete_genome.fasta -o genes.gbk
```

**Best for:**
- Complete bacterial genomes
- Complete archaeal genomes
- Large contigs or chromosomes
- Closed genomes

### Metagenomic Mode

Use for short contigs or mixed metagenomic assemblies. Uses pre-trained parameters instead of training on the input.

```bash
orphos -i metagenome_contigs.fasta -p meta -o genes.gff
```

**Best for:**
- Metagenomic assemblies
- Short contigs (<100kb)
- Mixed-species samples
- Fragmented sequences

## Advanced Examples

### Complete Circular Genome

For complete circular genomes (chromosomes, plasmids), use the `-c` flag to prevent genes from being called off the edges:

```bash
orphos -i circular_plasmid.fasta -c -o plasmid.gbk
```

### Custom Translation Table

Specify a custom genetic code (translation table):

```bash
# Use translation table 4 (Mycoplasma/Spiroplasma)
orphos -i mycoplasma.fasta -g 4 -o genes.gbk

# Use translation table 11 (Bacterial and Archaea)
orphos -i bacteria.fasta -g 11 -o genes.gbk
```

### Masking Low-Quality Regions

Mask runs of N's in low-quality sequences:

```bash
orphos -i draft_assembly.fasta -m -o genes.gff
```

### Batch Processing

Process multiple genomes:

```bash
for genome in genomes/*.fasta; do
    base=$(basename "$genome" .fasta)
    orphos -i "$genome" -f gff -o "results/${base}.gff"
done
```

### Pipeline Integration

Integrate with other bioinformatics tools:

```bash
# Find genes and extract protein sequences
orphos -i genome.fasta -f gff -o genes.gff
# ... then use genes.gff with other tools

# Combine with annotation pipelines
orphos -i assembly.fasta -p meta -f gff -o genes.gff
prokka --proteins genes.gff --outdir annotation genome.fasta
```

## Performance Tips

1. **Use multiple cores**: Orphos automatically uses all available CPU cores via Rayon
2. **Metagenomic mode for many small contigs**: Faster than single mode for fragmented assemblies
3. **Batch processing**: Process multiple files in parallel using shell scripting
4. **Large files**: Orphos handles multi-GB files efficiently

## Translation Tables

Orphos supports NCBI translation tables 1-25 (excluding 7, 8, 17-20). Common tables:

| Table | Name | Organisms |
|-------|------|-----------|
| 1 | Standard | Most eukaryotes |
| 4 | Mycoplasma/Spiroplasma | Mycoplasma, Spiroplasma |
| 11 | Bacterial, Archaeal, Plant Plastid | Most bacteria and archaea (default) |
| 25 | Candidate Division SR1, Gracilibacteria | Certain bacteria |


## Related Projects

- **[orphos-core](../orphos-core)**: Rust library for gene prediction
- **[orphos-python](../orphos-python)**: Python bindings
- **[orphos-wasm](../orphos-wasm)**: WebAssembly module for browser/Node.js

## Contributing

We welcome contributions! Please see the main [repository](https://github.com/FullHuman/orphos) for contribution guidelines.

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](../LICENSE) file for details.

## Citation

If you use Orphos in your research, please cite:

```bibtex
# TODO: Add citation information
```

## Acknowledgments

This project is inspired by the original [Prodigal](https://github.com/hyattpd/Prodigal) by Doug Hyatt. We thank the authors for their groundbreaking work in prokaryotic gene prediction.

## Support

- **Issues**: [GitHub Issues](https://github.com/FullHuman/orphos/issues)
- **Discussions**: [GitHub Discussions](https://github.com/FullHuman/orphos/discussions)
- **Documentation**: [docs.rs/orphos-cli](https://docs.rs/orphos-cli)
