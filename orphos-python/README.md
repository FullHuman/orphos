# Orphos Python Bindings

Python bindings for Orphos, a fast and accurate gene prediction tool for microbial genomes.

## Installation

### From Source

You'll need Rust and Python 3.8+ installed. Then:

```bash
# Install maturin (the build tool for Python+Rust projects)
pip install maturin

# Build and install in development mode
cd orphos-python
maturin develop --release

# Or build a wheel
maturin build --release
pip install target/wheels/orphos-*.whl
```

## Quick Start

```python
import orphos

# Analyze a FASTA file
result = orphos.analyze_file("genome.fasta")
print(f"Found {result.gene_count} genes in {result.sequence_count} sequences")
print(result.output)  # GenBank formatted output

# Analyze a sequence string
fasta_string = """>seq1
ATGCGATCGATCGATCGATCGATCGATCG...
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

## API Reference

### `analyze_sequence(fasta_content, options=None)`

Analyze DNA sequences from a FASTA-formatted string.

**Parameters:**
- `fasta_content` (str): FASTA-formatted sequence(s)
- `options` (OrphosOptions, optional): Configuration options

**Returns:** `OrphosResult`

### `analyze_file(file_path, options=None)`

Analyze DNA sequences from a FASTA file.

**Parameters:**
- `file_path` (str): Path to the FASTA file
- `options` (OrphosOptions, optional): Configuration options

**Returns:** `OrphosResult`

### `OrphosOptions`

Configuration options for gene prediction.

**Attributes:**
- `mode` (str): "single" for single genome mode, "meta" for metagenomic mode (default: "single")
- `format` (str): Output format - "gbk", "gff", "sco", or "gca" (default: "gbk")
- `closed_ends` (bool): Don't allow genes to run off edges (default: False)
- `mask_n_runs` (bool): Mask runs of N's in the sequence (default: False)
- `force_non_sd` (bool): Force non-Shine-Dalgarno model (default: False)
- `translation_table` (int, optional): Translation table 1-25 (excluding 7, 8, 17-20)
- `num_threads` (int, optional): Number of threads to use
- `quiet` (bool): Suppress informational output (default: True)

### `OrphosResult`

Result from gene prediction.

**Attributes:**
- `output` (str): Formatted output (GenBank, GFF, etc.)
- `gene_count` (int): Total number of genes predicted
- `sequence_count` (int): Number of sequences analyzed

## Examples

### Single Genome Mode (Default)

```python
import orphos

# Analyze with default settings
result = orphos.analyze_file("ecoli.fasta")
print(f"Found {result.gene_count} genes")

# Save GenBank output
with open("output.gbk", "w") as f:
    f.write(result.output)
```

### Metagenomic Mode

```python
import orphos

options = orphos.OrphosOptions(mode="meta")
result = orphos.analyze_file("metagenome.fasta", options)
```

### GFF Output Format

```python
import orphos

options = orphos.OrphosOptions(format="gff")
result = orphos.analyze_file("genome.fasta", options)

# Save GFF output
with open("output.gff", "w") as f:
    f.write(result.output)
```

### Custom Translation Table

```python
import orphos

# Use translation table 4 (Mycoplasma/Spiroplasma)
options = orphos.OrphosOptions(translation_table=4)
result = orphos.analyze_file("mycoplasma.fasta", options)
```

### Processing Multiple Files

```python
import orphos
import os

fasta_dir = "genomes"
output_dir = "predictions"
os.makedirs(output_dir, exist_ok=True)

options = orphos.OrphosOptions(format="gff")

for filename in os.listdir(fasta_dir):
    if filename.endswith(".fasta"):
        input_path = os.path.join(fasta_dir, filename)
        output_path = os.path.join(output_dir, filename.replace(".fasta", ".gff"))
        
        result = orphos.analyze_file(input_path, options)
        
        with open(output_path, "w") as f:
            f.write(result.output)
        
        print(f"{filename}: {result.gene_count} genes")
```

## Output Formats

### GenBank (gbk)
Standard GenBank feature table format with gene coordinates and annotations.

### GFF (gff)
GFF3 format with gene features and attributes.

### Simple Coordinates (sco)
Simple tab-delimited format with gene coordinates.

### Gene Calls (gca)
Detailed gene call information including scores and training data.

## Performance

The Python bindings have minimal overhead compared to the native Rust implementation. Large genomes and metagenomes can be processed efficiently.

## License

GPL-3.0-or-later

## Citation

If you use Orphos in your research, please cite:

Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Orphos: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119.
