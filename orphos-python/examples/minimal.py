"""
Minimal example showing how to use Orphos Python bindings.
This is the simplest possible usage - just 5 lines!
"""

import prodigal

# Analyze a genome file
result = prodigal.analyze_file("genome.fasta")

# Print results
print(f"Found {result.gene_count} genes in {result.sequence_count} sequence(s)")

# Output is in GenBank format by default
print(result.output[:200])  # Show first 200 characters
