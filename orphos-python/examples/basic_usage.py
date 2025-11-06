#!/usr/bin/env python3
"""
Basic usage example for Orphos Python bindings.
"""

import prodigal

def main():
    # Example FASTA sequence (E. coli sequence)
    fasta_sequence = """>NC_000913.3 Escherichia coli str. K-12 substr. MG1655
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC
ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG
CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA
GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC
"""

    print("=== Orphos Python Bindings Example ===\n")
    
    # Example 1: Analyze with default settings (single genome, GenBank format)
    print("1. Analyzing with default settings...")
    result = prodigal.analyze_sequence(fasta_sequence)
    print(f"   Genes found: {result.gene_count}")
    print(f"   Sequences analyzed: {result.sequence_count}")
    print(f"   Output length: {len(result.output)} characters\n")
    
    # Example 2: Use metagenomic mode
    print("2. Analyzing in metagenomic mode...")
    meta_options = prodigal.OrphosOptions(mode="meta")
    result = prodigal.analyze_sequence(fasta_sequence, meta_options)
    print(f"   Genes found: {result.gene_count}\n")
    
    # Example 3: GFF output format
    print("3. Generating GFF output...")
    gff_options = prodigal.OrphosOptions(format="gff")
    result = prodigal.analyze_sequence(fasta_sequence, gff_options)
    print(f"   First 200 characters of GFF output:")
    print(f"   {result.output[:200]}...\n")
    
    # Example 4: Custom options
    print("4. Using custom options...")
    custom_options = prodigal.OrphosOptions(
        mode="single",
        format="sco",
        closed_ends=True,
        translation_table=11
    )
    result = prodigal.analyze_sequence(fasta_sequence, custom_options)
    print(f"   Genes found: {result.gene_count}")
    print(f"   Format: Simple coordinates")
    print(f"   Output preview:")
    print(f"   {result.output[:200]}...\n")
    
    # Example 5: Show OrphosOptions representation
    print("5. OrphosOptions object:")
    print(f"   {custom_options}\n")
    
    print("=== Examples completed successfully! ===")

if __name__ == "__main__":
    main()
