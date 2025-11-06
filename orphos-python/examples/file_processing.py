#!/usr/bin/env python3
"""
Example showing how to process FASTA files with Orphos Python bindings.
"""

import prodigal
import os
import sys

def process_single_file(input_path, output_path=None, options=None):
    """Process a single FASTA file."""
    print(f"Processing: {input_path}")
    
    try:
        result = prodigal.analyze_file(input_path, options)
        print(f"  Genes found: {result.gene_count}")
        print(f"  Sequences: {result.sequence_count}")
        
        if output_path:
            with open(output_path, 'w') as f:
                f.write(result.output)
            print(f"  Output saved to: {output_path}")
        
        return result
    except Exception as e:
        print(f"  Error: {e}")
        return None

def process_directory(input_dir, output_dir, file_extension=".fasta", options=None):
    """Process all FASTA files in a directory."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    files = [f for f in os.listdir(input_dir) if f.endswith(file_extension)]
    
    if not files:
        print(f"No files with extension '{file_extension}' found in {input_dir}")
        return
    
    print(f"Found {len(files)} file(s) to process\n")
    
    total_genes = 0
    for filename in files:
        input_path = os.path.join(input_dir, filename)
        output_filename = filename.replace(file_extension, ".gff")
        output_path = os.path.join(output_dir, output_filename)
        
        result = process_single_file(input_path, output_path, options)
        if result:
            total_genes += result.gene_count
        print()
    
    print(f"Total genes found: {total_genes}")

def main():
    print("=== Orphos File Processing Example ===\n")
    
    # This example assumes you have test data in ../orphos-cli/tests/data/
    # Adjust the paths as needed for your setup
    test_data_dir = "../orphos-cli/tests/data"
    
    if not os.path.exists(test_data_dir):
        print(f"Test data directory not found: {test_data_dir}")
        print("This is just an example script.")
        print("\nTo use this script with your own data:")
        print("  python file_processing.py <input_file> [output_file]")
        return
    
    # Example 1: Process a single file with different formats
    ecoli_fasta = os.path.join(test_data_dir, "ecoli.fasta")
    
    if os.path.exists(ecoli_fasta):
        print("Example 1: Processing E. coli genome with different formats\n")
        
        # GenBank format
        process_single_file(
            ecoli_fasta,
            "ecoli_output.gbk",
            prodigal.OrphosOptions(format="gbk")
        )
        print()
        
        # GFF format
        process_single_file(
            ecoli_fasta,
            "ecoli_output.gff",
            prodigal.OrphosOptions(format="gff")
        )
        print()
        
        # Simple coordinates format
        process_single_file(
            ecoli_fasta,
            "ecoli_output.sco",
            prodigal.OrphosOptions(format="sco")
        )
        print()
    
    # Example 2: Batch process files
    single_dir = os.path.join(test_data_dir, "single")
    if os.path.exists(single_dir):
        print("\nExample 2: Batch processing genomes\n")
        process_directory(
            single_dir,
            "batch_output",
            file_extension=".fna",
            options=prodigal.OrphosOptions(format="gff")
        )

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Command-line usage
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else None
        
        options = prodigal.OrphosOptions(format="gff")
        process_single_file(input_file, output_file, options)
    else:
        # Run examples
        main()
