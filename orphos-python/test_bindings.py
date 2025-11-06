#!/usr/bin/env python3
"""
Quick test to verify the prodigal Python bindings work correctly.
"""

import sys

def test_import():
    """Test that the module can be imported."""
    try:
        import prodigal
        print("✓ Successfully imported prodigal module")
        print(f"  Version: {prodigal.__version__}")
        return True
    except ImportError as e:
        print(f"✗ Failed to import prodigal: {e}")
        return False

def test_basic_analysis():
    """Test basic sequence analysis."""
    try:
        import prodigal
        
        # Longer test sequence (needs to be >= 96 bp)
        fasta = """>test_sequence
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA
ACGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAA
"""
        
        result = prodigal.analyze_sequence(fasta)
        print(f"✓ Basic analysis works")
        print(f"  Genes found: {result.gene_count}")
        print(f"  Sequences: {result.sequence_count}")
        print(f"  Output length: {len(result.output)} chars")
        return True
    except Exception as e:
        print(f"✗ Basic analysis failed: {e}")
        return False

def test_options():
    """Test OrphosOptions."""
    try:
        import prodigal
        
        options = prodigal.OrphosOptions(
            mode="meta",
            format="gff",
            closed_ends=True
        )
        print(f"✓ OrphosOptions works")
        print(f"  {options}")
        return True
    except Exception as e:
        print(f"✗ OrphosOptions failed: {e}")
        return False

def test_metagenomic_mode():
    """Test metagenomic mode."""
    try:
        import prodigal
        
        fasta = """>test
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA
ACGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAA
"""
        
        options = prodigal.OrphosOptions(mode="meta", format="gff")
        result = prodigal.analyze_sequence(fasta, options)
        print(f"✓ Metagenomic mode works")
        print(f"  Genes found: {result.gene_count}")
        return True
    except Exception as e:
        print(f"✗ Metagenomic mode failed: {e}")
        return False

def main():
    print("=== Orphos Python Bindings Test Suite ===")
    print()
    
    tests = [
        ("Module import", test_import),
        ("Basic analysis", test_basic_analysis),
        ("OrphosOptions", test_options),
        ("Metagenomic mode", test_metagenomic_mode),
    ]
    
    results = []
    for name, test_func in tests:
        print(f"\nTest: {name}")
        print("-" * 40)
        results.append(test_func())
        print()
    
    print("=" * 40)
    passed = sum(results)
    total = len(results)
    print(f"\nResults: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n✓ All tests passed!")
        return 0
    else:
        print(f"\n✗ {total - passed} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
