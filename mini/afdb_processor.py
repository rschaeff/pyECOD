#!/usr/bin/env python3
"""
AFDB Processor - Fixed version with correct paths

Run from /home/rschaeff/dev/pyecod/mini directory
"""

import sys
import os
from pathlib import Path

# Add the mini directory to path
mini_dir = Path("/home/rschaeff/dev/pyecod/mini")
afdb_dir = Path("/home/rschaeff/work/afdb_pyecod_mini_bench")
sys.path.insert(0, str(mini_dir))
sys.path.insert(0, str(afdb_dir))

from afdb_pyecod_mini import partition_afdb_protein, process_afdb_batch, AFDBConfig

def create_fixed_config():
    """Create AFDB config with correct paths"""
    config = AFDBConfig()
    
    # Fix paths to use the actual location
    config.test_data_dir = mini_dir / "test_data"
    config.domain_lengths_file = config.test_data_dir / "domain_lengths.csv"
    config.protein_lengths_file = config.test_data_dir / "protein_lengths.csv"
    config.domain_definitions_file = config.test_data_dir / "domain_definitions.csv"
    config.reference_blacklist_file = config.test_data_dir / "reference_blacklist.csv"
    config.output_dir = Path("/tmp/afdb_domains")
    
    # Ensure output directory exists
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    return config

def test_single_protein(protein_id: str):
    """Test processing a single AFDB protein"""
    print(f"Testing protein: {protein_id}")
    print("=" * 50)
    
    config = create_fixed_config()
    result = partition_afdb_protein(protein_id, config, verbose=True)
    
    if result.success:
        print(f"✅ SUCCESS: {result.domain_count} domains found")
        if result.domain_count > 0:
            print(f"   Coverage: {result.coverage_percent:.1f}%")
            print(f"   Evidence: {result.evidence_counts}")
            
            # Show TSV output
            tsv_file = config.output_dir / f"{protein_id}_mini.domains.tsv"
            if tsv_file.exists():
                print(f"   Output: {tsv_file}")
                print("   First few lines:")
                with open(tsv_file, 'r') as f:
                    for i, line in enumerate(f):
                        if i < 5:  # Show first 5 lines
                            print(f"     {line.strip()}")
        else:
            print("   No domains assigned (valid result for some proteins)")
    else:
        print(f"❌ FAILED: {result.error}")
    
    return result.success

def test_small_batch(batch_name: str = "afdb_human_pyecod_batch_001", max_proteins: int = 10):
    """Test processing a small batch"""
    print(f"Testing batch: {batch_name} (first {max_proteins} proteins)")
    print("=" * 70)
    
    config = create_fixed_config()
    result = process_afdb_batch(batch_name, config, max_proteins=max_proteins, verbose=True)
    
    print(f"\nBatch Results:")
    print(f"  Success: {result['success_count']}/{result['total_proteins']}")
    print(f"  Failed: {result['failed_count']}")
    print(f"  Total domains: {result['total_domains']}")
    print(f"  Avg domains/protein: {result['avg_domains']:.2f}")
    print(f"  Output file: {result['output_file']}")
    
    return result['success_count'] > 0

def validate_setup():
    """Validate that everything is ready for AFDB processing"""
    print("Validating AFDB setup...")
    print("=" * 30)
    
    config = create_fixed_config()
    
    # Check reference files
    required_files = [
        ("Domain lengths", config.domain_lengths_file),
        ("Protein lengths", config.protein_lengths_file),
        ("Domain definitions", config.domain_definitions_file)
    ]
    
    all_good = True
    for name, file_path in required_files:
        if file_path.exists():
            size_mb = file_path.stat().st_size / (1024 * 1024)
            print(f"  ✅ {name}: {file_path} ({size_mb:.1f} MB)")
        else:
            print(f"  ❌ {name}: {file_path} (missing)")
            all_good = False
    
    # Check optional files
    if config.reference_blacklist_file.exists():
        print(f"  ✅ Blacklist: {config.reference_blacklist_file}")
    else:
        print(f"  ⚠️  Blacklist: {config.reference_blacklist_file} (optional)")
    
    # Check AFDB data
    afdb_base = Path("/data/ecod/pyecod_mini_bench/batches")
    if afdb_base.exists():
        batch_dirs = [d for d in afdb_base.iterdir() 
                     if d.is_dir() and d.name.startswith("afdb_human_pyecod_batch_")]
        print(f"  ✅ AFDB batches: {len(batch_dirs)} found")
        for batch_dir in batch_dirs[:3]:  # Show first 3
            summaries_dir = batch_dir / "summaries"
            if summaries_dir.exists():
                protein_count = len(list(summaries_dir.glob("*.develop291.domain_summary.xml")))
                print(f"    - {batch_dir.name}: {protein_count:,} proteins")
    else:
        print(f"  ❌ AFDB batches: {afdb_base} (missing)")
        all_good = False
    
    return all_good

def main():
    """Main test runner"""
    print("AFDB PyECOD Mini - Quick Test")
    print("=" * 60)
    
    # Step 1: Validate setup
    if not validate_setup():
        print("❌ Setup validation failed - fix issues before proceeding")
        return 1
    
    print("\n" + "=" * 60)
    
    # Step 2: Test single proteins
    test_proteins = [
        "P19784_F1",   # High evidence (10,000 items)
        "P0CZ25_F1",   # No evidence (0 items)
        "A8MT69_F1",   # Medium evidence (127 items)
    ]
    
    success_count = 0
    for protein_id in test_proteins:
        if test_single_protein(protein_id):
            success_count += 1
        print()
    
    print(f"Single protein tests: {success_count}/{len(test_proteins)} successful")
    
    # Step 3: Test small batch
    print("\n" + "=" * 60)
    if test_small_batch():
        print("✅ Small batch test successful")
    else:
        print("❌ Small batch test failed")
    
    # Step 4: Show what to do next
    print("\n" + "=" * 60)
    print("Next Steps:")
    if success_count == len(test_proteins):
        print("🎉 All tests passed! Ready for production processing.")
        print()
        print("To process full batches:")
        print("  python afdb_processor_fixed.py --batch afdb_human_pyecod_batch_001")
        print()
        print("To process all batches:")
        print("  for batch in afdb_human_pyecod_batch_{001..005}; do")
        print("    python afdb_processor_fixed.py --batch $batch")
        print("  done")
    else:
        print("❌ Some tests failed - investigate issues before full processing")
    
    return 0 if success_count == len(test_proteins) else 1

def run_full_batch(batch_name: str):
    """Run a full batch"""
    print(f"Processing full batch: {batch_name}")
    config = create_fixed_config()
    result = process_afdb_batch(batch_name, config, verbose=False)
    
    if result['success_count'] > 0:
        print(f"✅ Batch completed: {result['success_count']}/{result['total_proteins']} successful")
        print(f"   Total domains: {result['total_domains']}")
        print(f"   Output: {result['output_file']}")
    else:
        print(f"❌ Batch failed: {result['failed_count']} failures")
    
    return result['success_count'] > 0

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='AFDB processor with fixed paths')
    parser.add_argument('--batch', help='Process specific batch')
    parser.add_argument('--test', action='store_true', help='Run tests (default)')
    
    args = parser.parse_args()
    
    if args.batch:
        success = run_full_batch(args.batch)
        sys.exit(0 if success else 1)
    else:
        # Default: run tests
        sys.exit(main())
