# mini/run_test.py
"""Simple test runner for mini_pyecod"""

import sys
import os
from pathlib import Path

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths  # Add this import
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

def test_8ovp():
    """Test with proper reference lengths"""

    # Load reference lengths
    ref_lengths = {}
    ref_csv = "test_data/domain_lengths.csv"
    if os.path.exists(ref_csv):
        ref_lengths = load_reference_lengths(ref_csv)
    else:
        print("Warning: No reference lengths file found")

    # Parse with reference lengths
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    evidence = parse_domain_summary(xml_path, ref_lengths)

    # Show coverage stats
    coverage_stats = {}
    for ev in evidence:
        if ev.alignment_coverage > 0:
            bucket = f"{int(ev.alignment_coverage * 10) * 10}%"
            coverage_stats[bucket] = coverage_stats.get(bucket, 0) + 1

    print("\nAlignment coverage distribution:")
    for bucket, count in sorted(coverage_stats.items()):
        print(f"  {bucket}: {count} hits")

    # Partition with blocking
    domains = partition_domains(good_evidence, sequence_length=518)

    print(f"\nFound {len(domains)} domains:")
    for domain in domains:
        print(f"  Domain {domain.id}:")
        print(f"    Family: {domain.family}")
        print(f"    Range: {domain.range}")
        print(f"    Size: {domain.range.size} residues")

    # Should be close to 3!
    if len(domains) == 3:
        print("\n✅ SUCCESS: Found expected 3 domains!")
    else:
        print(f"\n⚠️  Found {len(domains)} domains (expected 3)")

if __name__ == "__main__":
    test_8ovp()
