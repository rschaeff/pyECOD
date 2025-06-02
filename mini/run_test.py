# mini/run_test.py
"""Simple test runner for mini_pyecod"""

import sys
import os
from pathlib import Path

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

def test_8ovp():
    """Test with residue blocking"""
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"

    evidence = parse_domain_summary(xml_path)
    print(f"Found {len(evidence)} evidence items")

    # Filter to high-quality evidence first
    good_evidence = [e for e in evidence if e.confidence > 0.6 or (e.evalue and e.evalue < 1e-5)]
    print(f"Filtered to {len(good_evidence)} high-quality evidence items")

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
