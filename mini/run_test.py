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
    """Test the GFP+PBP fusion case"""
    # You'll need to adjust this path
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    
    print(f"Testing: {xml_path}")
    
    # Parse evidence
    evidence = parse_domain_summary(xml_path)
    print(f"\nFound {len(evidence)} evidence items")
    
    # Show unique source families
    families = set()
    for ev in evidence:
        if ev.source_pdb:
            families.add(ev.source_pdb)
    print(f"Unique source PDBs: {sorted(families)[:10]}...")  # First 10
    
    # Partition
    domains = partition_domains(evidence, sequence_length=518)
    
    print(f"\nFound {len(domains)} domains:")
    for domain in domains[:5]:  # Show first 5
        print(f"  Domain {domain.id}:")
        print(f"    Family: {domain.family}")
        print(f"    Range: {domain.range}")
        print(f"    Evidence: {domain.evidence_count} items")
        print(f"    Source: {domain.source}")
    
    # Check for GFP and PBP
    gfp_domains = [d for d in domains if '6dgv' in d.family]
    pbp_domains = [d for d in domains if '2ia4' in d.family]
    
    print(f"\nGFP domains (6dgv): {len(gfp_domains)}")
    print(f"PBP domains (2ia4): {len(pbp_domains)}")
    
    # Write output
    output_path = "/tmp/8ovp_A_mini.domains.xml"
    write_domain_partition(domains, "8ovp", "A", output_path)
    print(f"\nWrote output to: {output_path}")

if __name__ == "__main__":
    test_8ovp()
