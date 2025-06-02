# mini/run_test.py
"""Test runner for mini_pyecod"""

import sys
import os
from pathlib import Path

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

def test_8ovp():
    """Test the GFP+PBP fusion case (expecting 3 domains)"""

    # Load reference lengths
    ref_lengths = {}
    ref_csv = "test_data/domain_lengths.csv"
    if os.path.exists(ref_csv):
        ref_lengths = load_reference_lengths(ref_csv)
    else:
        print("Warning: No reference lengths file found")

    # Parse evidence
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    evidence = parse_domain_summary(xml_path, ref_lengths)
    print(f"\nFound {len(evidence)} total evidence items")

    # Filter to high-quality evidence
    good_evidence = [e for e in evidence if e.confidence > 0.6 or (e.evalue and e.evalue < 1e-5)]
    print(f"Filtered to {len(good_evidence)} high-quality evidence items")

    # Show coverage statistics
    coverage_stats = {}
    for ev in good_evidence:
        if ev.alignment_coverage is not None and ev.alignment_coverage > 0:  # Check for None!
            bucket = f"{int(ev.alignment_coverage * 10) * 10}%"
            coverage_stats[bucket] = coverage_stats.get(bucket, 0) + 1

    if coverage_stats:
        print("\nAlignment coverage distribution:")
        for bucket, count in sorted(coverage_stats.items()):
            print(f"  {bucket}: {count} hits")

    # Show evidence distribution
    families = {}
    for ev in good_evidence:
        family = ev.source_pdb or "unknown"
        families[family] = families.get(family, 0) + 1
    print(f"\nTop families by evidence count:")
    for family, count in sorted(families.items(), key=lambda x: -x[1])[:10]:
        print(f"  {family}: {count} hits")

    # Partition with residue blocking
    print(f"\nPartitioning with residue blocking...")
    domains = partition_domains(good_evidence, sequence_length=518)

    print(f"\nFound {len(domains)} domains:")
    for domain in domains:
        print(f"  Domain {domain.id}:")
        print(f"    Family: {domain.family}")
        print(f"    Range: {domain.range}")
        print(f"    Size: {domain.range.total_length} residues")
        print(f"    Discontinuous: {domain.range.is_discontinuous}")

    # Check results
    if len(domains) == 3:
        print("\n✅ SUCCESS: Found expected 3 domains!")
    else:
        print(f"\n⚠️  Found {len(domains)} domains (expected 3)")
    
    # Write output
    output_path = "/tmp/8ovp_A_mini.domains.xml"
    write_domain_partition(domains, "8ovp", "A", output_path)
    print(f"\nWrote output to: {output_path}")

if __name__ == "__main__":
    test_8ovp()
