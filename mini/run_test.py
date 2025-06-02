# mini/run_test.py
"""Test runner for mini_pyecod with proper chain BLAST decomposition"""

import sys
import os
from pathlib import Path

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.decomposer import load_domain_definitions
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

def test_8ovp():
    """Test the GFP+PBP fusion case (expecting 3 domains)"""

    # Load reference lengths for domains
    ref_lengths = {}
    ref_csv = "test_data/domain_lengths.csv"
    if os.path.exists(ref_csv):
        ref_lengths = load_reference_lengths(ref_csv)
    else:
        print("Warning: No domain reference lengths file found")

    # Load protein/chain lengths
    protein_lengths = {}
    protein_csv = "test_data/protein_lengths.csv"
    if os.path.exists(protein_csv):
        protein_lengths = load_protein_lengths(protein_csv)
    else:
        print("Warning: No protein lengths file found")

    # Load domain definitions for decomposition
    domain_defs = {}
    domain_defs_csv = "test_data/domain_definitions.csv"
    if os.path.exists(domain_defs_csv):
        domain_defs = load_domain_definitions(domain_defs_csv)
    else:
        print("Warning: No domain definitions file found - chain BLAST decomposition disabled")

    # Parse evidence
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    evidence = parse_domain_summary(xml_path, ref_lengths, protein_lengths)
    print(f"\nFound {len(evidence)} total evidence items")

    # Show evidence type distribution
    type_counts = {}
    for ev in evidence:
        type_counts[ev.type] = type_counts.get(ev.type, 0) + 1
    print("\nEvidence types:")
    for ev_type, count in sorted(type_counts.items()):
        print(f"  {ev_type}: {count}")

    # Show chain BLAST hits with alignment data
    chain_blast_hits = [e for e in evidence if e.type == "chain_blast"]
    print(f"\nChain BLAST hits: {len(chain_blast_hits)}")
    for cb in chain_blast_hits:
        has_align = "with alignment" if cb.alignment else "no alignment"
        print(f"  {cb.source_pdb}: coverage={cb.alignment_coverage:.1%} ({has_align})")

    # Filter to high-quality evidence
    good_evidence = [e for e in evidence if e.confidence > 0.6 or (e.evalue and e.evalue < 1e-5)]
    print(f"\nFiltered to {len(good_evidence)} high-quality evidence items")

    # Show coverage statistics
    coverage_stats = {}
    for ev in good_evidence:
        if ev.alignment_coverage is not None and ev.alignment_coverage > 0:
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

    # Partition with residue blocking and decomposition
    print(f"\nPartitioning with residue blocking (chain BLAST highest precedence)...")
    domains = partition_domains(good_evidence, sequence_length=518,
                               domain_defs=domain_defs, use_precedence=True)

    print(f"\nFound {len(domains)} domains:")
    for domain in domains:
        print(f"  Domain {domain.id}:")
        print(f"    Family: {domain.family}")
        print(f"    Range: {domain.range}")
        print(f"    Size: {domain.range.total_length} residues")
        print(f"    Discontinuous: {domain.range.is_discontinuous}")
        print(f"    Source: {domain.source}")

    # Check expected domains
    print("\nExpected domain structure:")
    print("  GFP: 252-494 (243 residues)")
    print("  PBP1: 108-205 (98 residues)")
    print("  PBP2: 4-107,206-251,495-512 (168 residues, discontinuous)")

    # Check results
    if len(domains) == 3:
        print("\n✅ SUCCESS: Found expected 3 domains!")

        # Check for discontinuous domain
        discontinuous_domains = [d for d in domains if d.range.is_discontinuous]
        if discontinuous_domains:
            print(f"  Found {len(discontinuous_domains)} discontinuous domain(s)")
            for dd in discontinuous_domains:
                print(f"    {dd.family}: {dd.range}")
    else:
        print(f"\n⚠️  Found {len(domains)} domains (expected 3)")
    
    # Write output
    output_path = "/tmp/8ovp_A_mini.domains.xml"
    write_domain_partition(domains, "8ovp", "A", output_path)
    print(f"\nWrote output to: {output_path}")

if __name__ == "__main__":
    test_8ovp()
