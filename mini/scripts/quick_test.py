#!/usr/bin/env python3
"""Quick test runner for mini_pyecod"""

import sys
import os
from pathlib import Path

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

def run_test(protein_id, batch_dir=None, verbose=False):
    """Run domain partitioning test for a specific protein"""

    if batch_dir is None:
        batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"

    # Parse PDB and chain from protein_id
    parts = protein_id.split('_')
    if len(parts) >= 2:
        pdb_id = parts[0]
        chain_id = parts[1]
    else:
        print(f"Error: Invalid protein ID format: {protein_id}")
        print("Expected format: PDBID_CHAIN (e.g., 8ovp_A)")
        return

    # Construct file path
    xml_path = os.path.join(batch_dir, "domains", f"{protein_id}.develop291.domain_summary.xml")

    if not os.path.exists(xml_path):
        print(f"Error: Domain summary file not found: {xml_path}")
        return

    print(f"\n{'='*60}")
    print(f"Testing: {protein_id}")
    print(f"{'='*60}")

    try:
        # Parse evidence
        evidence = parse_domain_summary(xml_path, verbose=verbose)
        print(f"\nFound {len(evidence)} total evidence items")

        # Show evidence distribution
        evidence_types = {}
        families = {}
        for ev in evidence:
            # Count by type
            evidence_types[ev.type] = evidence_types.get(ev.type, 0) + 1

            # Count by family/source
            family = ev.source_pdb or "unknown"
            families[family] = families.get(family, 0) + 1

        print("\nEvidence by type:")
        for etype, count in sorted(evidence_types.items()):
            print(f"  {etype}: {count}")

        # Get sequence length from evidence ranges
        max_pos = 0
        for ev in evidence:
            for start, end in ev.query_range.segments:
                max_pos = max(max_pos, end)

        # Add 10% buffer to sequence length estimate
        sequence_length = int(max_pos * 1.1)
        print(f"\nEstimated sequence length: {sequence_length}")

        # Filter to high-quality evidence
        good_evidence = [e for e in evidence if e.confidence > 0.5 or (e.evalue and e.evalue < 1e-3)]
        print(f"Filtered to {len(good_evidence)} high-quality evidence items")

        # Show top families
        print(f"\nTop families by evidence count:")
        for family, count in sorted(families.items(), key=lambda x: -x[1])[:10]:
            print(f"  {family}: {count} hits")

        # Partition domains
        print(f"\n{'='*40}")
        print("DOMAIN PARTITIONING")
        print(f"{'='*40}")

        domains = partition_domains(good_evidence, sequence_length=sequence_length, verbose=verbose)

        # Results summary
        print(f"\n{'='*40}")
        print(f"RESULTS: Found {len(domains)} domains")
        print(f"{'='*40}")

        for domain in domains:
            print(f"\nDomain {domain.id}:")
            print(f"  Family: {domain.family}")
            print(f"  Range: {domain.range}")
            print(f"  Size: {domain.range.total_length} residues")
            print(f"  Source: {domain.source}")
            print(f"  Discontinuous: {domain.range.is_discontinuous}")

        # Write output
        output_path = f"/tmp/{protein_id}_mini.domains.xml"
        write_domain_partition(domains, pdb_id, chain_id, output_path)
        print(f"\n✓ Wrote output to: {output_path}")

        return domains

    except Exception as e:
        print(f"\nError processing {protein_id}: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(description='Quick test runner for mini_pyecod')
    parser.add_argument('protein_id', help='Protein ID to test (e.g., 8ovp_A)')
    parser.add_argument('--batch-dir', help='Batch directory containing domain files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')

    args = parser.parse_args()

    run_test(args.protein_id, args.batch_dir, verbose=args.verbose)

if __name__ == "__main__":
    main()
