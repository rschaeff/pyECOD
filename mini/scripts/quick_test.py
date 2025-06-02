#!/usr/bin/env python3
"""Quick test runner for mini_pyecod - Both domain and protein reference lengths are MANDATORY"""

import sys
import os
from pathlib import Path
from typing import Any

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

def run_test(protein_id, batch_dir=None, verbose=False, domain_lengths_file=None, protein_lengths_file=None):
    """
    Run domain partitioning test for a specific protein.

    Both domain and protein reference lengths are MANDATORY for accurate domain partitioning.
    """

    if batch_dir is None:
        batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"

    # MANDATORY: Both reference files are required
    if not domain_lengths_file:
        print("ERROR: Domain reference lengths are REQUIRED for domain partitioning")
        print("\nUsage: python quick_test.py PROTEIN_ID --domain-lengths DOMAIN_FILE --protein-lengths PROTEIN_FILE")
        print("\nExample:")
        print("  python quick_test.py 8ovp_A \\")
        print("    --domain-lengths test_data_range_cache/domain_lengths.csv \\")
        print("    --protein-lengths test_data_range_cache/protein_lengths.csv")
        return None

    if not protein_lengths_file:
        print("ERROR: Protein reference lengths are REQUIRED for chain BLAST evidence")
        print("\nUsage: python quick_test.py PROTEIN_ID --domain-lengths DOMAIN_FILE --protein-lengths PROTEIN_FILE")
        print("\nExample:")
        print("  python quick_test.py 8ovp_A \\")
        print("    --domain-lengths test_data_range_cache/domain_lengths.csv \\")
        print("    --protein-lengths test_data_range_cache/protein_lengths.csv")
        return None

    # Verify both files exist
    for file_path, file_type in [(domain_lengths_file, "domain lengths"), (protein_lengths_file, "protein lengths")]:
        if not os.path.exists(file_path):
            print(f"ERROR: {file_type.title()} file not found: {file_path}")
            print(f"\nPlease provide a valid {file_type} CSV file.")
            return None

    # Load both reference files
    print(f"Loading domain reference lengths from {domain_lengths_file}")
    domain_lengths = load_reference_lengths(domain_lengths_file)

    print(f"Loading protein reference lengths from {protein_lengths_file}")
    protein_lengths = load_protein_lengths(protein_lengths_file)

    if not domain_lengths:
        print("ERROR: No domain reference lengths loaded. Check file format.")
        print("Expected format: domain_id,length")
        print("Example: e6dgvA1,245")
        return None

    if not protein_lengths:
        print("ERROR: No protein reference lengths loaded. Check file format.")
        print("Expected format: pdb_id,chain_id,length")
        print("Example: 8ovp,A,518")
        return None

    print(f"Loaded {len(domain_lengths)} domain reference lengths")
    print(f"Loaded {len(protein_lengths)} protein reference lengths")

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
        # Parse evidence WITH BOTH reference length types
        evidence = parse_domain_summary(xml_path,
                                      reference_lengths=domain_lengths,
                                      protein_lengths=protein_lengths,
                                      require_reference_lengths=True,  # Skip evidence without refs
                                      verbose=verbose)

        print(f"\nFound {len(evidence)} evidence items with reference lengths")

        if len(evidence) == 0:
            print("\nERROR: No evidence with reference lengths found.")
            print("This means the reference files don't contain entries for the domains/proteins hitting this protein.")
            print("\nPossible solutions:")
            print("1. Use more complete reference files (e.g., full ECOD domain and protein lengths)")
            print("2. Check that the reference files have the expected domain/protein IDs")
            return None

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
            for segment in ev.query_range.segments:
                max_pos = max(max_pos, segment.end)

        # Add 10% buffer to sequence length estimate
        sequence_length = int(max_pos * 1.1)
        print(f"\nEstimated sequence length: {sequence_length}")

        # Count evidence with alignment coverage
        with_coverage = sum(1 for e in evidence if e.alignment_coverage is not None)
        print(f"Evidence with alignment coverage: {with_coverage}/{len(evidence)}")

        # Show top families
        print(f"\nTop families by evidence count:")
        for family, count in sorted(families.items(), key=lambda x: -x[1])[:10]:
            print(f"  {family}: {count} hits")

        # Partition domains - partitioner will handle evidence precedence
        print(f"\n{'='*40}")
        print("DOMAIN PARTITIONING")
        print(f"{'='*40}")

        domains = partition_domains(evidence, sequence_length=sequence_length, verbose=verbose)

        if not domains:
            print("\nNo domains found. This could mean:")
            print("1. The evidence doesn't meet coverage thresholds")
            print("2. All evidence was rejected due to overlaps")
            print("Run with --verbose to see rejection details")
            return None

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

            # Show evidence details
            if domain.evidence_items:
                ev = domain.evidence_items[0]
                if ev.reference_length:
                    print(f"  Reference length: {ev.reference_length} residues")
                if ev.alignment_coverage is not None:
                    print(f"  Alignment coverage: {ev.alignment_coverage:.1%}")

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

    parser = argparse.ArgumentParser(
        description='Quick test runner for mini_pyecod - Both domain and protein reference lengths are MANDATORY',
        epilog='Both domain and protein reference lengths are required for accurate domain partitioning.'
    )
    parser.add_argument('protein_id', help='Protein ID to test (e.g., 8ovp_A)')
    parser.add_argument('--batch-dir', help='Batch directory containing domain files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    parser.add_argument('--domain-lengths', required=True,
                        help='CSV file with domain reference lengths (REQUIRED)')
    parser.add_argument('--protein-lengths', required=True,
                        help='CSV file with protein reference lengths (REQUIRED)')

    args = parser.parse_args()

    run_test(args.protein_id, args.batch_dir,
             verbose=args.verbose,
             domain_lengths_file=args.domain_lengths,
             protein_lengths_file=args.protein_lengths)

if __name__ == "__main__":
    main()
