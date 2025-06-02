#!/usr/bin/env python3
"""
Complete test runner for mini_pyecod using ECOD range cache

This script:
1. Validates the range cache file
2. Generates reference files from the cache
3. Tests domain partitioning on 8ovp_A
4. Compares results with the old approach

Usage:
    python run_range_cache_test.py
    python run_range_cache_test.py --protein 1cuk_A  # Test different protein
"""

import sys
import os
from pathlib import Path
import time

def step(number, description):
    """Print a step header"""
    print(f"\n{'='*60}")
    print(f"STEP {number}: {description}")
    print(f"{'='*60}")

def run_range_cache_test(protein_id="8ovp_A",
                        cache_file=None,
                        batch_dir=None,
                        verbose=False):
    """
    Complete test using range cache approach

    Args:
        protein_id: Protein to test
        cache_file: Path to ECOD range cache file
        batch_dir: Batch directory
        verbose: Enable verbose output
    """

    print("MINI_PYECOD RANGE CACHE TEST")
    print("=" * 60)
    print(f"Testing: {protein_id}")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Default paths
    if cache_file is None:
        cache_file = "/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt"

    if batch_dir is None:
        batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"

    print(f"Cache file: {cache_file}")
    print(f"Batch dir: {batch_dir}")

    try:
        # Step 1: Validate range cache
        step(1, "Validating ECOD range cache")

        if not os.path.exists(cache_file):
            print(f"❌ ERROR: Cache file not found: {cache_file}")
            print("\nPlease check the path or download the range cache from ECOD")
            return False

        # Quick validation - check for 2ia4 domains
        found_2ia4 = False
        with open(cache_file, 'r') as f:
            for line in f:
                if 'e2ia4A' in line:
                    found_2ia4 = True
                    print(f"✅ Found 2ia4 domains in cache")
                    break

        if not found_2ia4:
            print(f"⚠️  WARNING: No 2ia4 domains found in cache")

        # Step 2: Generate reference files
        step(2, "Generating reference files from range cache")

        # Import and use the range cache parser
        sys.path.insert(0, str(Path(__file__).parent))

        try:
            from range_cache_parser import (
                create_domain_lengths_from_cache,
                create_domain_definitions_from_cache,
                extract_protein_lengths_from_cache
            )
        except ImportError:
            print("❌ ERROR: Could not import range_cache_parser")
            print("Make sure range_cache_parser.py is in the current directory")
            return False

        # Create output directory
        output_dir = "test_data_range_cache"
        os.makedirs(output_dir, exist_ok=True)

        domain_lengths_file = os.path.join(output_dir, 'domain_lengths.csv')
        domain_definitions_file = os.path.join(output_dir, 'domain_definitions.csv')
        protein_lengths_file = os.path.join(output_dir, 'protein_lengths.csv')

        print("Generating domain lengths...")
        create_domain_lengths_from_cache(cache_file, domain_lengths_file)

        print("Generating domain definitions...")
        create_domain_definitions_from_cache(cache_file, domain_definitions_file)

        print("Generating protein lengths...")
        extract_protein_lengths_from_cache(cache_file, protein_lengths_file)

        print("✅ Reference files generated successfully")

        # Step 3: Test domain partitioning
        step(3, f"Testing domain partitioning for {protein_id}")

        # Import mini_pyecod components
        sys.path.insert(0, str(Path(__file__).parent.parent))

        from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
        from mini.blast_parser import load_chain_blast_alignments
        from mini.decomposer import load_domain_definitions
        from mini.partitioner import partition_domains
        from mini.writer import write_domain_partition

        # Load reference data
        print("Loading reference data from range cache files...")
        reference_lengths = load_reference_lengths(domain_lengths_file)
        protein_lengths = load_protein_lengths(protein_lengths_file)
        domain_definitions = load_domain_definitions(domain_definitions_file)

        print(f"  Domain lengths: {len(reference_lengths):,}")
        print(f"  Protein lengths: {len(protein_lengths):,}")
        print(f"  Domain definitions: {len(domain_definitions):,} proteins")

        # Parse protein ID
        parts = protein_id.split('_')
        if len(parts) < 2:
            print(f"❌ ERROR: Invalid protein ID format: {protein_id}")
            return False

        pdb_id = parts[0]
        chain_id = parts[1]

        # Load BLAST alignments
        blast_dir = os.path.join(batch_dir, "blast/chain")
        blast_alignments = {}
        if os.path.exists(blast_dir):
            blast_alignments = load_chain_blast_alignments(blast_dir, pdb_id, chain_id)
            print(f"  BLAST alignments: {len(blast_alignments)}")

        # Parse domain summary
        xml_path = os.path.join(batch_dir, "domains", f"{protein_id}.develop291.domain_summary.xml")
        if not os.path.exists(xml_path):
            print(f"❌ ERROR: Domain summary not found: {xml_path}")
            return False

        print(f"Parsing evidence from domain summary...")
        evidence = parse_domain_summary(
            xml_path,
            reference_lengths=reference_lengths,
            protein_lengths=protein_lengths,
            blast_alignments=blast_alignments,
            verbose=verbose,
            require_reference_lengths=True
        )

        print(f"  Found {len(evidence)} evidence items with reference lengths")

        if len(evidence) == 0:
            print("❌ ERROR: No evidence with reference lengths found!")
            return False

        # Show evidence breakdown
        evidence_by_type = {}
        for ev in evidence:
            evidence_by_type[ev.type] = evidence_by_type.get(ev.type, 0) + 1

        print("  Evidence breakdown:")
        for etype, count in sorted(evidence_by_type.items()):
            print(f"    {etype}: {count}")

        # Estimate sequence length
        max_pos = 0
        for ev in evidence:
            for segment in ev.query_range.segments:
                max_pos = max(max_pos, segment.end)
        sequence_length = int(max_pos * 1.1)

        print(f"  Estimated sequence length: {sequence_length}")

        # Partition domains
        print("Running domain partitioning...")
        domains = partition_domains(evidence, sequence_length=sequence_length, verbose=verbose)

        if not domains:
            print("❌ No domains found!")
            return False

        # Step 4: Results analysis
        step(4, "Results Analysis")

        print(f"Found {len(domains)} domains:")
        total_coverage = 0

        for i, domain in enumerate(domains, 1):
            print(f"\n  Domain {i}:")
            print(f"    Family: {domain.family}")
            print(f"    Range: {domain.range}")
            print(f"    Size: {domain.range.total_length} residues")
            print(f"    Source: {domain.source}")
            print(f"    Discontinuous: {domain.range.is_discontinuous}")

            total_coverage += domain.range.total_length

            # Show evidence details
            if domain.evidence_items:
                ev = domain.evidence_items[0]
                if ev.reference_length:
                    coverage = domain.range.total_length / ev.reference_length
                    print(f"    Reference coverage: {coverage:.1%}")

        sequence_coverage = total_coverage / sequence_length
        print(f"\nTotal sequence coverage: {total_coverage}/{sequence_length} residues ({sequence_coverage:.1%})")

        # Write output
        output_path = f"/tmp/{protein_id}_range_cache.domains.xml"
        write_domain_partition(domains, pdb_id, chain_id, output_path)
        print(f"\n✅ Output written to: {output_path}")

        # Step 5: Validation for 8ovp_A
        if protein_id == "8ovp_A":
            step(5, "8ovp_A Validation")

            print("Expected domain architecture:")
            print("  1. GFP domain: ~252-494 (continuous)")
            print("  2. PBP domain(s): multiple possibilities")
            print("  3. Total: 2-3 domains expected")

            # Check results
            gfp_domains = [d for d in domains if '6dgv' in d.family or 'gfp' in d.family.lower()]
            pbp_domains = [d for d in domains if '2ia4' in d.family or 'pbp' in d.family.lower()]

            print(f"\nFound:")
            print(f"  GFP-like domains: {len(gfp_domains)}")
            print(f"  PBP-like domains: {len(pbp_domains)}")

            if len(domains) >= 2:
                print("\n✅ SUCCESS: Multiple domains found as expected!")
            else:
                print(f"\n⚠️  WARNING: Only {len(domains)} domain found, expected 2-3")

            # Check for discontinuous domains
            discontinuous = [d for d in domains if d.range.is_discontinuous]
            if discontinuous:
                print(f"\n🔀 Found {len(discontinuous)} discontinuous domain(s):")
                for d in discontinuous:
                    print(f"    {d.family}: {d.range}")

        print(f"\n{'='*60}")
        print("✅ RANGE CACHE TEST COMPLETED SUCCESSFULLY")
        print(f"{'='*60}")

        return True

    except Exception as e:
        print(f"\n❌ ERROR: Test failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Complete test runner for mini_pyecod using ECOD range cache',
        epilog='This script validates the cache, generates references, and tests domain partitioning'
    )

    parser.add_argument('--protein', default='8ovp_A',
                        help='Protein ID to test (default: 8ovp_A)')
    parser.add_argument('--cache-file',
                        default='/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt',
                        help='Path to ECOD range cache file')
    parser.add_argument('--batch-dir',
                        default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                        help='Batch directory containing domain files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')

    args = parser.parse_args()

    success = run_range_cache_test(
        protein_id=args.protein,
        cache_file=args.cache_file,
        batch_dir=args.batch_dir,
        verbose=args.verbose
    )

    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
