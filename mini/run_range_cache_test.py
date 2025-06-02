#!/usr/bin/env python3
"""
Test mini_pyecod using ECOD range cache as reference source

This script demonstrates using the authoritative ECOD range cache file
instead of database dumps for reference lengths and domain definitions.
"""

import sys
import os
from pathlib import Path

# Add parent directory to path so we can import ecod
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.blast_parser import load_chain_blast_alignments
from mini.decomposer import load_domain_definitions
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

def setup_range_cache_references(cache_file: str, force_regenerate: bool = False):
    """
    Set up reference files from range cache if they don't exist or if forced.
    
    Args:
        cache_file: Path to ECOD range cache file
        force_regenerate: Whether to regenerate files even if they exist
        
    Returns:
        Tuple of (domain_lengths_file, domain_definitions_file, protein_lengths_file)
    """
    # Import the range cache parser
    from range_cache_parser import (
        create_domain_lengths_from_cache,
        create_domain_definitions_from_cache, 
        extract_protein_lengths_from_cache
    )
    
    output_dir = "test_data_v291"
    os.makedirs(output_dir, exist_ok=True)
    
    domain_lengths_file = os.path.join(output_dir, 'domain_lengths.csv')
    domain_definitions_file = os.path.join(output_dir, 'domain_definitions.csv')
    protein_lengths_file = os.path.join(output_dir, 'protein_lengths.csv')
    
    # Check if files exist and are recent
    files_exist = all(os.path.exists(f) for f in [domain_lengths_file, domain_definitions_file, protein_lengths_file])
    
    if not files_exist or force_regenerate:
        print("Generating reference files from ECOD range cache...")
        print("This may take a few moments...")
        
        create_domain_lengths_from_cache(cache_file, domain_lengths_file)
        create_domain_definitions_from_cache(cache_file, domain_definitions_file)
        extract_protein_lengths_from_cache(cache_file, protein_lengths_file)
        
        print("Reference files generated successfully!")
    else:
        print("Using existing reference files from range cache")
    
    return domain_lengths_file, domain_definitions_file, protein_lengths_file

def test_protein_with_range_cache(protein_id: str, 
                                  batch_dir: str = None,
                                  cache_file: str = None,
                                  verbose: bool = False):
    """
    Test domain partitioning using range cache references.
    
    Args:
        protein_id: Protein to test (e.g., '8ovp_A')
        batch_dir: Batch directory containing domain files
        cache_file: Path to ECOD range cache file
        verbose: Whether to enable verbose output
    """
    if batch_dir is None:
        batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    
    if cache_file is None:
        cache_file = "/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt"
    
    print(f"\n{'='*60}")
    print(f"Testing {protein_id} with ECOD v291 range cache references")
    print(f"{'='*60}")
    
    # Verify cache file exists
    if not os.path.exists(cache_file):
        print(f"ERROR: Range cache file not found: {cache_file}")
        print("Please check the path or download from ECOD database")
        return None
    
    # Set up reference files from cache
    try:
        domain_lengths_file, domain_definitions_file, protein_lengths_file = setup_range_cache_references(cache_file)
    except Exception as e:
        print(f"ERROR: Failed to set up reference files: {e}")
        return None
    
    # Load reference data
    print("\nLoading reference data...")
    reference_lengths = load_reference_lengths(domain_lengths_file)
    protein_lengths = load_protein_lengths(protein_lengths_file)
    domain_definitions = load_domain_definitions(domain_definitions_file, verbose)
    
    print(f"Loaded {len(reference_lengths)} domain lengths")
    print(f"Loaded {len(protein_lengths)} protein lengths")
    print(f"Loaded domain definitions for {len(domain_definitions)} proteins")
    
    # Parse PDB and chain
    parts = protein_id.split('_')
    if len(parts) < 2:
        print(f"ERROR: Invalid protein ID format: {protein_id}")
        return None
    
    pdb_id = parts[0]
    chain_id = parts[1]
    
    # Load BLAST alignments
    blast_dir = os.path.join(batch_dir, "blast/chain")
    blast_alignments = {}
    if os.path.exists(blast_dir):
        blast_alignments = load_chain_blast_alignments(blast_dir, pdb_id, chain_id, verbose)
        print(f"Loaded {len(blast_alignments)} BLAST alignments")
    
    # Parse domain summary
    xml_path = os.path.join(batch_dir, "domains", f"{protein_id}.develop291.domain_summary.xml")
    if not os.path.exists(xml_path):
        print(f"ERROR: Domain summary not found: {xml_path}")
        return None
    
    print(f"\nParsing evidence from {xml_path}")
    evidence = parse_domain_summary(
        xml_path,
        reference_lengths=reference_lengths,
        protein_lengths=protein_lengths,
        blast_alignments=blast_alignments,
        verbose=verbose,
        require_reference_lengths=True  # Only use evidence with reference lengths
    )
    
    print(f"Found {len(evidence)} evidence items with reference lengths")
    
    if len(evidence) == 0:
        print("ERROR: No evidence with reference lengths found!")
        print("This suggests the domain summary doesn't contain hits to domains in the range cache.")
        return None
    
    # Show evidence breakdown
    evidence_by_type = {}
    for ev in evidence:
        evidence_by_type[ev.type] = evidence_by_type.get(ev.type, 0) + 1
    
    print("\nEvidence breakdown:")
    for etype, count in sorted(evidence_by_type.items()):
        print(f"  {etype}: {count}")
    
    # Check for key domains (for 8ovp_A, we expect 2ia4 hits)
    if pdb_id == "8ovp":
        pbp_evidence = [e for e in evidence if e.source_pdb == "2ia4"]
        print(f"\nPBP evidence (2ia4): {len(pbp_evidence)} hits")
        for ev in pbp_evidence:
            print(f"  {ev.type}: {ev.domain_id} @ {ev.query_range}")
    
    # Get sequence length estimate
    max_pos = 0
    for ev in evidence:
        for segment in ev.query_range.segments:
            max_pos = max(max_pos, segment.end)
    sequence_length = int(max_pos * 1.1)  # Add 10% buffer
    
    print(f"\nEstimated sequence length: {sequence_length}")
    
    # Partition domains
    print("\n" + "="*40)
    print("DOMAIN PARTITIONING")
    print("="*40)
    
    domains = partition_domains(evidence, sequence_length=sequence_length, verbose=verbose)
    
    if not domains:
        print("No domains found!")
        return None
    
    # Results
    print(f"\n{'='*40}")
    print(f"RESULTS: Found {len(domains)} domains")
    print(f"{'='*40}")
    
    total_coverage = 0
    for i, domain in enumerate(domains, 1):
        print(f"\nDomain {i}:")
        print(f"  Family: {domain.family}")
        print(f"  Range: {domain.range}")
        print(f"  Size: {domain.range.total_length} residues")
        print(f"  Source: {domain.source}")
        print(f"  Discontinuous: {domain.range.is_discontinuous}")
        
        total_coverage += domain.range.total_length
        
        # Show evidence details
        if domain.evidence_items:
            ev = domain.evidence_items[0]
            if ev.reference_length:
                coverage = domain.range.total_length / ev.reference_length
                print(f"  Reference coverage: {coverage:.1%} ({domain.range.total_length}/{ev.reference_length})")
    
    print(f"\nTotal domain coverage: {total_coverage}/{sequence_length} residues ({total_coverage/sequence_length:.1%})")
    
    # Write output
    output_path = f"/tmp/{protein_id}_range_cache.domains.xml"
    write_domain_partition(domains, pdb_id, chain_id, output_path)
    print(f"\n✓ Wrote output to: {output_path}")
    
    return domains

def test_2ia4_coverage():
    """Test that we can find 2ia4 domains in the range cache"""
    from range_cache_parser import parse_range_cache
    
    cache_file = "/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt"
    
    if not os.path.exists(cache_file):
        print(f"Range cache file not found: {cache_file}")
        return
    
    print("Checking 2ia4 domain coverage in range cache...")
    
    entries = parse_range_cache(cache_file)
    
    # Find all 2ia4 domains
    ia4_domains = {k: v for k, v in entries.items() if v.pdb_id == "2ia4"}
    
    print(f"Found {len(ia4_domains)} domains for 2ia4:")
    for domain_id, entry in sorted(ia4_domains.items()):
        print(f"  {domain_id}: {entry.range_spec} ({entry.length} residues)")
    
    return ia4_domains

def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Test mini_pyecod with ECOD range cache references'
    )
    parser.add_argument('protein_id', nargs='?', default='8ovp_A',
                        help='Protein ID to test (default: 8ovp_A)')
    parser.add_argument('--batch-dir', 
                        default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                        help='Batch directory containing domain files')
    parser.add_argument('--cache-file',
                        default='/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt',
                        help='Path to ECOD range cache file')
    parser.add_argument('--test-2ia4', action='store_true',
                        help='Test 2ia4 domain coverage in cache')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    
    if args.test_2ia4:
        test_2ia4_coverage()
        return
    
    # Test the protein
    domains = test_protein_with_range_cache(
        args.protein_id,
        args.batch_dir,
        args.cache_file,
        args.verbose
    )
    
    if domains and args.protein_id == "8ovp_A":
        print("\n" + "="*60)
        print("8ovp_A VALIDATION")
        print("="*60)
        print("Expected domain structure:")
        print("  1. GFP domain: ~252-494 (continuous)")
        print("  2. PBP domains: multiple segments, possibly discontinuous")
        print("  3. Total domains: 2-3 expected")
        
        if len(domains) >= 2:
            print("\n✅ SUCCESS: Found multiple domains as expected!")
        else:
            print(f"\n⚠️  Found only {len(domains)} domain(s), expected 2-3")

if __name__ == "__main__":
    main()
