#!/usr/bin/env python3
"""Quick test runner for individual test cases"""

import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from mini import (
    parse_domain_summary, 
    load_reference_lengths, 
    load_protein_lengths,
    load_chain_blast_alignments,
    load_domain_definitions,
    partition_domains,
    write_domain_partition
)

def quick_test(pdb_chain: str):
    """Quick test of a single PDB chain"""
    
    if '_' not in pdb_chain:
        print(f"Usage: {sys.argv[0]} PDB_CHAIN (e.g., 1cuk_A)")
        return
    
    pdb_id, chain_id = pdb_chain.split('_', 1)
    print(f"\nQuick test: {pdb_id}_{chain_id}")
    print("=" * 50)
    
    # Paths
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    xml_path = f"{batch_dir}/domains/{pdb_id}_{chain_id}.develop291.domain_summary.xml"
    
    if not os.path.exists(xml_path):
        print(f"❌ File not found: {xml_path}")
        return
    
    # Load minimal reference data
    domain_defs = {}
    if os.path.exists("test_data/domain_definitions.csv"):
        domain_defs = load_domain_definitions("test_data/domain_definitions.csv")
    
    # Load BLAST alignments
    blast_alignments = {}
    blast_dir = f"{batch_dir}/blast/chain"
    if os.path.exists(blast_dir):
        try:
            blast_alignments = load_chain_blast_alignments(blast_dir, pdb_id, chain_id)
            print(f"✓ Loaded {len(blast_alignments)} BLAST alignments")
        except:
            print("✗ Could not load BLAST alignments")
    
    # Parse evidence
    evidence = parse_domain_summary(xml_path, {}, {}, blast_alignments)
    print(f"✓ Found {len(evidence)} evidence items")
    
    # Show evidence summary
    evidence_types = {}
    for ev in evidence:
        evidence_types[ev.type] = evidence_types.get(ev.type, 0) + 1
    print(f"  Evidence types: {dict(evidence_types)}")
    
    # Get sequence length
    max_pos = 1
    for ev in evidence:
        try:
            positions = ev.query_range.to_positions_simple()
            if positions:
                max_pos = max(max_pos, max(positions))
        except:
            for pos, _ in ev.query_range.to_positions():
                max_pos = max(max_pos, pos)
    
    sequence_length = max_pos + 50
    print(f"  Estimated length: {sequence_length}")
    
    # Filter evidence
    good_evidence = [e for e in evidence if e.confidence > 0.6 or (e.evalue and e.evalue < 1e-5)]
    print(f"✓ Filtered to {len(good_evidence)} high-quality items")
    
    # Partition
    print(f"\nPartitioning...")
    domains = partition_domains(good_evidence, sequence_length=sequence_length,
                               domain_defs=domain_defs, use_precedence=True)
    
    # Results
    print(f"\n{'='*50}")
    print(f"RESULTS: Found {len(domains)} domains")
    print(f"{'='*50}")
    
    for i, domain in enumerate(domains, 1):
        disc = "discontinuous" if domain.range.is_discontinuous else "continuous"
        print(f"{i}. {domain.family}")
        print(f"   Range: {domain.range} ({domain.range.total_length} residues)")
        print(f"   Type: {disc}")
        print(f"   Source: {domain.source}")
    
    # Check for potential issues
    print(f"\nPotential issues:")
    
    # Check for very small domains
    small_domains = [d for d in domains if d.range.total_length < 30]
    if small_domains:
        print(f"  ⚠️  {len(small_domains)} very small domains (<30 residues)")
    
    # Check for gaps in coverage
    all_positions = set()
    for d in domains:
        all_positions.update(d.range.to_positions_simple())
    
    coverage = len(all_positions) / sequence_length * 100
    print(f"  Coverage: {coverage:.1f}% of sequence")
    
    if coverage < 80:
        print(f"  ⚠️  Low coverage - possible missing domains")
    
    # Save output
    output = f"/tmp/{pdb_id}_{chain_id}_quick.domains.xml"
    write_domain_partition(domains, pdb_id, chain_id, output)
    print(f"\nOutput: {output}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        quick_test(sys.argv[1])
    else:
        print(f"Usage: {sys.argv[0]} PDB_CHAIN")
        print(f"Example: {sys.argv[0]} 1cuk_A")
