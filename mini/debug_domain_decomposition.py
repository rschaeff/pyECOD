#!/usr/bin/env python3
"""Debug why chain BLAST decomposition isn't working"""

import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.blast_parser import load_chain_blast_alignments
from mini.decomposer import load_domain_definitions

def debug_decomposition():
    """Debug the decomposition process step by step"""
    
    # 1. Load BLAST alignments
    blast_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/blast/chain"
    blast_alignments = load_chain_blast_alignments(blast_dir, "8ovp", "A")
    print(f"1. Loaded {len(blast_alignments)} BLAST alignments")
    
    # Check for 2ia4
    for (pdb, chain), align in blast_alignments.items():
        if pdb == '2ia4':
            print(f"   Found alignment for {pdb}_{chain}")
    
    # 2. Load domain definitions
    domain_defs = load_domain_definitions("test_data/domain_definitions.csv")
    print(f"\n2. Loaded domain definitions for {len(domain_defs)} proteins")
    
    # Check for 2ia4
    for (pdb, chain) in domain_defs.keys():
        if pdb == '2ia4':
            print(f"   Found definitions for {pdb}_{chain}: {len(domain_defs[(pdb, chain)])} domains")
    
    # 3. Parse evidence
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    evidence = parse_domain_summary(xml_path, {}, {}, blast_alignments)
    
    print(f"\n3. Parsed {len(evidence)} evidence items")
    
    # 4. Check chain BLAST hits to 2ia4
    print("\n4. Chain BLAST hits to 2ia4:")
    chain_blast_2ia4 = []
    for ev in evidence:
        if ev.type == "chain_blast" and ev.source_pdb == "2ia4":
            chain_blast_2ia4.append(ev)
            print(f"   {ev.source_pdb}: {ev.query_range}")
            print(f"   - Has alignment: {ev.alignment is not None}")
            print(f"   - Domain ID: {ev.domain_id}")
            
            # Check hit key
            hit_key = (ev.source_pdb, ev.domain_id.split('_')[-1])
            print(f"   - Hit key: {hit_key}")
            print(f"   - Key in domain_defs: {hit_key in domain_defs}")
            
            if ev.alignment:
                print(f"   - Alignment query range: {ev.alignment.query_start}-{ev.alignment.query_end}")
                print(f"   - Alignment hit range: {ev.alignment.hit_start}-{ev.alignment.hit_end}")
    
    # 5. Test decomposition directly
    if chain_blast_2ia4:
        print("\n5. Testing decomposition directly:")
        from mini.decomposer import decompose_chain_blast_with_mapping
        
        ev = chain_blast_2ia4[0]
        hit_key = (ev.source_pdb, ev.domain_id.split('_')[-1])
        
        if ev.alignment and hit_key in domain_defs:
            print(f"   Attempting to decompose {ev.source_pdb}_{hit_key[1]}")
            decomposed = decompose_chain_blast_with_mapping(
                ev,
                ev.alignment.query_seq,
                ev.alignment.hit_seq,
                ev.alignment.query_start,
                ev.alignment.hit_start,
                domain_defs[hit_key]
            )
            print(f"   Result: {len(decomposed)} domains")
            for d in decomposed:
                print(f"     - {d.domain_id}: {d.query_range}")

if __name__ == "__main__":
    debug_decomposition()
