#!/usr/bin/env python3
"""Debug why certain evidence is selected first"""

import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.blast_parser import load_chain_blast_alignments
from mini.decomposer import load_domain_definitions
from mini.partitioner import partition_domains

def debug_evidence_ordering():
    """Check the ordering of evidence"""
    
    # Load all data
    ref_lengths = load_reference_lengths("test_data/domain_lengths.csv") if os.path.exists("test_data/domain_lengths.csv") else {}
    protein_lengths = load_protein_lengths("test_data/protein_lengths.csv") if os.path.exists("test_data/protein_lengths.csv") else {}
    domain_defs = load_domain_definitions("test_data/domain_definitions.csv") if os.path.exists("test_data/domain_definitions.csv") else {}
    
    blast_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/blast/chain"
    blast_alignments = load_chain_blast_alignments(blast_dir, "8ovp", "A")
    
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    evidence = parse_domain_summary(xml_path, ref_lengths, protein_lengths, blast_alignments)
    
    # Filter to good evidence
    good_evidence = [e for e in evidence if e.confidence > 0.6 or (e.evalue and e.evalue < 1e-5)]
    
    # Show PBP-related chain BLAST evidence
    print("PBP-related chain BLAST evidence (before decomposition):")
    pbp_evidence = []
    for ev in good_evidence:
        if ev.type == "chain_blast" and ev.source_pdb in ['2ia4', '2vha', '1pda']:
            pbp_evidence.append(ev)
            print(f"\n{ev.source_pdb}_{ev.domain_id.split('_')[-1]}:")
            print(f"  Range: {ev.query_range}")
            print(f"  E-value: {ev.evalue}")
            print(f"  Confidence: {ev.confidence}")
            print(f"  Coverage: {ev.alignment_coverage:.1%}")
            print(f"  Has alignment: {ev.alignment is not None}")
    
    # Now check what happens after decomposition
    print("\n\nChecking decomposition...")
    
    # Process chain BLAST hits like partitioner does
    processed_evidence = []
    for evidence in good_evidence:
        if evidence.type == "chain_blast":
            hit_key = (evidence.source_pdb, evidence.domain_id.split('_')[-1])
            
            if evidence.alignment is not None and hit_key in domain_defs:
                print(f"\nDecomposing {evidence.source_pdb}_{hit_key[1]}...")
                from mini.decomposer import decompose_chain_blast_with_mapping
                decomposed = decompose_chain_blast_with_mapping(
                    evidence,
                    evidence.alignment.query_seq,
                    evidence.alignment.hit_seq,
                    evidence.alignment.query_start,
                    evidence.alignment.hit_start,
                    domain_defs[hit_key]
                )
                for d in decomposed:
                    print(f"  -> {d.domain_id}: {d.query_range}")
                processed_evidence.extend(decomposed)
            else:
                processed_evidence.append(evidence)
        else:
            processed_evidence.append(evidence)
    
    # Show ordering after sorting
    print("\n\nEvidence ordering (top 20 by precedence):")
    precedence_map = {
        'chain_blast': 1,
        'chain_blast_decomposed': 1.5,
        'domain_blast': 2,
        'hhsearch': 3
    }
    
    sorted_evidence = sorted(processed_evidence,
                           key=lambda e: (precedence_map.get(e.type, 99), 
                                           e.evalue if e.evalue else 999))
    
    for i, ev in enumerate(sorted_evidence[:20]):
        if ev.source_pdb in ['2ia4', '2vha', '1pda', '6dgv', '5yr2'] or ev.type == 'chain_blast_decomposed':
            print(f"{i+1}. {ev.type}: {ev.source_pdb} @ {ev.query_range} (e={ev.evalue:.2e})")

if __name__ == "__main__":
    debug_evidence_ordering()
