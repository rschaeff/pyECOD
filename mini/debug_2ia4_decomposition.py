#!/usr/bin/env python3
"""
Debug 2ia4 decomposition step by step

We're still getting 2 domains instead of 3, which means 2ia4 isn't being decomposed.
Let's trace exactly what's happening in the decomposition pipeline.
"""

import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

def debug_2ia4_evidence():
    """Debug the 2ia4 evidence from parsing"""
    
    print("STEP 1: CHECKING 2ia4 EVIDENCE FROM PARSING")
    print("=" * 50)
    
    from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
    from mini.blast_parser import load_chain_blast_alignments
    from mini.decomposer import load_domain_definitions
    
    # Load all reference data
    reference_lengths = load_reference_lengths("test_data_v291/domain_lengths.csv")
    protein_lengths = load_protein_lengths("test_data_v291/protein_lengths.csv")
    domain_definitions = load_domain_definitions("test_data_v291/domain_definitions.csv")
    
    # Load BLAST alignments
    blast_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/blast/chain"
    blast_alignments = load_chain_blast_alignments(blast_dir, "8ovp", "A")
    
    print(f"Loaded {len(blast_alignments)} BLAST alignments")
    print(f"Loaded domain definitions for {len(domain_definitions)} proteins")
    
    # Parse evidence
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    evidence = parse_domain_summary(
        xml_path,
        reference_lengths=reference_lengths,
        protein_lengths=protein_lengths,
        blast_alignments=blast_alignments,
        require_reference_lengths=True
    )
    
    # Find 2ia4 evidence
    ia4_evidence = [e for e in evidence if e.source_pdb == "2ia4"]
    print(f"\nFound {len(ia4_evidence)} 2ia4 evidence items:")
    
    for i, ev in enumerate(ia4_evidence):
        print(f"  {i+1}. Type: {ev.type}")
        print(f"      Domain ID: {ev.domain_id}")
        print(f"      Range: {ev.query_range}")
        print(f"      Has alignment: {ev.alignment is not None}")
        if ev.alignment:
            print(f"      Alignment range: query {ev.alignment.query_start}-{ev.alignment.query_end}")
            print(f"                       hit {ev.alignment.hit_start}-{ev.alignment.hit_end}")
        print()
    
    return ia4_evidence, domain_definitions, blast_alignments

def debug_domain_definitions(domain_definitions):
    """Check domain definitions for 2ia4"""
    
    print("\nSTEP 2: CHECKING DOMAIN DEFINITIONS FOR 2ia4")
    print("=" * 50)
    
    # Check for 2ia4 domain definitions
    ia4_keys = [key for key in domain_definitions.keys() if key[0] == "2ia4"]
    print(f"Found domain definition keys for 2ia4: {ia4_keys}")
    
    for key in ia4_keys:
        domains = domain_definitions[key]
        print(f"\n2ia4_{key[1]} has {len(domains)} domain definitions:")
        for d in domains:
            print(f"  {d.domain_id}: {d.range} ({d.length} residues)")

def debug_alignment_data(blast_alignments):
    """Check BLAST alignment data for 2ia4"""
    
    print("\nSTEP 3: CHECKING BLAST ALIGNMENT DATA")
    print("=" * 50)
    
    # Check for 2ia4 alignments
    ia4_alignments = {k: v for k, v in blast_alignments.items() if k[0] == "2ia4"}
    print(f"Found {len(ia4_alignments)} 2ia4 alignments:")
    
    for (pdb, chain), align in ia4_alignments.items():
        print(f"\n{pdb}_{chain}:")
        print(f"  Query range: {align.query_start}-{align.query_end}")
        print(f"  Hit range: {align.hit_start}-{align.hit_end}")
        print(f"  E-value: {align.evalue}")
        print(f"  Query seq length: {len(align.query_seq)}")
        print(f"  Hit seq length: {len(align.hit_seq)}")

def test_decomposition_directly(ia4_evidence, domain_definitions):
    """Test decomposition function directly"""
    
    print("\nSTEP 4: TESTING DECOMPOSITION DIRECTLY")
    print("=" * 50)
    
    from mini.decomposer import decompose_chain_blast_with_mapping
    
    # Find chain BLAST evidence for 2ia4
    chain_blast_evidence = None
    for ev in ia4_evidence:
        if ev.type == "chain_blast":
            chain_blast_evidence = ev
            break
    
    if not chain_blast_evidence:
        print("✗ No chain BLAST evidence found for 2ia4")
        return []
    
    print(f"✓ Found chain BLAST evidence: {chain_blast_evidence.domain_id}")
    print(f"  Range: {chain_blast_evidence.query_range}")
    print(f"  Has alignment: {chain_blast_evidence.alignment is not None}")
    
    # Check for domain definitions
    hit_key = (chain_blast_evidence.source_pdb, chain_blast_evidence.domain_id.split('_')[-1])
    print(f"  Looking for domain definitions with key: {hit_key}")
    
    if hit_key not in domain_definitions:
        print(f"✗ No domain definitions found for key {hit_key}")
        print(f"  Available keys: {[k for k in domain_definitions.keys() if k[0] == '2ia4']}")
        return []
    
    print(f"✓ Found {len(domain_definitions[hit_key])} domain definitions")
    
    # Test decomposition
    if not chain_blast_evidence.alignment:
        print("✗ No alignment data available")
        return []
    
    print("✓ Testing decomposition with alignment data...")
    
    try:
        decomposed = decompose_chain_blast_with_mapping(
            chain_blast_evidence,
            chain_blast_evidence.alignment.query_seq,
            chain_blast_evidence.alignment.hit_seq,
            chain_blast_evidence.alignment.query_start,
            chain_blast_evidence.alignment.hit_start,
            domain_definitions[hit_key],
            verbose=True
        )
        
        print(f"\n✓ Decomposition successful: {len(decomposed)} evidence items")
        for i, d in enumerate(decomposed):
            print(f"  {i+1}. Type: {d.type}")
            print(f"      Domain ID: {d.domain_id}")
            print(f"      Range: {d.query_range}")
            print(f"      Reference length: {d.reference_length}")
        
        return decomposed
        
    except Exception as e:
        print(f"✗ Decomposition failed: {e}")
        import traceback
        traceback.print_exc()
        return []

def test_partitioner_integration(decomposed_evidence):
    """Test if the partitioner correctly processes decomposed evidence"""
    
    print("\nSTEP 5: TESTING PARTITIONER INTEGRATION")
    print("=" * 50)
    
    if not decomposed_evidence:
        print("✗ No decomposed evidence to test")
        return
    
    # Check evidence types
    print("Decomposed evidence types:")
    for d in decomposed_evidence:
        print(f"  {d.type}: {d.domain_id}")
    
    # Verify they have the right type
    decomposed_types = [d.type for d in decomposed_evidence]
    if "chain_blast_decomposed" in decomposed_types:
        print("✓ Found chain_blast_decomposed type")
    else:
        print("✗ Missing chain_blast_decomposed type - this is the bug!")
        print(f"  Actual types: {set(decomposed_types)}")

def main():
    """Main debugging function"""
    
    print("DEBUGGING 2ia4 DECOMPOSITION PIPELINE")
    print("=" * 60)
    print("Expected: 3 domains (GFP + PBP1 + PBP2)")
    print("Actual: 2 domains (GFP + PBP_combined)")
    print("=" * 60)
    
    # Step by step debugging
    ia4_evidence, domain_definitions, blast_alignments = debug_2ia4_evidence()
    debug_domain_definitions(domain_definitions)
    debug_alignment_data(blast_alignments)
    decomposed = test_decomposition_directly(ia4_evidence, domain_definitions)
    test_partitioner_integration(decomposed)
    
    print("\n" + "=" * 60)
    print("DIAGNOSIS:")
    print("=" * 60)
    if len(decomposed) > 1:
        print("✓ Decomposition works in isolation")
        print("✗ Integration with partitioner is the problem")
        print("  → Check partitioner.py for missing decomposition calls")
    else:
        print("✗ Decomposition itself is failing")
        print("  → Check alignment data and domain definitions")

if __name__ == "__main__":
    main()
