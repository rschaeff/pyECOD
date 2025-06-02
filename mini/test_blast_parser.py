#!/usr/bin/env python3
"""Test BLAST XML parser"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.blast_parser import parse_blast_xml, load_chain_blast_alignments

def test_blast_parser():
    """Test parsing BLAST XML for alignment data"""
    
    # Test with 8ovp_A BLAST results
    blast_file = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/blast/chain/8ovp_A.develop291.xml"
    
    print(f"Parsing BLAST XML: {blast_file}")
    alignments = parse_blast_xml(blast_file)
    
    print(f"\nFound {len(alignments)} alignments:")
    
    # Show first few alignments
    for i, ((pdb, chain), align) in enumerate(sorted(alignments.items())[:5]):
        print(f"\n{i+1}. {pdb}_{chain}:")
        print(f"   E-value: {align.evalue}")
        print(f"   Query range: {align.query_start}-{align.query_end}")
        print(f"   Hit range: {align.hit_start}-{align.hit_end}")
        print(f"   Query seq length: {len(align.query_seq)}")
        print(f"   Hit seq length: {len(align.hit_seq)}")
        
        # Show alignment preview
        if align.query_seq and align.hit_seq:
            print(f"   Query: {align.query_seq[:50]}...")
            print(f"   Hit:   {align.hit_seq[:50]}...")
    
    # Check specifically for 2ia4
    if ('2ia4', 'B') in alignments:
        print(f"\n2ia4_B alignment found!")
        align = alignments[('2ia4', 'B')]
        print(f"   Query positions: {align.query_start}-{align.query_end}")
        print(f"   This should help decompose the discontinuous range")
    
    # Test the loader function
    print("\n\nTesting load_chain_blast_alignments...")
    blast_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/blast/chain"
    alignments2 = load_chain_blast_alignments(blast_dir, "8ovp", "A")
    print(f"Loaded {len(alignments2)} alignments using loader function")

if __name__ == "__main__":
    test_blast_parser()
