#!/usr/bin/env python3
"""Debug XML structure for chain BLAST hits"""

import xml.etree.ElementTree as ET
import sys

def debug_chain_blast_xml(xml_path: str):
    """Examine the structure of chain BLAST hits in the XML"""
    
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except Exception as e:
        print(f"Error parsing XML: {e}")
        return
    
    # Find chain BLAST hits
    chain_blast_hits = root.findall(".//chain_blast_run/hits/hit")
    print(f"Found {len(chain_blast_hits)} chain BLAST hits")
    
    if chain_blast_hits:
        # Examine first few hits
        for i, hit in enumerate(chain_blast_hits[:3]):
            print(f"\n--- Chain BLAST Hit {i+1} ---")
            print(f"PDB: {hit.get('pdb_id')}")
            print(f"Chain: {hit.get('chain_id')}")
            
            # List all child elements
            print("Child elements:")
            for child in hit:
                text = child.text[:50] + "..." if child.text and len(child.text) > 50 else child.text
                print(f"  <{child.tag}>: {text}")
            
            # Check for alignment-related elements
            query_seq = hit.find("query_seq")
            hit_seq = hit.find("hit_seq")
            query_aln = hit.find("query_aln")
            hit_aln = hit.find("hit_aln")
            
            print("\nAlignment elements:")
            print(f"  query_seq: {'Found' if query_seq is not None else 'Not found'}")
            print(f"  hit_seq: {'Found' if hit_seq is not None else 'Not found'}")
            print(f"  query_aln: {'Found' if query_aln is not None else 'Not found'}")
            print(f"  hit_aln: {'Found' if hit_aln is not None else 'Not found'}")
            
            # Try other possible names
            for elem_name in ['qseq', 'hseq', 'query_alignment', 'hit_alignment', 'alignment']:
                elem = hit.find(elem_name)
                if elem is not None:
                    print(f"  {elem_name}: Found!")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        xml_path = sys.argv[1]
    else:
        xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    
    debug_chain_blast_xml(xml_path)
