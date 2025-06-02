#!/usr/bin/env python3
"""
Check what type of evidence is hitting 6dgv in the 8ovp_A domain summary

This will help us understand why 6dgv is being treated as chain BLAST
when it might should be domain-specific evidence.
"""

import sys
import os
import xml.etree.ElementTree as ET
from pathlib import Path

def analyze_6dgv_evidence(xml_path: str):
    """Analyze all 6dgv-related evidence in the domain summary XML"""
    
    print("Analyzing 6dgv evidence in domain summary XML")
    print("=" * 60)
    
    if not os.path.exists(xml_path):
        print(f"ERROR: XML file not found: {xml_path}")
        return
    
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except Exception as e:
        print(f"ERROR: Failed to parse XML: {e}")
        return
    
    # Check chain BLAST for 6dgv hits
    print("1. CHAIN BLAST HITS:")
    print("-" * 30)
    chain_blast_hits = root.findall(".//chain_blast_run/hits/hit")
    gfp_chain_hits = []
    
    for hit in chain_blast_hits:
        pdb_id = hit.get("pdb_id", "")
        if pdb_id == "6dgv":
            gfp_chain_hits.append(hit)
            chain_id = hit.get("chain_id", "")
            evalues = hit.get("evalues", "")
            query_reg = hit.find("query_reg")
            query_range = query_reg.text if query_reg is not None else "N/A"
            
            print(f"  6dgv_{chain_id}: {query_range} (e-value: {evalues})")
    
    print(f"Found {len(gfp_chain_hits)} chain BLAST hits to 6dgv")
    
    # Check domain BLAST for 6dgv domain hits
    print("\n2. DOMAIN BLAST HITS:")
    print("-" * 30)
    domain_blast_hits = root.findall(".//blast_run/hits/hit")
    gfp_domain_hits = []
    
    for hit in domain_blast_hits:
        domain_id = hit.get("domain_id", "")
        if "6dgv" in domain_id:
            gfp_domain_hits.append(hit)
            evalues = hit.get("evalues", "")
            query_reg = hit.find("query_reg")
            query_range = query_reg.text if query_reg is not None else "N/A"
            
            print(f"  {domain_id}: {query_range} (e-value: {evalues})")
    
    print(f"Found {len(gfp_domain_hits)} domain BLAST hits to 6dgv domains")
    
    # Check HHsearch for 6dgv hits
    print("\n3. HHSEARCH HITS:")
    print("-" * 30)
    hhsearch_hits = root.findall(".//hh_run/hits/hit")
    gfp_hh_hits = []
    
    for hit in hhsearch_hits:
        hit_id = hit.get("hit_id", "")
        domain_id = hit.get("domain_id", "")
        if "6dgv" in hit_id or "6dgv" in domain_id:
            gfp_hh_hits.append(hit)
            probability = hit.get("probability", "")
            query_reg = hit.find("query_reg")
            query_range = query_reg.text if query_reg is not None else "N/A"
            
            print(f"  {hit_id or domain_id}: {query_range} (prob: {probability})")
    
    print(f"Found {len(gfp_hh_hits)} HHsearch hits to 6dgv")
    
    # Analysis
    print("\n" + "=" * 60)
    print("ANALYSIS:")
    print("=" * 60)
    
    total_6dgv_evidence = len(gfp_chain_hits) + len(gfp_domain_hits) + len(gfp_hh_hits)
    print(f"Total 6dgv evidence: {total_6dgv_evidence}")
    
    if len(gfp_chain_hits) > 0:
        print(f"\nChain BLAST hits ({len(gfp_chain_hits)}):")
        print("  - These should use protein reference length (6dgv_A = 559 residues)")
        print("  - This is CORRECT behavior for chain BLAST")
        
    if len(gfp_domain_hits) > 0:
        print(f"\nDomain BLAST hits ({len(gfp_domain_hits)}):")
        print("  - These should use domain reference length (e6dgvA1 = 252 residues)")
        
    if len(gfp_hh_hits) > 0:
        print(f"\nHHsearch hits ({len(gfp_hh_hits)}):")
        print("  - These should use domain reference length (e6dgvA1 = 252 residues)")
    
    # Check which evidence is being selected
    print(f"\nThe algorithm reported 'Source: chain_blast' for the 6dgv domain.")
    print(f"This suggests that:")
    if len(gfp_chain_hits) > 0 and len(gfp_domain_hits) == 0 and len(gfp_hh_hits) == 0:
        print("  ✓ Only chain BLAST evidence exists - using protein length is CORRECT")
        print("  ✓ The 559 reference length is appropriate for chain BLAST")
        print("  ℹ️  This is actually the correct behavior!")
    elif len(gfp_chain_hits) > 0:
        print("  ⚠️  Chain BLAST was selected over domain-specific evidence")
        print("  ⚠️  Check evidence precedence/quality scoring")
    else:
        print("  ❌ Something is wrong with evidence detection")

def check_precedence_logic():
    """Explain the evidence precedence and why chain BLAST might be selected"""
    
    print("\n" + "=" * 60)
    print("EVIDENCE PRECEDENCE EXPLANATION:")
    print("=" * 60)
    
    print("mini_pyecod uses this precedence order:")
    print("1. chain_blast (precedence 0)")
    print("2. chain_blast_decomposed (precedence 0)")  
    print("3. domain_blast (precedence 1)")
    print("4. hhsearch (precedence 2)")
    print()
    print("If 6dgv has both chain BLAST and domain BLAST hits,")
    print("chain BLAST will be selected first due to higher precedence.")
    print()
    print("This means:")
    print("✓ Using protein length (559) for chain BLAST is CORRECT")
    print("✓ The algorithm is working as designed")
    print("ℹ️  The 44.2% coverage might be misleading but technically accurate")
    print()
    print("For comparison:")
    print("- Chain BLAST coverage: 247/559 = 44.2% (of full protein)")
    print("- Domain coverage: 247/252 = 98.0% (of GFP domain)")
    print("Both are correct for their respective reference types!")

if __name__ == "__main__":
    xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"
    analyze_6dgv_evidence(xml_path)
    check_precedence_logic()
