#!/usr/bin/env python3
"""Test simple cases that we verified exist"""

import os
import sys

def test_verified_cases():
    """Test cases we know exist"""
    
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    domains_dir = os.path.join(batch_dir, "domains")
    
    # Cases we know exist
    verified_cases = [
        ("8ovp_A", "GFP-PBP fusion - our working example"),
        ("8opd_Aa", "Two domain protein"),
    ]
    
    print("Testing verified cases:")
    print("=" * 60)
    
    for protein, description in verified_cases:
        print(f"\n{protein}: {description}")
        
        # Check exact filenames
        possible_files = [
            f"{protein}.develop291.domain_summary.xml",
            f"{protein}.domain_summary.xml",
        ]
        
        found = False
        for fname in possible_files:
            fpath = os.path.join(domains_dir, fname)
            if os.path.exists(fpath):
                print(f"  ✓ Found: {fname}")
                found = True
                break
        
        if found:
            print(f"  Run: python quick_test.py {protein}")
        else:
            print(f"  ✗ Not found")
    
    # Let's also check the exact structure of filenames
    print("\n" + "="*60)
    print("Sample of actual filenames in directory:")
    print("="*60)
    
    all_files = os.listdir(domains_dir)
    domain_files = [f for f in all_files if "domain_summary" in f]
    
    # Show first 10
    for f in sorted(domain_files)[:10]:
        print(f"  {f}")
    
    # Extract pattern
    print("\n" + "="*60)
    print("Filename patterns:")
    print("="*60)
    
    patterns = set()
    for f in domain_files[:100]:
        if f.endswith(".domain_summary.xml"):
            # Extract the pattern
            base = f.replace(".domain_summary.xml", "")
            if ".develop291" in base:
                patterns.add("NAME.develop291.domain_summary.xml")
            else:
                patterns.add("NAME.domain_summary.xml")
    
    for p in patterns:
        print(f"  {p}")
    
    # Find some simple test cases
    print("\n" + "="*60)
    print("Looking for simple test cases...")
    print("="*60)
    
    # Try to find proteins with simple names
    simple_proteins = []
    for f in domain_files[:200]:
        if f.endswith(".develop291.domain_summary.xml"):
            protein = f.replace(".develop291.domain_summary.xml", "")
            # Look for simple names (4 char PDB + chain)
            if len(protein) >= 5 and len(protein) <= 6 and '_' in protein:
                simple_proteins.append(protein)
    
    if simple_proteins:
        print(f"Found {len(simple_proteins)} proteins with standard names:")
        for p in simple_proteins[:10]:
            print(f"  python quick_test.py {p}")

if __name__ == "__main__":
    test_verified_cases()
