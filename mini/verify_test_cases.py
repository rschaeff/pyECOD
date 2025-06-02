#!/usr/bin/env python3
"""Verify which suggested test files actually exist"""

import os

def verify_test_files(batch_dir):
    """Check which test files actually exist"""
    
    # Test cases suggested by find_test_cases.py
    suggested_cases = [
        "8p7d_A",    # Single domain
        "8opl_Bp",   # Single domain
        "8opf_An",   # Single domain  
        "8oz8_D",    # 2 domains
        "8p4r_G",    # 3 domains
        "8olc_C",    # 2 domains
        "8p4r_a",    # 3 domains
        "8oll_Z",    # 8 domains
        "8p9a_C",    # 20 domains
        "8p0k_O",    # 6 domains
        "8ovp_A",    # Our known working case
    ]
    
    domains_dir = os.path.join(batch_dir, "domains")
    blast_dir = os.path.join(batch_dir, "blast/chain")
    
    print("Verifying test file availability...")
    print("=" * 70)
    print(f"{'Protein':<10} {'Domain Summary':<20} {'Chain BLAST':<20} {'Status':<20}")
    print("=" * 70)
    
    available_cases = []
    
    for case in suggested_cases:
        domain_file = os.path.join(domains_dir, f"{case}.develop291.domain_summary.xml")
        blast_file = os.path.join(blast_dir, f"{case}.develop291.xml")
        
        has_domain = os.path.exists(domain_file)
        has_blast = os.path.exists(blast_file)
        
        status = "Ready" if has_domain else "Missing domain summary"
        if has_domain:
            available_cases.append(case)
            
        print(f"{case:<10} {'✓' if has_domain else '✗':<20} {'✓' if has_blast else '✗':<20} {status:<20}")
    
    print("\n" + "=" * 70)
    print(f"Available for testing: {len(available_cases)} out of {len(suggested_cases)} cases")
    
    if available_cases:
        print("\nReady to test:")
        for case in available_cases[:5]:  # Show first 5
            print(f"  python quick_test.py {case}")
    
    # Also check for other files in the directory
    print("\n\nChecking for other available proteins...")
    all_files = os.listdir(domains_dir) if os.path.exists(domains_dir) else []
    domain_summary_files = [f for f in all_files if f.endswith(".domain_summary.xml")]
    
    print(f"Total domain summary files in batch: {len(domain_summary_files)}")
    
    # Sample a few random ones
    if len(domain_summary_files) > 10:
        print("\nRandom sample of available proteins:")
        import random
        random.seed(42)  # For reproducibility
        sample = random.sample(domain_summary_files, 10)
        for f in sample:
            protein = f.replace(".develop291.domain_summary.xml", "")
            print(f"  {protein}")
    
    return available_cases

if __name__ == "__main__":
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    verify_test_files(batch_dir)
