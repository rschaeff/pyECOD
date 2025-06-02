"""
Recommended test cases for mini_pyecod testing

These test cases cover different domain architectures and challenges:

1. SINGLE DOMAIN PROTEINS (baseline tests)
   - 1ubq_A: Ubiquitin - small, single domain
   - 2trx_A: Thioredoxin - single domain with clear boundaries
   - 1tim_A: TIM barrel - single domain, but large

2. TWO-DOMAIN PROTEINS (simple multi-domain)
   - 1cuk_A: Glycogen phosphorylase - two clear domains
   - 2src_A: Src kinase - SH3 + kinase domain
   - 1lba_A: Lactoferrin - two lobes

3. DISCONTINUOUS DOMAINS (challenging)
   - 8ovp_A: GFP-PBP fusion (already tested)
   - 1xkw_A: Ferredoxin-like fold with insertion
   - 3hhm_A: HUP domain with discontinuous topology

4. MANY DOMAINS (complex architectures)
   - 1fat_A: Focal adhesion kinase - 4 domains
   - 2v9v_A: Multi-domain signaling protein
   - 1vps_A: Clathrin heavy chain - many repeats

5. DOMAIN REPEATS (tandem repeats)
   - 1a0s_P: Ankyrin repeats
   - 1qgk_A: TPR repeats
   - 3lbx_A: WD40 repeats

6. DIFFICULT CASES
   - Small domains (<50 residues)
   - Intrinsically disordered regions
   - Domain interfaces that are ambiguous

To add these test cases:

1. Find them in your batch:
   ls /data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/ | grep -E '1ubq_A|2trx_A|1tim_A|1cuk_A|2src_A'

2. Generate reference data:
   python prepare_test_data.py 1ubq_A 2trx_A 1tim_A 1cuk_A 2src_A

3. Run the tests:
   python test_framework.py 1ubq_A
   python test_framework.py 1cuk_A
   etc.

4. Or add them to TEST_CASES in test_framework.py for batch testing
"""

# Script to check which recommended cases are available
import os
import sys

def check_available_test_cases(batch_dir):
    """Check which recommended test cases are available in the batch"""
    
    recommended = {
        "Single domain": ["1ubq_A", "2trx_A", "1tim_A", "1mbn_A", "2lyz_A"],
        "Two domain": ["1cuk_A", "2src_A", "1lba_A", "1gca_A", "3pgk_A"],
        "Discontinuous": ["8ovp_A", "1xkw_A", "3hhm_A", "2acy_A", "1mol_A"],
        "Many domains": ["1fat_A", "2v9v_A", "1vps_A", "1xwd_A", "2wss_A"],
        "Repeats": ["1a0s_P", "1qgk_A", "3lbx_A", "1bd8_A", "2bnh_A"]
    }
    
    domains_dir = os.path.join(batch_dir, "domains")
    available_files = set(os.listdir(domains_dir)) if os.path.exists(domains_dir) else set()
    
    print("Checking recommended test cases in batch...")
    print("=" * 60)
    
    found_cases = []
    
    for category, cases in recommended.items():
        print(f"\n{category}:")
        for case in cases:
            filename = f"{case}.develop291.domain_summary.xml"
            if filename in available_files:
                print(f"  ✓ {case} - Available")
                found_cases.append(case)
            else:
                print(f"  ✗ {case} - Not in batch")
    
    if found_cases:
        print(f"\n\nFound {len(found_cases)} test cases. Generate reference data with:")
        print(f"python prepare_test_data.py {' '.join(found_cases)}")

if __name__ == "__main__":
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    check_available_test_cases(batch_dir)
