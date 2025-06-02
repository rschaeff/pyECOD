#!/usr/bin/env python3
"""
Validate ECOD range cache for mini_pyecod test cases

This script checks that the range cache contains the domains we need
for testing, particularly for the 8ovp_A case with 2ia4 PBP domains.
"""

import sys
import os
from pathlib import Path

def check_range_cache_structure(cache_file: str):
    """Check the basic structure of the range cache file"""
    
    print(f"Checking range cache structure: {cache_file}")
    print("=" * 60)
    
    if not os.path.exists(cache_file):
        print(f"❌ File not found: {cache_file}")
        return False
    
    try:
        with open(cache_file, 'r') as f:
            lines = []
            for i, line in enumerate(f):
                lines.append(line.strip())
                if i >= 10:  # Just read first 10 lines
                    break
        
        print("First 10 lines:")
        for i, line in enumerate(lines):
            parts = line.split('\t')
            print(f"  {i+1}: {len(parts)} columns - {line[:80]}...")
        
        # Count total lines
        with open(cache_file, 'r') as f:
            total_lines = sum(1 for _ in f)
        
        print(f"\nTotal lines: {total_lines:,}")
        return True
        
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return False

def check_test_domains(cache_file: str):
    """Check for specific domains needed for testing"""
    
    print(f"\nChecking for test domains in range cache...")
    print("=" * 60)
    
    # Domains we expect from the grep output
    expected_2ia4_domains = [
        'e2ia4A1',  # A:110-209
        'e2ia4A2',  # A:3-109,A:210-279  (discontinuous!)
        'e2ia4B1',  # B:2-128
        'e2ia4B2'   # B:129-279
    ]
    
    # Other important test domains
    other_test_domains = [
        'e6dgvA1',  # GFP domain for 8ovp
        'e1gflA1',  # Alternative GFP
    ]
    
    found_domains = {}
    
    try:
        with open(cache_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) == 5:
                    cache_id, domain_id, range_spec, pdb_id, chain_id = parts
                    
                    # Check if this is a domain we're looking for
                    if domain_id in expected_2ia4_domains + other_test_domains:
                        found_domains[domain_id] = {
                            'range_spec': range_spec,
                            'pdb_id': pdb_id,
                            'chain_id': chain_id,
                            'line_num': line_num
                        }
        
        # Report findings
        print("2ia4 PBP domains (critical for 8ovp_A test):")
        for domain_id in expected_2ia4_domains:
            if domain_id in found_domains:
                entry = found_domains[domain_id]
                print(f"  ✅ {domain_id}: {entry['range_spec']} (line {entry['line_num']})")
                
                # Check if discontinuous
                if ',' in entry['range_spec']:
                    print(f"     🔀 DISCONTINUOUS - This is important for decomposition!")
            else:
                print(f"  ❌ {domain_id}: NOT FOUND")
        
        print(f"\nOther test domains:")
        for domain_id in other_test_domains:
            if domain_id in found_domains:
                entry = found_domains[domain_id]
                print(f"  ✅ {domain_id}: {entry['range_spec']} (line {entry['line_num']})")
            else:
                print(f"  ❌ {domain_id}: NOT FOUND")
        
        # Summary
        found_2ia4 = sum(1 for d in expected_2ia4_domains if d in found_domains)
        print(f"\nSummary: Found {found_2ia4}/{len(expected_2ia4_domains)} required 2ia4 domains")
        
        return found_2ia4 == len(expected_2ia4_domains)
        
    except Exception as e:
        print(f"❌ Error checking domains: {e}")
        return False

def check_domain_statistics(cache_file: str):
    """Get some statistics about domains in the cache"""
    
    print(f"\nDomain statistics from range cache...")
    print("=" * 60)
    
    pdb_counts = {}
    discontinuous_count = 0
    total_domains = 0
    length_sum = 0
    
    try:
        with open(cache_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) == 5:
                    cache_id, domain_id, range_spec, pdb_id, chain_id = parts
                    
                    total_domains += 1
                    
                    # Count by PDB
                    pdb_counts[pdb_id] = pdb_counts.get(pdb_id, 0) + 1
                    
                    # Check if discontinuous
                    if ',' in range_spec:
                        discontinuous_count += 1
                    
                    # Estimate length (rough)
                    try:
                        # Simple length estimation from range
                        segments = range_spec.split(',')
                        domain_length = 0
                        for seg in segments:
                            if ':' in seg:
                                range_part = seg.split(':', 1)[1]
                            else:
                                range_part = seg
                            
                            if '-' in range_part:
                                start, end = range_part.split('-')
                                domain_length += int(end) - int(start) + 1
                        
                        length_sum += domain_length
                    except:
                        pass  # Skip if can't parse
        
        print(f"Total domains: {total_domains:,}")
        print(f"Discontinuous domains: {discontinuous_count:,} ({discontinuous_count/total_domains:.1%})")
        print(f"Average domain length: {length_sum/total_domains:.1f} residues")
        print(f"Unique PDB structures: {len(pdb_counts):,}")
        
        # Show top PDBs by domain count
        top_pdbs = sorted(pdb_counts.items(), key=lambda x: -x[1])[:10]
        print(f"\nTop 10 PDBs by domain count:")
        for pdb_id, count in top_pdbs:
            print(f"  {pdb_id}: {count} domains")
        
        return True
        
    except Exception as e:
        print(f"❌ Error getting statistics: {e}")
        return False

def main():
    """Main validation function"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Validate ECOD range cache for mini_pyecod testing'
    )
    parser.add_argument('--cache-file',
                        default='/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt',
                        help='Path to ECOD range cache file')
    parser.add_argument('--full-stats', action='store_true',
                        help='Show full statistics (slower)')
    
    args = parser.parse_args()
    
    print("ECOD Range Cache Validation for mini_pyecod")
    print("=" * 60)
    
    # Check file structure
    if not check_range_cache_structure(args.cache_file):
        print("\n❌ Range cache validation FAILED")
        sys.exit(1)
    
    # Check for test domains
    test_domains_ok = check_test_domains(args.cache_file)
    
    # Optionally show statistics
    if args.full_stats:
        check_domain_statistics(args.cache_file)
    
    # Final verdict
    print("\n" + "=" * 60)
    if test_domains_ok:
        print("✅ VALIDATION PASSED: Range cache contains required test domains")
        print("\nNext steps:")
        print("1. Run: python range_cache_parser.py --cache-file YOUR_CACHE_FILE")
        print("2. Test: python test_with_range_cache.py 8ovp_A")
    else:
        print("❌ VALIDATION FAILED: Some required test domains are missing")
        print("\nThis may indicate:")
        print("1. Wrong cache file version")
        print("2. Cache file is corrupted")
        print("3. Domain IDs have changed between versions")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
