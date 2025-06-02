#!/usr/bin/env python3
"""Test mini_pyecod on proteins with known/expected domain counts"""

import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from quick_test import run_test

# Known test cases with expected domain counts
KNOWN_CASES = [
    # (protein_id, expected_domains, description)
    ('8ovp_A', 3, 'GFP + two PBP domains fusion'),
    ('9nlb_E', 2, 'Two domain protein from batch'),
    ('8opd_Aa', 2, 'Two domain protein with chain BLAST'),
    # Add more as you validate them
]

def test_known_cases(verbose=False):
    """Test proteins with known domain counts"""
    
    print("="*60)
    print("Testing proteins with known domain counts")
    print("="*60)
    
    results = []
    
    for protein_id, expected, description in KNOWN_CASES:
        print(f"\n{'='*60}")
        print(f"Test: {protein_id} - {description}")
        print(f"Expected: {expected} domains")
        print(f"{'='*60}")
        
        domains = run_test(protein_id, verbose=verbose, reference_lengths_file=None)
        
        if domains is not None:
            actual = len(domains)
            status = "✓ PASS" if actual == expected else "✗ FAIL"
            results.append({
                'protein': protein_id,
                'expected': expected,
                'actual': actual,
                'status': status,
                'description': description
            })
            
            print(f"\n{'='*40}")
            print(f"Result: {status} - Found {actual} domains (expected {expected})")
            print(f"{'='*40}")
        else:
            results.append({
                'protein': protein_id,
                'expected': expected,
                'actual': 'ERROR',
                'status': '✗ ERROR',
                'description': description
            })
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"{'Protein':<12} {'Expected':<10} {'Actual':<10} {'Status':<10} Description")
    print("-"*60)
    
    passed = 0
    for r in results:
        print(f"{r['protein']:<12} {r['expected']:<10} {str(r['actual']):<10} {r['status']:<10} {r['description']}")
        if r['status'] == '✓ PASS':
            passed += 1
    
    print("-"*60)
    print(f"Passed: {passed}/{len(results)}")
    
    # Recommendations
    if passed < len(results):
        print("\nRecommendations for failed cases:")
        print("1. Check evidence quality thresholds")
        print("2. Verify sequence length estimation")
        print("3. Review rejection statistics in detailed output")
        print("4. Consider adjusting NEW_COVERAGE_THRESHOLD or OLD_COVERAGE_THRESHOLD")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Test mini_pyecod on proteins with known domain counts')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    test_known_cases(verbose=args.verbose)
