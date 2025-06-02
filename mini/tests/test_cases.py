#!/usr/bin/env python3
"""
Formal test cases for mini_pyecod domain partitioning

This module defines the official test cases that validate
the domain partitioning algorithm works correctly.
"""

import os
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.blast_parser import load_chain_blast_alignments
from mini.decomposer import load_domain_definitions
from mini.partitioner import partition_domains

@dataclass
class ExpectedDomain:
    """Expected domain in test case"""
    family: str
    approximate_range: str
    min_size: int
    max_size: int
    discontinuous: bool = False
    notes: str = ""

@dataclass
class TestCase:
    """Formal test case definition"""
    protein_id: str
    description: str
    expected_domain_count: int
    expected_domains: List[ExpectedDomain]
    requires_decomposition: bool = False
    requires_blast_alignments: bool = False
    notes: str = ""

# Official test cases
TEST_CASES = [
    TestCase(
        protein_id="8ovp_A",
        description="GFP-PBP fusion with chain BLAST decomposition",
        expected_domain_count=3,
        expected_domains=[
            ExpectedDomain(
                family="6dgv",  # GFP
                approximate_range="252-494",
                min_size=240,
                max_size=250,
                discontinuous=False,
                notes="GFP domain should be continuous"
            ),
            ExpectedDomain(
                family="2ia4",  # PBP domain 1
                approximate_range="~100-200",
                min_size=90,
                max_size=110,
                discontinuous=False,
                notes="First PBP domain from decomposition"
            ),
            ExpectedDomain(
                family="2ia4",  # PBP domain 2 (discontinuous)
                approximate_range="~1-100,200-250,490-520",
                min_size=150,
                max_size=200,
                discontinuous=True,
                notes="Second PBP domain, should be discontinuous"
            )
        ],
        requires_decomposition=True,
        requires_blast_alignments=True,
        notes="This is the canonical test case - fusion protein with complex domain insertion"
    ),
    
    # Additional test cases can be added here as they are validated
    # TestCase(
    #     protein_id="1ubq_A",
    #     description="Single domain protein (ubiquitin)",
    #     expected_domain_count=1,
    #     expected_domains=[
    #         ExpectedDomain(
    #             family="ubiquitin",
    #             approximate_range="1-76",
    #             min_size=70,
    #             max_size=80,
    #             discontinuous=False
    #         )
    #     ],
    #     requires_decomposition=False
    # ),
]

class TestRunner:
    """Run formal test cases"""
    
    def __init__(self, batch_dir: str, test_data_dir: str = None):
        self.batch_dir = batch_dir
        self.test_data_dir = test_data_dir or str(Path(__file__).parent.parent / "test_data")
        
        # Load reference data once
        self.reference_lengths = self._load_reference_data("domain_lengths.csv")
        self.protein_lengths = self._load_reference_data("protein_lengths.csv", is_protein=True)
        self.domain_definitions = self._load_reference_data("domain_definitions.csv", is_definitions=True)
        
        print(f"Loaded reference data:")
        print(f"  Domain lengths: {len(self.reference_lengths)}")
        print(f"  Protein lengths: {len(self.protein_lengths)}")
        print(f"  Domain definitions: {len(self.domain_definitions)} proteins")
    
    def _load_reference_data(self, filename: str, is_protein: bool = False, is_definitions: bool = False):
        """Load reference data files"""
        filepath = os.path.join(self.test_data_dir, filename)
        
        if not os.path.exists(filepath):
            print(f"Warning: Reference file not found: {filepath}")
            return {}
        
        if is_definitions:
            from mini.decomposer import load_domain_definitions
            return load_domain_definitions(filepath)
        elif is_protein:
            return load_protein_lengths(filepath)
        else:
            return load_reference_lengths(filepath)
    
    def run_test_case(self, test_case: TestCase, verbose: bool = False) -> Dict:
        """Run a single test case"""
        
        print(f"\n{'='*60}")
        print(f"Testing: {test_case.protein_id}")
        print(f"Description: {test_case.description}")
        print(f"Expected: {test_case.expected_domain_count} domains")
        print(f"{'='*60}")
        
        # Parse protein ID
        parts = test_case.protein_id.split('_')
        pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else 'A'
        
        # Check required files
        xml_path = os.path.join(self.batch_dir, "domains", f"{test_case.protein_id}.develop291.domain_summary.xml")
        if not os.path.exists(xml_path):
            return {
                'status': 'missing_file',
                'error': f"Domain summary not found: {xml_path}",
                'test_case': test_case
            }
        
        try:
            # Load BLAST alignments if required
            blast_alignments = {}
            if test_case.requires_blast_alignments:
                blast_dir = os.path.join(self.batch_dir, "blast/chain")
                blast_alignments = load_chain_blast_alignments(blast_dir, pdb_id, chain_id, verbose=verbose)
                if verbose:
                    print(f"Loaded {len(blast_alignments)} BLAST alignments")
            
            # Parse evidence
            evidence = parse_domain_summary(
                xml_path,
                reference_lengths=self.reference_lengths,
                protein_lengths=self.protein_lengths,
                blast_alignments=blast_alignments,
                require_reference_lengths=True,
                verbose=verbose
            )
            
            if not evidence:
                return {
                    'status': 'no_evidence',
                    'error': 'No evidence with reference lengths found',
                    'test_case': test_case
                }
            
            # Estimate sequence length
            max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
            sequence_length = int(max_pos * 1.1)
            
            # Partition domains
            domain_definitions = self.domain_definitions if test_case.requires_decomposition else None
            domains = partition_domains(
                evidence, 
                sequence_length=sequence_length,
                domain_definitions=domain_definitions,
                verbose=verbose
            )
            
            # Analyze results
            result = self._analyze_results(test_case, domains, verbose)
            
            return result
            
        except Exception as e:
            import traceback
            return {
                'status': 'error',
                'error': str(e),
                'traceback': traceback.format_exc(),
                'test_case': test_case
            }
    
    def _analyze_results(self, test_case: TestCase, domains: List, verbose: bool = False) -> Dict:
        """Analyze test results against expectations"""
        
        result = {
            'status': 'completed',
            'test_case': test_case,
            'found_domains': len(domains),
            'expected_domains': test_case.expected_domain_count,
            'domains': [],
            'checks': {}
        }
        
        # Basic domain count check
        result['checks']['domain_count'] = len(domains) == test_case.expected_domain_count
        
        # Analyze each found domain
        for domain in domains:
            domain_info = {
                'id': domain.id,
                'family': domain.family,
                'range': str(domain.range),
                'size': domain.range.total_length,
                'discontinuous': domain.range.is_discontinuous,
                'source': domain.source
            }
            result['domains'].append(domain_info)
        
        # Check expected domain families
        found_families = [d.family for d in domains]
        expected_families = [ed.family for ed in test_case.expected_domains]
        result['checks']['families_found'] = all(family in found_families for family in expected_families)
        
        # Check discontinuous domains if expected
        has_discontinuous = any(d.range.is_discontinuous for d in domains)
        expects_discontinuous = any(ed.discontinuous for ed in test_case.expected_domains)
        result['checks']['discontinuous'] = has_discontinuous == expects_discontinuous
        
        # Overall pass/fail
        result['passed'] = all(result['checks'].values())
        
        # Print results
        print(f"\nResults:")
        print(f"  Found {len(domains)} domains (expected {test_case.expected_domain_count})")
        
        for i, domain in enumerate(domains, 1):
            disc = " (discontinuous)" if domain.range.is_discontinuous else ""
            print(f"  {i}. {domain.family}: {domain.range} ({domain.range.total_length} residues){disc}")
        
        print(f"\nChecks:")
        for check, passed in result['checks'].items():
            status = "✅" if passed else "❌"
            print(f"  {status} {check.replace('_', ' ').title()}")
        
        overall_status = "✅ PASSED" if result['passed'] else "❌ FAILED"
        print(f"\n{overall_status}")
        
        return result
    
    def run_all_tests(self, verbose: bool = False) -> Dict:
        """Run all test cases"""
        
        print(f"Running {len(TEST_CASES)} formal test cases")
        print("=" * 60)
        
        results = []
        passed = 0
        
        for test_case in TEST_CASES:
            result = self.run_test_case(test_case, verbose)
            results.append(result)
            
            if result.get('passed', False):
                passed += 1
        
        # Summary
        print(f"\n{'='*60}")
        print(f"SUMMARY: {passed}/{len(TEST_CASES)} tests passed")
        print(f"{'='*60}")
        
        return {
            'total': len(TEST_CASES),
            'passed': passed,
            'failed': len(TEST_CASES) - passed,
            'results': results
        }

def main():
    """Command line interface for test runner"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Run formal mini_pyecod test cases')
    parser.add_argument('--batch-dir',
                        default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                        help='Batch directory containing domain files')
    parser.add_argument('--test-data-dir',
                        help='Test data directory (default: mini/test_data)')
    parser.add_argument('--protein', help='Run specific protein test case')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    
    args = parser.parse_args()
    
    # Initialize test runner
    runner = TestRunner(args.batch_dir, args.test_data_dir)
    
    if args.protein:
        # Run specific test case
        test_case = None
        for tc in TEST_CASES:
            if tc.protein_id == args.protein:
                test_case = tc
                break
        
        if test_case:
            result = runner.run_test_case(test_case, args.verbose)
            if not result.get('passed', False):
                sys.exit(1)
        else:
            print(f"Test case not found: {args.protein}")
            print(f"Available test cases: {[tc.protein_id for tc in TEST_CASES]}")
            sys.exit(1)
    else:
        # Run all test cases
        summary = runner.run_all_tests(args.verbose)
        if summary['failed'] > 0:
            sys.exit(1)

if __name__ == "__main__":
    main()
