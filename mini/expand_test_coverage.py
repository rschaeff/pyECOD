#!/usr/bin/env python3
"""
Expand test coverage for mini_pyecod

Now that we've validated the core algorithm works with 8ovp_A,
let's systematically test different domain architectures and edge cases.
"""

import sys
import os
from pathlib import Path
import random
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class TestCase:
    """Test case definition"""
    protein_id: str
    expected_domains: int
    category: str
    description: str
    priority: int = 1  # 1=high, 2=medium, 3=low

def find_available_test_cases(batch_dir: str) -> List[str]:
    """Find all available protein IDs in the batch"""
    domains_dir = os.path.join(batch_dir, "domains")
    
    available = []
    if os.path.exists(domains_dir):
        for filename in os.listdir(domains_dir):
            if filename.endswith(".develop291.domain_summary.xml"):
                protein_id = filename.replace(".develop291.domain_summary.xml", "")
                available.append(protein_id)
    
    return sorted(available)

def create_systematic_test_suite(available_proteins: List[str]) -> List[TestCase]:
    """Create a systematic test suite covering different architectures"""
    
    # Known cases (if available)
    known_cases = [
        TestCase("8ovp_A", 2, "fusion", "GFP-PBP fusion (validated)", 1),
        TestCase("1ubq_A", 1, "single", "Ubiquitin - classic single domain", 1),
        TestCase("2trx_A", 1, "single", "Thioredoxin - single domain", 2),
        TestCase("1cuk_A", 2, "multi", "Two-domain kinase", 1),
        TestCase("2src_A", 2, "multi", "SH3 + kinase domain", 2),
    ]
    
    # Filter to only available cases
    test_suite = []
    for case in known_cases:
        if case.protein_id in available_proteins:
            test_suite.append(case)
    
    # Add systematic sampling of available proteins
    random.seed(42)  # Reproducible
    remaining = [p for p in available_proteins if p not in [c.protein_id for c in test_suite]]
    
    # Sample by PDB characteristics
    short_ids = [p for p in remaining if len(p) == 6 and '_' in p]  # Standard PDB_Chain format
    long_ids = [p for p in remaining if len(p) > 6]  # Complex naming
    
    # Add diverse samples
    if short_ids:
        sample_short = random.sample(short_ids, min(10, len(short_ids)))
        for protein_id in sample_short:
            test_suite.append(TestCase(protein_id, 0, "sample", f"Standard format sample", 2))
    
    if long_ids:
        sample_long = random.sample(long_ids, min(5, len(long_ids)))
        for protein_id in sample_long:
            test_suite.append(TestCase(protein_id, 0, "complex", f"Complex naming sample", 3))
    
    return test_suite

def run_test_case(test_case: TestCase, batch_dir: str, cache_file: str, verbose: bool = False) -> dict:
    """Run a single test case and return results"""
    
    print(f"\n{'='*60}")
    print(f"Testing {test_case.protein_id}: {test_case.description}")
    print(f"Category: {test_case.category}, Priority: {test_case.priority}")
    print(f"{'='*60}")
    
    # Import mini components
    sys.path.insert(0, str(Path(__file__).parent.parent))
    
    try:
        from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
        from mini.blast_parser import load_chain_blast_alignments
        from mini.decomposer import load_domain_definitions
        from mini.partitioner import partition_domains
        from mini.writer import write_domain_partition
    except ImportError as e:
        return {"status": "import_error", "error": str(e)}
    
    # Reference files (use range cache versions)
    ref_dir = "test_data_v291"
    domain_lengths_file = os.path.join(ref_dir, 'domain_lengths.csv')
    protein_lengths_file = os.path.join(ref_dir, 'protein_lengths.csv')
    domain_definitions_file = os.path.join(ref_dir, 'domain_definitions.csv')
    
    # Check if reference files exist
    if not all(os.path.exists(f) for f in [domain_lengths_file, protein_lengths_file]):
        return {"status": "missing_references", "error": "Run range cache setup first"}
    
    # Load reference data
    try:
        reference_lengths = load_reference_lengths(domain_lengths_file)
        protein_lengths = load_protein_lengths(protein_lengths_file)
        domain_definitions = load_domain_definitions(domain_definitions_file) if os.path.exists(domain_definitions_file) else {}
    except Exception as e:
        return {"status": "reference_error", "error": str(e)}
    
    # Parse protein ID
    parts = test_case.protein_id.split('_')
    if len(parts) < 2:
        return {"status": "invalid_id", "error": f"Cannot parse protein ID: {test_case.protein_id}"}
    
    pdb_id = parts[0]
    chain_id = parts[1]
    
    # Check if domain summary exists
    xml_path = os.path.join(batch_dir, "domains", f"{test_case.protein_id}.develop291.domain_summary.xml")
    if not os.path.exists(xml_path):
        return {"status": "missing_xml", "error": f"Domain summary not found"}
    
    try:
        # Load BLAST alignments
        blast_dir = os.path.join(batch_dir, "blast/chain")
        blast_alignments = {}
        if os.path.exists(blast_dir):
            blast_alignments = load_chain_blast_alignments(blast_dir, pdb_id, chain_id)
        
        # Parse evidence
        evidence = parse_domain_summary(
            xml_path,
            reference_lengths=reference_lengths,
            protein_lengths=protein_lengths,
            blast_alignments=blast_alignments,
            require_reference_lengths=True,
            verbose=verbose
        )
        
        if len(evidence) == 0:
            return {"status": "no_evidence", "error": "No evidence with reference lengths"}
        
        # Estimate sequence length
        max_pos = 0
        for ev in evidence:
            for segment in ev.query_range.segments:
                max_pos = max(max_pos, segment.end)
        sequence_length = int(max_pos * 1.1)
        
        # Partition domains
        domains = partition_domains(evidence, sequence_length=sequence_length, verbose=verbose)
        
        # Analyze results
        result = {
            "status": "success",
            "protein_id": test_case.protein_id,
            "category": test_case.category,
            "expected_domains": test_case.expected_domains,
            "found_domains": len(domains),
            "evidence_count": len(evidence),
            "sequence_length": sequence_length,
            "domains": []
        }
        
        # Domain details
        total_coverage = 0
        for domain in domains:
            domain_info = {
                "family": domain.family,
                "range": str(domain.range),
                "size": domain.range.total_length,
                "discontinuous": domain.range.is_discontinuous,
                "source": domain.source
            }
            result["domains"].append(domain_info)
            total_coverage += domain.range.total_length
        
        result["coverage"] = total_coverage / sequence_length if sequence_length > 0 else 0
        
        # Write output
        output_path = f"/tmp/test_coverage/{test_case.protein_id}.domains.xml"
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_domain_partition(domains, pdb_id, chain_id, output_path)
        
        # Print summary
        print(f"Results: {len(domains)} domains, {len(evidence)} evidence, {result['coverage']:.1%} coverage")
        for i, domain in enumerate(domains, 1):
            discontinuous = " (discontinuous)" if domain.range.is_discontinuous else ""
            print(f"  {i}. {domain.family}: {domain.range}{discontinuous}")
        
        return result
        
    except Exception as e:
        import traceback
        return {"status": "error", "error": str(e), "traceback": traceback.format_exc()}

def run_test_suite(test_cases: List[TestCase], batch_dir: str, cache_file: str, max_tests: int = 20) -> dict:
    """Run a suite of test cases"""
    
    print("MINI_PYECOD EXPANDED TEST COVERAGE")
    print("=" * 60)
    print(f"Running up to {max_tests} test cases")
    
    # Sort by priority and limit
    priority_cases = sorted(test_cases, key=lambda x: (x.priority, x.protein_id))[:max_tests]
    
    results = []
    categories = {}
    
    for i, test_case in enumerate(priority_cases, 1):
        print(f"\n[{i}/{len(priority_cases)}]", end=" ")
        
        result = run_test_case(test_case, batch_dir, cache_file)
        results.append(result)
        
        # Track by category
        category = test_case.category
        if category not in categories:
            categories[category] = {"total": 0, "success": 0, "domains": []}
        
        categories[category]["total"] += 1
        if result["status"] == "success":
            categories[category]["success"] += 1
            categories[category]["domains"].append(result["found_domains"])
    
    # Summary report
    print("\n" + "=" * 60)
    print("TEST SUITE SUMMARY")
    print("=" * 60)
    
    total_tests = len(results)
    successful = sum(1 for r in results if r["status"] == "success")
    
    print(f"Total tests: {total_tests}")
    print(f"Successful: {successful}")
    print(f"Failed: {total_tests - successful}")
    
    # Category breakdown
    print(f"\nBy category:")
    for category, stats in categories.items():
        success_rate = stats["success"] / stats["total"] if stats["total"] > 0 else 0
        avg_domains = sum(stats["domains"]) / len(stats["domains"]) if stats["domains"] else 0
        print(f"  {category}: {stats['success']}/{stats['total']} successful ({success_rate:.1%}), avg {avg_domains:.1f} domains")
    
    # Failure analysis
    failures = [r for r in results if r["status"] != "success"]
    if failures:
        print(f"\nFailure types:")
        failure_types = {}
        for f in failures:
            ftype = f["status"]
            failure_types[ftype] = failure_types.get(ftype, 0) + 1
        
        for ftype, count in sorted(failure_types.items()):
            print(f"  {ftype}: {count}")
    
    return {
        "results": results,
        "summary": {
            "total": total_tests,
            "successful": successful,
            "categories": categories
        }
    }

def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Expand mini_pyecod test coverage')
    parser.add_argument('--batch-dir', 
                        default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                        help='Batch directory')
    parser.add_argument('--cache-file',
                        default='/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt',
                        help='Range cache file')
    parser.add_argument('--max-tests', type=int, default=15,
                        help='Maximum number of tests to run')
    parser.add_argument('--setup-only', action='store_true',
                        help='Only set up reference files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    
    args = parser.parse_args()
    
    # Set up reference files if needed
    if not os.path.exists("test_data_v291/domain_lengths.csv") or args.setup_only:
        print("Setting up reference files from range cache...")
        try:
            sys.path.insert(0, str(Path(__file__).parent))
            from range_cache_parser import (
                create_domain_lengths_from_cache,
                create_domain_definitions_from_cache,
                extract_protein_lengths_from_cache
            )
            
            os.makedirs("test_data_v291", exist_ok=True)
            create_domain_lengths_from_cache(args.cache_file, "test_data_v291/domain_lengths.csv")
            create_domain_definitions_from_cache(args.cache_file, "test_data_v291/domain_definitions.csv")
            extract_protein_lengths_from_cache(args.cache_file, "test_data_v291/protein_lengths.csv")
            
            if args.setup_only:
                return
                
        except Exception as e:
            print(f"ERROR: Failed to set up reference files: {e}")
            return
    
    # Find available test cases
    available_proteins = find_available_test_cases(args.batch_dir)
    print(f"Found {len(available_proteins)} available proteins in batch")
    
    # Create test suite
    test_suite = create_systematic_test_suite(available_proteins)
    print(f"Created test suite with {len(test_suite)} cases")
    
    # Run tests
    results = run_test_suite(test_suite, args.batch_dir, args.cache_file, args.max_tests)
    
    print(f"\n✅ Test coverage expansion complete!")
    print(f"Detailed results in /tmp/test_coverage/")

if __name__ == "__main__":
    main()
