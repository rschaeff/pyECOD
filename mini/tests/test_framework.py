#!/usr/bin/env python3
"""Test framework for mini_pyecod with multiple test cases"""

import sys
import os
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.blast_parser import load_chain_blast_alignments
from mini.decomposer import load_domain_definitions
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition

@dataclass
class TestCase:
    """Test case definition"""
    pdb_id: str
    chain_id: str
    description: str
    expected_domains: int
    expected_structure: List[Dict[str, any]]
    notes: Optional[str] = None

# Define test cases covering different architectures
TEST_CASES = [
    TestCase(
        pdb_id="8ovp",
        chain_id="A", 
        description="GFP-PBP fusion with domain insertion",
        expected_domains=3,
        expected_structure=[
            {"name": "GFP", "range": "252-494", "discontinuous": False},
            {"name": "PBP1", "range": "108-205", "discontinuous": False},
            {"name": "PBP2", "range": "4-107,206-251,495-512", "discontinuous": True}
        ],
        notes="C-terminal segment may be missing due to alignment limitations"
    ),
    
    # Add more test cases here:
    # TestCase(
    #     pdb_id="1cuk",
    #     chain_id="A",
    #     description="Two-domain kinase",
    #     expected_domains=2,
    #     expected_structure=[
    #         {"name": "N-terminal", "range": "1-150", "discontinuous": False},
    #         {"name": "C-terminal", "range": "151-300", "discontinuous": False}
    #     ]
    # ),
]

def run_test_case(test_case: TestCase, batch_dir: str) -> Dict[str, any]:
    """Run a single test case"""
    
    print(f"\n{'='*60}")
    print(f"Testing {test_case.pdb_id}_{test_case.chain_id}: {test_case.description}")
    print(f"{'='*60}")
    
    # Paths
    xml_path = f"{batch_dir}/domains/{test_case.pdb_id}_{test_case.chain_id}.develop291.domain_summary.xml"
    blast_dir = f"{batch_dir}/blast/chain"
    
    if not os.path.exists(xml_path):
        print(f"❌ Domain summary not found: {xml_path}")
        return {"status": "missing_file", "test_case": test_case}
    
    # Load reference data
    ref_lengths = load_reference_lengths("test_data/domain_lengths.csv") if os.path.exists("test_data/domain_lengths.csv") else {}
    protein_lengths = load_protein_lengths("test_data/protein_lengths.csv") if os.path.exists("test_data/protein_lengths.csv") else {}
    domain_defs = load_domain_definitions("test_data/domain_definitions.csv") if os.path.exists("test_data/domain_definitions.csv") else {}
    
    # Load BLAST alignments
    blast_alignments = {}
    if os.path.exists(blast_dir):
        blast_alignments = load_chain_blast_alignments(blast_dir, test_case.pdb_id, test_case.chain_id)
        print(f"Loaded {len(blast_alignments)} BLAST alignments")
    
    # Parse evidence
    evidence = parse_domain_summary(xml_path, ref_lengths, protein_lengths, blast_alignments)
    print(f"Found {len(evidence)} total evidence items")
    
    # Filter to high-quality evidence
    good_evidence = [e for e in evidence if e.confidence > 0.6 or (e.evalue and e.evalue < 1e-5)]
    print(f"Filtered to {len(good_evidence)} high-quality evidence items")
    
    # Get sequence length (approximate from evidence)
    max_pos = 1
    for ev in evidence:
        for pos, _ in ev.query_range.to_positions():
            max_pos = max(max_pos, pos)
    sequence_length = max_pos + 50  # Add some buffer
    
    # Partition domains
    domains = partition_domains(good_evidence, sequence_length=sequence_length,
                               domain_defs=domain_defs, use_precedence=True)
    
    # Analyze results
    result = {
        "status": "completed",
        "test_case": test_case,
        "found_domains": len(domains),
        "expected_domains": test_case.expected_domains,
        "domains": domains,
        "passed": len(domains) == test_case.expected_domains
    }
    
    # Display results
    print(f"\nResults:")
    print(f"  Expected: {test_case.expected_domains} domains")
    print(f"  Found: {len(domains)} domains")
    
    for i, domain in enumerate(domains, 1):
        discontinuous = "discontinuous" if domain.range.is_discontinuous else "continuous"
        print(f"  {i}. {domain.family}: {domain.range} ({domain.range.total_length} residues, {discontinuous})")
    
    if test_case.notes:
        print(f"\nNotes: {test_case.notes}")
    
    if result["passed"]:
        print("\n✅ PASSED")
    else:
        print("\n❌ FAILED")
    
    # Write output
    output_path = f"/tmp/{test_case.pdb_id}_{test_case.chain_id}_mini.domains.xml"
    write_domain_partition(domains, test_case.pdb_id, test_case.chain_id, output_path)
    
    return result

def run_all_tests(batch_dir: str):
    """Run all test cases"""
    
    results = []
    passed = 0
    failed = 0
    
    print(f"Running {len(TEST_CASES)} test cases...")
    
    for test_case in TEST_CASES:
        result = run_test_case(test_case, batch_dir)
        results.append(result)
        
        if result["status"] == "completed":
            if result["passed"]:
                passed += 1
            else:
                failed += 1
    
    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY: {passed} passed, {failed} failed")
    print(f"{'='*60}")
    
    # Details of failures
    if failed > 0:
        print("\nFailed cases:")
        for result in results:
            if result["status"] == "completed" and not result["passed"]:
                tc = result["test_case"]
                print(f"  - {tc.pdb_id}_{tc.chain_id}: expected {tc.expected_domains}, found {result['found_domains']}")
    
    return results

if __name__ == "__main__":
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    
    if len(sys.argv) > 1:
        # Run specific test case
        pdb_chain = sys.argv[1]
        if '_' in pdb_chain:
            pdb_id, chain_id = pdb_chain.split('_')
            test_case = TestCase(
                pdb_id=pdb_id,
                chain_id=chain_id,
                description="Manual test case",
                expected_domains=0,  # Unknown
                expected_structure=[]
            )
            run_test_case(test_case, batch_dir)
        else:
            print(f"Usage: {sys.argv[0]} [PDB_CHAIN]")
            print(f"Example: {sys.argv[0]} 1cuk_A")
    else:
        # Run all defined test cases
        run_all_tests(batch_dir)
