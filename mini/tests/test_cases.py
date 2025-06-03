#!/usr/bin/env python3
"""
Official integration test cases for mini_pyecod

This module contains the formal test cases that validate the complete
domain partitioning pipeline works correctly on real proteins.

8ovp_A is the primary/canonical test case that all CI/CD should run against.
"""

import pytest
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional

from mini.parser import parse_domain_summary
from mini.partitioner import partition_domains
from mini.writer import write_domain_partition


@dataclass
class ExpectedDomain:
    """Expected domain characteristics for validation"""
    family: str
    approximate_range: str
    min_size: int
    max_size: int
    discontinuous: bool = False
    notes: str = ""


@dataclass
class TestCase:
    """Official test case definition"""
    protein_id: str
    description: str
    expected_domain_count: int
    expected_domains: List[ExpectedDomain]
    requires_decomposition: bool = False
    requires_blast_alignments: bool = False
    notes: str = ""


# Official test cases - 8ovp_A is the canonical/primary test
OFFICIAL_TEST_CASES = {
    "8ovp_A": TestCase(
        protein_id="8ovp_A",
        description="GFP-PBP fusion protein with domain insertion (PRIMARY TEST CASE)",
        expected_domain_count=3,
        expected_domains=[
            ExpectedDomain(
                family="6dgv",  # GFP domain
                approximate_range="252-494",
                min_size=240,
                max_size=250,
                discontinuous=False,
                notes="GFP domain should be continuous and well-defined"
            ),
            ExpectedDomain(
                family="e2ia4A1",  # PBP domain 1 (from decomposition)
                approximate_range="~108-205",
                min_size=90,
                max_size=110,
                discontinuous=False,
                notes="First PBP domain from chain BLAST decomposition"
            ),
            ExpectedDomain(
                family="e2ia4A2",  # PBP domain 2 (discontinuous from decomposition)
                approximate_range="~2-107,206-251,491-517",
                min_size=150,
                max_size=200,
                discontinuous=True,
                notes="Second PBP domain, should be discontinuous after GFP insertion"
            )
        ],
        requires_decomposition=True,
        requires_blast_alignments=True,
        notes="""
This is the canonical test case for mini_pyecod. It tests:
- Chain BLAST decomposition with alignment data
- Domain insertion architecture (GFP inserted into PBP)
- Discontinuous domain handling
- Multi-family partitioning

This test MUST pass for any release.
        """
    ),

    # Additional validated test cases can be added here as they are confirmed
    # "1ubq_A": TestCase(...),  # Single domain protein
    # "1cuk_A": TestCase(...),  # Two-domain protein
}


class TestOfficialCases:
    """Official integration test cases for mini_pyecod"""

    @pytest.mark.integration
    @pytest.mark.slow
    def test_8ovp_A_canonical(self, stable_batch_dir, real_reference_data,
                             blast_alignments, temp_output_dir):
        """
        PRIMARY TEST CASE: 8ovp_A with full decomposition

        This is the gold standard test that validates the complete pipeline:
        - Evidence parsing from real XML
        - Chain BLAST decomposition with alignment data
        - Proper domain family assignment
        - Discontinuous domain handling
        """
        test_case = OFFICIAL_TEST_CASES["8ovp_A"]
        result = self._run_test_case(
            test_case,
            stable_batch_dir,
            real_reference_data,
            blast_alignments,
            temp_output_dir
        )

        # Detailed validation for primary test case
        assert result['passed'], f"Primary test case failed: {result.get('error', 'Unknown error')}"
        assert result['found_domains'] == 3, f"Expected 3 domains, found {result['found_domains']}"

        # Validate specific domain characteristics
        domains = result['domains']

        # Should have GFP domain
        gfp_domains = [d for d in domains if '6dgv' in d['family']]
        assert len(gfp_domains) == 1, f"Expected 1 GFP domain, found {len(gfp_domains)}"

        # Should have PBP domains from decomposition
        pbp_domains = [d for d in domains if '2ia4' in d['family']]
        assert len(pbp_domains) == 2, f"Expected 2 PBP domains from decomposition, found {len(pbp_domains)}"

        # Should have one discontinuous domain
        discontinuous_domains = [d for d in domains if d['discontinuous']]
        assert len(discontinuous_domains) >= 1, "Expected at least one discontinuous domain"

        print(f"✅ PRIMARY TEST CASE PASSED: {test_case.protein_id}")
        for i, domain in enumerate(domains, 1):
            disc_note = " (discontinuous)" if domain['discontinuous'] else ""
            print(f"   {i}. {domain['family']}: {domain['range']} ({domain['size']} residues){disc_note}")

    @pytest.mark.integration
    def test_8ovp_A_without_decomposition(self, stable_batch_dir, real_reference_data, temp_output_dir):
        """
        Test 8ovp_A without chain BLAST decomposition

        This tests the basic partitioning algorithm without decomposition.
        Should produce 2 domains instead of 3.
        """
        test_case = OFFICIAL_TEST_CASES["8ovp_A"]

        # Run without domain definitions (disables decomposition)
        reference_data_no_decomp = real_reference_data.copy()
        reference_data_no_decomp['domain_definitions'] = {}

        result = self._run_test_case(
            test_case,
            stable_batch_dir,
            reference_data_no_decomp,
            blast_alignments={},  # No BLAST alignments
            output_dir=temp_output_dir,
            expect_decomposition=False
        )

        assert result['passed'], f"No-decomposition test failed: {result.get('error', 'Unknown error')}"

        # Without decomposition, should get fewer domains (GFP + single PBP chain)
        assert result['found_domains'] <= 2, f"Without decomposition, expected ≤2 domains, found {result['found_domains']}"

        print(f"✅ NO-DECOMPOSITION TEST PASSED: Found {result['found_domains']} domains (expected ≤2)")

    @pytest.mark.performance
    @pytest.mark.integration
    def test_8ovp_A_performance(self, stable_batch_dir, real_reference_data, blast_alignments, temp_output_dir):
        """
        Performance benchmark for 8ovp_A

        Ensures processing completes within reasonable time limits.
        """
        import time

        test_case = OFFICIAL_TEST_CASES["8ovp_A"]

        start_time = time.time()
        result = self._run_test_case(
            test_case,
            stable_batch_dir,
            real_reference_data,
            blast_alignments,
            temp_output_dir
        )
        processing_time = time.time() - start_time

        # Performance requirements
        MAX_PROCESSING_TIME = 60  # seconds

        assert result['passed'], "Performance test must use working algorithm"
        assert processing_time < MAX_PROCESSING_TIME, f"Processing took {processing_time:.1f}s (max: {MAX_PROCESSING_TIME}s)"

        print(f"✅ PERFORMANCE TEST PASSED: {processing_time:.2f}s (limit: {MAX_PROCESSING_TIME}s)")

    @pytest.mark.parametrize("protein_id", ["8ovp_A"])  # Extend as more cases are validated
    def test_output_file_generation(self, protein_id, stable_batch_dir, real_reference_data,
                                   blast_alignments, temp_output_dir):
        """
        Test that output XML files are generated correctly
        """
        test_case = OFFICIAL_TEST_CASES[protein_id]
        result = self._run_test_case(
            test_case,
            stable_batch_dir,
            real_reference_data,
            blast_alignments,
            temp_output_dir
        )

        assert result['passed'], f"Test case failed: {result.get('error', 'Unknown error')}"

        # Check output file was created
        output_file = Path(temp_output_dir) / f"{protein_id}_test.domains.xml"
        assert output_file.exists(), f"Output file not created: {output_file}"

        # Validate XML structure
        import xml.etree.ElementTree as ET
        tree = ET.parse(output_file)
        root = tree.getroot()

        assert root.tag == "domain_partition", "Root element should be domain_partition"
        assert root.get("pdb_id") == protein_id.split('_')[0], "PDB ID should be set correctly"

        domains_elem = root.find("domains")
        assert domains_elem is not None, "Should have domains element"

        domain_elems = domains_elem.findall("domain")
        assert len(domain_elems) == result['found_domains'], "XML should match found domain count"

        print(f"✅ OUTPUT VALIDATION PASSED: {output_file}")

    def _run_test_case(self, test_case: TestCase, batch_dir: str,
                      reference_data: Dict, blast_alignments: Dict,
                      output_dir: str, expect_decomposition: bool = True) -> Dict:
        """
        Run a single test case and return detailed results
        """
        import os

        # Parse protein ID
        parts = test_case.protein_id.split('_')
        pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else 'A'

        # Check domain summary file
        xml_path = os.path.join(batch_dir, "domains", f"{test_case.protein_id}.develop291.domain_summary.xml")
        if not os.path.exists(xml_path):
            return {
                'passed': False,
                'error': f"Domain summary not found: {xml_path}",
                'test_case': test_case
            }

        try:
            # Parse evidence
            evidence = parse_domain_summary(
                xml_path,
                reference_lengths=reference_data.get('domain_lengths', {}),
                protein_lengths=reference_data.get('protein_lengths', {}),
                blast_alignments=blast_alignments,
                require_reference_lengths=True,
                verbose=False
            )

            if not evidence:
                return {
                    'passed': False,
                    'error': 'No evidence with reference lengths found',
                    'test_case': test_case
                }

            # Estimate sequence length
            max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
            sequence_length = int(max_pos * 1.1)

            # Partition domains
            domain_definitions = reference_data.get('domain_definitions', {}) if expect_decomposition else {}
            domains = partition_domains(
                evidence,
                sequence_length=sequence_length,
                domain_definitions=domain_definitions,
                verbose=False
            )

            # Write output for validation
            output_file = os.path.join(output_dir, f"{test_case.protein_id}_test.domains.xml")
            write_domain_partition(domains, pdb_id, chain_id, output_file)

            # Analyze results
            result = self._analyze_results(test_case, domains)
            result['output_file'] = output_file

            return result

        except Exception as e:
            import traceback
            return {
                'passed': False,
                'error': str(e),
                'traceback': traceback.format_exc(),
                'test_case': test_case
            }

    def _analyze_results(self, test_case: TestCase, domains: List) -> Dict:
        """
        Analyze test results against expectations
        """
        result = {
            'passed': False,
            'test_case': test_case,
            'found_domains': len(domains),
            'expected_domains': test_case.expected_domain_count,
            'domains': [],
            'checks': {}
        }

        # Convert domains to analysis format
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

        # Run validation checks
        result['checks']['domain_count'] = len(domains) == test_case.expected_domain_count

        # Check that expected families are found
        found_families = [d.family for d in domains]
        expected_families = [ed.family for ed in test_case.expected_domains]

        # For decomposed domains, check partial matches (e.g., "2ia4" in "e2ia4A1")
        families_found = []
        for expected in expected_families:
            for found in found_families:
                if expected in found or found in expected:
                    families_found.append(expected)
                    break

        result['checks']['families_found'] = len(families_found) >= len(set(expected_families))

        # Check discontinuous domain expectation
        has_discontinuous = any(d.range.is_discontinuous for d in domains)
        expects_discontinuous = any(ed.discontinuous for ed in test_case.expected_domains)
        result['checks']['discontinuous'] = has_discontinuous == expects_discontinuous

        # Check domain sizes are reasonable
        size_checks = []
        for domain in domains:
            # Check that domain size is reasonable (not too small/large)
            size_ok = 20 <= domain.range.total_length <= 1000
            size_checks.append(size_ok)
        result['checks']['reasonable_sizes'] = all(size_checks)

        # Overall pass/fail
        result['passed'] = all(result['checks'].values())

        return result


class TestCaseValidation:
    """Validation utilities for test cases"""

    @pytest.mark.unit
    def test_test_case_definitions(self):
        """
        Validate that test case definitions are well-formed
        """
        for protein_id, test_case in OFFICIAL_TEST_CASES.items():
            # Basic validation
            assert test_case.protein_id == protein_id, f"Protein ID mismatch for {protein_id}"
            assert test_case.expected_domain_count > 0, f"Invalid domain count for {protein_id}"
            assert len(test_case.expected_domains) > 0, f"No expected domains defined for {protein_id}"

            # Validate expected domains
            for expected_domain in test_case.expected_domains:
                assert expected_domain.family, f"Missing family for domain in {protein_id}"
                assert expected_domain.min_size > 0, f"Invalid min_size for domain in {protein_id}"
                assert expected_domain.max_size >= expected_domain.min_size, f"Invalid size range for domain in {protein_id}"

        print(f"✅ TEST CASE DEFINITIONS VALIDATED: {len(OFFICIAL_TEST_CASES)} cases")

    @pytest.mark.unit
    def test_primary_test_case_properties(self):
        """
        Validate that 8ovp_A has all required properties for a primary test case
        """
        primary = OFFICIAL_TEST_CASES["8ovp_A"]

        # Primary test case requirements
        assert primary.requires_decomposition, "Primary test case must test decomposition"
        assert primary.requires_blast_alignments, "Primary test case must test BLAST alignments"
        assert primary.expected_domain_count >= 3, "Primary test case should be complex (≥3 domains)"
        assert any(ed.discontinuous for ed in primary.expected_domains), "Primary test case should test discontinuous domains"

        print("✅ PRIMARY TEST CASE PROPERTIES VALIDATED")


# Standalone runner for development/debugging
def main():
    """Run test cases directly (for development)"""
    import argparse

    parser = argparse.ArgumentParser(description='Run official mini_pyecod test cases')
    parser.add_argument('--protein', default='8ovp_A', help='Protein to test')
    parser.add_argument('--batch-dir', help='Batch directory override')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    if args.protein not in OFFICIAL_TEST_CASES:
        print(f"Unknown test case: {args.protein}")
        print(f"Available: {list(OFFICIAL_TEST_CASES.keys())}")
        return 1

    # This would run the test case directly (implementation depends on environment setup)
    print(f"Running test case: {args.protein}")
    print("Use 'pytest tests/test_cases.py::TestOfficialCases::test_8ovp_A_canonical -v' for full pytest integration")

    return 0


if __name__ == "__main__":
    exit(main())
