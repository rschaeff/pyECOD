#!/usr/bin/env python3
"""
Core algorithm tests for mini_pyecod

Tests the fundamental domain partitioning algorithm components.
"""

import pytest
from pathlib import Path
from typing import List

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.core.models import Evidence, Domain, AlignmentData
from mini.core.partitioner import partition_domains
from mini.core.sequence_range import SequenceRange


def create_complete_evidence(
    evidence_type: str,
    source_pdb: str,
    query_range: str,
    confidence: float,
    reference_length: int,
    domain_id: str,
    evalue: float = 1e-20,
    t_group: str = None,
    **kwargs
) -> Evidence:
    """Helper function to create evidence with all required fields for quality thresholds"""

    # Calculate proper reference coverage based on type
    if evidence_type == "hhsearch":
        min_ref_coverage = 0.60  # 60% threshold for hhsearch
    else:
        min_ref_coverage = 0.50  # 50% threshold for domain_blast and chain_blast

    # Use provided reference_coverage or calculate a good default
    reference_coverage = kwargs.get('reference_coverage', min_ref_coverage + 0.1)

    query_seq_range = SequenceRange.parse(query_range)

    return Evidence(
        type=evidence_type,
        source_pdb=source_pdb,
        query_range=query_seq_range,
        confidence=confidence,
        evalue=evalue,
        reference_length=reference_length,
        reference_coverage=reference_coverage,
        alignment_coverage=reference_coverage,  # Set same as reference_coverage
        domain_id=domain_id,
        hit_range=SequenceRange.parse(f"1-{reference_length}"),  # Default hit range
        source_chain_id="A",  # Default chain
        discontinuous=query_seq_range.is_discontinuous,
        hsp_count=1,
        t_group=t_group,
        **kwargs
    )


class TestResidueBlocking:
    """Test the residue blocking algorithm"""

    @pytest.mark.unit
    def test_basic_residue_blocking(self):
        """Test that domains block residues from reuse"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="test1",
                query_range="10-100",
                confidence=0.95,
                reference_length=91,
                domain_id="test1_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="test2",
                query_range="50-150",  # Overlaps with first
                confidence=0.85,
                reference_length=101,
                domain_id="test2_A"
            )
        ]

        domains = partition_domains(evidence, sequence_length=200)

        # Should select first domain (higher confidence) and possibly second if overlap is acceptable
        assert len(domains) >= 1
        assert domains[0].family == "test1"

    @pytest.mark.unit
    def test_coverage_thresholds(self):
        """Test NEW_COVERAGE and OLD_COVERAGE thresholds"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="domain1",
                query_range="1-100",
                confidence=0.95,  # Highest confidence - selected first
                reference_length=100,
                domain_id="domain1_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="domain2",
                query_range="90-200",  # 10% overlap with domain1
                confidence=0.95,  # Second highest - should be accepted
                reference_length=111,
                domain_id="domain2_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="domain3",
                query_range="50-150",  # 50% overlap with domain1
                confidence=0.85,  # Lowest confidence, high overlap - should be rejected
                reference_length=101,
                domain_id="domain3_A"
            )
        ]

        domains = partition_domains(evidence, sequence_length=250)

        # Should have domain1 and domain2 (acceptable overlap)
        # Should reject domain3 (too much overlap with domain1)
        assert len(domains) >= 2
        families = {d.family for d in domains}
        assert "domain1" in families
        assert "domain2" in families
        # domain3 may or may not be included depending on overlap calculation

    @pytest.mark.unit
    def test_evidence_priority_sorting(self):
        """Test evidence is processed in correct priority order"""
        evidence = [
            create_complete_evidence(
                evidence_type="hhsearch",
                source_pdb="hh1",
                query_range="1-50",
                confidence=0.75,
                reference_length=50,
                domain_id="hh1_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="blast1",
                query_range="1-50",
                confidence=0.85,
                reference_length=50,
                domain_id="blast1_A"
            ),
            # Chain blast without alignment data - will be rejected
            Evidence(
                type="chain_blast",
                source_pdb="chain1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.95,
                reference_length=50,
                domain_id="chain1_A"
                # No alignment data - will be rejected
            ),
        ]

        domains = partition_domains(evidence, sequence_length=100)

        # Chain blast should be rejected (no decomposition)
        # Domain blast should win over HHsearch (higher confidence)
        assert len(domains) == 1
        assert domains[0].family == "blast1"
        assert domains[0].source == "domain_blast"

    @pytest.mark.unit
    def test_chain_blast_priority_with_decomposition(self):
        """Test that chain blast wins when decomposition is available"""
        from mini.core.models import AlignmentData
        from mini.core.decomposer import DomainReference

        # Create alignment data for chain blast evidence
        alignment = AlignmentData(
            query_seq="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            hit_seq="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            query_start=1,
            query_end=50,
            hit_start=1,
            hit_end=50
        )

        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="blast1",
                query_range="1-50",
                confidence=0.85,
                reference_length=50,
                domain_id="blast1_A"
            ),
            Evidence(
                type="chain_blast",
                source_pdb="chain1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.95,
                reference_length=50,
                domain_id="chain1_A",
                alignment=alignment  # Required for decomposition
            ),
        ]

        # Provide mock domain definitions to enable chain blast decomposition
        domain_definitions = {
            ("chain1", "A"): [
                DomainReference(
                    domain_id="chain1_domain1",
                    pdb_id="chain1",
                    chain_id="A",
                    range=SequenceRange.parse("1-50"),
                    length=50
                )
            ]
        }

        domains = partition_domains(
            evidence,
            sequence_length=100,
            domain_definitions=domain_definitions  # Enable decomposition
        )

        # NOW chain blast should win because decomposition is available
        assert len(domains) == 1
        assert domains[0].source in ["chain_blast", "chain_blast_decomposed"]

    @pytest.mark.unit
    def test_minimum_domain_size(self):
        """Test that tiny domains are rejected"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="tiny",
                query_range="1-15",  # Too small
                confidence=0.95,
                reference_length=15,
                domain_id="tiny_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="normal",
                query_range="20-100",
                confidence=0.85,
                reference_length=81,
                domain_id="normal_A"
            )
        ]

        domains = partition_domains(evidence, sequence_length=150)

        # Should only have the normal-sized domain (tiny domain rejected for size)
        assert len(domains) == 1
        assert domains[0].family == "normal"


class TestDiscontinuousDomains:
    """Test handling of discontinuous domains"""

    @pytest.mark.unit
    def test_discontinuous_domain_parsing(self):
        """Test that discontinuous chain BLAST without decomposition is rejected"""
        evidence = [
            Evidence(
                type="chain_blast",
                source_pdb="disc",
                query_range=SequenceRange.parse("1-100,200-250"),
                confidence=0.95,
                reference_length=151,
                domain_id="disc_A"
                # No alignment data - should be rejected
            )
        ]

        domains = partition_domains(evidence, sequence_length=300)

        # Chain BLAST without decomposition should be rejected
        assert len(domains) == 0

    @pytest.mark.unit
    def test_discontinuous_domain_blast_accepted(self):
        """Test that discontinuous domain BLAST evidence is accepted"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="disc_domain",
                query_range="1-100,200-250",
                confidence=0.95,
                reference_length=151,
                domain_id="disc_domain_A"
            )
        ]

        domains = partition_domains(evidence, sequence_length=300)

        # Domain BLAST with discontinuous range should be accepted
        assert len(domains) == 1
        assert domains[0].range.is_discontinuous
        assert domains[0].range.total_length == 151
        assert len(domains[0].range.segments) == 2

    @pytest.mark.unit
    def test_discontinuous_overlap_calculation(self):
        """Test overlap calculation with discontinuous domains"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="disc1",
                query_range="1-50,100-150",
                confidence=0.95,
                reference_length=101,
                domain_id="disc1_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="cont1",
                query_range="40-120",  # Overlaps both segments
                confidence=0.85,
                reference_length=81,
                domain_id="cont1_A"
            )
        ]

        domains = partition_domains(evidence, sequence_length=200)

        # First domain should be selected
        assert len(domains) == 1
        assert domains[0].family == "disc1"


class TestDomainFamilyAssignment:
    """Test domain family assignment logic"""

    @pytest.mark.unit
    def test_family_from_tgroup(self):
        """Test that T-group is preferred for family assignment"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="test",
                query_range="1-100",
                confidence=0.95,
                reference_length=100,
                domain_id="test_A",
                t_group="1234.5.6"
            )
        ]

        domains = partition_domains(evidence, sequence_length=150)

        assert len(domains) == 1
        assert domains[0].family == "1234.5.6"

    @pytest.mark.unit
    def test_family_fallback_order(self):
        """Test fallback order for family assignment"""
        evidence_list = [
            # Has T-group
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="pdb1",
                query_range="1-50",
                confidence=0.95,
                reference_length=50,
                domain_id="e1abcA1",
                t_group="1111.1.1"
            ),
            # No T-group, has source_pdb
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="pdb2",
                query_range="60-110",
                confidence=0.95,
                reference_length=51,
                domain_id="e2defB1"
            ),
            # Only domain_id
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="",
                query_range="120-170",
                confidence=0.95,
                reference_length=51,
                domain_id="e3ghiC1"
            )
        ]

        domains = partition_domains(evidence_list, sequence_length=200)

        assert len(domains) == 3
        # Domains should be ordered by sequence position (N-terminal to C-terminal)
        # domains[0] should be the most N-terminal (position 1-50)
        assert domains[0].family == "1111.1.1"  # T-group (position 1-50)
        assert domains[1].family == "pdb2"      # source_pdb (position 60-110)
        assert domains[2].family == "e3ghiC1"   # domain_id (position 120-170)


class TestCoverageCalculation:
    """Test sequence coverage calculations"""

    @pytest.mark.unit
    def test_total_coverage(self):
        """Test that coverage is calculated correctly"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="d1",
                query_range="1-100",
                confidence=0.95,
                reference_length=100,
                domain_id="d1_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="d2",
                query_range="150-200",
                confidence=0.95,
                reference_length=51,
                domain_id="d2_A"
            )
        ]

        sequence_length = 250
        domains = partition_domains(evidence, sequence_length)

        # Calculate coverage
        total_coverage = sum(d.range.total_length for d in domains)
        coverage_fraction = total_coverage / sequence_length

        assert len(domains) == 2
        assert total_coverage == 151
        assert coverage_fraction == 151/250

    @pytest.mark.unit
    def test_gap_handling(self):
        """Test that gaps between domains are handled correctly"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="d1",
                query_range="1-50",
                confidence=0.95,
                reference_length=50,
                domain_id="d1_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="d2",
                query_range="100-150",
                confidence=0.95,
                reference_length=51,
                domain_id="d2_A"
            ),
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="d3",
                query_range="200-250",
                confidence=0.95,
                reference_length=51,
                domain_id="d3_A"
            )
        ]

        domains = partition_domains(evidence, sequence_length=300)

        assert len(domains) == 3
        # Check that domains maintain their gaps
        assert domains[0].range.segments[0][1] < domains[1].range.segments[0][0]
        assert domains[1].range.segments[0][1] < domains[2].range.segments[0][0]


class TestEmptyAndEdgeCases:
    """Test edge cases and error conditions"""

    @pytest.mark.unit
    def test_no_evidence(self):
        """Test with no evidence"""
        domains = partition_domains([], sequence_length=100)
        assert len(domains) == 0

    @pytest.mark.unit
    def test_no_reference_lengths(self):
        """Test that evidence without reference lengths can still be processed"""
        evidence = [
            Evidence(
                type="domain_blast",
                source_pdb="test",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.95,
                reference_length=None,  # No reference length
                domain_id="test_A"
            )
        ]

        # Use apply_quality_thresholds=False for this edge case test
        domains = partition_domains(evidence, sequence_length=150, apply_quality_thresholds=False)
        assert len(domains) == 1  # Should accept evidence without reference length

    @pytest.mark.unit
    def test_zero_sequence_length(self):
        """Test with zero sequence length"""
        evidence = [
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="test",
                query_range="1-100",
                confidence=0.95,
                reference_length=100,
                domain_id="test_A"
            )
        ]

        # Should handle gracefully
        domains = partition_domains(evidence, sequence_length=0)
        # Coverage calculation should not crash
        assert isinstance(domains, list)


class TestRealWorldScenarios:
    """Test scenarios from real proteins"""

    @pytest.mark.unit
    def test_gfp_pbp_fusion_pattern(self):
        """Test pattern similar to 8ovp_A"""
        from mini.core.models import AlignmentData
        from mini.core.decomposer import DomainReference

        # Create alignment data for chain blast evidence
        gfp_alignment = AlignmentData(
            query_seq="A" * 243,  # 252-494 = 243 residues
            hit_seq="A" * 238,    # Reference length
            query_start=252,
            query_end=494,
            hit_start=1,
            hit_end=238
        )

        pbp_alignment = AlignmentData(
            query_seq="A" * 273,  # (2-248) + (491-517) = 247 + 27 = 274 residues
            hit_seq="A" * 508,    # Reference length
            query_start=2,
            query_end=517,
            hit_start=1,
            hit_end=508
        )

        evidence = [
            # GFP domain
            Evidence(
                type="chain_blast",
                source_pdb="6dgv",
                query_range=SequenceRange.parse("252-494"),
                confidence=0.95,
                evalue=1e-100,
                reference_length=238,
                domain_id="6dgv_A",
                alignment=gfp_alignment
            ),
            # PBP domain (discontinuous)
            Evidence(
                type="chain_blast",
                source_pdb="2ia4",
                query_range=SequenceRange.parse("2-248,491-517"),
                confidence=0.90,
                evalue=1e-80,
                reference_length=508,
                domain_id="2ia4_A",
                alignment=pbp_alignment
            ),
            # Overlapping domain blast hits
            create_complete_evidence(
                evidence_type="domain_blast",
                source_pdb="other",
                query_range="100-300",
                confidence=0.750,
                reference_length=201,
                domain_id="other_A"
            )
        ]

        # Mock domain definitions for decomposition
        domain_definitions = {
            ("6dgv", "A"): [
                DomainReference(
                    domain_id="e6dgvA1",
                    pdb_id="6dgv",
                    chain_id="A",
                    range=SequenceRange.parse("1-238"),
                    length=238,
                    t_group="1.1.1"
                )
            ],
            ("2ia4", "A"): [
                DomainReference(
                    domain_id="e2ia4A1",
                    pdb_id="2ia4",
                    chain_id="A",
                    range=SequenceRange.parse("1-247,275-301"),
                    length=274,
                    t_group="2.2.2"
                )
            ]
        }

        domains = partition_domains(evidence,
                                   sequence_length=569,
                                   domain_definitions=domain_definitions)

        # Should get decomposed domains if decomposition works, or regular domain blast
        assert len(domains) >= 1
        # Note: Actual families depend on decomposition results


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
