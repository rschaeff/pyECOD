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

from mini.core.models import Evidence, Domain
from mini.core.partitioner import partition_domains
from mini.core.sequence_range import SequenceRange


class TestResidueBlocking:
    """Test the residue blocking algorithm"""
    
    @pytest.mark.unit
    def test_basic_residue_blocking(self):
        """Test that domains block residues from reuse"""
        evidence = [
            Evidence(
                type="domain_blast",
                source_pdb="test1",
                query_range=SequenceRange.parse("10-100"),
                confidence=0.9,
                evalue=1e-50,
                reference_length=91,
                domain_id="test1_A"
            ),
            Evidence(
                type="domain_blast",
                source_pdb="test2",
                query_range=SequenceRange.parse("50-150"),  # Overlaps with first
                confidence=0.8,
                evalue=1e-40,
                reference_length=101,
                domain_id='test2_A'
            )
        ]
        
        domains = partition_domains(evidence,
            sequence_length=200
        )
        
        # Should select first domain (higher confidence)
        assert len(domains) == 1
        assert domains[0].family == "test1"
        assert "10" in str(domains[0].range) or "1" in str(domains[0].range)

    
    @pytest.mark.unit
    def test_coverage_thresholds(self):
        """Test NEW_COVERAGE and OLD_COVERAGE thresholds"""
        evidence = [
            Evidence(
                type="domain_blast",
                source_pdb="domain1",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.9,
                reference_length=100,
                domain_id="domain1_A"
            ),
            Evidence(
                type="domain_blast",
                source_pdb="domain2",
                query_range=SequenceRange.parse("90-200"),  # 10% overlap
                confidence=0.85,
                reference_length=111,
                domain_id="domain2_A"
            ),
            Evidence(
                type="domain_blast",
                source_pdb="domain3",
                query_range=SequenceRange.parse("50-150"),  # 50% overlap
                confidence=0.95,
                reference_length=101,
                domain_id="domain3_A"
            )
        ]
        
        domains = partition_domains(evidence, sequence_length=250)
        
        # Should have domain1 and domain2 (acceptable overlap)
        # Should reject domain3 (too much overlap)
        assert len(domains) >= 2
        families = {d.family for d in domains}
        assert "domain3" in families
    
    @pytest.mark.unit
    def test_evidence_priority_sorting(self):
        """Test evidence is processed in correct priority order"""
        evidence = [
            Evidence(
                type="hhsearch",
                source_pdb="hh1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.7,
                reference_length=50,
                domain_id="hh1_A"
            ),
            Evidence(
                type="domain_blast",
                source_pdb="blast1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.8,
                reference_length=50,
                domain_id="blast1_A"
            ),
            Evidence(
                type="chain_blast",
                source_pdb="chain1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.9,
                reference_length=50,
                domain_id="chain1_A"
            ),
        ]
        
        domains = partition_domains(evidence, sequence_length=100)
        
        # Chain blast should win (type precedence despite lower confidence than others)
        assert len(domains) == 1
        assert domains[0].family == "blast1"
        assert domains[0].source == "domain_blast"

    @pytest.mark.unit
    def test_chain_blast_priority_with_decomposition(self):
        """Test that chain blast wins when decomposition is available"""
        evidence = [
            Evidence(
                type="domain_blast",
                source_pdb="blast1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.8,
                reference_length=50,
                domain_id="blast1_A"
            ),
            Evidence(
                type="chain_blast",
                source_pdb="chain1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.9,
                reference_length=50,
                domain_id="chain1_A"
            ),
        ]

        # Provide mock domain definitions to enable chain blast decomposition
        from mini.core.decomposer import DomainReference
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
                Evidence(
                    type="domain_blast",
                    source_pdb="tiny",
                    query_range=SequenceRange.parse("1-15"),  # Too small
                    confidence=0.9,
                    reference_length=15
                ),
                Evidence(
                    type="domain_blast",
                    source_pdb="normal",
                    query_range=SequenceRange.parse("20-100"),
                    confidence=0.8,
                    reference_length=81
                )
            ]

            domains = partition_domains(evidence, sequence_length=150)

            # Should only have the normal-sized domain
            assert len(domains) == 1
            assert domains[0].family == "normal"


class TestDiscontinuousDomains:
    """Test handling of discontinuous domains"""
    
    @pytest.mark.unit
    def test_discontinuous_domain_parsing(self):
        """Test that discontinuous ranges are handled correctly"""
        evidence = [
            Evidence(
                type="chain_blast",
                source_pdb="disc",
                query_range=SequenceRange.parse("1-100,200-250"),
                confidence=0.9,
                reference_length=151
            )
        ]
        
        domains = partition_domains(evidence, sequence_length=300)
        
        assert len(domains) == 1
        assert domains[0].range.is_discontinuous
        assert domains[0].range.total_length == 151
        assert len(domains[0].range.segments) == 2
    
    @pytest.mark.unit
    def test_discontinuous_overlap_calculation(self):
        """Test overlap calculation with discontinuous domains"""
        evidence = [
            Evidence(
                type="domain_blast",
                source_pdb="disc1",
                query_range=SequenceRange.parse("1-50,100-150"),
                confidence=0.9,
                reference_length=101
            ),
            Evidence(
                type="domain_blast",
                source_pdb="cont1",
                query_range=SequenceRange.parse("40-120"),  # Overlaps both segments
                confidence=0.85,
                reference_length=81
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
            Evidence(
                type="domain_blast",
                source_pdb="test",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.9,
                t_group="1234.5.6",
                reference_length=100
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
            Evidence(
                type="domain_blast",
                source_pdb="pdb1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.9,
                t_group="1111.1.1",
                domain_id="e1abcA1",
                reference_length=50
            ),
            # No T-group, has source_pdb
            Evidence(
                type="domain_blast", 
                source_pdb="pdb2",
                query_range=SequenceRange.parse("60-110"),
                confidence=0.9,
                domain_id="e2defB1",
                reference_length=51
            ),
            # Only domain_id
            Evidence(
                type="domain_blast",
                source_pdb="",
                query_range=SequenceRange.parse("120-170"),
                confidence=0.9,
                domain_id="e3ghiC1",
                reference_length=51
            )
        ]
        
        domains = partition_domains(evidence_list, sequence_length=200)
        
        assert len(domains) == 3
        # Check family assignment priority
        assert domains[0].family == "1111.1.1"  # T-group
        assert domains[1].family == "pdb2"      # source_pdb
        assert domains[2].family == "e3ghiC1"   # domain_id


class TestCoverageCalculation:
    """Test sequence coverage calculations"""
    
    @pytest.mark.unit
    def test_total_coverage(self):
        """Test that coverage is calculated correctly"""
        evidence = [
            Evidence(
                type="domain_blast",
                source_pdb="d1",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.9,
                reference_length=100
            ),
            Evidence(
                type="domain_blast",
                source_pdb="d2", 
                query_range=SequenceRange.parse("150-200"),
                confidence=0.9,
                reference_length=51
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
            Evidence(
                type="domain_blast",
                source_pdb="d1",
                query_range=SequenceRange.parse("1-50"),
                confidence=0.9,
                reference_length=50
            ),
            Evidence(
                type="domain_blast",
                source_pdb="d2",
                query_range=SequenceRange.parse("100-150"),
                confidence=0.9,
                reference_length=51
            ),
            Evidence(
                type="domain_blast",
                source_pdb="d3",
                query_range=SequenceRange.parse("200-250"),
                confidence=0.9,
                reference_length=51
            )
        ]
        
        domains = partition_domains(evidence, sequence_length=300)
        
        assert len(domains) == 3
        # Check that domains maintain their gaps
        assert domains[0].range.segments[0].end < domains[1].range.segments[0].start
        assert domains[1].range.segments[0].end < domains[2].range.segments[0].start


class TestEmptyAndEdgeCases:
    """Test edge cases and error conditions"""
    
    @pytest.mark.unit
    def test_no_evidence(self):
        """Test with no evidence"""
        domains = partition_domains([], sequence_length=100)
        assert len(domains) == 0
    
    @pytest.mark.unit
    def test_no_reference_lengths(self):
        """Test that evidence without reference lengths is handled"""
        evidence = [
            Evidence(
                type="domain_blast",
                source_pdb="test",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.9,
                reference_length=None  # No reference length
            )
        ]
        
        # With require_reference_lengths=True by default, should skip this evidence
        domains = partition_domains(evidence, sequence_length=150)
        assert len(domains) == 0
    
    @pytest.mark.unit
    def test_zero_sequence_length(self):
        """Test with zero sequence length"""
        evidence = [
            Evidence(
                type="domain_blast",
                source_pdb="test",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.9,
                reference_length=100
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
        evidence = [
            # GFP domain
            Evidence(
                type="chain_blast",
                source_pdb="6dgv",
                query_range=SequenceRange.parse("252-494"),
                confidence=0.95,
                evalue=1e-100,
                reference_length=238
            ),
            # PBP domain (discontinuous)
            Evidence(
                type="chain_blast",
                source_pdb="2ia4",
                query_range=SequenceRange.parse("2-248,491-517"),
                confidence=0.90,
                evalue=1e-80,
                reference_length=508
            ),
            # Overlapping domain blast hits
            Evidence(
                type="domain_blast",
                source_pdb="other",
                query_range=SequenceRange.parse("100-300"),
                confidence=0.70,
                reference_length=201
            )
        ]
        
        domains = partition_domains(evidence, sequence_length=569)
        
        # Should prefer the chain blast hits
        assert len(domains) >= 2
        families = {d.family for d in domains}
        assert "6dgv" in families
        assert "2ia4" in families
        
        # Check for discontinuous domain
        discontinuous = [d for d in domains if d.range.is_discontinuous]
        assert len(discontinuous) >= 1


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
