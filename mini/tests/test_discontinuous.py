#!/usr/bin/env python3
"""
Discontinuous domain decomposition tests for mini_pyecod

Tests the handling of discontinuous domains and chain BLAST decomposition.
"""

import pytest
from pathlib import Path
from typing import List

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.core.models import Evidence, AlignmentData
from mini.core.decomposer import (
    decompose_chain_blast_discontinuous,
    decompose_chain_blast_with_mapping,
    build_alignment_mapping,
    DomainReference
)
from mini.core.sequence_range import SequenceRange


class TestDiscontinuousDecomposition:
    """Test discontinuous chain BLAST decomposition"""
    
    @pytest.mark.unit
    def test_basic_discontinuous_decomposition(self):
        """Test that discontinuous hits are properly decomposed by segments"""
        # Create a discontinuous chain BLAST hit
        evidence = Evidence(
            type="chain_blast",
            source_pdb="2ia4",
            query_range=SequenceRange.parse("2-248,491-517"),
            confidence=0.95,
            evalue=1e-100,
            domain_id="2ia4_A",
            alignment_coverage=0.95
        )
        
        # Decompose without reference data
        decomposed = decompose_chain_blast_discontinuous(evidence)
        
        # Should split into 2 continuous segments
        assert len(decomposed) == 2
        
        # First segment
        assert str(decomposed[0].query_range) == "2-248"
        assert decomposed[0].type == "chain_blast_decomposed"
        assert decomposed[0].domain_id == "2ia4_A_seg1"
        assert decomposed[0].confidence < evidence.confidence  # Slightly reduced
        
        # Second segment
        assert str(decomposed[1].query_range) == "491-517"
        assert decomposed[1].type == "chain_blast_decomposed"
        assert decomposed[1].domain_id == "2ia4_A_seg2"
    
    @pytest.mark.unit
    def test_continuous_hit_not_decomposed(self):
        """Test that continuous hits are returned unchanged"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("252-494"),
            confidence=0.9,
            domain_id="6dgv_A"
        )
        
        decomposed = decompose_chain_blast_discontinuous(evidence)
        
        # Should return original evidence unchanged
        assert len(decomposed) == 1
        assert decomposed[0] is evidence
    
    @pytest.mark.unit
    def test_small_segment_filtering(self):
        """Test that tiny segments are filtered out"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="test",
            query_range=SequenceRange.parse("1-100,105-110,200-300"),  # Middle segment too small
            confidence=0.9,
            domain_id="test_A"
        )
        
        decomposed = decompose_chain_blast_discontinuous(evidence, min_domain=50)
        
        # Should only have 2 segments (tiny middle one filtered)
        assert len(decomposed) == 2
        assert str(decomposed[0].query_range) == "1-100"
        assert str(decomposed[1].query_range) == "200-300"
    
    @pytest.mark.unit
    def test_all_segments_too_small(self):
        """Test handling when all segments are too small"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="tiny",
            query_range=SequenceRange.parse("1-10,20-30,40-45"),
            confidence=0.9,
            domain_id="tiny_A"
        )
        
        decomposed = decompose_chain_blast_discontinuous(evidence, min_domain=50)
        
        # Should return original if no segments meet threshold
        assert len(decomposed) == 1
        assert decomposed[0] is evidence


class TestAlignmentMapping:
    """Test alignment-based position mapping"""
    
    @pytest.mark.unit
    def test_simple_alignment_mapping(self):
        """Test basic alignment mapping without gaps"""
        query_str = "ABCDEFG"
        hit_str = "ABCDEFG"
        query_start = 10
        hit_start = 20
        
        mapping = build_alignment_mapping(query_str, hit_str, query_start, hit_start)
        
        # Check mapping (0-based internally)
        assert mapping[9] == 19   # Position 10 -> 20
        assert mapping[10] == 20  # Position 11 -> 21
        assert mapping[15] == 25  # Position 16 -> 26
        assert len(mapping) == 7
    
    @pytest.mark.unit
    def test_alignment_with_query_gaps(self):
        """Test alignment mapping with gaps in query"""
        query_str = "A-CD-FG"
        hit_str = "ABCDEFG"
        query_start = 1
        hit_start = 1
        
        mapping = build_alignment_mapping(query_str, hit_str, query_start, hit_start)
        
        # Only non-gap positions are mapped
        assert len(mapping) == 5  # A, C, D, F, G
        assert mapping[0] == 0    # A -> A
        assert mapping[1] == 2    # C -> C
        assert mapping[2] == 3    # D -> D
        assert mapping[3] == 5    # F -> F
        assert mapping[4] == 6    # G -> G
    
    @pytest.mark.unit
    def test_alignment_with_hit_gaps(self):
        """Test alignment mapping with gaps in hit"""
        query_str = "ABCDEFG"
        hit_str = "AB--EFG"
        query_start = 1
        hit_start = 1
        
        mapping = build_alignment_mapping(query_str, hit_str, query_start, hit_start)
        
        # Only positions where both have residues are mapped
        assert len(mapping) == 5  # A, B, E, F, G
        assert mapping[0] == 0    # A -> A
        assert mapping[1] == 1    # B -> B
        assert mapping[4] == 2    # E -> E (position 3 in hit)
        assert mapping[5] == 3    # F -> F
        assert mapping[6] == 4    # G -> G
    
    @pytest.mark.unit
    def test_mismatched_alignment_lengths(self):
        """Test error handling for mismatched alignment lengths"""
        query_str = "ABC"
        hit_str = "ABCDEF"  # Different length
        
        with pytest.raises(ValueError, match="different lengths"):
            build_alignment_mapping(query_str, hit_str, 1, 1)


class TestMappingBasedDecomposition:
    """Test decomposition using alignment mapping and reference domains"""
    
    @pytest.mark.unit
    def test_single_domain_mapping(self):
        """Test mapping a single reference domain"""
        # Mock evidence with alignment
        evidence = Evidence(
            type="chain_blast",
            source_pdb="ref",
            query_range=SequenceRange.parse("10-60"),
            confidence=0.9,
            domain_id="ref_A"
        )
        
        # Mock alignment data (simple 1:1)
        query_str = "A" * 51
        hit_str = "A" * 51
        
        # Single reference domain
        ref_domains = [
            DomainReference(
                domain_id="eRefA1",
                pdb_id="ref",
                chain_id="A",
                range=SequenceRange.parse("5-30"),
                length=26,
                t_group="1234.1.1"
            )
        ]
        
        decomposed = decompose_chain_blast_with_mapping(
            evidence, query_str, hit_str, 10, 5, ref_domains
        )
        
        # Should map the reference domain to query
        assert len(decomposed) == 1
        assert decomposed[0].type == "chain_blast_decomposed"
        assert decomposed[0].domain_id == "eRefA1"
        assert decomposed[0].t_group == "1234.1.1"
        # Query positions 10-35 map to hit positions 5-30
        assert str(decomposed[0].query_range) == "10-35"
    
    @pytest.mark.unit
    def test_multi_domain_mapping(self):
        """Test mapping multiple reference domains"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="multi",
            query_range=SequenceRange.parse("1-200"),
            confidence=0.9,
            domain_id="multi_A"
        )
        
        # Simple alignment
        query_str = "A" * 200
        hit_str = "A" * 200
        
        # Two reference domains
        ref_domains = [
            DomainReference(
                domain_id="eMultiA1",
                pdb_id="multi",
                chain_id="A",
                range=SequenceRange.parse("10-90"),
                length=81,
                t_group="1111.1.1"
            ),
            DomainReference(
                domain_id="eMultiA2",
                pdb_id="multi",
                chain_id="A",
                range=SequenceRange.parse("100-180"),
                length=81,
                t_group="2222.2.2"
            )
        ]
        
        decomposed = decompose_chain_blast_with_mapping(
            evidence, query_str, hit_str, 1, 1, ref_domains
        )
        
        # Should get both domains
        assert len(decomposed) == 2
        
        # First domain
        assert decomposed[0].domain_id == "eMultiA1"
        assert decomposed[0].t_group == "1111.1.1"
        assert str(decomposed[0].query_range) == "10-90"
        
        # Second domain
        assert decomposed[1].domain_id == "eMultiA2"
        assert decomposed[1].t_group == "2222.2.2"
        assert str(decomposed[1].query_range) == "100-180"
    
    @pytest.mark.unit
    def test_partial_coverage_filtering(self):
        """Test that poorly covered domains are filtered"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="partial",
            query_range=SequenceRange.parse("1-100"),
            confidence=0.9
        )
        
        # Alignment covers only part of the reference
        query_str = "A" * 100
        hit_str = "A" * 100
        
        ref_domains = [
            DomainReference(
                domain_id="ePartialA1",
                pdb_id="partial",
                chain_id="A",
                range=SequenceRange.parse("1-200"),  # Much larger than alignment
                length=200
            )
        ]
        
        decomposed = decompose_chain_blast_with_mapping(
            evidence, query_str, hit_str, 1, 1, ref_domains
        )
        
        # Should filter out due to poor coverage
        # Or return with reduced confidence
        if decomposed:
            assert decomposed[0].alignment_coverage < 0.6  # Low coverage
            assert decomposed[0].confidence < evidence.confidence
    
    @pytest.mark.unit
    def test_no_reference_domains(self):
        """Test handling when no reference domains available"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="noref",
            query_range=SequenceRange.parse("1-100"),
            confidence=0.9
        )
        
        decomposed = decompose_chain_blast_with_mapping(
            evidence, "A"*100, "A"*100, 1, 1, []  # Empty reference list
        )
        
        # Should return original evidence
        assert len(decomposed) == 1
        assert decomposed[0] is evidence


class TestRealWorldDecomposition:
    """Test with patterns from real proteins"""
    
    @pytest.mark.unit
    def test_pbp_like_discontinuous_pattern(self):
        """Test pattern similar to PBP domains with insertions"""
        # Discontinuous hit like 2ia4 in 8ovp_A
        evidence = Evidence(
            type="chain_blast",
            source_pdb="2ia4",
            query_range=SequenceRange.parse("2-248,491-517"),
            confidence=0.95,
            evalue=1e-80,
            domain_id="2ia4_B"
        )
        
        # Mock alignment that shows the discontinuity
        alignment = AlignmentData(
            query_seq="M" * 247 + "X" * 27,  # Simplified
            hit_seq="M" * 247 + "X" * 27,
            query_start=2,
            query_end=517,
            hit_start=1,
            hit_end=274
        )
        
        evidence.alignment = alignment
        
        # Reference has 2 domains
        ref_domains = [
            DomainReference(
                domain_id="e2ia4B1",
                pdb_id="2ia4",
                chain_id="B",
                range=SequenceRange.parse("1-120"),
                length=120,
                t_group="7523.1.1.1"
            ),
            DomainReference(
                domain_id="e2ia4B2",
                pdb_id="2ia4",
                chain_id="B",
                range=SequenceRange.parse("121-274"),
                length=154,
                t_group="7523.1.1.2"
            )
        ]
        
        # First try discontinuous decomposition
        simple_decomposed = decompose_chain_blast_discontinuous(evidence)
        assert len(simple_decomposed) == 2
        
        # Then try mapping-based (would need proper alignment)
        # This is simplified - real case would map through alignment
        if hasattr(evidence, 'alignment'):
            mapped_decomposed = decompose_chain_blast_with_mapping(
                evidence,
                alignment.query_seq,
                alignment.hit_seq,
                alignment.query_start,
                alignment.hit_start,
                ref_domains
            )
            # Would produce mapped domains if alignment was real


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
