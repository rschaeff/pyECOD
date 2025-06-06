#!/usr/bin/env python3
"""
Model tests for mini_pyecod

Tests the data model classes.
"""

import pytest
from pathlib import Path

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.core.models import Evidence, Domain, AlignmentData
from mini.core.sequence_range import SequenceRange


class TestEvidenceModel:
    """Test the Evidence data model"""
    
    @pytest.mark.unit
    def test_evidence_creation_minimal(self):
        """Test creating evidence with minimal fields"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="1abc",
            query_range=SequenceRange.parse("10-100"),
        )
        
        assert evidence.type == "domain_blast"
        assert evidence.source_pdb == "1abc"
        assert str(evidence.query_range) == "10-100"
        assert evidence.confidence == 0.0  # Default
        assert evidence.evalue is None
        assert evidence.domain_id is None  # No auto-generation for minimal evidence

    @pytest.mark.unit
    def test_evidence_creation_full(self):
        """Test creating evidence with all fields"""
        evidence = Evidence(
            type="hhsearch",
            source_pdb="2def",
            query_range=SequenceRange.parse("1-50,100-150"),
            confidence=0.95,
            evalue=1e-50,
            domain_id="e2defA1",
            t_group="1234.5.6",
            h_group="1234.5",
            reference_length=101,
            alignment_coverage=0.85,
            alignment=AlignmentData(
                query_seq="ACDEF",
                hit_seq="ACDEF",
                query_start=1,
                query_end=50,
                hit_start=10,
                hit_end=60
            )
        )

        assert evidence.type == "hhsearch"
        assert evidence.confidence == 0.95
        assert evidence.evalue == 1e-50
        assert evidence.t_group == "1234.5.6"
        assert evidence.reference_length == 101
        assert evidence.alignment_coverage == 0.85
        assert evidence.alignment is not None
        assert evidence.alignment.query_seq == "ACDEF"

    @pytest.mark.unit
    def test_evidence_with_discontinuous_range(self):
        """Test evidence with discontinuous range - FIXED calculation"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="3ghi",
            query_range=SequenceRange.parse("10-50,60-100,150-200"),
            domain_id="e3ghiA1"  # Proper domain_id example
        )

        assert evidence.query_range.is_discontinuous
        assert len(evidence.query_range.segments) == 3
        # FIXED: (50-10+1) + (100-60+1) + (200-150+1) = 41 + 41 + 51 = 133
        assert evidence.query_range.total_length == 133


class TestDomainModel:
    """Test the Domain data model"""

    @pytest.mark.unit
    def test_domain_creation(self):
        """Test creating a domain"""
        evidence_items = [
            Evidence(
                type="domain_blast",
                source_pdb="test",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.9,
                domain_id="eTestA1"  # Proper domain_id
            )
        ]

        domain = Domain(
            id="d1",
            range=SequenceRange.parse("1-100"),
            family="test_family",
            evidence_count=1,
            source="domain_blast",
            evidence_items=evidence_items,
            t_group=None,
            h_group=None,
            x_group=None
        )

        assert domain.id == "d1"
        assert str(domain.range) == "1-100"
        assert domain.family == "test_family"
        assert domain.evidence_count == 1
        assert domain.source == "domain_blast"
        assert len(domain.evidence_items) == 1

    @pytest.mark.unit
    def test_domain_with_empty_evidence(self):
        """Test domain with no evidence items"""
        domain = Domain(
            id="d2",
            range=SequenceRange.parse("150-250"),
            family="another_family",
            evidence_count=0,
            source="unknown",
            evidence_items=[]
        )

        assert len(domain.evidence_items) == 0
        assert domain.evidence_count == 0

    @pytest.mark.unit
    def test_domain_discontinuous(self):
        """Test discontinuous domain"""
        domain = Domain(
            id="d3",
            range=SequenceRange.parse("1-50,100-150"),
            family="disc_family",
            evidence_count=2,
            source="chain_blast_decomposed",
            evidence_items=[]
        )

        assert domain.range.is_discontinuous
        assert domain.range.total_length == 101


class TestAlignmentDataModel:
    """Test the AlignmentData model"""

    @pytest.mark.unit
    def test_alignment_data_creation(self):
        """Test creating alignment data"""
        alignment = AlignmentData(
            query_seq="ACDEFGHIKLMNPQRSTVWY",
            hit_seq="ACDEFGHIKLMNPQRSTVWY",
            query_start=10,
            query_end=29,
            hit_start=1,
            hit_end=20
        )

        assert alignment.query_seq == "ACDEFGHIKLMNPQRSTVWY"
        assert alignment.hit_seq == "ACDEFGHIKLMNPQRSTVWY"
        assert alignment.query_start == 10
        assert alignment.query_end == 29
        assert alignment.hit_start == 1
        assert alignment.hit_end == 20

    @pytest.mark.unit
    def test_alignment_with_gaps(self):
        """Test alignment with gaps"""
        alignment = AlignmentData(
            query_seq="ACD-FGH",
            hit_seq="ACDEFGH",
            query_start=1,
            query_end=6,
            hit_start=1,
            hit_end=7
        )

        assert "-" in alignment.query_seq
        assert "-" not in alignment.hit_seq
        assert len(alignment.query_seq) == len(alignment.hit_seq)


class TestModelRelationships:
    """Test relationships between models"""

    @pytest.mark.unit
    def test_evidence_to_domain_relationship(self):
        """Test that evidence properly relates to domains"""
        # Create evidence
        evidence1 = Evidence(
            type="domain_blast",
            source_pdb="1abc",
            query_range=SequenceRange.parse("10-60"),
            confidence=0.8,
            t_group="1111.1.1",
            domain_id="e1abcA1"
        )

        evidence2 = Evidence(
            type="hhsearch",
            source_pdb="1abc",
            query_range=SequenceRange.parse("15-55"),
            confidence=0.9,
            t_group="1111.1.1",
            domain_id="e1abcA1"
        )

        # Create domain from evidence
        domain = Domain(
            id="d1",
            range=SequenceRange.parse("10-60"),  # Uses evidence1's range
            family=evidence1.t_group or evidence1.source_pdb,
            evidence_count=2,
            source=evidence1.type,
            evidence_items=[evidence1, evidence2]
        )

        # Verify relationships
        assert domain.family == "1111.1.1"
        assert domain.evidence_count == len(domain.evidence_items)
        assert domain.source == "domain_blast"
        assert all(e.t_group == "1111.1.1" for e in domain.evidence_items)

    @pytest.mark.unit
    def test_alignment_attached_to_evidence(self):
        """Test alignment data attached to evidence"""
        alignment = AlignmentData(
            query_seq="ABCDEFG",
            hit_seq="ABCDEFG",
            query_start=1,
            query_end=7,
            hit_start=10,
            hit_end=16
        )

        evidence = Evidence(
            type="chain_blast",
            source_pdb="2def",
            query_range=SequenceRange.parse("1-7"),
            alignment=alignment,
            domain_id="2def_A"  # Chain-level identifier for chain BLAST
        )

        assert evidence.alignment is not None
        assert evidence.alignment.query_start == 1
        assert evidence.alignment.hit_start == 10


class TestModelDefaults:
    """Test default values and optional fields"""

    @pytest.mark.unit
    def test_evidence_optional_fields(self):
        """Test that optional fields have proper defaults"""
        evidence = Evidence(
            type="test",
            source_pdb="test",
            query_range=SequenceRange.parse("1-10")
        )

        # Check defaults
        assert evidence.confidence == 0.0
        assert evidence.evalue is None
        assert evidence.domain_id is None
        assert evidence.t_group is None
        assert evidence.h_group is None
        assert evidence.reference_length is None
        assert evidence.alignment_coverage is None
        assert evidence.alignment is None

    @pytest.mark.unit
    def test_domain_evidence_list_initialization(self):
        """Test that evidence_items list is properly initialized"""
        # Domain with explicit empty list
        domain1 = Domain(
            id="d1",
            range=SequenceRange.parse("1-100"),
            family="test",
            evidence_count=0,
            source="test",
            evidence_items=[]
        )

        assert domain1.evidence_items == []
        assert isinstance(domain1.evidence_items, list)

        # Domain with evidence
        evidence = Evidence(
            type="test",
            source_pdb="test",
            query_range=SequenceRange.parse("1-100"),
            domain_id="eTestA1"
        )

        domain2 = Domain(
            id="d2",
            range=SequenceRange.parse("1-100"),
            family="test",
            evidence_count=1,
            source="test",
            evidence_items=[evidence]
        )

        assert len(domain2.evidence_items) == 1
        assert domain2.evidence_items[0] is evidence


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
