#!/usr/bin/env python3
"""
Comprehensive test coverage for ecod.models.pipeline.domain module

Tests the DomainModel class functionality including:
- Basic domain creation and properties
- Evidence handling and confidence calculation  
- Classification logic
- Overlap calculations
- Serialization (dict and XML)
- Domain operations (merge, split)
- Edge cases and error handling
"""

import pytest
import xml.etree.ElementTree as ET
from unittest.mock import Mock, patch

from ecod.models.pipeline.domain import DomainModel, Domain
from ecod.models.pipeline.evidence import Evidence


class TestDomainModelBasics:
    """Test basic DomainModel creation and properties"""
    
    def test_basic_domain_creation(self):
        """Test creating a basic domain"""
        domain = DomainModel(
            id="test_domain",
            start=10,
            end=100,
            range="10-100"
        )
        
        assert domain.id == "test_domain"
        assert domain.start == 10
        assert domain.end == 100
        assert domain.range == "10-100"
        assert domain.size == 91  # 100 - 10 + 1
        assert domain.length == 91  # Alias for size
        assert domain.confidence == 0.0  # Default when no evidence
    
    def test_domain_with_full_classification(self):
        """Test domain with complete ECOD classification"""
        domain = DomainModel(
            id="d1abcA1",
            start=1,
            end=150,
            range="1-150",
            t_group="2002.1.1.1",
            h_group="2002.1.1",
            x_group="2002.1",
            a_group="a.39",
            source="hhsearch",
            confidence=0.95,
            source_id="e4xfcA1"
        )
        
        assert domain.is_fully_classified()
        assert domain.is_classified()
        assert domain.get_classification_level() == "A-group"
        assert domain.confidence == 0.95
        assert domain.source == "hhsearch"
    
    def test_domain_with_partial_classification(self):
        """Test domain with partial classification"""
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100",
            t_group="2002.1.1.1",
            h_group="2002.1.1"
            # Missing x_group and a_group
        )
        
        assert not domain.is_fully_classified()
        assert domain.is_classified()
        assert domain.get_classification_level() == "H-group"
    
    def test_unclassified_domain(self):
        """Test completely unclassified domain"""
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100"
        )
        
        assert not domain.is_fully_classified()
        assert not domain.is_classified()
        assert domain.get_classification_level() == "Unclassified"
    
    def test_domain_properties(self):
        """Test domain size and position properties"""
        domain = DomainModel(
            id="test",
            start=20,
            end=30,
            range="20-30"
        )
        
        assert domain.size == 11
        assert domain.length == 11
        positions = domain.get_positions()
        assert positions == {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30}
        assert len(positions) == 11
    
    def test_range_segments_parsing(self):
        """Test parsing of multi-segment ranges"""
        # Simple range
        domain1 = DomainModel(id="test1", start=10, end=20, range="10-20")
        segments1 = domain1.get_range_segments()
        assert segments1 == [(10, 20)]
        
        # Multi-segment range
        domain2 = DomainModel(id="test2", start=10, end=40, range="10-20,30-40")
        segments2 = domain2.get_range_segments()
        assert segments2 == [(10, 20), (30, 40)]
        
        # Empty range falls back to start-end
        domain3 = DomainModel(id="test3", start=5, end=15, range="")
        segments3 = domain3.get_range_segments()
        assert segments3 == [(5, 15)]
        
        # Invalid range falls back to start-end
        domain4 = DomainModel(id="test4", start=5, end=15, range="invalid-range")
        segments4 = domain4.get_range_segments()
        assert segments4 == [(5, 15)]


class TestDomainModelEvidence:
    """Test evidence handling and confidence calculation"""
    
    def test_domain_with_dict_evidence(self):
        """Test domain creation with dictionary evidence"""
        evidence_dict = {
            "type": "hhsearch",
            "source_id": "d1abcA1",
            "probability": 95.0,
            "evalue": 1e-10,
            "confidence": 0.95
        }
        
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100",
            evidence=[evidence_dict]
        )
        
        # Should have converted dict to Evidence object
        assert len(domain.evidence) == 1
        assert isinstance(domain.evidence[0], Evidence)
        assert domain.evidence[0].type == "hhsearch"
        assert domain.evidence[0].probability == 95.0
    
    def test_domain_with_evidence_objects(self):
        """Test domain creation with Evidence objects"""
        evidence = Evidence(
            type="domain_blast",
            source_id="d2xyzB1",
            evalue=1e-8,
            confidence=0.85
        )
        
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100",
            evidence=[evidence]
        )
        
        assert len(domain.evidence) == 1
        assert isinstance(domain.evidence[0], Evidence)
        assert domain.evidence[0].type == "domain_blast"
        assert domain.evidence[0].evalue == 1e-8
    
    def test_confidence_calculation_from_evidence(self):
        """Test confidence calculation from multiple evidence sources"""
        evidence1 = Evidence(type="domain_blast", evalue=1e-10, confidence=0.9)
        evidence2 = Evidence(type="hhsearch", probability=85.0, confidence=0.8)
        evidence3 = Evidence(type="chain_blast", evalue=1e-5, confidence=0.6)
        
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100",
            evidence=[evidence1, evidence2, evidence3],
            confidence=0.0  # Will be recalculated from evidence
        )
        
        # Should calculate weighted average based on evidence types
        # domain_blast (weight 3.0): 0.9 * 3.0 = 2.7
        # hhsearch (weight 2.5): 0.8 * 2.5 = 2.0  
        # chain_blast (weight 2.0): 0.6 * 2.0 = 1.2
        # Total: 5.9, Total weight: 7.5, Average: 5.9/7.5 ≈ 0.787
        
        assert domain.confidence > 0.75
        assert domain.confidence < 0.85
    
    def test_add_evidence(self):
        """Test adding evidence to existing domain"""
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100"
        )
        
        initial_confidence = domain.confidence
        
        # Add evidence
        evidence = Evidence(type="hhsearch", probability=90.0, confidence=0.9)
        domain.add_evidence(evidence)
        
        assert len(domain.evidence) == 1
        assert domain.confidence > initial_confidence
        assert domain.confidence > 0.8
    
    def test_classification_update_from_evidence(self):
        """Test that classification is updated from evidence"""
        evidence_dict = {
            "type": "hhsearch",
            "source_id": "d1abcA1",
            "t_group": "2002.1.1.1",
            "h_group": "2002.1.1",
            "x_group": "2002.1",
            "a_group": "a.39"
        }
        
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100"
        )
        
        # Initially unclassified
        assert not domain.is_classified()
        
        # Add evidence with classification
        domain.add_evidence(evidence_dict)
        
        # Should now be classified
        assert domain.is_fully_classified()
        assert domain.t_group == "2002.1.1.1"
        assert domain.h_group == "2002.1.1"
        assert domain.x_group == "2002.1"
        assert domain.a_group == "a.39"
    
    def test_protected_flag_for_high_confidence(self):
        """Test that very high confidence domains are marked as protected"""
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100",
            confidence=0.99  # Very high confidence
        )
        
        assert domain.protected == True
        
        # Lower confidence should not be protected
        domain2 = DomainModel(
            id="test2",
            start=1,
            end=100,
            range="1-100",
            confidence=0.8
        )
        
        assert domain2.protected == False


class TestDomainModelOverlaps:
    """Test domain overlap calculations"""
    
    def test_overlapping_domains(self):
        """Test detection of overlapping domains"""
        domain1 = DomainModel(id="d1", start=10, end=30, range="10-30")
        domain2 = DomainModel(id="d2", start=25, end=45, range="25-45")
        
        assert domain1.overlaps(domain2)
        assert domain2.overlaps(domain1)
        
        # Overlap size should be 6 (25, 26, 27, 28, 29, 30)
        assert domain1.overlap_size(domain2) == 6
        assert domain2.overlap_size(domain1) == 6
    
    def test_non_overlapping_domains(self):
        """Test non-overlapping domains"""
        domain1 = DomainModel(id="d1", start=10, end=20, range="10-20")
        domain2 = DomainModel(id="d2", start=30, end=40, range="30-40")
        
        assert not domain1.overlaps(domain2)
        assert not domain2.overlaps(domain1)
        assert domain1.overlap_size(domain2) == 0
        assert domain2.overlap_size(domain1) == 0
    
    def test_adjacent_domains(self):
        """Test adjacent (touching) domains"""
        domain1 = DomainModel(id="d1", start=10, end=20, range="10-20")
        domain2 = DomainModel(id="d2", start=21, end=30, range="21-30")
        
        # Adjacent domains should not overlap
        assert not domain1.overlaps(domain2)
        assert domain1.overlap_size(domain2) == 0
    
    def test_contained_domains(self):
        """Test domain completely contained in another"""
        domain1 = DomainModel(id="d1", start=10, end=50, range="10-50")  # Large
        domain2 = DomainModel(id="d2", start=20, end=30, range="20-30")  # Small, inside
        
        assert domain1.overlaps(domain2)
        assert domain2.overlaps(domain1)
        
        # Overlap size should be size of smaller domain
        assert domain1.overlap_size(domain2) == 11  # Size of domain2
        assert domain2.overlap_size(domain1) == 11
    
    def test_overlap_percentage(self):
        """Test overlap percentage calculation"""
        domain1 = DomainModel(id="d1", start=10, end=20, range="10-20")  # Size 11
        domain2 = DomainModel(id="d2", start=15, end=25, range="15-25")  # Size 11
        
        # Overlap is 15-20 = 6 positions
        # Percentage = 6/11 ≈ 54.5%
        overlap_pct = domain1.overlap_percentage(domain2)
        assert 54.0 < overlap_pct < 55.0
    
    def test_identical_domains(self):
        """Test identical domains"""
        domain1 = DomainModel(id="d1", start=10, end=20, range="10-20")
        domain2 = DomainModel(id="d2", start=10, end=20, range="10-20")
        
        assert domain1.overlaps(domain2)
        assert domain1.overlap_size(domain2) == 11
        assert domain1.overlap_percentage(domain2) == 100.0


class TestDomainModelSerialization:
    """Test domain serialization to/from dict and XML"""
    
    def test_to_dict_basic(self):
        """Test basic domain to dictionary conversion"""
        domain = DomainModel(
            id="d1abcA1",
            start=10,
            end=100,
            range="10-100",
            source="hhsearch",
            confidence=0.9,
            t_group="2002.1.1.1",
            h_group="2002.1.1"
        )
        
        domain_dict = domain.to_dict()
        
        assert domain_dict["id"] == "d1abcA1"
        assert domain_dict["start"] == 10
        assert domain_dict["end"] == 100
        assert domain_dict["range"] == "10-100"
        assert domain_dict["size"] == 91
        assert domain_dict["source"] == "hhsearch"
        assert domain_dict["confidence"] == 0.9
        assert domain_dict["t_group"] == "2002.1.1.1"
        assert domain_dict["h_group"] == "2002.1.1"
    
    def test_to_dict_with_evidence(self):
        """Test domain to dict with evidence included"""
        evidence = Evidence(type="hhsearch", probability=90.0, confidence=0.9)
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100",
            evidence=[evidence]
        )
        
        domain_dict = domain.to_dict()
        
        assert "evidence" in domain_dict
        assert len(domain_dict["evidence"]) == 1
        assert domain_dict["evidence"][0]["type"] == "hhsearch"
        assert domain_dict["evidence"][0]["probability"] == 90.0
    
    def test_from_dict_basic(self):
        """Test creating domain from dictionary"""
        domain_dict = {
            "id": "test_domain",
            "start": 10,
            "end": 100,
            "range": "10-100",
            "source": "blast",
            "confidence": 0.8,
            "t_group": "2002.1.1.1"
        }
        
        domain = DomainModel.from_dict(domain_dict)
        
        assert domain.id == "test_domain"
        assert domain.start == 10
        assert domain.end == 100
        assert domain.range == "10-100"
        assert domain.source == "blast"
        assert domain.confidence == 0.8
        assert domain.t_group == "2002.1.1.1"
    
    def test_from_dict_with_evidence(self):
        """Test creating domain from dict with evidence"""
        domain_dict = {
            "id": "test",
            "start": 1,
            "end": 100,
            "range": "1-100",
            "evidence": [
                {
                    "type": "hhsearch",
                    "probability": 85.0,
                    "confidence": 0.85
                }
            ]
        }
        
        domain = DomainModel.from_dict(domain_dict)
        
        assert len(domain.evidence) == 1
        assert isinstance(domain.evidence[0], Evidence)
        assert domain.evidence[0].type == "hhsearch"
        assert domain.evidence[0].probability == 85.0
    
    def test_dict_round_trip(self):
        """Test domain dict serialization round-trip"""
        evidence = Evidence(type="blast", evalue=1e-8, confidence=0.8)
        original = DomainModel(
            id="test",
            start=10,
            end=100,
            range="10-100",
            source="blast",
            confidence=0.8,
            t_group="2002.1.1.1",
            evidence=[evidence]
        )
        
        domain_dict = original.to_dict()
        reconstructed = DomainModel.from_dict(domain_dict)
        
        assert reconstructed.id == original.id
        assert reconstructed.start == original.start
        assert reconstructed.end == original.end
        assert reconstructed.range == original.range
        assert reconstructed.source == original.source
        assert reconstructed.confidence == original.confidence
        assert reconstructed.t_group == original.t_group
        assert len(reconstructed.evidence) == len(original.evidence)
    
    def test_to_xml_basic(self):
        """Test domain to XML conversion"""
        domain = DomainModel(
            id="d1abcA1",
            start=10,
            end=100,
            range="10-100",
            source="hhsearch",
            confidence=0.9,
            t_group="2002.1.1.1"
        )
        
        xml_element = domain.to_xml()
        
        assert xml_element.tag == "domain"
        assert xml_element.get("id") == "d1abcA1"
        assert xml_element.get("start") == "10"
        assert xml_element.get("end") == "100"
        assert xml_element.get("range") == "10-100"
        assert xml_element.get("source") == "hhsearch"
        assert xml_element.get("confidence") == "0.9000"
        assert xml_element.get("t_group") == "2002.1.1.1"
    
    def test_to_xml_with_evidence(self):
        """Test domain to XML with evidence"""
        evidence = Evidence(type="hhsearch", probability=90.0, confidence=0.9)
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100",
            evidence=[evidence]
        )
        
        xml_element = domain.to_xml()
        
        # Should have evidence_list element
        evidence_list = xml_element.find("evidence_list")
        assert evidence_list is not None
        assert evidence_list.get("count") == "1"
        
        # Should have evidence elements
        evidence_elements = evidence_list.findall("evidence")
        assert len(evidence_elements) == 1
        assert evidence_elements[0].get("type") == "hhsearch"
    
    def test_from_xml_basic(self):
        """Test creating domain from XML"""
        xml_str = '''
        <domain id="test" start="10" end="100" range="10-100" 
                source="blast" confidence="0.8000" t_group="2002.1.1.1">
        </domain>
        '''
        xml_element = ET.fromstring(xml_str)
        
        domain = DomainModel.from_xml(xml_element)
        
        assert domain.id == "test"
        assert domain.start == 10
        assert domain.end == 100
        assert domain.range == "10-100"
        assert domain.source == "blast"
        assert domain.confidence == 0.8
        assert domain.t_group == "2002.1.1.1"
    
    def test_xml_round_trip(self):
        """Test domain XML serialization round-trip"""
        evidence = Evidence(type="hhsearch", probability=95.0, confidence=0.95)
        original = DomainModel(
            id="test",
            start=10,
            end=100,
            range="10-100",
            source="hhsearch",
            confidence=0.9,
            t_group="2002.1.1.1",
            evidence=[evidence]
        )
        
        xml_element = original.to_xml()
        reconstructed = DomainModel.from_xml(xml_element)
        
        assert reconstructed.id == original.id
        assert reconstructed.start == original.start
        assert reconstructed.end == original.end
        assert reconstructed.range == original.range
        assert reconstructed.source == original.source
        assert abs(reconstructed.confidence - original.confidence) < 0.001
        assert reconstructed.t_group == original.t_group
        assert len(reconstructed.evidence) == len(original.evidence)


class TestDomainModelOperations:
    """Test domain operations like merge and split"""
    
    def test_merge_domains(self):
        """Test merging two domains"""
        domain1 = DomainModel(
            id="d1",
            start=10,
            end=30,
            range="10-30",
            source="blast",
            confidence=0.7,
            t_group="2002.1.1.1"
        )
        
        domain2 = DomainModel(
            id="d2", 
            start=25,
            end=45,
            range="25-45",
            source="hhsearch",
            confidence=0.9,
            h_group="2002.1.1"
        )
        
        merged = domain1.merge_with(domain2)
        
        # Should use boundaries from both domains
        assert merged.start == 10  # min(10, 25)
        assert merged.end == 45    # max(30, 45)
        assert merged.range == "10-45"
        
        # Should use higher confidence and its source
        assert merged.confidence == 0.9
        assert merged.source == "hhsearch"
        
        # Should combine classification (primary wins, but secondary fills gaps)
        assert merged.t_group == "2002.1.1.1"  # From domain1
        assert merged.h_group == "2002.1.1"    # From domain2
        
        # Should combine evidence
        assert len(merged.evidence) == 0  # Original domains had no evidence
        
        # Should have merge note
        assert "Merged from" in merged.notes
    
    def test_merge_domains_with_evidence(self):
        """Test merging domains with evidence"""
        evidence1 = Evidence(type="blast", evalue=1e-5, confidence=0.7)
        evidence2 = Evidence(type="hhsearch", probability=90.0, confidence=0.9)
        
        domain1 = DomainModel(
            id="d1",
            start=10,
            end=30,
            range="10-30",
            evidence=[evidence1]
        )
        
        domain2 = DomainModel(
            id="d2",
            start=25,
            end=45,
            range="25-45",
            evidence=[evidence2]
        )
        
        merged = domain1.merge_with(domain2)
        
        # Should combine all evidence
        assert len(merged.evidence) == 2
        # Evidence from higher confidence domain should come first
        assert merged.evidence[0].type == "hhsearch"  # From domain2 (higher confidence)
        assert merged.evidence[1].type == "blast"     # From domain1
    
    def test_split_domain(self):
        """Test splitting a domain"""
        domain = DomainModel(
            id="original",
            start=10,
            end=50,
            range="10-50",
            source="blast",
            confidence=0.8,
            t_group="2002.1.1.1"
        )
        
        domain1, domain2 = domain.split_at(30)
        
        # Check first part
        assert domain1.id == "original_part1"
        assert domain1.start == 10
        assert domain1.end == 29
        assert domain1.range == "10-29"
        assert domain1.source == "blast"
        assert domain1.confidence == 0.72  # 0.8 * 0.9 (penalty for splitting)
        assert domain1.t_group == "2002.1.1.1"
        assert "Split from" in domain1.notes
        
        # Check second part
        assert domain2.id == "original_part2"
        assert domain2.start == 30
        assert domain2.end == 50
        assert domain2.range == "30-50"
        assert domain2.source == "blast"
        assert domain2.confidence == 0.72  # 0.8 * 0.9 (penalty for splitting)
        assert domain2.t_group == "2002.1.1.1"
        assert "Split from" in domain2.notes
    
    def test_split_domain_invalid_position(self):
        """Test splitting domain at invalid position"""
        domain = DomainModel(
            id="test",
            start=10,
            end=20,
            range="10-20"
        )
        
        # Split position outside range should raise error
        with pytest.raises(ValueError):
            domain.split_at(5)  # Before domain
        
        with pytest.raises(ValueError):
            domain.split_at(25)  # After domain
        
        with pytest.raises(ValueError):
            domain.split_at(10)  # At start
        
        with pytest.raises(ValueError):
            domain.split_at(20)  # At end


class TestDomainModelEdgeCases:
    """Test edge cases and error handling"""
    
    def test_domain_with_zero_size(self):
        """Test domain with single residue"""
        domain = DomainModel(
            id="single",
            start=10,
            end=10,
            range="10"
        )
        
        assert domain.size == 1
        assert domain.get_positions() == {10}
    
    def test_domain_with_invalid_range(self):
        """Test domain with invalid start/end"""
        # End before start should still work (model doesn't validate this)
        domain = DomainModel(
            id="invalid",
            start=20,
            end=10,
            range="20-10"  
        )
        
        # Size calculation will be negative, but that's not the model's job to validate
        assert domain.size == -9  # 10 - 20 + 1
    
    def test_domain_string_representations(self):
        """Test domain string representations"""
        domain = DomainModel(
            id="test_domain",
            start=10,
            end=50,
            range="10-50",
            source="blast",
            confidence=0.8,
            t_group="2002.1.1.1"
        )
        
        str_repr = str(domain)
        assert "test_domain" in str_repr
        assert "10-50" in str_repr
        assert "T-group" in str_repr  # Classification level
        
        repr_str = repr(domain)
        assert "DomainModel" in repr_str
        assert "test_domain" in repr_str
        assert "10-50" in repr_str
        assert "blast" in repr_str
        assert "0.800" in repr_str
    
    def test_domain_with_mixed_evidence_types(self):
        """Test domain with mix of Evidence objects and dictionaries"""
        evidence_obj = Evidence(type="blast", evalue=1e-5, confidence=0.7)
        evidence_dict = {
            "type": "hhsearch",
            "probability": 85.0,
            "confidence": 0.85
        }
        
        domain = DomainModel(
            id="mixed",
            start=1,
            end=100,
            range="1-100",
            evidence=[evidence_obj, evidence_dict]
        )
        
        # Both should be standardized to Evidence objects
        assert len(domain.evidence) == 2
        assert all(isinstance(ev, Evidence) for ev in domain.evidence)
        assert domain.evidence[0].type == "blast"
        assert domain.evidence[1].type == "hhsearch"
    
    def test_domain_with_malformed_evidence(self):
        """Test domain with evidence that can't be converted"""
        malformed_evidence = "not a dict or Evidence object"
        
        domain = DomainModel(
            id="test",
            start=1,
            end=100,
            range="1-100",
            evidence=[malformed_evidence]
        )
        
        # Should handle gracefully - might keep original or skip
        # The exact behavior depends on implementation
        assert len(domain.evidence) >= 0  # Shouldn't crash


class TestDomainAlias:
    """Test backward compatibility alias"""
    
    def test_domain_alias(self):
        """Test that Domain is an alias for DomainModel"""
        assert Domain is DomainModel
        
        # Should be able to create using alias
        domain = Domain(
            id="test",
            start=1,
            end=100,
            range="1-100"
        )
        
        assert isinstance(domain, DomainModel)
        assert domain.id == "test"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
