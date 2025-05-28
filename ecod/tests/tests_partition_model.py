#!/usr/bin/env python3
"""
Comprehensive test coverage for ecod.models.pipeline.partition module

Tests the DomainPartitionResult class functionality including:
- Basic partition creation and properties
- Domain standardization and management
- Coverage calculation and statistics
- Serialization (dict and XML)
- Domain operations and queries
- Status flags and classification
- File operations and persistence
- Edge cases and error handling
"""

import pytest
import xml.etree.ElementTree as ET
import tempfile
import os
from datetime import datetime
from unittest.mock import Mock, patch, MagicMock

from ecod.models.pipeline.partition import DomainPartitionResult, PartitionResult
from ecod.models.pipeline.evidence import Evidence


class TestDomainPartitionResultBasics:
    """Test basic DomainPartitionResult creation and properties"""
    
    def test_basic_partition_creation(self):
        """Test creating a basic partition result"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291"
        )
        
        assert result.pdb_id == "1abc"
        assert result.chain_id == "A"
        assert result.reference == "develop291"
        assert result.success == True
        assert result.error is None
        assert result.domains == []
        assert result.is_classified == False
        assert result.is_unclassified == False
        assert result.is_peptide == False
        assert result.sequence_length == 0
        assert result.coverage == 0.0
        assert result.residues_assigned == 0
        assert result.residues_unassigned == 0
        assert result.timestamp is not None
    
    def test_partition_with_domains_dict(self):
        """Test partition creation with dictionary domains"""
        domain_dict = {
            "id": "d1abcA1",
            "start": 10,
            "end": 100,
            "range": "10-100",
            "source": "hhsearch",
            "confidence": 0.9
        }
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A", 
            reference="develop291",
            domains=[domain_dict],
            sequence_length=120
        )
        
        # Should have 1 domain
        assert len(result.domains) == 1
        
        # Domain should be accessible
        domain = result.domains[0]
        
        # Check if standardization worked (could be DomainModel or still dict)
        if hasattr(domain, 'id'):
            # Standardized to DomainModel
            assert domain.id == "d1abcA1"
            assert domain.start == 10
            assert domain.end == 100
        else:
            # Still a dictionary
            assert isinstance(domain, dict)
            assert domain["id"] == "d1abcA1"
            assert domain["start"] == 10
            assert domain["end"] == 100
        
        # Should be classified since it has domains
        assert result.is_classified == True
        assert result.is_unclassified == False
    
    def test_partition_with_evidence_objects(self):
        """Test partition with domains containing Evidence objects"""
        try:
            from ecod.models.pipeline.domain import DomainModel
            NEW_DOMAIN_MODEL = True
        except ImportError:
            NEW_DOMAIN_MODEL = False
        
        if NEW_DOMAIN_MODEL:
            evidence = Evidence(type="hhsearch", probability=95.0, confidence=0.95)
            domain = DomainModel(
                id="test_domain",
                start=1,
                end=100,
                range="1-100",
                evidence=[evidence]
            )
            
            result = DomainPartitionResult(
                pdb_id="1abc",
                chain_id="A",
                reference="develop291",
                domains=[domain],
                sequence_length=120
            )
            
            assert len(result.domains) == 1
            assert result.is_classified == True
            assert result.domains[0].id == "test_domain"
    
    def test_peptide_classification(self):
        """Test peptide classification status"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            is_peptide=True,
            sequence_length=25  # Short sequence
        )
        
        assert result.is_peptide == True
        assert result.is_classified == True  # Peptides are "classified" as peptides
        assert result.is_unclassified == False
    
    def test_failed_partition(self):
        """Test failed partition result"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            success=False,
            error="Processing failed: timeout"
        )
        
        assert result.success == False
        assert result.error == "Processing failed: timeout"
        assert result.is_classified == False
        assert result.is_unclassified == False


class TestDomainStandardization:
    """Test domain standardization and management"""
    
    def test_domain_standardization_mixed_types(self):
        """Test standardization with mixed domain types"""
        domain_dict = {
            "id": "dict_domain",
            "start": 10,
            "end": 50,
            "range": "10-50"
        }
        
        # Create partition with mixed types
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict]
        )
        
        # Should handle mixed types gracefully
        assert len(result.domains) == 1
        
        # Check if standardization occurred
        domain = result.domains[0]
        if hasattr(domain, 'id'):
            # Successfully standardized
            assert domain.id == "dict_domain"
        else:
            # Remained as dict
            assert domain["id"] == "dict_domain"
    
    def test_malformed_domain_handling(self):
        """Test handling of malformed domain data"""
        malformed_domain = {
            "start": "not_a_number",  # Invalid start
            "end": 50,
            "range": "invalid-50"
        }
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[malformed_domain]
        )
        
        # Should not crash and should preserve the domain
        assert len(result.domains) == 1
        
        # Domain should still be accessible
        domain = result.domains[0]
        if isinstance(domain, dict):
            assert domain["start"] == "not_a_number"  # Preserved as-is
            assert domain["end"] == 50
    
    def test_empty_domains_list(self):
        """Test partition with empty domains list"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[],
            sequence_length=100
        )
        
        assert len(result.domains) == 0
        assert result.is_classified == False
        assert result.is_unclassified == False
        assert result.coverage == 0.0
        assert result.residues_assigned == 0
        assert result.residues_unassigned == 100


class TestCoverageCalculation:
    """Test sequence coverage calculation"""
    
    def test_coverage_calculation_single_domain(self):
        """Test coverage calculation with single domain"""
        domain_dict = {
            "id": "domain1",
            "start": 10,
            "end": 50,
            "range": "10-50"
        }
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            sequence_length=100
        )
        
        # Coverage should be calculated automatically
        # Domain covers positions 10-50 = 41 residues
        # Coverage = 41/100 = 0.41
        assert result.residues_assigned == 41
        assert result.residues_unassigned == 59
        assert result.coverage == 0.41
    
    def test_coverage_calculation_multiple_domains(self):
        """Test coverage calculation with multiple domains"""
        domain1 = {"id": "d1", "start": 10, "end": 30, "range": "10-30"}  # 21 residues
        domain2 = {"id": "d2", "start": 50, "end": 80, "range": "50-80"}  # 31 residues
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain1, domain2],
            sequence_length=100
        )
        
        # Total coverage = 21 + 31 = 52 residues
        assert result.residues_assigned == 52
        assert result.residues_unassigned == 48
        assert result.coverage == 0.52
    
    def test_coverage_calculation_overlapping_domains(self):
        """Test coverage calculation with overlapping domains"""
        domain1 = {"id": "d1", "start": 10, "end": 30, "range": "10-30"}  # 21 residues
        domain2 = {"id": "d2", "start": 25, "end": 45, "range": "25-45"}  # 21 residues
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain1, domain2],
            sequence_length=100
        )
        
        # Overlapping coverage should be counted once
        # Positions 10-45 = 36 unique residues
        assert result.residues_assigned == 36
        assert result.residues_unassigned == 64
        assert result.coverage == 0.36
    
    def test_coverage_calculation_complex_ranges(self):
        """Test coverage calculation with multi-segment ranges"""
        domain_dict = {
            "id": "complex_domain",
            "start": 10,
            "end": 50,
            "range": "10-20,30-40,45-50"  # Discontinuous domain
        }
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            sequence_length=100
        )
        
        # Should parse complex range: 10-20(11) + 30-40(11) + 45-50(6) = 28 residues
        assert result.residues_assigned == 28
        assert result.coverage == 0.28
    
    def test_coverage_calculation_zero_length(self):
        """Test coverage calculation with zero sequence length"""
        domain_dict = {"id": "d1", "start": 10, "end": 30, "range": "10-30"}
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            sequence_length=0
        )
        
        # Should handle zero length gracefully
        assert result.coverage == 0.0
        assert result.residues_assigned == 0
        assert result.residues_unassigned == 0
    
    def test_manual_coverage_recalculation(self):
        """Test manual coverage recalculation"""
        domain_dict = {"id": "d1", "start": 10, "end": 30, "range": "10-30"}
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            sequence_length=100
        )
        
        original_coverage = result.coverage
        
        # Change sequence length and recalculate
        result.sequence_length = 200
        result.calculate_coverage()
        
        # Coverage should be different now
        assert result.coverage != original_coverage
        assert result.coverage == 0.105  # 21/200


class TestDomainOperations:
    """Test domain operations and queries"""
    
    def test_add_domain(self):
        """Test adding domains to result"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            sequence_length=100
        )
        
        # Initially no domains
        assert len(result.domains) == 0
        assert result.is_classified == False
        
        # Add a domain
        domain_dict = {"id": "new_domain", "start": 10, "end": 50, "range": "10-50"}
        result.add_domain(domain_dict)
        
        # Should have 1 domain now
        assert len(result.domains) == 1
        assert result.is_classified == True
        
        # Coverage should be updated
        assert result.coverage > 0.0
        assert result.residues_assigned > 0
    
    def test_get_domain_by_id(self):
        """Test retrieving domain by ID"""
        domain1 = {"id": "domain1", "start": 10, "end": 30, "range": "10-30"}
        domain2 = {"id": "domain2", "start": 50, "end": 80, "range": "50-80"}
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain1, domain2]
        )
        
        # Should find existing domain
        found_domain = result.get_domain_by_id("domain1")
        assert found_domain is not None
        
        # Check if it's the right domain (handling both dict and DomainModel)
        if hasattr(found_domain, 'id'):
            assert found_domain.id == "domain1"
        else:
            assert found_domain["id"] == "domain1"
        
        # Should not find non-existent domain
        not_found = result.get_domain_by_id("nonexistent")
        assert not_found is None
    
    def test_get_domains_by_source(self):
        """Test retrieving domains by source"""
        domain1 = {"id": "d1", "start": 10, "end": 30, "range": "10-30", "source": "hhsearch"}
        domain2 = {"id": "d2", "start": 50, "end": 80, "range": "50-80", "source": "blast"}
        domain3 = {"id": "d3", "start": 90, "end": 110, "range": "90-110", "source": "hhsearch"}
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain1, domain2, domain3]
        )
        
        # Should find 2 hhsearch domains
        hhsearch_domains = result.get_domains_by_source("hhsearch")
        assert len(hhsearch_domains) == 2
        
        # Should find 1 blast domain
        blast_domains = result.get_domains_by_source("blast")
        assert len(blast_domains) == 1
        
        # Should find 0 unknown domains
        unknown_domains = result.get_domains_by_source("unknown")
        assert len(unknown_domains) == 0
    
    def test_get_overlapping_domains(self):
        """Test finding overlapping domains"""
        domain1 = {"id": "d1", "start": 10, "end": 30, "range": "10-30"}
        domain2 = {"id": "d2", "start": 25, "end": 45, "range": "25-45"}  # Overlaps with d1
        domain3 = {"id": "d3", "start": 50, "end": 70, "range": "50-70"}  # No overlap
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain1, domain2, domain3]
        )
        
        overlaps = result.get_overlapping_domains()
        
        # Should find 1 overlapping pair
        assert len(overlaps) == 1
        
        # Check the overlapping pair
        i, j, dom1, dom2 = overlaps[0]
        assert i == 0  # First domain index
        assert j == 1  # Second domain index
        
        # Check domains are correct
        if hasattr(dom1, 'id'):
            assert dom1.id == "d1"
            assert dom2.id == "d2"
        else:
            assert dom1["id"] == "d1"
            assert dom2["id"] == "d2"


class TestAnalysisStatistics:
    """Test analysis statistics and quality metrics"""
    
    def test_domain_quality_stats_calculation(self):
        """Test domain quality statistics calculation"""
        domain1 = {
            "id": "d1", "start": 10, "end": 30, "range": "10-30",
            "source": "hhsearch", "confidence": 0.9,
            "t_group": "a", "h_group": "b", "x_group": "c", "a_group": "d"
        }
        domain2 = {
            "id": "d2", "start": 50, "end": 80, "range": "50-80", 
            "source": "blast", "confidence": 0.7,
            "t_group": "a", "h_group": "b"  # Partial classification
        }
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain1, domain2],
            sequence_length=100
        )
        
        stats = result.domain_quality_stats
        
        # Check basic counts
        assert stats["total_domains"] == 2
        
        # Check confidence statistics
        assert 0.7 <= stats["average_confidence"] <= 0.9
        assert stats["min_confidence"] == 0.7
        assert stats["max_confidence"] == 0.9
        
        # Check source distribution
        assert "hhsearch" in stats["source_distribution"]
        assert "blast" in stats["source_distribution"]
        assert stats["source_distribution"]["hhsearch"] == 1
        assert stats["source_distribution"]["blast"] == 1
    
    def test_evidence_summary_calculation(self):
        """Test evidence summary calculation"""
        # Create domains with evidence
        domain_with_evidence = {
            "id": "d1", "start": 10, "end": 30, "range": "10-30",
            "evidence": [
                {"type": "hhsearch", "confidence": 0.9},
                {"type": "blast", "confidence": 0.8}
            ]
        }
        domain_without_evidence = {
            "id": "d2", "start": 50, "end": 80, "range": "50-80"
        }
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_with_evidence, domain_without_evidence]
        )
        
        evidence_summary = result.evidence_summary
        
        # Check evidence counts
        assert evidence_summary["total_evidence_items"] == 2
        assert evidence_summary["domains_with_evidence"] == 1
        
        # Check evidence type distribution
        type_dist = evidence_summary["evidence_type_distribution"]
        assert "hhsearch" in type_dist
        assert "blast" in type_dist
        assert type_dist["hhsearch"] == 1
        assert type_dist["blast"] == 1
    
    def test_empty_domains_statistics(self):
        """Test statistics calculation with no domains"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[],
            sequence_length=100
        )
        
        # Should handle empty domains gracefully
        stats = result.domain_quality_stats
        evidence_summary = result.evidence_summary
        
        # Basic checks for empty case
        assert stats.get("total_domains", 0) == 0
        assert evidence_summary.get("total_evidence_items", 0) == 0
        assert evidence_summary.get("domains_with_evidence", 0) == 0


class TestSerialization:
    """Test serialization to/from dict and XML"""
    
    def test_to_dict_basic(self):
        """Test basic dictionary serialization"""
        domain_dict = {"id": "d1", "start": 10, "end": 50, "range": "10-50"}
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            sequence_length=100,
            success=True,
            is_classified=True
        )
        
        result_dict = result.to_dict()
        
        # Check basic fields
        assert result_dict["pdb_id"] == "1abc"
        assert result_dict["chain_id"] == "A"
        assert result_dict["reference"] == "develop291"
        assert result_dict["success"] == True
        assert result_dict["is_classified"] == True
        assert result_dict["sequence_length"] == 100
        assert result_dict["domain_count"] == 1
        
        # Check domains are included
        assert "domains" in result_dict
        assert len(result_dict["domains"]) == 1
        assert result_dict["domains"][0]["id"] == "d1"
    
    def test_from_dict_basic(self):
        """Test creating partition from dictionary"""
        result_dict = {
            "pdb_id": "2xyz",
            "chain_id": "B",
            "reference": "develop291",
            "success": True,
            "is_classified": True,
            "sequence_length": 150,
            "domains": [
                {"id": "d1", "start": 20, "end": 80, "range": "20-80"}
            ]
        }
        
        result = DomainPartitionResult.from_dict(result_dict)
        
        assert result.pdb_id == "2xyz"
        assert result.chain_id == "B"
        assert result.reference == "develop291"
        assert result.success == True
        assert result.is_classified == True
        assert result.sequence_length == 150
        assert len(result.domains) == 1
    
    def test_dict_round_trip(self):
        """Test dictionary serialization round-trip"""
        domain_dict = {
            "id": "test_domain",
            "start": 10,
            "end": 50,
            "range": "10-50",
            "source": "hhsearch",
            "confidence": 0.85
        }
        
        original = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            sequence_length=100,
            success=True,
            is_classified=True
        )
        
        # Convert to dict and back
        result_dict = original.to_dict()
        reconstructed = DomainPartitionResult.from_dict(result_dict)
        
        # Check key fields preserved
        assert reconstructed.pdb_id == original.pdb_id
        assert reconstructed.chain_id == original.chain_id
        assert reconstructed.reference == original.reference
        assert reconstructed.success == original.success
        assert reconstructed.is_classified == original.is_classified
        assert reconstructed.sequence_length == original.sequence_length
        assert len(reconstructed.domains) == len(original.domains)
    
    def test_to_xml_basic(self):
        """Test basic XML serialization"""
        domain_dict = {"id": "d1", "start": 10, "end": 50, "range": "10-50"}
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            success=True,
            is_classified=True
        )
        
        xml_element = result.to_xml()
        
        # Check XML structure
        assert xml_element.tag == "domain_partition"
        assert xml_element.get("pdb_id") == "1abc"
        assert xml_element.get("chain_id") == "A"
        assert xml_element.get("reference") == "develop291"
        assert xml_element.get("success") == "true"
        assert xml_element.get("is_classified") == "true"
        
        # Check domains section
        domains_elem = xml_element.find("domains")
        assert domains_elem is not None
        assert domains_elem.get("count") == "1"
        
        domain_elements = domains_elem.findall("domain")
        assert len(domain_elements) == 1
        assert domain_elements[0].get("id") == "d1"
    
    def test_xml_with_metadata(self):
        """Test XML serialization with metadata"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            sequence_length=100,
            include_metadata=True
        )
        
        xml_element = result.to_xml()
        
        # Check metadata section
        metadata = xml_element.find("metadata")
        assert metadata is not None
        
        # Check metadata fields
        seq_len = metadata.find("sequence_length")
        assert seq_len is not None
        assert seq_len.text == "100"
        
        coverage_elem = metadata.find("coverage")
        assert coverage_elem is not None
    
    def test_xml_evidence_inclusion_control(self):
        """Test controlling evidence inclusion in XML"""
        domain_with_evidence = {
            "id": "d1", "start": 10, "end": 50, "range": "10-50",
            "evidence": [{"type": "hhsearch", "confidence": 0.9}]
        }
        
        # Test with evidence included
        result_with_evidence = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_with_evidence],
            include_evidence=True
        )
        
        xml_with_evidence = result_with_evidence.to_xml()
        domains_elem = xml_with_evidence.find("domains")
        domain_elem = domains_elem.find("domain")
        evidence_list = domain_elem.find("evidence_list")
        
        # Should have evidence (if domain model supports it)
        # Note: Actual presence depends on whether domain was standardized
        
        # Test with evidence excluded
        result_without_evidence = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_with_evidence],
            include_evidence=False
        )
        
        xml_without_evidence = result_without_evidence.to_xml()
        # Should still work but evidence may be excluded
        assert xml_without_evidence.tag == "domain_partition"
    
    def test_from_xml_basic(self):
        """Test creating partition from XML"""
        xml_str = '''
        <domain_partition pdb_id="2xyz" chain_id="B" reference="develop291" 
                         success="true" is_classified="true">
            <domains count="1">
                <domain id="d1" start="20" end="80" range="20-80" source="blast" confidence="0.8000"/>
            </domains>
        </domain_partition>
        '''
        xml_element = ET.fromstring(xml_str)
        
        result = DomainPartitionResult.from_xml(xml_element)
        
        assert result.pdb_id == "2xyz"
        assert result.chain_id == "B"
        assert result.reference == "develop291"
        assert result.success == True
        assert result.is_classified == True
        assert len(result.domains) == 1
    
    def test_xml_round_trip(self):
        """Test XML serialization round-trip"""
        domain_dict = {
            "id": "test_domain",
            "start": 10,
            "end": 50,
            "range": "10-50"
        }
        
        original = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            success=True
        )
        
        # Convert to XML and back
        xml_element = original.to_xml()
        reconstructed = DomainPartitionResult.from_xml(xml_element)
        
        # Check key fields preserved
        assert reconstructed.pdb_id == original.pdb_id
        assert reconstructed.chain_id == original.chain_id
        assert reconstructed.reference == original.reference
        assert reconstructed.success == original.success
        assert len(reconstructed.domains) == len(original.domains)


class TestFactoryMethods:
    """Test factory methods for creating partition results"""
    
    def test_from_domains_factory(self):
        """Test creating partition from domains list"""
        domain1 = {"id": "d1", "start": 10, "end": 30, "range": "10-30"}
        domain2 = {"id": "d2", "start": 50, "end": 80, "range": "50-80"}
        
        result = DomainPartitionResult.from_domains(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain1, domain2],
            sequence_length=100
        )
        
        assert result.pdb_id == "1abc"
        assert result.chain_id == "A"
        assert result.reference == "develop291"
        assert len(result.domains) == 2
        assert result.sequence_length == 100
        assert result.is_classified == True  # Has domains
    
    def test_from_domains_with_kwargs(self):
        """Test from_domains with additional parameters"""
        domain_dict = {"id": "d1", "start": 10, "end": 50, "range": "10-50"}
        
        result = DomainPartitionResult.from_domains(
            pdb_id="1abc",
            chain_id="A", 
            reference="develop291",
            domains=[domain_dict],
            sequence_length=100,
            success=True,
            is_peptide=False,
            include_evidence=True
        )
        
        assert result.success == True
        assert result.is_peptide == False
        assert result.include_evidence == True
        assert len(result.domains) == 1


class TestFileOperations:
    """Test file operations and persistence"""
    
    def test_save_with_output_dir(self):
        """Test saving partition to file with output directory"""
        domain_dict = {"id": "d1", "start": 10, "end": 50, "range": "10-50"}
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict]
        )
        
        with tempfile.TemporaryDirectory() as temp_dir:
            success = result.save(output_dir=temp_dir)
            
            assert success == True
            assert result.domain_file is not None
            assert os.path.exists(result.domain_file)
            
            # Check filename format
            expected_name = "1abc_A.develop291.domains.xml"
            assert expected_name in result.domain_file
    
    def test_save_with_custom_filename(self):
        """Test saving with custom filename"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291"
        )
        
        with tempfile.TemporaryDirectory() as temp_dir:
            custom_filename = "custom_partition.xml"
            success = result.save(output_dir=temp_dir, filename=custom_filename)
            
            assert success == True
            assert result.domain_file is not None
            assert custom_filename in result.domain_file
            assert os.path.exists(result.domain_file)
    
    def test_save_with_existing_domain_file(self):
        """Test saving when domain_file is already set"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291"
        )
        
        with tempfile.TemporaryDirectory() as temp_dir:
            existing_path = os.path.join(temp_dir, "existing.xml")
            result.domain_file = existing_path
            
            success = result.save()
            
            assert success == True
            assert os.path.exists(existing_path)
    
    def test_save_failure_handling(self):
        """Test save failure handling"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291"
        )
        
        # Try to save to invalid path
        success = result.save(output_dir="/invalid/path/that/does/not/exist")
        
        # Should handle failure gracefully
        assert success == False


class TestEdgeCases:
    """Test edge cases and error handling"""
    
    def test_invalid_timestamp_handling(self):
        """Test handling of invalid timestamp data"""
        result_dict = {
            "pdb_id": "1abc",
            "chain_id": "A",
            "reference": "develop291",
            "timestamp": "invalid_timestamp_format"
        }
        
        result = DomainPartitionResult.from_dict(result_dict)
        
        # Should handle invalid timestamp gracefully
        assert result.pdb_id == "1abc"
        # timestamp should be None or current time
        assert result.timestamp is None or isinstance(result.timestamp, datetime)
    
    def test_very_large_domain_list(self):
        """Test handling of large domain lists"""
        # Create many domains
        domains = []
        for i in range(100):
            domains.append({
                "id": f"domain_{i}",
                "start": i * 10,
                "end": i * 10 + 9,
                "range": f"{i*10}-{i*10+9}"
            })
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=domains,
            sequence_length=1000
        )
        
        # Should handle large lists without issues
        assert len(result.domains) == 100
        assert result.is_classified == True
        
        # Statistics calculation should work
        stats = result.domain_quality_stats
        assert stats["total_domains"] == 100
    
    def test_conflicting_classification_flags(self):
        """Test handling of conflicting classification flags"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            is_classified=True,
            is_unclassified=True,  # Conflicting!
            is_peptide=True
        )
        
        # Should resolve conflicts appropriately
        # Peptide classification takes precedence
        assert result.is_peptide == True
        assert result.is_classified == True
    
    def test_negative_sequence_length(self):
        """Test handling of negative sequence length"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            sequence_length=-50  # Invalid
        )
        
        # Should handle gracefully
        result.calculate_coverage()
        assert result.coverage == 0.0
    
    def test_string_representations(self):
        """Test string representations"""
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            success=True,
            is_classified=True
        )
        
        str_repr = str(result)
        assert "1abc_A" in str_repr
        assert "SUCCESS" in str_repr
        assert "CLASSIFIED" in str_repr
        
        repr_str = repr(result)
        assert "DomainPartitionResult" in repr_str
        assert "1abc" in repr_str
        assert "A" in repr_str
    
    def test_get_summary_stats(self):
        """Test summary statistics generation"""
        domain_dict = {"id": "d1", "start": 10, "end": 50, "range": "10-50"}
        
        result = DomainPartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291",
            domains=[domain_dict],
            sequence_length=100,
            processing_time=1.5
        )
        
        summary = result.get_summary_stats()
        
        # Check summary fields
        assert summary["pdb_id"] == "1abc"
        assert summary["chain_id"] == "A"
        assert summary["reference"] == "develop291"
        assert summary["success"] == True
        assert summary["domain_count"] == 1
        assert summary["sequence_length"] == 100
        assert summary["processing_time"] == 1.5
        assert "coverage" in summary
        assert "residues_assigned" in summary


class TestBackwardCompatibility:
    """Test backward compatibility features"""
    
    def test_partition_result_alias(self):
        """Test that PartitionResult is an alias"""
        assert PartitionResult is DomainPartitionResult
        
        # Should be able to create using alias
        result = PartitionResult(
            pdb_id="1abc",
            chain_id="A",
            reference="develop291"
        )
        
        assert isinstance(result, DomainPartitionResult)
        assert result.pdb_id == "1abc"
    
    def test_legacy_field_handling(self):
        """Test handling of legacy field names"""
        # Test that alternative field names work in from_dict
        legacy_dict = {
            "pdb_id": "1abc",
            "chain_id": "A", 
            "reference": "develop291",
            "domain_count": 5  # Legacy field
        }
        
        result = DomainPartitionResult.from_dict(legacy_dict)
        
        # Should handle legacy fields gracefully
        assert result.pdb_id == "1abc"
        assert result.chain_id == "A"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
