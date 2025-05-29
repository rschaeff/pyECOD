#!/usr/bin/env python3
"""
Unit tests for domain evidence integration
Tests the enhanced evidence handling during dictionary to model conversion
"""
import unittest
import os
import sys
import tempfile
from pathlib import Path
import xml.etree.ElementTree as ET

# Add parent directory to path for imports during testing
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.models.domain_analysis.domain_model import DomainModel
from ecod.models.domain_analysis.partition_result import DomainPartitionResult
from ecod.models.evidence import DomainEvidence
from ecod.utils.evidence_bridge import EvidenceBridge


class TestEvidenceIntegration(unittest.TestCase):
    """Test suite for evidence integration in domain analysis"""
    
    def setUp(self):
        """Set up test data"""
        # Create test evidence
        self.hhsearch_evidence = {
            "type": "hhsearch",
            "source_id": "e4xfcA1",
            "domain_id": "e4xfcA1",
            "query_range": "10-150",
            "hit_range": "5-145",
            "probability": 99.8,
            "evalue": 1e-45,
            "score": 425.2
        }
        
        self.blast_evidence = {
            "type": "domain_blast",
            "source_id": "e3a9qA1",
            "domain_id": "e3a9qA1",
            "query_range": "5-120",
            "hit_range": "1-115",
            "evalue": 1e-25,
            "pdb_id": "3a9q",
            "chain_id": "A"
        }
        
        # Create test domain dictionaries
        self.domain_dict1 = {
            "start": 10,
            "end": 150,
            "range": "10-150",
            "source": "hhsearch",
            "confidence": 0.998,
            "source_id": "e4xfcA1",
            "t_group": "2004.1.1",
            "h_group": "2004.1",
            "x_group": "2004",
            "a_group": "a.17",
            "is_manual_rep": False,
            "is_f70": False,
            "is_f40": False,
            "is_f99": False,
            "evidence": [self.hhsearch_evidence]
        }
        
        self.domain_dict2 = {
            "start": 5,
            "end": 120,
            "range": "5-120",
            "source": "domain_blast",
            "confidence": 0.8,
            "source_id": "e3a9qA1",
            "t_group": "2004.1.1",
            "h_group": "2004.1",
            "x_group": "2004",
            "a_group": "a.17",
            "is_manual_rep": False,
            "is_f70": True,
            "is_f40": False,
            "is_f99": False,
            "evidence": [self.blast_evidence]
        }
        
        # Create test DomainModel object
        self.domain_model = DomainModel(
            id="7xya_A_d1",
            start=10,
            end=150,
            range="10-150",
            source="hhsearch",
            confidence=0.998,
            source_id="e4xfcA1",
            t_group="2004.1.1",
            h_group="2004.1",
            x_group="2004",
            a_group="a.17",
            is_manual_rep=False,
            is_f70=False,
            is_f40=False,
            is_f99=False,
            evidence=[DomainEvidence.from_dict(self.hhsearch_evidence)]
        )
    
    def test_evidence_standardization(self):
        """Test evidence standardization from dictionary to model"""
        # Convert dictionary evidence to DomainEvidence
        std_evidence = EvidenceBridge.standardize_evidence(self.hhsearch_evidence)
        
        # Check type
        self.assertIsInstance(std_evidence, DomainEvidence)
        
        # Check attributes
        self.assertEqual(std_evidence.type, "hhsearch")
        self.assertEqual(std_evidence.source_id, "e4xfcA1")
        self.assertEqual(std_evidence.domain_id, "e4xfcA1")
        self.assertEqual(std_evidence.query_range, "10-150")
        self.assertEqual(std_evidence.hit_range, "5-145")
        
        # Check that confidence is normalized to 0-1 range
        self.assertGreaterEqual(std_evidence.confidence, 0.0)
        self.assertLessEqual(std_evidence.confidence, 1.0)
        
        # Check attributes contains other fields
        self.assertIn("probability", std_evidence.attributes)
        self.assertIn("evalue", std_evidence.attributes)
        self.assertIn("score", std_evidence.attributes)
    
    def test_domain_model_conversion(self):
        """Test conversion from dictionary domain to DomainModel"""
        # Convert dictionary domain to DomainModel
        domain_model = EvidenceBridge.ensure_domain_model(
            self.domain_dict1, "7xya", "A", 0
        )
        
        # Check type
        self.assertIsInstance(domain_model, DomainModel)
        
        # Check attributes
        self.assertEqual(domain_model.start, 10)
        self.assertEqual(domain_model.end, 150)
        self.assertEqual(domain_model.range, "10-150")
        self.assertEqual(domain_model.t_group, "2004.1.1")
        
        # Check evidence was converted correctly
        self.assertEqual(len(domain_model.evidence), 1)
        self.assertIsInstance(domain_model.evidence[0], DomainEvidence)
        self.assertEqual(domain_model.evidence[0].type, "hhsearch")
    
    def test_result_with_dict_domains(self):
        """Test DomainPartitionResult with dictionary domains"""
        # Create result with dictionary domains
        result = DomainPartitionResult(
            pdb_id="7xya",
            chain_id="A",
            reference="develop291",
            domains=[self.domain_dict1, self.domain_dict2],
            sequence_length=200
        )
        
        # Check domains were standardized to models
        self.assertEqual(len(result.domains), 2)
        for domain in result.domains:
            self.assertIsInstance(domain, DomainModel)
            
        # Check first domain evidence
        self.assertEqual(len(result.domains[0].evidence), 1)
        self.assertIsInstance(result.domains[0].evidence[0], DomainEvidence)
        
        # Check coverage calculation
        self.assertGreater(result.coverage, 0.0)
        self.assertEqual(result.sequence_length, 200)
    
    def test_result_with_mixed_domains(self):
        """Test DomainPartitionResult with mixed domain types"""
        # Create result with mixed domain types
        result = DomainPartitionResult(
            pdb_id="7xya",
            chain_id="A",
            reference="develop291",
            domains=[self.domain_dict1, self.domain_model],
            sequence_length=200
        )
        
        # Check domains were standardized to models
        self.assertEqual(len(result.domains), 2)
        for domain in result.domains:
            self.assertIsInstance(domain, DomainModel)
            
        # Check evidence in both domains
        for domain in result.domains:
            self.assertEqual(len(domain.evidence), 1)
            self.assertIsInstance(domain.evidence[0], DomainEvidence)
    
    def test_xml_generation_with_evidence(self):
        """Test XML generation with evidence included"""
        # Create result
        result = DomainPartitionResult(
            pdb_id="7xya",
            chain_id="A",
            reference="develop291",
            domains=[self.domain_dict1, self.domain_dict2],
            sequence_length=200,
            include_evidence=True
        )
        
        # Generate XML
        xml = result.to_xml()
        
        # Check root attributes
        self.assertEqual(xml.get("pdb_id"), "7xya")
        self.assertEqual(xml.get("chain_id"), "A")
        self.assertEqual(xml.get("reference"), "develop291")
        
        # Check domains
        domains = xml.findall(".//domain")
        self.assertEqual(len(domains), 2)
        
        # Check evidence included
        evidence_lists = xml.findall(".//evidence_list")
        self.assertEqual(len(evidence_lists), 2)
        
        # Check individual evidence items
        evidence_items = xml.findall(".//evidence_list/evidence")
        self.assertTrue(len(evidence_items) >= 2)
    
    def test_xml_generation_without_evidence(self):
        """Test XML generation with evidence excluded"""
        # Create result
        result = DomainPartitionResult(
            pdb_id="7xya",
            chain_id="A",
            reference="develop291",
            domains=[self.domain_dict1, self.domain_dict2],
            sequence_length=200,
            include_evidence=False
        )
        
        # Generate XML
        xml = result.to_xml()
        
        # Check domains
        domains = xml.findall(".//domain")
        self.assertEqual(len(domains), 2)
        
        # Check evidence excluded
        evidence_lists = xml.findall(".//evidence_list")
        self.assertEqual(len(evidence_lists), 0)
    
    def test_save_to_file(self):
        """Test saving DomainPartitionResult to a file"""
        # Create result
        result = DomainPartitionResult(
            pdb_id="7xya",
            chain_id="A",
            reference="develop291",
            domains=[self.domain_dict1, self.domain_dict2],
            sequence_length=200
        )
        
        # Create temp directory
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save to file
            success = result.save(temp_dir)
            self.assertTrue(success)
            
            # Check file exists
            domain_file = os.path.join(temp_dir, "domains", "7xya_A.develop291.domains.xml")
            self.assertTrue(os.path.exists(domain_file))
            
            # Check file is valid XML
            tree = ET.parse(domain_file)
            root = tree.getroot()
            
            # Check basic structure
            self.assertEqual(root.tag, "domain_partition")
            self.assertEqual(root.get("pdb_id"), "7xya")
            
            # Check domains
            domains = root.findall(".//domain")
            self.assertEqual(len(domains), 2)
            
            # Check evidence
            evidence_items = root.findall(".//evidence_list/evidence")
            self.assertTrue(len(evidence_items) >= 2)
    
    def test_from_xml(self):
        """Test creating DomainPartitionResult from XML"""
        # Create result and save to file
        original_result = DomainPartitionResult(
            pdb_id="7xya",
            chain_id="A",
            reference="develop291",
            domains=[self.domain_dict1, self.domain_dict2],
            sequence_length=200
        )
        
        # Create XML
        xml = original_result.to_xml()
        
        # Create new result from XML
        new_result = DomainPartitionResult.from_xml(xml)
        
        # Check basic attributes
        self.assertEqual(new_result.pdb_id, "7xya")
        self.assertEqual(new_result.chain_id, "A")
        self.assertEqual(new_result.reference, "develop291")
        self.assertEqual(new_result.sequence_length, 200)
        
        # Check domains
        self.assertEqual(len(new_result.domains), 2)
        for domain in new_result.domains:
            self.assertIsInstance(domain, DomainModel)


if __name__ == "__main__":
    unittest.main()
