#!/usr/bin/env python3
"""
Simple Integration Test for Domain Partition Pipeline

A simplified version that tests the integration without complex mocking.
"""

import pytest
import tempfile
import json
from pathlib import Path
from datetime import datetime

# Test that we can import the main components
def test_basic_integration_imports():
    """Test that all integration components can be imported"""
    try:
        from ecod.models.pipeline.evidence import Evidence
        from ecod.models.pipeline.domain import DomainModel  
        from ecod.models.pipeline.partition import DomainPartitionResult
        assert True, "All imports successful"
    except ImportError as e:
        pytest.fail(f"Import failed: {e}")

def test_evidence_creation_and_serialization():
    """Test creating evidence and serializing it"""
    from ecod.models.pipeline.evidence import Evidence
    
    # Create test evidence
    evidence = Evidence(
        type="hhsearch",
        source_id="test_source",
        domain_id="test_domain",
        probability=85.0,
        evalue=1e-10,
        query_range="1-100",
        hit_range="1-100"
    )
    
    assert evidence.type == "hhsearch"
    assert evidence.confidence > 0.0
    
    # Test XML serialization
    xml_element = evidence.to_xml()
    assert xml_element is not None

def test_domain_model_creation():
    """Test creating domain models"""
    from ecod.models.pipeline.domain import DomainModel
    from ecod.models.pipeline.evidence import Evidence
    
    # Create domain with evidence
    domain = DomainModel(
        id="test_domain_1",
        start=1,
        end=100,
        range="1-100",
        source="hhsearch"
    )
    
    # Add evidence
    evidence = Evidence(
        type="hhsearch",
        source_id="test_evidence",
        probability=90.0,
        evalue=1e-15
    )
    
    domain.add_evidence(evidence)
    assert len(domain.evidence) == 1
    assert domain.confidence > 0.0

def test_partition_result_creation():
    """Test creating partition results"""
    from ecod.models.pipeline.partition import DomainPartitionResult
    from ecod.models.pipeline.domain import DomainModel
    
    # Create partition result
    result = DomainPartitionResult(
        pdb_id="1test",
        chain_id="A",
        reference="develop291",
        sequence_length=200,
        is_classified=True
    )
    
    # Add domain
    domain = DomainModel(
        id="1test_A_d1",
        start=1,
        end=200,
        range="1-200",
        source="hhsearch",
        confidence=0.9
    )
    
    result.add_domain(domain)
    assert len(result.domains) == 1
    assert result.success

@pytest.mark.slow
def test_xml_serialization_integration():
    """Test full XML serialization workflow"""
    from ecod.models.pipeline.partition import DomainPartitionResult
    from ecod.models.pipeline.domain import DomainModel
    from ecod.models.pipeline.evidence import Evidence
    import xml.etree.ElementTree as ET
    
    # Create complete test result
    result = DomainPartitionResult(
        pdb_id="1xml",
        chain_id="A",
        reference="develop291",
        sequence_length=200,
        is_classified=True
    )
    
    # Create domain with evidence
    domain = DomainModel(
        id="1xml_A_d1",
        start=1,
        end=200,
        range="1-200",
        source="hhsearch",
        confidence=0.95
    )
    
    evidence = Evidence(
        type="hhsearch",
        source_id="test_hit",
        probability=98.0,
        evalue=1e-20,
        query_range="1-200",
        hit_range="1-200"
    )
    
    domain.add_evidence(evidence)
    result.add_domain(domain)
    
    # Test XML round-trip
    xml_element = result.to_xml()
    xml_str = ET.tostring(xml_element, encoding='unicode')
    
    # Parse back
    parsed_element = ET.fromstring(xml_str)
    reconstructed = DomainPartitionResult.from_xml(parsed_element)
    
    assert reconstructed.pdb_id == result.pdb_id
    assert reconstructed.chain_id == result.chain_id
    assert len(reconstructed.domains) == len(result.domains)

# Performance test
@pytest.mark.performance 
def test_bulk_evidence_processing():
    """Test processing many evidence objects"""
    from ecod.models.pipeline.evidence import Evidence
    
    evidences = []
    for i in range(100):
        evidence = Evidence(
            type="hhsearch",
            source_id=f"source_{i}",
            probability=90.0 - i*0.1,
            evalue=10**(-15 + i*0.1)
        )
        evidences.append(evidence)
    
    assert len(evidences) == 100
    assert all(e.confidence > 0 for e in evidences)
