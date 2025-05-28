#!/usr/bin/env python3
"""
Test configuration and fixtures for ECOD pipeline testing

This conftest.py provides clean, focused fixtures for testing the ECOD bioinformatics pipeline.
Built after comprehensive testing methodology development.
"""

import pytest
import tempfile
import shutil
import os
from pathlib import Path
import xml.etree.ElementTree as ET


@pytest.fixture(scope="session")
def project_root():
    """Get the project root directory"""
    return Path(__file__).parent.parent.parent


@pytest.fixture
def temp_test_dir():
    """Create a temporary directory for test files"""
    temp_dir = tempfile.mkdtemp(prefix="ecod_test_")
    yield temp_dir
    shutil.rmtree(temp_dir)


class TestDataFactory:
    """Factory for creating test data objects"""
    
    @staticmethod
    def create_blast_xml_element(**kwargs):
        """Create BLAST XML element for testing"""
        defaults = {
            "num": "1",
            "domain_id": "d1abcA1",
            "pdb_id": "1abc",
            "chain_id": "A",
            "evalues": "1e-5",
            "hsp_count": "1"
        }
        defaults.update(kwargs)
        
        element = ET.Element("hit")
        for key, value in defaults.items():
            if value is not None:
                element.set(key, str(value))
        
        # Add ranges if provided
        if "query_range" in kwargs:
            query_reg = ET.SubElement(element, "query_reg")
            query_reg.text = kwargs["query_range"]
        
        if "hit_range" in kwargs:
            hit_reg = ET.SubElement(element, "hit_reg")  
            hit_reg.text = kwargs["hit_range"]
            
        return element
    
    @staticmethod
    def create_hhsearch_xml_element(**kwargs):
        """Create HHSearch XML element for testing"""
        defaults = {
            "hit_id": "d1abcA1",
            "domain_id": "d1abcA1",
            "probability": "90.0",
            "evalue": "1e-6", 
            "score": "50.0",
            "num": "1"
        }
        defaults.update(kwargs)
        
        element = ET.Element("hit")
        for key, value in defaults.items():
            if value is not None:
                element.set(key, str(value))
        
        # Add ranges if provided
        if "query_range" in kwargs:
            query_reg = ET.SubElement(element, "query_reg")
            query_reg.text = kwargs["query_range"]
            
        if "hit_range" in kwargs:
            hit_reg = ET.SubElement(element, "hit_reg")
            hit_reg.text = kwargs["hit_range"]
            
        return element


@pytest.fixture
def test_data_factory():
    """Provide the test data factory"""
    return TestDataFactory


@pytest.fixture
def sample_evidence_data():
    """Sample evidence data for testing"""
    return {
        "blast_excellent": {
            "type": "domain_blast",
            "source_id": "d1abcA1",
            "domain_id": "d1abcA1",
            "evalue": 1e-10,
            "query_range": "10-50",
            "hit_range": "5-45",
            "identity": 85.0,
            "coverage": 90.0
        },
        "hhsearch_excellent": {
            "type": "hhsearch",
            "source_id": "d1abcA1", 
            "domain_id": "d1abcA1",
            "probability": 95.0,
            "evalue": 1e-8,
            "score": 55.0,
            "query_range": "10-50",
            "hit_range": "5-45"
        },
        "blast_poor": {
            "type": "domain_blast",
            "source_id": "d2xyzB1",
            "domain_id": "d2xyzB1", 
            "evalue": 1.0,
            "identity": 25.0,
            "coverage": 40.0
        },
        "chain_blast": {
            "type": "chain_blast",
            "source_id": "1xyz_B",
            "evalue": 1e-5,
            "hsp_count": 3
        }
    }


@pytest.fixture
def confidence_test_cases():
    """Standard test cases for confidence calculation validation"""
    return [
        # (evidence_params, expected_min_confidence, description)
        ({"type": "blast", "evalue": 1e-10}, 0.8, "Highly significant BLAST"),
        ({"type": "blast", "evalue": 1e-5}, 0.6, "Significant BLAST"),  
        ({"type": "blast", "evalue": 0.01}, 0.4, "Marginally significant BLAST"),
        ({"type": "blast", "evalue": 1.0}, 0.3, "Moderate BLAST"),
        
        ({"type": "hhsearch", "probability": 95.0}, 0.9, "Highly confident HHSearch"),
        ({"type": "hhsearch", "probability": 80.0}, 0.7, "Confident HHSearch"),
        ({"type": "hhsearch", "probability": 60.0}, 0.5, "Moderate HHSearch"),
    ]


@pytest.fixture
def serialization_test_cases():
    """Standard test cases for serialization validation"""
    return {
        "complete_evidence": {
            "type": "hhsearch",
            "source_id": "d1abcA1",
            "domain_id": "d1abcA1",
            "query_range": "10-50",
            "hit_range": "5-45", 
            "probability": 85.5,
            "evalue": 1e-8,
            "score": 42.3,
            "confidence": 0.9,
            "identity": 78.5,
            "coverage": 65.2,
            "hsp_count": 3,
            "t_group": "a",
            "h_group": "b",
            "extra_attributes": {"custom_field": "test_value"}
        },
        "minimal_evidence": {
            "type": "blast",
            "source_id": "test_hit"
        },
        "edge_case_evidence": {
            "type": "unknown",
            "source_id": "edge_case",
            "evalue": 0.0,  # Zero E-value
            "probability": 101.0,  # Invalid probability
            "confidence": None  # Auto-calculate
        }
    }


# Mock fixtures for external dependencies
@pytest.fixture
def mock_logging():
    """Mock logging to avoid log spam during tests"""
    import logging
    
    # Set logging level to ERROR to suppress debug/info messages
    logging.getLogger().setLevel(logging.ERROR)
    
    # Provide a logger for tests that need it
    return logging.getLogger('test')


@pytest.fixture
def isolated_imports():
    """Ensure tests don't interfere with each other's imports"""
    import sys
    
    # Save original modules
    original_modules = sys.modules.copy()
    
    yield
    
    # Clean up any modules added during test
    new_modules = set(sys.modules.keys()) - set(original_modules.keys())
    for module in new_modules:
        if module.startswith('ecod'):  # Only clean up our modules
            del sys.modules[module]


# Performance testing fixtures
@pytest.fixture
def performance_timer():
    """Timer for performance testing"""
    import time
    
    times = {}
    
    def timer(name):
        def decorator(func):
            start = time.time()
            result = func()
            end = time.time()
            times[name] = end - start
            return result
        return decorator
    
    timer.times = times
    return timer


# Markers for different test categories
def pytest_configure(config):
    """Configure custom markers"""
    config.addinivalue_line(
        "markers", "unit: mark test as a unit test"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test"
    )
    config.addinivalue_line(
        "markers", "performance: mark test as a performance test"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow (may be skipped)"
    )


# Skip slow tests by default unless explicitly requested
def pytest_collection_modifyitems(config, items):
    """Skip slow tests unless --runslow option is given"""
    if not config.getoption("--runslow"):
        skip_slow = pytest.mark.skip(reason="need --runslow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)


def pytest_addoption(parser):
    """Add custom command line options"""
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )
