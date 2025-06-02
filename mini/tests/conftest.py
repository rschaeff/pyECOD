#!/usr/bin/env python3
"""
Test configuration and fixtures for mini_pyecod tests

This file provides shared test fixtures and configuration
for all mini_pyecod tests.
"""

import os
import sys
import pytest
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

@pytest.fixture(scope="session")
def batch_dir():
    """Default batch directory for tests"""
    return "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"

@pytest.fixture(scope="session")
def test_data_dir():
    """Test data directory with reference files"""
    return str(Path(__file__).parent.parent / "test_data")

@pytest.fixture(scope="session")
def sample_protein():
    """Sample protein ID for testing"""
    return "8ovp_A"

@pytest.fixture(scope="session")
def reference_data(test_data_dir):
    """Load reference data once for all tests"""
    from mini.parser import load_reference_lengths, load_protein_lengths
    from mini.decomposer import load_domain_definitions
    
    domain_lengths_file = os.path.join(test_data_dir, "domain_lengths.csv")
    protein_lengths_file = os.path.join(test_data_dir, "protein_lengths.csv")
    domain_definitions_file = os.path.join(test_data_dir, "domain_definitions.csv")
    
    return {
        'domain_lengths': load_reference_lengths(domain_lengths_file) if os.path.exists(domain_lengths_file) else {},
        'protein_lengths': load_protein_lengths(protein_lengths_file) if os.path.exists(protein_lengths_file) else {},
        'domain_definitions': load_domain_definitions(domain_definitions_file) if os.path.exists(domain_definitions_file) else {}
    }

@pytest.fixture
def sample_evidence(batch_dir, sample_protein, reference_data):
    """Sample evidence for testing"""
    from mini.parser import parse_domain_summary
    from mini.blast_parser import load_chain_blast_alignments
    
    # Parse protein ID
    parts = sample_protein.split('_')
    pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else 'A'
    
    # Load BLAST alignments
    blast_dir = os.path.join(batch_dir, "blast/chain")
    blast_alignments = {}
    if os.path.exists(blast_dir):
        blast_alignments = load_chain_blast_alignments(blast_dir, pdb_id, chain_id)
    
    # Parse evidence
    xml_path = os.path.join(batch_dir, "domains", f"{sample_protein}.develop291.domain_summary.xml")
    
    if not os.path.exists(xml_path):
        pytest.skip(f"Domain summary file not found: {xml_path}")
    
    evidence = parse_domain_summary(
        xml_path,
        reference_lengths=reference_data['domain_lengths'],
        protein_lengths=reference_data['protein_lengths'],
        blast_alignments=blast_alignments,
        require_reference_lengths=True
    )
    
    if not evidence:
        pytest.skip("No evidence with reference lengths found")
    
    return evidence

@pytest.fixture
def temp_output_dir(tmp_path):
    """Temporary directory for test outputs"""
    return str(tmp_path)

# Test markers
def pytest_configure(config):
    """Configure custom test markers"""
    config.addinivalue_line("markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')")
    config.addinivalue_line("markers", "integration: marks tests as integration tests")
    config.addinivalue_line("markers", "visualization: marks tests that require visualization tools")

# Skip tests if required files are missing
def pytest_runtest_setup(item):
    """Setup function run before each test"""
    # Check if required test data exists
    test_data_dir = Path(__file__).parent.parent / "test_data"
    
    if not test_data_dir.exists():
        pytest.skip("Test data directory not found - run range_cache_parser.py first")
    
    required_files = [
        test_data_dir / "domain_lengths.csv",
        test_data_dir / "protein_lengths.csv"
    ]
    
    missing_files = [f for f in required_files if not f.exists()]
    if missing_files:
        pytest.skip(f"Required test files missing: {[str(f) for f in missing_files]}")

# Custom test collection
def pytest_collection_modifyitems(config, items):
    """Modify test collection"""
    # Add markers to tests based on their names
    for item in items:
        if "visualization" in item.name or "pymol" in item.name:
            item.add_marker(pytest.mark.visualization)
        
        if "integration" in item.name or "test_case" in item.name:
            item.add_marker(pytest.mark.integration)
        
        if "decomposition" in item.name and "chain_blast" in item.name:
            item.add_marker(pytest.mark.slow)
