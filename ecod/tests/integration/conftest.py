#!/usr/bin/env python3
"""
Integration Test Configuration and Fixtures

Specialized fixtures for integration testing of the ECOD domain partition pipeline.
Complements the main conftest.py with integration-specific functionality.
"""

import pytest
import os
import json
import yaml
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Any, Optional
from unittest.mock import Mock, patch
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT

# Add project root to path
import sys
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from ecod.core.context import ApplicationContext


@pytest.fixture(scope="session")
def integration_test_root():
    """Get the integration test root directory"""
    return Path(__file__).parent


@pytest.fixture(scope="session") 
def test_configs_dir(integration_test_root):
    """Get the test configurations directory"""
    return integration_test_root / "configs"


@pytest.fixture(scope="session")
def test_datasets_dir(integration_test_root):
    """Get the test datasets directory"""
    return integration_test_root / "datasets"


@pytest.fixture(scope="session")
def test_baselines_dir(integration_test_root):
    """Get the test baselines directory"""
    return integration_test_root / "baselines"


@pytest.fixture(scope="session")
def test_results_dir(integration_test_root):
    """Get the test results directory"""
    return integration_test_root / "results"


@pytest.fixture(scope="session")
def test_config(test_configs_dir, request):
    """Load test configuration"""
    # Get config name from command line or use default
    config_name = getattr(request.config.option, 'config', 'default')
    config_file = test_configs_dir / f"{config_name}.yml"
    
    if not config_file.exists():
        config_file = test_configs_dir / "default.yml"
    
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="session")
def test_database_config(test_config):
    """Get database configuration for tests"""
    db_config = test_config['database'].copy()
    
    # Override with environment variables if set
    db_config['host'] = os.getenv('PGHOST', db_config['host'])
    db_config['port'] = int(os.getenv('PGPORT', db_config['port']))
    db_config['user'] = os.getenv('PGUSER', db_config['user'])
    db_config['password'] = os.getenv('PGPASSWORD', db_config['password'])
    db_config['database'] = os.getenv('PGDATABASE', db_config['database'])
    
    return db_config


@pytest.fixture(scope="session")
def test_database_connection(test_database_config):
    """Create test database connection"""
    try:
        # Test connection
        conn = psycopg2.connect(**test_database_config)
        conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        
        # Verify we can query
        with conn.cursor() as cur:
            cur.execute("SELECT 1")
            result = cur.fetchone()
            assert result[0] == 1
        
        yield conn
        
    except psycopg2.Error as e:
        pytest.skip(f"Cannot connect to test database: {e}")
    finally:
        if 'conn' in locals():
            conn.close()


@pytest.fixture
def mock_application_context(test_config, temp_test_dir):
    """Create mock application context for integration tests"""
    context = Mock(spec=ApplicationContext)
    
    # Mock config manager
    context.config_manager = Mock()
    context.config_manager.config = test_config
    context.config_manager.get_db_config.return_value = test_config['database']
    context.config_manager.get_path.return_value = temp_test_dir
    
    # Mock database manager
    context.db_manager = Mock()
    
    # Mock logger
    context.logger = Mock()
    
    return context


@pytest.fixture
def golden_datasets(test_datasets_dir):
    """Load golden datasets"""
    catalog_file = test_datasets_dir / "catalog.json"
    
    if not catalog_file.exists():
        pytest.skip("Golden datasets not available")
    
    with open(catalog_file, 'r') as f:
        catalog = json.load(f)
    
    datasets = {}
    for category in catalog['datasets']:
        dataset_file = test_datasets_dir / f"{category}.json"
        if dataset_file.exists():
            with open(dataset_file, 'r') as f:
                datasets[category] = json.load(f)
    
    return datasets


@pytest.fixture
def baseline_results_dir(test_baselines_dir, test_config):
    """Get baseline results directory for current config"""
    config_name = test_config.get('config_name', 'default')
    baseline_dir = test_baselines_dir / config_name / "latest"
    
    if not baseline_dir.exists():
        pytest.skip(f"No baseline results found for config {config_name}")
    
    return baseline_dir


@pytest.fixture
def test_batch_data():
    """Sample batch data for testing"""
    return {
        'batch_id': 999,
        'batch_name': 'test_batch',
        'ref_version': 'develop291',
        'proteins': [
            {'pdb_id': '1test', 'chain_id': 'A', 'sequence_length': 150},
            {'pdb_id': '1test', 'chain_id': 'B', 'sequence_length': 200},
            {'pdb_id': '2test', 'chain_id': 'A', 'sequence_length': 300}
        ]
    }


@pytest.fixture
def mock_domain_summary_factory(temp_test_dir):
    """Factory for creating mock domain summary files"""
    import xml.etree.ElementTree as ET
    
    def create_summary(pdb_id: str, chain_id: str, evidence_data: Dict[str, Any]) -> Path:
        """Create a mock domain summary XML file"""
        root = ET.Element("blast_summ_doc")
        
        # Add metadata
        blast_summ = ET.SubElement(root, "blast_summ")
        blast_summ.set("pdb", pdb_id)
        blast_summ.set("chain", chain_id)
        
        # Add chain BLAST hits
        if 'chain_blast_hits' in evidence_data:
            chain_run = ET.SubElement(root, "chain_blast_run")
            chain_run.set("program", "blastp")
            hits_elem = ET.SubElement(chain_run, "hits")
            
            for hit in evidence_data['chain_blast_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    hit_elem.set(key, str(value))
        
        # Add domain BLAST hits
        if 'domain_blast_hits' in evidence_data:
            domain_run = ET.SubElement(root, "blast_run")
            domain_run.set("program", "blastp")
            hits_elem = ET.SubElement(domain_run, "hits")
            
            for hit in evidence_data['domain_blast_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    if key not in ['query_range', 'hit_range']:
                        hit_elem.set(key, str(value))
                
                # Add range elements
                if 'query_range' in hit:
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit['query_range']
                if 'hit_range' in hit:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit['hit_range']
        
        # Add HHSearch hits
        if 'hhsearch_hits' in evidence_data:
            hh_run = ET.SubElement(root, "hh_run")
            hh_run.set("program", "hhsearch")
            hits_elem = ET.SubElement(hh_run, "hits")
            
            for hit in evidence_data['hhsearch_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    if key not in ['query_range', 'hit_range']:
                        hit_elem.set(key, str(value))
                
                # Add range elements
                if 'query_range' in hit:
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit['query_range']
                if 'hit_range' in hit:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit['hit_range']
        
        # Save to file
        summary_file = Path(temp_test_dir) / f"{pdb_id}_{chain_id}.summary.xml"
        tree = ET.ElementTree(root)
        tree.write(str(summary_file), encoding='utf-8', xml_declaration=True)
        
        return summary_file
    
    return create_summary


@pytest.fixture
def evidence_data_factory():
    """Factory for creating test evidence data"""
    
    def create_evidence(protein_info: Dict[str, Any], quality: str = 'medium') -> Dict[str, Any]:
        """Create evidence data for a protein"""
        pdb_id = protein_info['pdb_id']
        chain_id = protein_info['chain_id']
        seq_len = protein_info.get('sequence_length', 200)
        expected_domains = protein_info.get('expected_domains', 1)
        
        # Quality settings
        quality_settings = {
            'high': {'evalue_exp': 20, 'probability': 95.0, 'score': 100.0},
            'medium': {'evalue_exp': 10, 'probability': 80.0, 'score': 60.0},
            'low': {'evalue_exp': 5, 'probability': 60.0, 'score': 30.0}
        }
        
        settings = quality_settings.get(quality, quality_settings['medium'])
        
        evidence_data = {
            'sequence_length': seq_len,
            'is_peptide': seq_len < 50,
            'chain_blast_hits': [],
            'domain_blast_hits': [],
            'hhsearch_hits': []
        }
        
        if not evidence_data['is_peptide']:
            # Create domain hits
            for i in range(expected_domains):
                start = i * (seq_len // expected_domains) + 1
                end = (i + 1) * (seq_len // expected_domains)
                
                # Domain BLAST hit
                evidence_data['domain_blast_hits'].append({
                    'domain_id': f'e{pdb_id}{chain_id}{i+1}',
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'evalues': f'1e-{settings["evalue_exp"] + i*2}',
                    'hsp_count': '1',
                    'query_range': f'{start}-{end}',
                    'hit_range': f'{start}-{end}'
                })
                
                # HHSearch hit
                evidence_data['hhsearch_hits'].append({
                    'hit_id': f'h{pdb_id}{chain_id}{i+1}',
                    'domain_id': f'e{pdb_id}{chain_id}{i+1}',
                    'probability': str(settings['probability'] - i*5.0),
                    'evalue': f'1e-{settings["evalue_exp"] + i*3}',
                    'score': str(settings['score'] - i*10.0),
                    'query_range': f'{start}-{end}',
                    'hit_range': f'{start}-{end}'
                })
            
            # Chain BLAST hit
            evidence_data['chain_blast_hits'].append({
                'num': '1',
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'evalues': f'1e-{settings["evalue_exp"] + 10}',
                'hsp_count': '1',
                'query_range': f'1-{seq_len}',
                'hit_range': f'1-{seq_len}'
            })
        
        return evidence_data
    
    return create_evidence


@pytest.fixture
def performance_monitor():
    """Monitor for performance testing"""
    import time
    import psutil
    import os
    
    class PerformanceMonitor:
        def __init__(self):
            self.process = psutil.Process(os.getpid())
            self.start_time = None
            self.start_memory = None
            self.measurements = {}
        
        def start_monitoring(self, test_name: str):
            self.start_time = time.time()
            self.start_memory = self.process.memory_info().rss
            self.measurements[test_name] = {'start_time': self.start_time}
        
        def stop_monitoring(self, test_name: str):
            if test_name in self.measurements:
                end_time = time.time()
                end_memory = self.process.memory_info().rss
                
                self.measurements[test_name].update({
                    'end_time': end_time,
                    'duration': end_time - self.start_time,
                    'memory_start_mb': self.start_memory / 1024 / 1024,
                    'memory_end_mb': end_memory / 1024 / 1024,
                    'memory_delta_mb': (end_memory - self.start_memory) / 1024 / 1024
                })
        
        def get_results(self):
            return self.measurements.copy()
    
    return PerformanceMonitor()


# Command line options for integration tests
def pytest_addoption(parser):
    """Add integration test specific command line options"""
    parser.addoption(
        "--config", action="store", default="default",
        help="Test configuration to use"
    )
    parser.addoption(
        "--establish-baseline", action="store_true", default=False,
        help="Establish baseline results instead of comparing"
    )
    parser.addoption(
        "--compare-baseline", action="store_true", default=False,
        help="Compare results against baseline"
    )
    parser.addoption(
        "--baseline-dir", action="store", default=None,
        help="Directory for baseline results"
    )
    parser.addoption(
        "--results-dir", action="store", default=None,
        help="Directory for test results"
    )
    parser.addoption(
        "--skip-db-tests", action="store_true", default=False,
        help="Skip tests that require database connection"
    )


def pytest_configure(config):
    """Configure integration test markers"""
    config.addinivalue_line(
        "markers", "database: mark test as requiring database connection"
    )
    config.addinivalue_line(
        "markers", "golden: mark test as golden dataset test"
    )
    config.addinivalue_line(
        "markers", "evidence: mark test as evidence processing test"
    )
    config.addinivalue_line(
        "markers", "weights: mark test as evidence weight test"
    )
    config.addinivalue_line(
        "markers", "regression: mark test as regression detection test"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection based on command line options"""
    skip_db = pytest.mark.skip(reason="--skip-db-tests option given")
    
    for item in items:
        # Skip database tests if requested
        if config.getoption("--skip-db-tests") and "database" in item.keywords:
            item.add_marker(skip_db)


# Cleanup fixture
@pytest.fixture(autouse=True)
def cleanup_temp_files():
    """Automatically cleanup temporary files after each test"""
    yield
    # Cleanup happens automatically with temp_test_dir fixture
    pass
