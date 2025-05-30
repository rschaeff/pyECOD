#!/usr/bin/env python3
"""
Integration Test Configuration and Fixtures

Synthesized configuration combining robust database setup, service architecture support,
golden datasets, baseline testing, and performance monitoring for comprehensive
integration testing of the ECOD domain partition pipeline.
"""

import pytest
import os
import sys
import json
import yaml
import tempfile
import shutil
import time
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from unittest.mock import Mock, patch
from contextlib import contextmanager

# Database imports
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# ECOD imports
from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.db.migration_manager import MigrationManager


class IntegrationTestConfig:
    """Enhanced configuration class for integration tests"""
    
    def __init__(self):
        # Database configuration with environment variable support
        self.test_db_config = {
            'host': os.getenv('PGHOST', os.getenv('TEST_DB_HOST', 'localhost')),
            'port': int(os.getenv('PGPORT', os.getenv('TEST_DB_PORT', '5432'))),
            'database': os.getenv('PGDATABASE', os.getenv('TEST_DB_NAME', 'ecod_test')),
            'user': os.getenv('PGUSER', os.getenv('TEST_DB_USER', 'test_user')),
            'password': os.getenv('PGPASSWORD', os.getenv('TEST_DB_PASSWORD', 'test_pass'))
        }
        
        # Script paths
        self.script_paths = {
            'domain_partition': PROJECT_ROOT / "scripts" / "domain_partition_run.py"
        }
        
        # Verify critical scripts exist
        for script_name, script_path in self.script_paths.items():
            if not script_path.exists():
                print(f"Warning: Script not found: {script_path}")


# ===== SESSION-SCOPED FIXTURES =====

@pytest.fixture(scope="session")
def integration_config():
    """Integration test configuration"""
    return IntegrationTestConfig()


@pytest.fixture(scope="session")
def integration_test_root():
    """Get the integration test root directory"""
    return Path(__file__).parent


@pytest.fixture(scope="session") 
def test_configs_dir(integration_test_root):
    """Get the test configurations directory"""
    configs_dir = integration_test_root / "configs"
    configs_dir.mkdir(exist_ok=True)
    return configs_dir


@pytest.fixture(scope="session")
def test_datasets_dir(integration_test_root):
    """Get the test datasets directory"""
    datasets_dir = integration_test_root / "datasets"
    datasets_dir.mkdir(exist_ok=True)
    return datasets_dir


@pytest.fixture(scope="session")
def test_baselines_dir(integration_test_root):
    """Get the test baselines directory"""
    baselines_dir = integration_test_root / "baselines"
    baselines_dir.mkdir(exist_ok=True)
    return baselines_dir


@pytest.fixture(scope="session")
def test_results_dir(integration_test_root):
    """Get the test results directory"""
    results_dir = integration_test_root / "results"
    results_dir.mkdir(exist_ok=True)
    return results_dir


@pytest.fixture(scope="session")
def test_config(test_configs_dir, request):
    """Load test configuration with fallback to default"""
    # Get config name from command line or use default
    config_name = getattr(request.config.option, 'config', 'default')
    config_file = test_configs_dir / f"{config_name}.yml"
    
    if not config_file.exists():
        config_file = test_configs_dir / "default.yml"
        
        # Create default config if it doesn't exist
        if not config_file.exists():
            default_config = create_default_test_config()
            with open(config_file, 'w') as f:
                yaml.dump(default_config, f, default_flow_style=False)
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Add config name for reference
    config['config_name'] = config_name
    return config


def create_default_test_config():
    """Create default test configuration"""
    return {
        'database': {
            'host': 'localhost',
            'port': 5432,
            'database': 'ecod_test',
            'user': 'test_user',
            'password': 'test_pass'
        },
        'reference': {
            'current_version': 'develop291'
        },
        'partition': {
            'confidence_thresholds': {'high': 0.9, 'medium': 0.7, 'low': 0.5},
            'evidence_weights': {'domain_blast': 3.0, 'hhsearch': 2.5, 'chain_blast': 2.0},
            'overlap_tolerance': 0.15,
            'min_domain_size': 20,
            'peptide_threshold': 50
        },
        'job_manager': {'type': 'local'},
        'services': {
            'domain_partition': {
                'max_workers': 2,
                'use_multiprocessing': False,
                'batch_size': 10,
                'save_intermediate': True,
                'track_status': True
            }
        },
        'tools': {
            'blast_path': '/usr/bin/blastp',
            'hhblits_path': '/usr/bin/hhblits',
            'hhsearch_path': '/usr/bin/hhsearch'
        },
        'logging': {
            'level': 'DEBUG',
            'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        }
    }


@pytest.fixture(scope="session")
def test_database_config(integration_config, test_config):
    """Get database configuration for tests with environment override"""
    # Start with config file
    db_config = test_config.get('database', {}).copy()
    
    # Override with integration config (which includes env vars)
    db_config.update(integration_config.test_db_config)
    
    return db_config


@pytest.fixture(scope="session")
def test_database_connection(test_database_config):
    """Create and validate test database connection"""
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


@pytest.fixture(scope="session")
def test_database(test_database_config):
    """Set up test database with schema and migrations"""
    try:
        db_manager = DBManager(test_database_config)
        
        # Apply migrations to ensure schema exists
        migration_manager = MigrationManager(
            test_database_config, 
            str(PROJECT_ROOT / "ecod" / "db" / "migrations")
        )
        try:
            migration_manager.apply_migrations()
        except Exception as e:
            print(f"Migration warning: {e}")
        
        yield db_manager
        
    except Exception as e:
        pytest.skip(f"Could not set up test database: {e}")


@pytest.fixture(scope="session")
def script_paths(integration_config):
    """Provide paths to scripts"""
    return integration_config.script_paths


# ===== FUNCTION-SCOPED FIXTURES =====

@pytest.fixture
def clean_database(test_database):
    """Provide clean database for each test"""
    # Clean up before test
    with test_database.get_connection() as conn:
        with conn.cursor() as cursor:
            # Support both schema types
            tables = [
                # pdb_analysis schema
                'pdb_analysis.domains',
                'pdb_analysis.domain_sequence', 
                'pdb_analysis.domain_dssp_detail',
                'pdb_analysis.chain_sequences',
                'pdb_analysis.pdb_chains',
                'pdb_analysis.protein_structure',
                'pdb_analysis.protein_sequence',
                'pdb_analysis.protein',
                'pdb_analysis.pdb_entries',
                # ecod_schema
                'ecod_schema.job_item',
                'ecod_schema.job', 
                'ecod_schema.process_file',
                'ecod_schema.process_status',
                'ecod_schema.protein_sequence',
                'ecod_schema.protein',
                'ecod_schema.batch',
            ]
            for table in tables:
                try:
                    cursor.execute(f"TRUNCATE TABLE {table} CASCADE")
                except Exception:
                    # Table might not exist, continue
                    pass
    
    yield test_database
    
    # Clean up after test
    with test_database.get_connection() as conn:
        with conn.cursor() as cursor:
            for table in tables:
                try:
                    cursor.execute(f"TRUNCATE TABLE {table} CASCADE")
                except Exception:
                    pass


@pytest.fixture
def temp_test_dir():
    """Create temporary directory for each test"""
    with tempfile.TemporaryDirectory(prefix="ecod_test_") as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create standard subdirectories
        (temp_path / "domains").mkdir()
        (temp_path / "sequences").mkdir()
        (temp_path / "logs").mkdir()
        (temp_path / "temp").mkdir()
        
        yield temp_path


@pytest.fixture
def test_config_file(test_config, temp_test_dir, test_database_config):
    """Create comprehensive test configuration file"""
    config = test_config.copy()
    
    # Update database config
    config['database'] = test_database_config
    
    # Update paths
    config['paths'] = {
        'output_dir': str(temp_test_dir / "output"),
        'temp_dir': str(temp_test_dir / "temp"),
        'data_dir': str(temp_test_dir / "data")
    }
    
    # Create necessary directories
    for path_key, path_value in config['paths'].items():
        os.makedirs(path_value, exist_ok=True)
    
    config_file = temp_test_dir / "test_config.yml"
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    return str(config_file)


@pytest.fixture
def mock_application_context(test_config, temp_test_dir):
    """Create mock application context for integration tests"""
    context = Mock(spec=ApplicationContext)
    
    # Mock config manager
    context.config_manager = Mock()
    context.config_manager.config = test_config
    context.config_manager.get_db_config.return_value = test_config['database']
    context.config_manager.get_path.return_value = str(temp_test_dir)
    
    # Mock database manager
    context.db_manager = Mock()
    
    # Mock logger
    context.logger = Mock()
    
    return context


# ===== GOLDEN DATASETS AND BASELINES =====

@pytest.fixture
def golden_datasets(test_datasets_dir):
    """Load golden datasets"""
    catalog_file = test_datasets_dir / "catalog.json"
    
    if not catalog_file.exists():
        # Create default catalog if none exists
        default_catalog = {
            'datasets': ['single_domain', 'multi_domain', 'peptides', 'complex'],
            'version': '1.0',
            'description': 'Test datasets for domain partition validation'
        }
        with open(catalog_file, 'w') as f:
            json.dump(default_catalog, f, indent=2)
        
        pytest.skip("Golden datasets not available - created default catalog")
    
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
def baseline_results_dir(test_baselines_dir, test_config, request):
    """Get baseline results directory for current config"""
    config_name = test_config.get('config_name', 'default')
    baseline_dir = test_baselines_dir / config_name / "latest"
    
    # Check if we should establish baseline
    if getattr(request.config.option, 'establish_baseline', False):
        baseline_dir.mkdir(parents=True, exist_ok=True)
        return baseline_dir
    
    if not baseline_dir.exists():
        pytest.skip(f"No baseline results found for config {config_name}")
    
    return baseline_dir


# ===== TEST DATA FACTORIES =====

@pytest.fixture
def sample_protein_data():
    """Sample protein data for testing"""
    return [
        {
            "pdb_id": "1TST", 
            "chain_id": "A", 
            "length": 150, 
            "is_rep": True,
            "sequence": "M" + "AKVLTKSPG" * 16 + "AKVLT",  # 150 chars
            "expected_domains": 1
        },
        {
            "pdb_id": "2TST", 
            "chain_id": "A", 
            "length": 300, 
            "is_rep": False,
            "sequence": "M" + "AKVLTKSPG" * 33 + "AKVLTK",  # 300 chars
            "expected_domains": 2
        }, 
        {
            "pdb_id": "3TST", 
            "chain_id": "B", 
            "length": 45, 
            "is_rep": True,  # Peptide
            "sequence": "M" + "AKVLTKSPG" * 4 + "AKVL",  # 45 chars
            "expected_domains": 0
        }
    ]


@pytest.fixture
def test_batch_data():
    """Enhanced sample batch data for testing"""
    return {
        'batch_id': 999,
        'batch_name': 'integration_test_batch',
        'ref_version': 'develop291',
        'proteins': [
            {'pdb_id': '1test', 'chain_id': 'A', 'sequence_length': 150, 'expected_domains': 1},
            {'pdb_id': '1test', 'chain_id': 'B', 'sequence_length': 200, 'expected_domains': 1},
            {'pdb_id': '2test', 'chain_id': 'A', 'sequence_length': 300, 'expected_domains': 2},
            {'pdb_id': '3test', 'chain_id': 'A', 'sequence_length': 45, 'expected_domains': 0}  # Peptide
        ]
    }


@pytest.fixture
def evidence_data_factory():
    """Factory for creating test evidence data matching script parsing expectations"""

    def create_evidence(protein_info: Dict[str, Any], quality: str = 'medium') -> Dict[str, Any]:
        """Create evidence data for a protein"""
        pdb_id = protein_info['pdb_id']
        chain_id = protein_info['chain_id']
        seq_len = protein_info.get('sequence_length', 200)
        expected_domains = protein_info.get('expected_domains', 1)

        # Quality settings affect e-values, probabilities, and scores
        quality_settings = {
            'high': {'evalue_exp': 25, 'probability': 98.0, 'score': 120.0},
            'medium': {'evalue_exp': 15, 'probability': 85.0, 'score': 75.0},
            'low': {'evalue_exp': 8, 'probability': 65.0, 'score': 35.0}
        }

        settings = quality_settings.get(quality, quality_settings['medium'])

        evidence_data = {
            'sequence_length': seq_len,
            'is_peptide': seq_len < 50,
            'chain_blast_hits': [],
            'domain_blast_hits': [],
            'hhsearch_hits': []
        }

        if not evidence_data['is_peptide'] and expected_domains > 0:
            # Create domain hits based on expected domain count
            domain_size = seq_len // expected_domains if expected_domains > 0 else seq_len

            for i in range(expected_domains):
                start = i * domain_size + 1
                end = min((i + 1) * domain_size, seq_len)

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
                    'probability': str(max(50.0, settings['probability'] - i*8.0)),
                    'evalue': f'1e-{settings["evalue_exp"] + i*3}',
                    'score': str(max(20.0, settings['score'] - i*15.0)),
                    'query_range': f'{start}-{end}',
                    'hit_range': f'{start}-{end}'
                })

            # Chain BLAST hit (covers full sequence)
            evidence_data['chain_blast_hits'].append({
                'num': '1',
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'evalues': f'1e-{settings["evalue_exp"] + 5}',
                'hsp_count': str(expected_domains),
                'query_range': f'1-{seq_len}',
                'hit_range': f'1-{seq_len}'
            })
        
        return evidence_data
    
    return create_evidence


@pytest.fixture
def mock_domain_summary_factory(temp_test_dir):
    """Factory for creating mock domain summary files that match real format"""
    import xml.etree.ElementTree as ET

    def create_summary(pdb_id: str, chain_id: str, evidence_data: Dict[str, Any],
                      reference: str = "develop291") -> Path:
        """Create a mock domain summary XML file matching expected format"""
        root = ET.Element("blast_summ_doc")

        # Add metadata
        blast_summ = ET.SubElement(root, "blast_summ")
        blast_summ.set("pdb", pdb_id)
        blast_summ.set("chain", chain_id)

        # Add chain BLAST hits
        if 'chain_blast_hits' in evidence_data and evidence_data['chain_blast_hits']:
            chain_run = ET.SubElement(root, "chain_blast_run")
            chain_run.set("program", "blastp")
            hits_elem = ET.SubElement(chain_run, "hits")

            for hit in evidence_data['chain_blast_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    if key not in ['query_range', 'hit_range']:
                        hit_elem.set(key, str(value))

                if 'query_range' in hit:
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit['query_range']
                if 'hit_range' in hit:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit['hit_range']

        # Add domain BLAST hits
        if 'domain_blast_hits' in evidence_data and evidence_data['domain_blast_hits']:
            domain_run = ET.SubElement(root, "blast_run")
            domain_run.set("program", "blastp")
            hits_elem = ET.SubElement(domain_run, "hits")

            for hit in evidence_data['domain_blast_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    if key not in ['query_range', 'hit_range']:
                        hit_elem.set(key, str(value))

                if 'query_range' in hit:
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit['query_range']
                if 'hit_range' in hit:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit['hit_range']

        # Add HHSearch hits
        if 'hhsearch_hits' in evidence_data and evidence_data['hhsearch_hits']:
            hh_run = ET.SubElement(root, "hh_run")
            hh_run.set("program", "hhsearch")
            hits_elem = ET.SubElement(hh_run, "hits")

            for hit in evidence_data['hhsearch_hits']:
                hit_elem = ET.SubElement(hits_elem, "hit")
                for key, value in hit.items():
                    if key not in ['query_range', 'hit_range']:
                        hit_elem.set(key, str(value))

                if 'query_range' in hit:
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = hit['query_range']
                if 'hit_range' in hit:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit['hit_range']

        # Save to file with proper naming convention expected by scripts
        domains_dir = Path(temp_test_dir) / "domains"
        domains_dir.mkdir(exist_ok=True)

        summary_file = domains_dir / f"{pdb_id}_{chain_id}.{reference}.domains_v14.xml"
        tree = ET.ElementTree(root)
        tree.write(str(summary_file), encoding='utf-8', xml_declaration=True)

        return summary_file

    return create_summary



# ===== PERFORMANCE MONITORING =====

@pytest.fixture
def performance_monitor():
    """Enhanced monitor for performance testing"""
    import psutil
    
    class PerformanceMonitor:
        def __init__(self):
            self.process = psutil.Process(os.getpid())
            self.start_time = None
            self.start_memory = None
            self.start_cpu = None
            self.measurements = {}
        
        def start_monitoring(self, test_name: str):
            self.start_time = time.time()
            self.start_memory = self.process.memory_info().rss
            self.start_cpu = self.process.cpu_percent()
            self.measurements[test_name] = {
                'start_time': self.start_time,
                'start_memory_mb': self.start_memory / 1024 / 1024
            }
        
        def stop_monitoring(self, test_name: str):
            if test_name in self.measurements:
                end_time = time.time()
                end_memory = self.process.memory_info().rss
                end_cpu = self.process.cpu_percent()
                
                self.measurements[test_name].update({
                    'end_time': end_time,
                    'duration': end_time - self.start_time,
                    'memory_end_mb': end_memory / 1024 / 1024,
                    'memory_delta_mb': (end_memory - self.start_memory) / 1024 / 1024,
                    'cpu_percent': end_cpu,
                    'proteins_per_second': 0  # Will be calculated by caller
                })
        
        def get_results(self) -> Dict[str, Any]:
            return self.measurements.copy()
        
        def save_results(self, output_file: Path):
            """Save performance results to file"""
            with open(output_file, 'w') as f:
                json.dump(self.get_results(), f, indent=2)
    
    return PerformanceMonitor()


# ===== HELPER CLASSES AND UTILITIES =====

class TestDataCreator:
    """Test data creator using raw SQL queries to match actual codebase patterns"""

    def __init__(self, db_manager):
        self.db = db_manager

    def create_test_batch(self, batch_dir: Path, proteins_data: List[Dict[str, Any]],
                         ref_version: str = "develop291") -> Dict[str, Any]:
        """Create test batch using raw SQL queries like domain_partition_run.py"""
        batch_dir.mkdir(exist_ok=True)

        # Create batch using raw SQL INSERT - matches get_batch_info() expectations
        batch_insert_query = """
        INSERT INTO ecod_schema.batch (batch_name, base_path, ref_version, status, total_items, type, created_at)
        VALUES (%s, %s, %s, %s, %s, %s, CURRENT_TIMESTAMP)
        RETURNING id
        """

        batch_result = self.db.execute_query(
            batch_insert_query,
            ("integration_test_batch", str(batch_dir), ref_version,
             "processing", len(proteins_data), "domain_analysis")
        )
        batch_id = batch_result[0][0]

        proteins_created = []

        for prot_data in proteins_data:
            # Create protein using raw SQL - matches find_protein_in_database() expectations
            protein_insert_query = """
            INSERT INTO ecod_schema.protein (pdb_id, chain_id, sequence_length, created_at)
            VALUES (%s, %s, %s, CURRENT_TIMESTAMP)
            RETURNING id
            """

            protein_result = self.db.execute_query(
                protein_insert_query,
                (prot_data["pdb_id"], prot_data["chain_id"], prot_data["length"])
            )
            protein_id = protein_result[0][0]

            # Create protein sequence if provided (skip if table doesn't exist)
            if "sequence" in prot_data:
                try:
                    sequence_insert_query = """
                    INSERT INTO ecod_schema.protein_sequence (protein_id, sequence, sequence_md5, created_at)
                    VALUES (%s, %s, %s, CURRENT_TIMESTAMP)
                    """

                    sequence = prot_data["sequence"]
                    import hashlib
                    sequence_md5 = hashlib.md5(sequence.encode()).hexdigest()

                    self.db.execute_query(
                        sequence_insert_query,
                        (protein_id, sequence, sequence_md5)
                    )
                except Exception:
                    # Table might not exist in minimal test schema, skip
                    pass

            # Create process status using raw SQL - matches script expectations
            process_insert_query = """
            INSERT INTO ecod_schema.process_status
            (protein_id, batch_id, current_stage, status, is_representative, created_at, updated_at)
            VALUES (%s, %s, %s, %s, %s, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)
            RETURNING id
            """

            process_result = self.db.execute_query(
                process_insert_query,
                (protein_id, batch_id, "domain_summary_complete", "ready",
                 prot_data.get("is_rep", False))
            )
            process_id = process_result[0][0]

            proteins_created.append({
                "protein_id": protein_id,
                "process_id": process_id,
                **prot_data
            })

        return {
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
            "proteins": proteins_created
        }

    def create_test_files(self, process_id: int, file_mappings: Dict[str, str]):
        """Create test file records using raw SQL - matches script file checking"""
        file_insert_query = """
        INSERT INTO ecod_schema.process_file (process_id, file_type, file_path, file_exists, file_size, created_at)
        VALUES (%s, %s, %s, %s, %s, CURRENT_TIMESTAMP)
        """

        for file_type, file_path in file_mappings.items():
            # Check if file actually exists
            full_path = Path(file_path)
            file_exists = full_path.exists()
            file_size = full_path.stat().st_size if file_exists else 0

            self.db.execute_query(
                file_insert_query,
                (process_id, file_type, file_path, file_exists, file_size)
            )


@pytest.fixture
def test_data_creator(clean_database):
    """Provide enhanced test data creator"""
    return TestDataCreator(clean_database)


# ===== UTILITY FUNCTIONS =====

def run_script_subprocess(script_path: Union[str, Path], config_file: str, args: List[str], 
                         timeout: int = 60, cwd: Optional[str] = None):
    """Helper function to run scripts as subprocess"""
    import subprocess
    
    cmd = [
        sys.executable, str(script_path),
        "--config", config_file,
        "--verbose"
    ] + args
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=cwd or str(PROJECT_ROOT)
    )
    
    return result


@contextmanager
def temporary_environment(**env_vars):
    """Context manager to temporarily set environment variables"""
    old_environ = dict(os.environ)
    os.environ.update(env_vars)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(old_environ)


# ===== PYTEST CONFIGURATION =====

#def pytest_addoption(parser):
    """Add integration test specific command line options"""
#    parser.addoption(
#        "--config", action="store", default="default",
#        help="Test configuration to use"
#    )
#    parser.addoption(
#        "--establish-baseline", action="store_true", default=False,
#        help="Establish baseline results instead of comparing"
#    )
#    parser.addoption(
#        "--compare-baseline", action="store_true", default=False,
#        help="Compare results against baseline"
#    )
#    parser.addoption(
#        "--baseline-dir", action="store", default=None,
#        help="Directory for baseline results"
#    )
#    parser.addoption(
#        "--results-dir", action="store", default=None,
#        help="Directory for test results"
#    )
#    parser.addoption(
#        "--skip-db-tests", action="store_true", default=False,
#       help="Skip tests that require database connection"
#    )
#    parser.addoption(
#        "--performance-tests", action="store_true", default=False,
#        help="Enable performance testing mode"
#    )


def pytest_configure(config):
    """Configure integration test markers"""
    # Integration-specific markers
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
    config.addinivalue_line(
        "markers", "baseline: mark test as baseline comparison test"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection based on command line options"""
    skip_db = pytest.mark.skip(reason="--skip-db-tests option given")
    
    for item in items:
        # Skip database tests if requested
        if config.getoption("--skip-db-tests") and "database" in item.keywords:
            item.add_marker(skip_db)
        
        # Auto-mark tests based on location and content
        if "integration" in str(item.fspath):
            item.add_marker(pytest.mark.integration)
        
        if "database" in item.name or "test_db" in item.name:
            item.add_marker(pytest.mark.database)

# Add this to your integration/conftest.py to fix database connection issues

# ===== ENHANCED DATABASE CONFIGURATION =====

@pytest.fixture(scope="session")
def test_database_config_fixed():
    """Get database configuration with proper environment variable handling"""

    # Default test database configuration
    default_config = {
        'host': 'localhost',
        'port': 5432,
        'database': 'ecod_test',
        'user': 'test_user',
        'password': 'test_pass'
    }

    # Override with environment variables if they exist
    # This handles both PG* and TEST_DB_* prefixes
    config = {}

    # Check for PostgreSQL standard environment variables first
    config['host'] = os.getenv('PGHOST') or os.getenv('TEST_DB_HOST') or default_config['host']
    config['port'] = int(os.getenv('PGPORT') or os.getenv('TEST_DB_PORT') or default_config['port'])
    config['database'] = os.getenv('PGDATABASE') or os.getenv('TEST_DB_NAME') or default_config['database']
    config['user'] = os.getenv('PGUSER') or os.getenv('TEST_DB_USER') or default_config['user']
    config['password'] = os.getenv('PGPASSWORD') or os.getenv('TEST_DB_PASSWORD') or default_config['password']

    # Debug output to help troubleshoot
    print(f"\nTest Database Configuration:")
    print(f"  Host: {config['host']}")
    print(f"  Port: {config['port']}")
    print(f"  Database: {config['database']}")
    print(f"  User: {config['user']}")
    print(f"  Password: {'*' * len(str(config['password']))}")

    return config


@pytest.fixture(scope="session")
def validate_test_database(test_database_config_fixed):
    """Validate test database connection before running tests"""
    try:
        import psycopg2

        # Test the connection
        print(f"Testing database connection to {test_database_config_fixed['database']}...")

        conn = psycopg2.connect(**test_database_config_fixed)
        conn.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)

        # Test basic query
        with conn.cursor() as cur:
            cur.execute("SELECT version()")
            version = cur.fetchone()[0]
            print(f"  Connected successfully! PostgreSQL version: {version[:50]}...")

            # Test if we can create/drop a test table
            try:
                cur.execute("CREATE TEMP TABLE connection_test (id INTEGER)")
                cur.execute("DROP TABLE connection_test")
                print("  Database permissions: OK")
            except Exception as e:
                print(f"  Warning: Limited database permissions: {e}")

        conn.close()
        return test_database_config_fixed

    except psycopg2.Error as e:
        error_msg = f"""
Database connection failed!

Configuration used:
  Host: {test_database_config_fixed['host']}
  Port: {test_database_config_fixed['port']}
  Database: {test_database_config_fixed['database']}
  User: {test_database_config_fixed['user']}

Error: {e}

To fix this, either:
1. Set environment variables:
   export TEST_DB_HOST=localhost
   export TEST_DB_PORT=5432
   export TEST_DB_NAME=ecod_test
   export TEST_DB_USER=test_user
   export TEST_DB_PASSWORD=test_pass

2. Or use PostgreSQL standard variables:
   export PGHOST=localhost
   export PGPORT=5432
   export PGDATABASE=ecod_test
   export PGUSER=test_user
   export PGPASSWORD=test_pass

3. Ensure the database exists:
   psql -U test_user -h localhost -c "CREATE DATABASE ecod_test;" postgres

4. Test the connection manually:
   psql -U test_user -d ecod_test -h localhost
"""
        pytest.skip(error_msg)


# Override the existing fixtures to use the fixed version
@pytest.fixture(scope="session")
def test_database_config(validate_test_database):
    """Provide validated test database configuration"""
    return validate_test_database


@pytest.fixture(scope="session")
def test_database(test_database_config):
    """Set up test database with schema and migrations using fixed config"""
    try:
        from ecod.db import DBManager
        from ecod.db.migration_manager import MigrationManager

        print(f"Setting up test database with DBManager...")
        db_manager = DBManager(test_database_config)

        # Test the DBManager connection
        with db_manager.get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute("SELECT current_database(), current_user")
                db_name, user = cur.fetchone()
                print(f"  DBManager connected to database '{db_name}' as user '{user}'")

        # Try to apply migrations
        try:
            migration_dir = str(PROJECT_ROOT / "ecod" / "db" / "migrations")
            if Path(migration_dir).exists():
                print(f"Applying migrations from {migration_dir}...")
                migration_manager = MigrationManager(test_database_config, migration_dir)
                migration_manager.apply_migrations()
                print("  Migrations applied successfully")
            else:
                print(f"  No migrations directory found at {migration_dir}")
        except Exception as e:
            print(f"  Migration warning (continuing anyway): {e}")

        yield db_manager

    except Exception as e:
        pytest.skip(f"Could not set up test database with DBManager: {e}")


# ===== DATABASE CONNECTION HELPERS =====

def test_database_connection_manual(config):
    """Test database connection manually with detailed error reporting"""
    try:
        import psycopg2

        # Try to connect
        conn = psycopg2.connect(
            host=config['host'],
            port=config['port'],
            database=config['database'],
            user=config['user'],
            password=config['password']
        )

        # Test basic operations
        with conn.cursor() as cur:
            cur.execute("SELECT 1")
            result = cur.fetchone()

        conn.close()
        return True, "Connection successful"

    except psycopg2.OperationalError as e:
        if "does not exist" in str(e):
            return False, f"Database '{config['database']}' does not exist. Create it with: createdb -U {config['user']} -h {config['host']} {config['database']}"
        elif "authentication failed" in str(e):
            return False, f"Authentication failed for user '{config['user']}'. Check username/password."
        elif "could not connect" in str(e):
            return False, f"Could not connect to server at {config['host']}:{config['port']}. Is PostgreSQL running?"
        else:
            return False, f"Connection error: {e}"
    except Exception as e:
        return False, f"Unexpected error: {e}"


@pytest.fixture
def db_connection_test():
    """Provide database connection testing utility"""
    return test_database_connection_manual


# ===== CONFIGURATION FILE OVERRIDE =====

def create_test_config_with_db(db_config: Dict[str, Any], temp_dir: Path) -> Dict[str, Any]:
    """Create test configuration with proper database settings"""
    config = {
        'database': db_config,  # Use the validated database config
        'reference': {'current_version': 'develop291'},
        'partition': {
            'confidence_thresholds': {'high': 0.9, 'medium': 0.7, 'low': 0.5},
            'evidence_weights': {'domain_blast': 3.0, 'hhsearch': 2.5, 'chain_blast': 2.0},
            'overlap_tolerance': 0.15,
            'min_domain_size': 20,
            'peptide_threshold': 50
        },
        'job_manager': {'type': 'local'},
        'paths': {
            'output_dir': str(temp_dir / "output"),
            'temp_dir': str(temp_dir / "temp"),
            'data_dir': str(temp_dir / "data")
        },
        'services': {
            'domain_partition': {
                'max_workers': 1,  # Keep low for tests
                'use_multiprocessing': False,
                'batch_size': 5,
                'save_intermediate': True,
                'track_status': True
            }
        },
        'tools': {
            'blast_path': '/usr/bin/blastp',
            'hhblits_path': '/usr/bin/hhblits',
            'hhsearch_path': '/usr/bin/hhsearch'
        },
        'logging': {
            'level': 'INFO',
            'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        }
    }

    return config


@pytest.fixture
def test_config_file(test_database_config, temp_test_dir):
    """Create test configuration file with validated database config"""
    config = create_test_config_with_db(test_database_config, temp_test_dir)

    # Create necessary directories
    for path_key, path_value in config['paths'].items():
        os.makedirs(path_value, exist_ok=True)

    config_file = temp_test_dir / "test_config.yml"
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    print(f"Created test config file: {config_file}")
    print(f"Database config in file: {config['database']}")

    return str(config_file)

# ===== CLEANUP =====

@pytest.fixture(autouse=True)
def cleanup_temp_files():
    """Automatically cleanup temporary files after each test"""
    yield
    # Cleanup happens automatically with temp_test_dir fixture
    pass
