#!/usr/bin/env python3
"""
Refactored Integration Test Configuration and Fixtures

MAJOR CHANGE: Removed all repository patterns and replaced with raw SQL queries
to match the actual codebase patterns used in domain_partition_run.py
"""

import pytest
import os
import sys
import json
import yaml
import tempfile
import shutil
import time
import hashlib
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from unittest.mock import Mock, patch
from contextlib import contextmanager

# Database imports
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# ECOD imports - REMOVED REPOSITORY IMPORTS
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
    config = {}
    config['host'] = os.getenv('PGHOST') or os.getenv('TEST_DB_HOST') or default_config['host']
    config['port'] = int(os.getenv('PGPORT') or os.getenv('TEST_DB_PORT') or default_config['port'])
    config['database'] = os.getenv('PGDATABASE') or os.getenv('TEST_DB_NAME') or default_config['database']
    config['user'] = os.getenv('PGUSER') or os.getenv('TEST_DB_USER') or default_config['user']
    config['password'] = os.getenv('PGPASSWORD') or os.getenv('TEST_DB_PASSWORD') or default_config['password']

    return config


@pytest.fixture(scope="session")
def validate_test_database(test_database_config_fixed):
    """Validate test database connection before running tests"""
    try:
        import psycopg2

        print(f"Testing database connection to {test_database_config_fixed['database']}...")

        conn = psycopg2.connect(**test_database_config_fixed)
        conn.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)

        with conn.cursor() as cur:
            cur.execute("SELECT version()")
            version = cur.fetchone()[0]
            print(f"  Connected successfully! PostgreSQL version: {version[:50]}...")

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

To fix this, set environment variables:
  export PGHOST=localhost
  export PGPORT=5432
  export PGDATABASE=ecod_test
  export PGUSER=test_user
  export PGPASSWORD=test_pass
"""
        pytest.skip(error_msg)


@pytest.fixture(scope="session")
def test_database_config(validate_test_database):
    """Provide validated test database configuration"""
    return validate_test_database


@pytest.fixture(scope="session")
def test_database(test_database_config):
    """Set up test database with schema and migrations using fixed config"""
    try:
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


# ===== REFACTORED TEST DATA CREATOR - RAW SQL BASED =====

class RawSQLTestDataCreator:
    """
    REFACTORED: Test data creator using raw SQL queries to match
    actual codebase patterns from domain_partition_run.py

    REMOVED: All repository patterns and imports
    ADDED: Direct SQL execution matching production code
    """

    def __init__(self, db_manager: DBManager):
        self.db = db_manager

    def create_test_batch(self, batch_dir: Path, proteins_data: List[Dict[str, Any]],
                         ref_version: str = "develop291") -> Dict[str, Any]:
        """Create test batch using raw SQL queries like the actual script"""
        batch_dir.mkdir(exist_ok=True)

        # Create batch using raw SQL INSERT - matches production patterns
        batch_insert_query = """
        INSERT INTO ecod_schema.batch (batch_name, base_path, type, ref_version, total_items, status)
        VALUES (%s, %s, %s, %s, %s, %s)
        RETURNING id
        """

        with self.db.get_connection() as conn:
            with conn.cursor() as cursor:
                cursor.execute(
                    batch_insert_query,
                    ("integration_test_batch", str(batch_dir), "domain_analysis",
                     ref_version, len(proteins_data), "processing")
                )
                batch_id = cursor.fetchone()[0]

        proteins_created = []

        for prot_data in proteins_data:
            # Create protein using raw SQL - matches get_protein_for_process() pattern
            protein_insert_query = """
            INSERT INTO ecod_schema.protein (pdb_id, chain_id, source_id, sequence_length)
            VALUES (%s, %s, %s, %s)
            RETURNING id
            """

            with self.db.get_connection() as conn:
                with conn.cursor() as cursor:
                    cursor.execute(
                        protein_insert_query,
                        (prot_data["pdb_id"], prot_data["chain_id"],
                         f"{prot_data['pdb_id']}_{prot_data['chain_id']}", prot_data["length"])
                    )
                    protein_id = cursor.fetchone()[0]

            # Create protein sequence using raw SQL if provided
            if "sequence" in prot_data:
                sequence_insert_query = """
                INSERT INTO ecod_schema.protein_sequence (protein_id, sequence, sequence_md5)
                VALUES (%s, %s, %s)
                """

                sequence = prot_data["sequence"]
                sequence_md5 = hashlib.md5(sequence.encode()).hexdigest()

                with self.db.get_connection() as conn:
                    with conn.cursor() as cursor:
                        cursor.execute(
                            sequence_insert_query,
                            (protein_id, sequence, sequence_md5)
                        )

            # Create process status using raw SQL - matches find_protein_in_database() pattern
            process_insert_query = """
            INSERT INTO ecod_schema.process_status
            (protein_id, batch_id, current_stage, status, is_representative)
            VALUES (%s, %s, %s, %s, %s)
            RETURNING id
            """

            with self.db.get_connection() as conn:
                with conn.cursor() as cursor:
                    cursor.execute(
                        process_insert_query,
                        (protein_id, batch_id, "domain_summary_complete", "ready",
                         prot_data.get("is_rep", False))
                    )
                    process_id = cursor.fetchone()[0]

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
        """Create test file records using raw SQL - matches script's file checking"""
        file_insert_query = """
        INSERT INTO ecod_schema.process_file (process_id, file_type, file_path, file_exists, file_size)
        VALUES (%s, %s, %s, %s, %s)
        """

        with self.db.get_connection() as conn:
            with conn.cursor() as cursor:
                for file_type, file_path in file_mappings.items():
                    # Check if file actually exists
                    full_path = Path(file_path)
                    file_exists = full_path.exists()
                    file_size = full_path.stat().st_size if file_exists else 0

                    cursor.execute(
                        file_insert_query,
                        (process_id, file_type, file_path, file_exists, file_size)
                    )

    def create_domain_summary_file(self, pdb_id: str, chain_id: str,
                                  batch_dir: Path, ref_version: str = "develop291",
                                  evidence_data: Optional[Dict] = None) -> Path:
        """Create mock domain summary XML file"""
        import xml.etree.ElementTree as ET

        # Create XML structure
        root = ET.Element("blast_summ_doc")

        # Add metadata
        blast_summ = ET.SubElement(root, "blast_summ")
        blast_summ.set("pdb", pdb_id)
        blast_summ.set("chain", chain_id)

        # Add basic evidence if provided
        if evidence_data:
            # Add chain BLAST hits
            if 'chain_blast_hits' in evidence_data:
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
            if 'domain_blast_hits' in evidence_data:
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
            if 'hhsearch_hits' in evidence_data:
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

        # Save to file with proper naming convention
        domains_dir = batch_dir / "domains"
        domains_dir.mkdir(exist_ok=True)

        summary_file = domains_dir / f"{pdb_id}_{chain_id}.{ref_version}.domains_v14.xml"
        tree = ET.ElementTree(root)
        tree.write(str(summary_file), encoding='utf-8', xml_declaration=True)

        return summary_file

    def get_batch_info_like_script(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information using the same query as the actual script"""
        # EXACT SAME QUERY as get_batch_info() in domain_partition_run.py
        query = """
        SELECT id, batch_name, base_path, ref_version, status
        FROM ecod_schema.batch
        WHERE id = %s
        """

        with self.db.get_connection() as conn:
            with conn.cursor() as cursor:
                cursor.execute(query, (batch_id,))
                row = cursor.fetchone()

                if row:
                    return {
                        'id': row[0],
                        'batch_name': row[1],
                        'base_path': row[2],
                        'ref_version': row[3],
                        'status': row[4]
                    }
        return None

    def find_protein_like_script(self, pdb_id: str, chain_id: str,
                                batch_id: Optional[int] = None) -> List[Dict[str, Any]]:
        """Find protein using the same query as the actual script"""
        # EXACT SAME QUERY as find_protein_in_database() in domain_partition_run.py
        query = """
        SELECT ps.id as process_id, ps.batch_id, ps.current_stage, ps.status,
               ps.error_message, ps.is_representative
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        """

        params = [pdb_id, chain_id]

        if batch_id is not None:
            query += " AND ps.batch_id = %s"
            params.append(batch_id)

        query += " ORDER BY ps.batch_id DESC"

        with self.db.get_connection() as conn:
            with conn.cursor() as cursor:
                cursor.execute(query, tuple(params))
                rows = cursor.fetchall()

                results = []
                for row in rows:
                    result = {
                        'process_id': row[0],
                        'batch_id': row[1],
                        'current_stage': row[2],
                        'status': row[3],
                        'error_message': row[4],
                        'is_representative': row[5]
                    }

                    # Add file information like the script does
                    file_query = """
                    SELECT file_type, file_path, file_exists, file_size
                    FROM ecod_schema.process_file
                    WHERE process_id = %s
                    """

                    cursor.execute(file_query, (row[0],))
                    file_rows = cursor.fetchall()

                    files = {}
                    for file_row in file_rows:
                        files[file_row[0]] = {
                            'file_type': file_row[0],
                            'file_path': file_row[1],
                            'file_exists': file_row[2],
                            'file_size': file_row[3]
                        }

                    result['files'] = files
                    results.append(result)

        return results


@pytest.fixture
def raw_sql_test_data_creator(clean_database):
    """Provide refactored test data creator using raw SQL"""
    return RawSQLTestDataCreator(clean_database)


# ===== GOLDEN DATASETS =====

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


# ===== SCRIPT EXECUTION UTILITIES =====

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


# ===== EVIDENCE DATA FACTORIES =====

@pytest.fixture
def evidence_data_factory():
    """Factory for creating test evidence data"""

    def create_evidence(protein_info: Dict[str, Any], quality: str = 'medium') -> Dict[str, Any]:
        """Create evidence data for a protein"""
        pdb_id = protein_info['pdb_id']
        chain_id = protein_info['chain_id']
        seq_len = protein_info.get('sequence_length', protein_info.get('length', 200))
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
            'hhsearch_hits': [],
            'expected_domains': expected_domains
        }

        if not evidence_data['is_peptide'] and expected_domains > 0:
            # Create domain hits based on expected domain count
            for i in range(expected_domains):
                if expected_domains == 1:
                    start, end = 1, seq_len
                else:
                    domain_size = seq_len // expected_domains
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
                'hsp_count': '1',
                'query_range': f'1-{seq_len}',
                'hit_range': f'1-{seq_len}'
            })

        return evidence_data

    return create_evidence


# ===== SAMPLE PROTEIN DATA =====

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


# ===== CLEANUP =====

@pytest.fixture(autouse=True)
def cleanup_temp_files():
    """Automatically cleanup temporary files after each test"""
    yield
    # Cleanup happens automatically with temp_test_dir fixture
    pass
