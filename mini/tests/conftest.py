#!/usr/bin/env python3
"""
Enhanced test configuration and fixtures for mini_pyecod tests

This file provides improved test fixtures with better error handling,
environment detection, and flexible test data loading.
"""

import os
import sys
import pytest
from pathlib import Path
from typing import Dict, List, Optional, Any

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.core.parser import load_reference_lengths, load_protein_lengths
from mini.core.decomposer import load_domain_definitions
from mini.core.blast_parser import load_chain_blast_alignments
from mini.core.models import Evidence
from ecod.core.sequence_range import SequenceRange


class TestEnvironment:
    """Test environment detection and validation"""

    def __init__(self):
        self.base_dir = Path("/data/ecod/pdb_updates/batches")
        self.stable_batch = "ecod_batch_036_20250406_1424"
        self._available_batches = None
        self._test_data_dir = Path(__file__).parent.parent / "test_data"

    @property
    def available_batches(self) -> List[str]:
        """Get list of available batch directories"""
        if self._available_batches is None:
            if not self.base_dir.exists():
                self._available_batches = []
            else:
                self._available_batches = [
                    d.name for d in self.base_dir.iterdir()
                    if d.is_dir() and d.name.startswith("ecod_batch_")
                ]
        return self._available_batches

    def get_stable_batch_dir(self) -> Optional[str]:
        """Get stable batch directory for consistent testing"""
        # Try stable batch first
        if self.stable_batch in self.available_batches:
            return str(self.base_dir / self.stable_batch)

        # Fallback to any available batch
        if self.available_batches:
            return str(self.base_dir / self.available_batches[-1])

        return None

    def validate_test_data(self) -> Dict[str, bool]:
        """Validate test data files exist"""
        required_files = {
            "domain_lengths": self._test_data_dir / "domain_lengths.csv",
            "protein_lengths": self._test_data_dir / "protein_lengths.csv",
            "domain_definitions": self._test_data_dir / "domain_definitions.csv"
        }

        return {name: path.exists() for name, path in required_files.items()}


class ReferenceDataLoader:
    """Lazy-loading reference data to improve test performance"""

    def __init__(self, test_data_dir: Path):
        self.test_data_dir = test_data_dir
        self._cache = {}

    def get_domain_lengths(self) -> Dict[str, int]:
        """Get domain lengths with caching"""
        if 'domain_lengths' not in self._cache:
            path = self.test_data_dir / "domain_lengths.csv"
            if path.exists():
                self._cache['domain_lengths'] = load_reference_lengths(str(path))
            else:
                self._cache['domain_lengths'] = {}
        return self._cache['domain_lengths']

    def get_protein_lengths(self) -> Dict[tuple, int]:
        """Get protein lengths with caching"""
        if 'protein_lengths' not in self._cache:
            path = self.test_data_dir / "protein_lengths.csv"
            if path.exists():
                self._cache['protein_lengths'] = load_protein_lengths(str(path))
            else:
                self._cache['protein_lengths'] = {}
        return self._cache['protein_lengths']

    def get_domain_definitions(self, blacklist_path: Optional[str] = None) -> Dict[tuple, List]:
        """Get domain definitions with caching"""
        cache_key = f'domain_definitions_{blacklist_path or "none"}'
        if cache_key not in self._cache:
            path = self.test_data_dir / "domain_definitions.csv"
            if path.exists():
                self._cache[cache_key] = load_domain_definitions(
                    str(path),
                    blacklist_path=blacklist_path
                )
            else:
                self._cache[cache_key] = {}
        return self._cache[cache_key]


# Session-scoped fixtures for expensive setup
@pytest.fixture(scope="session", autouse=True)
def test_environment():
    """Detect and validate test environment"""
    env = TestEnvironment()

    # Check if we have any test environment at all
    if not env.available_batches and not env._test_data_dir.exists():
        pytest.fail("""
No test environment detected. Please set up test data:

For test data files:
    cd mini
    python range_cache_parser.py --output-dir test_data

For ECOD batch data:
    Ensure /data/ecod/pdb_updates/batches/ exists with batch directories
    Or set ECOD_BATCH_DIR environment variable
        """)

    return env


@pytest.fixture(scope="session")
def test_data_dir():
    """More flexible test data directory detection"""
    possible_paths = [
        str(Path(__file__).parent.parent / "test_data"),
        "/tmp/test_data",  # fallback
        os.environ.get('MINI_TEST_DATA', str(Path(__file__).parent.parent / "test_data"))
    ]

    for path in possible_paths:
        if os.path.exists(path):
            return path

    # Create minimal test data if none exists
    test_data = str(Path(__file__).parent.parent / "test_data")
    os.makedirs(test_data, exist_ok=True)
    return test_data


@pytest.fixture(scope="session")
def reference_data_loader(test_data_dir):
    """Lazy-loading reference data factory"""
    return ReferenceDataLoader(Path(test_data_dir))


@pytest.fixture(scope="session")
def stable_batch_dir(test_environment):
    """Stable batch directory for consistent testing"""
    batch_dir = test_environment.get_stable_batch_dir()
    if batch_dir is None:
        # Try environment variable fallback
        batch_dir = os.environ.get('ECOD_BATCH_DIR')
        if not batch_dir or not os.path.exists(batch_dir):
            pytest.skip("No ECOD batch directories available - set ECOD_BATCH_DIR environment variable")
    return batch_dir


@pytest.fixture(scope="session")
def blacklist_file(test_data_dir):
    """Reference blacklist file if it exists"""
    blacklist_path = Path(test_data_dir) / "reference_blacklist.csv"
    return str(blacklist_path) if blacklist_path.exists() else None


# Primary test case fixtures
@pytest.fixture
def primary_test_protein():
    """Primary test protein for canonical tests"""
    return "8ovp_A"


@pytest.fixture
def sample_protein():
    """Alias for primary_test_protein for backward compatibility"""
    return "8ovp_A"


# Mock data fixtures for unit tests
@pytest.fixture
def mock_evidence():
    """Mock evidence data for unit tests (no file I/O)"""
    return [
        Evidence(
            type="chain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("252-494"),
            confidence=0.9,
            evalue=1e-50,
            domain_id="6dgv_A",
            reference_length=238
        ),
        Evidence(
            type="chain_blast",
            source_pdb="2ia4",
            query_range=SequenceRange.parse("2-248,491-517"),
            confidence=0.8,
            evalue=1e-40,
            domain_id="2ia4_A",
            reference_length=508
        ),
        Evidence(
            type="domain_blast",
            source_pdb="2ia4",
            query_range=SequenceRange.parse("108-205"),
            confidence=0.75,
            evalue=1e-20,
            domain_id="e2ia4A1",
            reference_length=98
        ),
        Evidence(
            type="hhsearch",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("260-480"),
            confidence=0.85,
            evalue=1e-15,
            domain_id="6dgv_fold"
        )
    ]


@pytest.fixture
def mock_reference_data():
    """Mock reference data for unit tests"""
    return {
        'domain_lengths': {
            "e6dgvA1": 238,
            "e2ia4A1": 98,
            "e2ia4A2": 120
        },
        'protein_lengths': {
            ("6dgv", "A"): 238,
            ("2ia4", "A"): 508,
            ("8ovp", "A"): 569
        },
        'domain_definitions': {
            ("2ia4", "A"): [
                # Mock domain references would go here
            ]
        }
    }


# Real data fixtures for integration tests (lazy-loaded)
@pytest.fixture
def real_reference_data(reference_data_loader, blacklist_file):
    """Real reference data loaded on demand"""
    return {
        'domain_lengths': reference_data_loader.get_domain_lengths(),
        'protein_lengths': reference_data_loader.get_protein_lengths(),
        'domain_definitions': reference_data_loader.get_domain_definitions(blacklist_file)
    }


@pytest.fixture
def blast_alignments(primary_test_protein, stable_batch_dir):
    """BLAST alignments for the primary test protein"""
    parts = primary_test_protein.split('_')
    pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else 'A'

    blast_dir = os.path.join(stable_batch_dir, "blast", "chain")
    if os.path.exists(blast_dir):
        return load_chain_blast_alignments(blast_dir, pdb_id, chain_id)
    else:
        pytest.skip(f"BLAST directory not found: {blast_dir}")


@pytest.fixture
def domain_summary_path(primary_test_protein, stable_batch_dir):
    """Path to domain summary XML for primary test protein"""
    xml_path = os.path.join(
        stable_batch_dir,
        "domains",
        f"{primary_test_protein}.develop291.domain_summary.xml"
    )

    if not os.path.exists(xml_path):
        pytest.skip(f"Domain summary not found: {xml_path}")

    return xml_path


# Output fixtures
@pytest.fixture
def temp_output_dir(tmp_path):
    """Clean temporary directory for each test"""
    return str(tmp_path)


@pytest.fixture
def output_file_path(temp_output_dir, primary_test_protein):
    """Standard output file path for test results"""
    return os.path.join(temp_output_dir, f"{primary_test_protein}_test.domains.xml")


# Parametrized fixtures for comprehensive testing
@pytest.fixture(params=[True, False])
def with_decomposition(request):
    """Test both with and without decomposition"""
    return request.param


@pytest.fixture(params=["chain_blast", "domain_blast", "hhsearch"])
def evidence_type(request):
    """Test different evidence types"""
    return request.param


# Test markers configuration
def pytest_configure(config):
    """Configure custom test markers"""
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers",
        "integration: marks tests as integration tests requiring real data"
    )
    config.addinivalue_line(
        "markers",
        "unit: marks tests as fast unit tests with mock data"
    )
    config.addinivalue_line(
        "markers",
        "visualization: marks tests that require PyMOL or visualization tools"
    )
    config.addinivalue_line(
        "markers",
        "performance: marks tests that measure performance/timing"
    )
    config.addinivalue_line(
        "markers",
        "regression: marks tests that compare mini vs current engine performance"
    )


# Test setup and validation
def pytest_runtest_setup(item):
    """Enhanced setup with better error reporting"""
    # Skip integration tests if no real data available
    if "integration" in [mark.name for mark in item.iter_markers()]:
        test_data_dir = Path(__file__).parent.parent / "test_data"
        env = TestEnvironment()

        if not test_data_dir.exists() and not env.available_batches:
            pytest.skip("Integration tests require test data or ECOD batches")

        # Check for specific requirements
        if not env.available_batches:
            pytest.skip("Integration tests require ECOD batch directories")


# Smart test categorization
def pytest_collection_modifyitems(config, items):
    """Automatically categorize tests based on naming and dependencies"""
    for item in items:
        # Add markers based on test names
        if "visualization" in item.name or "pymol" in item.name:
            item.add_marker(pytest.mark.visualization)

        if "performance" in item.name or "benchmark" in item.name:
            item.add_marker(pytest.mark.performance)

        # Add markers based on fixtures used
        fixture_names = getattr(item, 'fixturenames', [])

        if any(name.startswith('real_') for name in fixture_names):
            item.add_marker(pytest.mark.integration)
        elif any(name.startswith('mock_') for name in fixture_names):
            item.add_marker(pytest.mark.unit)

        # Mark slow tests
        if any(name in fixture_names for name in ['blast_alignments', 'domain_summary_path']):
            item.add_marker(pytest.mark.slow)


# Test session hooks for reporting
def pytest_sessionstart(session):
    """Report test environment at start"""
    env = TestEnvironment()
    print(f"\n=== Mini PyECOD Test Environment ===")
    print(f"Available batches: {len(env.available_batches)}")
    if env.available_batches:
        print(f"Using batch: {env.get_stable_batch_dir()}")

    test_data_status = env.validate_test_data()
    print(f"Test data files: {sum(test_data_status.values())}/{len(test_data_status)} available")

    for name, exists in test_data_status.items():
        status = "✓" if exists else "✗"
        print(f"  {status} {name}")
    print("=" * 40)


def pytest_sessionfinish(session, exitstatus):
    """Report test session results"""
    if hasattr(session, 'testscollected'):
        print(f"\n=== Test Session Complete ===")
        print(f"Exit status: {exitstatus}")
        print("=" * 30)
