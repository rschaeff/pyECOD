#!/usr/bin/env python3
"""
Fixed Comprehensive Test Suite for Algorithm Versioning
"""

import pytest
import tempfile
import json
import os
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path

# Import modules under test
from ecod.core.context import ApplicationContext
from ecod.evaluation.algorithm_versions.manager import (
    AlgorithmVersionManager,
    AlgorithmVersion,
    AlgorithmStatus
)

class TestAlgorithmVersioning:
    """Comprehensive test suite for algorithm versioning"""

    @pytest.fixture
    def mock_context(self):
        """Create mock application context"""
        context = Mock(spec=ApplicationContext)
        context.db = Mock()
        context.config = Mock()
        # Mock database responses
        context.db.execute_query.return_value = [[1]]  # Mock ID return
        context.db.execute_dict_query.return_value = []  # Empty list for queries
        return context

    @pytest.fixture
    def sample_algorithm_data(self):
        """Sample algorithm data for testing"""
        return {
            'version_id': 'test_v1.0',
            'name': 'Test Algorithm',
            'description': 'Test algorithm for unit testing',
            'status': 'development',
            'partition_config': {'min_domain_size': 20},
            'evidence_weights': {'hhsearch': 3.0},
            'coverage_thresholds': {'min_reference_coverage': 0.7},
            'behavioral_flags': {'prefer_hhsearch_classification': True}
        }

    @pytest.fixture
    def temp_config_file(self, sample_algorithm_data):
        """Create temporary algorithm config file"""
        config_data = {
            'algorithm': {
                'version_id': sample_algorithm_data['version_id'],
                'name': sample_algorithm_data['name'],
                'description': sample_algorithm_data['description'],
                'status': sample_algorithm_data['status'],
                'created_by': 'test_suite'
            },
            'domain_analysis': {
                'partition': sample_algorithm_data['partition_config'],
                'evidence_weights': sample_algorithm_data['evidence_weights'],
                'coverage_thresholds': sample_algorithm_data['coverage_thresholds'],
                'behavioral_flags': sample_algorithm_data['behavioral_flags']
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
            import yaml
            yaml.dump(config_data, f)
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    # AlgorithmStatus tests
    def test_algorithm_status_enum_values(self):
        """Test AlgorithmStatus enum has correct values"""
        assert AlgorithmStatus.DEVELOPMENT.value == "development"
        assert AlgorithmStatus.TESTING.value == "testing"
        assert AlgorithmStatus.PRODUCTION.value == "production"
        assert AlgorithmStatus.DEPRECATED.value == "deprecated"

    def test_algorithm_status_from_string(self):
        """Test creating AlgorithmStatus from string"""
        status = AlgorithmStatus("development")
        assert status == AlgorithmStatus.DEVELOPMENT

        with pytest.raises(ValueError):
            AlgorithmStatus("invalid_status")

    # AlgorithmVersion tests
    def test_algorithm_version_creation(self, sample_algorithm_data):
        """Test AlgorithmVersion creation with required parameters"""
        algorithm = AlgorithmVersion(
            version_id=sample_algorithm_data['version_id'],
            name=sample_algorithm_data['name'],
            description=sample_algorithm_data['description']
        )

        assert algorithm.version_id == sample_algorithm_data['version_id']
        assert algorithm.name == sample_algorithm_data['name']
        assert algorithm.description == sample_algorithm_data['description']
        assert algorithm.status == AlgorithmStatus.DEVELOPMENT  # Default

    def test_algorithm_version_with_all_parameters(self, sample_algorithm_data):
        """Test AlgorithmVersion creation with all parameters"""
        algorithm = AlgorithmVersion(
            version_id=sample_algorithm_data['version_id'],
            name=sample_algorithm_data['name'],
            description=sample_algorithm_data['description'],
            status=AlgorithmStatus.TESTING,
            partition_config=sample_algorithm_data['partition_config'],
            evidence_weights=sample_algorithm_data['evidence_weights'],
            coverage_thresholds=sample_algorithm_data['coverage_thresholds'],
            behavioral_flags=sample_algorithm_data['behavioral_flags']
        )

        assert algorithm.status == AlgorithmStatus.TESTING
        assert algorithm.partition_config == sample_algorithm_data['partition_config']
        assert algorithm.evidence_weights == sample_algorithm_data['evidence_weights']

    def test_algorithm_version_validation(self):
        """Test AlgorithmVersion validation"""
        # Test missing version_id
        with pytest.raises(Exception):  # ValidationError
            AlgorithmVersion(version_id="", name="Test", description="Test")

        # Test missing name
        with pytest.raises(Exception):  # ValidationError
            AlgorithmVersion(version_id="test", name="", description="Test")

    def test_algorithm_version_to_config_dict(self, sample_algorithm_data):
        """Test AlgorithmVersion.to_config_dict()"""
        algorithm = AlgorithmVersion(
            version_id=sample_algorithm_data['version_id'],
            name=sample_algorithm_data['name'],
            description=sample_algorithm_data['description'],
            partition_config=sample_algorithm_data['partition_config'],
            evidence_weights=sample_algorithm_data['evidence_weights']
        )

        config_dict = algorithm.to_config_dict()

        assert config_dict['version_id'] == sample_algorithm_data['version_id']
        assert 'domain_analysis' in config_dict
        assert config_dict['domain_analysis']['partition'] == sample_algorithm_data['partition_config']

    def test_algorithm_version_from_config_file(self, temp_config_file):
        """Test AlgorithmVersion.from_config_file()"""
        algorithm = AlgorithmVersion.from_config_file(temp_config_file)

        assert algorithm.version_id == 'test_v1.0'
        assert algorithm.name == 'Test Algorithm'
        assert algorithm.description == 'Test algorithm for unit testing'
        assert algorithm.status == AlgorithmStatus.DEVELOPMENT

    def test_algorithm_version_from_config_file_missing_file(self):
        """Test AlgorithmVersion.from_config_file() with missing file"""
        with pytest.raises(Exception):  # ConfigurationError
            AlgorithmVersion.from_config_file('/nonexistent/file.yml')

    # AlgorithmVersionManager tests
    def test_algorithm_version_manager_creation(self, mock_context):
        """Test AlgorithmVersionManager creation"""
        manager = AlgorithmVersionManager(mock_context)

        assert manager.context == mock_context
        assert manager.db == mock_context.db
        assert hasattr(manager, 'logger')

    def test_algorithm_version_manager_list_versions_empty(self, mock_context):
        """Test listing versions when database is empty"""
        mock_context.db.execute_dict_query.return_value = []

        manager = AlgorithmVersionManager(mock_context)
        versions = manager.list_versions()

        assert isinstance(versions, list)
        assert len(versions) == 0

    def test_algorithm_version_manager_list_versions_with_data(self, mock_context):
        """Test listing versions with mock data"""
        mock_data = [{
            'id': 1,
            'version_id': 'test_v1.0',
            'name': 'Test Algorithm',
            'description': 'Test description',
            'status': 'development',
            'config_data': '{"test": true}',
            'created_at': '2024-01-01T00:00:00',
            'created_by': 'test_user',
            'parent_version_id_str': None
        }]
        mock_context.db.execute_dict_query.return_value = mock_data

        manager = AlgorithmVersionManager(mock_context)
        versions = manager.list_versions()

        assert len(versions) == 1
        assert versions[0].version_id == 'test_v1.0'
        assert versions[0].name == 'Test Algorithm'

    def test_algorithm_version_manager_register_version(self, mock_context, sample_algorithm_data):
        """Test registering a new algorithm version"""
        # Mock the get_version to return None (doesn't exist)
        mock_context.db.execute_dict_query.return_value = []
        # Mock the insert to return an ID
        mock_context.db.execute_query.return_value = [[42]]

        manager = AlgorithmVersionManager(mock_context)
        algorithm = AlgorithmVersion(
            version_id=sample_algorithm_data['version_id'],
            name=sample_algorithm_data['name'],
            description=sample_algorithm_data['description']
        )

        algorithm_id = manager.register_version(algorithm)

        assert algorithm_id == 42
        assert algorithm.database_id == 42
        # Verify the database was called
        mock_context.db.execute_query.assert_called()

    def test_algorithm_version_manager_get_version_not_found(self, mock_context):
        """Test getting a version that doesn't exist"""
        mock_context.db.execute_dict_query.return_value = []

        manager = AlgorithmVersionManager(mock_context)
        result = manager.get_version('nonexistent')

        assert result is None

    def test_algorithm_version_manager_promote_version(self, mock_context):
        """Test promoting an algorithm version"""
        # Mock getting the existing version
        mock_data = [{
            'id': 1,
            'version_id': 'test_v1.0',
            'name': 'Test Algorithm',
            'description': 'Test description',
            'status': 'development',
            'config_data': '{"test": true}',
            'created_at': '2024-01-01T00:00:00',
            'created_by': 'test_user',
            'parent_version_id_str': None
        }]
        mock_context.db.execute_dict_query.return_value = mock_data

        manager = AlgorithmVersionManager(mock_context)

        result = manager.promote_version('test_v1.0', AlgorithmStatus.TESTING)

        assert result is True
        # Verify update query was called
        mock_context.db.execute_query.assert_called()

class TestAlgorithmVersioningIntegration:
    """Integration tests for algorithm versioning"""

    @pytest.fixture
    def temp_config_file(self):
        """Create temporary algorithm config file"""
        config_data = {
            'algorithm': {
                'version_id': 'integration_test_v1.0',
                'name': 'Integration Test Algorithm',
                'description': 'Algorithm for integration testing',
                'status': 'development',
                'created_by': 'test_suite'
            },
            'domain_analysis': {
                'partition': {'min_domain_size': 20},
                'evidence_weights': {'hhsearch': 3.0},
                'coverage_thresholds': {'min_reference_coverage': 0.7},
                'behavioral_flags': {'prefer_hhsearch_classification': True}
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
            import yaml
            yaml.dump(config_data, f)
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    def test_config_file_round_trip(self, temp_config_file):
        """Test loading config file and converting back to config dict"""
        # Load from file
        algorithm = AlgorithmVersion.from_config_file(temp_config_file)

        # Convert back to config dict
        config_dict = algorithm.to_config_dict()

        # Verify round trip
        assert config_dict['version_id'] == 'integration_test_v1.0'
        assert 'domain_analysis' in config_dict
        assert config_dict['domain_analysis']['partition']['min_domain_size'] == 20

class TestAlgorithmVersioningValidation:
    """Tests for validation and error handling"""

    def test_invalid_status_string(self):
        """Test handling of invalid status strings"""
        with pytest.raises(ValueError):
            AlgorithmStatus("invalid_status")

    def test_algorithm_version_empty_required_fields(self):
        """Test validation of required fields"""
        with pytest.raises(Exception):  # ValidationError
            AlgorithmVersion(version_id="", name="Test", description="Test")

        with pytest.raises(Exception):  # ValidationError
            AlgorithmVersion(version_id="test", name="", description="Test")

class TestAlgorithmVersioningEdgeCases:
    """Tests for edge cases and boundary conditions"""

    def test_algorithm_version_with_special_characters(self):
        """Test algorithm versions with special characters"""
        algorithm = AlgorithmVersion(
            version_id="test-v1.0_special.chars",
            name="Test Algorithm (Special)",
            description="Test with special characters: !@#$%^&*()"
        )

        assert algorithm.version_id == "test-v1.0_special.chars"
        assert "special characters" in algorithm.description

    def test_algorithm_version_unicode_support(self):
        """Test unicode support in algorithm versions"""
        algorithm = AlgorithmVersion(
            version_id="test_unicode",
            name="Test Algorithμ",  # Greek mu
            description="Test with unicode: αβγδε"
        )

        assert "Algorithμ" in algorithm.name
        assert "αβγδε" in algorithm.description

# Utility functions for test helpers
def create_test_algorithm(version_id: str = "test_v1.0") -> AlgorithmVersion:
    """Helper function to create test algorithm instances"""
    return AlgorithmVersion(
        version_id=version_id,
        name=f"Test Algorithm {version_id}",
        description="Test algorithm for unit testing",
        partition_config={'min_domain_size': 20},
        evidence_weights={'hhsearch': 3.0},
        coverage_thresholds={'min_reference_coverage': 0.7},
        behavioral_flags={'prefer_hhsearch_classification': True}
    )

def create_mock_database_response(algorithm_data: dict) -> dict:
    """Helper function to create mock database responses"""
    return {
        'id': 1,
        'version_id': algorithm_data.get('version_id', 'test_v1.0'),
        'name': algorithm_data.get('name', 'Test Algorithm'),
        'description': algorithm_data.get('description', 'Test description'),
        'status': algorithm_data.get('status', 'development'),
        'config_data': json.dumps(algorithm_data.get('config', {})),
        'created_at': '2024-01-01T00:00:00',
        'created_by': 'test_user',
        'parent_version_id_str': None
    }

if __name__ == "__main__":
    # Run specific test groups
    import sys

    if len(sys.argv) > 1:
        test_group = sys.argv[1]
        if test_group == "basic":
            pytest.main([__file__ + "::TestAlgorithmVersioning", "-v"])
        elif test_group == "integration":
            pytest.main([__file__ + "::TestAlgorithmVersioningIntegration", "-v"])
        elif test_group == "validation":
            pytest.main([__file__ + "::TestAlgorithmVersioningValidation", "-v"])
        elif test_group == "edge_cases":
            pytest.main([__file__ + "::TestAlgorithmVersioningEdgeCases", "-v"])
        else:
            print(f"Unknown test group: {test_group}")
            sys.exit(1)
    else:
        # Run all tests
        pytest.main([__file__, "-v"])
