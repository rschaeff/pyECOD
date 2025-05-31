#!/usr/bin/env python3
"""
Helper script to extract algorithm versioning classes from curation_test_suite.py
and properly structure them in the ecod package
"""

import os
import re
import ast
from pathlib import Path

def extract_class_from_file(file_path, class_name):
    """Extract a class definition from a Python file"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    tree = ast.parse(content)
    
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            # Get the source code for this class
            lines = content.split('\n')
            start_line = node.lineno - 1
            
            # Find the end of the class
            end_line = start_line + 1
            indent_level = len(lines[start_line]) - len(lines[start_line].lstrip())
            
            for i in range(start_line + 1, len(lines)):
                line = lines[i]
                if line.strip() == "":
                    continue
                current_indent = len(line) - len(line.lstrip())
                if current_indent <= indent_level and line.strip():
                    break
                end_line = i + 1
            
            return '\n'.join(lines[start_line:end_line])
    
    return None

def extract_function_from_file(file_path, function_name):
    """Extract a function definition from a Python file"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    tree = ast.parse(content)
    
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == function_name:
            lines = content.split('\n')
            start_line = node.lineno - 1
            
            # Find the end of the function
            end_line = start_line + 1
            indent_level = len(lines[start_line]) - len(lines[start_line].lstrip())
            
            for i in range(start_line + 1, len(lines)):
                line = lines[i]
                if line.strip() == "":
                    continue
                current_indent = len(line) - len(line.lstrip())
                if current_indent <= indent_level and line.strip():
                    break
                end_line = i + 1
            
            return '\n'.join(lines[start_line:end_line])
    
    return None

def create_proper_manager_module():
    """Create the proper algorithm versions manager module"""
    
    # Extract classes from curation_test_suite.py
    curation_file = "scripts/curation_test_suite.py"
    
    if not Path(curation_file).exists():
        print(f"Warning: {curation_file} not found")
        return
    
    # Create the manager.py content
    manager_content = '''"""
Algorithm Version Management for pyECOD

This module provides version management for domain partitioning algorithms,
supporting systematic evaluation and deployment of algorithm improvements.
"""

import os
import json
import yaml
import logging
from enum import Enum
from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

'''
    
    # Extract the enums and classes we need
    classes_to_extract = [
        "AlgorithmIterationType",
        "AlgorithmVersion", 
        "AlgorithmVersionManager"
    ]
    
    functions_to_extract = [
        "create_baseline_algorithm",
        "create_coverage_focused_algorithm",
        "create_chain_blast_priority_algorithm", 
        "create_chain_blast_only_algorithm"
    ]
    
    for class_name in classes_to_extract:
        class_code = extract_class_from_file(curation_file, class_name)
        if class_code:
            manager_content += "\n" + class_code + "\n"
        else:
            print(f"Warning: Could not extract {class_name}")
    
    for func_name in functions_to_extract:
        func_code = extract_function_from_file(curation_file, func_name)
        if func_code:
            manager_content += "\n" + func_code + "\n"
        else:
            print(f"Warning: Could not extract {func_name}")
    
    # Create the directory and file
    os.makedirs("ecod/evaluation/algorithm_versions", exist_ok=True)
    
    with open("ecod/evaluation/algorithm_versions/manager.py", "w") as f:
        f.write(manager_content)
    
    print("✓ Created ecod/evaluation/algorithm_versions/manager.py")

def create_models_module():
    """Create a separate models module for algorithm data structures"""
    
    models_content = '''"""
Algorithm Version Data Models

Core data structures for algorithm versioning and evaluation.
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Any, Optional
from enum import Enum

@dataclass
class CurationDecision:
    """Manual curation decision for a protein"""
    protein_id: int
    source_id: str
    pdb_id: str
    chain_id: str
    has_domain: bool
    domain_assigned_correctly: Optional[bool]
    boundaries_correct: Optional[bool]
    is_fragment: bool
    is_repeat_protein: bool
    confidence_level: int
    primary_evidence_type: Optional[str]
    reference_domain_id: Optional[str] 
    reference_pdb_id: Optional[str]
    reference_chain_id: Optional[str]
    evidence_confidence: Optional[float]
    evidence_evalue: Optional[float]
    curator_name: str
    session_id: int
    created_at: datetime

@dataclass 
class PartitionResult:
    """Automated domain partition result"""
    protein_id: int
    source_id: str
    pdb_id: str
    chain_id: str
    is_classified: bool
    is_peptide: bool
    domain_count: int
    domains: List[Dict[str, Any]]
    coverage: float
    sequence_length: int
    confidence_scores: List[float]
    algorithm_version: str
    processing_timestamp: datetime

@dataclass
class TestSet:
    """A test set of curated proteins"""
    test_set_id: int
    name: str
    description: str
    created_at: datetime
    protein_count: int
    curator_breakdown: Dict[str, int]
    decision_breakdown: Dict[str, int]
    proteins: List[CurationDecision]

@dataclass
class ComparisonMetrics:
    """Metrics comparing automated results to manual curation"""
    test_set_id: int
    algorithm_version: str
    
    # Primary metrics (focus areas)
    domain_presence_accuracy: float
    fragment_detection_accuracy: float
    boundary_agreement_rate: float
    
    # Discontinuous domain metrics
    discontinuous_domain_detection: float
    discontinuous_boundary_accuracy: float
    coverage_cutoff_effectiveness: float
    
    # Boundary-specific metrics
    exact_boundary_matches: float
    boundary_tolerance_5: float
    boundary_tolerance_10: float
    boundary_over_segmentation: float
    boundary_under_segmentation: float
    
    # Fragment/peptide detection
    peptide_precision: float
    peptide_recall: float
    fragment_vs_domain_accuracy: float
    
    # Secondary metrics
    domain_count_accuracy: float
    classification_agreement_rate: float
    
    # Detailed breakdowns
    confusion_matrix: Dict[str, int]
    boundary_error_distribution: Dict[str, List[int]]
    improvement_cases: List[Dict[str, Any]]
    regression_cases: List[Dict[str, Any]]
    discontinuous_cases: List[Dict[str, Any]]
'''
    
    os.makedirs("ecod/evaluation/algorithm_versions", exist_ok=True)
    
    with open("ecod/evaluation/algorithm_versions/models.py", "w") as f:
        f.write(models_content)
    
    print("✓ Created ecod/evaluation/algorithm_versions/models.py")

def create_database_schema():
    """Create database schema for algorithm versioning"""
    
    schema_content = '''-- Algorithm Versioning Tables for pyECOD

-- Algorithm versions table
CREATE TABLE IF NOT EXISTS pdb_analysis.algorithm_versions (
    id SERIAL PRIMARY KEY,
    version_id VARCHAR(100) UNIQUE NOT NULL,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    iteration_type VARCHAR(50) NOT NULL,
    parent_version VARCHAR(100),
    created_by VARCHAR(100),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    deployment_status VARCHAR(20) DEFAULT 'development',
    
    -- Configuration as JSON
    partition_options JSONB,
    evidence_weights JSONB,
    evidence_filters JSONB,
    coverage_thresholds JSONB,
    boundary_settings JSONB,
    performance_settings JSONB,
    behavioral_flags JSONB,
    
    -- Metadata
    changes_from_parent TEXT[],
    improvement_targets TEXT[],
    known_limitations TEXT[],
    test_results JSONB,
    notes TEXT
);

-- Algorithm test runs
CREATE TABLE IF NOT EXISTS pdb_analysis.algorithm_test_runs (
    id SERIAL PRIMARY KEY,
    version_id VARCHAR(100) REFERENCES pdb_analysis.algorithm_versions(version_id),
    test_set_id INTEGER,
    run_type VARCHAR(50) NOT NULL, -- 'quick_test', 'full_evaluation', etc.
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP,
    proteins_tested INTEGER,
    successful_runs INTEGER,
    status VARCHAR(20) DEFAULT 'running', -- 'running', 'completed', 'failed'
    test_parameters JSONB,
    results_summary JSONB
);

-- Indexes for performance
CREATE INDEX IF NOT EXISTS idx_algorithm_versions_status ON pdb_analysis.algorithm_versions(deployment_status);
CREATE INDEX IF NOT EXISTS idx_algorithm_versions_created ON pdb_analysis.algorithm_versions(created_at);
CREATE INDEX IF NOT EXISTS idx_algorithm_test_runs_version ON pdb_analysis.algorithm_test_runs(version_id);
CREATE INDEX IF NOT EXISTS idx_algorithm_test_runs_test_set ON pdb_analysis.algorithm_test_runs(test_set_id);
'''
    
    with open("sql/algorithm_versioning_schema.sql", "w") as f:
        f.write(schema_content)
    
    print("✓ Created sql/algorithm_versioning_schema.sql")

def create_enhanced_tests():
    """Create enhanced unit tests with better coverage"""
    
    test_content = '''#!/usr/bin/env python3
"""
Enhanced unit tests for algorithm version management
"""

import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from ecod.core.context import ApplicationContext
from ecod.evaluation.algorithm_versions.manager import (
    AlgorithmVersionManager, 
    AlgorithmVersion, 
    AlgorithmStatus,
    AlgorithmIterationType,
    create_baseline_algorithm
)

class TestAlgorithmVersion:
    """Test AlgorithmVersion class"""
    
    def test_algorithm_version_creation(self):
        """Test basic algorithm version creation"""
        algorithm = AlgorithmVersion(
            version_id="test_v1.0",
            iteration_type=AlgorithmIterationType.BASELINE,
            name="Test Algorithm",
            description="Test algorithm for unit testing"
        )
        
        assert algorithm.version_id == "test_v1.0"
        assert algorithm.iteration_type == AlgorithmIterationType.BASELINE
        assert algorithm.name == "Test Algorithm"
        assert algorithm.deployment_status == "development"
    
    def test_baseline_algorithm_creation(self):
        """Test baseline algorithm factory function"""
        baseline = create_baseline_algorithm()
        
        assert baseline.version_id == "v1.0_baseline"
        assert baseline.iteration_type == AlgorithmIterationType.BASELINE
        assert "Original Baseline" in baseline.name
        assert baseline.deployment_status == "development"
        assert len(baseline.known_limitations) > 0
    
    @patch('builtins.open')
    @patch('yaml.safe_load')
    def test_from_config_file(self, mock_yaml_load, mock_open):
        """Test loading algorithm from config file"""
        mock_config = {
            'algorithm': {
                'version_id': 'test_config',
                'name': 'Config Test',
                'description': 'Test from config'
            }
        }
        mock_yaml_load.return_value = mock_config
        
        # This will fail until from_config_file is implemented
        with pytest.raises(NotImplementedError):
            AlgorithmVersion.from_config_file("test.yml")

class TestAlgorithmVersionManager:
    """Test AlgorithmVersionManager class"""
    
    @pytest.fixture
    def mock_context(self):
        """Create mock application context"""
        context = Mock()
        context.db = Mock()
        context.config = Mock()
        return context
    
    @pytest.fixture
    def manager(self, mock_context):
        """Create algorithm version manager"""
        return AlgorithmVersionManager(mock_context)
    
    def test_manager_creation(self, manager):
        """Test manager creation"""
        assert manager is not None
        assert hasattr(manager, 'context')
        assert hasattr(manager, 'db')
    
    def test_list_versions_empty(self, manager):
        """Test listing versions when none exist"""
        # This will fail until list_versions is implemented
        versions = manager.list_versions()
        assert isinstance(versions, list)
        assert len(versions) == 0
    
    def test_register_version_not_implemented(self, manager):
        """Test that register_version raises NotImplementedError"""
        algorithm = create_baseline_algorithm()
        
        with pytest.raises(NotImplementedError):
            manager.register_version(algorithm)
    
    def test_get_version_not_implemented(self, manager):
        """Test that get_version returns None (not implemented)"""
        result = manager.get_version("nonexistent")
        assert result is None

class TestAlgorithmVersionIntegration:
    """Integration tests for algorithm versioning"""
    
    @pytest.fixture
    def temp_config(self):
        """Create temporary config file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
            f.write("""
database:
  host: localhost
  port: 5432
  database: test_db
  user: test_user
  password: test_pass

logging:
  level: INFO
""")
            temp_path = f.name
        
        yield temp_path
        os.unlink(temp_path)
    
    @patch('ecod.core.context.ApplicationContext')
    def test_integration_with_context(self, mock_context_class, temp_config):
        """Test integration with ApplicationContext"""
        mock_context = Mock()
        mock_context.db = Mock()
        mock_context_class.return_value = mock_context
        
        # This test will help identify if ApplicationContext integration works
        context = ApplicationContext(temp_config)
        manager = AlgorithmVersionManager(context)
        
        assert manager.context == context
        assert manager.db == context.db

class TestDatabaseIntegration:
    """Test database integration for algorithm versioning"""
    
    @pytest.fixture
    def mock_db(self):
        """Create mock database connection"""
        db = Mock()
        db.execute_query = Mock()
        db.execute_dict_query = Mock()
        return db
    
    def test_algorithm_version_storage(self, mock_db):
        """Test storing algorithm version in database"""
        # This test will be implemented once database methods are added
        pass
    
    def test_algorithm_version_retrieval(self, mock_db):
        """Test retrieving algorithm version from database"""
        # This test will be implemented once database methods are added
        pass

# Helper functions for running specific test groups
def run_basic_tests():
    """Run just the basic functionality tests"""
    pytest.main([
        __file__ + "::TestAlgorithmVersion::test_algorithm_version_creation",
        __file__ + "::TestAlgorithmVersion::test_baseline_algorithm_creation",
        __file__ + "::TestAlgorithmVersionManager::test_manager_creation",
        "-v"
    ])

def run_implementation_tests():
    """Run tests that check for NotImplementedError (expected to fail initially)"""
    pytest.main([
        __file__ + "::TestAlgorithmVersion::test_from_config_file",
        __file__ + "::TestAlgorithmVersionManager::test_register_version_not_implemented",
        "-v"
    ])

if __name__ == "__main__":
    # Run all tests
    pytest.main([__file__, "-v"])
'''
    
    os.makedirs("ecod/tests/unit", exist_ok=True)
    
    with open("ecod/tests/unit/test_algorithm_versioning_enhanced.py", "w") as f:
        f.write(test_content)
    
    print("✓ Created ecod/tests/unit/test_algorithm_versioning_enhanced.py")

def update_imports_in_scripts():
    """Update imports in existing scripts to use the new package structure"""
    
    # Update algorithm_manager.py
    manager_script = "scripts/algorithm_manager.py"
    if Path(manager_script).exists():
        with open(manager_script, 'r') as f:
            content = f.read()
        
        # Update imports
        updated_content = content.replace(
            "from ecod.evaluation.algorithm_versions.manager import (",
            "from ecod.evaluation.algorithm_versions.manager import ("
        )
        
        # Add proper imports if they're missing
        if "from ecod.evaluation.algorithm_versions.manager import" not in updated_content:
            # Add the import after the sys.path.insert line
            import_line = "from ecod.evaluation.algorithm_versions.manager import (\n    AlgorithmVersionManager, \n    AlgorithmVersion, \n    AlgorithmStatus\n)\n"
            lines = updated_content.split('\n')
            
            for i, line in enumerate(lines):
                if "sys.path.insert" in line:
                    lines.insert(i + 2, import_line)
                    break
            
            updated_content = '\n'.join(lines)
        
        with open(manager_script, 'w') as f:
            f.write(updated_content)
        
        print(f"✓ Updated imports in {manager_script}")

def main():
    """Run the refactoring process"""
    print("PYECOD ALGORITHM VERSIONING REFACTORING")
    print("=" * 50)
    
    print("\n1. Creating proper module structure...")
    create_proper_manager_module()
    
    print("\n2. Creating models module...")
    create_models_module()
    
    print("\n3. Creating database schema...")
    os.makedirs("sql", exist_ok=True)
    create_database_schema()
    
    print("\n4. Creating enhanced tests...")
    create_enhanced_tests()
    
    print("\n5. Updating script imports...")
    update_imports_in_scripts()
    
    print("\n6. Creating __init__.py files...")
    init_files = [
        "ecod/__init__.py",
        "ecod/evaluation/__init__.py", 
        "ecod/evaluation/algorithm_versions/__init__.py"
    ]
    
    for init_file in init_files:
        os.makedirs(os.path.dirname(init_file), exist_ok=True)
        if not Path(init_file).exists():
            with open(init_file, 'w') as f:
                f.write(f'"""{"Algorithm evaluation" if "evaluation" in init_file else "pyECOD package"}"""\n')
            print(f"✓ Created {init_file}")
    
    print("\nREFACTORING COMPLETE!")
    print("=" * 30)
    print("\nNext steps:")
    print("1. Run: python diagnostic_script.py")
    print("2. Run: pytest ecod/tests/unit/test_algorithm_versioning_enhanced.py -v")
    print("3. Implement the NotImplementedError methods in manager.py")
    print("4. Add database schema to your PostgreSQL database")
    print("5. Test the algorithm_manager.py script")

if __name__ == "__main__":
    main()
