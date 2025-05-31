#!/usr/bin/env python3
"""
Diagnostic script to identify issues with algorithm versioning refactoring
"""

import os
import sys
import importlib.util
from pathlib import Path

def check_file_exists(path):
    """Check if a file exists and report status"""
    if Path(path).exists():
        print(f"✓ {path} exists")
        return True
    else:
        print(f"✗ {path} missing")
        return False

def check_import(module_path, class_names=None):
    """Try to import a module and optionally check for specific classes"""
    try:
        module = __import__(module_path, fromlist=[''])
        print(f"✓ Successfully imported {module_path}")
        
        if class_names:
            for class_name in class_names:
                if hasattr(module, class_name):
                    print(f"  ✓ {class_name} found in {module_path}")
                else:
                    print(f"  ✗ {class_name} missing from {module_path}")
        return True
    except ImportError as e:
        print(f"✗ Failed to import {module_path}: {e}")
        return False

def check_directory_structure():
    """Check if the expected directory structure exists"""
    print("DIRECTORY STRUCTURE CHECK")
    print("=" * 50)
    
    # Expected directories
    directories = [
        "ecod",
        "ecod/core",
        "ecod/evaluation", 
        "ecod/evaluation/algorithm_versions",
        "ecod/models",
        "ecod/models/pipeline",
        "ecod/pipelines",
        "ecod/pipelines/domain_analysis",
        "ecod/tests",
        "ecod/tests/unit",
        "scripts"
    ]
    
    for directory in directories:
        check_file_exists(directory)
    
    # Expected files
    files = [
        "ecod/__init__.py",
        "ecod/core/__init__.py",
        "ecod/core/context.py",
        "ecod/evaluation/__init__.py",
        "ecod/evaluation/algorithm_versions/__init__.py",
        "ecod/evaluation/algorithm_versions/manager.py",
        "ecod/models/__init__.py",
        "ecod/models/pipeline/__init__.py",
        "ecod/models/pipeline/partition.py",
        "ecod/config.py",
        "ecod/db.py",
        "scripts/algorithm_manager.py"
    ]
    
    for file_path in files:
        check_file_exists(file_path)

def check_imports():
    """Check if all expected imports work"""
    print("\nIMPORT CHECK")
    print("=" * 50)
    
    # Add current directory to path for testing
    sys.path.insert(0, os.getcwd())
    
    # Test imports that the unit test expects
    imports_to_test = [
        ("ecod.core.context", ["ApplicationContext"]),
        ("ecod.evaluation.algorithm_versions.manager", ["AlgorithmVersionManager", "AlgorithmVersion", "AlgorithmStatus"]),
        ("ecod.config", ["ConfigManager"]),
        ("ecod.db", ["DBManager"]),
        ("ecod.models.pipeline.partition", ["DomainPartitionResult"]),
        ("ecod.pipelines.domain_analysis.partition", ["DomainPartitionService"])
    ]
    
    for module_path, expected_classes in imports_to_test:
        check_import(module_path, expected_classes)

def suggest_fixes():
    """Suggest fixes based on what we found"""
    print("\nSUGGESTED FIXES")
    print("=" * 50)
    
    fixes = [
        "1. CREATE MISSING MODULE STRUCTURE:",
        "   mkdir -p ecod/evaluation/algorithm_versions",
        "   touch ecod/evaluation/__init__.py",
        "   touch ecod/evaluation/algorithm_versions/__init__.py",
        "",
        "2. EXTRACT CLASSES FROM curation_test_suite.py:",
        "   Move AlgorithmVersion, AlgorithmVersionManager, AlgorithmStatus",
        "   from curation_test_suite.py to ecod/evaluation/algorithm_versions/manager.py",
        "",
        "3. CREATE PROPER MODELS:",
        "   Move data classes to ecod/models/algorithm/ or similar",
        "",
        "4. FIX IMPORTS:",
        "   Update all imports to use the new module paths",
        "",
        "5. CREATE DATABASE TABLES:",
        "   Add algorithm versioning tables to your schema",
        "",
        "6. UPDATE CONTEXT:",
        "   Ensure ApplicationContext can handle algorithm versioning config"
    ]
    
    for fix in fixes:
        print(fix)

def create_missing_files():
    """Create basic missing files with minimal content"""
    print("\nCREATING MISSING FILES")
    print("=" * 50)
    
    # Create directories
    os.makedirs("ecod/evaluation/algorithm_versions", exist_ok=True)
    
    # Create __init__.py files
    init_files = [
        "ecod/evaluation/__init__.py",
        "ecod/evaluation/algorithm_versions/__init__.py"
    ]
    
    for init_file in init_files:
        if not Path(init_file).exists():
            with open(init_file, 'w') as f:
                f.write('"""Algorithm evaluation package"""\n')
            print(f"✓ Created {init_file}")
    
    # Create basic manager.py if it doesn't exist
    manager_file = "ecod/evaluation/algorithm_versions/manager.py"
    if not Path(manager_file).exists():
        with open(manager_file, 'w') as f:
            f.write('''"""
Algorithm Version Manager

This module should contain:
- AlgorithmVersion class
- AlgorithmVersionManager class  
- AlgorithmStatus enum

TODO: Move classes from curation_test_suite.py here
"""

from enum import Enum
from dataclasses import dataclass
from typing import Optional, List, Dict, Any
from datetime import datetime

class AlgorithmStatus(Enum):
    """Algorithm deployment status"""
    DEVELOPMENT = "development"
    TESTING = "testing"
    PRODUCTION = "production"
    DEPRECATED = "deprecated"

@dataclass
class AlgorithmVersion:
    """Basic algorithm version - TODO: implement fully"""
    version_id: str
    name: str = ""
    description: str = ""
    status: AlgorithmStatus = AlgorithmStatus.DEVELOPMENT
    created_at: Optional[datetime] = None
    
    @classmethod
    def from_config_file(cls, config_path: str):
        """TODO: Implement config file loading"""
        raise NotImplementedError("Config file loading not implemented yet")
    
    def to_config_dict(self) -> Dict[str, Any]:
        """TODO: Implement config dict conversion"""
        raise NotImplementedError("Config dict conversion not implemented yet")

class AlgorithmVersionManager:
    """Basic algorithm version manager - TODO: implement fully"""
    
    def __init__(self, context):
        self.context = context
        self.db = context.db if hasattr(context, 'db') else None
    
    def list_versions(self, status=None) -> List[AlgorithmVersion]:
        """TODO: Implement version listing"""
        return []
    
    def register_version(self, algorithm: AlgorithmVersion) -> int:
        """TODO: Implement version registration"""
        raise NotImplementedError("Version registration not implemented yet")
    
    def get_version(self, version_id: str) -> Optional[AlgorithmVersion]:
        """TODO: Implement version retrieval"""
        return None
    
    def promote_version(self, version_id: str, new_status: AlgorithmStatus) -> bool:
        """TODO: Implement version promotion"""
        raise NotImplementedError("Version promotion not implemented yet")
    
    def export_version(self, version_id: str, output_path: str):
        """TODO: Implement version export"""
        raise NotImplementedError("Version export not implemented yet")
    
    def import_version(self, config_path: str) -> int:
        """TODO: Implement version import"""
        raise NotImplementedError("Version import not implemented yet")
''')
        print(f"✓ Created {manager_file}")

def main():
    """Run full diagnostic"""
    print("PYECOD ALGORITHM VERSIONING DIAGNOSTIC")
    print("=" * 60)
    
    check_directory_structure()
    check_imports()
    
    print("\nCREATE MISSING FILES? (y/n): ", end="")
    try:
        response = input().lower().strip()
        if response == 'y':
            create_missing_files()
            print("\nRe-running import check after creating files...")
            check_imports()
    except KeyboardInterrupt:
        print("\nSkipping file creation")
    
    suggest_fixes()

if __name__ == "__main__":
    main()
