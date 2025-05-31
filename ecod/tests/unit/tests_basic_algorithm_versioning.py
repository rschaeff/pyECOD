#!/usr/bin/env python3
"""
Validation script for algorithm version management

This script tests the basic functionality to ensure everything is working correctly.
"""

import os
import sys
import tempfile
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

def test_imports():
    """Test that all required modules can be imported"""
    print("Testing imports...")
    
    try:
        from ecod.core.context import ApplicationContext
        print("  ✓ ApplicationContext imported successfully")
    except ImportError as e:
        print(f"  ✗ Failed to import ApplicationContext: {e}")
        return False
    
    try:
        from ecod.evaluation.algorithm_versions.manager import (
            AlgorithmVersionManager, 
            AlgorithmVersion, 
            AlgorithmStatus
        )
        print("  ✓ Algorithm management classes imported successfully")
    except ImportError as e:
        print(f"  ✗ Failed to import algorithm management: {e}")
        return False
    
    return True

def test_context_creation():
    """Test context creation"""
    print("\nTesting context creation...")
    
    try:
        context = ApplicationContext("config/config.yml")
        print("  ✓ Context created successfully")
        return context
    except Exception as e:
        print(f"  ✗ Failed to create context: {e}")
        return None

def test_manager_creation(context):
    """Test manager creation"""
    print("\nTesting manager creation...")
    
    try:
        from ecod.evaluation.algorithm_versions.manager import AlgorithmVersionManager
        manager = AlgorithmVersionManager(context)
        print("  ✓ Algorithm manager created successfully")
        return manager
    except Exception as e:
        print(f"  ✗ Failed to create manager: {e}")
        return None

def test_algorithm_operations(manager):
    """Test basic algorithm operations"""
    print("\nTesting algorithm operations...")
    
    try:
        from ecod.evaluation.algorithm_versions.manager import AlgorithmVersion, AlgorithmStatus
        
        # Test listing (should work even if empty)
        algorithms = manager.list_versions()
        print(f"  ✓ Listed algorithms: {len(algorithms)} found")
        
        # Create a test algorithm
        test_config = """
algorithm:
  version_id: "test_validation"
  name: "Test Algorithm"
  description: "Test algorithm for validation"
  status: "development"
  created_by: "validation_script"

domain_analysis:
  partition:
    min_domain_size: 20
  evidence_weights:
    hhsearch: 3.0
  coverage_thresholds:
    min_reference_coverage: 0.7
  behavioral_flags:
    prefer_hhsearch_classification: true
"""
        
        # Write to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
            f.write(test_config)
            temp_path = f.name
        
        try:
            # Test loading from file
            algorithm = AlgorithmVersion.from_config_file(temp_path)
            print("  ✓ Algorithm loaded from config file")
            
            # Test registration
            algorithm_id = manager.register_version(algorithm)
            print(f"  ✓ Algorithm registered with ID: {algorithm_id}")
            
            # Test retrieval
            retrieved = manager.get_version("test_validation")
            if retrieved:
                print("  ✓ Algorithm retrieved successfully")
            else:
                print("  ✗ Failed to retrieve algorithm")
                return False
            
            # Test promotion
            success = manager.promote_version("test_validation", AlgorithmStatus.TESTING)
            if success:
                print("  ✓ Algorithm promoted successfully")
            else:
                print("  ✗ Failed to promote algorithm")
                return False
            
            # Test export
            with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
                export_path = f.name
            
            manager.export_version("test_validation", export_path)
            if Path(export_path).exists():
                print("  ✓ Algorithm exported successfully")
                os.unlink(export_path)  # Clean up
            else:
                print("  ✗ Failed to export algorithm")
                return False
            
            return True
            
        finally:
            # Clean up temp file
            os.unlink(temp_path)
            
    except Exception as e:
        print(f"  ✗ Algorithm operations failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_standalone_script():
    """Test the standalone script"""
    print("\nTesting standalone script...")
    
    script_path = "scripts/algorithm_manager.py"
    if not Path(script_path).exists():
        print(f"  ✗ Standalone script not found at {script_path}")
        return False
    
    try:
        # Test help command
        result = os.system(f"python {script_path} --help > /dev/null 2>&1")
        if result == 0:
            print("  ✓ Standalone script help works")
        else:
            print("  ✗ Standalone script help failed")
            return False
        
        # Test list command
        result = os.system(f"python {script_path} list > /dev/null 2>&1")
        if result == 0:
            print("  ✓ Standalone script list works")
        else:
            print("  ✗ Standalone script list failed")
            return False
        
        return True
        
    except Exception as e:
        print(f"  ✗ Standalone script test failed: {e}")
        return False

def main():
    """Run all validation tests"""
    print("Algorithm Version Management Validation")
    print("=" * 50)
    
    success_count = 0
    total_tests = 5
    
    # Test 1: Imports
    if test_imports():
        success_count += 1
    
    # Test 2: Context creation
    context = test_context_creation()
    if context:
        success_count += 1
    else:
        print("\n❌ Cannot continue without valid context")
        return 1
    
    # Test 3: Manager creation
    manager = test_manager_creation(context)
    if manager:
        success_count += 1
    else:
        print("\n❌ Cannot continue without valid manager")
        return 1
    
    # Test 4: Algorithm operations
    if test_algorithm_operations(manager):
        success_count += 1
    
    # Test 5: Standalone script
    if test_standalone_script():
        success_count += 1
    
    # Summary
    print(f"\nValidation Summary")
    print("=" * 20)
    print(f"Passed: {success_count}/{total_tests} tests")
    
    if success_count == total_tests:
        print("\n🎉 All tests passed! Algorithm version management is working correctly.")
        print("\nNext steps:")
        print("1. Register your baseline algorithm:")
        print("   python scripts/algorithm_manager.py register config/algorithms/v1_baseline.yml")
        print("2. Create and register your V2 algorithm")
        print("3. Start using algorithm tracking in your domain partition scripts")
        return 0
    else:
        print(f"\n⚠️  {total_tests - success_count} tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
