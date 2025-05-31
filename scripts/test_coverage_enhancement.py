#!/usr/bin/env python3
"""
Test Coverage Enhancement Script

This script helps improve test coverage for the algorithm versioning system
by identifying untested functionality and creating comprehensive test cases.
"""

import ast
import os
import sys
from pathlib import Path
from typing import List, Dict, Set, Any

def analyze_module_coverage(module_path: str) -> Dict[str, Any]:
    """Analyze a Python module to identify classes, methods, and functions"""
    
    with open(module_path, 'r') as f:
        content = f.read()
    
    tree = ast.parse(content)
    
    analysis = {
        'classes': {},
        'functions': [],
        'imports': [],
        'complexity_score': 0
    }
    
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            # Analyze class
            class_info = {
                'name': node.name,
                'methods': [],
                'properties': [],
                'line_start': node.lineno,
                'docstring': ast.get_docstring(node)
            }
            
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    method_info = {
                        'name': item.name,
                        'args': [arg.arg for arg in item.args.args],
                        'is_private': item.name.startswith('_'),
                        'is_property': any(isinstance(d, ast.Name) and d.id == 'property' 
                                         for d in item.decorator_list),
                        'line_start': item.lineno,
                        'docstring': ast.get_docstring(item)
                    }
                    
                    if method_info['is_property']:
                        class_info['properties'].append(method_info)
                    else:
                        class_info['methods'].append(method_info)
            
            analysis['classes'][node.name] = class_info
            
        elif isinstance(node, ast.FunctionDef) and not any(isinstance(parent, ast.ClassDef) for parent in ast.walk(tree)):
            # Top-level function
            func_info = {
                'name': node.name,
                'args': [arg.arg for arg in node.args.args],
                'line_start': node.lineno,
                'docstring': ast.get_docstring(node)
            }
            analysis['functions'].append(func_info)
            
        elif isinstance(node, ast.Import):
            for alias in node.names:
                analysis['imports'].append(alias.name)
                
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                for alias in node.names:
                    analysis['imports'].append(f"{node.module}.{alias.name}")
    
    # Calculate complexity score
    analysis['complexity_score'] = (
        len(analysis['classes']) * 3 +
        sum(len(cls['methods']) for cls in analysis['classes'].values()) * 2 +
        len(analysis['functions']) * 1
    )
    
    return analysis

def identify_missing_tests(module_analysis: Dict[str, Any], test_file_path: str = None) -> Dict[str, List[str]]:
    """Identify functions/methods that don't have corresponding tests"""
    
    tested_items = set()
    
    # If test file exists, analyze it to see what's already tested
    if test_file_path and Path(test_file_path).exists():
        test_analysis = analyze_module_coverage(test_file_path)
        
        for func in test_analysis['functions']:
            if func['name'].startswith('test_'):
                # Extract what this test is testing
                tested_name = func['name'].replace('test_', '').replace('_', '')
                tested_items.add(tested_name.lower())
    
    missing_tests = {
        'classes': [],
        'methods': [],
        'functions': []
    }
    
    # Check classes
    for class_name, class_info in module_analysis['classes'].items():
        test_name = f"test_{class_name.lower()}"
        if test_name.replace('_', '') not in tested_items:
            missing_tests['classes'].append(class_name)
        
        # Check methods
        for method in class_info['methods']:
            if not method['is_private']:  # Don't test private methods by default
                method_test_name = f"test_{class_name.lower()}_{method['name']}"
                if method_test_name.replace('_', '') not in tested_items:
                    missing_tests['methods'].append(f"{class_name}.{method['name']}")
    
    # Check functions
    for func in module_analysis['functions']:
        test_name = f"test_{func['name']}"
        if test_name.replace('_', '') not in tested_items:
            missing_tests['functions'].append(func['name'])
    
    return missing_tests

def generate_test_template(module_analysis: Dict[str, Any], missing_tests: Dict[str, List[str]]) -> str:
    """Generate pytest test template for missing test cases"""
    
    template = '''#!/usr/bin/env python3
"""
Enhanced Test Suite for Algorithm Versioning
Auto-generated test template - customize as needed
"""

import pytest
import tempfile
import json
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
        return context
    
    @pytest.fixture
    def mock_db_response(self):
        """Mock database responses"""
        return {
            'empty_list': [],
            'sample_algorithm': {
                'id': 1,
                'version_id': 'test_v1.0',
                'name': 'Test Algorithm',
                'description': 'Test description',
                'status': 'development',
                'config_data': '{"test": true}',
                'created_at': '2024-01-01T00:00:00',
                'created_by': 'test_user',
                'notes': 'Test notes'
            }
        }

'''
    
    # Generate class tests
    for class_name in missing_tests['classes']:
        class_info = module_analysis['classes'][class_name]
        
        template += f'''
    def test_{class_name.lower()}_creation(self, mock_context):
        """Test {class_name} creation and basic properties"""
        # TODO: Implement test for {class_name} creation
        instance = {class_name}(mock_context)
        assert instance is not None
        # Add more specific assertions

    def test_{class_name.lower()}_initialization(self, mock_context):
        """Test {class_name} initialization with different parameters"""
        # TODO: Test different initialization scenarios
        pass

'''
        
        # Generate method tests
        for method in class_info['methods']:
            if not method['is_private']:
                method_name = method['name']
                template += f'''
    def test_{class_name.lower()}_{method_name}(self, mock_context):
        """Test {class_name}.{method_name} method"""
        instance = {class_name}(mock_context)
        
        # TODO: Implement test for {method_name}
        # Mock dependencies as needed
        # Test various input scenarios
        # Assert expected outcomes
        pass

    def test_{class_name.lower()}_{method_name}_error_handling(self, mock_context):
        """Test {class_name}.{method_name} error handling"""
        instance = {class_name}(mock_context)
        
        # TODO: Test error cases for {method_name}
        # Test invalid inputs
        # Test exception handling
        pass

'''
    
    # Generate function tests
    for func_name in missing_tests['functions']:
        template += f'''
    def test_{func_name}(self):
        """Test {func_name} function"""
        # TODO: Implement test for {func_name}
        pass

    def test_{func_name}_edge_cases(self):
        """Test {func_name} edge cases"""
        # TODO: Test edge cases for {func_name}
        pass

'''
    
    # Add integration tests
    template += '''
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
    
    def test_full_algorithm_lifecycle(self, temp_config_file):
        """Test complete algorithm lifecycle: create -> register -> promote -> export"""
        # This is a key integration test
        # TODO: Implement full lifecycle test
        pass
    
    def test_algorithm_inheritance(self):
        """Test algorithm version inheritance (parent-child relationships)"""
        # TODO: Test parent-child algorithm relationships
        pass
    
    def test_concurrent_algorithm_operations(self):
        """Test concurrent access to algorithm versions"""
        # TODO: Test thread safety and concurrent operations
        pass

class TestAlgorithmVersioningValidation:
    """Tests for validation and error handling"""
    
    def test_invalid_algorithm_config(self):
        """Test handling of invalid algorithm configurations"""
        # TODO: Test various invalid configurations
        pass
    
    def test_algorithm_version_conflicts(self):
        """Test handling of version ID conflicts"""
        # TODO: Test duplicate version IDs
        pass
    
    def test_database_connection_errors(self):
        """Test handling of database connection issues"""
        # TODO: Test database error scenarios
        pass

class TestAlgorithmVersioningPerformance:
    """Performance tests for algorithm versioning"""
    
    def test_large_algorithm_list_performance(self):
        """Test performance with large number of algorithm versions"""
        # TODO: Test performance with many algorithms
        pass
    
    def test_complex_algorithm_config_performance(self):
        """Test performance with complex algorithm configurations"""
        # TODO: Test performance with large config files
        pass

class TestAlgorithmVersioningEdgeCases:
    """Tests for edge cases and boundary conditions"""
    
    def test_empty_algorithm_database(self):
        """Test behavior with empty algorithm database"""
        # TODO: Test empty database scenarios
        pass
    
    def test_corrupted_algorithm_config(self):
        """Test handling of corrupted configuration data"""
        # TODO: Test corrupted config handling
        pass
    
    def test_version_id_edge_cases(self):
        """Test edge cases for version ID formats"""
        # TODO: Test various version ID formats
        pass

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
        'notes': 'Test notes'
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
        elif test_group == "performance":
            pytest.main([__file__ + "::TestAlgorithmVersioningPerformance", "-v"])
        elif test_group == "edge_cases":
            pytest.main([__file__ + "::TestAlgorithmVersioningEdgeCases", "-v"])
        else:
            print(f"Unknown test group: {test_group}")
            sys.exit(1)
    else:
        # Run all tests
        pytest.main([__file__, "-v"])
'''
    
    return template

def create_coverage_report() -> str:
    """Create a comprehensive coverage analysis report"""
    
    modules_to_analyze = [
        "ecod/evaluation/algorithm_versions/manager.py",
        "ecod/cli/evaluation.py",
        "scripts/algorithm_manager.py"
    ]
    
    report = "# Algorithm Versioning Test Coverage Report\n\n"
    
    total_missing = 0
    
    for module_path in modules_to_analyze:
        if not Path(module_path).exists():
            report += f"## {module_path} - NOT FOUND\n\n"
            continue
            
        analysis = analyze_module_coverage(module_path)
        
        # Look for corresponding test file
        test_file = module_path.replace('.py', '_test.py').replace('ecod/', 'ecod/tests/unit/')
        if not Path(test_file).exists():
            test_file = f"ecod/tests/unit/test_{Path(module_path).stem}.py"
        
        missing = identify_missing_tests(analysis, test_file)
        
        module_missing = sum(len(items) for items in missing.values())
        total_missing += module_missing
        
        report += f"## {module_path}\n\n"
        report += f"**Complexity Score:** {analysis['complexity_score']}\n\n"
        report += f"**Classes:** {len(analysis['classes'])}\n"
        report += f"**Functions:** {len(analysis['functions'])}\n"
        report += f"**Missing Tests:** {module_missing}\n\n"
        
        if missing['classes']:
            report += "### Missing Class Tests\n"
            for cls in missing['classes']:
                report += f"- [ ] {cls}\n"
            report += "\n"
        
        if missing['methods']:
            report += "### Missing Method Tests\n"
            for method in missing['methods']:
                report += f"- [ ] {method}\n"
            report += "\n"
        
        if missing['functions']:
            report += "### Missing Function Tests\n"
            for func in missing['functions']:
                report += f"- [ ] {func}\n"
            report += "\n"
        
        # Class details
        if analysis['classes']:
            report += "### Class Details\n"
            for class_name, class_info in analysis['classes'].items():
                report += f"#### {class_name}\n"
                report += f"- Methods: {len(class_info['methods'])}\n"
                report += f"- Properties: {len(class_info['properties'])}\n"
                if class_info['docstring']:
                    report += f"- Documented: ✅\n"
                else:
                    report += f"- Documented: ❌\n"
                report += "\n"
    
    report += f"## Summary\n\n"
    report += f"**Total Missing Tests:** {total_missing}\n\n"
    
    if total_missing > 0:
        report += "### Recommendations\n\n"
        report += "1. **High Priority:** Implement tests for core algorithm operations\n"
        report += "2. **Medium Priority:** Add integration tests for CLI commands\n"
        report += "3. **Low Priority:** Add edge case and performance tests\n\n"
        report += "Use the generated test template to implement missing tests.\n"
    else:
        report += "✅ **Excellent!** All major functionality appears to be tested.\n"
    
    return report

def main():
    """Main entry point for test coverage analysis"""
    print("Algorithm Versioning Test Coverage Analysis")
    print("=" * 50)
    
    # Analyze main manager module
    manager_path = "ecod/evaluation/algorithm_versions/manager.py"
    
    if not Path(manager_path).exists():
        print(f"❌ Manager module not found: {manager_path}")
        print("Run the refactoring script first!")
        return 1
    
    print(f"📊 Analyzing {manager_path}...")
    analysis = analyze_module_coverage(manager_path)
    
    print(f"Found:")
    print(f"  - {len(analysis['classes'])} classes")
    print(f"  - {len(analysis['functions'])} functions")
    print(f"  - Complexity score: {analysis['complexity_score']}")
    
    # Check for missing tests
    test_file = "ecod/tests/unit/tests_basic_algorithm_versioning.py"
    missing = identify_missing_tests(analysis, test_file)
    
    total_missing = sum(len(items) for items in missing.values())
    
    print(f"\n🔍 Test Coverage Analysis:")
    print(f"  - Missing class tests: {len(missing['classes'])}")
    print(f"  - Missing method tests: {len(missing['methods'])}")
    print(f"  - Missing function tests: {len(missing['functions'])}")
    print(f"  - Total missing tests: {total_missing}")
    
    if total_missing > 0:
        print(f"\n📝 Generating enhanced test template...")
        
        # Generate comprehensive test template
        test_template = generate_test_template(analysis, missing)
        
        output_file = "ecod/tests/unit/test_algorithm_versioning_comprehensive.py"
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        with open(output_file, 'w') as f:
            f.write(test_template)
        
        print(f"✅ Created comprehensive test template: {output_file}")
        
        # Generate coverage report
        print(f"\n📊 Generating coverage report...")
        
        report = create_coverage_report()
        report_file = "test_coverage_report.md"
        
        with open(report_file, 'w') as f:
            f.write(report)
        
        print(f"✅ Created coverage report: {report_file}")
        
        print(f"\n🎯 Next Steps:")
        print(f"1. Review the comprehensive test template")
        print(f"2. Implement the high-priority missing tests")
        print(f"3. Run: pytest {output_file} -v")
        print(f"4. Check the coverage report for detailed analysis")
        
    else:
        print(f"\n🎉 Great! Test coverage appears to be comprehensive.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
