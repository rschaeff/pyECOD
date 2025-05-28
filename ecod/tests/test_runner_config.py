#!/usr/bin/env python3
"""
Test runner configuration for partition module tests.

Provides organized test execution, coverage reporting, and integration
with existing test infrastructure.
"""

import pytest
import sys
import os
from pathlib import Path


class PartitionTestRunner:
    """Organized test runner for partition module tests"""
    
    def __init__(self, test_dir=None):
        """Initialize test runner
        
        Args:
            test_dir: Directory containing test files (defaults to current directory)
        """
        self.test_dir = Path(test_dir) if test_dir else Path(__file__).parent
        
        # Test suite organization
        self.test_suites = {
            'evidence_analyzer': 'test_evidence_analyzer.py',
            'evidence_confidence': 'test_evidence_confidence.py', 
            'batch_processing': 'test_batch_processing.py',
            'error_handling': 'test_error_handling.py',
            'existing_models': 'test_partition_models.py',  # Your existing tests
            'existing_service': 'test_partition_service.py'  # Your existing tests
        }
        
        # Test categories for targeted execution
        self.categories = {
            'high_priority': [
                'evidence_analyzer',
                'evidence_confidence', 
                'batch_processing',
                'error_handling'
            ],
            'core_models': [
                'existing_models',
                'evidence_confidence'
            ],
            'integration': [
                'batch_processing',
                'existing_service'
            ],
            'reliability': [
                'error_handling',
                'batch_processing'
            ],
            'all': list(self.test_suites.keys())
        }
    
    def run_category(self, category='high_priority', coverage=True, verbose=True):
        """Run tests for a specific category
        
        Args:
            category: Test category to run ('high_priority', 'core_models', etc.)
            coverage: Whether to generate coverage report
            verbose: Whether to run in verbose mode
            
        Returns:
            pytest exit code
        """
        if category not in self.categories:
            print(f"Unknown category: {category}")
            print(f"Available categories: {list(self.categories.keys())}")
            return 1
        
        test_files = []
        for suite_name in self.categories[category]:
            if suite_name in self.test_suites:
                test_file = self.test_dir / self.test_suites[suite_name]
                if test_file.exists():
                    test_files.append(str(test_file))
                else:
                    print(f"Warning: Test file not found: {test_file}")
        
        if not test_files:
            print(f"No test files found for category: {category}")
            return 1
        
        # Build pytest arguments
        pytest_args = []
        
        if verbose:
            pytest_args.extend(['-v', '--tb=short'])
        
        if coverage:
            pytest_args.extend([
                '--cov=ecod.pipelines.domain_analysis.partition',
                '--cov-report=html:htmlcov',
                '--cov-report=term-missing',
                '--cov-fail-under=80'
            ])
        
        # Add specific markers and filters
        pytest_args.extend([
            '--strict-markers',
            '--disable-warnings',  # Reduce noise during development
            '-x'  # Stop on first failure for faster iteration
        ])
        
        pytest_args.extend(test_files)
        
        print(f"\nRunning {category} tests:")
        print(f"Test files: {[Path(f).name for f in test_files]}")
        print(f"Pytest args: {' '.join(pytest_args)}")
        print("-" * 60)
        
        return pytest.main(pytest_args)
    
    def run_single_suite(self, suite_name, **kwargs):
        """Run a single test suite
        
        Args:
            suite_name: Name of test suite to run
            **kwargs: Additional arguments passed to pytest
        """
        if suite_name not in self.test_suites:
            print(f"Unknown test suite: {suite_name}")
            print(f"Available suites: {list(self.test_suites.keys())}")
            return 1
        
        test_file = self.test_dir / self.test_suites[suite_name]
        if not test_file.exists():
            print(f"Test file not found: {test_file}")
            return 1
        
        pytest_args = ['-v', str(test_file)]
        
        # Add any additional pytest arguments
        for key, value in kwargs.items():
            if key.startswith('pytest_'):
                arg_name = key.replace('pytest_', '--')
                if value is True:
                    pytest_args.append(arg_name)
                elif value is not False:
                    pytest_args.extend([arg_name, str(value)])
        
        return pytest.main(pytest_args)
    
    def run_focused_tests(self, pattern, **kwargs):
        """Run tests matching a specific pattern
        
        Args:
            pattern: Test name pattern to match
            **kwargs: Additional pytest arguments
        """
        all_test_files = [
            str(self.test_dir / filename) 
            for filename in self.test_suites.values()
            if (self.test_dir / filename).exists()
        ]
        
        pytest_args = ['-v', '-k', pattern] + all_test_files
        
        return pytest.main(pytest_args)
    
    def run_performance_tests(self):
        """Run performance-focused tests"""
        return self.run_focused_tests('performance or large or concurrent', 
                                    pytest_durations=10)
    
    def run_quick_tests(self):
        """Run quick smoke tests for fast feedback"""
        pytest_args = [
            '-v', '-x',  # Stop on first failure
            '--tb=line',  # Minimal traceback
            '-k', 'not (large or performance or concurrent)',  # Skip slow tests
        ]
        
        # Add core test files
        core_files = [
            'test_evidence_analyzer.py',
            'test_evidence_confidence.py'
        ]
        
        for filename in core_files:
            test_file = self.test_dir / filename
            if test_file.exists():
                pytest_args.append(str(test_file))
        
        return pytest.main(pytest_args)
    
    def generate_test_report(self, output_file='test_report.html'):
        """Generate comprehensive test report"""
        pytest_args = [
            '--html=' + output_file,
            '--self-contained-html',
            '--cov=ecod.pipelines.domain_analysis.partition',
            '--cov-report=html:htmlcov',
            '--junitxml=test_results.xml'
        ]
        
        # Add all test files
        for filename in self.test_suites.values():
            test_file = self.test_dir / filename
            if test_file.exists():
                pytest_args.append(str(test_file))
        
        return pytest.main(pytest_args)


def main():
    """Command line interface for test runner"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Run partition module tests')
    parser.add_argument('--category', '-c', 
                       choices=['high_priority', 'core_models', 'integration', 'reliability', 'all'],
                       default='high_priority',
                       help='Test category to run')
    parser.add_argument('--suite', '-s',
                       help='Run specific test suite')
    parser.add_argument('--pattern', '-k',
                       help='Run tests matching pattern')
    parser.add_argument('--no-coverage', action='store_true',
                       help='Disable coverage reporting')
    parser.add_argument('--quick', action='store_true',
                       help='Run quick smoke tests only')
    parser.add_argument('--performance', action='store_true',
                       help='Run performance tests only')
    parser.add_argument('--report', action='store_true',
                       help='Generate comprehensive test report')
    parser.add_argument('--test-dir',
                       help='Directory containing test files')
    
    args = parser.parse_args()
    
    runner = PartitionTestRunner(args.test_dir)
    
    if args.quick:
        return runner.run_quick_tests()
    elif args.performance:
        return runner.run_performance_tests()
    elif args.report:
        return runner.generate_test_report()
    elif args.suite:
        return runner.run_single_suite(args.suite)
    elif args.pattern:
        return runner.run_focused_tests(args.pattern)
    else:
        return runner.run_category(args.category, 
                                 coverage=not args.no_coverage)


if __name__ == '__main__':
    sys.exit(main())
