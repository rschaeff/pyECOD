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
            'basic_imports': 'test_basic_imports.py',
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

    def _check_dependencies(self):
        """Check for optional test dependencies"""
        missing_deps = []
        available_features = {}

        # Check for pytest-cov
        try:
            import pytest_cov
            available_features['coverage'] = True
        except ImportError:
            available_features['coverage'] = False
            missing_deps.append('pytest-cov')

        # Check for pytest-html
        try:
            import pytest_html
            available_features['html_reports'] = True
        except ImportError:
            available_features['html_reports'] = False
            missing_deps.append('pytest-html')

        return available_features, missing_deps

    def run_category(self, category='high_priority', coverage=False, verbose=True):
        """Run tests for a specific category

        Args:
            category: Test category to run ('high_priority', 'core_models', etc.)
            coverage: Whether to generate coverage report (requires pytest-cov)
            verbose: Whether to run in verbose mode

        Returns:
            pytest exit code
        """
        if category not in self.categories:
            print(f"Unknown category: {category}")
            print(f"Available categories: {list(self.categories.keys())}")
            return 1

        # Check dependencies
        available_features, missing_deps = self._check_dependencies()

        if coverage and not available_features['coverage']:
            print("⚠️  Coverage requested but pytest-cov not installed.")
            print("   Install with: pip install pytest-cov")
            print("   Running without coverage...\n")
            coverage = False

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

        if coverage and available_features['coverage']:
            pytest_args.extend([
                '--cov=ecod.pipelines.domain_analysis.partition',
                '--cov-report=html:htmlcov',
                '--cov-report=term-missing',
                '--cov-fail-under=70'  # Reduced threshold for initial runs
            ])

        # Add specific markers and filters
        pytest_args.extend([
            '--tb=short',  # Short traceback format
            '-x'  # Stop on first failure for faster iteration
        ])

        # Only add strict-markers if not causing issues
        try:
            pytest_args.append('--strict-markers')
        except:
            pass

        pytest_args.extend(test_files)

        print(f"\nRunning {category} tests:")
        print(f"Test files: {[Path(f).name for f in test_files]}")
        if coverage and available_features['coverage']:
            print("📊 Coverage reporting enabled")
        print(f"Command: pytest {' '.join(pytest_args[:-len(test_files)])}")
        print("-" * 60)

        return pytest.main(pytest_args)

    def run_single_suite(self, suite_name, coverage=False, **kwargs):
        """Run a single test suite

        Args:
            suite_name: Name of test suite to run
            coverage: Whether to generate coverage
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

        available_features, _ = self._check_dependencies()

        pytest_args = ['-v', str(test_file)]

        if coverage and available_features['coverage']:
            pytest_args.extend([
                '--cov=ecod.pipelines.domain_analysis.partition',
                '--cov-report=term-missing'
            ])

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

    def run_basic_fallback(self):
        """Run basic tests using pure Python (no pytest required)"""
        basic_test_file = self.test_dir / 'test_basic_imports.py'

        if not basic_test_file.exists():
            print(f"❌ Basic test file not found: {basic_test_file}")
            return 1

        import subprocess
        import sys

        try:
            result = subprocess.run([sys.executable, str(basic_test_file)],
                                  capture_output=True, text=True, cwd=self.test_dir.parent)

            print(result.stdout)
            if result.stderr:
                print("STDERR:", result.stderr)

            return result.returncode
        except Exception as e:
            print(f"❌ Error running basic tests: {e}")
            return 1

    def run_quick_tests(self):
        """Run quick smoke tests for fast feedback"""
        pytest_args = [
            '-v', '-x',  # Stop on first failure
            '--tb=line',  # Minimal traceback
            '-k', 'not (large or performance or concurrent)',  # Skip slow tests
        ]

        # Add core test files (prioritize basic imports for quick feedback)
        core_files = [
            'test_basic_imports.py',
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
        available_features, missing_deps = self._check_dependencies()

        pytest_args = []

        if available_features['html_reports']:
            pytest_args.extend([
                '--html=' + output_file,
                '--self-contained-html'
            ])
        else:
            print("⚠️  HTML reports requested but pytest-html not installed.")
            print("   Install with: pip install pytest-html")

        if available_features['coverage']:
            pytest_args.extend([
                '--cov=ecod.pipelines.domain_analysis.partition',
                '--cov-report=html:htmlcov',
                '--junitxml=test_results.xml'
            ])
        else:
            print("⚠️  Coverage requested but pytest-cov not installed.")
            print("   Install with: pip install pytest-cov")

        # Add all test files
        for filename in self.test_suites.values():
            test_file = self.test_dir / filename
            if test_file.exists():
                pytest_args.append(str(test_file))

        if not pytest_args or all(arg.startswith('-') for arg in pytest_args[:-len([f for f in self.test_suites.values() if (self.test_dir / f).exists()])]):
            # No special reporting features available, run basic tests
            print("Running basic test execution...")
            pytest_args = ['-v'] + [str(self.test_dir / f) for f in self.test_suites.values() if (self.test_dir / f).exists()]

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
    parser.add_argument('--coverage', action='store_true',
                       help='Enable coverage reporting (requires pytest-cov)')
    parser.add_argument('--quick', action='store_true',
                       help='Run quick smoke tests only')
    parser.add_argument('--performance', action='store_true',
                       help='Run performance tests only')
    parser.add_argument('--report', action='store_true',
                       help='Generate comprehensive test report')
    parser.add_argument('--test-dir',
                       help='Directory containing test files')
    parser.add_argument('--basic', action='store_true',
                       help='Run basic import tests (no pytest required)')
    parser.add_argument('--check-deps', action='store_true',
                       help='Check test dependencies and exit')

    args = parser.parse_args()

    runner = PartitionTestRunner(args.test_dir)

    # Check dependencies if requested
    if args.check_deps:
        available_features, missing_deps = runner._check_dependencies()
        print("Test Dependencies Status:")
        print("=" * 40)

        if available_features['coverage']:
            print("✅ pytest-cov (coverage reporting)")
        else:
            print("❌ pytest-cov (coverage reporting)")

        if available_features['html_reports']:
            print("✅ pytest-html (HTML reports)")
        else:
            print("❌ pytest-html (HTML reports)")

        if missing_deps:
            print(f"\nMissing dependencies: {', '.join(missing_deps)}")
            print(f"Install with: pip install {' '.join(missing_deps)}")
            return 1
        else:
            print("\n✅ All optional dependencies available!")
            return 0

    if args.basic:
        return runner.run_basic_fallback()
    elif args.quick:
        return runner.run_quick_tests()
    elif args.performance:
        return runner.run_performance_tests()
    elif args.report:
        return runner.generate_test_report()
    elif args.suite:
        return runner.run_single_suite(args.suite, coverage=args.coverage)
    elif args.pattern:
        return runner.run_focused_tests(args.pattern)
    else:
        return runner.run_category(args.category, coverage=args.coverage)


if __name__ == '__main__':
    sys.exit(main())
