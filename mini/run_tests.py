#!/usr/bin/env python3
"""
Test runner for mini_pyecod with coverage reporting

This script provides convenient commands for running different test suites
and generating coverage reports.

Test Structure:
- Unit tests: Fast tests with mock data (test_core.py, test_models.py, etc.)
- Integration tests: Tests with real data (test_cases.py with real files)
- Primary test: The canonical 8ovp_A test case
- Performance tests: Benchmark tests (marked with @pytest.mark.performance)

Usage:
    python run_tests.py --help
    python run_tests.py unit              # Fast unit tests only
    python run_tests.py integration       # Integration tests with real data
    python run_tests.py primary           # Just the primary 8ovp_A test
    python run_tests.py all               # All tests
    python run_tests.py coverage          # All tests with coverage report
    python run_tests.py specific tests/test_parser.py  # Run specific test file
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
from typing import List, Optional


class TestRunner:
    """Test runner with coverage and reporting"""

    def __init__(self):
        self.mini_dir = Path(__file__).parent
        self.tests_dir = self.mini_dir / "tests"
        self.coverage_dir = self.mini_dir / "htmlcov"

        # Test file categories
        self.unit_test_files = [
            "test_core.py",          # Core algorithm tests
            "test_models.py",        # Data model tests
            "test_parser.py",        # XML parser tests
            "test_writer.py",        # Output writer tests
            "test_decomposer.py",    # Decomposer tests
            "test_blast_parser.py",  # BLAST parser tests
            "test_discontinuous.py", # Discontinuous domain tests
            "test_range_cache_parser.py",  # Cache parser tests
            "test_cli.py",          # CLI tests
        ]

        self.integration_test_files = [
            "test_cases.py",         # Full pipeline tests
            "test_ecod_regression.py",  # Regression tests
            "test_ecod_tgroup.py",   # T-group validation
        ]

    def run_command(self, cmd: List[str], description: str = None) -> int:
        """Run a command and return exit code"""
        if description:
            print(f"\n{'='*70}")
            print(f"🔧 {description}")
            print(f"{'='*70}")

        print(f"📍 Running: {' '.join(cmd)}")
        print(f"📂 Working directory: {self.mini_dir}")
        print("-" * 70)

        try:
            result = subprocess.run(cmd, cwd=self.mini_dir, check=False)

            if result.returncode == 0:
                print(f"\n✅ {description or 'Command'} completed successfully")
            else:
                print(f"\n❌ {description or 'Command'} failed with exit code {result.returncode}")

            return result.returncode

        except KeyboardInterrupt:
            print("\n\n⚠️  Interrupted by user")
            return 130
        except Exception as e:
            print(f"\n❌ Error running command: {e}")
            return 1

    def check_dependencies(self) -> bool:
        """Check that required test dependencies are installed"""
        required_packages = ['pytest', 'pytest-cov']
        missing = []

        for package in required_packages:
            try:
                __import__(package.replace('-', '_'))
            except ImportError:
                missing.append(package)

        if missing:
            print("❌ Missing required packages:")
            for pkg in missing:
                print(f"   pip install {pkg}")
            print("\nInstall all test dependencies:")
            print("   pip install pytest pytest-cov")
            return False

        return True

    def validate_environment(self) -> bool:
        """Validate test environment"""
        print("\n🔍 Validating test environment...")

        # Check test data
        test_data_dir = self.mini_dir / "test_data"
        required_files = [
            "domain_lengths.csv",
            "protein_lengths.csv",
            "domain_definitions.csv"
        ]

        missing_files = []
        for file in required_files:
            if not (test_data_dir / file).exists():
                missing_files.append(file)

        if missing_files:
            print("\n⚠️  Missing test data files:")
            for file in missing_files:
                print(f"   {test_data_dir / file}")
            print("\nTo generate test data:")
            print("   cd mini")
            print("   python range_cache_parser.py --output-dir test_data")
            return False

        # Check ECOD batch directory (optional for unit tests)
        batch_dir = Path("/data/ecod/pdb_updates/batches")
        if not batch_dir.exists():
            print("\n⚠️  ECOD batch directory not found (integration tests will be skipped)")
            print(f"   Expected: {batch_dir}")
        else:
            print(f"✅ ECOD batch directory found: {batch_dir}")

        # Check for ecod.core module
        try:
            import ecod.core.sequence_range
            print("✅ ecod.core module available")
        except ImportError:
            print("❌ ecod.core module not found - ensure PYTHONPATH is set correctly")
            return False

        print("\n✅ Test environment validated")
        return True

    def run_unit_tests(self, verbose: bool = False) -> int:
        """Run fast unit tests with mock data"""
        cmd = [
            "python", "-m", "pytest",
            "tests/",
            "-m", "unit",
            "--tb=short"
        ]

        if verbose:
            cmd.append("-v")

        # Add coverage for unit tests
        cmd.extend(["--cov=mini", "--cov-report=term-missing:skip-covered"])

        return self.run_command(cmd, "Running Unit Tests")

    def run_integration_tests(self, verbose: bool = False) -> int:
        """Run integration tests with real data"""
        cmd = [
            "python", "-m", "pytest",
            "tests/",
            "-m", "integration",
            "--tb=short"
        ]

        if verbose:
            cmd.append("-v")

        return self.run_command(cmd, "Running Integration Tests")

    def run_primary_test(self, verbose: bool = False) -> int:
        """Run just the primary 8ovp_A test case"""
        cmd = [
            "python", "-m", "pytest",
            "tests/test_cases.py::TestOfficialCases::test_8ovp_A_canonical",
            "--tb=short"
        ]

        if verbose:
            cmd.extend(["-v", "-s"])  # -s to see print statements

        return self.run_command(cmd, "Running Primary Test Case (8ovp_A)")

    def run_all_tests(self, verbose: bool = False, skip_slow: bool = False) -> int:
        """Run all tests"""
        cmd = [
            "python", "-m", "pytest",
            "tests/",
            "--tb=short"
        ]

        if skip_slow:
            cmd.extend(["-m", "not slow"])

        if verbose:
            cmd.append("-v")

        return self.run_command(cmd, "Running All Tests")

    def run_with_coverage(self, verbose: bool = False,
                         html_report: bool = True,
                         term_report: bool = True) -> int:
        """Run tests with coverage reporting"""
        cmd = [
            "python", "-m", "pytest",
            "tests/",
            "--cov=mini",
            "--cov-branch",
            "--tb=short"
        ]

        if html_report:
            cmd.append("--cov-report=html")

        if term_report:
            cmd.append("--cov-report=term-missing")

        if verbose:
            cmd.append("-v")

        result = self.run_command(cmd, "Running Tests with Coverage")

        if result == 0 and html_report:
            coverage_index = self.coverage_dir / "index.html"
            if coverage_index.exists():
                print(f"\n📊 Coverage report generated: {coverage_index}")
                print(f"   Open in browser: file://{coverage_index.absolute()}")

        return result

    def run_performance_tests(self, verbose: bool = False) -> int:
        """Run performance/benchmark tests"""
        cmd = [
            "python", "-m", "pytest",
            "tests/",
            "-m", "performance",
            "--tb=short"
        ]

        if verbose:
            cmd.extend(["-v", "--durations=10"])  # Show slowest tests

        return self.run_command(cmd, "Running Performance Tests")

    def clean_coverage(self):
        """Clean coverage artifacts"""
        import shutil

        print("\n🧹 Cleaning test artifacts...")

        artifacts = [
            self.coverage_dir,
            self.mini_dir / ".coverage",
            self.mini_dir / ".pytest_cache"
        ]

        for artifact in artifacts:
            if artifact.exists():
                if artifact.is_dir():
                    shutil.rmtree(artifact)
                else:
                    artifact.unlink()
                print(f"   Cleaned: {artifact.name}")

        print("✅ Test artifacts cleaned")

    def list_tests(self) -> int:
        """List available tests"""
        cmd = [
            "python", "-m", "pytest",
            "tests/",
            "--collect-only",
            "-q"
        ]

        return self.run_command(cmd, "Available Tests")

    def test_specific(self, test_pattern: str, verbose: bool = False) -> int:
        """Run specific test pattern"""
        cmd = [
            "python", "-m", "pytest",
            test_pattern,
            "--tb=short"
        ]

        if verbose:
            cmd.extend(["-v", "-s"])

        return self.run_command(cmd, f"Running Tests: {test_pattern}")

    def show_test_summary(self):
        """Show summary of test organization"""
        print("\n📋 TEST SUITE ORGANIZATION")
        print("=" * 70)

        print("\n📦 Unit Test Files (fast, mock data):")
        for test_file in self.unit_test_files:
            test_path = self.tests_dir / test_file
            if test_path.exists():
                print(f"   ✅ {test_file}")
            else:
                print(f"   ❌ {test_file} (missing)")

        print("\n🔗 Integration Test Files (slower, real data):")
        for test_file in self.integration_test_files:
            test_path = self.tests_dir / test_file
            if test_path.exists():
                print(f"   ✅ {test_file}")
            else:
                print(f"   ❌ {test_file} (missing)")

        print("\n🏷️  Test Markers:")
        print("   @pytest.mark.unit         - Fast unit tests")
        print("   @pytest.mark.integration  - Integration tests")
        print("   @pytest.mark.slow         - Slow tests")
        print("   @pytest.mark.performance  - Performance benchmarks")
        print("   @pytest.mark.visualization - Visualization tests")

        print("\n🎯 Primary Test Case:")
        print("   8ovp_A - GFP-PBP fusion with chain BLAST decomposition")
        print("   Location: tests/test_cases.py::TestOfficialCases::test_8ovp_A_canonical")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Mini PyECOD Test Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Test Categories:
  unit          Fast unit tests with mock data (~5 seconds)
  integration   Integration tests with real data (~30 seconds)
  primary       Just the primary 8ovp_A test case (~10 seconds)
  performance   Performance/benchmark tests
  all           All tests (~60 seconds)
  coverage      All tests with coverage report

Examples:
  python run_tests.py unit -v               # Verbose unit tests
  python run_tests.py integration           # Integration tests
  python run_tests.py primary               # Quick validation with 8ovp_A
  python run_tests.py coverage              # Full coverage report
  python run_tests.py specific tests/test_parser.py  # Run specific file
  python run_tests.py --list                # Show available tests
  python run_tests.py --summary             # Show test organization
        """
    )

    # Test categories
    parser.add_argument('category', nargs='?',
                       choices=['unit', 'integration', 'primary', 'performance', 'all', 'coverage', 'specific'],
                       help='Test category to run')

    # Options
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose test output')
    parser.add_argument('--skip-slow', action='store_true',
                       help='Skip slow tests (for all/coverage)')
    parser.add_argument('--no-html', action='store_true',
                       help='Skip HTML coverage report')
    parser.add_argument('--clean', action='store_true',
                       help='Clean coverage artifacts and exit')
    parser.add_argument('--list', action='store_true',
                       help='List available tests and exit')
    parser.add_argument('--summary', action='store_true',
                       help='Show test suite organization')
    parser.add_argument('--check-deps', action='store_true',
                       help='Check dependencies and exit')
    parser.add_argument('test_path', nargs='?',
                       help='Specific test path for "specific" category')

    args = parser.parse_args()

    runner = TestRunner()

    # Handle special commands
    if args.clean:
        runner.clean_coverage()
        return 0

    if args.check_deps:
        if runner.check_dependencies():
            print("✅ All dependencies available")
            return 0
        else:
            return 1

    if args.summary:
        runner.show_test_summary()
        return 0

    if args.list:
        return runner.list_tests()

    # Check dependencies before running tests
    if not runner.check_dependencies():
        return 1

    # Validate environment (warn but don't fail)
    runner.validate_environment()

    # Handle specific test pattern
    if args.category == 'specific':
        if not args.test_path:
            print("❌ Please specify a test path for 'specific' category")
            print("   Example: python run_tests.py specific tests/test_parser.py")
            return 1
        return runner.test_specific(args.test_path, args.verbose)

    # Run requested test category
    if not args.category:
        parser.print_help()
        print("\n❌ Please specify a test category")
        print("   Try: python run_tests.py unit")
        return 1

    print(f"\n🚀 Mini PyECOD Test Runner")
    print(f"📂 Working directory: {runner.mini_dir}")
    print(f"🎯 Test category: {args.category}")

    if args.category == 'unit':
        result = runner.run_unit_tests(args.verbose)
    elif args.category == 'integration':
        result = runner.run_integration_tests(args.verbose)
    elif args.category == 'primary':
        result = runner.run_primary_test(args.verbose)
    elif args.category == 'performance':
        result = runner.run_performance_tests(args.verbose)
    elif args.category == 'all':
        result = runner.run_all_tests(args.verbose, args.skip_slow)
    elif args.category == 'coverage':
        result = runner.run_with_coverage(
            args.verbose,
            html_report=not args.no_html
        )
    else:
        print(f"❌ Unknown category: {args.category}")
        return 1

    # Summary
    print("\n" + "=" * 70)
    if result == 0:
        print("✅ TESTS PASSED")
    else:
        print(f"❌ TESTS FAILED (exit code: {result})")

    return result


if __name__ == "__main__":
    sys.exit(main())
