#!/usr/bin/env python3
"""
Production regression test runner for mini_pyecod

This script runs comprehensive regression tests to validate that mini
can replace the current pyecod partitioning engine.

Usage:
    python run_production_tests.py --primary           # Just 8ovp_A
    python run_production_tests.py --suite             # All test proteins  
    python run_production_tests.py --validate-setup    # Check test environment
    python run_production_tests.py --baseline          # Establish performance baseline
"""

import sys
import subprocess
from pathlib import Path

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Mini PyECOD Production Testing')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--primary', action='store_true', 
                      help='Run primary test case (8ovp_A)')
    group.add_argument('--suite', action='store_true',
                      help='Run full test suite')
    group.add_argument('--validate-setup', action='store_true',
                      help='Validate test environment setup')
    group.add_argument('--baseline', action='store_true',
                      help='Establish performance baseline')
    
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    if args.validate_setup:
        # Validate test environment
        cmd = ["python", "-c", "from setup_production_testing import ProductionTestSetup; setup = ProductionTestSetup(); setup.validate_environment()"]
        return subprocess.call(cmd)
    
    elif args.primary:
        # Run primary regression test
        cmd = ["python", "-m", "pytest", 
               "tests/test_ecod_regression.py::TestRegressionSuite::test_mini_vs_current_primary",
                  "-v" if args.verbose else ""]
        return subprocess.call([c for c in cmd if c])
      
    elif args.suite:
        # Run full regression suite
        cmd = ["python", "-m", "pytest",
                 "tests/test_ecod_regression.py::TestRegressionSuite::test_mini_vs_current_suite", 
                "-v" if args.verbose else ""]
        return subprocess.call([c for c in cmd if c])
    
    elif args.baseline:
        # Run baseline establishment
        cmd = ["python", "tests/test_ecod_regression.py", "--proteins"] + ['8ovp_A']
        return subprocess.call(cmd)

if __name__ == "__main__":
    exit(main())
