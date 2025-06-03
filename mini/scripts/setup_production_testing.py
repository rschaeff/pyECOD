#!/usr/bin/env python3
"""
Complete setup script for mini_pyecod production testing

This script sets up all required data files for comprehensive regression testing:
1. Range cache parsing (domain/protein lengths)
2. ECOD domains parsing (T-group classifications)
3. Test environment validation
4. Performance baseline establishment

Usage:
    python setup_production_testing.py --full-setup
    python setup_production_testing.py --test-proteins-only
    python setup_production_testing.py --validate
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Tuple

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

# Import our parsers
from mini.range_cache_parser import (
    create_domain_lengths_from_cache,
    create_domain_definitions_from_cache,
    extract_protein_lengths_from_cache,
    validate_cache_data
)
from mini.ecod_domains_parser import (
    create_ecod_classifications_file,
    validate_ecod_data_for_protein,
    get_ecod_summary,
    parse_ecod_domains_file
)

class ProductionTestSetup:
    """Setup coordinator for production testing environment"""
    
    def __init__(self, output_dir: str = "test_data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Default file paths
        self.range_cache_file = "/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt"
        self.domains_file = "/data/ecod/database_versions/v291/ecod.develop291.domains.txt"
        
        # Production test proteins (start with validated cases, expand to 10-15)
        self.production_test_proteins = [
            "8ovp_A",   # Primary: GFP-PBP fusion with decomposition
            # Add more as they are validated:
            # "1abc_A",   # Single domain
            # "2def_A",   # Two domain
            # "3ghi_B",   # Multi-domain
            # etc.
        ]
    
    def validate_environment(self) -> Dict[str, bool]:
        """Validate that required ECOD files are available"""
        
        validation = {
            'range_cache_exists': os.path.exists(self.range_cache_file),
            'domains_file_exists': os.path.exists(self.domains_file),
            'output_dir_writable': os.access(self.output_dir, os.W_OK),
            'ecod_batch_dir_exists': os.path.exists("/data/ecod/pdb_updates/batches")
        }
        
        print("=== ENVIRONMENT VALIDATION ===")
        for check, status in validation.items():
            icon = "✓" if status else "✗"
            print(f"{icon} {check.replace('_', ' ').title()}")
        
        if not validation['range_cache_exists']:
            print(f"   Range cache not found: {self.range_cache_file}")
        
        if not validation['domains_file_exists']:
            print(f"   Domains file not found: {self.domains_file}")
        
        return validation
    
    def setup_range_cache_data(self, verbose: bool = False) -> bool:
        """Setup reference data from range cache"""
        
        print("\n=== SETTING UP RANGE CACHE DATA ===")
        
        if not os.path.exists(self.range_cache_file):
            print(f"ERROR: Range cache file not found: {self.range_cache_file}")
            return False
        
        try:
            # Generate domain lengths
            domain_lengths_file = self.output_dir / "domain_lengths.csv"
            print(f"Generating domain lengths: {domain_lengths_file}")
            create_domain_lengths_from_cache(self.range_cache_file, str(domain_lengths_file), verbose)
            
            # Generate domain definitions  
            domain_definitions_file = self.output_dir / "domain_definitions.csv"
            print(f"Generating domain definitions: {domain_definitions_file}")
            create_domain_definitions_from_cache(self.range_cache_file, str(domain_definitions_file), verbose)
            
            # Generate protein lengths
            protein_lengths_file = self.output_dir / "protein_lengths.csv"
            print(f"Generating protein lengths: {protein_lengths_file}")
            extract_protein_lengths_from_cache(self.range_cache_file, str(protein_lengths_file), verbose)
            
            # Validate key test proteins
            print(f"Validating test proteins in range cache...")
            validate_cache_data(self.range_cache_file, ['e8ovpA1', 'e8ovpA2'])  # Expected 8ovp domains
            
            print("✅ Range cache data setup complete")
            return True
            
        except Exception as e:
            print(f"❌ Range cache setup failed: {e}")
            return False
    
    def setup_ecod_classifications(self, verbose: bool = False) -> bool:
        """Setup ECOD classifications for production testing"""
        
        print("\n=== SETTING UP ECOD CLASSIFICATIONS ===")
        
        if not os.path.exists(self.domains_file):
            print(f"ERROR: ECOD domains file not found: {self.domains_file}")
            return False
        
        try:
            # Create ECOD classifications for test proteins
            classifications_file = self.output_dir / "ecod_classifications.csv"
            print(f"Generating ECOD classifications: {classifications_file}")
            
            create_ecod_classifications_file(
                self.domains_file,
                str(classifications_file),
                self.production_test_proteins,
                verbose
            )
            
            # Validate that our test proteins have ECOD data
            print(f"Validating ECOD data for test proteins...")
            classifications = parse_ecod_domains_file(self.domains_file, verbose=False)
            
            missing_proteins = []
            for protein_id in self.production_test_proteins:
                has_data = validate_ecod_data_for_protein(protein_id, classifications)
                status = "✓" if has_data else "✗"
                print(f"  {status} {protein_id}")
                if not has_data:
                    missing_proteins.append(protein_id)
            
            if missing_proteins:
                print(f"⚠️  Warning: {len(missing_proteins)} test proteins missing ECOD data: {missing_proteins}")
                print("   This may indicate they are new structures not yet in ECOD")
            
            # Print summary
            summary = get_ecod_summary(classifications)
            print(f"\nECOD Data Summary:")
            print(f"  Total protein chains: {summary['total_chains']:,}")
            print(f"  Total domains: {summary['total_domains']:,}")
            print(f"  Unique T-groups: {summary['unique_t_groups']:,}")
            
            print("✅ ECOD classifications setup complete")
            return True
            
        except Exception as e:
            print(f"❌ ECOD classifications setup failed: {e}")
            return False
    
    def create_test_configuration(self) -> bool:
        """Create test configuration file"""
        
        print("\n=== CREATING TEST CONFIGURATION ===")
        
        config_content = f'''# Mini PyECOD Production Test Configuration
# Generated by setup_production_testing.py

# Test Data Files
DOMAIN_LENGTHS_FILE = "{self.output_dir}/domain_lengths.csv"
PROTEIN_LENGTHS_FILE = "{self.output_dir}/protein_lengths.csv"
DOMAIN_DEFINITIONS_FILE = "{self.output_dir}/domain_definitions.csv"
ECOD_CLASSIFICATIONS_FILE = "{self.output_dir}/ecod_classifications.csv"

# Production Test Proteins
PRODUCTION_TEST_PROTEINS = {self.production_test_proteins}

# ECOD Source Files (for reference)
RANGE_CACHE_FILE = "{self.range_cache_file}"
DOMAINS_FILE = "{self.domains_file}"

# Test Thresholds
MIN_ACCEPTABLE_QUALITY_SCORE = -0.2    # Allow up to 20% degradation from current
PREFERRED_QUALITY_SCORE = 0.1          # Prefer 10% improvement
MIN_ACCEPTABLE_RATE = 0.8               # 80% of tests must be acceptable

# Performance Benchmarks
MAX_PROCESSING_TIME_8OVP_A = 60.0       # seconds
TARGET_PERFORMANCE_IMPROVEMENT = 0.2    # 20% faster than current

# Coverage Requirements
MIN_COVERAGE_8OVP_A = 0.85              # 85% sequence coverage
MIN_DOMAIN_COUNT_8OVP_A = 2             # At least 2 domains
EXPECTED_DOMAIN_COUNT_8OVP_A = 3        # Preferably 3 domains

# ECOD Classification Requirements
MIN_ECOD_ACCURACY = 0.7                 # 70% classification accuracy
TARGET_ECOD_ACCURACY = 0.9              # 90% classification accuracy
'''
        
        config_file = self.output_dir / "test_config.py"
        with open(config_file, 'w') as f:
            f.write(config_content)
        
        print(f"Created test configuration: {config_file}")
        return True
    
    def create_regression_test_script(self) -> bool:
        """Create convenient regression test script"""
        
        print("\n=== CREATING REGRESSION TEST SCRIPT ===")
        
        script_content = f'''#!/usr/bin/env python3
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
        cmd = ["python", "tests/test_ecod_regression.py", "--proteins"] + {self.production_test_proteins}
        return subprocess.call(cmd)

if __name__ == "__main__":
    exit(main())
'''
        
        script_file = Path("run_production_tests.py")
        with open(script_file, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_file, 0o755)  # Make executable
        print(f"Created regression test script: {script_file}")
        return True
    
    def create_blacklist_file(self) -> bool:
        """Create reference blacklist file for problematic structures"""
        
        print("\n=== CREATING REFERENCE BLACKLIST ===")
        
        # Start with known problematic cases (expand as needed)
        blacklist_content = '''pdb_id,chain_id,reason,date_added,added_by
# Known problematic reference structures for domain decomposition
# Format: pdb_id,chain_id,reason,date_added,added_by

# Example entries (add real problematic cases as discovered):
# 1abc,A,incomplete_structure,2025-01-01,setup_script
# 2def,B,non_standard_residues,2025-01-01,setup_script
'''
        
        blacklist_file = self.output_dir / "reference_blacklist.csv"
        with open(blacklist_file, 'w') as f:
            f.write(blacklist_content)
        
        print(f"Created reference blacklist: {blacklist_file}")
        print("   (Currently empty - add problematic structures as discovered)")
        return True
    
    def full_setup(self, verbose: bool = False) -> bool:
        """Run complete production testing setup"""
        
        print("MINI PYECOD PRODUCTION TESTING SETUP")
        print("=" * 60)
        
        success_count = 0
        total_steps = 6
        
        # 1. Validate environment
        validation = self.validate_environment()
        if not all(validation.values()):
            print("❌ Environment validation failed - please fix issues above")
            return False
        success_count += 1
        
        # 2. Setup range cache data
        if self.setup_range_cache_data(verbose):
            success_count += 1
        
        # 3. Setup ECOD classifications
        if self.setup_ecod_classifications(verbose):
            success_count += 1
        
        # 4. Create test configuration
        if self.create_test_configuration():
            success_count += 1
        
        # 5. Create regression test script
        if self.create_regression_test_script():
            success_count += 1
        
        # 6. Create blacklist file
        if self.create_blacklist_file():
            success_count += 1
        
        print("\n" + "=" * 60)
        print(f"SETUP COMPLETE: {success_count}/{total_steps} steps successful")
        print("=" * 60)
        
        if success_count == total_steps:
            print("✅ Production testing environment ready!")
            print("\nNext steps:")
            print("1. Validate setup: python run_production_tests.py --validate-setup")
            print("2. Run primary test: python run_production_tests.py --primary -v")
            print("3. Run full suite: python run_production_tests.py --suite -v")
            print("4. Check results in: regression_output/regression_results.json")
            return True
        else:
            print("❌ Setup incomplete - please fix issues above")
            return False
    
    def test_proteins_only_setup(self, verbose: bool = False) -> bool:
        """Quick setup for just the test proteins (faster for development)"""
        
        print("MINI PYECOD TEST PROTEINS SETUP")
        print("=" * 50)
        
        success = True
        
        # Just create essential files for test proteins
        if os.path.exists(self.range_cache_file):
            success &= self.setup_range_cache_data(verbose)
        else:
            print("⚠️  Range cache not available - skipping range cache data")
        
        if os.path.exists(self.domains_file):
            success &= self.setup_ecod_classifications(verbose)
        else:
            print("⚠️  ECOD domains file not available - skipping ECOD classifications")
        
        success &= self.create_test_configuration()
        
        return success

def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Setup mini_pyecod production testing environment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Setup Options:
  --full-setup         Complete setup with all data files (recommended)
    --test-proteins-only Quick setup for test proteins only (development)
  --validate           Just validate environment (no file generation)
    
Examples:
  python setup_production_testing.py --full-setup
    python setup_production_testing.py --test-proteins-only -v
  python setup_production_testing.py --validate
          '''
      )
    
    parser.add_argument('--full-setup', action='store_true',
                       help='Complete production testing setup')
    parser.add_argument('--test-proteins-only', action='store_true',
                       help='Quick setup for test proteins only')
    parser.add_argument('--validate', action='store_true',
                       help='Validate environment only')
    parser.add_argument('--output-dir', default='test_data',
                       help='Output directory for generated files')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    if not any([args.full_setup, args.test_proteins_only, args.validate]):
        parser.print_help()
        print("\n❌ Please specify a setup option")
        return 1
    
    setup = ProductionTestSetup(args.output_dir)
    
    if args.validate:
        validation = setup.validate_environment()
        return 0 if all(validation.values()) else 1
    
    elif args.full_setup:
        success = setup.full_setup(args.verbose)
        return 0 if success else 1
    
    elif args.test_proteins_only:
        success = setup.test_proteins_only_setup(args.verbose)
        return 0 if success else 1

if __name__ == "__main__":
    exit(main())
