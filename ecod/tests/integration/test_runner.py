#!/usr/bin/env python3
"""
Domain Partition Test Runner and Configuration

This script provides the infrastructure to run domain partition integration tests
and manage test configurations, baselines, and regression detection.
"""

import os
import sys
import json
import yaml
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
import logging

# Test configuration templates
TEST_CONFIGS = {
    'default': {
        'database': {
            'host': 'localhost',  # PostgreSQL is listening on 127.0.0.1
            'port': 5432,         # Confirmed port from ss output
            'database': 'ecod_test',
            'user': 'test_user',
            'password': 'test_pass'
        },
        'reference': {
            'current_version': 'develop291'
        },
        'partition': {
            'confidence_thresholds': {
                'high': 0.9,
                'medium': 0.7,
                'low': 0.5
            },
            'evidence_weights': {
                'domain_blast': 3.0,
                'hhsearch': 2.5,
                'chain_blast': 2.0,
                'blast': 1.5
            },
            'overlap_tolerance': 0.15,
            'min_domain_size': 20,
            'peptide_threshold': 50
        }
    },

    'hhsearch_heavy': {
        'partition': {
            'evidence_weights': {
                'domain_blast': 2.0,
                'hhsearch': 4.0,
                'chain_blast': 1.5,
                'blast': 1.0
            }
        }
    },

    'blast_heavy': {
        'partition': {
            'evidence_weights': {
                'domain_blast': 4.0,
                'hhsearch': 1.5,
                'chain_blast': 3.0,
                'blast': 2.0
            }
        }
    },

    'strict_thresholds': {
        'partition': {
            'confidence_thresholds': {
                'high': 0.95,
                'medium': 0.8,
                'low': 0.6
            }
        }
    },

    'lenient_thresholds': {
        'partition': {
            'confidence_thresholds': {
                'high': 0.8,
                'medium': 0.6,
                'low': 0.4
            }
        }
    }
}

# Golden dataset definitions
GOLDEN_DATASETS = {
    'single_domain': [
        {
            'pdb_id': '1cbs',
            'chain_id': 'A',
            'sequence_length': 141,
            'expected_domains': 1,
            'expected_classification': 'a.39.1.5',
            'description': 'Small calcium-binding domain protein'
        },
        {
            'pdb_id': '1ubq',
            'chain_id': 'A',
            'sequence_length': 76,
            'expected_domains': 1,
            'expected_classification': 'a.5.2.1',
            'description': 'Ubiquitin - well-characterized single domain'
        }
    ],

    'multi_domain': [
        {
            'pdb_id': '2pth',
            'chain_id': 'A',
            'sequence_length': 393,
            'expected_domains': 2,
            'expected_classifications': ['c.37.1.13', 'c.37.1.13'],
            'description': 'Two-domain protein with clear boundaries'
        },
        {
            'pdb_id': '1igd',
            'chain_id': 'A',
            'sequence_length': 425,
            'expected_domains': 3,
            'expected_classifications': ['c.37.1.1', 'c.37.1.1', 'c.37.1.1'],
            'description': 'Immunoglobulin with multiple domains'
        }
    ],

    'complex_cases': [
        {
            'pdb_id': '3hhp',
            'chain_id': 'A',
            'sequence_length': 506,
            'expected_domains': 3,
            'expected_overlaps': True,
            'description': 'Complex architecture with potential overlaps'
        }
    ],

    'peptides': [
        {
            'pdb_id': '1pep',
            'chain_id': 'A',
            'sequence_length': 25,
            'expected_peptide': True,
            'description': 'Short peptide sequence'
        }
    ],

    'edge_cases': [
        {
            'pdb_id': '1edge',
            'chain_id': 'A',
            'sequence_length': 55,
            'expected_domains': 0,
            'expected_unclassified': True,
            'description': 'Borderline peptide/protein size'
        }
    ]
}


class TestConfigManager:
    """Manages test configurations and golden datasets"""

    def __init__(self, test_root_dir: str):
        self.test_root = Path(test_root_dir)
        self.config_dir = self.test_root / "configs"
        self.datasets_dir = self.test_root / "datasets"
        self.baselines_dir = self.test_root / "baselines"
        self.results_dir = self.test_root / "results"

        # Create directories
        for dir_path in [self.config_dir, self.datasets_dir, self.baselines_dir, self.results_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

        self.logger = logging.getLogger(__name__)

    def setup_test_environment(self):
        """Set up the complete test environment"""
        self.logger.info("Setting up test environment...")

        # Create test configurations
        self._create_test_configs()

        # Create golden datasets
        self._create_golden_datasets()

        # Create pytest configuration
        self._create_pytest_config()

        # Create test database setup script
        self._create_db_setup_script()

        self.logger.info("Test environment setup complete")

    def _create_test_configs(self):
        """Create test configuration files"""
        for config_name, config_data in TEST_CONFIGS.items():
            config_file = self.config_dir / f"{config_name}.yml"

            # Merge with default if not default
            if config_name != 'default':
                merged_config = self._deep_merge(TEST_CONFIGS['default'].copy(), config_data)
            else:
                merged_config = config_data

            with open(config_file, 'w') as f:
                yaml.dump(merged_config, f, default_flow_style=False, indent=2)

            self.logger.info(f"Created config: {config_file}")

    def _create_golden_datasets(self):
        """Create golden dataset files"""
        # Create dataset catalog
        catalog = {
            'created': datetime.now().isoformat(),
            'datasets': GOLDEN_DATASETS,
            'total_proteins': sum(len(proteins) for proteins in GOLDEN_DATASETS.values())
        }

        catalog_file = self.datasets_dir / "catalog.json"
        with open(catalog_file, 'w') as f:
            json.dump(catalog, f, indent=2)

        # Create individual dataset files for each category
        for category, proteins in GOLDEN_DATASETS.items():
            dataset_file = self.datasets_dir / f"{category}.json"
            with open(dataset_file, 'w') as f:
                json.dump(proteins, f, indent=2)

            self.logger.info(f"Created dataset: {dataset_file} ({len(proteins)} proteins)")

    def _create_pytest_config(self):
        """Create pytest configuration"""
        pytest_ini_content = """
[tool:pytest]
testpaths = tests/integration
python_files = test_*.py *_test.py
python_classes = Test* *Tests
python_functions = test_*
addopts =
    -v
    --tb=short
    --strict-markers
    --disable-warnings
markers =
    golden: Golden dataset regression tests
    evidence: Evidence processing tests
    weights: Evidence weight sensitivity tests
    service: Service integration tests
    performance: Performance regression tests
    slow: Slow running tests
    integration: Integration tests
filterwarnings =
    ignore::DeprecationWarning
    ignore::PendingDeprecationWarning
"""

        pytest_file = self.test_root / "pytest.ini"
        with open(pytest_file, 'w') as f:
            f.write(pytest_ini_content.strip())

    def _create_db_setup_script(self):
        """Create database setup script for tests"""
        db_setup_content = '''#!/bin/bash
# Database setup script for domain partition tests

set -e

# PostgreSQL paths - adjust if needed
PG_BIN_PATH="/sw/apps/postgresql-17.4/bin"
if [ -d "$PG_BIN_PATH" ]; then
    export PATH="$PG_BIN_PATH:$PATH"
fi

# Database configuration - adjust as needed
DB_NAME="ecod_test"
DB_USER="test_user"
DB_PASS="test_pass"
DB_HOST="${PGHOST:-lotta}"     # Use PGHOST env var or default to lotta
DB_PORT="${PGPORT:-5432}"      # Use PGPORT env var or default to 5432
DB_SUPERUSER="${PGUSER:-rschaeff}"  # Use PGUSER env var or default to rschaeff

echo "Setting up test database on $DB_HOST:$DB_PORT as user $DB_SUPERUSER..."

# Check if PostgreSQL utilities are available
if ! command -v createdb &> /dev/null; then
    echo "Error: PostgreSQL utilities not found in PATH"
    echo "Please add PostgreSQL bin directory to PATH:"
    echo "export PATH=\"/sw/apps/postgresql-17.4/bin:\$PATH\""
    exit 1
fi

# Test connection first
echo "Testing PostgreSQL connection..."
if ! pg_isready -p $DB_PORT -h $DB_HOST; then
    echo "Error: Cannot connect to PostgreSQL on $DB_HOST:$DB_PORT"
    echo "Please check:"
    echo "1. PostgreSQL is running"
    echo "2. Correct host (set PGHOST=<host> if not lotta)"
    echo "3. Correct port (set PGPORT=<port> if not 5432)"
    echo "4. PostgreSQL is accepting connections"
    echo ""
    echo "Try testing connection manually:"
    echo "  psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d postgres -c '\\l'"
    exit 1
fi

# Check if we can connect as superuser
echo "Testing superuser connection..."
if ! psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d postgres -c "SELECT version();" > /dev/null 2>&1; then
    echo "Error: Cannot connect as user $DB_SUPERUSER"
    echo "You may need to:"
    echo "1. Use a different superuser (set PGUSER=<username>)"
    echo "2. Set up password authentication"
    echo "3. Check pg_hba.conf settings"
    echo ""
    echo "Try manual connection:"
    echo "  psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d postgres"
    exit 1
fi

# Create test database
echo "Creating test database..."
createdb -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT $DB_NAME || echo "Database '$DB_NAME' may already exist"

# Check if test database was created or already exists
if ! psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d $DB_NAME -c "SELECT 1;" > /dev/null 2>&1; then
    echo "Error: Cannot access test database $DB_NAME"
    exit 1
fi

# Create test user
echo "Creating test user..."
psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d $DB_NAME -c "
    DO \$\$
    BEGIN
        IF NOT EXISTS (SELECT FROM pg_catalog.pg_roles WHERE rolname = '$DB_USER') THEN
            CREATE USER $DB_USER WITH PASSWORD '$DB_PASS';
            RAISE NOTICE 'User $DB_USER created';
        ELSE
            RAISE NOTICE 'User $DB_USER already exists';
        END IF;
    END
    \$\$;
"

# Grant permissions
echo "Granting permissions..."
psql -U $DB_SUPERUSER -h $DB_HOST -p $DB_PORT -d $DB_NAME -c "
    GRANT ALL PRIVILEGES ON DATABASE $DB_NAME TO $DB_USER;
    GRANT ALL PRIVILEGES ON SCHEMA public TO $DB_USER;
    GRANT CREATE ON SCHEMA public TO $DB_USER;
"

# Create test schema (minimal version of production schema)
psql -U postgres -d $DB_NAME -c "
CREATE SCHEMA IF NOT EXISTS ecod_schema;
CREATE SCHEMA IF NOT EXISTS pdb_analysis;

-- Basic tables for testing
CREATE TABLE IF NOT EXISTS ecod_schema.batch (
    id SERIAL PRIMARY KEY,
    batch_name VARCHAR(255),
    base_path TEXT,
    ref_version VARCHAR(50),
    status VARCHAR(50),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ecod_schema.protein (
    id SERIAL PRIMARY KEY,
    pdb_id VARCHAR(4),
    chain_id VARCHAR(10),
    sequence_length INTEGER,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ecod_schema.process_status (
    id SERIAL PRIMARY KEY,
    batch_id INTEGER REFERENCES ecod_schema.batch(id),
    protein_id INTEGER REFERENCES ecod_schema.protein(id),
    current_stage VARCHAR(100),
    status VARCHAR(50),
    error_message TEXT,
    is_representative BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ecod_schema.process_file (
    id SERIAL PRIMARY KEY,
    process_id INTEGER REFERENCES ecod_schema.process_status(id),
    file_type VARCHAR(100),
    file_path TEXT,
    file_exists BOOLEAN DEFAULT FALSE,
    file_size BIGINT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Grant permissions on all tables
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA ecod_schema TO $DB_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA pdb_analysis TO $DB_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA ecod_schema TO $DB_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA pdb_analysis TO $DB_USER;
"

echo "Test database setup complete"
'''
        
        db_script_file = self.test_root / "setup_test_db.sh"
        with open(db_script_file, 'w') as f:
            f.write(db_setup_content)
        
        # Make executable
        os.chmod(db_script_file, 0o755)

    def _deep_merge(self, base_dict: Dict, update_dict: Dict) -> Dict:
        """Deep merge two dictionaries"""
        result = base_dict.copy()
        
        for key, value in update_dict.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = self._deep_merge(result[key], value)
            else:
                result[key] = value
        
        return result


class RegressionTestRunner:
    """Runs regression tests and manages baselines"""
    
    def __init__(self, test_root_dir: str):
        self.test_root = Path(test_root_dir)
        self.baselines_dir = self.test_root / "baselines"
        self.results_dir = self.test_root / "results"
        self.logger = logging.getLogger(__name__)

    def establish_baselines(self, config_name: str = 'default'):
        """Establish baseline results for regression testing"""
        self.logger.info(f"Establishing baselines with config: {config_name}")
        
        # Run golden dataset tests to establish baselines
        baseline_run_id = datetime.now().strftime('%Y%m%d_%H%M%S')
        baseline_dir = self.baselines_dir / config_name / baseline_run_id
        baseline_dir.mkdir(parents=True, exist_ok=True)
        
        # Run pytest with baseline establishment mode
        cmd = [
            'python', '-m', 'pytest',
            'tests/integration/test_domain_partition_integration.py::GoldenDatasetTests',
            '-v',
            '--tb=short',
            f'--config={config_name}',
            f'--baseline-dir={baseline_dir}',
            '--establish-baseline'
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.test_root)
            
            if result.returncode == 0:
                self.logger.info(f"Baselines established successfully in {baseline_dir}")
                
                # Create baseline metadata
                metadata = {
                    'created': datetime.now().isoformat(),
                    'config': config_name,
                    'run_id': baseline_run_id,
                    'command': ' '.join(cmd),
                    'status': 'success'
                }
                
                with open(baseline_dir / 'metadata.json', 'w') as f:
                    json.dump(metadata, f, indent=2)
                
                # Update latest baseline link
                latest_link = self.baselines_dir / config_name / 'latest'
                if latest_link.exists():
                    latest_link.unlink()
                latest_link.symlink_to(baseline_run_id)
                
            else:
                self.logger.error(f"Baseline establishment failed: {result.stderr}")
                return False
                
        except Exception as e:
            self.logger.error(f"Error establishing baselines: {e}")
            return False
        
        return True

    def run_regression_tests(self, config_name: str = 'default', 
                           test_categories: List[str] = None) -> Dict[str, Any]:
        """Run regression tests against established baselines"""
        self.logger.info(f"Running regression tests with config: {config_name}")
        
        if test_categories is None:
            test_categories = ['golden', 'evidence', 'weights', 'service', 'performance']
        
        # Create results directory
        run_id = datetime.now().strftime('%Y%m%d_%H%M%S')
        results_dir = self.results_dir / config_name / run_id
        results_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        
        for category in test_categories:
            self.logger.info(f"Running {category} tests...")
            
            # Build pytest command for category
            if category == 'golden':
                test_class = 'GoldenDatasetTests'
            elif category == 'evidence':
                test_class = 'EvidenceVariationTests'
            elif category == 'weights':
                test_class = 'ConfidenceWeightTests'
            elif category == 'service':
                test_class = 'ServiceIntegrationTests'
            elif category == 'performance':
                test_class = 'PerformanceRegressionTests'
            else:
                continue
            
            cmd = [
                'python', '-m', 'pytest',
                f'tests/integration/test_domain_partition_integration.py::{test_class}',
                '-v',
                '--tb=short',
                f'--config={config_name}',
                f'--results-dir={results_dir}',
                '--compare-baseline'
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.test_root)
                
                # Parse test results
                category_results = self._parse_test_output(result.stdout, result.stderr)
                category_results['return_code'] = result.returncode
                results[category] = category_results
                
            except Exception as e:
                self.logger.error(f"Error running {category} tests: {e}")
                results[category] = {'error': str(e), 'return_code': -1}
        
        # Generate comprehensive report
        report = self._generate_regression_report(results, config_name, run_id)
        
        # Save report
        report_file = results_dir / 'regression_report.json'
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        self.logger.info(f"Regression test results saved to {report_file}")
        
        return report

    def compare_configurations(self, config_names: List[str]) -> Dict[str, Any]:
        """Compare results across different configurations"""
        self.logger.info(f"Comparing configurations: {config_names}")
        
        comparison_results = {}
        
        for config_name in config_names:
            results = self.run_regression_tests(config_name, ['golden', 'evidence'])
            comparison_results[config_name] = results
        
        # Generate comparison report
        comparison_report = {
            'timestamp': datetime.now().isoformat(),
            'configurations': config_names,
            'comparison_results': comparison_results,
            'summary': self._summarize_config_comparison(comparison_results)
        }
        
        # Save comparison report
        comparison_file = self.results_dir / f'config_comparison_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
        with open(comparison_file, 'w') as f:
            json.dump(comparison_report, f, indent=2)
        
        return comparison_report

    def _parse_test_output(self, stdout: str, stderr: str) -> Dict[str, Any]:
        """Parse pytest output to extract test results"""
        lines = stdout.split('\n')
        
        results = {
            'tests_run': 0,
            'tests_passed': 0,
            'tests_failed': 0,
            'tests_skipped': 0,
            'failures': [],
            'errors': []
        }
        
        # Simple parsing - in real implementation, you'd want more sophisticated parsing
        for line in lines:
            if '::' in line and ('PASSED' in line or 'FAILED' in line or 'SKIPPED' in line):
                results['tests_run'] += 1
                if 'PASSED' in line:
                    results['tests_passed'] += 1
                elif 'FAILED' in line:
                    results['tests_failed'] += 1
                    results['failures'].append(line)
                elif 'SKIPPED' in line:
                    results['tests_skipped'] += 1
        
        if stderr:
            results['errors'].append(stderr)
        
        return results

    def _generate_regression_report(self, results: Dict[str, Any], 
                                  config_name: str, run_id: str) -> Dict[str, Any]:
        """Generate comprehensive regression report"""
        total_tests = sum(r.get('tests_run', 0) for r in results.values() if isinstance(r, dict))
        total_passed = sum(r.get('tests_passed', 0) for r in results.values() if isinstance(r, dict))
        total_failed = sum(r.get('tests_failed', 0) for r in results.values() if isinstance(r, dict))
        
        report = {
            'timestamp': datetime.now().isoformat(),
            'config_name': config_name,
            'run_id': run_id,
            'summary': {
                'total_tests': total_tests,
                'passed': total_passed,
                'failed': total_failed,
                'success_rate': total_passed / total_tests if total_tests > 0 else 0
            },
            'category_results': results,
            'regressions_detected': self._detect_regressions(results),
            'recommendations': self._generate_recommendations(results)
        }
        
        return report

    def _detect_regressions(self, results: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Detect regressions from test results"""
        regressions = []
        
        for category, category_results in results.items():
            if isinstance(category_results, dict) and category_results.get('tests_failed', 0) > 0:
                regressions.append({
                    'category': category,
                    'type': 'test_failure',
                    'count': category_results['tests_failed'],
                    'details': category_results.get('failures', [])
                })
        
        return regressions

    def _generate_recommendations(self, results: Dict[str, Any]) -> List[str]:
        """Generate recommendations based on test results"""
        recommendations = []
        
        # Check for golden dataset failures
        golden_results = results.get('golden', {})
        if isinstance(golden_results, dict) and golden_results.get('tests_failed', 0) > 0:
            recommendations.append(
                "Golden dataset tests failed - core domain partitioning logic may be broken. "
                "Review recent changes to evidence processing and confidence calculation."
            )
        
        # Check for evidence processing issues
        evidence_results = results.get('evidence', {})
        if isinstance(evidence_results, dict) and evidence_results.get('tests_failed', 0) > 0:
            recommendations.append(
                "Evidence processing tests failed - check evidence parsing and integration logic."
            )
        
        # Check for weight sensitivity issues
        weights_results = results.get('weights', {})
        if isinstance(weights_results, dict) and weights_results.get('tests_failed', 0) > 0:
            recommendations.append(
                "Evidence weight tests failed - confidence calculation changes may have broken "
                "expected behavior. Review weight application and normalization."
            )
        
        # Check for performance regressions
        perf_results = results.get('performance', {})
        if isinstance(perf_results, dict) and perf_results.get('tests_failed', 0) > 0:
            recommendations.append(
                "Performance tests failed - processing time or memory usage may have regressed. "
                "Profile recent changes for performance impact."
            )
        
        if not recommendations:
            recommendations.append("All tests passed - no immediate issues detected.")
        
        return recommendations

    def _summarize_config_comparison(self, comparison_results: Dict[str, Any]) -> Dict[str, Any]:
        """Summarize comparison across configurations"""
        summary = {
            'best_performing': None,
            'most_reliable': None,
            'performance_differences': {},
            'reliability_differences': {}
        }
        
        # Find best performing configuration
        best_success_rate = 0
        for config_name, results in comparison_results.items():
            if isinstance(results, dict) and 'summary' in results:
                success_rate = results['summary'].get('success_rate', 0)
                if success_rate > best_success_rate:
                    best_success_rate = success_rate
                    summary['best_performing'] = config_name
        
        return summary


def main():
    """Main entry point for test runner"""
    parser = argparse.ArgumentParser(description="Domain Partition Test Runner")
    
    parser.add_argument('--test-root', type=str, default='tests/integration',
                       help='Root directory for integration tests')
    parser.add_argument('--config', type=str, default='default',
                       help='Test configuration to use')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Setup command
    setup_parser = subparsers.add_parser('setup', help='Set up test environment')
    
    # Baseline command
    baseline_parser = subparsers.add_parser('baseline', help='Establish test baselines')
    baseline_parser.add_argument('--config', type=str, default='default',
                                help='Configuration for baseline')
    
    # Test command
    test_parser = subparsers.add_parser('test', help='Run regression tests')
    test_parser.add_argument('--config', type=str, default='default',
                            help='Configuration for tests')
    test_parser.add_argument('--categories', nargs='+',
                            choices=['golden', 'evidence', 'weights', 'service', 'performance'],
                            help='Test categories to run')
    
    # Compare command
    compare_parser = subparsers.add_parser('compare', help='Compare configurations')
    compare_parser.add_argument('--configs', nargs='+', required=True,
                               help='Configurations to compare')
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Create test manager and runner
    config_manager = TestConfigManager(args.test_root)
    test_runner = RegressionTestRunner(args.test_root)
    
    if args.command == 'setup':
        config_manager.setup_test_environment()
        print("Test environment setup complete!")
        print(f"Test root: {args.test_root}")
        print("Next steps:")
        print("1. Run './setup_test_db.sh' to set up test database")
        print("2. Run 'python test_runner.py baseline' to establish baselines")
        print("3. Run 'python test_runner.py test' to run regression tests")
        
    elif args.command == 'baseline':
        success = test_runner.establish_baselines(args.config)
        if success:
            print(f"Baselines established for config: {args.config}")
        else:
            print("Failed to establish baselines")
            return 1
            
    elif args.command == 'test':
        report = test_runner.run_regression_tests(args.config, args.categories)
        
        print(f"\nRegression Test Results (Config: {args.config})")
        print("=" * 50)
        print(f"Total Tests: {report['summary']['total_tests']}")
        print(f"Passed: {report['summary']['passed']}")
        print(f"Failed: {report['summary']['failed']}")
        print(f"Success Rate: {report['summary']['success_rate']:.1%}")
        
        if report['regressions_detected']:
            print(f"\nRegressions Detected: {len(report['regressions_detected'])}")
            for regression in report['regressions_detected']:
                print(f"  - {regression['category']}: {regression['type']}")
        
        print(f"\nRecommendations:")
        for rec in report['recommendations']:
            print(f"  - {rec}")
            
    elif args.command == 'compare':
        report = test_runner.compare_configurations(args.configs)
        
        print(f"\nConfiguration Comparison")
        print("=" * 50)
        
        for config_name in args.configs:
            if config_name in report['comparison_results']:
                results = report['comparison_results'][config_name]
                if 'summary' in results:
                    summary = results['summary']
                    print(f"{config_name:15s}: {summary['success_rate']:.1%} success rate "
                          f"({summary['passed']}/{summary['total_tests']} tests)")
        
        if report['summary']['best_performing']:
            print(f"\nBest performing: {report['summary']['best_performing']}")
    
    else:
        parser.print_help()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
