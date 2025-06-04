#!/usr/bin/env python3
"""
Mini PyECOD Integration Testing Suite

Comprehensive testing of import, comparison, and assessment scripts
on the 6% subset to identify bugs before full pipeline completion.

Usage:
    python test_integration.py --full-test          # Run all tests
    python test_integration.py --test-import        # Test import only
    python test_integration.py --test-comparison    # Test comparison only  
    python test_integration.py --test-assessment    # Test quality assessment
    python test_integration.py --validate-setup     # Validate environment
"""

import os
import sys
import subprocess
import argparse
import yaml
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MiniIntegrationTester:
    """Comprehensive testing suite for mini PyECOD integration"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        self.config = self._load_config(config_path)
        self.test_results = {}
        
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration"""
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            logger.error(f"Config file not found: {config_path}")
            raise
    
    def validate_environment(self) -> Dict[str, bool]:
        """Validate that all required components are available"""
        
        logger.info("🔍 Validating environment setup...")
        validation = {}
        
        # Check config file
        validation['config_file'] = os.path.exists("config.local.yml")
        logger.info(f"  Config file: {'✓' if validation['config_file'] else '✗'}")
        
        # Check database connectivity
        try:
            import psycopg2
            conn = psycopg2.connect(**self.config['database'])
            conn.close()
            validation['database_connection'] = True
            logger.info("  Database connection: ✓")
        except Exception as e:
            validation['database_connection'] = False
            logger.error(f"  Database connection: ✗ ({e})")
        
        # Check batch directories
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        validation['batch_directories'] = batch_base.exists()
        logger.info(f"  Batch directories: {'✓' if validation['batch_directories'] else '✗'}")
        
        # Check for mini results
        mini_results_found = False
        if batch_base.exists():
            for batch_dir in batch_base.iterdir():
                mini_domains_dir = batch_dir / "mini_domains"
                if mini_domains_dir.exists() and list(mini_domains_dir.glob("*.mini.domains.xml")):
                    mini_results_found = True
                    break
        validation['mini_results_available'] = mini_results_found
        logger.info(f"  Mini results available: {'✓' if mini_results_found else '✗'}")
        
        # Check script files
        scripts = [
            "mini_production/import_results.py",
            "mini_production/compare_results.py", 
            "mini_production/assess_quality.py",
            "mini_production/monitor_progress.py"
        ]
        
        for script in scripts:
            exists = os.path.exists(script)
            validation[f'script_{Path(script).stem}'] = exists
            logger.info(f"  {script}: {'✓' if exists else '✗'}")
        
        # Summary
        passed = sum(validation.values())
        total = len(validation)
        logger.info(f"\n📊 Environment validation: {passed}/{total} checks passed")
        
        if passed == total:
            logger.info("✅ Environment ready for testing")
        else:
            logger.warning("⚠️  Some environment issues detected")
        
        return validation
    
    def test_import_functionality(self, test_limit: int = 10) -> Dict[str, any]:
        """Test import functionality with collision detection"""
        
        logger.info(f"🧪 Testing import functionality (limit: {test_limit})...")
        test_result = {
            'success': False,
            'collision_check_passed': False,
            'import_test_passed': False,
            'verification_passed': False,
            'errors': []
        }
        
        try:
            # 1. Test collision detection first (safe - no data changes)
            logger.info("  Step 1: Testing collision detection...")
            result = subprocess.run([
                "python", "mini_production/import_results.py",
                "--check-collisions", 
                "--limit", str(test_limit)
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                test_result['collision_check_passed'] = True
                logger.info("    ✓ Collision detection working")
                
                # Parse output for collision info
                collision_info = result.stdout
                if "potential collisions" in collision_info:
                    logger.info(f"    Found collision info in output")
                
            else:
                test_result['errors'].append(f"Collision check failed: {result.stderr}")
                logger.error(f"    ✗ Collision check failed: {result.stderr}")
                return test_result
            
            # 2. Test actual import with separate strategy (safe)
            logger.info("  Step 2: Testing import with 'separate' strategy...")
            result = subprocess.run([
                "python", "mini_production/import_results.py",
                "--import-all",
                "--collision-strategy", "separate",
                "--limit", str(test_limit)
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                test_result['import_test_passed'] = True
                logger.info("    ✓ Import test working")
                
                # Parse import results
                import_info = result.stdout
                if "imported" in import_info.lower():
                    logger.info(f"    Import completed successfully")
                
            else:
                test_result['errors'].append(f"Import test failed: {result.stderr}")
                logger.error(f"    ✗ Import test failed: {result.stderr}")
                return test_result
            
            # 3. Test verification
            logger.info("  Step 3: Testing import verification...")
            result = subprocess.run([
                "python", "mini_production/import_results.py",
                "--verify"
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                test_result['verification_passed'] = True
                logger.info("    ✓ Verification working")
            else:
                test_result['errors'].append(f"Verification failed: {result.stderr}")
                logger.error(f"    ✗ Verification failed: {result.stderr}")
            
            test_result['success'] = (test_result['collision_check_passed'] and 
                                    test_result['import_test_passed'] and 
                                    test_result['verification_passed'])
            
        except Exception as e:
            test_result['errors'].append(f"Import test exception: {str(e)}")
            logger.error(f"  ✗ Import test exception: {e}")
        
        status = "✅ PASSED" if test_result['success'] else "❌ FAILED"
        logger.info(f"📊 Import functionality test: {status}")
        
        return test_result
    
    def test_comparison_functionality(self, sample_size: int = 20) -> Dict[str, any]:
        """Test comparison functionality"""
        
        logger.info(f"🧪 Testing comparison functionality (sample: {sample_size})...")
        test_result = {
            'success': False,
            'comparison_completed': False,
            'results_generated': False,
            'errors': []
        }
        
        try:
            # Test sample comparison
            logger.info("  Step 1: Testing sample comparison...")
            result = subprocess.run([
                "python", "mini_production/compare_results.py",
                "--sample-comparison", str(sample_size)
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                test_result['comparison_completed'] = True
                logger.info("    ✓ Comparison completed")
                
                # Check for expected output patterns
                output = result.stdout
                if any(keyword in output.lower() for keyword in ["mini wins", "comparison", "improvement"]):
                    test_result['results_generated'] = True
                    logger.info("    ✓ Results generated")
                    
                    # Extract some key metrics from output
                    lines = output.split('\n')
                    for line in lines:
                        if "mini wins" in line.lower() or "comparison" in line.lower():
                            logger.info(f"    Sample result: {line.strip()}")
                            break
                else:
                    logger.warning("    ⚠️  Unexpected output format")
                
            else:
                test_result['errors'].append(f"Comparison failed: {result.stderr}")
                logger.error(f"    ✗ Comparison failed: {result.stderr}")
        
        except Exception as e:
            test_result['errors'].append(f"Comparison test exception: {str(e)}")
            logger.error(f"  ✗ Comparison test exception: {e}")
        
        test_result['success'] = test_result['comparison_completed'] and test_result['results_generated']
        
        status = "✅ PASSED" if test_result['success'] else "❌ FAILED"
        logger.info(f"📊 Comparison functionality test: {status}")
        
        return test_result
    
    def test_assessment_functionality(self) -> Dict[str, any]:
        """Test quality assessment functionality"""
        
        logger.info("🧪 Testing quality assessment functionality...")
        test_result = {
            'success': False,
            'scan_completed': False,
            'assessment_completed': False,
            'export_completed': False,
            'errors': []
        }
        
        try:
            # Test scanning for results
            logger.info("  Step 1: Testing result scanning...")
            result = subprocess.run([
                "python", "mini_production/assess_quality.py",
                "--scan-all"
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                test_result['scan_completed'] = True
                logger.info("    ✓ Scanning completed")
                
                # Check for quality metrics in output
                output = result.stdout
                if any(keyword in output.lower() for keyword in ["quality", "assessment", "tier", "production"]):
                    test_result['assessment_completed'] = True
                    logger.info("    ✓ Quality assessment completed")
                    
                    # Extract quality summary
                    lines = output.split('\n')
                    for line in lines:
                        if "production ready" in line.lower() or "quality" in line.lower():
                            logger.info(f"    Quality metric: {line.strip()}")
                
            else:
                test_result['errors'].append(f"Assessment failed: {result.stderr}")
                logger.error(f"    ✗ Assessment failed: {result.stderr}")
                return test_result
            
            # Test export functionality
            logger.info("  Step 2: Testing export functionality...")
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as tmp_file:
                tmp_export_path = tmp_file.name
            
            try:
                result = subprocess.run([
                    "python", "mini_production/assess_quality.py",
                    "--scan-all",
                    "--export-production",
                    "--output-file", tmp_export_path
                ], capture_output=True, text=True)
                
                if result.returncode == 0 and os.path.exists(tmp_export_path):
                    test_result['export_completed'] = True
                    logger.info("    ✓ Export completed")
                    
                    # Check export file content
                    try:
                        import json
                        with open(tmp_export_path, 'r') as f:
                            export_data = json.load(f)
                        
                        if 'proteins' in export_data and 'production_ready' in export_data:
                            logger.info(f"    Export contains {export_data.get('production_ready', 0)} production-ready results")
                    except Exception as e:
                        logger.warning(f"    Could not parse export file: {e}")
                
                else:
                    test_result['errors'].append(f"Export failed: {result.stderr}")
                    logger.error(f"    ✗ Export failed: {result.stderr}")
            
            finally:
                # Clean up temp file
                if os.path.exists(tmp_export_path):
                    os.unlink(tmp_export_path)
            
        except Exception as e:
            test_result['errors'].append(f"Assessment test exception: {str(e)}")
            logger.error(f"  ✗ Assessment test exception: {e}")
        
        test_result['success'] = (test_result['scan_completed'] and 
                                test_result['assessment_completed'] and 
                                test_result['export_completed'])
        
        status = "✅ PASSED" if test_result['success'] else "❌ FAILED"
        logger.info(f"📊 Assessment functionality test: {status}")
        
        return test_result
    
    def test_monitoring_functionality(self) -> Dict[str, any]:
        """Test monitoring functionality"""
        
        logger.info("🧪 Testing monitoring functionality...")
        test_result = {
            'success': False,
            'status_check_passed': False,
            'progress_scan_passed': False,
            'import_readiness_passed': False,
            'errors': []
        }
        
        try:
            # Test basic status
            logger.info("  Step 1: Testing status check...")
            result = subprocess.run([
                "python", "mini_production/monitor_progress.py"
            ], capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0:
                test_result['status_check_passed'] = True
                logger.info("    ✓ Status check completed")
                
                # Look for expected output patterns
                output = result.stdout
                if any(keyword in output.lower() for keyword in ["progress", "completed", "slurm"]):
                    test_result['progress_scan_passed'] = True
                    logger.info("    ✓ Progress scanning working")
            else:
                test_result['errors'].append(f"Status check failed: {result.stderr}")
                logger.error(f"    ✗ Status check failed: {result.stderr}")
                return test_result
            
            # Test import readiness
            logger.info("  Step 2: Testing import readiness check...")
            result = subprocess.run([
                "python", "mini_production/monitor_progress.py",
                "--ready-import"
            ], capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0:
                test_result['import_readiness_passed'] = True
                logger.info("    ✓ Import readiness check completed")
            else:
                test_result['errors'].append(f"Import readiness check failed: {result.stderr}")
                logger.error(f"    ✗ Import readiness check failed: {result.stderr}")
            
        except subprocess.TimeoutExpired:
            test_result['errors'].append("Monitoring test timed out")
            logger.error("  ✗ Monitoring test timed out")
        except Exception as e:
            test_result['errors'].append(f"Monitoring test exception: {str(e)}")
            logger.error(f"  ✗ Monitoring test exception: {e}")
        
        test_result['success'] = (test_result['status_check_passed'] and 
                                test_result['progress_scan_passed'] and 
                                test_result['import_readiness_passed'])
        
        status = "✅ PASSED" if test_result['success'] else "❌ FAILED"
        logger.info(f"📊 Monitoring functionality test: {status}")
        
        return test_result
    
    def generate_test_report(self) -> str:
        """Generate comprehensive test report"""
        
        report = []
        report.append("🧪 Mini PyECOD Integration Test Report")
        report.append("=" * 60)
        report.append(f"Test Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append(f"Data Subset: ~6% of aggregate results\n")
        
        # Environment validation
        if 'environment' in self.test_results:
            env_results = self.test_results['environment']
            passed = sum(env_results.values())
            total = len(env_results)
            report.append(f"🔍 Environment Validation: {passed}/{total} checks passed")
            
            for check, result in env_results.items():
                status = "✓" if result else "✗"
                report.append(f"  {check}: {status}")
            report.append("")
        
        # Test results summary
        test_categories = ['import', 'comparison', 'assessment', 'monitoring']
        passed_tests = 0
        total_tests = 0
        
        for category in test_categories:
            if category in self.test_results:
                result = self.test_results[category]
                total_tests += 1
                if result.get('success', False):
                    passed_tests += 1
                
                status = "✅ PASSED" if result.get('success', False) else "❌ FAILED"
                report.append(f"📊 {category.title()} Test: {status}")
                
                # Show errors if any
                if result.get('errors'):
                    report.append("  Errors:")
                    for error in result['errors']:
                        report.append(f"    • {error}")
                report.append("")
        
        # Overall assessment
        report.append("🎯 Overall Assessment:")
        report.append(f"  Tests passed: {passed_tests}/{total_tests}")
        
        if passed_tests == total_tests:
            report.append("  Status: ✅ ALL TESTS PASSED - Ready for full pipeline")
            report.append("  Recommendation: Continue with full import when pipeline completes")
        elif passed_tests >= total_tests * 0.75:
            report.append("  Status: ⚠️  MOSTLY WORKING - Minor issues detected")
            report.append("  Recommendation: Address issues but pipeline can continue")
        else:
            report.append("  Status: ❌ SIGNIFICANT ISSUES - Pipeline may need restart")
            report.append("  Recommendation: Fix issues before full pipeline completion")
        
        report.append("")
        report.append("🚀 Next Steps:")
        
        if passed_tests == total_tests:
            report.append("  1. Monitor ongoing pipeline completion")
            report.append("  2. Prepare for full-scale import")
            report.append("  3. Plan research analysis of results")
        else:
            report.append("  1. Review and fix identified issues")
            report.append("  2. Re-run integration tests")
            report.append("  3. Consider if pipeline restart is needed")
        
        return "\n".join(report)
    
    def run_full_test_suite(self) -> Dict[str, any]:
        """Run all integration tests"""
        
        logger.info("🚀 Starting full integration test suite...")
        
        # 1. Environment validation
        self.test_results['environment'] = self.validate_environment()
        
        # 2. Import functionality
        self.test_results['import'] = self.test_import_functionality(test_limit=5)
        
        # 3. Comparison functionality  
        self.test_results['comparison'] = self.test_comparison_functionality(sample_size=10)
        
        # 4. Assessment functionality
        self.test_results['assessment'] = self.test_assessment_functionality()
        
        # 5. Monitoring functionality
        self.test_results['monitoring'] = self.test_monitoring_functionality()
        
        # Generate and print report
        report = self.generate_test_report()
        print("\n" + report)
        
        return self.test_results


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Mini PyECOD Integration Testing Suite'
    )
    
    parser.add_argument('--full-test', action='store_true',
                       help='Run complete test suite')
    parser.add_argument('--test-import', action='store_true',
                       help='Test import functionality only')
    parser.add_argument('--test-comparison', action='store_true',
                       help='Test comparison functionality only')
    parser.add_argument('--test-assessment', action='store_true',
                       help='Test assessment functionality only')
    parser.add_argument('--test-monitoring', action='store_true',
                       help='Test monitoring functionality only')
    parser.add_argument('--validate-setup', action='store_true',
                       help='Validate environment setup only')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize tester
    tester = MiniIntegrationTester(args.config)
    
    if args.validate_setup:
        results = tester.validate_environment()
        return
    
    if args.test_import:
        results = tester.test_import_functionality()
        return
    
    if args.test_comparison:
        results = tester.test_comparison_functionality()
        return
    
    if args.test_assessment:
        results = tester.test_assessment_functionality()
        return
    
    if args.test_monitoring:
        results = tester.test_monitoring_functionality()
        return
    
    if args.full_test:
        results = tester.run_full_test_suite()
        return
    
    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
