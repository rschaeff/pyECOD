#!/usr/bin/env python3
"""
Mini PyECOD Production Workflow Guide

Step-by-step workflow for importing high-quality mini results and
propagating classifications to redundant sequences.

Usage:
    python mini_production_workflow.py --status                    # Check current status
    python mini_production_workflow.py --run-workflow             # Run complete workflow
    python mini_production_workflow.py --step quality-assessment  # Run specific step
    python mini_production_workflow.py --validate-results         # Validate final results
"""

import os
import sys
import argparse
import yaml
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MiniProductionWorkflow:
    """Orchestrates the complete mini production workflow"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        self.config_path = config_path
        self.script_dir = Path(__file__).parent
        
    def check_prerequisites(self) -> Dict[str, bool]:
        """Check if all required components are available"""
        
        logger.info("🔍 Checking prerequisites...")
        
        checks = {}
        
        # Check required scripts
        required_scripts = [
            'assess_quality_v2.py',
            'detect_collisions.py', 
            'import_results.py',
            'quality_filtered_importer.py'
        ]
        
        for script in required_scripts:
            script_path = self.script_dir / script
            checks[f"script_{script}"] = script_path.exists()
            if checks[f"script_{script}"]:
                logger.info(f"✓ Found {script}")
            else:
                logger.error(f"✗ Missing {script}")
        
        # Check config file
        config_path = Path(self.config_path)
        checks["config_file"] = config_path.exists()
        if checks["config_file"]:
            logger.info(f"✓ Found config: {config_path}")
        else:
            logger.error(f"✗ Missing config: {config_path}")
        
        # Check for mini results
        try:
            if checks["config_file"]:
                with open(config_path, 'r') as f:
                    config = yaml.safe_load(f)
                
                batch_base = Path(config["paths"]["batch_base_dir"])
                mini_results_found = False
                
                for batch_dir in batch_base.iterdir():
                    if (batch_dir / "mini_domains").exists():
                        xml_files = list((batch_dir / "mini_domains").glob("*.mini.domains.xml"))
                        if xml_files:
                            mini_results_found = True
                            break
                
                checks["mini_results"] = mini_results_found
                if mini_results_found:
                    logger.info("✓ Found mini domain results")
                else:
                    logger.warning("⚠️  No mini domain results found")
            else:
                checks["mini_results"] = False
        except Exception as e:
            logger.error(f"Error checking mini results: {e}")
            checks["mini_results"] = False
        
        # Check database connectivity
        try:
            result = subprocess.run([
                'python', str(self.script_dir / 'detect_collisions.py'), '--summary', '--config', self.config_path
            ], capture_output=True, text=True, timeout=30)
            
            checks["database_connection"] = result.returncode == 0
            if checks["database_connection"]:
                logger.info("✓ Database connection working")
            else:
                logger.error("✗ Database connection failed")
        except Exception as e:
            logger.error(f"✗ Database connection test failed: {e}")
            checks["database_connection"] = False
        
        return checks
    
    def run_quality_assessment(self) -> Dict[str, any]:
        """Step 1: Assess quality of all mini results"""
        
        logger.info("📊 Step 1: Quality Assessment")
        
        cmd = [
            'python', str(self.script_dir / 'assess_quality_v2.py'),
            '--scan-all',
            '--config', self.config_path
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info("✓ Quality assessment completed")
            # Parse key metrics from output
            lines = result.stdout.split('\n')
            excellent_count = 0
            good_count = 0
            total_count = 0
            
            for line in lines:
                if "Excellent" in line and "(" in line:
                    try:
                        excellent_count = int(line.split()[1].replace(',', ''))
                    except:
                        pass
                elif "Good" in line and "(" in line:
                    try:
                        good_count = int(line.split()[1].replace(',', ''))
                    except:
                        pass
                elif "Total results:" in line:
                    try:
                        total_count = int(line.split()[-1].replace(',', ''))
                    except:
                        pass
            
            return {
                'success': True,
                'excellent_count': excellent_count,
                'good_count': good_count,
                'total_count': total_count,
                'production_ready': excellent_count + good_count,
                'output': result.stdout
            }
        else:
            logger.error("✗ Quality assessment failed")
            logger.error(result.stderr)
            return {'success': False, 'error': result.stderr}
    
    def run_collision_detection(self) -> Dict[str, any]:
        """Step 2: Check for potential database collisions"""
        
        logger.info("🚨 Step 2: Collision Detection")
        
        cmd = [
            'python', str(self.script_dir / 'detect_collisions.py'),
            '--config', self.config_path
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info("✓ Collision detection completed")
            
            # Parse collision info
            lines = result.stdout.split('\n')
            collision_count = 0
            mini_proteins = 0
            
            for line in lines:
                if "Mini proteins to import:" in line:
                    try:
                        mini_proteins = int(line.split()[-1].replace(',', ''))
                    except:
                        pass
                elif "Potential collisions:" in line:
                    try:
                        collision_count = int(line.split()[-1].replace(',', ''))
                    except:
                        pass
            
            return {
                'success': True,
                'mini_proteins': mini_proteins,
                'collision_count': collision_count,
                'collision_rate': collision_count / max(1, mini_proteins) * 100,
                'output': result.stdout
            }
        else:
            logger.error("✗ Collision detection failed")
            logger.error(result.stderr)
            return {'success': False, 'error': result.stderr}
    
    def run_quality_filtered_import(self, limit: Optional[int] = None) -> Dict[str, any]:
        """Step 3: Import high-quality results with collision avoidance"""
        
        logger.info("🚀 Step 3: Quality-Filtered Import")
        
        cmd = [
            'python', str(self.script_dir / 'quality_filtered_importer.py'),
            '--assess-and-import',
            '--tier-filter', 'excellent,good',
            '--collision-strategy', 'separate',
            '--config', self.config_path
        ]
        
        if limit:
            cmd.extend(['--limit', str(limit)])
        
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info("✓ Quality-filtered import completed")
            
            # Parse import results
            lines = result.stdout.split('\n')
            imported_count = 0
            failed_count = 0
            
            for line in lines:
                if "Successfully imported:" in line:
                    try:
                        imported_count = int(line.split()[-1])
                    except:
                        pass
                elif "Failed:" in line:
                    try:
                        failed_count = int(line.split()[-1])
                    except:
                        pass
            
            return {
                'success': True,
                'imported_count': imported_count,
                'failed_count': failed_count,
                'output': result.stdout
            }
        else:
            logger.error("✗ Quality-filtered import failed")
            logger.error(result.stderr)
            return {'success': False, 'error': result.stderr}
    
    def run_sequence_propagation(self, dry_run: bool = False) -> Dict[str, any]:
        """Step 4: Propagate classifications to sequence-identical chains"""
        
        action = "DRY RUN" if dry_run else "EXECUTION"
        logger.info(f"🔄 Step 4: Sequence Propagation ({action})")
        
        cmd = [
            'python', str(self.script_dir / 'quality_filtered_importer.py'),
            '--propagate-sequences',
            '--limit-per-sequence', '10',
            '--config', self.config_path
        ]
        
        if dry_run:
            cmd.append('--dry-run')
        
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info(f"✓ Sequence propagation completed ({action})")
            
            # Parse propagation results
            lines = result.stdout.split('\n')
            candidates_count = 0
            successful_count = 0
            domains_propagated = 0
            
            for line in lines:
                if "Candidates processed:" in line:
                    try:
                        candidates_count = int(line.split()[-1])
                    except:
                        pass
                elif "Successful:" in line:
                    try:
                        successful_count = int(line.split()[-1])
                    except:
                        pass
                elif "Total domains propagated:" in line:
                    try:
                        domains_propagated = int(line.split()[-1])
                    except:
                        pass
            
            return {
                'success': True,
                'candidates_count': candidates_count,
                'successful_count': successful_count,
                'domains_propagated': domains_propagated,
                'dry_run': dry_run,
                'output': result.stdout
            }
        else:
            logger.error(f"✗ Sequence propagation failed ({action})")
            logger.error(result.stderr)
            return {'success': False, 'error': result.stderr}
    
    def run_comparative_analysis(self) -> Dict[str, any]:
        """Step 5: Generate comparative analysis"""
        
        logger.info("📊 Step 5: Comparative Analysis")
        
        cmd = [
            'python', str(self.script_dir / 'quality_filtered_importer.py'),
            '--comparative-analysis',
            '--config', self.config_path
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info("✓ Comparative analysis completed")
            return {
                'success': True,
                'output': result.stdout
            }
        else:
            logger.error("✗ Comparative analysis failed")
            logger.error(result.stderr)
            return {'success': False, 'error': result.stderr}
    
    def run_complete_workflow(self, limit: Optional[int] = None, 
                            propagation_dry_run: bool = False) -> Dict[str, any]:
        """Run the complete workflow"""
        
        logger.info("🚀 Starting Complete Mini Production Workflow")
        logger.info(f"Timestamp: {datetime.now().isoformat()}")
        
        workflow_results = {}
        
        # Step 1: Quality Assessment
        step1 = self.run_quality_assessment()
        workflow_results['quality_assessment'] = step1
        
        if not step1['success']:
            logger.error("❌ Workflow stopped - Quality assessment failed")
            return workflow_results
        
        production_ready = step1.get('production_ready', 0)
        if production_ready == 0:
            logger.warning("⚠️  No production-ready results found")
            return workflow_results
        
        logger.info(f"📊 Quality Assessment: {production_ready} production-ready results")
        
        # Step 2: Collision Detection
        step2 = self.run_collision_detection()
        workflow_results['collision_detection'] = step2
        
        if not step2['success']:
            logger.error("❌ Workflow stopped - Collision detection failed")
            return workflow_results
        
        collision_rate = step2.get('collision_rate', 0)
        logger.info(f"🚨 Collision Detection: {collision_rate:.1f}% collision rate")
        
        # Step 3: Quality-Filtered Import
        step3 = self.run_quality_filtered_import(limit)
        workflow_results['quality_import'] = step3
        
        if not step3['success']:
            logger.error("❌ Workflow stopped - Import failed")
            return workflow_results
        
        imported = step3.get('imported_count', 0)
        logger.info(f"🚀 Quality Import: {imported} proteins imported successfully")
        
        if imported == 0:
            logger.warning("⚠️  No proteins imported - skipping propagation")
            return workflow_results
        
        # Step 4: Sequence Propagation (with dry run option)
        step4 = self.run_sequence_propagation(propagation_dry_run)
        workflow_results['sequence_propagation'] = step4
        
        if not step4['success']:
            logger.error("❌ Sequence propagation failed")
            return workflow_results
        
        propagated = step4.get('successful_count', 0)
        action = "analyzed" if propagation_dry_run else "propagated"
        logger.info(f"🔄 Sequence Propagation: {propagated} chains {action}")
        
        # Step 5: Comparative Analysis
        step5 = self.run_comparative_analysis()
        workflow_results['comparative_analysis'] = step5
        
        if step5['success']:
            logger.info("📊 Comparative Analysis: Report generated")
        
        # Summary
        workflow_results['summary'] = {
            'timestamp': datetime.now().isoformat(),
            'success': True,
            'production_ready_count': production_ready,
            'imported_count': imported,
            'propagated_count': propagated,
            'propagation_dry_run': propagation_dry_run
        }
        
        logger.info("🎉 Complete workflow finished successfully!")
        
        return workflow_results
    
    def print_workflow_status(self):
        """Print current workflow status"""
        
        print("\n🔍 Mini PyECOD Production Workflow Status")
        print("=" * 60)
        print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print()
        
        # Check prerequisites
        checks = self.check_prerequisites()
        
        print("📋 Prerequisites:")
        all_good = True
        for check, passed in checks.items():
            status = "✓" if passed else "✗"
            print(f"  {status} {check.replace('_', ' ').title()}")
            if not passed:
                all_good = False
        
        if not all_good:
            print("\n❌ Prerequisites not met - fix issues before proceeding")
            return
        
        print("\n✅ All prerequisites met - ready to run workflow")
        
        # Quick quality check
        try:
            result = subprocess.run([
                'python', str(self.script_dir / 'assess_quality_v2.py'),
                '--scan-all', '--config', self.config_path
            ], capture_output=True, text=True, timeout=60)
            
            if result.returncode == 0:
                lines = result.stdout.split('\n')
                for line in lines:
                    if "Production ready:" in line:
                        print(f"\n📊 Quick Status: {line}")
                        break
            
        except Exception as e:
            logger.warning(f"Could not get quick status: {e}")
        
        print("\n🚀 Next Steps:")
        print("  1. Run complete workflow:")
        print("     python mini_production_workflow.py --run-workflow")
        print()
        print("  2. Or run individual steps:")
        print("     python mini_production_workflow.py --step quality-assessment")
        print("     python mini_production_workflow.py --step collision-detection")
        print("     python mini_production_workflow.py --step quality-import")
        print("     python mini_production_workflow.py --step sequence-propagation")
        print("     python mini_production_workflow.py --step comparative-analysis")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Mini PyECOD Production Workflow Orchestrator'
    )
    
    parser.add_argument('--status', action='store_true',
                       help='Check current workflow status')
    parser.add_argument('--run-workflow', action='store_true',
                       help='Run complete workflow')
    parser.add_argument('--step', type=str,
                       choices=['quality-assessment', 'collision-detection', 
                              'quality-import', 'sequence-propagation', 'comparative-analysis'],
                       help='Run specific workflow step')
    parser.add_argument('--limit', type=int,
                       help='Limit number of results to import (for testing)')
    parser.add_argument('--propagation-dry-run', action='store_true',
                       help='Run propagation in dry-run mode (analyze only)')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize workflow
    workflow = MiniProductionWorkflow(args.config)
    
    if args.status:
        workflow.print_workflow_status()
        return
    
    if args.run_workflow:
        results = workflow.run_complete_workflow(args.limit, args.propagation_dry_run)
        
        if results.get('summary', {}).get('success'):
            summary = results['summary']
            print(f"\n🎉 Workflow Complete!")
            print(f"  Production ready: {summary['production_ready_count']}")
            print(f"  Imported: {summary['imported_count']}")
            action = "analyzed" if summary['propagation_dry_run'] else "propagated"
            print(f"  Propagated: {summary['propagated_count']} ({action})")
        else:
            print(f"\n❌ Workflow failed - check logs for details")
        return
    
    if args.step:
        if args.step == 'quality-assessment':
            result = workflow.run_quality_assessment()
        elif args.step == 'collision-detection':
            result = workflow.run_collision_detection()
        elif args.step == 'quality-import':
            result = workflow.run_quality_filtered_import(args.limit)
        elif args.step == 'sequence-propagation':
            result = workflow.run_sequence_propagation(args.propagation_dry_run)
        elif args.step == 'comparative-analysis':
            result = workflow.run_comparative_analysis()
        
        if result['success']:
            print(f"✅ Step '{args.step}' completed successfully")
            print(result['output'])
        else:
            print(f"❌ Step '{args.step}' failed")
            print(result.get('error', 'Unknown error'))
        return
    
    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
