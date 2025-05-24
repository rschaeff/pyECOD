#!/usr/bin/env python3
"""
Integrated Batch Repair Workflow

This script provides a complete workflow for identifying and repairing problematic batches.
It combines problematic batch identification with domain partition repair and import processes.

Usage:
    python repair_batch_workflow.py --config config.yml [options]

Workflow Steps:
    1. Identify problematic batches
    2. Run import script in repair mode to capture existing partition files
    3. Identify remaining missing partitions
    4. Optionally trigger domain partition pipeline for missing proteins
    5. Re-import newly generated partitions

Options:
    --config CONFIG          Path to configuration file
    --batch-id BATCH_ID      Process only specific batch (skip identification)
    --dry-run                Show what would be done without executing
    --severity-threshold     Minimum severity for auto-repair: medium, high, critical
    --skip-import            Skip import step (only identify issues)
    --skip-partition         Skip domain partition step
    --include-non-rep        Include non-representative proteins
    --max-batches N          Limit number of batches to repair
    --verbose                Enable verbose output
"""

import os
import sys
import logging
import argparse
import subprocess
import json
import yaml
from datetime import datetime
from typing import List, Dict, Any, Optional
from dataclasses import dataclass


@dataclass
class RepairAction:
    """Represents a repair action to be taken."""
    batch_id: int
    batch_name: str
    action_type: str
    command: str
    description: str
    required: bool = True


def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format_str)
    return logging.getLogger(__name__)


def run_command(command: str, description: str, dry_run: bool = False, logger=None) -> bool:
    """Run a command with proper logging and error handling."""
    if not logger:
        logger = logging.getLogger(__name__)
    
    logger.info(f"{'[DRY RUN] ' if dry_run else ''}{description}")
    logger.debug(f"Command: {command}")
    
    if dry_run:
        return True
    
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            capture_output=True, 
            text=True,
            timeout=3600  # 1 hour timeout
        )
        
        if result.returncode == 0:
            logger.info(f"✅ {description} - completed successfully")
            if result.stdout.strip():
                logger.debug(f"Output: {result.stdout.strip()}")
            return True
        else:
            logger.error(f"❌ {description} - failed with return code {result.returncode}")
            if result.stderr.strip():
                logger.error(f"Error: {result.stderr.strip()}")
            if result.stdout.strip():
                logger.debug(f"Output: {result.stdout.strip()}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error(f"❌ {description} - timed out after 1 hour")
        return False
    except Exception as e:
        logger.error(f"❌ {description} - exception: {str(e)}")
        return False


def identify_problematic_batches(config_path: str, severity: str, include_non_rep: bool, 
                               max_batches: Optional[int] = None, logger=None) -> List[Dict[str, Any]]:
    """Identify problematic batches using the identification script."""
    if not logger:
        logger = logging.getLogger(__name__)
    
    logger.info("Identifying problematic batches...")
    
    # Build command
    script_dir = os.path.dirname(os.path.abspath(__file__))
    identify_script = os.path.join(script_dir, "identify_problematic_batches.py")
    
    command = f"python {identify_script} --config {config_path} --output-format json --severity-threshold {severity}"
    
    if include_non_rep:
        command += " --include-non-rep"
    
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            logger.error(f"Failed to identify problematic batches: {result.stderr}")
            return []
        
        # Parse JSON output
        data = json.loads(result.stdout)
        issues = data.get('issues', [])
        
        # Group by batch and filter auto-fixable
        batch_issues = {}
        for issue in issues:
            if issue.get('auto_fixable', False):
                batch_id = issue['batch_id']
                if batch_id not in batch_issues:
                    batch_issues[batch_id] = {
                        'batch_id': batch_id,
                        'batch_name': issue['batch_name'],
                        'issues': []
                    }
                batch_issues[batch_id]['issues'].append(issue)
        
        # Sort by severity and limit if requested
        problematic_batches = list(batch_issues.values())
        problematic_batches.sort(key=lambda x: len([i for i in x['issues'] if i['severity'] == 'critical']), reverse=True)
        
        if max_batches:
            problematic_batches = problematic_batches[:max_batches]
        
        logger.info(f"Found {len(problematic_batches)} problematic batches with auto-fixable issues")
        
        return problematic_batches
        
    except Exception as e:
        logger.error(f"Error identifying problematic batches: {str(e)}")
        return []


def create_repair_plan(problematic_batches: List[Dict[str, Any]], config_path: str, 
                      skip_import: bool, skip_partition: bool, logger=None) -> List[RepairAction]:
    """Create a repair plan for the problematic batches."""
    if not logger:
        logger = logging.getLogger(__name__)
    
    repair_actions = []
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    for batch_data in problematic_batches:
        batch_id = batch_data['batch_id']
        batch_name = batch_data['batch_name']
        issues = batch_data['issues']
        
        logger.info(f"Creating repair plan for batch {batch_id}: {batch_name}")
        
        # Analyze issues to determine required actions
        has_missing_partitions = any(i['issue_type'] == 'missing_partitions' for i in issues)
        has_partition_gap = any(i['issue_type'] == 'partition_gap' for i in issues)
        
        # Step 1: Import existing partition files in repair mode
        if not skip_import and (has_missing_partitions or has_partition_gap):
            import_script = os.path.join(script_dir, "import_domain_partitions_with_repair.py")
            import_command = (f"python {import_script} --config {config_path} "
                            f"--batch-id {batch_id} --repair-mode --conflict-strategy update "
                            f"--validate-before-import")
            
            repair_actions.append(RepairAction(
                batch_id=batch_id,
                batch_name=batch_name,
                action_type="import_repair",
                command=import_command,
                description=f"Import existing partition files for batch {batch_id}",
                required=True
            ))
        
        # Step 2: Re-analyze to see what's still missing after import
        # This would require running the identification script again, but for now we'll assume
        # that if there were missing partitions, we might need to regenerate some
        
        # Step 3: Generate missing partitions if needed
        if not skip_partition and has_missing_partitions:
            partition_script = os.path.join(script_dir, "domain_partition_run.py")
            
            # First try repair of failed processes
            repair_command = (f"python {partition_script} --config {config_path} "
                            f"repair failed --batch-id {batch_id} --rerun")
            
            repair_actions.append(RepairAction(
                batch_id=batch_id,
                batch_name=batch_name,
                action_type="repair_failed",
                command=repair_command,
                description=f"Repair failed processes for batch {batch_id}",
                required=False  # This might not be needed for all batches
            ))
            
            # Then run partition for any remaining missing files
            partition_command = (f"python {partition_script} --config {config_path} "
                               f"process batch --batch-id {batch_id} --force")
            
            repair_actions.append(RepairAction(
                batch_id=batch_id,
                batch_name=batch_name,
                action_type="generate_partitions",
                command=partition_command,
                description=f"Generate missing domain partitions for batch {batch_id}",
                required=False  # Only if import doesn't solve the issue
            ))
        
        # Step 4: Final import of newly generated partitions
        if not skip_import and not skip_partition and has_missing_partitions:
            final_import_script = os.path.join(script_dir, "import_domain_partitions_with_repair.py")
            final_import_command = (f"python {final_import_script} --config {config_path} "
                                  f"--batch-id {batch_id} --conflict-strategy skip")
            
            repair_actions.append(RepairAction(
                batch_id=batch_id,
                batch_name=batch_name,
                action_type="final_import",
                command=final_import_command,
                description=f"Import newly generated partitions for batch {batch_id}",
                required=False
            ))
    
    return repair_actions


def execute_repair_plan(repair_actions: List[RepairAction], dry_run: bool, logger=None) -> Dict[str, int]:
    """Execute the repair plan."""
    if not logger:
        logger = logging.getLogger(__name__)
    
    stats = {
        'total_actions': len(repair_actions),
        'successful_actions': 0,
        'failed_actions': 0,
        'skipped_actions': 0
    }
    
    # Group actions by batch for better organization
    batch_actions = {}
    for action in repair_actions:
        if action.batch_id not in batch_actions:
            batch_actions[action.batch_id] = []
        batch_actions[action.batch_id].append(action)
    
    logger.info(f"Executing repair plan for {len(batch_actions)} batches...")
    
    for batch_id in sorted(batch_actions.keys()):
        actions = batch_actions[batch_id]
        batch_name = actions[0].batch_name
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing Batch {batch_id}: {batch_name}")
        logger.info(f"{'='*60}")
        
        batch_success = True
        
        for action in actions:
            logger.info(f"\nStep: {action.action_type}")
            
            success = run_command(action.command, action.description, dry_run, logger)
            
            if success:
                stats['successful_actions'] += 1
            else:
                stats['failed_actions'] += 1
                if action.required:
                    logger.error(f"Required action failed for batch {batch_id}, skipping remaining actions")
                    batch_success = False
                    break
                else:
                    logger.warning(f"Optional action failed for batch {batch_id}, continuing")
        
        if not batch_success:
            # Mark remaining actions as skipped
            remaining_actions = actions[actions.index(action) + 1:]
            stats['skipped_actions'] += len(remaining_actions)
            logger.warning(f"Skipped {len(remaining_actions)} remaining actions for batch {batch_id}")
    
    return stats


def generate_repair_report(problematic_batches: List[Dict[str, Any]], 
                         repair_actions: List[RepairAction], 
                         execution_stats: Dict[str, int], 
                         logger=None):
    """Generate a summary report of the repair process."""
    if not logger:
        logger = logging.getLogger(__name__)
    
    logger.info(f"\n{'='*80}")
    logger.info("BATCH REPAIR WORKFLOW SUMMARY")
    logger.info(f"{'='*80}")
    
    # Batch analysis summary
    logger.info(f"Batches analyzed: {len(problematic_batches)}")
    
    issue_summary = {}
    for batch_data in problematic_batches:
        for issue in batch_data['issues']:
            issue_type = issue['issue_type']
            severity = issue['severity']
            key = f"{issue_type} ({severity})"
            issue_summary[key] = issue_summary.get(key, 0) + 1
    
    if issue_summary:
        logger.info("Issues found:")
        for issue_type, count in sorted(issue_summary.items()):
            logger.info(f"  {issue_type}: {count}")
    
    # Execution summary
    logger.info(f"\nRepair actions executed:")
    logger.info(f"  Total actions: {execution_stats['total_actions']}")
    logger.info(f"  Successful: {execution_stats['successful_actions']}")
    logger.info(f"  Failed: {execution_stats['failed_actions']}")
    logger.info(f"  Skipped: {execution_stats['skipped_actions']}")
    
    success_rate = (execution_stats['successful_actions'] / execution_stats['total_actions'] * 100 
                   if execution_stats['total_actions'] > 0 else 0)
    logger.info(f"  Success rate: {success_rate:.1f}%")
    
    # Recommendations
    logger.info(f"\nRecommendations:")
    if execution_stats['failed_actions'] > 0:
        logger.info("  - Review failed actions and address underlying issues")
        logger.info("  - Check log files for detailed error information")
    
    if execution_stats['successful_actions'] > 0:
        logger.info("  - Run audit again to verify repair effectiveness")
        logger.info("  - Monitor batch processing status")
    
    logger.info("  - Consider updating batch metadata if definition mismatches persist")
    logger.info(f"{'='*80}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Integrated batch repair workflow')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, help='Process only specific batch')
    parser.add_argument('--dry-run', action='store_true', help='Show actions without executing')
    parser.add_argument('--severity-threshold', choices=['medium', 'high', 'critical'], 
                       default='high', help='Minimum severity for auto-repair')
    parser.add_argument('--skip-import', action='store_true', help='Skip import steps')
    parser.add_argument('--skip-partition', action='store_true', help='Skip domain partition steps')
    parser.add_argument('--include-non-rep', action='store_true', 
                       help='Include non-representative proteins')
    parser.add_argument('--max-batches', type=int, default=10, 
                       help='Maximum number of batches to repair')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.verbose)

    logger.info(f"Starting batch repair workflow at {datetime.now()}")
    logger.info(f"Configuration: {args.config}")
    
    if args.dry_run:
        logger.info("DRY RUN MODE - no actual changes will be made")

    try:
        # Step 1: Identify problematic batches (unless specific batch provided)
        if args.batch_id:
            # Process specific batch
            problematic_batches = [{
                'batch_id': args.batch_id,
                'batch_name': f'Batch {args.batch_id}',
                'issues': [{'issue_type': 'manual_repair', 'severity': 'high', 'auto_fixable': True}]
            }]
            logger.info(f"Processing specific batch: {args.batch_id}")
        else:
            # Identify problematic batches
            problematic_batches = identify_problematic_batches(
                args.config, args.severity_threshold, args.include_non_rep, 
                args.max_batches, logger
            )
        
        if not problematic_batches:
            logger.info("No problematic batches found or specified. Exiting.")
            return 0
        
        # Step 2: Create repair plan
        logger.info("Creating repair plan...")
        repair_actions = create_repair_plan(
            problematic_batches, args.config, args.skip_import, args.skip_partition, logger
        )
        
        if not repair_actions:
            logger.info("No repair actions needed. Exiting.")
            return 0
        
        logger.info(f"Created repair plan with {len(repair_actions)} actions")
        
        # Step 3: Execute repair plan
        logger.info("Executing repair plan...")
        execution_stats = execute_repair_plan(repair_actions, args.dry_run, logger)
        
        # Step 4: Generate summary report
        generate_repair_report(problematic_batches, repair_actions, execution_stats, logger)
        
        # Determine exit code based on results
        if execution_stats['failed_actions'] > 0:
            logger.warning("Some repair actions failed")
            return 1
        else:
            logger.info("Batch repair workflow completed successfully")
            return 0

    except KeyboardInterrupt:
        logger.info("Repair workflow interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Error during repair workflow: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
