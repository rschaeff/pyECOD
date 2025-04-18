#!/usr/bin/env python3
"""
submit_repopulate_jobs.py - Submit SLURM jobs for repopulating ECOD schema

This script identifies batches in the 'created' state and submits SLURM jobs
to run the repopulate_ecod_schema.py script for each batch.
"""

import os
import sys
import argparse
import logging
import subprocess
from typing import List, Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.core.slurm_job_manager import SlurmJobManager

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )

def get_created_batches(context: ApplicationContext) -> List[Dict[str, Any]]:
    """
    Get list of batches in 'created' state
    
    Args:
        context: Application context
        
    Returns:
        List of batch information dictionaries
    """
    query = """
    SELECT id, batch_name, base_path, type, ref_version, total_items
    FROM ecod_schema.batch
    WHERE status = 'created'
    ORDER BY batch_name
    """
    
    result = context.db.execute_query(query)
    
    batches = []
    for row in result:
        if isinstance(row, dict):
            batches.append(row)
        elif isinstance(row, tuple):
            # Convert tuple to dict
            batches.append({
                'id': row[0],
                'batch_name': row[1],
                'base_path': row[2],
                'type': row[3],
                'ref_version': row[4],
                'total_items': row[5]
            })
    
    return batches

def submit_repopulate_job(job_manager: SlurmJobManager, batch_info: Dict[str, Any], 
                         config_path: str, script_path: str, threads: int, memory: str, time: str, 
                         dry_run: bool = False) -> Optional[str]:
    """
    Submit a SLURM job for repopulating a batch
    
    Args:
        job_manager: SlurmJobManager instance
        batch_info: Batch information dictionary
        config_path: Path to configuration file
        script_path: Path to repopulate_ecod_schema.py script
        threads: Number of threads to allocate
        memory: Memory to allocate (e.g. "8G")
        time: Time limit (e.g. "24:00:00")
        dry_run: If True, don't actually submit the job
        
    Returns:
        Job ID if successful, None otherwise
    """
    batch_name = batch_info['batch_name']
    
    # Create commands for the job
    commands = [
        f"cd {os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}", # Go to project root
        f"python {script_path} --config {config_path} --batch-name {batch_name}"
    ]
    
    # Create job script
    job_name = f"repop_{batch_name.replace('ecod_batch_', '')}"
    output_dir = f"/tmp/ecod_jobs/{batch_info['id']}"
    
    os.makedirs(output_dir, exist_ok=True)
    
    script_path = job_manager.create_job_script(
        commands=commands,
        job_name=job_name,
        output_dir=output_dir,
        threads=threads,
        memory=memory,
        time=time
    )
    
    if dry_run:
        logging.info(f"DRY RUN: Would submit job for batch {batch_name} (ID: {batch_info['id']})")
        logging.info(f"Script path: {script_path}")
        return None
    
    # Submit the job
    try:
        job_id = job_manager.submit_job(script_path)
        logging.info(f"Submitted job {job_id} for batch {batch_name} (ID: {batch_info['id']})")
        return job_id
    except Exception as e:
        logging.error(f"Failed to submit job for batch {batch_name} (ID: {batch_info['id']}): {str(e)}")
        return None

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Submit SLURM jobs for repopulating ECOD schema')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--script', type=str, default='scripts/core/repopulate_ecod_schema.py',
                      help='Path to repopulate_ecod_schema.py script')
    parser.add_argument('--threads', type=int, default=4,
                      help='Number of threads to allocate')
    parser.add_argument('--memory', type=str, default='8G',
                      help='Memory to allocate (e.g. "8G")')
    parser.add_argument('--time', type=str, default='4:00:00',
                      help='Time limit (e.g. "4:00:00")')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show what would be done without making changes')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.submit_repopulate_jobs")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Initialize job manager
    job_manager = SlurmJobManager(args.config)
    
    # Get batches in 'created' state
    logger.info("Retrieving batches in 'created' state")
    batches = get_created_batches(context)
    
    if not batches:
        logger.info("No batches in 'created' state found")
        return 0
    
    logger.info(f"Found {len(batches)} batches in 'created' state")
    
    # Submit jobs for each batch
    submitted_jobs = 0
    failed_jobs = 0
    
    for batch_info in batches:
        job_id = submit_repopulate_job(
            job_manager, 
            batch_info, 
            args.config,
            args.script,
            args.threads, 
            args.memory, 
            args.time, 
            args.dry_run
        )
        
        if job_id:
            submitted_jobs += 1
        else:
            failed_jobs += 1
    
    logger.info(f"Submitted {submitted_jobs} jobs, {failed_jobs} failed")
    
    if args.dry_run:
        logger.info("This was a dry run. No actual jobs were submitted.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())