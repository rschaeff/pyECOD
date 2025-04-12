#!/usr/bin/env python3
"""
submit_summaries_to_slurm.py - Submit domain summary generation jobs to SLURM cluster
"""

import os
import sys
import time
import logging
import argparse
import subprocess
from datetime import datetime
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

def get_indexed_batches(context: ApplicationContext) -> List[Dict[str, Any]]:
    """Get list of indexed batches directly from the database"""
    logger = logging.getLogger("ecod.submit_summaries")
    
    query = """
    SELECT 
        id, batch_name, base_path, ref_version, total_items
    FROM 
        ecod_schema.batch
    WHERE 
        total_items = 5000 AND status = 'indexed'
    ORDER BY 
        id
    """
    
    result = context.db.execute_query(query)
    
    batches = []
    for row in result:
        batches.append({
            'id': row[0],
            'name': row[1],
            'path': row[2],
            'reference': row[3],
            'total_items': row[4]
        })
    
    logger.info(f"Found {len(batches)} indexed batches")
    
    return batches

def check_batch_summaries(batch_id: int, context: ApplicationContext) -> Dict[str, int]:
    """Check how many domain summaries exist for this batch"""
    query = """
    SELECT 
        COUNT(DISTINCT pf.id) as summary_count
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    LEFT JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'domain_summary' AND pf.file_exists = TRUE
    WHERE 
        b.id = %s
    """
    
    result = context.db.execute_query(query, (batch_id,))
    
    summary_count = result[0][0] if result else 0
    
    return {
        'batch_id': batch_id,
        'summary_count': summary_count
    }

def submit_domain_summary_job(batch_id: int, job_manager: SlurmJobManager, 
                           config_path: str, output_dir: Optional[str] = None,
                           blast_only: bool = True, threads: int = 8, 
                           memory: str = "16G", time: str = "8:00:00",
                           force: bool = False, limit: Optional[int] = None) -> str:
    """Submit a job to generate domain summaries for a batch"""
    logger = logging.getLogger("ecod.submit_summaries")
    
    # Build command for the job
    cmd = [
        "python", "scripts/generate_batch_domain_summaries.py",
        "--config", config_path,
        "--batch-id", str(batch_id),
        "--threads", str(threads)
    ]
    
    if blast_only:
        cmd.append("--blast-only")
    
    if output_dir:
        cmd.extend(["--output-dir", output_dir])
    
    if force:
        cmd.append("--force")
    
    if limit:
        cmd.extend(["--limit", str(limit)])
    
    # Create log directory
    log_dir = os.path.join("logs", "domain_summaries", f"batch_{batch_id}")
    os.makedirs(log_dir, exist_ok=True)
    
    # Add log file to command
    log_file = os.path.join(log_dir, f"domain_summary.log")
    cmd.extend(["--log-file", log_file])
    
    # Create a job script
    job_name = f"summary_batch_{batch_id}"
    script_path = job_manager.create_job_script(
        commands=[" ".join(cmd)],
        job_name=job_name,
        output_dir=log_dir,
        threads=threads,
        memory=memory,
        time=time
    )
    
    # Submit the job
    job_id = job_manager.submit_job(script_path)
    logger.info(f"Submitted job {job_id} for batch {batch_id}")
    
    # Update job tracking in database
    try:
        query = """
        INSERT INTO ecod_schema.job_tracking
        (job_id, batch_id, job_type, status, submitted_at)
        VALUES (%s, %s, %s, %s, NOW())
        """
        job_manager.context.db.execute_query(query, (job_id, batch_id, "domain_summary", "submitted"))
    except Exception as e:
        logger.warning(f"Failed to record job in database: {e}")
    
    return job_id

def main():
    """Main function to submit domain summary jobs to SLURM"""
    parser = argparse.ArgumentParser(description='Submit domain summary generation jobs to SLURM')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--output-dir', type=str,
                      help='Override output directory (default: from batch path)')
    parser.add_argument('--batch-id', type=int,
                      help='Process specific batch ID only')
    parser.add_argument('--start-id', type=int,
                      help='Start processing from this batch ID')
    parser.add_argument('--end-id', type=int,
                      help='Stop processing at this batch ID')
    parser.add_argument('--blast-only', action='store_true', default=True,
                      help='Generate blast-only summaries (no HHSearch)')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads per job')
    parser.add_argument('--memory', type=str, default="16G",
                      help='Memory request per job (e.g., 16G)')
    parser.add_argument('--time', type=str, default="8:00:00",
                      help='Time limit per job (e.g., 8:00:00)')
    parser.add_argument('--max-parallel', type=int, default=5,
                      help='Maximum number of parallel jobs to submit')
    parser.add_argument('--sleep', type=int, default=5,
                      help='Seconds to sleep between job submissions')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show batches but do not submit jobs')
    parser.add_argument('--force', action='store_true',
                      help='Force regeneration of existing summaries')
    parser.add_argument('--limit', type=int,
                      help='Limit number of proteins per batch to process')
    parser.add_argument('--skip-complete', action='store_true',
                      help='Skip batches that already have all summaries')
    parser.add_argument('--check-interval', type=int, default=60,
                      help='Interval (seconds) to check job status')
    parser.add_argument('--wait', action='store_true',
                      help='Wait for all jobs to complete')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    log_dir = os.path.join("logs", "domain_summaries")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f'submit_jobs_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    setup_logging(args.verbose, log_file)
    logger = logging.getLogger("ecod.submit_summaries")
    
    # Initialize application context and job manager
    context = ApplicationContext(args.config)
    job_manager = SlurmJobManager(args.config)
    
    # Get indexed batches
    if args.batch_id:
        batches = [{'id': args.batch_id}]
        logger.info(f"Processing single batch ID: {args.batch_id}")
    else:
        batches = get_indexed_batches(context)
        
        # Filter by start/end ID if specified
        if args.start_id:
            batches = [b for b in batches if b['id'] >= args.start_id]
        if args.end_id:
            batches = [b for b in batches if b['id'] <= args.end_id]
        
    if not batches:
        logger.warning(f"No indexed batches found to process")
        return 0
    
    # Check current summary status for each batch
    for i, batch in enumerate(batches):
        batch_id = batch['id']
        summary_status = check_batch_summaries(batch_id, context)
        batch['summary_count'] = summary_status['summary_count']
        
        # Skip batches that already have all summaries
        if args.skip_complete and batch.get('total_items', 5000) == batch['summary_count']:
            logger.info(f"Batch {batch_id} already has all summaries ({batch['summary_count']}/{batch.get('total_items', 5000)})")
            batch['skip'] = True
        else:
            batch['skip'] = False
    
    # Display batch info
    logger.info(f"Found {len(batches)} batches to process:")
    for i, batch in enumerate(batches):
        batch_id = batch['id']
        if 'name' in batch:
            skip_status = " [SKIP - COMPLETE]" if batch.get('skip', False) else ""
            logger.info(f"  {i+1}. Batch ID: {batch_id} - {batch['name']}{skip_status}")
            logger.info(f"     Summary files: {batch['summary_count']}/{batch.get('total_items', 5000)}")
        else:
            logger.info(f"  {i+1}. Batch ID: {batch_id}")
    
    if args.dry_run:
        logger.info("Dry run completed. No jobs submitted.")
        return 0
    
    # Submit jobs for each batch
    submitted_jobs = []
    active_jobs = []
    
    for batch in batches:
        batch_id = batch['id']
        
        # Skip if marked for skipping
        if batch.get('skip', False):
            logger.info(f"Skipping batch {batch_id} (already complete)")
            continue
        
        # Check if we need to wait for active jobs to finish
        while len(active_jobs) >= args.max_parallel:
            logger.info(f"Waiting for job slots ({len(active_jobs)}/{args.max_parallel} active)...")
            
            # Check active jobs
            still_active = []
            for job_id in active_jobs:
                status = job_manager.check_job_status(job_id)
                if status and status['state'] in ['PENDING', 'RUNNING', 'CONFIGURING']:
                    still_active.append(job_id)
            
            active_jobs = still_active
            
            if len(active_jobs) >= args.max_parallel:
                # Still waiting for slots
                time.sleep(args.check_interval)
            
        # Submit job
        job_id = submit_domain_summary_job(
            batch_id=batch_id,
            job_manager=job_manager,
            config_path=args.config,
            output_dir=args.output_dir,
            blast_only=args.blast_only,
            threads=args.threads,
            memory=args.memory,
            time=args.time,
            force=args.force,
            limit=args.limit
        )
        
        submitted_jobs.append(job_id)
        active_jobs.append(job_id)
        
        # Sleep between submissions
        if args.sleep > 0:
            time.sleep(args.sleep)
    
    logger.info(f"Submitted {len(submitted_jobs)} jobs")
    
    # Wait for all jobs to complete if requested
    if args.wait and submitted_jobs:
        logger.info("Waiting for all jobs to complete...")
        
        while True:
            active_count = 0
            completed_count = 0
            failed_count = 0
            
            for job_id in submitted_jobs:
                status = job_manager.check_job_status(job_id)
                if status:
                    if status['state'] in ['PENDING', 'RUNNING', 'CONFIGURING']:
                        active_count += 1
                    elif status['state'] == 'COMPLETED':
                        completed_count += 1
                    else:
                        failed_count += 1
            
            logger.info(f"Job status: {active_count} active, {completed_count} completed, {failed_count} failed")
            
            if active_count == 0:
                break
                
            time.sleep(args.check_interval)
        
        logger.info(f"All jobs completed: {completed_count} successful, {failed_count} failed")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())