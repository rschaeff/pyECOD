#!/usr/bin/env python3
"""
run_domain_summaries_slurm.py - Submit domain summary jobs to SLURM for indexed batches
"""

import os
import sys
import logging
import argparse
from datetime import datetime
from typing import List, Dict, Any, Optional

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.jobs import SlurmJobManager

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
    """Get list of indexed batches from the database"""
    logger = logging.getLogger("ecod.slurm_summaries")
    
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

def check_summary_status(batch_id: int, context: ApplicationContext) -> Dict[str, int]:
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

def create_summary_job(batch_id: int, job_manager: SlurmJobManager, config_path: str, 
                     blast_only: bool = True, threads: int = 8, memory: str = "16G", 
                     time: str = "8:00:00", log_dir: str = "logs") -> Dict[str, Any]:
    """Create a job script for generating domain summaries"""
    logger = logging.getLogger("ecod.slurm_summaries")
    
    # Create batch-specific log directory
    batch_log_dir = os.path.join(log_dir, f"batch_{batch_id}")
    os.makedirs(batch_log_dir, exist_ok=True)
    
    # Build command
    command = f"python scripts/generate_batch_domain_summaries.py --config {config_path} --batch-id {batch_id} --threads {threads}"
    
    if blast_only:
        command += " --blast-only"
    
    # Add log file to command
    log_file = os.path.join(batch_log_dir, "domain_summary.log")
    command += f" --log-file {log_file}"
    
    # Job name
    job_name = f"summary_batch_{batch_id}"
    
    # Create job script
    script_path = job_manager.create_job_script(
        commands=[command],
        job_name=job_name,
        output_dir=batch_log_dir,
        threads=threads,
        memory=memory,
        time=time
    )
    
    return {
        'batch_id': batch_id,
        'script_path': script_path,
        'job_name': job_name,
        'log_dir': batch_log_dir
    }

def submit_job(job_info: Dict[str, Any], job_manager: SlurmJobManager, context: ApplicationContext) -> str:
    """Submit a job and record in database"""
    logger = logging.getLogger("ecod.slurm_summaries")
    
    # Submit the job
    job_id = job_manager.submit_job(job_info['script_path'])
    
    if not job_id:
        logger.error(f"Failed to submit job for batch {job_info['batch_id']}")
        return None
    
    logger.info(f"Submitted job {job_id} for batch {job_info['batch_id']}")
    
    # Record in database
    try:
        query = """
        INSERT INTO ecod_schema.job
        (batch_id, job_type, slurm_job_id, status, submission_time)
        VALUES (%s, %s, %s, %s, NOW())
        RETURNING id
        """
        
        result = context.db.execute_query(
            query, 
            (job_info['batch_id'], 'domain_summary', job_id, 'submitted')
        )
        
        db_job_id = result[0][0] if result else None
        logger.info(f"Recorded job in database with ID {db_job_id}")
        
    except Exception as e:
        logger.warning(f"Failed to record job in database: {e}")
    
    return job_id

def main():
    """Main function to submit domain summary jobs to SLURM"""
    parser = argparse.ArgumentParser(description='Submit domain summary jobs to SLURM for indexed batches')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--log-dir', type=str, default='logs/domain_summaries',
                      help='Directory for log files')
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
    parser.add_argument('--dry-run', action='store_true',
                      help='Show batches but do not submit jobs')
    parser.add_argument('--force', action='store_true',
                      help='Force regeneration of existing summaries')
    parser.add_argument('--skip-complete', action='store_true',
                      help='Skip batches that already have all summaries')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    os.makedirs(args.log_dir, exist_ok=True)
    log_file = os.path.join(args.log_dir, f'submit_jobs_{timestamp}.log')
    setup_logging(args.verbose, log_file)
    logger = logging.getLogger("ecod.slurm_summaries")
    
    # Initialize context and job manager
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
        logger.warning("No indexed batches found to process")
        return 1
    
    # Check current summary status for each batch
    for batch in batches:
        batch_id = batch['id']
        summary_status = check_summary_status(batch_id, context)
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
    for batch in batches:
        batch_id = batch['id']
        
        # Skip if marked for skipping
        if batch.get('skip', False):
            logger.info(f"Skipping batch {batch_id} (already complete)")
            continue
        
        # Create job script
        job_info = create_summary_job(
            batch_id=batch_id,
            job_manager=job_manager,
            config_path=args.config,
            blast_only=args.blast_only,
            threads=args.threads,
            memory=args.memory,
            time=args.time,
            log_dir=args.log_dir
        )
        
        # Add force flag if requested
        if args.force:
            command = job_info['script_path'].replace("--batch-id", "--force --batch-id")
            # Update the script file
            with open(job_info['script_path'], 'r') as f:
                content = f.read()
            
            content = content.replace("--batch-id", "--force --batch-id")
            
            with open(job_info['script_path'], 'w') as f:
                f.write(content)
        
        # Submit job
        job_id = submit_job(job_info, job_manager, context)
        
        if job_id:
            submitted_jobs.append({
                'batch_id': batch_id,
                'job_id': job_id
            })
    
    logger.info(f"Submitted {len(submitted_jobs)} jobs")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())