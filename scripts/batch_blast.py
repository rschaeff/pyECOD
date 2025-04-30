#!/usr/bin/env python3
import os
import sys
import argparse
import logging
from typing import List, Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
# Import required modules from pyECOD
from ecod.core.context import ApplicationContext
from ecod.core.logging_config import LoggingManager
from ecod.pipelines.blast_pipeline import BlastPipeline
from ecod.exceptions import ConfigurationError, PipelineError

def setup_logging(verbose: bool, log_file: Optional[str]) -> logging.Logger:
    """Set up logging configuration"""
    return LoggingManager.configure(
        verbose=verbose,
        log_file=log_file,
        component="batch_blast",
    )

def get_pending_batches(context: ApplicationContext) -> List[Dict[str, Any]]:
    """Get batches that need blast processing"""
    query = """
    SELECT 
        b.id, b.batch_name, b.base_path, b.type, b.ref_version
    FROM 
        ecod_schema.batch b
    WHERE 
        b.status = 'created'
        OR EXISTS (
            SELECT 1 FROM ecod_schema.process_status ps
            WHERE ps.batch_id = b.id
            AND ps.status IN ('pending', 'processing')
        )
    ORDER BY 
        b.id
    """
    
    return context.db.execute_dict_query(query)

def get_batch_statistics(context: ApplicationContext, batch_id: int) -> Dict[str, int]:
    """Get statistics for a batch"""
    query = """
    SELECT
        COUNT(*) as total,
        SUM(CASE WHEN ps.status = 'pending' THEN 1 ELSE 0 END) as pending,
        SUM(CASE WHEN ps.status = 'processing' THEN 1 ELSE 0 END) as processing,
        SUM(CASE WHEN ps.status = 'success' THEN 1 ELSE 0 END) as completed,
        SUM(CASE WHEN ps.status = 'error' THEN 1 ELSE 0 END) as failed
    FROM
        ecod_schema.process_status ps
    WHERE
        ps.batch_id = %s
    """
    
    result = context.db.execute_dict_query(query, (batch_id,))
    return result[0] if result else {}

def run_chain_blast(pipeline: BlastPipeline, batch_id: int, threads: int, 
                  memory: str, batch_size: int, force: bool) -> List[str]:
    """Run chain blast for a batch"""
    try:
        return pipeline.run_chain_blast(
            batch_id=batch_id,
            batch_size=batch_size
        )
    except (ConfigurationError, PipelineError) as e:
        logger.error(f"Error running chain blast: {str(e)}")
        return []

def run_domain_blast(pipeline: BlastPipeline, batch_id: int, threads: int, 
                   memory: str, batch_size: int, force: bool) -> List[str]:
    """Run domain blast for a batch"""
    try:
        return pipeline.run_domain_blast(
            batch_id=batch_id,
            batch_size=batch_size
        )
    except (ConfigurationError, PipelineError) as e:
        logger.error(f"Error running domain blast: {str(e)}")
        return []

def process_batch(context: ApplicationContext, batch_id: int, threads: int, 
                memory: str, batch_size: int, force: bool) -> None:
    """Process a single batch"""
    logger.info(f"Processing batch {batch_id}")
    
    # Get batch statistics before processing
    stats_before = get_batch_statistics(context, batch_id)
    logger.info(f"Batch {batch_id} statistics before processing: {stats_before}")
    
    # Create pipeline
    pipeline = BlastPipeline(context)
    
    # Run chain blast
    logger.info(f"Running chain blast for batch {batch_id}")
    chain_jobs = run_chain_blast(pipeline, batch_id, threads, memory, batch_size, force)
    logger.info(f"Submitted {len(chain_jobs)} chain blast jobs")
    
    # Run domain blast
    logger.info(f"Running domain blast for batch {batch_id}")
    domain_jobs = run_domain_blast(pipeline, batch_id, threads, memory, batch_size, force)
    logger.info(f"Submitted {len(domain_jobs)} domain blast jobs")
    
    # Update batch status if no jobs were submitted
    if not chain_jobs and not domain_jobs:
        # Get current batch status
        query = "SELECT status FROM ecod_schema.batch WHERE id = %s"
        result = context.db.execute_query(query, (batch_id,))
        current_status = result[0][0] if result else None
        
        # Update status to 'processed' if currently 'created'
        if current_status == 'created':
            context.db.update(
                "ecod_schema.batch",
                {"status": "processed"},
                "id = %s",
                (batch_id,)
            )
            logger.info(f"Updated batch {batch_id} status to 'processed'")
    
    # Wait for jobs to complete if requested
    if chain_jobs or domain_jobs:
        logger.info(f"Checking job status for batch {batch_id}")
        pipeline.check_job_status(batch_id)

def main():
    parser = argparse.ArgumentParser(description="Run BLAST for ECOD batches")
    parser.add_argument("--config", default="config/config.yml", help="Path to config file")
    parser.add_argument("--batch-id", type=int, help="Process specific batch ID")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--memory", default="4G", help="Memory allocation")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for job creation")
    parser.add_argument("--force", action="store_true", help="Force reprocessing")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--log-file", help="Log file path")
    args = parser.parse_args()
    
    # Setup logging
    global logger
    logger = setup_logging(args.verbose, args.log_file)
    
    # Create application context
    context = ApplicationContext(args.config)
    
    # Process specific batch or all pending batches
    if args.batch_id:
        process_batch(
            context=context,
            batch_id=args.batch_id,
            threads=args.threads,
            memory=args.memory,
            batch_size=args.batch_size,
            force=args.force
        )
    else:
        # Get all pending batches
        batches = get_pending_batches(context)
        logger.info(f"Found {len(batches)} pending batches")
        
        for batch in batches:
            process_batch(
                context=context,
                batch_id=batch['id'],
                threads=args.threads,
                memory=args.memory,
                batch_size=args.batch_size,
                force=args.force
            )

if __name__ == "__main__":
    main()
