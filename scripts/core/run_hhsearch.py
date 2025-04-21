#!/usr/bin/env python3
"""
run_hhsearch.py - Run HHSearch pipeline on proteins with BLAST results

This script runs the HHSearch pipeline for protein chains that already have
BLAST results but need HHblits profiles and HHSearch results.
"""

import os
import sys
import argparse
import logging
import time
from typing import Dict, Any, Optional, List

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

from ecod.core.context import ApplicationContext
from ecod.exceptions import ECODError, PipelineError, JobSubmissionError
from ecod.error_handlers import handle_exceptions

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Configure logging
    
    Args:
        verbose: Enable debug logging if True
        log_file: Optional log file path
    """
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

def get_chains_with_blast_results(context: ApplicationContext, batch_id: int) -> List[Dict[str, Any]]:
    """Get chains that have BLAST results but not HHSearch profiles"""
    logger = logging.getLogger("ecod.hhsearch")
    
    # Query to find chains with BLAST results but no HHSearch profiles
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, p.source_id, p.length, 
        ps.id as process_id, ps.relative_path
    FROM 
        ecod_schema.protein p
    JOIN
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN
        ecod_schema.process_file pf ON ps.id = pf.process_id 
    WHERE 
        ps.batch_id = %s
        AND ps.status = 'success'
        AND ps.is_representative = TRUE  -- Filter for representatives
        AND pf.file_type IN ('chain_blast_result', 'domain_blast_result')
        AND pf.file_exists = TRUE
        AND NOT EXISTS (
            SELECT 1 FROM ecod_schema.process_file 
            WHERE process_id = ps.id AND file_type = 'a3m'
        )
    GROUP BY 
        p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.id, ps.relative_path
    HAVING 
        COUNT(DISTINCT pf.file_type) >= 1
    """
    
    try:
        rows = context.db.execute_dict_query(query, (batch_id,))
        logger.info(f"Found {len(rows)} representative chains with BLAST results but no HHSearch profiles")
        
        # Add sequence to each chain - Modified this part
        result = []
        for row in rows:
            seq_query = """
            SELECT sequence FROM ecod_schema.protein_sequence WHERE protein_id = %s
            """
            seq_result = context.db.execute_query(seq_query, (row['id'],))
            if seq_result and len(seq_result) > 0:
                row_copy = row.copy()  # Create copy to avoid modifying the original
                row_copy['sequence'] = seq_result[0][0]  # Access as tuple
                result.append(row_copy)
        
        logger.info(f"Found {len(result)} representative chains with sequences ready for HHSearch")
        return result
    except Exception as e:
        logger.error(f"Error querying representative chains with BLAST results: {str(e)}")
        return []

@handle_exceptions(exit_on_error=True)
def run_hhsearch_pipeline(context: ApplicationContext, batch_id: int, threads: int = 8, 
                        memory: str = "16G", check_interval: int = 300, 
                        wait: bool = True, force: bool = False) -> bool:
    """Run HHSearch pipeline for a batch
    
    Args:
        context: Application context
        batch_id: Batch ID to process
        threads: Number of threads for HHblits/HHsearch
        memory: Memory allocation for jobs
        check_interval: Interval in seconds to check job status
        wait: Whether to wait for jobs to complete
        force: Force regeneration of profiles even if they exist
        
    Returns:
        True if successful
    """
    logger = logging.getLogger("ecod.hhsearch")
    
    try:
        # Import the necessary modules
        from ecod.jobs import SlurmJobManager
        from ecod.pipelines.hhsearch_pipeline import HHSearchPipeline
        
        # Initialize job manager
        slurm_manager = create_job_manager(context.config_manager.config, manager_type="slurm")
        
        # Initialize HHSearch pipeline
        hhsearch_pipeline = HHSearchPipeline(context)
        
        # Get chains ready for processing
        chains = get_chains_with_blast_results(context, batch_id)
        
        if not chains:
            logger.warning(f"No chains found with BLAST results but without HHSearch profiles in batch {batch_id}")
            return False
        
        logger.info(f"Found {len(chains)} chains ready for HHSearch processing")
        
        # Check if batch exists and update if needed
        existing_batch_query = """
        SELECT id, type FROM ecod_schema.batch 
        WHERE id = %s
        """
        existing_batch = context.db.execute_dict_query(existing_batch_query, (batch_id,))
        
        if existing_batch:
            current_type = existing_batch[0]['type']
            if 'hhsearch' not in current_type:
                # Update batch type
                new_type = f"{current_type}_hhsearch"
                context.db.update(
                    "ecod_schema.batch",
                    {"type": new_type},
                    "id = %s",
                    (batch_id,)
                )
                logger.info(f"Updated batch type from '{current_type}' to '{new_type}'")
        else:
            logger.error(f"Batch {batch_id} does not exist")
            return False
        
        # Run profile generation
        logger.info("Generating HHblits profiles")
        profile_job_ids = hhsearch_pipeline.generate_profiles(batch_id, threads, memory, force)
        
        if not profile_job_ids or len(profile_job_ids) == 0:
            logger.error("Failed to submit profile generation jobs")
            return False
        
        logger.info(f"Submitted {len(profile_job_ids)} profile generation jobs")
        
        if wait:
            # Wait for profile generation to complete
            logger.info(f"Waiting for profile generation jobs to complete (checking every {check_interval} seconds)")
            
            while True:
                # Check job status
                completed, failed, running = slurm_manager.check_all_jobs(batch_id)
                
                logger.info(f"Profile generation: {completed} completed, {failed} failed, {running} running")
                
                if running == 0:
                    if failed > 0:
                        logger.warning(f"{failed} profile generation jobs failed")
                    break
                
                # Wait before checking again
                time.sleep(check_interval)
            
            # Run HHSearch
            logger.info("Running HHSearch")
            search_job_ids = hhsearch_pipeline.run_hhsearch(batch_id, threads, memory)
            
            if not search_job_ids or len(search_job_ids) == 0:
                logger.error("Failed to submit HHSearch jobs")
                return False
            
            logger.info(f"Submitted {len(search_job_ids)} HHSearch jobs")
            
            # Wait for HHSearch to complete
            logger.info(f"Waiting for HHSearch jobs to complete (checking every {check_interval} seconds)")
            
            while True:
                # Check job status
                completed, failed, running = slurm_manager.check_all_jobs(batch_id)
                
                logger.info(f"HHSearch: {completed} completed, {failed} failed, {running} running")
                
                if running == 0:
                    if failed > 0:
                        logger.warning(f"{failed} HHSearch jobs failed")
                    break
                
                # Wait before checking again
                time.sleep(check_interval)
        
        return True
    
    except ImportError as e:
        logger.error(f"Error importing required modules: {str(e)}")
        raise PipelineError(f"Error importing required modules: {str(e)}")
    except JobSubmissionError as e:
        logger.error(f"Job submission error: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in HHSearch pipeline: {str(e)}", exc_info=True)
        raise PipelineError(f"HHSearch pipeline error: {str(e)}")

@handle_exceptions(exit_on_error=True)
def main():
    """Main entry point for HHSearch pipeline script"""
    parser = argparse.ArgumentParser(description='Run ECOD HHSearch Pipeline')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads for HHblits/HHsearch')
    parser.add_argument('--memory', type=str, default="16G",
                      help='Memory allocation for jobs')
    parser.add_argument('--check-interval', type=int, default=300,
                      help='Interval in seconds to check job status')
    parser.add_argument('--no-wait', action='store_true',
                      help='Do not wait for jobs to complete')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--force', action='store_true',
                  help='Force regeneration of profiles even if they exist')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    logger = logging.getLogger("ecod.hhsearch")
    
    # Check if batch exists
    query = "SELECT id FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_query(query, (args.batch_id,))
    
    if not result:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
    
    # Run HHSearch pipeline
    try:
        logger.info(f"Starting HHSearch pipeline for batch {args.batch_id}")
        
        success = run_hhsearch_pipeline(
            context,
            args.batch_id,
            args.threads,
            args.memory,
            args.check_interval,
            not args.no_wait,
            args.force
        )
        
        if success:
            logger.info(f"HHSearch pipeline started successfully for batch {args.batch_id}")
            return 0
        else:
            logger.error(f"HHSearch pipeline failed for batch {args.batch_id}")
            return 1
    
    except PipelineError as e:
        logger.error(f"Pipeline error: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}", exc_info=True)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())