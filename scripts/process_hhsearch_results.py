#!/usr/bin/env python3
"""
process_hhsearch_results.py - Process HHSearch results and integrate with BLAST evidence

This script processes HHSearch results for protein chains, converts them to XML,
and integrates them with BLAST evidence to create domain summaries.
"""

import os
import sys
import argparse
import logging
from typing import Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.exceptions import PipelineError
from ecod.error_handlers import handle_exceptions
from ecod.pipelines.hhsearch.processor import HHSearchProcessor  # Corrected import path

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

@handle_exceptions(exit_on_error=True)
def main():
    """Main entry point for HHSearch results processing script"""
    parser = argparse.ArgumentParser(description='Process ECOD HHSearch Results')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--protein-id', type=int,
                      help='Process a specific protein (optional)')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    logger = logging.getLogger("ecod.hhsearch_processor")
    
    # Check if batch exists
    query = "SELECT id, batch_name, ref_version FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_dict_query(query, (args.batch_id,))
    
    if not result:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
    
    batch_info = result[0]
    logger.info(f"Processing batch {args.batch_id} ({batch_info['batch_name']})")
    
    # Process HHSearch results
    try:
        logger.info(f"Starting HHSearch results processing for batch {args.batch_id}")
        
        processor = HHSearchProcessor(context)
        
        if args.protein_id:
            # Process a specific protein
            logger.info(f"Processing protein ID {args.protein_id}")
            success = processor.process_protein(args.protein_id, args.force)
            processed_count = 1 if success else 0
        else:
            # Process all proteins in the batch
            processed_count = processor.process_batch(args.batch_id, args.force)
        
        if processed_count > 0:
            logger.info(f"Successfully processed HHSearch results for {processed_count} chains in batch {args.batch_id}")
            return 0
        else:
            logger.warning(f"No chains were processed in batch {args.batch_id}")
            
            # Check for chains that need processing
            chains_query = """
            SELECT COUNT(*) 
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE ps.batch_id = %s
            AND pf.file_type = 'hhr'
            AND pf.file_exists = TRUE
            AND NOT EXISTS (
                SELECT 1 FROM ecod_schema.process_file 
                WHERE process_id = ps.id AND file_type = 'domain_summary'
            )
            """
            pending = context.db.execute_query(chains_query, (args.batch_id,))
            
            if pending and pending[0][0] > 0:
                logger.info(f"Found {pending[0][0]} chains with HHR files that need processing")
                logger.info("Consider using the --force option to reprocess these chains")
            
            return 0
    
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