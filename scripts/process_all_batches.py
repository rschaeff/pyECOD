#!/usr/bin/env python3
"""
Process HHSearch results for multiple batches

This script processes HHSearch result files for all specified batches using the
improved hybrid column/regex parsing approach.
"""

import os
import sys
import logging
import argparse
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Tuple, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import from pipeline module for consistency
try:
    from ecod.core.context import ApplicationContext
    from ecod.pipelines.domain_analysis.hhresult_registrar import HHResultRegistrar
except ImportError:
    print("Could not import required modules. Make sure you're running this script from the correct directory.")
    sys.exit(1)


def setup_logging(verbose=False, log_file=None):
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


def process_batch(batch_id, config_path, force=False):
    """Process a single batch"""
    logger = logging.getLogger(f"process_batch_{batch_id}")
    logger.info(f"Processing HHSearch results for batch {batch_id}")
    
    # Initialize application context
    context = ApplicationContext(config_path)
    
    # Create registrar
    registrar = HHResultRegistrar(context)
    
    try:
        # Process batch
        result = registrar.register_batch_results(batch_id, force)
        
        logger.info(f"Successfully registered {result} HHR files for batch {batch_id}")
        return batch_id, result
    except Exception as e:
        logger.error(f"Error processing batch {batch_id}: {e}")
        logger.exception("Stack trace:")
        return batch_id, -1


def process_all_batches(batch_ids, config_path, force=False, max_workers=4):
    """Process multiple batches in parallel"""
    logger = logging.getLogger("process_all_batches")
    logger.info(f"Processing {len(batch_ids)} batches with {max_workers} workers")
    
    results = {}
    failed_batches = []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_batch = {
            executor.submit(process_batch, batch_id, config_path, force): batch_id
            for batch_id in batch_ids
        }
        
        for future in as_completed(future_to_batch):
            batch_id = future_to_batch[future]
            try:
                batch_id, result = future.result()
                results[batch_id] = result
                
                if result < 0:
                    failed_batches.append(batch_id)
                    logger.error(f"Batch {batch_id} failed processing")
                else:
                    logger.info(f"Completed batch {batch_id}: registered {result} files")
            except Exception as e:
                failed_batches.append(batch_id)
                logger.error(f"Batch {batch_id} failed with exception: {e}")
    
    return results, failed_batches


def get_all_batch_ids(context):
    """Get all batch IDs from the database"""
    query = "SELECT id, name FROM ecod_schema.batch ORDER BY id DESC"
    batches = context.db.execute_query(query)
    return [(row[0], row[1]) for row in batches]


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Process HHSearch results for multiple batches')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-ids', nargs='+', type=int, default=None,
                      help='Specific batch IDs to process')
    parser.add_argument('--exclude-batch-ids', nargs='+', type=int, default=[62],
                      help='Batch IDs to exclude')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')
    parser.add_argument('--max-workers', type=int, default=4,
                      help='Maximum number of worker threads')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)

    logger = logging.getLogger("process_all_batches")
    logger.info("Starting HHSearch registration for multiple batches")

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    # Initialize application context
    context = ApplicationContext(args.config)

    # Get batch IDs to process
    if args.batch_ids:
        # Use specified batch IDs
        batch_ids = [b_id for b_id in args.batch_ids if b_id not in args.exclude_batch_ids]
    else:
        # Get all batch IDs from database
        all_batches = get_all_batch_ids(context)
        batch_ids = [b_id for b_id, _ in all_batches if b_id not in args.exclude_batch_ids]
    
    logger.info(f"Will process {len(batch_ids)} batches: {batch_ids}")

    # Process all batches
    results, failed_batches = process_all_batches(
        batch_ids,
        args.config,
        args.force,
        args.max_workers
    )
    
    # Print summary
    total_processed = sum(count for count in results.values() if count > 0)
    logger.info(f"Completed processing all batches")
    logger.info(f"Total files processed: {total_processed}")
    
    if failed_batches:
        logger.error(f"Failed batches: {failed_batches}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
