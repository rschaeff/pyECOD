#!/usr/bin/env python3
"""
update_batch_reference.py - Update batch reference version
"""

import os
import sys
import logging
import argparse
from typing import Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

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

def main():
    """Update batch reference version"""
    parser = argparse.ArgumentParser(description='Update batch reference version')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to update')
    parser.add_argument('--reference', type=str, required=True,
                      help='New reference version (e.g., v291)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.update_reference")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get current batch information
    query = """
    SELECT batch_name, ref_version
    FROM ecod_schema.batch
    WHERE id = %s
    """
    
    result = context.db.execute_query(query, (args.batch_id,))
    
    if not result:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
    
    batch_name = result[0][0]
    current_ref = result[0][1]
    
    logger.info(f"Updating reference version for batch {batch_name} (ID: {args.batch_id})")
    logger.info(f"Current reference: {current_ref}")
    logger.info(f"New reference: {args.reference}")
    
    # Update batch reference version
    update_query = """
    UPDATE ecod_schema.batch
    SET ref_version = %s
    WHERE id = %s
    """
    
    context.db.execute_query(update_query, (args.reference, args.batch_id))
    logger.info(f"Successfully updated batch reference version to {args.reference}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())