#!/usr/bin/env python3
"""
simple_check_summaries.py - Simple check of domain summary completion
"""

import os
import sys
import logging
import argparse
from typing import Dict, List, Any, Optional

# Add parent directory to path for imports
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
    """Main function to check domain summary quality"""
    parser = argparse.ArgumentParser(description='Check domain summary completion')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger("ecod.summary_check")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get all batches with summaries
    batch_query = """
    SELECT 
        b.id, b.batch_name, b.total_items
    FROM 
        ecod_schema.batch b
    WHERE 
        b.status = 'indexed' AND b.total_items = 5000
    ORDER BY 
        b.id
    """
    
    batches = context.db.execute_query(batch_query)
    
    if not batches:
        logger.warning("No batches found")
        return 1
    
    logger.info(f"Found {len(batches)} indexed batches")
    
    # For each batch, check summary completion
    for batch in batches:
        batch_id = batch[0]
        batch_name = batch[1]
        total_items = batch[2]
        
        # Get summary count for this batch
        summary_query = """
        SELECT 
            COUNT(DISTINCT pf.id) as summary_count
        FROM 
            ecod_schema.process_status ps
        JOIN 
            ecod_schema.batch b ON ps.batch_id = b.id
        LEFT JOIN 
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 
            b.id = %s
            AND pf.file_type = 'domain_summary'
            AND pf.file_exists = TRUE
        """
        
        summary_result = context.db.execute_query(summary_query, (batch_id,))
        summary_count = summary_result[0][0] if summary_result else 0
        
        # Get some basic stats about the summaries
        size_query = """
        SELECT 
            MIN(pf.file_size) as min_size,
            MAX(pf.file_size) as max_size,
            AVG(pf.file_size) as avg_size
        FROM 
            ecod_schema.process_file pf
        JOIN 
            ecod_schema.process_status ps ON pf.process_id = ps.id
        WHERE 
            ps.batch_id = %s
            AND pf.file_type = 'domain_summary'
            AND pf.file_exists = TRUE
        """
        
        size_result = context.db.execute_query(size_query, (batch_id,))
        
        # Calculate completion percentage
        completion_pct = (summary_count / total_items) * 100 if total_items > 0 else 0
        
        # Report results
        logger.info(f"Batch {batch_id} ({batch_name}): {summary_count}/{total_items} summaries ({completion_pct:.1f}%)")
        
        if size_result and size_result[0][0] is not None:
            min_size = size_result[0][0]
            max_size = size_result[0][1]
            avg_size = size_result[0][2]
            logger.info(f"  File sizes: min={min_size}, max={max_size}, avg={avg_size:.1f}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())