#!/usr/bin/env python3
"""
run_domain_analysis.py - Run domain summary and partition on BLAST results

This script runs the domain analysis pipeline on a specified batch,
with an option to use only BLAST results (--blast-only flag).
"""

import os
import sys
import argparse
import logging
from typing import Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.exceptions import ECODError, PipelineError
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
    
def verify_database_records(context, batch_id):
    """Verify that required database records exist for the batch"""
    logger = logging.getLogger("ecod.domain_analysis.debug")
    
    # Check process_status records
    status_query = """
    SELECT 
        COUNT(*) as total,
        SUM(CASE WHEN status = 'success' THEN 1 ELSE 0 END) as success,
        SUM(CASE WHEN status = 'error' THEN 1 ELSE 0 END) as error,
        SUM(CASE WHEN status = 'processing' THEN 1 ELSE 0 END) as processing,
        SUM(CASE WHEN status = 'pending' THEN 1 ELSE 0 END) as pending
    FROM 
        ecod_schema.process_status
    WHERE 
        batch_id = %s
    """
    
    status_results = context.db.execute_dict_query(status_query, (batch_id,))[0]
    logger.info(f"Process status counts: {status_results}")
    
    # Check process_file records
    file_query = """
    SELECT 
        file_type, COUNT(*) as count,
        SUM(CASE WHEN file_exists = TRUE THEN 1 ELSE 0 END) as exist_count
    FROM 
        ecod_schema.process_file pf
    JOIN
        ecod_schema.process_status ps ON pf.process_id = ps.id
    WHERE 
        ps.batch_id = %s
    GROUP BY 
        file_type
    """
    
    file_results = context.db.execute_dict_query(file_query, (batch_id,))
    logger.info("Process file counts:")
    for row in file_results:
        logger.info(f"  - {row['file_type']}: {row['count']} records, {row['exist_count']} exist")
    
    # Check for domain_summary files specifically
    summary_query = """
    SELECT COUNT(*) 
    FROM ecod_schema.process_file pf
    JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
    WHERE ps.batch_id = %s AND pf.file_type = 'domain_summary'
    """
    
    summary_count = context.db.execute_query(summary_query, (batch_id,))[0][0]
    logger.info(f"Domain summary files: {summary_count}")
    
    # Check for proteins with both BLAST results but no summary
    ready_query = """
    SELECT COUNT(*)
    FROM ecod_schema.process_status ps
    WHERE ps.batch_id = %s
    AND EXISTS (
        SELECT 1 FROM ecod_schema.process_file pf1
        WHERE pf1.process_id = ps.id AND pf1.file_type = 'chain_blast_result' AND pf1.file_exists = TRUE
    )
    AND EXISTS (
        SELECT 1 FROM ecod_schema.process_file pf2
        WHERE pf2.process_id = ps.id AND pf2.file_type = 'domain_blast_result' AND pf2.file_exists = TRUE
    )
    AND NOT EXISTS (
        SELECT 1 FROM ecod_schema.process_file pf3
        WHERE pf3.process_id = ps.id AND pf3.file_type = 'domain_summary'
    )
    """
    
    ready_count = context.db.execute_query(ready_query, (batch_id,))[0][0]
    logger.info(f"Proteins ready for domain summary: {ready_count}")
    
    return {
        'status_counts': status_results,
        'file_counts': {row['file_type']: row['count'] for row in file_results},
        'summary_count': summary_count,
        'ready_count': ready_count
    }

@handle_exceptions(exit_on_error=True)
def main():
    """Main entry point for domain analysis script"""
    parser = argparse.ArgumentParser(description='Run ECOD Domain Analysis')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--blast-only', action='store_true',
                      help='Use only BLAST results (no HHSearch)')
    parser.add_argument('--limit', type=int, default=None,
                      help='Maximum number of proteins to process')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    logger = logging.getLogger("ecod.domain_analysis")
    
    # Check if batch exists
    query = "SELECT id FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_query(query, (args.batch_id,))
    
    if not result:
        logger.error(f"Batch {args.batch_id} not found")
        print(f"Error: Batch {args.batch_id} not found")
        return 1
    
    # Initialize domain analysis pipeline
    try:
        # Import the pipeline module here to avoid circular imports
        from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
        
        # Use the same config path for consistency
        pipeline = DomainAnalysisPipeline(context.config_manager.config_path)
        
        logger.info(f"Starting domain analysis for batch {args.batch_id} (blast_only={args.blast_only})")
        print(f"Starting domain analysis for batch {args.batch_id}...")
        
        # Run the domain analysis pipeline
        success = pipeline.run_pipeline(
            args.batch_id,
            args.blast_only,
            args.limit
        )
        
        if success:
            logger.info(f"Domain analysis completed successfully for batch {args.batch_id}")
            print(f"Domain analysis completed successfully for batch {args.batch_id}")
            return 0
        else:
            logger.error(f"Domain analysis failed for batch {args.batch_id}")
            print(f"Domain analysis failed for batch {args.batch_id}")
            return 1
    except ImportError as e:
        logger.error(f"Error importing pipeline module: {str(e)}")
        print(f"Error: Could not import the domain analysis pipeline module: {str(e)}")
        print("Make sure the ecod package is properly installed and accessible.")
        return 1
    except PipelineError as e:
        logger.error(f"Pipeline error: {str(e)}", exc_info=True)
        print(f"Pipeline error: {str(e)}")
        print("Check the log file for more details.")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}", exc_info=True)
        print(f"Unexpected error: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())