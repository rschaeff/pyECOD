#!/usr/bin/env python3
"""
run_domain_analysis_v2.py - Run domain summary and partition on BLAST results

This script runs the domain analysis pipeline on a specified batch,
with options for targeted testing of the partition model.
"""

import os
import sys
import argparse
import logging
from typing import Dict, Any, Optional, List

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.exceptions import ECODError, PipelineError
from ecod.error_handlers import handle_exceptions
from ecod.core.validation import validate_file_path

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


def verify_summary_completion(context, batch_id):
    """Verify that domain summaries are complete for all proteins in batch"""
    logger = logging.getLogger("ecod.domain_analysis.debug")
    
    # Query to check summary completion
    query = """
    SELECT 
        COUNT(*) as total,
        SUM(CASE WHEN EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf
            WHERE pf.process_id = ps.id
            AND pf.file_type = 'domain_summary'
            AND pf.file_exists = TRUE
        ) THEN 1 ELSE 0 END) as complete_count
    FROM 
        ecod_schema.process_status ps
    WHERE 
        ps.batch_id = %s
    """
    
    results = context.db.execute_dict_query(query, (batch_id,))[0]
    total = results.get('total', 0)
    complete = results.get('complete_count', 0)
    is_complete = total > 0 and total == complete
    
    logger.info(f"Summary completion: {complete}/{total} ({is_complete})")
    
    return {
        'total_count': total,
        'complete_count': complete,
        'complete': is_complete
    }

def validate_partition_inputs(context, batch_id, process_ids=None):
    """Validate inputs required for partition"""
    logger = logging.getLogger("ecod.domain_analysis.debug")
    
    # Base query to get relevant process information
    query = """
    SELECT 
        ps.id as process_id, 
        p.pdb_id, 
        p.chain_id,
        pf_summary.file_path as summary_path,
        pf_blast.file_path as blast_path,
        b.base_path
    FROM 
        ecod_schema.process_status ps
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    JOIN
        ecod_schema.batch b ON ps.batch_id = b.id
    LEFT JOIN
        ecod_schema.process_file pf_summary ON ps.id = pf_summary.process_id AND pf_summary.file_type = 'domain_summary'
    LEFT JOIN
        ecod_schema.process_file pf_blast ON ps.id = pf_blast.process_id AND pf_blast.file_type = 'domain_blast_result'
    WHERE 
        ps.batch_id = %s
    """
    
    # Add process ID filter if specified
    if process_ids:
        query += f" AND ps.id IN ({','.join(map(str, process_ids))})"
    
    # Execute query
    rows = context.db.execute_dict_query(query, (batch_id,))
    
    # Check results
    valid_processes = []
    issues = []
    
    for row in rows:
        process_id = row['process_id']
        pdb_id = row['pdb_id']
        chain_id = row['chain_id']
        base_path = row['base_path']
        
        # Build full paths
        summary_path = os.path.join(base_path, row['summary_path']) if row['summary_path'] else None
        blast_path = os.path.join(base_path, row['blast_path']) if row['blast_path'] else None
        
        # Check files exist
        if not summary_path:
            issues.append(f"Process {process_id} ({pdb_id}_{chain_id}): Missing domain summary")
            continue
            
        if not os.path.exists(summary_path):
            issues.append(f"Process {process_id} ({pdb_id}_{chain_id}): Domain summary file not found at {summary_path}")
            continue
            
        if not blast_path:
            issues.append(f"Process {process_id} ({pdb_id}_{chain_id}): Missing BLAST result")
            continue
            
        if not os.path.exists(blast_path):
            issues.append(f"Process {process_id} ({pdb_id}_{chain_id}): BLAST result file not found at {blast_path}")
            continue
            
        # Process is valid
        valid_processes.append(process_id)
    
    # Return validation results
    result = {
        'valid': len(issues) == 0,
        'total_count': len(rows),
        'valid_count': len(valid_processes),
        'valid_processes': valid_processes,
        'issues': issues
    }
    
    logger.info(f"Validation results: {result['valid_count']}/{result['total_count']} valid")
    if issues:
        logger.warning(f"Validation issues: {len(issues)}")
        for issue in issues[:5]:  # Log first 5 issues
            logger.warning(f"  - {issue}")
        if len(issues) > 5:
            logger.warning(f"  ... and {len(issues) - 5} more issues")
    
    return result

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
                      
    # New options for targeted testing
    parser.add_argument('--partition-only', action='store_true',
                      help='Run only the partition step on batches with complete summaries')
    parser.add_argument('--process-ids', type=str,
                      help='Comma-separated list of process IDs to analyze')
    parser.add_argument('--validate-only', action='store_true',
                      help='Only validate database and filesystem structures without running partition')
 
    parser.add_argument('--force', action='store_true',
                       help="Force overwrite of existing domain files")

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

    if args.force:
        context.update_config('pipeline', 'force_overwrite', True)
        context.logger.info("Force overwrite enabled - will regenerate all files")
    
    # Parse process IDs if provided
    process_ids = None
    if args.process_ids:
        process_ids = [int(pid.strip()) for pid in args.process_ids.split(',')]
        logger.info(f"Processing specific process IDs: {process_ids}")
    
    # Handle validation-only mode
    if args.validate_only:
        logger.info(f"Validating partition inputs for batch {args.batch_id}")
        validation_result = validate_partition_inputs(context, args.batch_id, process_ids)
        
        if validation_result['valid']:
            print("Validation successful. All required structures are in place.")
            print(f"Found {validation_result['valid_count']}/{validation_result['total_count']} valid processes.")
            return 0
        else:
            print(f"Validation failed. {validation_result['valid_count']}/{validation_result['total_count']} valid processes.")
            print("Issues:")
            for issue in validation_result['issues'][:10]:  # Show first 10 issues
                print(f"  - {issue}")
            if len(validation_result['issues']) > 10:
                print(f"  ... and {len(validation_result['issues']) - 10} more issues")
            return 1
    
    # Handle partition-only mode
    if args.partition_only:
        logger.info(f"Running partition only for batch {args.batch_id}")
        
        # Check if summaries are complete
        summary_status = verify_summary_completion(context, args.batch_id)
        if not summary_status['complete']:
            logger.error(f"Cannot run partition: Summaries incomplete ({summary_status['complete_count']}/{summary_status['total_count']} complete)")
            print(f"Error: Summaries are not complete for batch {args.batch_id}")
            return 1
        
        # Initialize partition component directly
        try:
            from ecod.pipelines.domain_analysis.partition import DomainPartition
            
            # Create partition object
            partition = DomainPartition(context.config_manager.config_path)

            if args.force:
                partition.config['force_overwrite'] = True

            # Get base path and reference
            query = "SELECT base_path, ref_version FROM ecod_schema.batch WHERE id = %s"
            batch_info = context.db.execute_dict_query(query, (args.batch_id,))[0]
            base_path = batch_info['base_path']
            reference = batch_info['ref_version']
            
            # Run partition
            if process_ids:
                logger.info(f"Running partition for {len(process_ids)} specific processes")
                # A method would need to be implemented in DomainPartition to support specific process IDs
                success = partition.process_specific_ids(args.batch_id, process_ids, base_path, reference, args.blast_only, args.force)
            else:
                logger.info(f"Running partition for batch {args.batch_id}")
                success = partition.process_batch(args.batch_id, base_path, reference, args.blast_only, args.limit, args.force)
            
            if success:
                logger.info(f"Partition completed successfully for batch {args.batch_id}")
                print(f"Partition completed successfully for batch {args.batch_id}")
                return 0
            else:
                logger.error(f"Partition failed for batch {args.batch_id}")
                print(f"Partition failed for batch {args.batch_id}")
                return 1
                
        except ImportError as e:
            logger.error(f"Error importing partition module: {str(e)}")
            print(f"Error: Could not import the domain partition module: {str(e)}")
            return 1
        except Exception as e:
            logger.error(f"Error in partition: {str(e)}", exc_info=True)
            print(f"Error in partition: {str(e)}")
            return 1
    
    # Standard domain analysis pipeline
    try:
        # Import the pipeline module
        from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
        
        # Initialize pipeline
        pipeline = DomainAnalysisPipeline(context.config_manager.config_path)
        
        logger.info(f"Starting domain analysis for batch {args.batch_id} (blast_only={args.blast_only})")
        print(f"Starting domain analysis for batch {args.batch_id}...")
        
        # Run pipeline with process IDs if specified
        if process_ids:
            success = pipeline.process_proteins(args.batch_id, process_ids, args.blast_only)
        else:
            success = pipeline.run_pipeline(args.batch_id, args.blast_only, args.limit)
        
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