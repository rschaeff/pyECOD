#!/usr/bin/env python3
"""
generate_batch_domain_summaries.py - Generate domain summaries for all proteins in a batch

This script processes proteins in a batch to generate domain summaries,
with proper handling of ApplicationContext and force_overwrite flag.
"""

import os
import sys
import time
import logging
import argparse
from typing import Dict, Any, Optional, List
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.summary import DomainSummary
from ecod.exceptions import PipelineError

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

def process_protein(protein_info: Dict[str, Any], context: ApplicationContext, output_dir: str, 
                   blast_only: bool, force: bool, verbose: bool) -> Dict[str, Any]:
    """Process a single protein for domain summary generation"""
    pdb_id = protein_info['pdb_id']
    chain_id = protein_info['chain_id']
    protein_id = protein_info['id']
    batch_id = protein_info['batch_id']
    batch_path = protein_info['batch_path']
    reference = protein_info['reference']
    
    # Create logger for this protein
    logger = logging.getLogger(f"ecod.generate_summary.{pdb_id}_{chain_id}")
    
    result = {
        'protein_id': protein_id,
        'pdb_id': pdb_id,
        'chain_id': chain_id,
        'success': False,
        'error': None,
        'output_file': None,
        'duration': 0
    }
    
    start_time = time.time()
    
    try:
        logger.info(f"Generating domain summary for protein: {pdb_id}_{chain_id} (ID: {protein_id})")
        
        # If force flag is set, update the context
        if force:
            context.set_force_overwrite(True)
            logger.debug("Force overwrite enabled")
        
        # Initialize domain summary processor with proper context object
        summary_processor = DomainSummary(context)
        
        # Use provided output dir or batch path
        effective_output_dir = output_dir or batch_path
        
        # Generate domain summary
        summary_file = summary_processor.create_summary(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            job_dump_dir=effective_output_dir,
            blast_only=blast_only
        )
        
        if summary_file:
            logger.info(f"Successfully generated domain summary: {summary_file}")
            result['success'] = True
            result['output_file'] = summary_file
            
            # Register the summary file in the database
            update_database_for_summary(context, protein_id, batch_id, summary_file, batch_path, logger)
        else:
            error_msg = "Failed to generate domain summary - no file returned"
            logger.error(error_msg)
            result['error'] = error_msg
            
    except Exception as e:
        error_msg = f"Error generating domain summary: {str(e)}"
        logger.error(error_msg, exc_info=verbose)
        result['error'] = error_msg
        
        # Update error in database
        try:
            update_process_error(context, protein_id, batch_id, error_msg, logger)
        except Exception as db_error:
            logger.error(f"Failed to update error in database: {str(db_error)}")
    
    end_time = time.time()
    result['duration'] = end_time - start_time
    
    return result

def update_database_for_summary(context, protein_id, batch_id, summary_file, batch_path, logger):
    """Update database with summary file information"""
    try:
        # First get the process ID
        process_query = """
        SELECT id FROM ecod_schema.process_status
        WHERE protein_id = %s AND batch_id = %s
        """
        process_result = context.db.execute_query(process_query, (protein_id, batch_id))
        
        if not process_result:
            logger.warning(f"No process found for protein {protein_id} in batch {batch_id}")
            return
            
        process_id = process_result[0][0]
        
        # Check if summary file record already exists
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'domain_summary'
        """
        
        existing = context.db.execute_query(check_query, (process_id,))
        
        # Get relative path and file size
        rel_path = os.path.relpath(summary_file, batch_path)
        file_size = os.path.getsize(summary_file) if os.path.exists(summary_file) else 0
        
        if existing:
            # Update existing record
            update_query = """
            UPDATE ecod_schema.process_file
            SET file_path = %s, file_exists = TRUE, file_size = %s, last_checked = NOW()
            WHERE id = %s
            """
            
            context.db.execute_query(update_query, (rel_path, file_size, existing[0][0]))
            logger.debug(f"Updated domain_summary file record (ID: {existing[0][0]})")
        else:
            # Create new record
            insert_query = """
            INSERT INTO ecod_schema.process_file
            (process_id, file_type, file_path, file_exists, file_size, last_checked)
            VALUES (%s, %s, %s, %s, %s, NOW())
            """
            
            context.db.execute_query(insert_query, (process_id, 'domain_summary', rel_path, True, file_size))
            logger.debug("Added new domain_summary file record")
        
        # Update process status
        status_query = """
        UPDATE ecod_schema.process_status
        SET current_stage = 'domain_summary', status = 'success', error_message = NULL, updated_at = NOW()
        WHERE id = %s
        """
        
        context.db.execute_query(status_query, (process_id,))
        logger.debug("Updated process status to domain_summary:success")
        
    except Exception as e:
        logger.error(f"Database update error: {str(e)}")
        raise

def update_process_error(context, protein_id, batch_id, error_message, logger):
    """Update process status with error message"""
    try:
        # Get process ID
        process_query = """
        SELECT id FROM ecod_schema.process_status
        WHERE protein_id = %s AND batch_id = %s
        """
        process_result = context.db.execute_query(process_query, (protein_id, batch_id))
        
        if not process_result:
            logger.warning(f"No process found for protein {protein_id} in batch {batch_id}")
            return
            
        process_id = process_result[0][0]
        
        # Update process status
        status_query = """
        UPDATE ecod_schema.process_status
        SET current_stage = 'domain_summary_failed', status = 'error', error_message = %s, updated_at = NOW()
        WHERE id = %s
        """
        
        context.db.execute_query(status_query, (error_message[:500], process_id))  # Truncate to avoid DB field size issues
        logger.debug("Updated process status to error")
        
    except Exception as e:
        logger.error(f"Error updating process status: {str(e)}")
        raise

def get_proteins_to_process(context, batch_id, limit=None, retry_failed=False, force=False, representative_only=False):
    """Get list of proteins to process from the database"""
    logger = logging.getLogger("ecod.batch_summary")
    
    # First get batch information
    batch_query = """
    SELECT 
        id, batch_name, base_path, ref_version
    FROM 
        ecod_schema.batch
    WHERE 
        id = %s
    """
    
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return None, None
    
    batch_info = {
        'id': batch_result[0][0],
        'name': batch_result[0][1],
        'path': batch_result[0][2],
        'reference': batch_result[0][3]
    }
    
    # Construct protein query based on arguments
    protein_query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id,
        ps.id as process_id, b.id as batch_id, b.base_path, b.ref_version
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    LEFT JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'domain_summary'
    WHERE 
        b.id = %s
    """
    
    # Add filters
    conditions = []
    
    if not retry_failed:
        conditions.append("ps.status != 'error'")
    
    if not force:
        # Only process proteins without a domain summary or where file doesn't exist
        conditions.append("(pf.id IS NULL OR pf.file_exists = FALSE)")
    
    if representative_only:
        conditions.append("ps.is_representative = TRUE")

    if conditions:
        protein_query += " AND " + " AND ".join(conditions)
    
    # Always sort by protein ID
    protein_query += " ORDER BY p.id"
    
    # Add limit if specified
    if limit:
        protein_query += f" LIMIT {limit}"


    
    # Execute query
    protein_results = context.db.execute_query(protein_query, (batch_id,))
    
    if not protein_results:
        logger.info(f"No proteins found to process in batch {batch_id}")
        return batch_info, []
    
    # Prepare protein info list
    proteins = []
    for row in protein_results:
        proteins.append({
            'id': row[0],
            'pdb_id': row[1],
            'chain_id': row[2],
            'process_id': row[3],
            'batch_id': row[4],
            'batch_path': row[5],
            'reference': row[6]
        })
    
    return batch_info, proteins

def update_batch_completion(context, batch_id, successful_count, logger):
    """Update batch completion status in database"""
    if successful_count <= 0:
        return
        
    try:
        update_query = """
        UPDATE ecod_schema.batch
        SET completed_items = completed_items + %s,
            status = CASE 
                WHEN completed_items + %s >= total_items THEN 'completed' 
                ELSE status 
            END,
            completed_at = CASE 
                WHEN completed_items + %s >= total_items THEN NOW() 
                ELSE completed_at 
            END
        WHERE id = %s
        """
        
        context.db.execute_query(update_query, (successful_count, successful_count, successful_count, batch_id))
        logger.info(f"Updated batch {batch_id} completion status (+{successful_count} completed items)")
    except Exception as e:
        logger.error(f"Error updating batch completion: {str(e)}")

def main():
    """Main function to generate domain summaries for a batch"""
    parser = argparse.ArgumentParser(description='Generate domain summaries for all proteins in a batch')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--blast-only', action='store_true',
                      help='Generate blast-only summaries (no HHSearch)')
    parser.add_argument('--output-dir', type=str,
                      help='Override output directory (default: from batch path)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('--limit', type=int,
                      help='Limit number of proteins to process')
    parser.add_argument('--threads', type=int, default=1,
                      help='Number of threads to use for parallel processing')
    parser.add_argument('--retry-failed', action='store_true',
                      help='Retry proteins that previously failed')
    parser.add_argument('--force', action='store_true',
                      help='Force regeneration of existing summaries')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--representative-only', action='store_true',
                  help='Only process representative proteins')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.batch_summary")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get proteins to process
    batch_info, proteins = get_proteins_to_process(
        context, 
        args.batch_id, 
        args.limit, 
        args.retry_failed,
        args.force,
        args.representative_only
    )
    
    if not batch_info:
        return 1
        
    if not proteins:
        return 0
    
    logger.info(f"Processing batch {batch_info['id']} ({batch_info['name']})")
    logger.info(f"Found {len(proteins)} proteins to process")
    
    # Process proteins
    results = []
    
    if args.threads > 1 and len(proteins) > 1:
        logger.info(f"Processing with {args.threads} threads")
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            # Create separate context for each thread to avoid concurrent access issues
            futures = {}
            for protein in proteins:
                # Create a new context for each protein to avoid thread safety issues
                thread_context = ApplicationContext(args.config)
                if args.force:
                    thread_context.set_force_overwrite(True)
                
                future = executor.submit(
                    process_protein, 
                    protein, 
                    thread_context, 
                    args.output_dir, 
                    args.blast_only,
                    args.force,
                    args.verbose
                )
                futures[future] = protein
            
            for future in as_completed(futures):
                protein = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    status = 'SUCCESS' if result['success'] else 'FAILED'
                    logger.info(f"Completed {result['pdb_id']}_{result['chain_id']}: {status} ({result['duration']:.2f}s)")
                except Exception as e:
                    logger.error(f"Error processing {protein['pdb_id']}_{protein['chain_id']}: {str(e)}")
    else:
        logger.info("Processing proteins sequentially")
        # Set force flag on context if specified
        if args.force:
            context.set_force_overwrite(True)
            
        for protein in proteins:
            result = process_protein(
                protein, 
                context, 
                args.output_dir, 
                args.blast_only,
                args.force,
                args.verbose
            )
            results.append(result)
            status = 'SUCCESS' if result['success'] else 'FAILED'
            logger.info(f"Completed {result['pdb_id']}_{result['chain_id']}: {status} ({result['duration']:.2f}s)")
    
    # Summarize results
    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]
    
    logger.info(f"Batch processing complete: {len(successful)} successful, {len(failed)} failed")
    
    # Update batch completion in database
    if successful:
        update_batch_completion(context, args.batch_id, len(successful), logger)
    
    # List failed proteins if any
    if failed:
        logger.warning("Failed proteins:")
        for f in failed:
            logger.warning(f"  {f['pdb_id']}_{f['chain_id']}: {f['error']}")
    
    return 0 if not failed else 1

if __name__ == "__main__":
    sys.exit(main())