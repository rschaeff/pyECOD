#!/usr/bin/env python3
"""
generate_batch_domain_summaries.py - Generate domain summaries for all proteins in a batch
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

def process_protein(protein_info: Dict[str, Any], context: ApplicationContext, output_dir: str, blast_only: bool, verbose: bool) -> Dict[str, Any]:
    """Process a single protein for domain summary generation"""
    logger = logging.getLogger(f"ecod.generate_summary.{protein_info['pdb_id']}_{protein_info['chain_id']}")
    
    protein_id = protein_info['id']
    pdb_id = protein_info['pdb_id']
    chain_id = protein_info['chain_id']
    batch_id = protein_info['batch_id']
    batch_path = protein_info['batch_path']
    reference = protein_info['reference']
    
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
        
        # Initialize domain summary processor
        summary_processor = DomainSummary(context.config_manager.config_path)
        
        # Use provided output dir or batch path
        effective_output_dir = output_dir if output_dir else batch_path
        
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
            process_query = """
            SELECT id FROM ecod_schema.process_status
            WHERE protein_id = %s AND batch_id = %s
            """
            process_result = context.db.execute_query(process_query, (protein_id, batch_id))
            
            if process_result:
                process_id = process_result[0][0]
                
                # Check if summary file record already exists
                check_query = """
                SELECT id FROM ecod_schema.process_file
                WHERE process_id = %s AND file_type = 'domain_summary'
                """
                
                existing = context.db.execute_query(check_query, (process_id,))
                
                if existing:
                    # Update existing record
                    update_query = """
                    UPDATE ecod_schema.process_file
                    SET file_path = %s, file_exists = TRUE, file_size = %s, last_checked = NOW()
                    WHERE id = %s
                    """
                    
                    rel_path = os.path.relpath(summary_file, batch_path)
                    file_size = os.path.getsize(summary_file)
                    
                    context.db.execute_query(update_query, (rel_path, file_size, existing[0][0]))
                    logger.info(f"Updated domain_summary file record in database")
                else:
                    # Create new record
                    insert_query = """
                    INSERT INTO ecod_schema.process_file
                    (process_id, file_type, file_path, file_exists, file_size, last_checked)
                    VALUES (%s, %s, %s, %s, %s, NOW())
                    """
                    
                    rel_path = os.path.relpath(summary_file, batch_path)
                    file_size = os.path.getsize(summary_file)
                    
                    context.db.execute_query(insert_query, (process_id, 'domain_summary', rel_path, True, file_size))
                    logger.info(f"Added domain_summary file record to database")
                
                # Update process status
                status_query = """
                UPDATE ecod_schema.process_status
                SET current_stage = 'domain_summary', status = 'success', updated_at = NOW()
                WHERE id = %s
                """
                
                context.db.execute_query(status_query, (process_id,))
                logger.info(f"Updated process status to domain_summary:success")
        else:
            error_msg = "Failed to generate domain summary"
            logger.error(error_msg)
            result['error'] = error_msg
            
    except Exception as e:
        error_msg = f"Error generating domain summary: {str(e)}"
        logger.error(error_msg, exc_info=True)
        result['error'] = error_msg
    
    end_time = time.time()
    result['duration'] = end_time - start_time
    
    return result

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
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.batch_summary")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch information
    batch_query = """
    SELECT 
        id, batch_name, base_path, ref_version
    FROM 
        ecod_schema.batch
    WHERE 
        id = %s
    """
    
    batch_result = context.db.execute_query(batch_query, (args.batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
    
    batch_info = {
        'id': batch_result[0][0],
        'name': batch_result[0][1],
        'path': batch_result[0][2],
        'reference': batch_result[0][3]
    }
    
    logger.info(f"Processing batch {batch_info['id']} ({batch_info['name']})")
    
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
    
    if not args.retry_failed:
        protein_query += " AND ps.status != 'error'"
    
    if not args.force:
        protein_query += " AND (pf.id IS NULL OR pf.file_exists = FALSE)"
    
    protein_query += " ORDER BY p.id"
    
    if args.limit:
        protein_query += f" LIMIT {args.limit}"
    
    protein_results = context.db.execute_query(protein_query, (args.batch_id,))
    
    if not protein_results:
        logger.info(f"No proteins found to process in batch {args.batch_id}")
        return 0
    
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
    
    logger.info(f"Found {len(proteins)} proteins to process")
    
    # Process proteins (in parallel if threads > 1)
    results = []
    
    if args.threads > 1:
        logger.info(f"Processing with {args.threads} threads")
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = {
                executor.submit(
                    process_protein, 
                    protein, 
                    context, 
                    args.output_dir, 
                    args.blast_only, 
                    args.verbose
                ): protein for protein in proteins
            }
            
            for future in as_completed(futures):
                protein = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                    logger.info(f"Completed {result['pdb_id']}_{result['chain_id']}: {'SUCCESS' if result['success'] else 'FAILED'} ({result['duration']:.2f}s)")
                except Exception as e:
                    logger.error(f"Error processing {protein['pdb_id']}_{protein['chain_id']}: {str(e)}")
    else:
        logger.info("Processing proteins sequentially")
        for protein in proteins:
            result = process_protein(protein, context, args.output_dir, args.blast_only, args.verbose)
            results.append(result)
            logger.info(f"Completed {result['pdb_id']}_{result['chain_id']}: {'SUCCESS' if result['success'] else 'FAILED'} ({result['duration']:.2f}s)")
    
    # Summarize results
    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]
    
    logger.info(f"Batch processing complete: {len(successful)} successful, {len(failed)} failed")
    
    # Update batch completion in database
    if len(successful) > 0:
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
        
        context.db.execute_query(update_query, (len(successful), len(successful), len(successful), args.batch_id))
        logger.info(f"Updated batch {args.batch_id} completion status")
    
    # List failed proteins if any
    if failed:
        logger.warning("Failed proteins:")
        for f in failed:
            logger.warning(f"  {f['pdb_id']}_{f['chain_id']}: {f['error']}")
    
    return 0 if not failed else 1

if __name__ == "__main__":
    sys.exit(main())