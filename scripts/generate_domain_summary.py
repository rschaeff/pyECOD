#!/usr/bin/env python3
"""
generate_domain_summary.py - Generate domain summary for a single protein
"""

import os
import sys
import logging
import argparse
from typing import Dict, Any, Optional

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

def patch_domain_summary(domain_summary_obj):
    """
    Patch the DomainSummary class to fix file type mismatches and path issues
    """
    original_method = domain_summary_obj.simplified_file_path_resolution
    
    def patched_method(pdb_id, chain_id, file_type, job_dump_dir):
        # Fix file type mismatch
        logger = logging.getLogger("ecod.pipelines.domain_analysis.summary")
        logger.debug(f"Original file_type requested: {file_type}")
        
        # Map file types
        file_type_map = {
            'domain_blast_results': 'domain_blast_result',
            'chain_blast_result': 'blast_result'  # Map 'chain_blast_result' to 'blast_result'
        }
        
        if file_type in file_type_map:
            mapped_type = file_type_map[file_type]
            logger.debug(f"Remapping file type from '{file_type}' to '{mapped_type}'")
            file_type = mapped_type
        
        # Call original method with remapped file type
        paths = original_method(pdb_id, chain_id, file_type, job_dump_dir)
        
        # Normalize paths if needed
        fixed_paths = []
        for path in paths:
            if '..' in path:
                normalized = os.path.normpath(path)
                logger.debug(f"Normalizing path: {path} -> {normalized}")
                fixed_paths.append(normalized)
            else:
                fixed_paths.append(path)
                
        logger.debug(f"Resolved paths for {pdb_id}_{chain_id}, type '{file_type}': {fixed_paths}")
        return fixed_paths
    
    # Replace the method with our patched version
    domain_summary_obj.simplified_file_path_resolution = patched_method
    
def main():
    """Main function to generate domain summary for a single protein"""
    parser = argparse.ArgumentParser(description='Generate domain summary for a single protein')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID containing the protein')
    parser.add_argument('--protein-id', type=int, required=True,
                      help='Protein ID to process')
    parser.add_argument('--blast-only', action='store_true',
                      help='Generate blast-only summary (no HHSearch)')
    parser.add_argument('--output-dir', type=str,
                      help='Override output directory (default: from batch path)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.generate_summary")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get protein and batch information
    query = """
    SELECT 
        p.id, p.source_id, p.pdb_id, p.chain_id,
        b.id as batch_id, b.batch_name, b.base_path, b.ref_version
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    WHERE 
        p.id = %s AND b.id = %s
    """
    
    result = context.db.execute_query(query, (args.protein_id, args.batch_id))
    
    if not result:
        logger.error(f"Protein {args.protein_id} not found in batch {args.batch_id}")
        return 1
    
    protein_info = result[0]
    pdb_id = protein_info[2]
    chain_id = protein_info[3]
    batch_path = protein_info[6]
    reference = protein_info[7]
    
    logger.info(f"Generating domain summary for protein: {pdb_id}_{chain_id} (ID: {args.protein_id})")
    
    # Initialize domain summary processor
    summary_processor = DomainSummary(args.config)
    
    # Apply patch to fix file type mismatch
    patch_domain_summary(summary_processor)
    logger.info("Applied patch to fix domain_blast_result/domain_blast_results file type mismatch")
    
    # Determine output directory
    output_dir = args.output_dir if args.output_dir else batch_path
    
    # Generate domain summary
    try:
        summary_file = summary_processor.create_summary(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            job_dump_dir=output_dir,
            blast_only=args.blast_only
        )
        
        if summary_file:
            logger.info(f"Successfully generated domain summary: {summary_file}")
            
            # Register the summary file in the database
            process_query = """
            SELECT id FROM ecod_schema.process_status
            WHERE protein_id = %s AND batch_id = %s
            """
            process_result = context.db.execute_query(process_query, (args.protein_id, args.batch_id))
            
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
            
            return 0
        else:
            logger.error(f"Failed to generate domain summary")
            return 1
            
    except Exception as e:
        logger.error(f"Error generating domain summary: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())