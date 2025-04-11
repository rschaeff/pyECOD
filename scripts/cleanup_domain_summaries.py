#!/usr/bin/env python3
"""
cleanup_domain_summaries.py - Find and relocate domain summaries to standardized locations
"""

import os
import sys
import logging
import shutil
import argparse
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

def find_all_domain_summaries(base_dir):
    """Find all domain summary files in any location"""
    # Look for any XML files that could be domain summaries
    summaries = []
    
    # Walk the entire directory tree
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            # Check for various possible domain summary naming patterns
            if file.endswith('.domain_summary.xml') or '.blast_summ.' in file:
                full_path = os.path.join(root, file)
                # Parse the PDB and chain ID from the filename
                if '_' in file:
                    parts = file.split('_', 1)
                    pdb_id = parts[0]
                    
                    # Extract chain ID (handle various naming patterns)
                    chain_part = parts[1]
                    if '.' in chain_part:
                        chain_id = chain_part.split('.')[0]
                    else:
                        chain_id = chain_part
                    
                    summaries.append({
                        'path': full_path,
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'filename': file
                    })
    
    return summaries

def standardize_summary_locations(batch_id, config_path, dry_run=True):
    """Find and standardize domain summary locations"""
    context = ApplicationContext(config_path)
    logger = logging.getLogger("ecod.standardize_summaries")
    
    # Get batch path
    batch_query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return 1
        
    batch_path = batch_result[0][0]
    domains_dir = os.path.join(batch_path, "domains")
    os.makedirs(domains_dir, exist_ok=True)
    
    # Find all domain summaries
    all_summaries = find_all_domain_summaries(batch_path)
    logger.info(f"Found {len(all_summaries)} domain summary files in total")
    
    # Count how many are already in the correct location
    correct_location = 0
    for summary in all_summaries:
        if os.path.dirname(summary['path']) == domains_dir:
            correct_location += 1
    
    logger.info(f"{correct_location} files already in correct location")
    logger.info(f"{len(all_summaries) - correct_location} files need to be moved")
    
    # Process each file that needs to be moved
    moved_count = 0
    errors = 0
    
    for summary in all_summaries:
        if os.path.dirname(summary['path']) == domains_dir:
            continue  # Already in the right place
        
        # Define the standardized filename and path
        std_filename = f"{summary['pdb_id']}_{summary['chain_id']}.domain_summary.xml"
        std_path = os.path.join(domains_dir, std_filename)
        
        # Check if destination already exists
        if os.path.exists(std_path):
            logger.warning(f"Destination already exists: {std_path}, skipping")
            continue
        
        logger.info(f"Moving {summary['path']} -> {std_path}")
        
        if not dry_run:
            try:
                # Copy file to new location
                shutil.copy2(summary['path'], std_path)
                
                # Update database record
                process_query = """
                SELECT ps.id 
                FROM ecod_schema.process_status ps
                JOIN ecod_schema.protein p ON ps.protein_id = p.id
                WHERE ps.batch_id = %s AND p.pdb_id = %s AND p.chain_id = %s
                """
                process_result = context.db.execute_query(process_query, (batch_id, summary['pdb_id'], summary['chain_id']))
                
                if process_result:
                    process_id = process_result[0][0]
                    
                    # Check if a domain_summary record exists
                    file_query = """
                    SELECT id FROM ecod_schema.process_file
                    WHERE process_id = %s AND file_type = 'domain_summary'
                    """
                    file_result = context.db.execute_query(file_query, (process_id,))
                    
                    rel_path = os.path.relpath(std_path, batch_path)
                    
                    if file_result:
                        # Update existing record
                        file_id = file_result[0][0]
                        update_query = """
                        UPDATE ecod_schema.process_file
                        SET file_path = %s, file_exists = TRUE,
                            file_size = %s, last_checked = CURRENT_TIMESTAMP
                        WHERE id = %s
                        """
                        context.db.execute_query(update_query, (rel_path, os.path.getsize(std_path), file_id))
                    else:
                        # Insert new record
                        insert_query = """
                        INSERT INTO ecod_schema.process_file
                        (process_id, file_type, file_path, file_exists, file_size)
                        VALUES (%s, 'domain_summary', %s, TRUE, %s)
                        """
                        context.db.execute_query(insert_query, (process_id, rel_path, os.path.getsize(std_path)))
                
                moved_count += 1
            except Exception as e:
                logger.error(f"Error processing {summary['path']}: {str(e)}")
                errors += 1
    
    logger.info(f"Summary: {moved_count} files moved, {errors} errors")
    
    if dry_run:
        logger.info("DRY RUN - No files were actually moved")
    
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Standardize domain summary file locations')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--apply', action='store_true',
                      help='Actually move files (default is dry run)')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    sys.exit(standardize_summary_locations(
        args.batch_id, 
        args.config, 
        dry_run=not args.apply
    ))