#!/usr/bin/env python3
"""
migrate_domain_summaries.py - Migrate domain summary files to new location
"""

import os
import sys
import argparse
import shutil
import logging
from tqdm import tqdm

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.error_handlers import handle_exceptions

def setup_logging(verbose=False, log_file=None):
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
    parser = argparse.ArgumentParser(description='Migrate domain summary files to new location')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, default=None,
                      help='Specific batch ID to migrate (default: all batches)')
    parser.add_argument('--dry-run', action='store_true',
                      help='Perform a dry run without making changes')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    logger = logging.getLogger("ecod.migration")
    
    # Initialize context
    context = ApplicationContext(args.config)
    
    # Get batches to process
    if args.batch_id:
        batch_query = "SELECT id, base_path FROM ecod_schema.batch WHERE id = %s"
        batches = context.db.execute_dict_query(batch_query, (args.batch_id,))
    else:
        batch_query = "SELECT id, base_path FROM ecod_schema.batch"
        batches = context.db.execute_dict_query(batch_query)
    
    if not batches:
        logger.error(f"No batches found{' with ID ' + str(args.batch_id) if args.batch_id else ''}")
        return 1
    
    total_files = 0
    moved_files = 0
    updated_records = 0
    
    for batch in batches:
        batch_id = batch['id']
        base_path = batch['base_path']
        
        logger.info(f"Processing batch {batch_id} at {base_path}")
        
        # Get domain summary files for this batch
        query = """
        SELECT pf.id, pf.process_id, pf.file_path, p.pdb_id, p.chain_id
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.batch_id = %s AND pf.file_type = 'domain_summary'
        """
        
        summary_files = context.db.execute_dict_query(query, (batch_id,))
        
        logger.info(f"Found {len(summary_files)} domain summary files in batch {batch_id}")
        
        # Create domains directory if it doesn't exist
        domains_dir = os.path.join(base_path, "domains")
        os.makedirs(domains_dir, exist_ok=True)
        
        # Process each file
        for file_info in tqdm(summary_files):
            total_files += 1
            
            file_id = file_info['id']
            process_id = file_info['process_id']
            file_path = file_info['file_path']
            pdb_id = file_info['pdb_id']
            chain_id = file_info['chain_id']
            
            # Skip if already in domains directory
            if file_path.startswith("domains/"):
                logger.debug(f"File already in domains directory: {file_path}")
                continue
            
            # Get full path to current file
            current_full_path = os.path.join(base_path, file_path)
            
            # Check if file exists at old location
            if not os.path.exists(current_full_path):
                logger.warning(f"File not found at {current_full_path}, skipping")
                continue
            
            # Extract filename from path
            filename = os.path.basename(file_path)
            
            # Construct new path
            new_relative_path = os.path.join("domains", filename)
            new_full_path = os.path.join(base_path, new_relative_path)
            
            # Move file if not dry run
            if not args.dry_run:
                try:
                    # Copy file to new location
                    shutil.copy2(current_full_path, new_full_path)
                    
                    # Verify copy succeeded
                    if os.path.exists(new_full_path):
                        # Update database record
                        context.db.update(
                            "ecod_schema.process_file",
                            {"file_path": new_relative_path},
                            "id = %s",
                            (file_id,)
                        )
                        
                        # Only remove original after successful copy and DB update
                        os.remove(current_full_path)
                        
                        moved_files += 1
                        updated_records += 1
                        logger.debug(f"Moved {current_full_path} to {new_full_path}")
                    else:
                        logger.error(f"Failed to copy file to {new_full_path}")
                except Exception as e:
                    logger.error(f"Error processing {current_full_path}: {str(e)}")
            else:
                logger.info(f"[DRY RUN] Would move {current_full_path} to {new_full_path}")
                moved_files += 1
                updated_records += 1
    
    logger.info(f"Migration {'' if not args.dry_run else 'would have'} processed {total_files} files")
    logger.info(f"Moved {moved_files} files to new location")
    logger.info(f"Updated {updated_records} database records")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())