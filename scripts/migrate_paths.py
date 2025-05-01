#!/usr/bin/env python3
"""
migrate_paths.py - Migrate files to standardized path structure
"""

import os
import sys
import argparse
import logging
from typing import Dict, List, Tuple

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.utils.path_utils import (
    get_standardized_paths,
    find_files_with_legacy_paths,
    migrate_file_to_standard_path,
    update_db_file_paths,
    get_file_db_path
)

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

def get_batch_info(db, batch_id):
    """Get batch information"""
    query = """
    SELECT id, batch_name, base_path, ref_version 
    FROM ecod_schema.batch 
    WHERE id = %s
    """
    
    batch_info = db.execute_dict_query(query, (batch_id,))
    return batch_info[0] if batch_info else None

def get_protein_chains(db, batch_id, limit=None):
    """Get protein chains for a batch"""
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, ps.id as process_id
    FROM 
        ecod_schema.protein p
    JOIN
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        ps.batch_id = %s
    ORDER BY
        p.pdb_id, p.chain_id
    """
    
    if limit:
        query += f" LIMIT {limit}"
    
    return db.execute_dict_query(query, (batch_id,))

def get_process_files(db, process_id):
    """Get process files"""
    query = """
    SELECT 
        id, file_type, file_path, file_exists, file_size
    FROM 
        ecod_schema.process_file
    WHERE 
        process_id = %s
    """
    
    return db.execute_dict_query(query, (process_id,))

def migrate_batch_files(db, batch_id, dry_run=False, limit=None):
    """Migrate files for a batch to standardized paths"""
    logger = logging.getLogger("migration")
    
    # Get batch info
    batch_info = get_batch_info(db, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return 0
    
    batch_path = batch_info['base_path']
    ref_version = batch_info['ref_version']
    
    logger.info(f"Migrating files for batch {batch_id} ({batch_info['batch_name']})")
    logger.info(f"Base path: {batch_path}")
    logger.info(f"Reference version: {ref_version}")
    
    # Get protein chains
    chains = get_protein_chains(db, batch_id, limit)
    logger.info(f"Found {len(chains)} protein chains")
    
    # Statistics
    stats = {
        'total_files': 0,
        'already_standard': 0,
        'migrated': 0,
        'failed': 0,
        'no_file': 0
    }
    
    # Process each chain
    for chain in chains:
        pdb_id = chain['pdb_id']
        chain_id = chain['chain_id']
        process_id = chain['process_id']
        
        logger.info(f"Processing {pdb_id}_{chain_id} (process ID: {process_id})")
        
        # Get standard paths
        standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=not dry_run)
        
        # Get existing files
        existing_files = find_files_with_legacy_paths(batch_path, pdb_id, chain_id, ref_version)
        
        # Get process files from database
        process_files = get_process_files(db, process_id)
        process_files_by_type = {file['file_type']: file for file in process_files}
        
        # Process each file type
        for file_type, paths in existing_files.items():
            standard_path = paths['standard']
            legacy_path = paths['legacy']
            exists_at = paths['exists_at']
            
            stats['total_files'] += 1
            
            # Skip if file doesn't exist
            if not exists_at:
                logger.debug(f"File not found for {pdb_id}_{chain_id}: {file_type}")
                stats['no_file'] += 1
                continue
            
            # Already at standard path
            if exists_at == standard_path:
                logger.debug(f"Already at standard path: {standard_path}")
                stats['already_standard'] += 1
                continue
            
            # Need to migrate
            logger.info(f"Need to migrate: {exists_at} -> {standard_path}")
            
            # In dry run mode, just report
            if dry_run:
                logger.info(f"Would migrate: {exists_at} -> {standard_path}")
                stats['migrated'] += 1
                continue
            
            # Migrate file
            success = migrate_file_to_standard_path(exists_at, standard_path)
            
            if success:
                stats['migrated'] += 1
                
                # Update database if this file type is in the database
                if file_type in process_files_by_type:
                    db_file = process_files_by_type[file_type]
                    new_rel_path = get_file_db_path(batch_path, standard_path)
                    
                    if not dry_run:
                        update_db_file_path(db, db_file['id'], new_rel_path)
            else:
                stats['failed'] += 1
    
    # Print statistics
    logger.info("Migration Statistics:")
    logger.info(f"Total files processed: {stats['total_files']}")
    logger.info(f"Already at standard path: {stats['already_standard']}")
    logger.info(f"Files migrated: {stats['migrated']}")
    logger.info(f"Migration failed: {stats['failed']}")
    logger.info(f"Files not found: {stats['no_file']}")
    
    return stats['migrated']

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Migrate files to standardized path structure')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--dry-run', action='store_true',
                      help='Perform a dry run (no actual changes)')
    parser.add_argument('--limit', type=int,
                      help='Limit the number of chains to process')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    logger = logging.getLogger("main")
    mode = "dry run" if args.dry_run else "actual migration"
    logger.info(f"Starting path migration ({mode}) for batch {args.batch_id}")
    
    # Initialize configuration and database
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    migrated_count = migrate_batch_files(db, args.batch_id, args.dry_run, args.limit)
    
    if args.dry_run:
        logger.info(f"Dry run complete. {migrated_count} files would be migrated.")
    else:
        logger.info(f"Migration complete. {migrated_count} files were migrated.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
