#!/usr/bin/env python3
"""
migrate_domain_files.py - Standardize domain file paths and update database records

This script specifically handles the domain file standardization:
1. Finds all .domains.xml files and converts them to the standard .domain_partition.xml format
2. Updates database records to reflect the new paths
3. Updates process_status for chains that have domain files but aren't marked complete

Usage:
  python migrate_domain_files.py --batch-id 68 --config config/config.yml
  python migrate_domain_files.py --batch-id 68 --config config/config.yml --dry-run
"""

import os
import sys
import glob
import shutil
import argparse
import logging
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.utils.path_utils import get_standardized_paths, get_file_db_path

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

def find_domain_files(batch_path, pattern="*.domains.xml"):
    """Find all domain files in the batch directory"""
    domains_dir = os.path.join(batch_path, "domains")
    if not os.path.exists(domains_dir):
        return []
    
    return glob.glob(os.path.join(domains_dir, pattern))

def parse_filename(filename):
    """Parse domain filename to extract PDB ID and chain ID"""
    basename = os.path.basename(filename)
    parts = basename.split(".")
    
    if len(parts) < 3:
        return None, None, None
    
    pdb_chain = parts[0].split("_")
    ref_version = parts[1] if len(parts) > 1 else None
    
    if len(pdb_chain) != 2:
        return None, None, None
    
    return pdb_chain[0], pdb_chain[1], ref_version

def get_standard_domain_path(batch_path, pdb_id, chain_id, ref_version):
    """Get the standardized path for a domain partition file"""
    domains_dir = os.path.join(batch_path, "domains")
    return os.path.join(domains_dir, f"{pdb_id}_{chain_id}.{ref_version}.domain_partition.xml")

def get_protein_info(db, pdb_id, chain_id, batch_id):
    """Get protein and process information from the database"""
    query = """
    SELECT 
        p.id as protein_id, ps.id as process_id, ps.current_stage
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE 
        p.pdb_id = %s AND p.chain_id = %s AND ps.batch_id = %s
    """
    
    results = db.execute_dict_query(query, (pdb_id, chain_id, batch_id))
    return results[0] if results else None

def get_domain_partition_file(db, process_id):
    """Get domain partition file information from the database"""
    query = """
    SELECT 
        id, file_path, file_exists
    FROM 
        ecod_schema.process_file
    WHERE 
        process_id = %s AND file_type = 'domain_partition'
    """
    
    results = db.execute_dict_query(query, (process_id,))
    return results[0] if results else None

def update_process_status(db, process_id, new_stage, error_message=None, dry_run=False):
    """Update process status in the database"""
    if dry_run:
        return True
        
    update_data = {"current_stage": new_stage}
    if error_message is not None:
        update_data["error_message"] = error_message
    
    try:
        db.update(
            "ecod_schema.process_status",
            update_data,
            "id = %s",
            (process_id,)
        )
        return True
    except Exception as e:
        logging.error(f"Failed to update process status: {str(e)}")
        return False

def update_file_record(db, file_id, new_path, file_exists=True, dry_run=False):
    """Update file record in the database"""
    if dry_run:
        return True
        
    try:
        db.update(
            "ecod_schema.process_file",
            {"file_path": new_path, "file_exists": file_exists},
            "id = %s",
            (file_id,)
        )
        return True
    except Exception as e:
        logging.error(f"Failed to update file record: {str(e)}")
        return False

def create_file_record(db, process_id, file_type, file_path, file_exists=True, dry_run=False):
    """Create a new file record in the database"""
    if dry_run:
        return None
        
    try:
        file_id = db.insert(
            "ecod_schema.process_file",
            {
                "process_id": process_id,
                "file_type": file_type,
                "file_path": file_path,
                "file_exists": file_exists
            },
            returning="id"
        )
        return file_id
    except Exception as e:
        logging.error(f"Failed to create file record: {str(e)}")
        return None

def validate_domain_file(file_path):
    """Validate that the domain file has the expected structure"""
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Check for domain_doc root
        if root.tag != "domain_doc":
            return False, "Root element is not domain_doc"
        
        # Check for domain_list and at least one domain
        domain_list = root.find("domain_list")
        if domain_list is None:
            return False, "No domain_list element found"
        
        domains = domain_list.findall("domain")
        if not domains:
            return False, "No domains found in domain_list"
        
        return True, f"Valid domain file with {len(domains)} domains"
    except ET.ParseError as e:
        return False, f"XML parsing error: {str(e)}"
    except Exception as e:
        return False, f"Validation error: {str(e)}"

def migrate_batch_domain_files(db, batch_id, dry_run=False):
    """Migrate domain files for a batch to standardized paths and update database records"""
    logger = logging.getLogger("migration")
    
    # Get batch info
    batch_info = get_batch_info(db, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return 0
    
    batch_path = batch_info['base_path']
    ref_version = batch_info['ref_version']
    
    logger.info(f"Migrating domain files for batch {batch_id} ({batch_info['batch_name']})")
    logger.info(f"Base path: {batch_path}")
    logger.info(f"Reference version: {ref_version}")
    
    # Find all domain files
    domain_files = find_domain_files(batch_path)
    logger.info(f"Found {len(domain_files)} domain files")
    
    # Statistics
    stats = {
        'total_files': len(domain_files),
        'migrated': 0,
        'already_standard': 0,
        'failed': 0,
        'db_updated': 0,
        'status_updated': 0,
        'errors': []
    }
    
    # Process each file
    for domain_file in domain_files:
        # Parse filename
        pdb_id, chain_id, file_ref_version = parse_filename(domain_file)
        
        if not pdb_id or not chain_id:
            logger.warning(f"Could not parse filename: {domain_file}")
            stats['failed'] += 1
            stats['errors'].append(f"Failed to parse: {os.path.basename(domain_file)}")
            continue
        
        # Use file-specific ref version if available, otherwise batch ref version
        current_ref_version = file_ref_version or ref_version
        
        logger.info(f"Processing domain file for {pdb_id}_{chain_id}")
        
        # Get standard path
        standard_path = get_standard_domain_path(batch_path, pdb_id, chain_id, current_ref_version)
        
        # Check if already at standard path
        if domain_file == standard_path:
            logger.info(f"Already at standard path: {domain_file}")
            stats['already_standard'] += 1
            continue
        
        # Validate domain file
        valid, message = validate_domain_file(domain_file)
        if not valid:
            logger.warning(f"Invalid domain file {domain_file}: {message}")
            stats['failed'] += 1
            stats['errors'].append(f"Invalid domain file: {os.path.basename(domain_file)} - {message}")
            continue
        else:
            logger.debug(f"Validated domain file: {message}")
        
        # Migrate file
        if not dry_run:
            try:
                # Create domains directory if it doesn't exist
                os.makedirs(os.path.dirname(standard_path), exist_ok=True)
                
                # Copy file to standard path
                shutil.copy2(domain_file, standard_path)
                logger.info(f"Migrated file: {domain_file} -> {standard_path}")
                
                # Optionally remove old file if different location
                # if os.path.normpath(domain_file) != os.path.normpath(standard_path):
                #     os.remove(domain_file)
                #     logger.info(f"Removed old file: {domain_file}")
            except Exception as e:
                logger.error(f"Failed to migrate file {domain_file}: {str(e)}")
                stats['failed'] += 1
                stats['errors'].append(f"Migration failed: {os.path.basename(domain_file)} - {str(e)}")
                continue
        else:
            logger.info(f"Would migrate: {domain_file} -> {standard_path}")
        
        stats['migrated'] += 1
        
        # Update database
        # Get protein info
        protein_info = get_protein_info(db, pdb_id, chain_id, batch_id)
        
        if not protein_info:
            logger.warning(f"Protein {pdb_id}_{chain_id} not found in database for batch {batch_id}")
            continue
        
        process_id = protein_info['process_id']
        current_stage = protein_info['current_stage']
        
        # Update domain partition file record
        rel_path = get_file_db_path(batch_path, standard_path)
        
        domain_file_record = get_domain_partition_file(db, process_id)
        
        if domain_file_record:
            # Update existing record
            file_updated = update_file_record(
                db, domain_file_record['id'], rel_path, True, dry_run
            )
            if file_updated:
                logger.info(f"Updated domain partition file record for {pdb_id}_{chain_id}")
                stats['db_updated'] += 1
        else:
            # Create new record
            file_id = create_file_record(
                db, process_id, 'domain_partition', rel_path, True, dry_run
            )
            if file_id:
                logger.info(f"Created domain partition file record for {pdb_id}_{chain_id}")
                stats['db_updated'] += 1
        
        # Update process status if not already complete
        if current_stage != 'domain_partition_complete':
            status_updated = update_process_status(
                db, process_id, 'domain_partition_complete', None, dry_run
            )
            if status_updated:
                logger.info(f"Updated process status for {pdb_id}_{chain_id} to domain_partition_complete")
                stats['status_updated'] += 1
    
    # Print statistics
    logger.info("\nMigration Statistics:")
    logger.info(f"Total domain files: {stats['total_files']}")
    logger.info(f"Already at standard path: {stats['already_standard']}")
    logger.info(f"Files migrated: {stats['migrated']}")
    logger.info(f"Migration failed: {stats['failed']}")
    logger.info(f"Database records updated: {stats['db_updated']}")
    logger.info(f"Process status updated: {stats['status_updated']}")
    
    if stats['errors']:
        logger.info(f"\nEncountered {len(stats['errors'])} errors:")
        for error in stats['errors'][:10]:  # Show first 10 errors
            logger.info(f"  - {error}")
        if len(stats['errors']) > 10:
            logger.info(f"  ... and {len(stats['errors']) - 10} more errors")
    
    return stats

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Migrate domain files to standardized paths and update database records')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--dry-run', action='store_true',
                      help='Perform a dry run (no actual changes)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    logger = logging.getLogger("main")
    mode = "dry run" if args.dry_run else "actual migration"
    logger.info(f"Starting domain file migration ({mode}) for batch {args.batch_id}")
    
    # Initialize configuration and database
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    stats = migrate_batch_domain_files(db, args.batch_id, args.dry_run)
    
    if args.dry_run:
        logger.info(f"Dry run complete. {stats['migrated']} files would be migrated.")
    else:
        logger.info(f"Migration complete. {stats['migrated']} files were migrated.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
