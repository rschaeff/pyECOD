#!/usr/bin/env python3
"""
register_batch_files.py - Find and register existing files that aren't properly
registered in the database.

This script addresses validation issues where files exist on disk but aren't
properly registered in the database, causing validation scripts to report them
as missing. It scans the batch directory, identifies existing files, and registers
them correctly in the database.
"""

import os
import sys
import logging
import argparse
import glob
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Set, Tuple
import xml.etree.ElementTree as ET

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    find_files_with_legacy_paths,
    extract_pdb_chain_from_path,
    scan_batch_directory
)

def setup_logging(verbose: bool = False, log_file: Optional[str] = None) -> None:
    """Configure logging with appropriate handlers and format"""
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

def get_batch_info(db, batch_id: int) -> Dict[str, Any]:
    """Get information about a batch"""
    logger = logging.getLogger("register_files")

    query = """
    SELECT id, batch_name, base_path, ref_version 
    FROM ecod_schema.batch 
    WHERE id = %s
    """
    
    results = db.execute_dict_query(query, (batch_id,))
    
    if not results:
        logger.error(f"Batch {batch_id} not found")
        return {}
        
    return results[0]

def get_batch_proteins(db, batch_id: int) -> List[Dict[str, Any]]:
    """Get all proteins in a batch"""
    logger = logging.getLogger("register_files")

    query = """
    SELECT 
        ps.id as process_id,
        p.id as protein_id,
        p.pdb_id,
        p.chain_id,
        ps.current_stage,
        ps.is_representative
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
    ORDER BY 
        p.pdb_id, p.chain_id
    """
    
    proteins = db.execute_dict_query(query, (batch_id,))
    
    logger.info(f"Found {len(proteins)} proteins in batch {batch_id}")
    
    return proteins

def get_existing_file_registrations(db, process_ids: List[int]) -> Dict[int, Dict[str, Dict[str, Any]]]:
    """Get existing file registrations for a list of process_ids"""
    logger = logging.getLogger("register_files")
    
    if not process_ids:
        return {}
        
    query = """
    SELECT 
        process_id, file_type, file_path, file_exists
    FROM 
        ecod_schema.process_file 
    WHERE 
        process_id = ANY(%s)
    """
    
    files = db.execute_dict_query(query, (process_ids,))
    
    logger.info(f"Found {len(files)} existing file registrations")
    
    # Organize by process_id and file_type
    result = {}
    for file in files:
        process_id = file['process_id']
        file_type = file['file_type']
        
        if process_id not in result:
            result[process_id] = {}
            
        result[process_id][file_type] = file
        
    return result

def find_blast_files(batch_path: str, pdb_id: str, chain_id: str, ref_version: str) -> Dict[str, List[str]]:
    """Find BLAST files for a protein using various path patterns"""
    logger = logging.getLogger("register_files")
    
    # Common patterns for BLAST files
    chain_blast_patterns = [
        f"{batch_path}/blast/chain/{pdb_id}_{chain_id}.{ref_version}.xml",
        f"{batch_path}/blast/chain/{pdb_id}_{chain_id}*.xml",
        f"{batch_path}/blast/{pdb_id}_{chain_id}.{ref_version}.chainwise_blast.xml",
        f"{batch_path}/blast/{pdb_id}_{chain_id}.{ref_version}.chain_blast.xml",
        f"{batch_path}/blast/chain/batch_*/{pdb_id}_{chain_id}*.xml"
    ]
    
    domain_blast_patterns = [
        f"{batch_path}/blast/domain/{pdb_id}_{chain_id}.{ref_version}.xml",
        f"{batch_path}/blast/domain/{pdb_id}_{chain_id}*.xml",
        f"{batch_path}/blast/{pdb_id}_{chain_id}.{ref_version}.blast.xml",
        f"{batch_path}/blast/{pdb_id}_{chain_id}.{ref_version}.domain_blast.xml",
        f"{batch_path}/blast/domain/batch_*/{pdb_id}_{chain_id}*.xml"
    ]
    
    # Find files matching patterns
    chain_blast_files = []
    for pattern in chain_blast_patterns:
        matches = glob.glob(pattern)
        chain_blast_files.extend(matches)
    
    domain_blast_files = []
    for pattern in domain_blast_patterns:
        matches = glob.glob(pattern)
        domain_blast_files.extend(matches)
    
    return {
        'chain_blast_result': chain_blast_files,
        'domain_blast_result': domain_blast_files
    }

def register_file_in_db(db, process_id: int, file_type: str, 
                      file_path: str, batch_path: str, dry_run: bool = False) -> bool:
    """Register a file in the database"""
    logger = logging.getLogger("register_files")
    
    # Verify file exists
    if not os.path.exists(file_path):
        logger.warning(f"File does not exist: {file_path}")
        return False
    
    # Get file size
    file_size = os.path.getsize(file_path)
    
    # Convert to relative path for database storage
    rel_path = os.path.relpath(file_path, batch_path)
    
    # Check if there's an existing registration for this file type
    query = """
    SELECT id, file_path, file_exists FROM ecod_schema.process_file
    WHERE process_id = %s AND file_type = %s
    """
    
    existing = db.execute_dict_query(query, (process_id, file_type))
    
    if dry_run:
        if existing:
            logger.info(f"[DRY RUN] Would update {file_type} registration for process {process_id}")
        else:
            logger.info(f"[DRY RUN] Would register {file_type} for process {process_id}: {rel_path}")
        return True
    
    try:
        if existing:
            # Update existing registration
            db.update(
                "ecod_schema.process_file",
                {
                    "file_path": rel_path,
                    "file_exists": True,
                    "file_size": file_size,
                    "last_checked": "NOW()"
                },
                "id = %s",
                (existing[0]['id'],)
            )
            logger.debug(f"Updated {file_type} registration for process {process_id}")
        else:
            # Create new registration
            db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": file_type,
                    "file_path": rel_path,
                    "file_exists": True,
                    "file_size": file_size,
                    "last_checked": "NOW()"
                }
            )
            logger.debug(f"Registered {file_type} for process {process_id}")
        
        return True
    except Exception as e:
        logger.error(f"Error registering file: {str(e)}")
        return False

def scan_and_register_files(db, batch_id: int, file_types: List[str], 
                          dry_run: bool = False) -> Dict[str, int]:
    """Scan batch directory for files and register them in the database"""
    logger = logging.getLogger("register_files")
    
    # Get batch info
    batch_info = get_batch_info(db, batch_id)
    if not batch_info:
        return {}
    
    batch_path = batch_info['base_path']
    ref_version = batch_info['ref_version']
    batch_name = batch_info['batch_name']
    
    logger.info(f"Scanning files for batch {batch_id} ({batch_name})")
    logger.info(f"Batch path: {batch_path}")
    logger.info(f"Reference version: {ref_version}")
    
    # Get proteins in batch
    proteins = get_batch_proteins(db, batch_id)
    
    # Get existing file registrations
    process_ids = [p['process_id'] for p in proteins]
    existing_files = get_existing_file_registrations(db, process_ids)
    
    # Initialize statistics
    stats = {
        'batch_id': batch_id,
        'batch_name': batch_name,
        'total_proteins': len(proteins),
        'registered': {ft: 0 for ft in file_types},
        'already_registered': {ft: 0 for ft in file_types},
        'not_found': {ft: 0 for ft in file_types}
    }
    
    # Process each protein
    for protein in proteins:
        process_id = protein['process_id']
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        
        # Use path_utils to get standardized paths
        paths = get_standardized_paths(
            batch_path=batch_path,
            pdb_id=pdb_id,
            chain_id=chain_id,
            ref_version=ref_version,
            create_dirs=False
        )
        
        # Override legacy path finding with explicit BLAST file search
        blast_files = find_blast_files(batch_path, pdb_id, chain_id, ref_version)
        
        # Process each requested file type
        for file_type in file_types:
            # Skip if already registered
            if process_id in existing_files and file_type in existing_files[process_id]:
                file_info = existing_files[process_id][file_type]
                if file_info['file_exists']:
                    stats['already_registered'][file_type] += 1
                    continue
            
            # Determine file path based on type
            file_path = None
            
            if file_type == 'chain_blast_result' and blast_files['chain_blast_result']:
                # Use the first found chain blast file
                file_path = blast_files['chain_blast_result'][0]
            elif file_type == 'domain_blast_result' and blast_files['domain_blast_result']:
                # Use the first found domain blast file
                file_path = blast_files['domain_blast_result'][0]
            elif file_type in paths:
                # Use standardized path for other file types
                file_path = paths[file_type]
                if not os.path.exists(file_path):
                    # Try alternate path patterns
                    legacy_files = find_files_with_legacy_paths(
                        batch_path, pdb_id, chain_id, ref_version
                    )
                    if file_type in legacy_files and legacy_files[file_type]['exists_at']:
                        file_path = legacy_files[file_type]['exists_at']
            
            # Register file if found
            if file_path and os.path.exists(file_path):
                if register_file_in_db(db, process_id, file_type, file_path, batch_path, dry_run):
                    stats['registered'][file_type] += 1
            else:
                stats['not_found'][file_type] += 1
                
                if file_type.endswith('blast_result'):
                    logger.debug(f"No {file_type} found for {pdb_id}_{chain_id}")
    
    return stats

def get_missing_file_counts(db, batch_id: int, file_types: List[str]) -> Dict[str, int]:
    """Get counts of missing file registrations for a batch"""
    logger = logging.getLogger("register_files")
    
    missing_counts = {}
    
    # Get total protein count in batch
    count_query = """
    SELECT COUNT(*) as count 
    FROM ecod_schema.process_status 
    WHERE batch_id = %s
    """
    
    result = db.execute_dict_query(count_query, (batch_id,))
    total = result[0]['count'] if result else 0
    
    for file_type in file_types:
        # Count registrations of this file type
        reg_query = """
        SELECT COUNT(*) as count 
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        WHERE ps.batch_id = %s AND pf.file_type = %s
        """
        
        result = db.execute_dict_query(reg_query, (batch_id, file_type))
        registered = result[0]['count'] if result else 0
        
        missing_counts[file_type] = total - registered
    
    return missing_counts

def main():
    parser = argparse.ArgumentParser(description='Register batch files in the database')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--blast', action='store_true',
                      help='Register BLAST result files')
    parser.add_argument('--hhsearch', action='store_true',
                      help='Register HHSearch result files')
    parser.add_argument('--all', action='store_true',
                      help='Register all file types')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show what would be done without making changes')
    parser.add_argument('--log-file', type=str,
                      help='Path to log file')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()

    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("register_files")

    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)

    # Determine which file types to process
    file_types = []
    
    if args.blast or args.all:
        file_types.extend(['chain_blast_result', 'domain_blast_result'])
        
    if args.hhsearch or args.all:
        file_types.extend(['hhr', 'hh_xml'])
        
    if not file_types:
        logger.error("No file types specified. Use --blast, --hhsearch, or --all")
        return 1
    
    # Get batch info
    batch_info = get_batch_info(db, args.batch_id)
    if not batch_info:
        return 1
    
    # Get counts of missing registrations
    missing_counts = get_missing_file_counts(db, args.batch_id, file_types)
    
    logger.info(f"Missing file registrations for batch {args.batch_id}:")
    for file_type, count in missing_counts.items():
        logger.info(f"  {file_type}: {count}")
    
    # Skip if no missing registrations
    if all(count == 0 for count in missing_counts.values()):
        logger.info("No missing registrations found. Nothing to do.")
        return 0
    
    # Scan and register files
    stats = scan_and_register_files(db, args.batch_id, file_types, args.dry_run)
    
    # Print summary
    logger.info("=" * 60)
    logger.info(f"Summary for batch {args.batch_id} ({stats['batch_name']}):")
    logger.info(f"Total proteins: {stats['total_proteins']}")
    
    for file_type in file_types:
        registered = stats['registered'][file_type]
        already = stats['already_registered'][file_type]
        not_found = stats['not_found'][file_type]
        
        logger.info(f"\n{file_type}:")
        if args.dry_run:
            logger.info(f"  Would register: {registered}")
        else:
            logger.info(f"  Registered: {registered}")
        logger.info(f"  Already registered: {already}")
        logger.info(f"  Not found: {not_found}")
    
    if args.dry_run:
        logger.info("\nThis was a dry run. No changes were made to the database.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
