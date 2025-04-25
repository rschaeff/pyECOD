#!/usr/bin/env python3
"""
register_alt_rep_hmms.py - Find, copy and register alternative representative HHMs in the database
"""

import os
import sys
import logging
import argparse
import shutil
import glob
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager

def setup_logging(verbose: bool = False, log_file: str = None):
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

def get_batch_alt_reps(db, batch_id: int) -> List[Dict[str, Any]]:
    """Get list of alternative representative proteins from a specific batch"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    # First check if there are any alt reps with is_representative = FALSE
    check_query = """
    SELECT 
        COUNT(*) as count
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    JOIN 
        pdb_analysis.alt_representative_proteins a ON p.pdb_id = a.pdb_id AND p.chain_id = a.chain_id
    JOIN 
        pdb_analysis.rep_classification rc ON a.id = rc.alt_rep_id
    WHERE 
        ps.batch_id = %s AND
        rc.classification = 'novel' AND
        ps.is_representative = FALSE
    """
    
    result = db.execute_dict_query(check_query, (batch_id,))
    if result and result[0]['count'] > 0:
        logger.warning(f"Found {result[0]['count']} alt reps with is_representative = FALSE. Consider fixing this!")
    
    # Main query to get alt reps, enforcing is_representative = TRUE for safety
    query = """
    SELECT 
        ps.id as process_id,
        p.id as protein_id,
        p.pdb_id,
        p.chain_id,
        p.source_id,
        b.base_path
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    JOIN 
        pdb_analysis.alt_representative_proteins a ON p.pdb_id = a.pdb_id AND p.chain_id = a.chain_id
    JOIN 
        pdb_analysis.rep_classification rc ON a.id = rc.alt_rep_id
    WHERE 
        ps.batch_id = %s AND
        rc.classification = 'novel' AND
        ps.is_representative = TRUE
    ORDER BY 
        ps.id
    """
    
    proteins = db.execute_dict_query(query, (batch_id,))
    logger.info(f"Found {len(proteins)} novel alternative representatives with is_representative = TRUE in batch {batch_id}")
    
    return proteins

def fix_representative_flags(db, batch_id: int, dry_run: bool = False) -> int:
    """Fix is_representative flags for alternative representatives"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    query = """
    UPDATE ecod_schema.process_status ps
    SET is_representative = TRUE
    FROM ecod_schema.protein p
    JOIN pdb_analysis.alt_representative_proteins a ON p.pdb_id = a.pdb_id AND p.chain_id = a.chain_id
    JOIN pdb_analysis.rep_classification rc ON a.id = rc.alt_rep_id
    WHERE 
        ps.protein_id = p.id AND
        ps.batch_id = %s AND
        rc.classification = 'novel' AND
        ps.is_representative = FALSE
    RETURNING ps.id
    """
    
    if dry_run:
        # In dry run mode, just count how many would be updated
        count_query = """
        SELECT COUNT(*) as count
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        JOIN pdb_analysis.alt_representative_proteins a ON p.pdb_id = a.pdb_id AND p.chain_id = a.chain_id
        JOIN pdb_analysis.rep_classification rc ON a.id = rc.alt_rep_id
        WHERE 
            ps.batch_id = %s AND
            rc.classification = 'novel' AND
            ps.is_representative = FALSE
        """
        result = db.execute_dict_query(count_query, (batch_id,))
        count = result[0]['count'] if result else 0
        logger.info(f"[DRY RUN] Would update is_representative flag for {count} alt reps")
        return count
    else:
        result = db.execute_query(query, (batch_id,))
        count = len(result) if result else 0
        logger.info(f"Updated is_representative flag for {count} alt reps")
        return count

def find_hmm_file(source_dirs: List[str], protein_id: str) -> Optional[str]:
    """Search for HMM file in source directories"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    # Try looking for an exact match
    for source_dir in source_dirs:
        # Check for HMM files (both .hmm and .hhm extensions)
        for ext in ['.hmm', '.hhm']:
            path = os.path.join(source_dir, f"{protein_id}{ext}")
            if os.path.exists(path):
                logger.debug(f"Found exact match for {protein_id}: {path}")
                return path
    
    # If no exact match, try more flexible search through subdirectories
    pdb_id, chain_id = protein_id.split('_', 1)
    possible_matches = []
    
    for source_dir in source_dirs:
        # Check numbered subdirectories if they exist
        subdirs = [d for d in os.listdir(source_dir) 
                  if os.path.isdir(os.path.join(source_dir, d)) and d.isdigit()]
        
        # Add the base directory and all numbered subdirectories to search paths
        search_dirs = [source_dir] + [os.path.join(source_dir, d) for d in subdirs]
        
        for search_dir in search_dirs:
            # Look for files with matching PDB ID and chain ID
            for ext in ['.hmm', '.hhm']:
                pattern = os.path.join(search_dir, f"{pdb_id}*{chain_id}*{ext}")
                matches = glob.glob(pattern)
                
                if matches:
                    # Sort by length of filename (shorter is likely more exact)
                    matches.sort(key=lambda x: len(os.path.basename(x)))
                    logger.debug(f"Found possible match for {protein_id}: {matches[0]}")
                    possible_matches.extend(matches)
    
    if possible_matches:
        logger.info(f"Found best match for {protein_id}: {possible_matches[0]}")
        return possible_matches[0]
    
    logger.warning(f"No HMM file found for {protein_id}")
    return None

def copy_hmm_to_batch(source_path: str, batch_path: str, protein_id: str, dry_run: bool = False) -> Optional[str]:
    """Copy HMM file to batch directory, converting to .hhm if needed"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    # Ensure hhsearch directory exists
    dest_dir = os.path.join(batch_path, "hhsearch")
    if not os.path.exists(dest_dir) and not dry_run:
        os.makedirs(dest_dir, exist_ok=True)
        logger.info(f"Created directory: {dest_dir}")
    
    # Determine destination path (always use .hhm extension)
    dest_path = os.path.join(dest_dir, f"{protein_id}.hhm")
    
    if dry_run:
        logger.info(f"[DRY RUN] Would copy {source_path} to {dest_path}")
        return dest_path
    
    try:
        shutil.copy2(source_path, dest_path)
        logger.info(f"Copied {source_path} to {dest_path}")
        return dest_path
    except Exception as e:
        logger.error(f"Failed to copy {source_path} to {dest_path}: {str(e)}")
        return None

def register_hhm_in_db(db, protein_id: int, hhm_path: str, process_id: int, dry_run: bool = False) -> bool:
    """Register HHM file in the database"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    try:
        # First check if this HHM is already registered
        check_query = """
        SELECT id FROM ecod_schema.protein_profile
        WHERE protein_id = %s AND profile_type = 'hhm'
        """
        
        existing = db.execute_query(check_query, (protein_id,))
        
        if dry_run:
            if existing:
                logger.info(f"[DRY RUN] Would update existing HHM record for protein {protein_id}")
            else:
                logger.info(f"[DRY RUN] Would create new HHM record for protein {protein_id}")
            return True
        
        if existing:
            # Update existing record
            db.update(
                "ecod_schema.protein_profile",
                {
                    "file_path": hhm_path,
                    "updated_at": "NOW()"
                },
                "id = %s",
                (existing[0][0],)
            )
            logger.debug(f"Updated existing HHM record for protein {protein_id}")
        else:
            # Insert new record
            db.insert(
                "ecod_schema.protein_profile",
                {
                    "protein_id": protein_id,
                    "profile_type": "hhm",
                    "file_path": hhm_path,
                    "process_id": process_id
                }
            )
            logger.debug(f"Created new HHM record for protein {protein_id}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error registering HHM in database: {str(e)}")
        return False

def update_process_status(db, process_id: int, success: bool, error_message: Optional[str] = None, dry_run: bool = False) -> None:
    """Update process status after HHM registration"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    try:
        status = "completed" if success else "error"
        
        if dry_run:
            logger.info(f"[DRY RUN] Would update process_status {process_id} to {status}")
            return
        
        db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "hhm_registered",
                "status": status,
                "error_message": error_message
            },
            "id = %s",
            (process_id,)
        )
    except Exception as e:
        logger.error(f"Error updating process status: {str(e)}")

def process_batch_hmms(db, batch_id: int, source_dirs: List[str], dry_run: bool = False) -> Dict[str, int]:
    """Process HMMs for alternative representatives in a batch"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    proteins = get_batch_alt_reps(db, batch_id)
    
    if not proteins:
        logger.warning(f"No alternative representatives found in batch {batch_id}")
        return {"total": 0, "success": 0, "failed": 0, "not_found": 0}
    
    stats = {"total": len(proteins), "success": 0, "failed": 0, "not_found": 0}
    
    for protein in proteins:
        protein_id = protein['protein_id']
        process_id = protein['process_id']
        source_id = protein.get('source_id') or f"{protein['pdb_id']}_{protein['chain_id']}"
        base_path = protein['base_path']
        
        logger.debug(f"Processing HMM for alt rep {source_id} (protein_id={protein_id})")
        
        # Find HMM file in source directories
        source_path = find_hmm_file(source_dirs, source_id)
        
        if not source_path:
            logger.error(f"HMM file not found for {source_id}")
            stats["not_found"] += 1
            if not dry_run:
                update_process_status(db, process_id, False, f"HMM file not found in source directories")
            continue
        
        # Copy HMM file to batch directory
        dest_path = copy_hmm_to_batch(source_path, base_path, source_id, dry_run)
        
        if not dest_path:
            stats["failed"] += 1
            if not dry_run:
                update_process_status(db, process_id, False, f"Failed to copy HMM file to batch directory")
            continue
        
        # Register HHM in database
        db_success = register_hhm_in_db(db, protein_id, dest_path, process_id, dry_run)
        
        if db_success:
            stats["success"] += 1
            logger.info(f"Successfully registered HHM for alt rep {source_id}")
            if not dry_run:
                update_process_status(db, process_id, True)
        else:
            stats["failed"] += 1
            logger.error(f"Database registration failed for alt rep {source_id}")
            if not dry_run:
                update_process_status(db, process_id, False, "Failed to register HHM in database")
    
    return stats

def main():
    parser = argparse.ArgumentParser(description='Find, copy and register alternative representative HHMs')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--source-dir', type=str, action='append', required=True,
                      help='Source directory to search for HMM files (can be used multiple times)')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show what would be done without making changes')
    parser.add_argument('--fix-flags', action='store_true',
                      help='Fix is_representative flags for alt reps')
    parser.add_argument('--log-file', type=str,
                      help='Path to log file')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Option to fix representative flags
    if args.fix_flags:
        fixed = fix_representative_flags(db, args.batch_id, args.dry_run)
        if not args.dry_run and fixed > 0:
            logger.info(f"Fixed {fixed} alternative representatives with is_representative = FALSE")
    
    if args.dry_run:
        logger.info(f"Dry run: checking batch {args.batch_id}")
    
    stats = process_batch_hmms(db, args.batch_id, args.source_dir, args.dry_run)
    
    if args.dry_run:
        logger.info(f"Dry run completed for batch {args.batch_id}")
        logger.info(f"Would process {stats['total']} HHMs")
        logger.info(f"Found: {stats['success']}, Not found: {stats['not_found']}, Failed: {stats['failed']}")
    else:
        logger.info(f"Completed processing HHMs for batch {args.batch_id}")
        logger.info(f"Total: {stats['total']}, Success: {stats['success']}, Not found: {stats['not_found']}, Failed: {stats['failed']}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())