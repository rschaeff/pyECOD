#!/usr/bin/env python3
"""
register_alt_rep_profiles.py - Find, copy and register alternative representative profile files in the database
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

def get_batch_ids(db, batch_id: Optional[int] = None, all_batches: bool = False) -> List[int]:
    """Get list of batch IDs to process"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    if batch_id is not None and not all_batches:
        logger.info(f"Processing single batch: {batch_id}")
        return [batch_id]
    
    # Query to get all batches with alt representatives
    query = """
    SELECT DISTINCT ps.batch_id
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    JOIN pdb_analysis.alt_representative_proteins a ON p.pdb_id = a.pdb_id AND p.chain_id = a.chain_id
    JOIN ecod_schema.batch b ON ps.batch_id = b.id
    WHERE ps.is_representative = TRUE
    ORDER BY ps.batch_id
    """
    
    result = db.execute_query(query)
    batch_ids = [row[0] for row in result]
    
    # If a starting batch_id was provided, start from there
    if batch_id is not None:
        batch_ids = [b for b in batch_ids if b >= batch_id]
    
    logger.info(f"Found {len(batch_ids)} batches with alternative representatives")
    
    return batch_ids

def get_batch_alt_reps(db, batch_id: int) -> List[Dict[str, Any]]:
    """Get list of alternative representative proteins from a specific batch"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
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
        ps.is_representative = FALSE
    """
    
    result = db.execute_dict_query(check_query, (batch_id,))
    if result and result[0]['count'] > 0:
        logger.warning(f"Found {result[0]['count']} alt reps with is_representative = FALSE in batch {batch_id}. Consider fixing this!")
    
    # Main query to get alt reps, enforcing is_representative = TRUE for safety
    # Now includes both novel and existing classifications
    query = """
    SELECT 
        ps.id as process_id,
        p.id as protein_id,
        p.pdb_id,
        p.chain_id,
        p.source_id,
        b.base_path,
        rc.classification
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
        ps.is_representative = TRUE
    ORDER BY 
        ps.id
    """
    
    proteins = db.execute_dict_query(query, (batch_id,))
    
    # Count novel and existing
    novel_count = sum(1 for p in proteins if p['classification'] == 'novel')
    existing_count = sum(1 for p in proteins if p['classification'] == 'existing')
    
    logger.info(f"Found {len(proteins)} alternative representatives with is_representative = TRUE in batch {batch_id}")
    logger.info(f"   Novel: {novel_count}, Existing: {existing_count}")
    
    return proteins

def fix_representative_flags(db, batch_id: int, dry_run: bool = False) -> int:
    """Fix is_representative flags for alternative representatives"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    # Modified to include all alt reps (both novel and existing)
    query = """
    UPDATE ecod_schema.process_status ps
    SET is_representative = TRUE
    FROM ecod_schema.protein p
    JOIN pdb_analysis.alt_representative_proteins a ON p.pdb_id = a.pdb_id AND p.chain_id = a.chain_id
    WHERE 
        ps.protein_id = p.id AND
        ps.batch_id = %s AND
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
        WHERE 
            ps.batch_id = %s AND
            ps.is_representative = FALSE
        """
        result = db.execute_dict_query(count_query, (batch_id,))
        count = result[0]['count'] if result else 0
        logger.info(f"[DRY RUN] Would update is_representative flag for {count} alt reps in batch {batch_id}")
        return count
    else:
        result = db.execute_query(query, (batch_id,))
        count = len(result) if result else 0
        logger.info(f"Updated is_representative flag for {count} alt reps in batch {batch_id}")
        return count

def find_profile_file(source_dirs: List[str], protein_id: str, file_ext: str) -> Optional[str]:
    """Search for profile file (HMM or A3M) in source directories"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    # Try looking for an exact match
    for source_dir in source_dirs:
        # Check for profile files with the specified extension
        path = os.path.join(source_dir, f"{protein_id}{file_ext}")
        if os.path.exists(path):
            logger.debug(f"Found exact match for {protein_id}{file_ext}: {path}")
            return path
    
    # If no exact match, try more flexible search through subdirectories
    pdb_id, chain_id = protein_id.split('_', 1)
    possible_matches = []
    
    for source_dir in source_dirs:
        # Check numbered subdirectories if they exist
        try:
            subdirs = [d for d in os.listdir(source_dir) 
                    if os.path.isdir(os.path.join(source_dir, d)) and d.isdigit()]
        except (FileNotFoundError, PermissionError):
            logger.warning(f"Cannot access directory: {source_dir}")
            continue
        
        # Add the base directory and all numbered subdirectories to search paths
        search_dirs = [source_dir] + [os.path.join(source_dir, d) for d in subdirs]
        
        for search_dir in search_dirs:
            # Look for files with matching PDB ID and chain ID
            pattern = os.path.join(search_dir, f"{pdb_id}*{chain_id}*{file_ext}")
            matches = glob.glob(pattern)
            
            if matches:
                # Sort by length of filename (shorter is likely more exact)
                matches.sort(key=lambda x: len(os.path.basename(x)))
                logger.debug(f"Found possible match for {protein_id}{file_ext}: {matches[0]}")
                possible_matches.extend(matches)
    
    if possible_matches:
        logger.info(f"Found best match for {protein_id}{file_ext}: {possible_matches[0]}")
        return possible_matches[0]
    
    logger.warning(f"No {file_ext} file found for {protein_id}")
    return None

def copy_profile_to_batch(source_path: str, batch_path: str, protein_id: str, file_type: str, dry_run: bool = False) -> Optional[str]:
    """Copy profile file to batch directory"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    # All profile files (both HMM/HHM and A3M) go to the hhsearch directory
    dest_dir = os.path.join(batch_path, "hhsearch")
    
    # Determine file extension based on file type
    if file_type in ['hhm', 'hmm']:
        dest_ext = '.hhm'  # Always use .hhm extension for hmm/hhm files
    elif file_type == 'a3m':
        dest_ext = '.a3m'
    else:
        logger.error(f"Unknown file type: {file_type}")
        return None
    
    # Ensure directory exists
    if not os.path.exists(dest_dir) and not dry_run:
        os.makedirs(dest_dir, exist_ok=True)
        logger.info(f"Created directory: {dest_dir}")
    
    # Determine destination path
    dest_path = os.path.join(dest_dir, f"{protein_id}{dest_ext}")
    
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

def register_profile_in_db(db, process_id: int, file_path: str, file_type: str, dry_run: bool = False) -> bool:
    """Register profile file in the process_file table"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    try:
        # First check if this file is already registered
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = %s
        """
        
        existing = db.execute_query(check_query, (process_id, file_type))
        
        if dry_run:
            if existing:
                logger.info(f"[DRY RUN] Would update existing {file_type} record for process {process_id}")
            else:
                logger.info(f"[DRY RUN] Would create new {file_type} record for process {process_id}")
            return True
        
        file_exists = os.path.exists(file_path)
        file_size = os.path.getsize(file_path) if file_exists else 0
        
        if existing:
            # Update existing record
            db.update(
                "ecod_schema.process_file",
                {
                    "file_path": file_path,
                    "file_exists": file_exists,
                    "file_size": file_size,
                    "last_checked": "NOW()"
                },
                "id = %s",
                (existing[0][0],)
            )
            logger.debug(f"Updated existing {file_type} record for process {process_id}")
        else:
            # Insert new record
            db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": file_type,
                    "file_path": file_path,
                    "file_exists": file_exists,
                    "file_size": file_size,
                    "last_checked": "NOW()"
                }
            )
            logger.debug(f"Created new {file_type} record for process {process_id}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error registering {file_type} in database: {str(e)}")
        return False

def update_process_status(db, process_id: int, success: bool, file_type: str, error_message: Optional[str] = None, dry_run: bool = False) -> None:
    """Update process status after file registration"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    try:
        stage = f"{file_type}_registered"
        status = "completed" if success else "error"
        
        if dry_run:
            logger.info(f"[DRY RUN] Would update process_status {process_id} to {stage}/{status}")
            return
        
        db.update(
            "ecod_schema.process_status",
            {
                "current_stage": stage,
                "status": status,
                "error_message": error_message,
                "updated_at": "NOW()"
            },
            "id = %s",
            (process_id,)
        )
    except Exception as e:
        logger.error(f"Error updating process status: {str(e)}")

def process_profiles_for_protein(db, protein: Dict[str, Any], source_dirs: List[str], file_types: List[str], dry_run: bool = False) -> Dict[str, bool]:
    """Process profile files for a single protein"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    protein_id = protein['protein_id']
    process_id = protein['process_id']
    source_id = protein.get('source_id') or f"{protein['pdb_id']}_{protein['chain_id']}"
    base_path = protein['base_path']
    classification = protein.get('classification', 'unknown')
    
    results = {}
    
    for file_type in file_types:
        logger.debug(f"Processing {file_type} for alt rep {source_id} (protein_id={protein_id}, process_id={process_id}, classification={classification})")
        
        # Determine file extension
        file_ext = '.hhm' if file_type in ['hhm', 'hmm'] else f'.{file_type}'
        
        # Find profile file in source directories
        source_path = find_profile_file(source_dirs, source_id, file_ext)
        
        if not source_path:
            logger.error(f"{file_type.upper()} file not found for {source_id}")
            results[file_type] = False
            if not dry_run:
                update_process_status(db, process_id, False, file_type, f"{file_type.upper()} file not found in source directories", dry_run)
            continue
        
        # Copy profile file to batch directory
        dest_path = copy_profile_to_batch(source_path, base_path, source_id, file_type, dry_run)
        
        if not dest_path:
            results[file_type] = False
            if not dry_run:
                update_process_status(db, process_id, False, file_type, f"Failed to copy {file_type.upper()} file to batch directory", dry_run)
            continue
        
        # Register profile in database
        db_success = register_profile_in_db(db, process_id, dest_path, file_type, dry_run)
        
        if db_success:
            results[file_type] = True
            logger.info(f"Successfully registered {file_type.upper()} for alt rep {source_id}")
            if not dry_run:
                update_process_status(db, process_id, True, file_type, None, dry_run)
        else:
            results[file_type] = False
            logger.error(f"Database registration failed for {file_type.upper()} alt rep {source_id}")
            if not dry_run:
                update_process_status(db, process_id, False, file_type, f"Failed to register {file_type.upper()} in database", dry_run)
    
    return results

def process_batch_profiles(db, batch_id: int, source_dirs: List[str], file_types: List[str], dry_run: bool = False) -> Dict[str, Any]:
    """Process profile files for alternative representatives in a batch"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    proteins = get_batch_alt_reps(db, batch_id)
    
    if not proteins:
        logger.warning(f"No alternative representatives found in batch {batch_id}")
        return {
            "batch_id": batch_id,
            "total": 0,
            "file_types": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types},
            "classifications": {
                "novel": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types},
                "existing": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types}
            }
        }
    
    stats = {
        "batch_id": batch_id,
        "total": len(proteins),
        "file_types": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types},
        "classifications": {
            "novel": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types},
            "existing": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types}
        }
    }
    
    for protein in proteins:
        classification = protein.get('classification', 'unknown')
        results = process_profiles_for_protein(db, protein, source_dirs, file_types, dry_run)
        
        for file_type, success in results.items():
            if success:
                stats["file_types"][file_type]["success"] += 1
                if classification in stats["classifications"]:
                    stats["classifications"][classification][file_type]["success"] += 1
            else:
                stats["file_types"][file_type]["failed"] += 1
                if classification in stats["classifications"]:
                    stats["classifications"][classification][file_type]["failed"] += 1
    
    return stats

def print_stats(stats: Dict[str, Any], file_types: List[str], dry_run: bool = False):
    """Print statistics in a nice format"""
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    prefix = "[DRY RUN] Would process" if dry_run else "Processed"
    
    logger.info(f"{prefix} batch {stats['batch_id']}: Total proteins: {stats['total']}")
    
    for file_type in file_types:
        ft_stats = stats["file_types"][file_type]
        logger.info(f"  {file_type.upper()}: Success: {ft_stats['success']}, Failed: {ft_stats['failed']}")
    
    logger.info("By classification:")
    for classification in ["novel", "existing"]:
        logger.info(f"  {classification.capitalize()}:")
        for file_type in file_types:
            ft_stats = stats["classifications"][classification][file_type]
            logger.info(f"    {file_type.upper()}: Success: {ft_stats['success']}, Failed: {ft_stats['failed']}")

def main():
    parser = argparse.ArgumentParser(description='Find, copy and register alternative representative profile files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Batch ID to process (if not specified, all batches will be processed)')
    parser.add_argument('--all-batches', action='store_true',
                      help='Process all batches with alt representatives')
    parser.add_argument('--source-dir', type=str, action='append', required=True,
                      help='Source directory to search for profile files (can be used multiple times)')
    parser.add_argument('--no-hmm', action='store_true',
                      help='Skip HMM/HHM files')
    parser.add_argument('--no-a3m', action='store_true',
                      help='Skip A3M files')
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
    logger = logging.getLogger("ecod.alt_rep_profiles")
    
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Determine which file types to process
    file_types = []
    if not args.no_hmm:
        file_types.append('hhm')
    if not args.no_a3m:
        file_types.append('a3m')
    
    if not file_types:
        logger.error("No file types selected for processing (both --no-hmm and --no-a3m specified)")
        return 1
    
    # Get list of batches to process
    batch_ids = get_batch_ids(db, args.batch_id, args.all_batches)
    
    if not batch_ids:
        logger.error("No batches found to process")
        return 1
    
    # Process each batch
    overall_stats = {
        "total_batches": len(batch_ids),
        "total_proteins": 0,
        "file_types": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types},
        "classifications": {
            "novel": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types},
            "existing": {ft: {"success": 0, "failed": 0, "not_found": 0} for ft in file_types}
        }
    }
    
    for batch_id in batch_ids:
        logger.info(f"Processing batch {batch_id}")
        
        # Option to fix representative flags
        if args.fix_flags:
            fixed = fix_representative_flags(db, batch_id, args.dry_run)
            if not args.dry_run and fixed > 0:
                logger.info(f"Fixed {fixed} alternative representatives with is_representative = FALSE")
        
        # Process profiles for the batch
        stats = process_batch_profiles(db, batch_id, args.source_dir, file_types, args.dry_run)
        
        # Print stats for this batch
        print_stats(stats, file_types, args.dry_run)
        
        # Update overall stats
        overall_stats["total_proteins"] += stats["total"]
        for file_type in file_types:
            for stat_key in ["success", "failed", "not_found"]:
                overall_stats["file_types"][file_type][stat_key] += stats["file_types"][file_type][stat_key]
                for classification in ["novel", "existing"]:
                    overall_stats["classifications"][classification][file_type][stat_key] += stats["classifications"][classification][file_type][stat_key]
    
    # Print overall stats
    logger.info("=" * 50)
    logger.info(f"Overall statistics for {overall_stats['total_batches']} batches:")
    logger.info(f"Total proteins: {overall_stats['total_proteins']}")
    
    for file_type in file_types:
        ft_stats = overall_stats["file_types"][file_type]
        logger.info(f"  {file_type.upper()}: Success: {ft_stats['success']}, Failed: {ft_stats['failed']}")
    
    logger.info("By classification:")
    for classification in ["novel", "existing"]:
        logger.info(f"  {classification.capitalize()}:")
        for file_type in file_types:
            ft_stats = overall_stats["classifications"][classification][file_type]
            logger.info(f"    {file_type.upper()}: Success: {ft_stats['success']}, Failed: {ft_stats['failed']}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())