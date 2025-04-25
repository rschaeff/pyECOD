#!/usr/bin/env python3
"""
register_alt_rep_hmms.py - Register alternative representative HHMs in the database
"""

import os
import sys
import logging
import argparse
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

def get_hhm_path(base_path: str, protein_data: Dict[str, Any]) -> str:
    """Get the HHM file path for a protein in a batch"""
    # Use source_id if available, otherwise construct from pdb_id and chain_id
    if protein_data.get('source_id'):
        identifier = protein_data['source_id']
    else:
        identifier = f"{protein_data['pdb_id']}_{protein_data['chain_id']}"
    
    # HHM should be in the hhsearch directory of the batch
    hhm_path = os.path.join(base_path, "hhsearch", f"{identifier}.hhm")
    
    return hhm_path

def verify_hhm_exists(hhm_path: str) -> Tuple[bool, Optional[str]]:
    """Verify that the HHM file exists and is valid"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    try:
        if not os.path.exists(hhm_path):
            return False, f"HHM file not found: {hhm_path}"
        
        # Check if file is non-empty
        if os.path.getsize(hhm_path) == 0:
            return False, f"HHM file is empty: {hhm_path}"
        
        # You could add more validation here if needed
        
        return True, None
    
    except Exception as e:
        return False, f"Error verifying HHM file: {str(e)}"

def register_hhm_in_db(db, protein_id: int, hhm_path: str, process_id: int) -> bool:
    """Register HHM file in the database"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    try:
        # First check if this HHM is already registered
        check_query = """
        SELECT id FROM ecod_schema.protein_profile
        WHERE protein_id = %s AND profile_type = 'hhm'
        """
        
        existing = db.execute_query(check_query, (protein_id,))
        
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

def update_process_status(db, process_id: int, success: bool, error_message: Optional[str] = None) -> None:
    """Update process status after HHM registration"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    try:
        status = "completed" if success else "error"
        
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

def process_batch_hmms(db, batch_id: int, dry_run: bool = False) -> Dict[str, int]:
    """Process HHMs for alternative representatives in a batch"""
    logger = logging.getLogger("ecod.alt_rep_hmms")
    
    proteins = get_batch_alt_reps(db, batch_id)
    
    if not proteins:
        logger.warning(f"No alternative representatives found in batch {batch_id}")
        return {"total": 0, "success": 0, "failed": 0}
    
    stats = {"total": len(proteins), "success": 0, "failed": 0}
    
    for protein in proteins:
        protein_id = protein['protein_id']
        process_id = protein['process_id']
        source_id = protein.get('source_id') or f"{protein['pdb_id']}_{protein['chain_id']}"
        base_path = protein['base_path']
        
        # Get HHM path
        hhm_path = get_hhm_path(base_path, protein)
        
        logger.debug(f"Processing HHM for alt rep {source_id} (protein_id={protein_id})")
        logger.debug(f"  HHM path: {hhm_path}")
        
        if dry_run:
            logger.info(f"[DRY RUN] Would register {hhm_path} in database")
            continue
        
        # Verify HHM exists
        hhm_exists, error_message = verify_hhm_exists(hhm_path)
        
        if hhm_exists:
            # Register in database
            db_success = register_hhm_in_db(db, protein_id, hhm_path, process_id)
            
            if db_success:
                stats["success"] += 1
                logger.info(f"Successfully registered HHM for alt rep {source_id}")
                update_process_status(db, process_id, True)
            else:
                stats["failed"] += 1
                error_message = "Failed to register HHM in database"
                logger.error(f"Database registration failed for alt rep {source_id}")
                update_process_status(db, process_id, False, error_message)
        else:
            stats["failed"] += 1
            logger.error(f"Failed to find valid HHM for alt rep {source_id}: {error_message}")
            update_process_status(db, process_id, False, error_message)
    
    return stats

def main():
    parser = argparse.ArgumentParser(description='Register alternative representative HHMs')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
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
    
    stats = process_batch_hmms(db, args.batch_id, args.dry_run)
    
    if args.dry_run:
        logger.info(f"Dry run completed for batch {args.batch_id}")
        logger.info(f"Would process {stats['total']} HHMs")
    else:
        logger.info(f"Completed processing HHMs for batch {args.batch_id}")
        logger.info(f"Total: {stats['total']}, Success: {stats['success']}, Failed: {stats['failed']}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())