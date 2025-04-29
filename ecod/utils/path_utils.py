# ecod/utils/path_utils.py
"""
Standardized file path utilities for pyECOD

This module provides consistent path construction for all file types
used in the pyECOD processing pipeline, ensuring standardization across
different components.
"""

import os
import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

def get_standardized_paths(batch_path: str, pdb_id: str, chain_id: str, 
                          ref_version: str, create_dirs: bool = True) -> Dict[str, str]:
    """Get standardized file paths for a protein chain
    
    Args:
        batch_path: Path to batch directory
        pdb_id: PDB ID
        chain_id: Chain ID
        ref_version: Reference version
        create_dirs: Whether to create directories if they don't exist
        
    Returns:
        Dictionary of standardized paths
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    # Define all subdirectories following the established structure
    dirs = {
        'fastas': os.path.join(batch_path, "fastas"),
        'hhsearch': os.path.join(batch_path, "hhsearch"),
        'hhsearch_profiles': os.path.join(batch_path, "hhsearch", "profiles"),
        'blast_chain': os.path.join(batch_path, "blast", "chain"),
        'blast_domain': os.path.join(batch_path, "blast", "domain"),
        'domains': os.path.join(batch_path, "domains"),
    }
    
    # Create directories if requested
    if create_dirs:
        for directory in dirs.values():
            try:
                os.makedirs(directory, exist_ok=True)
                logger.debug(f"Ensured directory exists: {directory}")
            except Exception as e:
                logger.warning(f"Failed to create directory {directory}: {str(e)}")
    
    # Define paths dictionary with standardized naming
    paths = {
        # FASTA files
        'fasta': os.path.join(dirs['fastas'], f"{pdb_chain}.fa"),
        
        # HHblits profile files
        'a3m': os.path.join(dirs['hhsearch_profiles'], f"{pdb_chain}.a3m"),
        'hhm': os.path.join(dirs['hhsearch_profiles'], f"{pdb_chain}.hhm"),
        
        # HHSearch result files
        'hhr': os.path.join(dirs['hhsearch'], f"{pdb_chain}.{ref_version}.hhr"),
        'hh_xml': os.path.join(dirs['hhsearch'], f"{pdb_chain}.{ref_version}.xml"),
        
        # BLAST result files
        'chain_blast': os.path.join(dirs['blast_chain'], f"{pdb_chain}.{ref_version}.xml"),
        'domain_blast': os.path.join(dirs['blast_domain'], f"{pdb_chain}.{ref_version}.xml"),
        
        # Domain summary file
        'domain_summary': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.domain_summary.xml"),
        'blast_only_summary': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.blast_only.domain_summary.xml"),
    }
    
    return paths

def get_batch_path(base_dir: str, batch_name: str) -> str:
    """Get standardized batch path
    
    Args:
        base_dir: Base directory for ECOD data
        batch_name: Name of the batch
        
    Returns:
        Standardized batch path
    """
    return os.path.join(base_dir, "batches", batch_name)

def get_file_db_path(batch_path: str, file_path: str) -> str:
    """Convert absolute file path to relative path for database storage
    
    Args:
        batch_path: Base path of the batch
        file_path: Absolute file path
        
    Returns:
        Relative path for database storage
    """
    try:
        return os.path.relpath(file_path, batch_path)
    except Exception as e:
        logger.warning(f"Failed to get relative path for {file_path}: {str(e)}")
        return file_path  # Return full path if relative path fails

def resolve_file_path(batch_path: str, db_path: str) -> str:
    """Resolve database-stored relative path to absolute path
    
    Args:
        batch_path: Base path of the batch
        db_path: Relative path from database
        
    Returns:
        Absolute file path
    """
    if os.path.isabs(db_path):
        return db_path
    return os.path.join(batch_path, db_path)

def find_files_with_legacy_paths(batch_path: str, pdb_id: str, chain_id: str,
                              ref_version: str) -> Dict[str, Dict[str, str]]:
    """Find files using both standard and legacy path patterns
    
    Args:
        batch_path: Path to batch directory
        pdb_id: PDB ID
        chain_id: Chain ID
        ref_version: Reference version
        
    Returns:
        Dictionary with file types and paths for both standard and legacy formats
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    results = {}
    
    # Define legacy patterns
    legacy_patterns = {
        'hhr': [
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.hhsearch.{ref_version}.hhr"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}_{ref_version}.hhr"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.{ref_version}.hhr")
        ],
        'hh_xml': [
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.hhsearch.{ref_version}.xml"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.{ref_version}.hhsearch.xml"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}_{ref_version}.xml")
        ],
        'chain_blast': [
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.chainwise_blast.xml"),
            os.path.join(batch_path, "blast", "chain", f"{pdb_chain}.chainwise.{ref_version}.xml"),
        ],
        'domain_blast': [
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.blast.xml"),
            os.path.join(batch_path, "blast", "domain", f"{pdb_chain}.{ref_version}.blast.xml"),
        ],
        'domain_summary': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.domains.xml"),
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.domain_summary.xml"),
        ]
    }
    
    # Get standard paths
    standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)
    
    # Check both standard and legacy paths
    for file_type, standard_path in standard_paths.items():
        results[file_type] = {
            'standard': standard_path,
            'legacy': None,
            'exists_at': None
        }
        
        # Check standard path
        if os.path.exists(standard_path):
            results[file_type]['exists_at'] = standard_path
            continue
            
        # Check legacy paths if file type has legacy patterns
        if file_type in legacy_patterns:
            for legacy_path in legacy_patterns[file_type]:
                if os.path.exists(legacy_path):
                    results[file_type]['legacy'] = legacy_path
                    results[file_type]['exists_at'] = legacy_path
                    break
    
    return results

def migrate_file_to_standard_path(src_path: str, dst_path: str) -> bool:
    """Migrate a file from legacy path to standard path
    
    Args:
        src_path: Source (legacy) path
        dst_path: Destination (standard) path
        
    Returns:
        True if migration was successful
    """
    if not os.path.exists(src_path):
        logger.warning(f"Source file does not exist: {src_path}")
        return False
        
    if os.path.exists(dst_path):
        logger.warning(f"Destination file already exists: {dst_path}")
        return False
        
    try:
        # Create destination directory if it doesn't exist
        dst_dir = os.path.dirname(dst_path)
        os.makedirs(dst_dir, exist_ok=True)
        
        # Copy the file (don't move to avoid data loss)
        import shutil
        shutil.copy2(src_path, dst_path)
        logger.info(f"Migrated file: {src_path} -> {dst_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to migrate file {src_path} -> {dst_path}: {str(e)}")
        return False

def update_db_file_path(db_connection, process_file_id: int, new_path: str) -> bool:
    """Update file path in the database
    
    Args:
        db_connection: Database connection
        process_file_id: ID of the process_file record
        new_path: New file path
        
    Returns:
        True if update was successful
    """
    try:
        db_connection.update(
            "ecod_schema.process_file",
            {
                "file_path": new_path,
                "last_checked": "CURRENT_TIMESTAMP"
            },
            "id = %s",
            (process_file_id,)
        )
        logger.info(f"Updated database path for file ID {process_file_id} to {new_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to update database path for file ID {process_file_id}: {str(e)}")
        return False
