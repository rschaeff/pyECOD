# ecod/utils/path_utils.py
"""
Standardized file path utilities for pyECOD

This module provides consistent path construction for all file types
used in the pyECOD processing pipeline, ensuring standardization across
different components.
"""

import os
import logging
import re
import shutil
from typing import Dict, Optional, List, Tuple, Any

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

        # Domain files - clear distinction between summary and partition
        'domain_summary': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.domain_summary.xml"),
        'domain_partition': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.domains.xml"),

        # Blast-only variants
        'blast_only_summary': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.blast_only.domain_summary.xml"),
        'blast_only_partition': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.blast_only.domains.xml"),
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

def get_file_type_from_path(file_path: str) -> Optional[str]:
    """Determine file type from file path

    Args:
        file_path: Path to file

    Returns:
        File type identifier or None if unknown
    """
    basename = os.path.basename(file_path).lower()

    # Determine type by extension and naming pattern
    if basename.endswith('.fa'):
        return 'fasta'
    elif basename.endswith('.a3m'):
        return 'a3m'
    elif basename.endswith('.hhm'):
        return 'hhm'
    elif basename.endswith('.hhr'):
        return 'hhr'
    elif basename.endswith('.domains.xml'):
        if 'blast_only' in basename:
            return 'blast_only_partition'
        else:
            return 'domain_partition'
    elif basename.endswith('.domain_summary.xml'):
        if 'blast_only' in basename:
            return 'blast_only_summary'
        else:
            return 'domain_summary'
    elif basename.endswith('.hhsearch.xml') or basename.endswith('.hh.xml'):
        return 'hh_xml'
    elif basename.endswith('.xml') and '/blast/chain/' in file_path.replace('\\', '/'):
        return 'chain_blast'
    elif basename.endswith('.xml') and '/blast/domain/' in file_path.replace('\\', '/'):
        return 'domain_blast'
    elif basename.endswith('.chainwise_blast.xml'):
        return 'chain_blast'
    elif basename.endswith('.blast.xml'):
        return 'domain_blast'

    # Default case - unknown file type
    return None

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

    # Define legacy patterns with clear distinction between domain_summary and domain_partition
    legacy_patterns = {
        'hhr': [
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.hhsearch.{ref_version}.hhr"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}_{ref_version}.hhr"),
            os.path.join(batch_path, pdb_chain, f"{pdb_chain}.{ref_version}.hhr"),
            os.path.join(batch_path, "ecod_dump", pdb_chain, f"{pdb_chain}.{ref_version}.hhr"),
            os.path.join(batch_path, "ecod_dump", f"{pdb_chain}.{ref_version}.hhr")
        ],
        'a3m': [
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.a3m")
        ],
        'hhm': [
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.hhm")
        ],
        'hh_xml': [
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.hhsearch.{ref_version}.xml"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.{ref_version}.hhsearch.xml"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}_{ref_version}.xml")
        ],
        'chain_blast': [
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.chainwise_blast.xml"),
            os.path.join(batch_path, "blast", "chain", f"{pdb_chain}.chainwise.{ref_version}.xml"),
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.chain_blast.xml"),
            os.path.join(batch_path, "blast", "chain", "batch_0", f"{pdb_chain}.chainwise_blast.xml"),
            os.path.join(batch_path, "blast", "chain", "batch_1", f"{pdb_chain}.chainwise_blast.xml"),
            os.path.join(batch_path, "blast", "chain", "batch_2", f"{pdb_chain}.chainwise_blast.xml")
        ],
        'domain_blast': [
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.blast.xml"),
            os.path.join(batch_path, "blast", "domain", f"{pdb_chain}.{ref_version}.blast.xml"),
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.domain_blast.xml"),
            os.path.join(batch_path, "blast", "domain", "batch_0", f"{pdb_chain}.domain_blast.xml"),
            os.path.join(batch_path, "blast", "domain", "batch_1", f"{pdb_chain}.domain_blast.xml"),
            os.path.join(batch_path, "blast", "domain", "batch_2", f"{pdb_chain}.domain_blast.xml")
        ],
        'domain_summary': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.domain_summary.xml"),
        ],
        'domain_partition': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.domains.xml"),
        ],
        'blast_only_summary': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.blast_only.domain_summary.xml"),
        ],
        'blast_only_partition': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.blast_only.domains.xml"),
        ]
    }

    # Get standard paths
    standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)

    # In find_files_with_legacy_paths function in path_utils.py
    logger.info(f"Looking for legacy paths for {pdb_id}_{chain_id}")
    logger.info(f"Legacy pattern for chain_blast: {legacy_patterns['chain_blast']}")
    logger.info(f"Legacy pattern for domain_blast: {legacy_patterns['domain_blast']}")

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
        # If files are the same, return success immediately
        if os.path.samefile(src_path, dst_path):
            logger.debug(f"Source and destination are the same file: {src_path}")
            return True

        # Check if the destination is newer
        if os.path.getmtime(dst_path) > os.path.getmtime(src_path):
            logger.info(f"Destination file is newer than source, keeping destination: {dst_path}")
            return True

        logger.warning(f"Destination file already exists: {dst_path}")

        # If destination already exists but is older, we'll replace it
        # But first check file sizes to ensure we're not replacing with empty file
        if os.path.getsize(src_path) < os.path.getsize(dst_path):
            logger.warning(f"Source file is smaller than destination, keeping destination: {dst_path}")
            return False

    try:
        # Create destination directory if it doesn't exist
        dst_dir = os.path.dirname(dst_path)
        os.makedirs(dst_dir, exist_ok=True)

        # Copy the file with metadata preserved
        shutil.copy2(src_path, dst_path)
        logger.info(f"Migrated file: {src_path} -> {dst_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to migrate file {src_path} -> {dst_path}: {str(e)}")
        return False

def update_db_file_paths(db_connection, batch_id: int, dry_run: bool = True) -> Dict[str, int]:
    """Update file paths in the database to match standardized paths

    Args:
        db_connection: Database connection
        batch_id: Batch ID
        dry_run: Whether to perform a dry run without updating the database

    Returns:
        Dictionary with statistics
    """
    logger.info(f"Updating database file paths for batch {batch_id} (dry_run={dry_run})")

    # Get batch information
    query = """
    SELECT id, batch_name, base_path, ref_version
    FROM ecod_schema.batch
    WHERE id = %s
    """

    batch_info = db_connection.execute_dict_query(query, (batch_id,))
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return {"error": 1}

    batch_path = batch_info[0]['base_path']
    ref_version = batch_info[0]['ref_version']

    # Get all process files for this batch
    query = """
    SELECT
        pf.id, pf.process_id, pf.file_type, pf.file_path, pf.file_exists,
        p.pdb_id, p.chain_id
    FROM
        ecod_schema.process_file pf
    JOIN
        ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE
        ps.batch_id = %s
    """

    process_files = db_connection.execute_dict_query(query, (batch_id,))
    logger.info(f"Found {len(process_files)} files in database for batch {batch_id}")

    # Statistics
    stats = {
        'total': len(process_files),
        'updated': 0,
        'already_standard': 0,
        'errors': 0,
        'file_missing': 0
    }

    # Process each file
    for file in process_files:
        try:
            pdb_id = file['pdb_id']
            chain_id = file['chain_id']
            file_type = file['file_type']
            current_path = file['file_path']

            # Resolve current absolute path
            current_abs_path = resolve_file_path(batch_path, current_path)

            # Get standardized paths
            paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)

            # Skip if file type not in standardized paths
            if file_type not in paths:
                logger.warning(f"Unknown file type '{file_type}' for {pdb_id}_{chain_id}")
                stats['errors'] += 1
                continue

            # Get standard path
            standard_path = paths[file_type]
            standard_rel_path = get_file_db_path(batch_path, standard_path)

            # Skip if already using standard path
            if current_path == standard_rel_path:
                logger.debug(f"Already using standard path: {current_path}")
                stats['already_standard'] += 1
                continue

            # Check if file exists at current path
            if not os.path.exists(current_abs_path):
                logger.warning(f"File does not exist at current path: {current_abs_path}")

                # Check if file exists at standard path
                if os.path.exists(standard_path):
                    logger.info(f"File exists at standard path: {standard_path}")
                else:
                    logger.warning(f"File missing at both current and standard paths")
                    stats['file_missing'] += 1
                    continue
            else:
                # Migrate file to standard path if needed
                if not os.path.exists(standard_path):
                    if not dry_run:
                        migrate_file_to_standard_path(current_abs_path, standard_path)

            # Update database path
            if not dry_run:
                db_connection.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": standard_rel_path
                    },
                    "id = %s",
                    (file['id'],)
                )
                logger.info(f"Updated path in database: {current_path} -> {standard_rel_path}")
            else:
                logger.info(f"Would update path: {current_path} -> {standard_rel_path}")

            stats['updated'] += 1

        except Exception as e:
            logger.error(f"Error processing file {file['id']}: {str(e)}")
            stats['errors'] += 1

    # Log statistics
    logger.info("Update Statistics:")
    logger.info(f"Total files processed: {stats['total']}")
    logger.info(f"Files already using standard paths: {stats['already_standard']}")
    logger.info(f"Files updated: {stats['updated']}")
    logger.info(f"Files missing: {stats['file_missing']}")
    logger.info(f"Errors: {stats['errors']}")

    return stats

def list_files_for_chain(batch_path: str, pdb_id: str, chain_id: str,
                       ref_version: str) -> Dict[str, Dict[str, Any]]:
    """List all files for a chain with their paths and existence status

    Args:
        batch_path: Path to batch directory
        pdb_id: PDB ID
        chain_id: Chain ID
        ref_version: Reference version

    Returns:
        Dictionary with file types, paths, and existence status
    """
    # Get standardized paths
    standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)

    # Find files with legacy paths
    file_info = find_files_with_legacy_paths(batch_path, pdb_id, chain_id, ref_version)

    # Build result with all information
    result = {}
    for file_type, paths in file_info.items():
        standard_path = paths['standard']
        legacy_path = paths['legacy']
        exists_at = paths['exists_at']

        result[file_type] = {
            'standard_path': standard_path,
            'legacy_path': legacy_path,
            'exists_at': exists_at,
            'exists': exists_at is not None,
            'size': os.path.getsize(exists_at) if exists_at and os.path.exists(exists_at) else 0,
            'status': 'missing' if not exists_at else 'standard' if exists_at == standard_path else 'legacy'
        }

    return result

def scan_batch_directory(batch_path: str, ref_version: str) -> Dict[str, List[str]]:
    """Scan batch directory for files and classify them

    Args:
        batch_path: Path to batch directory
        ref_version: Reference version

    Returns:
        Dictionary with file types and lists of paths
    """
    result = {
        'fasta': [],
        'a3m': [],
        'hhm': [],
        'hhr': [],
        'hh_xml': [],
        'chain_blast': [],
        'domain_blast': [],
        'domain_summary': [],
        'domain_partition': [],
        'blast_only_summary': [],
        'blast_only_partition': [],
        'unknown': []
    }

    # Walk through the batch directory
    for root, dirs, files in os.walk(batch_path):
        for file in files:
            file_path = os.path.join(root, file)

            # Determine file type
            file_type = get_file_type_from_path(file_path)

            if file_type:
                result[file_type].append(file_path)
            else:
                result['unknown'].append(file_path)

    return result

def extract_pdb_chain_from_path(file_path: str) -> Optional[Tuple[str, str]]:
    """Extract PDB ID and chain ID from file path

    Args:
        file_path: Path to file

    Returns:
        Tuple of (pdb_id, chain_id) or None if not found
    """
    basename = os.path.basename(file_path)

    # Pattern: pdbid_chainid.* (e.g., 1abc_A.fa, 1abc_A.develop291.hhr)
    match = re.match(r'([a-zA-Z0-9]{4})_([a-zA-Z0-9])\.', basename)
    if match:
        return match.group(1), match.group(2)

    return None

# Enhance ecod/utils/path_utils.py

def get_all_evidence_paths(batch_path: str, pdb_id: str, chain_id: str, ref_version: str) -> Dict[str, Dict[str, str]]:
    """
    Get all paths related to evidence processing for a protein chain

    Returns:
        Dictionary with file types and paths info including:
        - standard_path: Path using standard naming
        - legacy_path: Path using legacy naming (if applicable)
        - exists_at: Path where file exists (or None)
    """
    # Get standard paths
    standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)

    # Check for legacy paths
    legacy_files = find_files_with_legacy_paths(batch_path, pdb_id, chain_id, ref_version)

    # Combine information
    result = {}
    for file_type, standard_path in standard_paths.items():
        result[file_type] = {
            'standard_path': standard_path,
            'legacy_path': None,
            'exists_at': None
        }

        # Check if file exists at standard path
        if os.path.exists(standard_path):
            result[file_type]['exists_at'] = standard_path
        elif file_type in legacy_files and legacy_files[file_type].get('exists_at'):
            # Use legacy path if standard doesn't exist
            # Safely access 'legacy_path' with .get() method to avoid KeyError
            result[file_type]['legacy_path'] = legacy_files[file_type].get('legacy_path')
            result[file_type]['exists_at'] = legacy_files[file_type]['exists_at']

    return result
