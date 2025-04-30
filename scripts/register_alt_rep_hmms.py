#!/usr/bin/env python3
"""
register_alt_rep_profiles.py - Find, copy and register alternative representative profile files in the database
using standardized path handling
"""

import os
import sys
import logging
import argparse
import re
import shutil
import glob
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import xml.etree.ElementTree as ET
from xml.dom import minidom

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.pipelines.hhsearch.processor import HHRToXMLConverter
from ecod.utils.hhsearch_utils import HHRParser
# Import standardized path utilities
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    find_files_with_legacy_paths,
    extract_pdb_chain_from_path
)

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
        b.ref_version,
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
    """Search for profile file (HMM, A3M, or HHSEARCH) in source directories"""
    logger = logging.getLogger("ecod.alt_rep_profiles")

    # Try looking for an exact match first
    for source_dir in source_dirs:
        path = os.path.join(source_dir, f"{protein_id}{file_ext}")
        if os.path.exists(path):
            logger.debug(f"Found exact match for {protein_id}{file_ext}: {path}")
            return path

    # If no exact match, try more flexible pattern matching
    pdb_id, chain_id = protein_id.split('_', 1)
    possible_matches = []

    # Check in numbered subdirectories (0-9)
    for source_dir in source_dirs:
        for subdir in range(10):  # Check subdirectories 0-9
            search_dir = os.path.join(source_dir, str(subdir))
            if os.path.exists(search_dir):
                # Pattern for the specific file type
                pattern = os.path.join(search_dir, f"{pdb_id}_{chain_id}*{file_ext}")
                matches = glob.glob(pattern)

                if matches:
                    # Sort by filename length - shorter names likely more precise matches
                    matches.sort(key=lambda x: len(os.path.basename(x)))
                    possible_matches.extend(matches)

    if possible_matches:
        best_match = possible_matches[0]
        logger.info(f"Found best match for {protein_id}{file_ext}: {best_match}")
        return best_match

    logger.warning(f"No {file_ext} file found for {protein_id}")
    return None

def copy_profile_to_batch(source_path: str, protein: Dict[str, Any], file_type: str, dry_run: bool = False) -> Optional[str]:
    """Copy profile file to batch directory using standardized paths"""
    logger = logging.getLogger("ecod.alt_rep_profiles")

    # Get standardized paths for this protein
    batch_path = protein['base_path']
    pdb_id = protein['pdb_id']
    chain_id = protein['chain_id']
    ref_version = protein['ref_version']

    # Get standardized paths
    paths = get_standardized_paths(
        batch_path=batch_path,
        pdb_id=pdb_id,
        chain_id=chain_id,
        ref_version=ref_version,
        create_dirs=not dry_run
    )

    # Determine which path to use based on file_type
    if file_type == 'hhm' or file_type == 'hmm':
        dest_path = paths['hhm']
    elif file_type == 'a3m':
        dest_path = paths['a3m']
    elif file_type == 'hhsearch':
        dest_path = paths['hhr']
    else:
        logger.error(f"Unknown file type: {file_type}")
        return None

    if dry_run:
        logger.info(f"[DRY RUN] Would copy {source_path} to {dest_path}")
        return dest_path

    try:
        # Check if source and destination are the same file
        if os.path.exists(dest_path) and os.path.samefile(source_path, dest_path):
            logger.info(f"Source and destination are the same file: {dest_path}")
            return dest_path

        # Ensure directory exists (should be created by get_standardized_paths)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)

        shutil.copy2(source_path, dest_path)
        logger.info(f"Copied {source_path} to {dest_path}")
        return dest_path
    except Exception as e:
        logger.error(f"Failed to copy {source_path} to {dest_path}: {str(e)}")
        return None

def register_profile_in_db(db, process_id: int, batch_path: str, file_path: str, file_type: str, dry_run: bool = False) -> bool:
    """Register profile file in the process_file table using relative paths"""
    logger = logging.getLogger("ecod.alt_rep_profiles")

    try:
        # Make sure we use the correct file_type string for database registration
        # This addresses the inconsistency between file extensions and db types
        db_file_type = file_type
        if file_type == 'hmm':
            db_file_type = 'hhm'  # Use 'hhm' in the database for both .hmm and .hhm files

        # Convert to relative path for database storage
        rel_path = get_file_db_path(batch_path, file_path)

        # Check if file already registered
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = %s
        """

        existing = db.execute_query(check_query, (process_id, db_file_type))

        if dry_run:
            if existing:
                logger.info(f"[DRY RUN] Would update existing {db_file_type} record for process {process_id}")
            else:
                logger.info(f"[DRY RUN] Would create new {db_file_type} record for process {process_id}")
            return True

        # Verify file exists before registering
        file_exists = os.path.exists(file_path)
        file_size = os.path.getsize(file_path) if file_exists else 0

        if not file_exists:
            logger.warning(f"Attempting to register non-existent file: {file_path}")

        if existing:
            # Update existing record
            db.update(
                "ecod_schema.process_file",
                {
                    "file_path": rel_path,
                    "file_exists": file_exists,
                    "file_size": file_size,
                    "last_checked": "NOW()"
                },
                "id = %s",
                (existing[0][0],)
            )
            logger.debug(f"Updated existing {db_file_type} record for process {process_id}")
        else:
            # Insert new record
            db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": db_file_type,
                    "file_path": rel_path,
                    "file_exists": file_exists,
                    "file_size": file_size,
                    "last_checked": "NOW()"
                }
            )
            logger.debug(f"Created new {db_file_type} record for process {process_id}")

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

def process_profiles_for_protein(db, protein: Dict[str, Any], source_dirs: List[str],
                                file_types: List[str], dry_run: bool = False) -> Dict[str, bool]:
    """Process profile files for a single protein"""
    logger = logging.getLogger("ecod.alt_rep_profiles")

    protein_id = protein['protein_id']
    process_id = protein['process_id']
    source_id = f"{protein['pdb_id']}_{protein['chain_id']}"
    base_path = protein['base_path']
    batch_id = None
    classification = protein.get('classification', 'unknown')

    results = {}

    for file_type in file_types:
        logger.debug(f"Processing {file_type} for alt rep {source_id} (protein_id={protein_id}, process_id={process_id}, classification={classification})")

        # Determine file extension to search for
        if file_type == 'hhm':
            # Search for both .hhm and .hmm files
            source_path = find_profile_file(source_dirs, source_id, '.hhm')
            if not source_path:
                source_path = find_profile_file(source_dirs, source_id, '.hmm')
        elif file_type == 'hhsearch':
            # Search for .hhsearch files
            source_path = find_profile_file(source_dirs, source_id, '.hhsearch')
        else:
            source_path = find_profile_file(source_dirs, source_id, f'.{file_type}')

        if not source_path:
            logger.error(f"{file_type.upper()} file not found for {source_id}")
            results[file_type] = False
            continue

        # Special processing for hhsearch files
        if file_type == 'hhsearch':
            success = process_hhsearch_file(db, protein, source_path, dry_run)
            results[file_type] = success
            continue

        # Copy profile file to batch directory
        dest_path = copy_profile_to_batch(source_path, protein, file_type, dry_run)

        if not dest_path:
            results[file_type] = False
            if not dry_run:
                update_process_status(db, process_id, False, file_type, f"Failed to copy {file_type.upper()} file to batch directory", dry_run)
            continue

        # Register profile in database
        db_success = register_profile_in_db(db, process_id, base_path, dest_path, file_type, dry_run)

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

# Process HHSearch file function
def process_hhsearch_file(db, protein: Dict[str, Any], source_path: str, dry_run: bool = False) -> bool:
    """Process HHSearch result file, convert to XML, and register in database"""
    logger = logging.getLogger("ecod.alt_rep_profiles")

    process_id = protein['process_id']
    pdb_id = protein['pdb_id']
    chain_id = protein['chain_id']
    source_id = f"{pdb_id}_{chain_id}"
    base_path = protein['base_path']
    ref_version = protein['ref_version']

    # Get standardized paths
    paths = get_standardized_paths(
        batch_path=base_path,
        pdb_id=pdb_id,
        chain_id=chain_id,
        ref_version=ref_version,
        create_dirs=not dry_run
    )

    hhr_file = paths['hhr']
    xml_file = paths['hh_xml']

    if dry_run:
        logger.info(f"[DRY RUN] Would copy {source_path} to {hhr_file}")
        logger.info(f"[DRY RUN] Would register HHR file without XML conversion")
        return True

    # Copy HHSearch file to standard location with standardized name
    try:
        # Ensure directory exists
        os.makedirs(os.path.dirname(hhr_file), exist_ok=True)

        shutil.copy2(source_path, hhr_file)
        logger.info(f"Copied {source_path} to {hhr_file}")
    except Exception as e:
        logger.error(f"Failed to copy HHSearch file: {str(e)}")
        return False

    # Register HHR file in database
    if not register_profile_in_db(db, process_id, base_path, hhr_file, 'hhr', dry_run):
        logger.error(f"Failed to register HHR file in database: {hhr_file}")
        return False

    # Generate a simple XML version of the hhsearch file manually
    try:
        logger.info(f"Converting HHR data to XML for {source_id}")

        # Parse the hhsearch file
        with open(hhr_file, 'r') as f:
            content = f.read()

        # Create XML structure
        root = ET.Element("hh_summ_doc")

        # Add metadata
        metadata = ET.SubElement(root, "metadata")
        ET.SubElement(metadata, "pdb_id").text = pdb_id
        ET.SubElement(metadata, "chain_id").text = chain_id
        ET.SubElement(metadata, "reference").text = ref_version
        ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Add hits element
        hits_elem = ET.SubElement(root, "hh_hit_list")

        # Process the file line by line to extract hits
        lines = content.split('\n')

        # Find the hit table header (No Hit...)
        table_start = None
        for i, line in enumerate(lines):
            if line.startswith(' No Hit'):
                table_start = i + 1
                break

        if table_start:
            # Process hit lines
            i = table_start
            while i < len(lines):
                line = lines[i].strip()
                if not line:
                    break

                # Parse hit line
                match = re.match(r'^\s*(\d+)\s+(\S+)\s+(\S+)\s+(.*)$', line)
                if match:
                    hit_num = match.group(1)
                    hit_id = match.group(2)
                    domain_info = match.group(3)
                    rest = match.group(4)

                    # Create hit element
                    hit_elem = ET.SubElement(hits_elem, "hh_hit")
                    hit_elem.set("hit_num", hit_num)
                    hit_elem.set("hit_id", hit_id)
                    hit_elem.set("domain_info", domain_info)

                    # Parse values from rest of line
                    values = rest.split()
                    if len(values) >= 6:
                        try:
                            prob = float(values[0])
                            e_value = values[1]
                            p_value = values[2]
                            score = values[3]
                            ss = values[4]
                            cols = values[5]

                            hit_elem.set("probability", str(prob))
                            hit_elem.set("e_value", e_value)
                            hit_elem.set("p_value", p_value)
                            hit_elem.set("score", score)
                            hit_elem.set("ss", ss)
                            hit_elem.set("cols", cols)

                            # Extract query and template ranges if available
                            if len(values) >= 8:
                                query_range = values[6]
                                template_range = values[7]

                                query_range_elem = ET.SubElement(hit_elem, "query_range")
                                query_range_elem.text = query_range

                                template_range_elem = ET.SubElement(hit_elem, "template_range")
                                template_range_elem.text = template_range
                        except Exception as e:
                            logger.warning(f"Error parsing hit values: {str(e)}")

                i += 1

        # Convert to string
        rough_string = ET.tostring(root, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")

        # Save XML file
        with open(xml_file, 'w') as f:
            f.write(pretty_xml)

        logger.info(f"Successfully created XML file: {xml_file}")

        # Register XML file in database
        if not register_profile_in_db(db, process_id, base_path, xml_file, 'hh_xml', dry_run):
            logger.error(f"Failed to register XML file in database: {xml_file}")
            return False

        logger.info(f"Successfully processed HHSearch file for {source_id}")
        update_process_status(db, process_id, True, 'hhsearch', None, dry_run)

        return True
    except Exception as e:
        logger.error(f"Error converting HHSearch file to XML: {str(e)}")
        logger.exception("Full traceback:")

        # We still successfully registered the HHR file, so return true
        # but with a warning in the status
        update_process_status(db, process_id, True, 'hhsearch',
                            f"Registered HHR file but XML conversion failed: {str(e)}", dry_run)
        return True

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
    parser.add_argument('--no-hhsearch', action="store_true",
                      help="Skip .hhsearch files")
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
    if not args.no_hhsearch:
        file_types.append('hhsearch')

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
