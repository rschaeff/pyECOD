#!/usr/bin/env python3
"""
regenerate_missing_summaries.py - Create missing domain summary files and fix database status

This script identifies proteins with missing domain summary files, determines if they
are no-hits cases or peptides (very short), generates appropriate stub XML files,
and updates the database to properly reflect their status.
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
import glob
import re
from datetime import datetime
from typing import Dict, Any, Optional, List, Tuple

# Add parent directory to path to allow imports from ecod modules
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Now we can import ecod modules
from ecod.core.context import ApplicationContext
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_file_db_path,
    resolve_file_path,
    find_files_with_legacy_paths,
    migrate_file_to_standard_path
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

def get_batch_info(context: ApplicationContext, batch_id: int) -> Dict[str, Any]:
    """
    Get batch information from database

    Args:
        context: Application context
        batch_id: Batch ID to retrieve

    Returns:
        Dictionary with batch information or empty dict if not found
    """
    query = """
    SELECT
        id, batch_name, base_path, ref_version,
        total_items, completed_items, status
    FROM
        ecod_schema.batch
    WHERE
        id = %s
    """

    result = context.db.execute_query(query, (batch_id,))
    if not result:
        return {}

    return {
        'id': result[0][0],
        'name': result[0][1],
        'path': result[0][2],
        'reference': result[0][3],
        'total_items': result[0][4],
        'completed_items': result[0][5],
        'status': result[0][6]
    }

def get_missing_domain_summaries(context: ApplicationContext, batch_id: int) -> List[Dict[str, Any]]:
    """
    Get information about proteins with missing domain summary files

    Args:
        context: Application context
        batch_id: Batch ID to process

    Returns:
        List of dictionaries with details about proteins with missing files
    """
    query = """
    SELECT
        pf.id AS file_id,
        pf.process_id,
        pf.file_path,
        ps.current_stage,
        ps.status,
        p.pdb_id,
        p.chain_id,
        p.length,
        ps.id AS process_status_id
    FROM
        ecod_schema.process_file pf
    JOIN
        ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE
        ps.batch_id = %s
        AND pf.file_type = 'domain_summary'
        AND (pf.file_exists = FALSE OR pf.file_exists IS NULL)
    ORDER BY p.pdb_id, p.chain_id
    """

    results = context.db.execute_dict_query(query, (batch_id,))
    if not results:
        return []

    return results

def check_blast_hits_in_filesystem(batch_path: str, pdb_id: str, chain_id: str,
                                 ref_version: str) -> bool:
    """
    Check if a protein has any BLAST hits by looking at the blast result files

    Args:
        batch_path: Base path of the batch
        pdb_id: PDB identifier
        chain_id: Chain identifier
        ref_version: Reference version

    Returns:
        True if blast hit files exist and contain hits, False otherwise
    """
    logger = logging.getLogger("ecod.regenerate.blast_check")

    # Get standardized paths
    paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)

    # Check for chain and domain blast files
    chain_blast_path = paths['chain_blast']
    domain_blast_path = paths['domain_blast']

    # Find files using both standard and legacy paths
    chain_blast_info = find_files_with_legacy_paths(batch_path, pdb_id, chain_id, ref_version)
    chain_blast_exists = chain_blast_info['chain_blast']['exists_at'] is not None
    domain_blast_exists = chain_blast_info['domain_blast']['exists_at'] is not None

    if not chain_blast_exists and not domain_blast_exists:
        logger.debug(f"No BLAST files found for {pdb_id}_{chain_id}")
        return False

    # Get actual file paths if they exist
    actual_chain_blast = chain_blast_info['chain_blast']['exists_at']
    actual_domain_blast = chain_blast_info['domain_blast']['exists_at']

    # Check content of chain blast file for hits
    has_chain_hits = False
    if actual_chain_blast and os.path.exists(actual_chain_blast):
        try:
            with open(actual_chain_blast, 'r') as f:
                content = f.read()
                # Check if file contains hit elements
                if "<hit" in content and "<hits" in content:
                    # Parse XML to count hits
                    tree = ET.parse(actual_chain_blast)
                    root = tree.getroot()
                    hits = root.findall(".//hit")
                    has_chain_hits = len(hits) > 0
        except Exception as e:
            logger.warning(f"Error checking chain blast file {actual_chain_blast}: {str(e)}")

    # Check content of domain blast file for hits
    has_domain_hits = False
    if actual_domain_blast and os.path.exists(actual_domain_blast):
        try:
            with open(actual_domain_blast, 'r') as f:
                content = f.read()
                # Check if file contains hit elements
                if "<hit" in content and "<hits" in content:
                    # Parse XML to count hits
                    tree = ET.parse(actual_domain_blast)
                    root = tree.getroot()
                    hits = root.findall(".//hit")
                    has_domain_hits = len(hits) > 0
        except Exception as e:
            logger.warning(f"Error checking domain blast file {actual_domain_blast}: {str(e)}")

    # Return True if either chain or domain blast has hits
    return has_chain_hits or has_domain_hits

def create_no_hits_xml(pdb_id: str, chain_id: str, is_peptide: bool = False,
                    ref_version: str = "develop") -> str:
    """
    Create XML content for a no-hits domain summary

    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        is_peptide: Whether this is a peptide (very short protein)
        ref_version: Reference version to use

    Returns:
        XML content as string
    """
    # Create structured XML document with proper namespaces
    root = ET.Element("domain_summ_doc")

    # Add metadata
    metadata = ET.SubElement(root, "metadata")
    ET.SubElement(metadata, "pdb_id").text = pdb_id
    ET.SubElement(metadata, "chain_id").text = chain_id
    ET.SubElement(metadata, "reference").text = ref_version
    ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Add empty evidence sections
    chain_blast_elem = ET.SubElement(root, "chain_blast_evidence")
    domain_blast_elem = ET.SubElement(root, "domain_blast_evidence")
    hhsearch_elem = ET.SubElement(root, "hhsearch_evidence")

    # Add status information
    status_elem = ET.SubElement(root, "status")
    status_elem.set("has_chain_blast", "false")
    status_elem.set("has_domain_blast", "false")
    status_elem.set("has_hhsearch", "false")

    if is_peptide:
        status_elem.set("is_peptide", "true")
        status_elem.set("reason", "peptide")
    else:
        status_elem.set("no_hits", "true")
        status_elem.set("reason", "no_hits")

    # Add empty domain suggestions
    ET.SubElement(root, "domain_suggestions")

    # Convert to properly formatted XML
    rough_string = ET.tostring(root, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="  ")

    return pretty_xml

def write_xml_file(file_path: str, content: str, base_path: str) -> bool:
    """
    Write XML content to file, ensuring directory exists

    Args:
        file_path: Relative or absolute path where to write the file
        content: XML content to write
        base_path: Base path to resolve relative paths

    Returns:
        True if file was written successfully
    """
    logger = logging.getLogger("ecod.regenerate.file_writer")

    # Handle absolute and relative paths
    if os.path.isabs(file_path):
        full_path = file_path
    else:
        full_path = resolve_file_path(base_path, file_path)

    # Ensure directory exists
    try:
        os.makedirs(os.path.dirname(full_path), exist_ok=True)
        logger.debug(f"Ensured directory exists: {os.path.dirname(full_path)}")
    except Exception as e:
        logger.error(f"Failed to create directory for {full_path}: {str(e)}")
        return False

    try:
        with open(full_path, 'w', encoding='utf-8') as f:
            f.write(content)
        logger.debug(f"Successfully wrote file: {full_path}")
        return True
    except Exception as e:
        logger.error(f"Error writing file {full_path}: {str(e)}")
        return False

def update_database_status(context: ApplicationContext, protein_info: Dict[str, Any],
                         file_exists: bool = True, is_complete: bool = True) -> bool:
    """
    Update database status for a protein

    Args:
        context: Application context
        protein_info: Dictionary with protein information
        file_exists: Whether the file exists
        is_complete: Whether processing is complete

    Returns:
        True if database was updated successfully
    """
    logger = logging.getLogger("ecod.regenerate.db_updater")

    file_id = protein_info['file_id']
    process_id = protein_info['process_id']

    try:
        # Update the process_file record
        if file_exists:
            full_path = resolve_file_path(protein_info.get('base_path', ''), protein_info['file_path'])
            if os.path.exists(full_path):
                file_size = os.path.getsize(full_path)
            else:
                file_size = 0

            update_file_query = """
            UPDATE ecod_schema.process_file
            SET file_exists = TRUE, file_size = %s, last_checked = NOW()
            WHERE id = %s
            """
            context.db.execute_query(update_file_query, (file_size, file_id))
            logger.debug(f"Updated process_file record {file_id} (size: {file_size})")
        else:
            update_file_query = """
            UPDATE ecod_schema.process_file
            SET file_exists = FALSE, file_size = 0, last_checked = NOW()
            WHERE id = %s
            """
            context.db.execute_query(update_file_query, (file_id,))
            logger.debug(f"Updated process_file record {file_id} (missing file)")

        # Update the process_status record based on whether it's complete
        if is_complete:
            # For no-hits or peptides that are properly processed
            status_query = """
            UPDATE ecod_schema.process_status
            SET status = 'success', current_stage = 'domain_summary_complete',
                updated_at = NOW(), error_message = NULL
            WHERE id = %s
            """
            stage = "domain_summary_complete"
        else:
            # For cases that still need processing
            status_query = """
            UPDATE ecod_schema.process_status
            SET status = 'pending', current_stage = 'domain_summary',
                updated_at = NOW(), error_message = NULL
            WHERE id = %s
            """
            stage = "domain_summary"

        context.db.execute_query(status_query, (process_id,))
        logger.debug(f"Updated process_status record {process_id} (stage: {stage})")

        return True

    except Exception as e:
        logger.error(f"Error updating database: {str(e)}")
        return False

def update_batch_completion(context: ApplicationContext, batch_id: int) -> bool:
    """
    Update the batch completion count and status

    Args:
        context: Application context
        batch_id: Batch ID to update

    Returns:
        True if batch was updated successfully
    """
    logger = logging.getLogger("ecod.regenerate.batch_updater")

    try:
        # Get the current count of completed items
        count_query = """
        SELECT COUNT(*)
        FROM ecod_schema.process_status
        WHERE batch_id = %s
        AND status = 'success'
        AND current_stage IN ('domain_summary_complete', 'domain_partition_complete')
        """

        count_result = context.db.execute_query(count_query, (batch_id,))
        if not count_result:
            logger.warning(f"Failed to get completed count for batch {batch_id}")
            return False

        completed_count = count_result[0][0]

        # Get total items for the batch
        batch_query = """
        SELECT total_items FROM ecod_schema.batch WHERE id = %s
        """

        batch_result = context.db.execute_query(batch_query, (batch_id,))
        if not batch_result:
            logger.warning(f"Failed to get total items for batch {batch_id}")
            return False

        total_items = batch_result[0][0]

        # Update batch status
        update_query = """
        UPDATE ecod_schema.batch
        SET completed_items = %s,
            status = CASE
                WHEN %s >= total_items THEN 'completed'
                ELSE 'processing'
            END,
            updated_at = NOW()
        WHERE id = %s
        """

        context.db.execute_query(update_query, (completed_count, completed_count, batch_id))

        logger.info(f"Updated batch status: {completed_count}/{total_items} completed")
        return True

    except Exception as e:
        logger.error(f"Error updating batch status: {str(e)}")
        return False

def process_missing_files(context: ApplicationContext, batch_id: int, dry_run: bool = False,
                       reference_version: str = None) -> Tuple[int, int, int]:
    """
    Process missing domain summary files

    Args:
        context: Application context
        batch_id: Batch ID to process
        dry_run: If True, don't make any changes
        reference_version: Reference version to use for XML files

    Returns:
        Tuple of (total processed, files created, database updates)
    """
    logger = logging.getLogger("ecod.regenerate")

    # Get batch information
    batch_info = get_batch_info(context, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return 0, 0, 0

    base_path = batch_info['path']
    logger.info(f"Processing batch {batch_id} ({batch_info['name']}) with base path: {base_path}")

    # Use provided reference version or get from batch
    ref_version = reference_version or batch_info['reference']

    # Get missing domain summaries
    missing_files = get_missing_domain_summaries(context, batch_id)

    if not missing_files:
        logger.info("No missing domain summary files found")
        return 0, 0, 0

    logger.info(f"Found {len(missing_files)} missing domain summary files")

    # Process each missing file
    total_processed = 0
    files_created = 0
    db_updates = 0

    # Categorize proteins
    peptides = []
    no_hits = []
    has_hits = []

    for protein_info in missing_files:
        total_processed += 1

        pdb_id = protein_info['pdb_id']
        chain_id = protein_info['chain_id']
        length = protein_info['length']
        file_path = protein_info['file_path']

        # Add base path to the protein info for use in other functions
        protein_info['base_path'] = base_path

        logger.debug(f"Processing {pdb_id}_{chain_id} (length: {length})")

        # Get standardized paths
        paths = get_standardized_paths(base_path, pdb_id, chain_id, ref_version)

        # Update file path to use standard paths if needed
        if file_path != get_file_db_path(base_path, paths['domain_summary']):
            standard_path = get_file_db_path(base_path, paths['domain_summary'])
            logger.debug(f"Updating file path from {file_path} to {standard_path}")
            protein_info['file_path'] = standard_path

        # Determine if this is a peptide or no-hits case
        is_peptide = length is not None and length < 25  # Adjust threshold as needed

        if is_peptide:
            peptides.append(protein_info)
            protein_info['category'] = 'peptide'
            protein_info['has_hits'] = False
        else:
            # Check for blast hits
            has_blast_hits = check_blast_hits_in_filesystem(base_path, pdb_id, chain_id, ref_version)
            protein_info['has_hits'] = has_blast_hits

            if has_blast_hits:
                has_hits.append(protein_info)
                protein_info['category'] = 'has_hits'
            else:
                no_hits.append(protein_info)
                protein_info['category'] = 'no_hits'

        # Print details for the first few proteins in each category
        if (len(peptides) <= 3 and protein_info in peptides) or \
           (len(no_hits) <= 3 and protein_info in no_hits) or \
           (len(has_hits) <= 3 and protein_info in has_hits) or \
           (total_processed % 100 == 0):
            logger.info(f"Processing {pdb_id}_{chain_id}: length={length}, category={protein_info['category']}, has_hits={protein_info['has_hits']}")

    # Log summary of categorization
    logger.info(f"Categorization summary:")
    logger.info(f"  Total proteins: {total_processed}")
    logger.info(f"  Peptides (too short): {len(peptides)}")
    logger.info(f"  No BLAST hits: {len(no_hits)}")
    logger.info(f"  Has BLAST hits: {len(has_hits)}")

    # Process peptides and no-hits proteins
    proteins_to_process = peptides + no_hits
    logger.info(f"Will create stub XML files for {len(proteins_to_process)} proteins (peptides and no-hits)")

    for protein_info in proteins_to_process:
        pdb_id = protein_info['pdb_id']
        chain_id = protein_info['chain_id']
        file_path = protein_info['file_path']
        is_peptide = protein_info['category'] == 'peptide'

        # Generate appropriate XML content
        xml_content = create_no_hits_xml(pdb_id, chain_id, is_peptide, ref_version)

        # Create the file if not in dry run mode
        if not dry_run:
            if write_xml_file(file_path, xml_content, base_path):
                files_created += 1
                logger.debug(f"Created XML file for {pdb_id}_{chain_id}")

            # Update database status - mark as complete for peptides and no-hits
            if update_database_status(context, protein_info, True, True):
                db_updates += 1
                logger.debug(f"Updated database for {pdb_id}_{chain_id}")

    # Log information about proteins with hits
    if has_hits:
        logger.warning(f"Found {len(has_hits)} proteins with BLAST hits but missing domain summary files")
        logger.warning("These proteins should have domain summary files. They may need further investigation.")
        logger.warning("First 5 examples:")
        for protein_info in has_hits[:5]:
            logger.warning(f"  {protein_info['pdb_id']}_{protein_info['chain_id']} (length: {protein_info['length']})")

    # Update batch completion status if files were created
    if files_created > 0 and not dry_run:
        update_batch_completion(context, batch_id)

    # Log summary
    logger.info(f"Processing summary:")
    logger.info(f"  Total proteins processed: {total_processed}")

    if dry_run:
        logger.info("  This was a dry run - no changes were made")
        logger.info("  Run without --dry-run to create files and update database")
    else:
        logger.info(f"  Files created: {files_created}")
        logger.info(f"  Database records updated: {db_updates}")

    return total_processed, files_created, db_updates

def main():
    """Main function to regenerate missing domain summary files"""
    parser = argparse.ArgumentParser(description='Regenerate missing domain summary files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--dry-run', action='store_true',
                      help='Check files but don\'t make changes')
    parser.add_argument('--reference-version', type=str,
                      help='Reference version to use (defaults to batch reference)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.regenerate")

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    # Initialize application context with the config file path
    context = ApplicationContext(args.config)

    logger.info(f"Starting regeneration of missing domain summary files for batch {args.batch_id}")

    total, created, updated = process_missing_files(
        context,
        args.batch_id,
        args.dry_run,
        args.reference_version
    )

    if created > 0 or args.dry_run:
        logger.info(f"Successfully completed regeneration process")
    else:
        logger.warning(f"No files were created during regeneration process")

    # Exit with success
    return 0

if __name__ == "__main__":
    sys.exit(main())
