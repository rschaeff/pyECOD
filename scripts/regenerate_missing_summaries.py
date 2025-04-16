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
import glob
import re
from typing import Dict, Any, Optional, List, Tuple

# Add parent directory to path to allow imports from ecod modules
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Now we can import ecod modules
from ecod.core.context import ApplicationContext

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
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

def get_batch_info(context, batch_id: int) -> Dict[str, Any]:
    """Get batch information from database"""
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

def get_missing_domain_summaries(context, batch_id: int) -> List[Dict[str, Any]]:
    """
    Get information about proteins with missing domain summary files
    
    Returns a list of dictionaries with details about proteins with missing files
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
    
    results = context.db.execute_query(query, (batch_id,))
    
    missing_files = []
    for row in results:
        missing_files.append({
            'file_id': row[0],
            'process_id': row[1],
            'file_path': row[2],
            'current_stage': row[3],
            'status': row[4],
            'pdb_id': row[5],
            'chain_id': row[6],
            'length': row[7],
            'process_status_id': row[8]
        })
    
    return missing_files

def check_blast_hits_in_filesystem(batch_path: str, pdb_id: str, chain_id: str) -> bool:
    """
    Check if a protein has any BLAST hits by looking at the blast result files in the filesystem
    
    Args:
        batch_path: Base path of the batch
        pdb_id: PDB identifier
        chain_id: Chain identifier
        
    Returns:
        True if blast hit files exist and contain hits, False otherwise
    """
    # Define paths for chain blast and domain blast results
    chain_blast_pattern = os.path.join(batch_path, 'blast', 'chain', 'batch_*', f'{pdb_id}_{chain_id}*.out')
    domain_blast_pattern = os.path.join(batch_path, 'blast', 'domain', 'batch_*', f'{pdb_id}_{chain_id}*.out')
    
    # Check if any files match the patterns
    chain_blast_files = glob.glob(chain_blast_pattern)
    domain_blast_files = glob.glob(domain_blast_pattern)
    
    # If no files found, there are no hits
    if not chain_blast_files and not domain_blast_files:
        return False
    
    # Check content of chain blast files for hits
    has_chain_hits = False
    for blast_file in chain_blast_files:
        try:
            with open(blast_file, 'r') as f:
                content = f.read()
                # Check if file contains hits - typically a blast hit would have alignment blocks
                # and scores. Here we look for standard BLAST output patterns.
                if re.search(r'Score\s*=', content) and re.search(r'Identities\s*=', content):
                    has_chain_hits = True
                    break
        except Exception as e:
            logging.warning(f"Error reading chain blast file {blast_file}: {str(e)}")
    
    # Check content of domain blast files for hits
    has_domain_hits = False
    for blast_file in domain_blast_files:
        try:
            with open(blast_file, 'r') as f:
                content = f.read()
                if re.search(r'Score\s*=', content) and re.search(r'Identities\s*=', content):
                    has_domain_hits = True
                    break
        except Exception as e:
            logging.warning(f"Error reading domain blast file {blast_file}: {str(e)}")
    
    # Return True if either chain or domain blast has hits
    return has_chain_hits or has_domain_hits

def create_no_hits_xml(pdb_id: str, chain_id: str, is_peptide: bool = False, ref_version: str = "develop") -> str:
    """
    Create XML content for a no-hits domain summary
    
    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        is_peptide: Whether this is a peptide (very short protein)
        ref_version: Reference version to use in file naming
        
    Returns:
        XML content as string
    """
    root = ET.Element("blast_summ_doc")
    blast_summ = ET.SubElement(root, "blast_summ")
    
    blast_summ.set("pdb", pdb_id)
    blast_summ.set("chain", chain_id)
    blast_summ.set("no_selfcomp", "true")
    blast_summ.set("chain_blast_no_hits", "true")
    blast_summ.set("domain_blast_no_hits", "true")
    
    if is_peptide:
        blast_summ.set("is_peptide", "true")
    
    # Format the XML with proper indentation
    rough_xml = ET.tostring(root, encoding='utf-8')
    
    # Manual pretty printing (simple version)
    xml_str = '<?xml version=\'1.0\' encoding=\'utf-8\'?>\n' + rough_xml.decode('utf-8')
    
    return xml_str

def write_xml_file(file_path: str, content: str, base_path: str) -> bool:
    """
    Write XML content to file
    
    Args:
        file_path: Relative or absolute path where to write the file
        content: XML content to write
        base_path: Base path to resolve relative paths
        
    Returns:
        True if file was written successfully
    """
    # Handle absolute and relative paths
    if os.path.isabs(file_path):
        full_path = file_path
    else:
        full_path = os.path.join(base_path, file_path)
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(full_path), exist_ok=True)
    
    try:
        with open(full_path, 'w') as f:
            f.write(content)
        return True
    except Exception as e:
        logging.error(f"Error writing file {full_path}: {str(e)}")
        return False

def update_database_status(context, protein_info: Dict[str, Any], file_exists: bool = True, 
                         is_complete: bool = True) -> bool:
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
    file_id = protein_info['file_id']
    process_id = protein_info['process_id']
    
    try:
        # Start a transaction
        cursor = context.db.connection.cursor()
        
        # Update the process_file record
        if file_exists:
            full_path = os.path.join(protein_info.get('base_path', ''), protein_info['file_path'])
            if os.path.exists(full_path):
                file_size = os.path.getsize(full_path)
            else:
                file_size = 0
                
            update_file_query = """
            UPDATE ecod_schema.process_file
            SET file_exists = TRUE, file_size = %s, last_checked = NOW()
            WHERE id = %s
            """
            cursor.execute(update_file_query, (file_size, file_id))
        else:
            update_file_query = """
            UPDATE ecod_schema.process_file
            SET file_exists = FALSE, file_size = 0, last_checked = NOW()
            WHERE id = %s
            """
            cursor.execute(update_file_query, (file_id,))
        
        # Update the process_status record based on whether it's complete
        if is_complete:
            # For no-hits or peptides that are properly processed
            status_query = """
            UPDATE ecod_schema.process_status
            SET status = 'success', current_stage = 'domain_partition_complete'
            WHERE id = %s
            """
        else:
            # For cases that still need processing
            status_query = """
            UPDATE ecod_schema.process_status
            SET status = 'pending', current_stage = 'domain_summary'
            WHERE id = %s
            """
        
        cursor.execute(status_query, (process_id,))
        
        # Commit the transaction
        context.db.connection.commit()
        cursor.close()
        return True
        
    except Exception as e:
        logging.error(f"Error updating database for {protein_info['pdb_id']}_{protein_info['chain_id']}: {str(e)}")
        if 'cursor' in locals():
            cursor.close()
        return False

def update_batch_completion(context, batch_id: int) -> bool:
    """
    Update the batch completion count and status
    
    Args:
        context: Application context
        batch_id: Batch ID to update
        
    Returns:
        True if batch was updated successfully
    """
    try:
        # Get the current count of completed items
        count_query = """
        SELECT COUNT(*)
        FROM ecod_schema.process_status
        WHERE batch_id = %s AND status = 'success' AND current_stage = 'domain_partition_complete'
        """
        
        count_result = context.db.execute_query(count_query, (batch_id,))
        if not count_result:
            return False
        
        completed_count = count_result[0][0]
        
        # Get total items for the batch
        batch_query = """
        SELECT total_items FROM ecod_schema.batch WHERE id = %s
        """
        
        batch_result = context.db.execute_query(batch_query, (batch_id,))
        if not batch_result:
            return False
        
        total_items = batch_result[0][0]
        
        # Update batch status
        update_query = """
        UPDATE ecod_schema.batch
        SET completed_items = %s,
            status = CASE 
                WHEN %s >= total_items THEN 'completed' 
                ELSE 'processing' 
            END
        WHERE id = %s
        """
        
        cursor = context.db.connection.cursor()
        cursor.execute(update_query, (completed_count, completed_count, batch_id))
        context.db.connection.commit()
        cursor.close()
        
        logging.info(f"Updated batch status: {completed_count}/{total_items} completed")
        return True
        
    except Exception as e:
        logging.error(f"Error updating batch status: {str(e)}")
        return False

def process_missing_files(context, batch_id: int, dry_run: bool = False,
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
        
        # Determine if this is a peptide or no-hits case
        is_peptide = length < 25  # Adjust threshold as needed
        
        if is_peptide:
            peptides.append(protein_info)
            protein_info['category'] = 'peptide'
            protein_info['has_hits'] = False
        else:
            # Check for blast hits
            has_blast_hits = check_blast_hits_in_filesystem(base_path, pdb_id, chain_id)
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
    
    # Initialize application context with the config file path
    # This assumes ApplicationContext handles loading and merging config files internally
    context = ApplicationContext(args.config)
    
    total, created, updated = process_missing_files(
        context, 
        args.batch_id, 
        args.dry_run,
        args.reference_version
    )
    
    # Exit with success if changes were made or it was a dry run
    return 0

if __name__ == "__main__":
    sys.exit(main())