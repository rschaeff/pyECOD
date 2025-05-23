#!/usr/bin/env python3
"""
repopulate_ecod_schema.py - Repopulate the ECOD schema tables based on files in the new directory structure

This script scans the file system and rebuilds the database tables from scratch.
"""

import os
import sys
import json
import re
import hashlib
import logging
import argparse
from datetime import datetime
from typing import Dict, Any, Optional, List, Tuple

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

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

def calculate_md5(file_path: str) -> str:
    """Calculate MD5 hash of file content"""
    md5_hash = hashlib.md5()
    
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
            
    return md5_hash.hexdigest()

def parse_fasta_length(fasta_path: str) -> int:
    """Parse protein length from FASTA file"""
    try:
        with open(fasta_path, 'r') as f:
            # Skip header line
            f.readline()
            sequence = ''
            for line in f:
                if not line.startswith('>'):
                    sequence += line.strip()
        return len(sequence)
    except Exception as e:
        logging.error(f"Error parsing FASTA file {fasta_path}: {str(e)}")
        return 0

def get_fasta_sequence(fasta_path: str) -> str:
    """Get protein sequence from FASTA file"""
    try:
        with open(fasta_path, 'r') as f:
            # Skip header line
            f.readline()
            sequence = ''
            for line in f:
                if not line.startswith('>'):
                    sequence += line.strip()
        return sequence
    except Exception as e:
        logging.error(f"Error reading FASTA file {fasta_path}: {str(e)}")
        return ""

def parse_source_id(filename: str) -> Tuple[str, str, str]:
    """
    Parse source_id from filename
    
    Returns:
        Tuple of (source_id, pdb_id, chain_id)
    """
    # Extract base filename without extension
    base = os.path.splitext(os.path.basename(filename))[0]
    
    # Remove common suffixes
    base = re.sub(r'_domains$|_domain_summary$|_blast$', '', base)
    
    # Split by underscore
    parts = base.split('_')
    
    if len(parts) >= 2:
        pdb_id = parts[0].lower()
        chain_id = parts[1]
        return f"{pdb_id}_{chain_id}", pdb_id, chain_id
    else:
        # If only one part, try to extract 4-letter PDB ID and rest as chain
        if len(base) >= 5:
            pdb_id = base[:4].lower()
            chain_id = base[4:]
            return f"{pdb_id}_{chain_id}", pdb_id, chain_id
        else:
            return base, base, ""

def discover_batches(base_path: str) -> List[Dict[str, Any]]:
    """
    Discover batch directories and their metadata
    
    Args:
        base_path: Base path to search (/data/ecod/pdb_updates)
        
    Returns:
        List of batch information dictionaries
    """
    batches = []
    batches_dir = os.path.join(base_path, 'batches')
    
    if not os.path.exists(batches_dir):
        logging.error(f"Batches directory not found: {batches_dir}")
        return batches
    
    for batch_name in os.listdir(batches_dir):
        batch_path = os.path.join(batches_dir, batch_name)
        
        if not os.path.isdir(batch_path):
            continue
        
        # Initialize batch info
        batch_info = {
            'name': batch_name,
            'path': batch_path,
            'type': 'pdb',  # default
            'ref_version': 'unknown',
            'total_items': 0,
            'metadata': {}
        }
        
        # Check for metadata file
        metadata_path = os.path.join(batch_path, 'metadata.json')
        if os.path.exists(metadata_path):
            try:
                with open(metadata_path, 'r') as f:
                    metadata = json.load(f)
                batch_info['metadata'] = metadata
                
                if 'type' in metadata:
                    batch_info['type'] = metadata['type']
                    
                if 'ref_version' in metadata:
                    batch_info['ref_version'] = metadata['ref_version']
            except Exception as e:
                logging.warning(f"Error reading metadata file {metadata_path}: {str(e)}")
        
        # Count proteins in the batch
        fasta_dir = os.path.join(batch_path, 'fastas')
        if os.path.exists(fasta_dir):
            count = count_fasta_files(fasta_dir)
            batch_info['total_items'] = count
        
        batches.append(batch_info)
    
    return batches

def count_fasta_files(directory: str) -> int:
    """Count all FASTA files in a directory and its subdirectories"""
    count = 0
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.fasta') or file.endswith('.fa'):
                count += 1
    return count

def discover_proteins_in_batch(batch_path: str) -> List[Dict[str, Any]]:
    """
    Discover proteins in a batch and their associated files,
    recursively searching sub-batch directories
    
    Args:
        batch_path: Path to batch directory
        
    Returns:
        List of protein information dictionaries
    """
    proteins = {}  # Use a dictionary to prevent duplicates by source_id
    logger = logging.getLogger("ecod.repopulate_schema")
    
    # Check for FASTA directory
    fasta_dir = os.path.join(batch_path, 'fastas')
    if not os.path.exists(fasta_dir):
        logger.warning(f"FASTA directory not found: {fasta_dir}")
        return []
    
    # Function to scan FASTA files and subdirectories
    def scan_fasta_directory(directory):
        for item in os.listdir(directory):
            item_path = os.path.join(directory, item)
            
            # If it's a directory, scan it recursively
            if os.path.isdir(item_path):
                scan_fasta_directory(item_path)
            # If it's a FASTA file, process it
            elif item.endswith('.fasta') or item.endswith('.fa'):
                source_id, pdb_id, chain_id = parse_source_id(item)
                
                # Skip if we've already found this protein
                if source_id in proteins:
                    continue
                
                # Initialize protein info
                proteins[source_id] = {
                    'source_id': source_id,
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'fasta_path': item_path,
                    'length': parse_fasta_length(item_path),
                    'files': {
                        'fasta': item_path
                    }
                }
    
    # Scan FASTA directory and subdirectories
    logger.info(f"Scanning FASTA directory: {fasta_dir}")
    scan_fasta_directory(fasta_dir)
    
    # Find associated files for each protein
    for source_id, protein_info in proteins.items():
        # Look for blast results in all subdirectories
        blast_chain_dir = os.path.join(batch_path, 'blast', 'chain')
        if os.path.exists(blast_chain_dir):
            # Function to recursively look for blast files
            def find_blast_file(directory, source_id):
                for root, dirs, files in os.walk(directory):
                    for file in files:
                        if (file.startswith(source_id) and 
                            (file.endswith('chainwise_blast.xml'))):
                            return os.path.join(root, file)
                return None
            
            blast_path = find_blast_file(blast_chain_dir, source_id)
            if blast_path:
                # Use chain_blast_result instead of blast_result
                protein_info['files']['chain_blast_result'] = blast_path
        
        # Look for domain blast results
        blast_domain_dir = os.path.join(batch_path, 'blast', 'domain')
        if os.path.exists(blast_domain_dir):
            def find_domain_blast_file(directory, source_id):
                for root, dirs, files in os.walk(directory):
                    for file in files:
                        if (file.startswith(source_id) and 
                            (file.endswith('domain_blast.xml'))):
                            return os.path.join(root, file)
                return None
            
            domain_blast_path = find_domain_blast_file(blast_domain_dir, source_id)
            if domain_blast_path:
                protein_info['files']['domain_blast_result'] = domain_blast_path
        
        # Look for HHSearch results
        hhsearch_dir = os.path.join(batch_path, 'hhsearch')
        if os.path.exists(hhsearch_dir):
            def find_hh_file(directory, patterns):
                for root, dirs, files in os.walk(directory):
                    for file in files:
                        for pattern in patterns:
                            if file == pattern:
                                return os.path.join(root, file)
                return None
            
            # HHSearch profile
            profile_patterns = [f"{source_id}.hhm", f"{protein_info['pdb_id']}_{protein_info['chain_id']}.hhm"]
            profile_path = find_hh_file(hhsearch_dir, profile_patterns)
            if profile_path:
                protein_info['files']['hhblits_profile'] = profile_path
            
            # HHSearch results
            result_patterns = [f"{source_id}.hhr", f"{protein_info['pdb_id']}_{protein_info['chain_id']}.hhr"]
            result_path = find_hh_file(hhsearch_dir, result_patterns)
            if result_path:
                protein_info['files']['hhsearch_result'] = result_path
        
        # Look for domain summary files
        domains_dir = os.path.join(batch_path, 'domains')
        if os.path.exists(domains_dir):
            def find_domain_summary(directory, patterns):
                for root, dirs, files in os.walk(directory):
                    for file in files:
                        for pattern in patterns:
                            if file == pattern:
                                return os.path.join(root, file)
                return None
            
            summary_patterns = [
                f"{source_id}_domains.txt",
                f"{source_id}_domain_summary.json",
                f"{source_id}_domain_summary.txt",
                f"{protein_info['pdb_id']}_{protein_info['chain_id']}_domains.txt",
                f"{protein_info['pdb_id']}_{protein_info['chain_id']}_domain_summary.json"
            ]
            
            summary_path = find_domain_summary(domains_dir, summary_patterns)
            if summary_path:
                protein_info['files']['domain_summary'] = summary_path
    
    return list(proteins.values())

def get_batch_id(context: Any, batch_name: str, batch_info: Dict[str, Any], dry_run: bool = False) -> int:
    """
    Get batch ID, creating the batch if it doesn't exist
    
    Args:
        context: Application context
        batch_name: Batch name
        batch_info: Batch information dictionary
        dry_run: If True, don't actually create batch
        
    Returns:
        Batch ID or 0 if batch doesn't exist and this is a dry run
    """
    query = "SELECT id FROM ecod_schema.batch WHERE batch_name = %s"
    result = context.db.execute_query(query, (batch_name,))
    
    if result and len(result) > 0:
        # Extract ID from result
        if isinstance(result[0], dict) and 'id' in result[0]:
            return result[0]['id']
        elif isinstance(result[0], tuple) and len(result[0]) > 0:
            return result[0][0]
    
    # Batch doesn't exist, create it
    if not dry_run:
        insert_query = """
        INSERT INTO ecod_schema.batch
        (batch_name, base_path, type, ref_version, total_items, status)
        VALUES (%s, %s, %s, %s, %s, %s)
        RETURNING id
        """
        
        insert_result = context.db.execute_query(
            insert_query, 
            (
                batch_info['name'],
                batch_info['path'],
                batch_info['type'],
                batch_info['ref_version'],
                batch_info['total_items'],
                'created'
            )
        )
        
        if insert_result and len(insert_result) > 0:
            # Extract ID from result
            if isinstance(insert_result[0], dict) and 'id' in insert_result[0]:
                return insert_result[0]['id']
            elif isinstance(insert_result[0], tuple) and len(insert_result[0]) > 0:
                return insert_result[0][0]
    
    return 0  # Return 0 for dry run or if insertion failed

def get_protein_id(context: Any, source_id: str, protein_info: Dict[str, Any], dry_run: bool = False) -> int:
    """
    Get protein ID, creating the protein if it doesn't exist
    
    Args:
        context: Application context
        source_id: Protein source ID
        protein_info: Protein information dictionary
        dry_run: If True, don't actually create protein
        
    Returns:
        Protein ID or 0 if protein doesn't exist and this is a dry run
    """
    query = "SELECT id FROM ecod_schema.protein WHERE source_id = %s"
    result = context.db.execute_query(query, (source_id,))
    
    if result and len(result) > 0:
        # Extract ID from result
        if isinstance(result[0], dict) and 'id' in result[0]:
            return result[0]['id']
        elif isinstance(result[0], tuple) and len(result[0]) > 0:
            return result[0][0]
    
    # Protein doesn't exist, create it
    if not dry_run:
        insert_query = """
        INSERT INTO ecod_schema.protein
        (pdb_id, chain_id, source_id, length)
        VALUES (%s, %s, %s, %s)
        RETURNING id
        """
        
        insert_result = context.db.execute_query(
            insert_query, 
            (
                protein_info['pdb_id'],
                protein_info['chain_id'],
                protein_info['source_id'],
                protein_info['length']
            )
        )
        
        if insert_result and len(insert_result) > 0:
            # Extract ID from result
            if isinstance(insert_result[0], dict) and 'id' in insert_result[0]:
                return insert_result[0]['id']
            elif isinstance(insert_result[0], tuple) and len(insert_result[0]) > 0:
                return insert_result[0][0]
                
            # If we get here but ID is not easily extractable, try to look it up
            lookup_result = context.db.execute_query(query, (source_id,))
            if lookup_result and len(lookup_result) > 0:
                if isinstance(lookup_result[0], dict) and 'id' in lookup_result[0]:
                    return lookup_result[0]['id']
                elif isinstance(lookup_result[0], tuple) and len(lookup_result[0]) > 0:
                    return lookup_result[0][0]
    
    return 0  # Return 0 for dry run or if insertion failed

def ensure_protein_sequence(context: Any, protein_id: int, fasta_path: str, dry_run: bool = False) -> int:
    """
    Ensure protein sequence exists, creating it if necessary
    
    Args:
        context: Application context
        protein_id: Protein ID
        fasta_path: Path to FASTA file
        dry_run: If True, don't actually create sequence
        
    Returns:
        Sequence ID or 0 if sequence doesn't exist and this is a dry run
    """
    # Check if sequence already exists for this protein
    query = "SELECT id FROM ecod_schema.protein_sequence WHERE protein_id = %s"
    result = context.db.execute_query(query, (protein_id,))
    
    if result and len(result) > 0:
        # Extract ID from result
        if isinstance(result[0], dict) and 'id' in result[0]:
            return result[0]['id']
        elif isinstance(result[0], tuple) and len(result[0]) > 0:
            return result[0][0]
    
    # Sequence doesn't exist, create it if not dry run
    if not dry_run:
        sequence = get_fasta_sequence(fasta_path)
        if not sequence:
            return 0
        
        md5_hash = hashlib.md5(sequence.encode()).hexdigest()
        
        insert_query = """
        INSERT INTO ecod_schema.protein_sequence
        (protein_id, sequence, md5_hash)
        VALUES (%s, %s, %s)
        RETURNING id
        """
        
        insert_result = context.db.execute_query(insert_query, (protein_id, sequence, md5_hash))
        
        if insert_result and len(insert_result) > 0:
            # Extract ID from result
            if isinstance(insert_result[0], dict) and 'id' in insert_result[0]:
                return insert_result[0]['id']
            elif isinstance(insert_result[0], tuple) and len(insert_result[0]) > 0:
                return insert_result[0][0]
    
    return 0  # Return 0 for dry run or if insertion failed

def get_process_id(context: Any, protein_id: int, batch_id: int, dry_run: bool = False) -> int:
    """
    Get process ID, creating the process if it doesn't exist
    
    Args:
        context: Application context
        protein_id: Protein ID
        batch_id: Batch ID
        dry_run: If True, don't actually create process
        
    Returns:
        Process ID or 0 if process doesn't exist and this is a dry run
    """
    query = """
    SELECT id FROM ecod_schema.process_status
    WHERE protein_id = %s AND batch_id = %s
    """
    result = context.db.execute_query(query, (protein_id, batch_id))
    
    if result and len(result) > 0:
        # Extract ID from result
        if isinstance(result[0], dict) and 'id' in result[0]:
            return result[0]['id']
        elif isinstance(result[0], tuple) and len(result[0]) > 0:
            return result[0][0]
    
    # Process doesn't exist, create it if not dry run
    if not dry_run:
        insert_query = """
        INSERT INTO ecod_schema.process_status
        (protein_id, batch_id, current_stage, status)
        VALUES (%s, %s, %s, %s)
        RETURNING id
        """
        
        insert_result = context.db.execute_query(
            insert_query, 
            (protein_id, batch_id, 'initial', 'pending')
        )
        
        if insert_result and len(insert_result) > 0:
            # Extract ID from result
            if isinstance(insert_result[0], dict) and 'id' in insert_result[0]:
                return insert_result[0]['id']
            elif isinstance(insert_result[0], tuple) and len(insert_result[0]) > 0:
                return insert_result[0][0]
    
    return 0  # Return 0 for dry run or if insertion failed

def ensure_process_file(context: Any, process_id: int, file_type: str, file_path: str, dry_run: bool = False) -> int:
    """
    Ensure process file exists, creating or updating it as necessary
    
    Args:
        context: Application context
        process_id: Process ID
        file_type: File type
        file_path: File path
        dry_run: If True, don't actually create or update file
        
    Returns:
        File ID or 0 if file doesn't exist and this is a dry run
    """
    # Check if file already exists
    query = """
    SELECT id FROM ecod_schema.process_file
    WHERE process_id = %s AND file_type = %s
    """
    result = context.db.execute_query(query, (process_id, file_type))
    
    file_exists = os.path.exists(file_path)
    file_size = os.path.getsize(file_path) if file_exists else 0
    
    if result and len(result) > 0:
        # File exists, update it
        file_id = result[0]['id'] if isinstance(result[0], dict) else result[0][0]
        
        if not dry_run:
            update_query = """
            UPDATE ecod_schema.process_file
            SET file_path = %s, file_exists = %s, file_size = %s, last_checked = NOW()
            WHERE id = %s
            """
            context.db.execute_query(update_query, (file_path, file_exists, file_size, file_id))
        
        return file_id
    
    # File doesn't exist, create it if not dry run
    if not dry_run:
        insert_query = """
        INSERT INTO ecod_schema.process_file
        (process_id, file_type, file_path, file_exists, file_size, last_checked)
        VALUES (%s, %s, %s, %s, %s, %s)
        RETURNING id
        """
        
        insert_result = context.db.execute_query(
            insert_query, 
            (process_id, file_type, file_path, file_exists, file_size, datetime.now())
        )
        
        if insert_result and len(insert_result) > 0:
            # Extract ID from result
            if isinstance(insert_result[0], dict) and 'id' in insert_result[0]:
                return insert_result[0]['id']
            elif isinstance(insert_result[0], tuple) and len(insert_result[0]) > 0:
                return insert_result[0][0]
    
    return 0  # Return 0 for dry run or if insertion failed

def update_process_status(context: Any, process_id: int, protein_info: Dict[str, Any], dry_run: bool = False) -> None:
    """
    Update process status based on available files
    
    Args:
        context: Application context
        process_id: Process ID
        protein_info: Protein information dictionary
        dry_run: If True, don't actually update status
    """
    # Determine current stage and status based on available files
    # Determine current stage and status based on available files
    files = protein_info['files']

    if 'domain_summary' in files:
        current_stage = 'domain_summary'
        status = 'success'  # Using 'success' to match pipeline expectations
    elif 'hhsearch_result' in files:
        current_stage = 'hhsearch'
        status = 'success'  # Using 'success' to match pipeline expectations
    elif 'hhblits_profile' in files:
        current_stage = 'hhblits'
        status = 'success'  # Using 'success' to match pipeline expectations
    elif 'chain_blast_result' in files and 'domain_blast_result' in files:
        current_stage = 'blast'
        status = 'success'  # Using 'success' to match pipeline expectations
    else:
        current_stage = 'initial'
        status = 'pending'
    
    if not dry_run:
        query = """
        UPDATE ecod_schema.process_status
        SET current_stage = %s, status = %s, updated_at = NOW()
        WHERE id = %s
        """
        
        context.db.execute_query(query, (current_stage, status, process_id))

def update_batch_status(context: Any, batch_id: int, total_proteins: int, dry_run: bool = False) -> None:
    """
    Update batch status
    
    Args:
        context: Application context
        batch_id: Batch ID
        total_proteins: Total number of proteins in batch
        dry_run: If True, don't actually update status
    """
    if not dry_run:
        query = """
        UPDATE ecod_schema.batch
        SET completed_items = %s, status = 'indexed'
        WHERE id = %s
        """
        
        context.db.execute_query(query, (total_proteins, batch_id))

def main():
    """Main function to repopulate ECOD schema"""
    parser = argparse.ArgumentParser(description='Repopulate ECOD schema tables from file system')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--base-path', type=str, default='/data/ecod/pdb_updates',
                      help='Base path for ECOD data structure')
    parser.add_argument('--batch-name', type=str,
                      help='Process only this specific batch (omit to process all)')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show what would be done without making changes')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.repopulate_schema")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Discover batches
    logger.info(f"Discovering batches in {args.base_path}")
    batches = discover_batches(args.base_path)
    
    if args.batch_name:
        batches = [b for b in batches if b['name'] == args.batch_name]
        if not batches:
            logger.error(f"Batch '{args.batch_name}' not found")
            return 1
    
    logger.info(f"Found {len(batches)} batches")
    
    total_batches = 0
    total_proteins = 0
    total_files = 0
    
    # Process each batch
    for batch_info in batches:
        logger.info(f"Processing batch: {batch_info['name']}")
        
        # Get or create batch
        batch_id = get_batch_id(context, batch_info['name'], batch_info, args.dry_run)
        if not batch_id and not args.dry_run:
            logger.error(f"Failed to get or create batch {batch_info['name']}")
            continue
        elif args.dry_run and not batch_id:
            # Use dummy ID for dry run
            batch_id = -1
            
        logger.info(f"Using batch ID: {batch_id}")
        
        # Discover proteins in batch
        proteins = discover_proteins_in_batch(batch_info['path'])
        logger.info(f"Found {len(proteins)} proteins in batch {batch_info['name']}")
        
        batch_proteins = 0
        batch_files = 0
        
        # Process each protein
        for protein_info in proteins:
            source_id = protein_info['source_id']
            logger.info(f"Processing protein: {source_id}")
            
            try:
                # Get or create protein
                protein_id = get_protein_id(context, source_id, protein_info, args.dry_run)
                if not protein_id and not args.dry_run:
                    logger.error(f"Failed to get or create protein {source_id}")
                    continue
                elif args.dry_run and not protein_id:
                    # Use dummy ID for dry run
                    protein_id = -1
                
                # Ensure protein sequence exists
                sequence_id = ensure_protein_sequence(context, protein_id, protein_info['fasta_path'], args.dry_run)
                if not sequence_id and not args.dry_run:
                    logger.warning(f"Failed to ensure sequence for protein {source_id}")
                
                # Get or create process status
                process_id = get_process_id(context, protein_id, batch_id, args.dry_run)
                if not process_id and not args.dry_run:
                    logger.error(f"Failed to get or create process status for protein {source_id}")
                    continue
                elif args.dry_run and not process_id:
                    # Use dummy ID for dry run
                    process_id = -1
                
                # Process files
                for file_type, file_path in protein_info['files'].items():
                    file_id = ensure_process_file(context, process_id, file_type, file_path, args.dry_run)
                    if not file_id and not args.dry_run:
                        logger.warning(f"Failed to ensure process file {file_type} for protein {source_id}")
                    
                # Update process status
                update_process_status(context, process_id, protein_info, args.dry_run)
                
                batch_proteins += 1
                batch_files += len(protein_info['files'])
                
            except Exception as e:
                logger.error(f"Error processing protein {source_id}: {str(e)}", exc_info=True)
                continue
        
        # Update batch status
        update_batch_status(context, batch_id, batch_proteins, args.dry_run)
        
        total_batches += 1
        total_proteins += batch_proteins
        total_files += batch_files
        
        logger.info(f"Completed batch {batch_info['name']}: {batch_proteins} proteins, {batch_files} files")
    
    logger.info(f"Processing complete: {total_batches} batches, {total_proteins} proteins, {total_files} files")
    
    if args.dry_run:
        logger.info("This was a dry run. No changes were made to the database.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())