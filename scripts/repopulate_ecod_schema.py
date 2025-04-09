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
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

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

def extract_id_from_query_result(result):
    """
    Extract ID from query result, handling different return formats
    
    Args:
        result: Query result in various formats
        
    Returns:
        ID value or 0 if not found
    """
    if not result:
        return 0
        
    # Handle list of dictionaries
    if isinstance(result, list) and result:
        if isinstance(result[0], dict) and 'id' in result[0]:
            return result[0]['id']
        # Handle list of tuples
        elif isinstance(result[0], tuple) and len(result[0]) > 0:
            return result[0][0]
    
    # Handle direct object with id attribute
    elif hasattr(result, 'id'):
        return result.id
    
    # Handle direct tuple
    elif isinstance(result, tuple) and len(result) > 0:
        return result[0]
    
    return 0

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
        fasta_dir = os.path.join(batch_path, 'fastas/batch_0')
        if os.path.exists(fasta_dir):
            fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta') or f.endswith('.fa')]
            batch_info['total_items'] = len(fasta_files)
        
        batches.append(batch_info)
    
    return batches

def discover_proteins_in_batch(batch_path: str) -> List[Dict[str, Any]]:
    """
    Discover proteins in a batch and their associated files
    
    Args:
        batch_path: Path to batch directory
        
    Returns:
        List of protein information dictionaries
    """
    proteins = []
    
    # Check for FASTA directory
    fasta_dir = os.path.join(batch_path, 'fastas')
    if not os.path.exists(fasta_dir):
        logging.warning(f"FASTA directory not found: {fasta_dir}")
        return proteins
    
    # Scan FASTA files to find proteins
    fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    
    for fasta_file in fasta_files:
        fasta_path = os.path.join(fasta_dir, fasta_file)
        source_id, pdb_id, chain_id = parse_source_id(fasta_file)
        
        # Initialize protein info
        protein_info = {
            'source_id': source_id,
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'fasta_path': fasta_path,
            'length': parse_fasta_length(fasta_path),
            'files': {
                'fasta': fasta_path
            }
        }
        
        # Look for associated files
        
        # BLAST results
        blast_chain_dir = os.path.join(batch_path, 'blast', 'chain')
        if os.path.exists(blast_chain_dir):
            blast_patterns = [
                f"{source_id}.blast",
                f"{source_id}_blast.txt",
                f"{source_id}.xml"
            ]
            
            for pattern in blast_patterns:
                blast_path = os.path.join(blast_chain_dir, pattern)
                if os.path.exists(blast_path):
                    protein_info['files']['blast_result'] = blast_path
                    break
        
        # Domain BLAST results
        blast_domain_dir = os.path.join(batch_path, 'blast', 'domain')
        if os.path.exists(blast_domain_dir):
            domain_blast_files = []
            for file in os.listdir(blast_domain_dir):
                if file.startswith(source_id) and (file.endswith('.blast') or file.endswith('_blast.txt')):
                    domain_blast_files.append(os.path.join(blast_domain_dir, file))
            
            if domain_blast_files:
                protein_info['files']['domain_blast_result'] = domain_blast_files[0]  # Take the first one
        
        # HHSearch results
        hhsearch_dir = os.path.join(batch_path, 'hhsearch')
        if os.path.exists(hhsearch_dir):
            # HHSearch profile
            profile_patterns = [
                f"{source_id}.hhm",
                f"{pdb_id}_{chain_id}.hhm"
            ]
            
            for pattern in profile_patterns:
                profile_path = os.path.join(hhsearch_dir, pattern)
                if os.path.exists(profile_path):
                    protein_info['files']['hhblits_profile'] = profile_path
                    break
            
            # HHSearch results
            result_patterns = [
                f"{source_id}.hhr",
                f"{pdb_id}_{chain_id}.hhr"
            ]
            
            for pattern in result_patterns:
                result_path = os.path.join(hhsearch_dir, pattern)
                if os.path.exists(result_path):
                    protein_info['files']['hhsearch_result'] = result_path
                    break
        
        # Domain partition/summary files
        domains_dir = os.path.join(batch_path, 'domains')
        if os.path.exists(domains_dir):
            summary_patterns = [
                f"{source_id}_domains.txt",
                f"{source_id}_domain_summary.json",
                f"{source_id}_domain_summary.txt",
                f"{pdb_id}_{chain_id}_domains.txt",
                f"{pdb_id}_{chain_id}_domain_summary.json"
            ]
            
            for pattern in summary_patterns:
                summary_path = os.path.join(domains_dir, pattern)
                if os.path.exists(summary_path):
                    protein_info['files']['domain_summary'] = summary_path
                    break
        
        proteins.append(protein_info)
    
    return proteins

def insert_batch(context: Any, batch_info: Dict[str, Any]) -> int:
    """
    Insert batch into database
    
    Args:
        context: Application context
        batch_info: Batch information dictionary
        
    Returns:
        Batch ID
    """
    query = """
    INSERT INTO ecod_schema.batch
    (batch_name, base_path, type, ref_version, total_items, status)
    VALUES (%s, %s, %s, %s, %s, %s)
    RETURNING id
    """
    
    result = context.db.execute_query(
        query, 
        (
            batch_info['name'],
            batch_info['path'],
            batch_info['type'],
            batch_info['ref_version'],
            batch_info['total_items'],
            'created'
        )
    )
    
    return extract_id_from_query_result(result)

def insert_protein(context: Any, protein_info: Dict[str, Any]) -> int:
    """
    Insert protein into database
    
    Args:
        context: Application context
        protein_info: Protein information dictionary
        
    Returns:
        Protein ID
    """
    query = """
    INSERT INTO ecod_schema.protein
    (pdb_id, chain_id, source_id, length)
    VALUES (%s, %s, %s, %s)
    RETURNING id
    """
    
    result = context.db.execute_query(
        query, 
        (
            protein_info['pdb_id'],
            protein_info['chain_id'],
            protein_info['source_id'],
            protein_info['length']
        )
    )
    
    return extract_id_from_query_result(result)

def insert_protein_sequence(context: Any, protein_id: int, fasta_path: str) -> int:
    """
    Insert protein sequence into database
    
    Args:
        context: Application context
        protein_id: Protein ID
        fasta_path: Path to FASTA file
        
    Returns:
        Sequence ID
    """
    sequence = get_fasta_sequence(fasta_path)
    if not sequence:
        return 0
    
    md5_hash = hashlib.md5(sequence.encode()).hexdigest()
    
    query = """
    INSERT INTO ecod_schema.protein_sequence
    (protein_id, sequence, md5_hash)
    VALUES (%s, %s, %s)
    RETURNING id
    """
    
    result = context.db.execute_query(query, (protein_id, sequence, md5_hash))
    
    return extract_id_from_query_result(result)

def insert_process_status(context: Any, protein_id: int, batch_id: int) -> int:
    """
    Insert process status into database
    
    Args:
        context: Application context
        protein_id: Protein ID
        batch_id: Batch ID
        
    Returns:
        Process ID
    """
    query = """
    INSERT INTO ecod_schema.process_status
    (protein_id, batch_id, current_stage, status)
    VALUES (%s, %s, %s, %s)
    RETURNING id
    """
    
    result = context.db.execute_query(
        query, 
        (protein_id, batch_id, 'initial', 'pending')
    )
    
    return extract_id_from_query_result(result)

def insert_process_file(context: Any, process_id: int, file_type: str, file_path: str) -> int:
    """
    Insert process file into database
    
    Args:
        context: Application context
        process_id: Process ID
        file_type: File type
        file_path: File path
        
    Returns:
        File ID
    """
    # Check if file exists
    file_exists = os.path.exists(file_path)
    file_size = os.path.getsize(file_path) if file_exists else 0
    
    query = """
    INSERT INTO ecod_schema.process_file
    (process_id, file_type, file_path, file_exists, file_size, last_checked)
    VALUES (%s, %s, %s, %s, %s, %s)
    RETURNING id
    """
    
    result = context.db.execute_query(
        query, 
        (process_id, file_type, file_path, file_exists, file_size, datetime.now())
    )
    
    return extract_id_from_query_result(result)

def update_process_status(context: Any, process_id: int, batch_id: int, protein_info: Dict[str, Any]) -> None:
    """
    Update process status based on available files
    
    Args:
        context: Application context
        process_id: Process ID
        batch_id: Batch ID
        protein_info: Protein information dictionary
    """
    # Determine current stage and status based on available files
    files = protein_info['files']
    
    if 'domain_summary' in files:
        current_stage = 'domain_summary'
        status = 'completed'
    elif 'hhsearch_result' in files:
        current_stage = 'hhsearch'
        status = 'completed'
    elif 'hhblits_profile' in files:
        current_stage = 'hhblits'
        status = 'completed'
    elif 'blast_result' in files:
        current_stage = 'blast'
        status = 'completed'
    else:
        current_stage = 'initial'
        status = 'pending'
    
    query = """
    UPDATE ecod_schema.process_status
    SET current_stage = %s, status = %s, updated_at = NOW()
    WHERE id = %s
    """
    
    context.db.execute(query, (current_stage, status, process_id))

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
        
        # Insert batch
        if not args.dry_run:
            batch_id = insert_batch(context, batch_info)
            if not batch_id:
                logger.error(f"Failed to insert batch {batch_info['name']}")
                continue
        else:
            # Dummy batch ID for dry run
            batch_id = -1
        
        # Discover proteins in batch
        proteins = discover_proteins_in_batch(batch_info['path'])
        logger.info(f"Found {len(proteins)} proteins in batch {batch_info['name']}")
        
        # Process each protein
        for protein_info in proteins:
            source_id = protein_info['source_id']
            logger.info(f"Processing protein: {source_id}")
            
            if not args.dry_run:
                # Insert protein
                protein_id = insert_protein(context, protein_info)
                if not protein_id:
                    logger.error(f"Failed to insert protein {source_id}")
                    continue
                
                # Insert protein sequence
                sequence_id = insert_protein_sequence(context, protein_id, protein_info['fasta_path'])
                if not sequence_id:
                    logger.warning(f"Failed to insert sequence for protein {source_id}")
                
                # Insert process status
                process_id = insert_process_status(context, protein_id, batch_id)
                if not process_id:
                    logger.error(f"Failed to insert process status for protein {source_id}")
                    continue
                
                # Insert process files
                for file_type, file_path in protein_info['files'].items():
                    file_id = insert_process_file(context, process_id, file_type, file_path)
                    if not file_id:
                        logger.warning(f"Failed to insert process file {file_type} for protein {source_id}")
                
                # Update process status based on available files
                update_process_status(context, process_id, batch_id, protein_info)
                
                total_proteins += 1
                total_files += len(protein_info['files'])
            else:
                # For dry run, just log what would be done
                logger.info(f"  Would insert protein: {source_id} with {len(protein_info['files'])} files")
                for file_type, file_path in protein_info['files'].items():
                    logger.debug(f"    {file_type}: {file_path}")
        
        total_batches += 1
        
        # If not dry run, update batch completion status
        if not args.dry_run:
            update_query = """
            UPDATE ecod_schema.batch
            SET completed_items = %s, status = 'indexed'
            WHERE id = %s
            """
            context.db.execute(update_query, (len(proteins), batch_id))
    
    logger.info(f"Processing complete: {total_batches} batches, {total_proteins} proteins, {total_files} files")
    
    if args.dry_run:
        logger.info("This was a dry run. No changes were made to the database.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())