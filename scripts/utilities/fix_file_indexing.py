#!/usr/bin/env python3
"""
fix_file_indexing.py - Update process_file records to reflect the correct file paths
in the new ECOD file organization system
"""

import os
import sys
import json
import logging
import argparse
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

def get_batch_path(base_path: str, batch_name: str) -> str:
    """Get the full path to a batch directory"""
    return os.path.join(base_path, 'batches', batch_name)

def find_files_in_batch(batch_path: str, source_id: str) -> Dict[str, str]:
    """
    Find files for a specific protein in the batch directory structure
    
    Args:
        batch_path: Path to the batch directory
        source_id: Protein source ID (typically PDB_ID_CHAIN_ID format)
        
    Returns:
        Dictionary mapping file types to file paths
    """
    # Parse source_id to get pdb_id and chain_id
    parts = source_id.split('_')
    if len(parts) >= 2:
        pdb_id = parts[0].lower()
        chain_id = parts[1]
    else:
        pdb_id = source_id
        chain_id = ""
    
    # Initialize result dictionary
    found_files = {}
    
    # Check for FASTA file
    fasta_dir = os.path.join(batch_path, 'fastas')
    if os.path.exists(fasta_dir):
        fasta_patterns = [
            f"{source_id}.fasta",
            f"{source_id}.fa",
            f"{pdb_id}_{chain_id}.fasta",
            f"{pdb_id}_{chain_id}.fa"
        ]
        
        for pattern in fasta_patterns:
            fasta_path = os.path.join(fasta_dir, pattern)
            if os.path.exists(fasta_path):
                found_files['fasta'] = fasta_path
                break
    
    # Check for BLAST result files
    blast_chain_dir = os.path.join(batch_path, 'blast', 'chain')
    if os.path.exists(blast_chain_dir):
        blast_patterns = [
            f"{source_id}.blast",
            f"{source_id}_blast.txt",
            f"{source_id}.xml",
            f"{pdb_id}_{chain_id}.blast",
            f"{pdb_id}_{chain_id}_blast.txt"
        ]
        
        for pattern in blast_patterns:
            blast_path = os.path.join(blast_chain_dir, pattern)
            if os.path.exists(blast_path):
                found_files['blast_result'] = blast_path
                break
    
    # Check for domain BLAST results
    blast_domain_dir = os.path.join(batch_path, 'blast', 'domain')
    if os.path.exists(blast_domain_dir):
        # These could be multiple files for different domains
        domain_blast_files = []
        for file in os.listdir(blast_domain_dir):
            if file.startswith(source_id) and (file.endswith('.blast') or file.endswith('_blast.txt')):
                domain_blast_files.append(os.path.join(blast_domain_dir, file))
        
        if domain_blast_files:
            found_files['domain_blast_result'] = domain_blast_files[0]  # Take the first one
    
    # Check for HHSearch results
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
                found_files['hhblits_profile'] = profile_path
                break
        
        # HHSearch results
        result_patterns = [
            f"{source_id}.hhr",
            f"{pdb_id}_{chain_id}.hhr"
        ]
        
        for pattern in result_patterns:
            result_path = os.path.join(hhsearch_dir, pattern)
            if os.path.exists(result_path):
                found_files['hhsearch_result'] = result_path
                break
    
    # Check for domain partition/summary files
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
                found_files['domain_summary'] = summary_path
                break
    
    return found_files

def update_process_file_records(
    context: Any, 
    process_id: int, 
    file_paths: Dict[str, str],
    dry_run: bool = False
) -> Tuple[int, int]:
    """
    Update process_file records with correct file paths
    
    Args:
        context: Application context with DB connection
        process_id: Process ID to update
        file_paths: Dictionary of file types to file paths
        dry_run: If True, don't actually update database
        
    Returns:
        Tuple of (updated_count, created_count)
    """
    updated_count = 0
    created_count = 0
    logger = logging.getLogger("ecod.fix_file_indexing")
    
    # Get existing records
    query = """
    SELECT id, file_type, file_path 
    FROM ecod_schema.process_file 
    WHERE process_id = %s
    """
    existing_records = context.db.execute_query(query, (process_id,))
    
    # Create a lookup by file type
    existing_by_type = {}
    for record in existing_records:
        file_type = record['file_type']
        existing_by_type[file_type] = record
    
    # Update or create records for each file type
    for file_type, file_path in file_paths.items():
        if file_type in existing_by_type:
            # Update existing record
            record = existing_by_type[file_type]
            record_id = record['id']
            old_path = record['file_path']
            
            logger.info(f"Updating {file_type} record: {record_id}")
            logger.info(f"  Old path: {old_path}")
            logger.info(f"  New path: {file_path}")
            
            if not dry_run:
                update_query = """
                UPDATE ecod_schema.process_file 
                SET file_path = %s, file_exists = %s, last_checked = NOW()
                WHERE id = %s
                """
                context.db.execute(update_query, (file_path, True, record_id))
            
            updated_count += 1
        else:
            # Create new record
            logger.info(f"Creating new {file_type} record")
            logger.info(f"  Path: {file_path}")
            
            if not dry_run:
                insert_query = """
                INSERT INTO ecod_schema.process_file 
                (process_id, file_type, file_path, file_exists, last_checked)
                VALUES (%s, %s, %s, %s, NOW())
                """
                context.db.execute(insert_query, (process_id, file_type, file_path, True))
            
            created_count += 1
    
    return updated_count, created_count

def get_batch_name(context: Any, batch_id: int) -> Optional[str]:
    """Get batch name from batch ID"""
    query = "SELECT batch_name FROM ecod_schema.batch WHERE id = %s"
    result = context.db.execute_query(query, (batch_id,))
    
    if result:
        return result[0][0]
    return None

def get_proteins_in_batch(context: Any, batch_id: int, protein_id: Optional[int] = None) -> List[Dict[str, Any]]:
    """Get proteins in a batch, optionally filtered by protein ID"""
    if protein_id:
        query = """
        SELECT p.id, p.source_id, ps.id as process_id
        FROM ecod_schema.protein p
        JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
        WHERE p.id = %s AND ps.batch_id = %s
        """
        return context.db.execute_query(query, (protein_id, batch_id))
    else:
        query = """
        SELECT p.id, p.source_id, ps.id as process_id
        FROM ecod_schema.protein p
        JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
        WHERE ps.batch_id = %s
        """
        return context.db.execute_query(query, (batch_id,))

def main():
    """Main function to fix file indexing"""
    parser = argparse.ArgumentParser(description='Fix file indexing in process_file table')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--protein-id', type=int,
                      help='Specific protein ID to fix (omit to fix all in batch)')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show what would be updated without making changes')
    parser.add_argument('--base-path', type=str, default='/data/ecod/pdb_updates',
                      help='Base path for ECOD data structure')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.fix_file_indexing")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get batch name
    batch_name = get_batch_name(context, args.batch_id)
    if not batch_name:
        logger.error(f"Batch {args.batch_id} not found")
        return 1
    
    logger.info(f"Processing batch {args.batch_id} ({batch_name})")
    
    # Get batch directory path
    batch_path = get_batch_path(args.base_path, batch_name)
    if not os.path.exists(batch_path):
        logger.error(f"Batch directory not found: {batch_path}")
        return 1
    
    logger.info(f"Using batch path: {batch_path}")
    
    # Check for metadata
    metadata_path = os.path.join(batch_path, 'metadata.json')
    if os.path.exists(metadata_path):
        try:
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
            logger.info(f"Found batch metadata: {metadata.get('description', 'No description')}")
        except Exception as e:
            logger.warning(f"Could not read metadata file: {str(e)}")
    
    # Get proteins to process
    proteins = get_proteins_in_batch(context, args.batch_id, args.protein_id)
    
    if not proteins:
        if args.protein_id:
            logger.error(f"Protein {args.protein_id} not found in batch {args.batch_id}")
        else:
            logger.error(f"No proteins found in batch {args.batch_id}")
        return 1
    
    logger.info(f"Found {len(proteins)} proteins to process")
    
    total_updated = 0
    total_created = 0
    
    # Process each protein
    for protein in proteins:
        protein_id = protein['id']
        source_id = protein['source_id']
        process_id = protein['process_id']
        
        logger.info(f"Processing protein {protein_id} ({source_id})")
        
        # Find files in the batch directory
        file_paths = find_files_in_batch(batch_path, source_id)
        
        # Count found files
        file_count = len(file_paths)
        logger.info(f"Found {file_count} files for protein {source_id}")
        
        if file_count == 0:
            logger.warning(f"No files found for protein {source_id}")
            continue
        
        # Log found files
        for file_type, file_path in file_paths.items():
            logger.debug(f"  {file_type}: {file_path}")
        
        # Update process_file records
        updated, created = update_process_file_records(
            context, process_id, file_paths, args.dry_run
        )
        
        total_updated += updated
        total_created += created
        
        logger.info(f"Updated {updated} records, created {created} records for protein {source_id}")
    
    logger.info(f"Finished processing. Total: {total_updated} updated, {total_created} created")
    
    if args.dry_run:
        logger.info("This was a dry run. No changes were made to the database.")
    else:
        # If we're successfully fixing file indexing, let's try to update the process status
        try:
            # Update process status for proteins that now have properly indexed files
            update_query = """
            UPDATE ecod_schema.process_status ps
            SET status = 'ready'
            FROM ecod_schema.process_file pf
            WHERE ps.id = pf.process_id
            AND ps.batch_id = %s
            AND ps.status = 'pending'
            AND pf.file_exists = true
            """
            if args.protein_id:
                update_query += " AND ps.protein_id = %s"
                context.db.execute(update_query, (args.batch_id, args.protein_id))
            else:
                context.db.execute(update_query, (args.batch_id,))
                
            logger.info("Updated process status for proteins with indexed files")
        except Exception as e:
            logger.error(f"Error updating process status: {str(e)}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())