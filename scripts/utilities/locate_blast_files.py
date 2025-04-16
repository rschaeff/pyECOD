#!/usr/bin/env python3
"""
locate_blast_files.py - Locate BLAST files for a specific protein
"""

import os
import sys
import argparse
import logging
import xml.etree.ElementTree as ET
from typing import List, Dict, Any, Optional
import subprocess

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

def search_for_blast_files(pdb_id: str, chain_id: str, batch_path: str) -> List[Dict[str, Any]]:
    """Search for BLAST files in batch directory"""
    logger = logging.getLogger("ecod.locate_blast")
    
    # Common patterns for blast files
    patterns = [
        f"{pdb_id}_{chain_id}.chainwise_blast.xml",
        f"{pdb_id}_{chain_id}.chain_blast.xml",
        f"{pdb_id}_{chain_id}.blast.xml",
        f"{pdb_id}_{chain_id}_blast.xml",
        f"{pdb_id}{chain_id}.blast.xml"
    ]
    
    # Common directories
    directories = [
        os.path.join(batch_path),
        os.path.join(batch_path, pdb_id, chain_id),
        os.path.join(batch_path, f"{pdb_id}_{chain_id}"),
        os.path.join(batch_path, "blast", "chain"),
        os.path.join(batch_path, "blast", "chain", "batch_0"),
        os.path.join(batch_path, "blast", "chain", "batch_1"),
        os.path.join(batch_path, "blast", "chain", "batch_2"),
        os.path.join(batch_path, "chain_blast"),
        os.path.join(batch_path, "chain_blast_results")
    ]
    
    found_files = []
    
    for directory in directories:
        if not os.path.exists(directory):
            continue
            
        logger.info(f"Checking directory: {directory}")
        
        try:
            files = os.listdir(directory)
        except PermissionError:
            logger.warning(f"Permission denied: {directory}")
            continue
            
        # Check for pattern matches
        for pattern in patterns:
            matching_files = [f for f in files if pattern.lower() in f.lower()]
            for file in matching_files:
                full_path = os.path.join(directory, file)
                size = os.path.getsize(full_path) if os.path.exists(full_path) else 0
                logger.info(f"FOUND: {full_path} (size: {size} bytes)")
                
                # Check if file has valid XML and BLAST hits
                valid, hit_count = check_blast_xml(full_path)
                
                found_files.append({
                    'file_path': full_path,
                    'relative_path': os.path.relpath(full_path, batch_path),
                    'size': size,
                    'is_valid_xml': valid,
                    'hit_count': hit_count
                })
        
        # Check for any XML files with PDB ID and chain ID in the name
        all_xml = [f for f in files if f.endswith('.xml') and 
                  pdb_id.lower() in f.lower() and 
                  chain_id.lower() in f.lower()]
        
        for file in all_xml:
            if not any(pattern.lower() in file.lower() for pattern in patterns):
                full_path = os.path.join(directory, file)
                size = os.path.getsize(full_path) if os.path.exists(full_path) else 0
                logger.info(f"OTHER XML: {full_path} (size: {size} bytes)")
                
                # Check if file has valid XML and BLAST hits
                valid, hit_count = check_blast_xml(full_path)
                
                found_files.append({
                    'file_path': full_path,
                    'relative_path': os.path.relpath(full_path, batch_path),
                    'size': size,
                    'is_valid_xml': valid,
                    'hit_count': hit_count
                })
    
    # Use find command for more thorough search
    try:
        cmd = f"find {batch_path} -name '*{pdb_id}*{chain_id}*blast*.xml' -type f 2>/dev/null"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.stdout:
            logger.info("Find command results:")
            for line in result.stdout.splitlines():
                if os.path.exists(line) and not any(f['file_path'] == line for f in found_files):
                    size = os.path.getsize(line)
                    logger.info(f"FOUND: {line} (size: {size} bytes)")
                    
                    # Check if file has valid XML and BLAST hits
                    valid, hit_count = check_blast_xml(line)
                    
                    found_files.append({
                        'file_path': line,
                        'relative_path': os.path.relpath(line, batch_path),
                        'size': size,
                        'is_valid_xml': valid,
                        'hit_count': hit_count
                    })
    except Exception as e:
        logger.error(f"Error running find command: {e}")
    
    return found_files

def check_blast_xml(file_path: str) -> tuple[bool, int]:
    """Check if file is valid XML and contains BLAST hits"""
    logger = logging.getLogger("ecod.locate_blast")
    
    if not os.path.exists(file_path):
        return False, 0
    
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Count hits in BLAST XML format
        hits = root.findall(".//Hit")
        hit_count = len(hits)
        
        logger.debug(f"File {file_path} is valid XML with {hit_count} hits")
        return True, hit_count
    except Exception as e:
        logger.debug(f"Error parsing XML file {file_path}: {e}")
        return False, 0

def check_db_file_paths(context: ApplicationContext, process_id: int) -> List[Dict[str, Any]]:
    """Check file paths in the database"""
    logger = logging.getLogger("ecod.locate_blast")
    
    # Get batch path
    batch_query = """
    SELECT b.base_path
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.batch b ON ps.batch_id = b.id
    WHERE ps.id = %s
    """
    
    batch_result = context.db.execute_query(batch_query, (process_id,))
    
    if not batch_result:
        logger.error(f"Could not find batch path for process {process_id}")
        return []
    
    batch_path = batch_result[0][0]
    
    # Get file records
    files_query = """
    SELECT id, file_type, file_path, file_exists, file_size
    FROM ecod_schema.process_file
    WHERE process_id = %s
    """
    
    files = context.db.execute_query(files_query, (process_id,))
    
    # Process file records
    file_records = []
    for file in files:
        file_id = file[0]
        file_type = file[1]
        file_path = file[2]
        file_exists_db = file[3]
        file_size_db = file[4]
        
        # Construct full path
        if os.path.isabs(file_path):
            full_path = file_path
        else:
            full_path = os.path.join(batch_path, file_path)
        
        # Check if file actually exists
        file_exists = os.path.exists(full_path)
        file_size = os.path.getsize(full_path) if file_exists else 0
        
        # Check for '..' in path that might cause normalization issues
        has_path_issue = '..' in file_path
        if has_path_issue:
            normalized_path = os.path.normpath(full_path)
            normalized_exists = os.path.exists(normalized_path)
        else:
            normalized_path = None
            normalized_exists = None
        
        file_records.append({
            'id': file_id,
            'file_type': file_type,
            'file_path': file_path,
            'full_path': full_path,
            'db_exists': file_exists_db,
            'actual_exists': file_exists,
            'db_size': file_size_db,
            'actual_size': file_size,
            'has_path_issue': has_path_issue,
            'normalized_path': normalized_path,
            'normalized_exists': normalized_exists
        })
    
    return file_records

def update_file_record(context: ApplicationContext, file_id: int, file_path: str, file_exists: bool, file_size: int) -> bool:
    """Update a file record in the database"""
    logger = logging.getLogger("ecod.locate_blast")
    
    # Update file record
    update_query = """
    UPDATE ecod_schema.process_file
    SET file_path = %s, file_exists = %s, file_size = %s, last_checked = NOW()
    WHERE id = %s
    """
    
    try:
        context.db.execute_query(update_query, (file_path, file_exists, file_size, file_id))
        logger.info(f"Updated file record {file_id} with path {file_path}")
        return True
    except Exception as e:
        logger.error(f"Error updating file record: {e}")
        return False

def create_file_record(context: ApplicationContext, process_id: int, file_type: str, file_path: str, file_exists: bool, file_size: int) -> bool:
    """Create a new file record in the database"""
    logger = logging.getLogger("ecod.locate_blast")
    
    # Create file record
    insert_query = """
    INSERT INTO ecod_schema.process_file
    (process_id, file_type, file_path, file_exists, file_size, last_checked)
    VALUES (%s, %s, %s, %s, %s, NOW())
    RETURNING id
    """
    
    try:
        result = context.db.execute_query(insert_query, (process_id, file_type, file_path, file_exists, file_size))
        file_id = result[0][0]
        logger.info(f"Created new file record {file_id} for {file_type} with path {file_path}")
        return True
    except Exception as e:
        logger.error(f"Error creating file record: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Locate BLAST files for a protein')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID')
    parser.add_argument('--protein-id', type=int, required=True,
                      help='Protein ID')
    parser.add_argument('--fix', action='store_true',
                      help='Fix database records if needed')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.locate_blast")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get protein and batch information
    query = """
    SELECT 
        p.id, p.source_id, p.pdb_id, p.chain_id,
        ps.id as process_id, b.base_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    WHERE 
        p.id = %s AND ps.batch_id = %s
    """
    
    result = context.db.execute_query(query, (args.protein_id, args.batch_id))
    
    if not result:
        logger.error(f"Protein {args.protein_id} not found in batch {args.batch_id}")
        return 1
    
    protein_info = result[0]
    pdb_id = protein_info[2]
    chain_id = protein_info[3]
    process_id = protein_info[4]
    batch_path = protein_info[5]
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    logger.info(f"Looking for BLAST files for protein: {pdb_chain} (ID: {args.protein_id})")
    logger.info(f"Process ID: {process_id}, Batch path: {batch_path}")
    
    # Check file records in database
    logger.info("Checking file records in database...")
    file_records = check_db_file_paths(context, process_id)
    
    if not file_records:
        logger.warning("No file records found in database")
    else:
        logger.info(f"Found {len(file_records)} file records in database")
        
        # Check each file record
        for record in file_records:
            logger.info(f"File record: {record['file_type']}")
            logger.info(f"  Path: {record['file_path']}")
            logger.info(f"  Full path: {record['full_path']}")
            logger.info(f"  DB says exists: {record['db_exists']}, Actually exists: {record['actual_exists']}")
            logger.info(f"  DB size: {record['db_size']}, Actual size: {record['actual_size']}")
            
            if record['has_path_issue']:
                logger.warning(f"  Path contains '..' which may cause resolution issues")
                logger.info(f"  Normalized path: {record['normalized_path']}")
                logger.info(f"  Normalized path exists: {record['normalized_exists']}")
    
    # Search for BLAST files
    logger.info("Searching for BLAST files...")
    found_files = search_for_blast_files(pdb_id, chain_id, batch_path)
    
    if not found_files:
        logger.warning("No BLAST files found")
    else:
        logger.info(f"Found {len(found_files)} potential BLAST files")
        
        for file in found_files:
            logger.info(f"Found file: {file['file_path']}")
            logger.info(f"  Relative path: {file['relative_path']}")
            logger.info(f"  Size: {file['size']} bytes")
            logger.info(f"  Valid XML: {file['is_valid_xml']}")
            logger.info(f"  Hit count: {file['hit_count']}")
    
    # Fix database records if needed
    if args.fix and found_files:
        # Find files with the most hits
        chain_blast_files = [f for f in found_files if 'chain' in f['file_path'].lower() and f['is_valid_xml'] and f['hit_count'] > 0]
        domain_blast_files = [f for f in found_files if 'domain' in f['file_path'].lower() and f['is_valid_xml'] and f['hit_count'] > 0]
        
        # Sort by hit count
        chain_blast_files.sort(key=lambda x: x['hit_count'], reverse=True)
        domain_blast_files.sort(key=lambda x: x['hit_count'], reverse=True)
        
        # Check if we need to update chain_blast_result record
        chain_blast_record = next((r for r in file_records if r['file_type'] == 'chain_blast_result'), None)
        
        if chain_blast_files:
            best_chain_file = chain_blast_files[0]
            
            if chain_blast_record:
                if not chain_blast_record['actual_exists'] or chain_blast_record['actual_size'] == 0:
                    # Update record with the found file
                    if args.fix:
                        logger.info(f"Updating chain_blast_result record with new path: {best_chain_file['relative_path']}")
                        update_file_record(
                            context, 
                            chain_blast_record['id'], 
                            best_chain_file['relative_path'], 
                            True, 
                            best_chain_file['size']
                        )
            else:
                # Create new record
                if args.fix:
                    logger.info(f"Creating new chain_blast_result record: {best_chain_file['relative_path']}")
                    create_file_record(
                        context,
                        process_id,
                        'chain_blast_result',
                        best_chain_file['relative_path'],
                        True,
                        best_chain_file['size']
                    )
        
        # Check if we need to update domain_blast_result record
        domain_blast_record = next((r for r in file_records if r['file_type'] == 'domain_blast_result'), None)
        
        if domain_blast_files:
            best_domain_file = domain_blast_files[0]
            
            if domain_blast_record:
                if not domain_blast_record['actual_exists'] or domain_blast_record['actual_size'] == 0:
                    # Update record with the found file
                    if args.fix:
                        logger.info(f"Updating domain_blast_result record with new path: {best_domain_file['relative_path']}")
                        update_file_record(
                            context, 
                            domain_blast_record['id'], 
                            best_domain_file['relative_path'], 
                            True, 
                            best_domain_file['size']
                        )
            else:
                # Create new record
                if args.fix:
                    logger.info(f"Creating new domain_blast_result record: {best_domain_file['relative_path']}")
                    create_file_record(
                        context,
                        process_id,
                        'domain_blast_result',
                        best_domain_file['relative_path'],
                        True,
                        best_domain_file['size']
                    )
    
    return 0

if __name__ == "__main__":
    sys.exit(main())