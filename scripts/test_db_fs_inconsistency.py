#!/usr/bin/env python3
"""
Test script to verify database-filesystem inconsistency for a single protein.
This checks whether the database says files exist when they actually don't.
"""

import os
import sys
import logging
import argparse
from typing import Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def check_protein_files(context, protein_id: int, batch_id: int) -> Dict[str, Any]:
    """
    Detailed check of a single protein's files and database records
    
    Args:
        context: Application context
        protein_id: ID of the protein to check
        batch_id: ID of the batch
        
    Returns:
        Dictionary with diagnostic information
    """
    logger = logging.getLogger("test_db_fs_inconsistency")
    
    # Get protein information
    protein_query = """
    SELECT p.id, p.pdb_id, p.chain_id, p.length, ps.id as process_id, ps.status, ps.current_stage
    FROM ecod_schema.protein p
    JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE p.id = %s AND ps.batch_id = %s
    """
    
    result = context.db.execute_query(protein_query, (protein_id, batch_id))
    if not result:
        logger.error(f"Protein {protein_id} not found in batch {batch_id}")
        return {"error": "Protein not found"}
    
    protein_info = {
        "id": result[0][0],
        "pdb_id": result[0][1],
        "chain_id": result[0][2],
        "length": result[0][3],
        "process_id": result[0][4],
        "status": result[0][5],
        "current_stage": result[0][6]
    }
    
    logger.info(f"Checking protein {protein_info['pdb_id']}_{protein_info['chain_id']} (ID: {protein_info['id']})")
    
    # Get batch information to determine base path
    batch_query = """
    SELECT base_path, ref_version 
    FROM ecod_schema.batch 
    WHERE id = %s
    """
    
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return {**protein_info, "error": "Batch not found"}
    
    batch_info = {
        "base_path": batch_result[0][0],
        "ref_version": batch_result[0][1]
    }
    
    # Get file records from the database
    file_query = """
    SELECT id, file_type, file_path, file_exists, file_size, last_checked
    FROM ecod_schema.process_file
    WHERE process_id = %s
    """
    
    file_results = context.db.execute_query(file_query, (protein_info["process_id"],))
    
    file_records = []
    for row in file_results:
        file_record = {
            "id": row[0],
            "file_type": row[1],
            "file_path": row[2],
            "file_exists_db": row[3],
            "file_size": row[4],
            "last_checked": row[5],
            "full_path": os.path.join(os.path.dirname(batch_info["base_path"]), row[2]) if row[2] else None,
            "file_exists_fs": False
        }
        
        # Check if file actually exists in the filesystem
        if file_record["full_path"] and os.path.exists(file_record["full_path"]):
            file_record["file_exists_fs"] = True
            file_record["actual_file_size"] = os.path.getsize(file_record["full_path"])
        else:
            file_record["file_exists_fs"] = False
            file_record["actual_file_size"] = None
        
        file_records.append(file_record)
    
    # Check for expected domain summary file path if not in records
    pdb_id = protein_info["pdb_id"]
    chain_id = protein_info["chain_id"]
    ref_version = batch_info["ref_version"]
    
    # Test multiple potential path formats for domain summary
    potential_paths = [
        os.path.join(batch_info["base_path"], "domains", f"{pdb_id}_{chain_id}.domain_summary.xml"),
        os.path.join(batch_info["base_path"], "domains", f"{pdb_id}_{chain_id}.{ref_version}.domain_summary.xml"),
        os.path.join(batch_info["base_path"], "domains", f"{pdb_id}_{chain_id}.{ref_version}.domains.xml")
    ]
    
    # Also check for non-standard chain IDs (could be issue with unusual chain IDs)
    # Try looking for files that might match our PDB ID
    domain_dir = os.path.join(batch_info["base_path"], "domains")
    matching_files = []
    
    if os.path.exists(domain_dir):
        for filename in os.listdir(domain_dir):
            if filename.startswith(f"{pdb_id}_") and filename.endswith(".xml"):
                matching_files.append(os.path.join(domain_dir, filename))
    
    # Include these in potential paths
    potential_paths.extend(matching_files)
    
    # Check each potential path
    found_files = []
    for path in potential_paths:
        if os.path.exists(path):
            found_files.append({
                "path": path,
                "size": os.path.getsize(path),
                "relative_path": os.path.relpath(path, os.path.dirname(batch_info["base_path"]))
            })
    
    # Determine if there's a database-filesystem inconsistency
    has_inconsistency = False
    domain_summary_in_db = False
    domain_summary_in_fs = len(found_files) > 0
    
    for record in file_records:
        if record["file_type"] == "domain_summary":
            domain_summary_in_db = True
            if record["file_exists_db"] != record["file_exists_fs"]:
                has_inconsistency = True
    
    # Prepare diagnosis
    diagnosis = {
        "protein": protein_info,
        "batch": batch_info,
        "file_records": file_records,
        "found_files": found_files,
        "potential_paths": potential_paths,
        "has_inconsistency": has_inconsistency,
        "domain_summary_in_db": domain_summary_in_db,
        "domain_summary_in_fs": domain_summary_in_fs
    }
    
    # Print detailed diagnosis
    logger.info(f"Database says domain summary exists: {domain_summary_in_db}")
    logger.info(f"Filesystem has matching files: {domain_summary_in_fs}")
    logger.info(f"Database-filesystem inconsistency: {has_inconsistency}")
    
    if domain_summary_in_db:
        for record in file_records:
            if record["file_type"] == "domain_summary":
                logger.info(f"DB record: path={record['file_path']}, exists={record['file_exists_db']}")
                logger.info(f"Filesystem check: exists={record['file_exists_fs']}")
    
    if found_files:
        logger.info(f"Found {len(found_files)} potential domain summary files:")
        for file in found_files:
            logger.info(f"  {file['path']} (size: {file['size']} bytes)")
    
    return diagnosis

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Test database-filesystem inconsistency')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--protein-id', type=int, required=True,
                      help='Protein ID to check')
    parser.add_argument('--batch-id', type=int, default=31,
                      help='Batch ID (default: 31)')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Run diagnosis
    diagnosis = check_protein_files(context, args.protein_id, args.batch_id)
    
    # If inconsistency found, suggest fix
    if diagnosis.get('has_inconsistency', False):
        print("\nSuggested fix for database inconsistency:")
        print("1. Update the process_file record to match reality:")
        
        for record in diagnosis['file_records']:
            if record['file_type'] == 'domain_summary':
                if record['file_exists_db'] and not record['file_exists_fs']:
                    print(f"   - Mark file_exists=FALSE for file ID {record['id']}")
                    print(f"   - SQL: UPDATE ecod_schema.process_file SET file_exists=FALSE WHERE id={record['id']};")
                elif not record['file_exists_db'] and record['file_exists_fs']:
                    print(f"   - Mark file_exists=TRUE for file ID {record['id']}")
                    print(f"   - SQL: UPDATE ecod_schema.process_file SET file_exists=TRUE WHERE id={record['id']};")
        
        # If domain summary exists in DB but not FS, and we found potential files
        domain_summary_in_db = diagnosis.get('domain_summary_in_db', False)
        domain_summary_in_fs = diagnosis.get('domain_summary_in_fs', False)
        
        if domain_summary_in_db and not domain_summary_in_fs and diagnosis.get('found_files'):
            print("\n2. Update the file path to point to an existing file:")
            for file in diagnosis['found_files']:
                print(f"   - Update file_path to '{file['relative_path']}' for the domain_summary record")
                for record in diagnosis['file_records']:
                    if record['file_type'] == 'domain_summary':
                        print(f"   - SQL: UPDATE ecod_schema.process_file SET file_path='{file['relative_path']}', file_exists=TRUE WHERE id={record['id']};")
                        break
                break  # Only suggest the first found file
    
    return 0

if __name__ == "__main__":
    sys.exit(main())