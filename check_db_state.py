#!/usr/bin/env python3
"""
check_db_state.py - Check the current state of the database
"""

import sys
import os
import argparse
import logging
from typing import Dict, Any, List

# Add parent directory to path if needed
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from ecod.core.config import ConfigManager
from ecod.core.db_manager import DBManager

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

def check_proteins():
    """Check proteins table"""
    config_manager = ConfigManager('config/config.yml')
    db = DBManager(config_manager.get_db_config())
    
    # Check proteins
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, p.source_id, p.length, 
        CASE WHEN ps.id IS NOT NULL THEN true ELSE false END AS has_sequence
    FROM 
        ecod_schema.protein p
    LEFT JOIN
        ecod_schema.protein_sequence ps ON p.id = ps.protein_id
    ORDER BY 
        p.id
    """
    rows = db.execute_dict_query(query)
    
    print(f"Found {len(rows)} proteins:")
    for row in rows:
        sequence_status = "with sequence" if row['has_sequence'] else "no sequence"
        print(f"  - ID: {row['id']}, PDB: {row['pdb_id']}, Chain: {row['chain_id']}, " 
              f"Length: {row['length']}, {sequence_status}")
    
    # Check process status
    query = """
    SELECT 
        ps.id, p.pdb_id, p.chain_id, ps.batch_id, ps.current_stage, ps.status
    FROM 
        ecod_schema.process_status ps
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    ORDER BY 
        ps.id
    """
    rows = db.execute_dict_query(query)
    
    print(f"\nFound {len(rows)} process status records:")
    for row in rows:
        print(f"  - Process ID: {row['id']}, PDB: {row['pdb_id']}, Chain: {row['chain_id']}, "
              f"Batch: {row['batch_id']}, Stage: {row['current_stage']}, Status: {row['status']}")
    
    # Check batches
    query = """
    SELECT 
        id, batch_name, type, total_items, completed_items, status
    FROM 
        ecod_schema.batch
    ORDER BY 
        id
    """
    rows = db.execute_dict_query(query)
    
    print(f"\nFound {len(rows)} batches:")
    for row in rows:
        print(f"  - Batch ID: {row['id']}, Name: {row['batch_name']}, Type: {row['type']}, "
              f"Progress: {row['completed_items']}/{row['total_items']}, Status: {row['status']}")
    
    # Check files
    query = """
    SELECT 
        pf.id, ps.id AS process_id, p.pdb_id, p.chain_id, 
        pf.file_type, pf.file_path, pf.file_exists
    FROM 
        ecod_schema.process_file pf
    JOIN
        ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    ORDER BY 
        pf.id
    """
    rows = db.execute_dict_query(query)
    
    print(f"\nFound {len(rows)} file records:")
    for row in rows:
        exists_status = "exists" if row['file_exists'] else "missing"
        print(f"  - File ID: {row['id']}, Process: {row['process_id']}, "
              f"PDB: {row['pdb_id']}, Chain: {row['chain_id']}, "
              f"Type: {row['file_type']}, Path: {row['file_path']}, {exists_status}")

def main():
    parser = argparse.ArgumentParser(description='Check PyECOD Database State')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    check_proteins()

if __name__ == "__main__":
    main()