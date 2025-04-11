#!/usr/bin/env python3
"""
check_protein_files.py - Check file records for a specific protein
"""

import os
import sys
import argparse
# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from ecod.core.context import ApplicationContext

def main():
    parser = argparse.ArgumentParser(description='Check file records for a protein')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--protein-id', type=int, required=True,
                      help='Protein ID to check')
    
    args = parser.parse_args()
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get protein information
    protein_query = """
    SELECT id, pdb_id, chain_id, source_id
    FROM ecod_schema.protein
    WHERE id = %s
    """
    
    protein = context.db.execute_query(protein_query, (args.protein_id,))
    if not protein:
        print(f"Protein ID {args.protein_id} not found in database")
        return 1
    
    pdb_id = protein[0][1]
    chain_id = protein[0][2]
    print(f"Checking files for protein: {pdb_id}_{chain_id} (ID: {args.protein_id})")
    
    # Get process_id
    process_query = """
    SELECT ps.id, b.id, b.batch_name, b.base_path
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.batch b ON ps.batch_id = b.id
    WHERE ps.protein_id = %s
    """
    
    process = context.db.execute_query(process_query, (args.protein_id,))
    if not process:
        print(f"No process record found for protein ID {args.protein_id}")
        return 1
    
    process_id = process[0][0]
    batch_id = process[0][1]
    batch_name = process[0][2]
    batch_path = process[0][3]
    print(f"Process ID: {process_id}, Batch: {batch_name} (ID: {batch_id})")
    print(f"Batch path: {batch_path}")
    
    # Get file records
    files_query = """
    SELECT file_type, file_path, file_exists, file_size, last_checked
    FROM ecod_schema.process_file
    WHERE process_id = %s
    """
    
    files = context.db.execute_query(files_query, (process_id,))
    print(f"Found {len(files)} file records:")
    
    for file in files:
        file_type = file[0]
        file_path = file[1]
        file_exists = file[2]
        file_size = file[3]
        
        # Construct full path
        full_path = os.path.join(batch_path, file_path)
        actual_exists = os.path.exists(full_path)
        
        print(f"  {file_type}:")
        print(f"    Path: {file_path}")
        print(f"    Full path: {full_path}")
        print(f"    DB says exists: {file_exists}, Actually exists: {actual_exists}")
        print(f"    Size: {file_size} bytes")
        
        # Print potential path normalization issues
        if ".." in file_path:
            print(f"    WARNING: Path contains '..' which may cause resolution issues")
            normalized_path = os.path.normpath(full_path)
            print(f"    Normalized path: {normalized_path}")
            norm_exists = os.path.exists(normalized_path)
            print(f"    Normalized path exists: {norm_exists}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())