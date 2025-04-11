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

# Add at the end of the script

# Specifically look for chain blast files in common locations
def search_for_blast_files(pdb_id, chain_id, batch_path):
    print("\nSearching for blast files for {}_{}:".format(pdb_id, chain_id))
    
    # Common patterns for blast files
    patterns = [
        f"{pdb_id}_{chain_id}.chainwise_blast.xml",
        f"{pdb_id}_{chain_id}.chain_blast.xml",
        f"{pdb_id}_{chain_id}.blast.xml",
        f"{pdb_id}_{chain_id}_blast.xml"
    ]
    
    # Common directories
    directories = [
        os.path.join(batch_path, "blast", "chain"),
        os.path.join(batch_path, "blast", "chain", "batch_0"),
        os.path.join(batch_path, "blast", "chain", "batch_1"),
        os.path.join(batch_path, "blast", "chain", "batch_2"),
        os.path.join(batch_path, "chain_blast"),
        os.path.join(batch_path, "chain_blast_results")
    ]
    
    for directory in directories:
        if not os.path.exists(directory):
            continue
            
        print(f"  Checking directory: {directory}")
        
        # Check for pattern matches
        for pattern in patterns:
            matching_files = [f for f in os.listdir(directory) if pattern.lower() in f.lower()]
            for file in matching_files:
                full_path = os.path.join(directory, file)
                size = os.path.getsize(full_path) if os.path.exists(full_path) else 0
                print(f"  FOUND: {full_path} (size: {size} bytes)")
        
        # Check for any XML files with PDB ID and chain ID in the name
        all_xml = [f for f in os.listdir(directory) if f.endswith('.xml') and pdb_id.lower() in f.lower() and chain_id.lower() in f.lower()]
        for file in all_xml:
            if not any(pattern.lower() in file.lower() for pattern in patterns):
                full_path = os.path.join(directory, file)
                size = os.path.getsize(full_path) if os.path.exists(full_path) else 0
                print(f"  OTHER XML: {full_path} (size: {size} bytes)")
    
    # Use find command for more thorough search
    import subprocess
    try:
        cmd = f"find {batch_path} -name '*{pdb_id}*{chain_id}*blast*.xml' -type f"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.stdout:
            print("\nFind command results:")
            for line in result.stdout.splitlines():
                size = os.path.getsize(line) if os.path.exists(line) else 0
                print(f"  {line} (size: {size} bytes)")
    except Exception as e:
        print(f"Error running find command: {e}")

if __name__ == "__main__":
    # Run the main checks
    process_id = main()
    
    # Additional blast file search if process_id is available
    if process_id and len(sys.argv) > 2:
        protein_id = int(sys.argv[2])
        batch_path = None
        pdb_id = None
        chain_id = None
        
        # Get batch path and PDB/chain ID
        db = DBManager(ConfigManager().get_db_config())
        query = """
        SELECT b.base_path, p.pdb_id, p.chain_id
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.batch b ON ps.batch_id = b.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.protein_id = %s
        """
        result = db.execute_query(query, (protein_id,))
        
        if result:
            batch_path = result[0][0]
            pdb_id = result[0][1]
            chain_id = result[0][2]
            
            search_for_blast_files(pdb_id, chain_id, batch_path)