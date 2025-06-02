#!/usr/bin/env python3
"""
Example script to prepare reference length files for mini_pyecod

This shows how to create the CSV files needed for accurate coverage calculations.
"""

import csv
import os
from typing import Dict, Tuple

def create_domain_lengths_from_ecod(ecod_domains_file: str, output_file: str):
    """
    Create domain lengths CSV from ECOD domain definitions.
    
    Expected input format (tab-separated):
    ecod_domain_id  pdb  chain  range  length  ...
    e1a0aA1         1a0a A      1-100  100     ...
    """
    domain_lengths = {}
    
    print(f"Reading ECOD domains from {ecod_domains_file}")
    
    with open(ecod_domains_file, 'r') as f:
        # Skip header if present
        header = f.readline()
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                domain_id = parts[0]
                try:
                    length = int(parts[4])
                    domain_lengths[domain_id] = length
                except ValueError:
                    continue
    
    print(f"Found {len(domain_lengths)} domain lengths")
    
    # Write to CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['domain_id', 'length'])
        for domain_id, length in sorted(domain_lengths.items()):
            writer.writerow([domain_id, length])
    
    print(f"Wrote domain lengths to {output_file}")
    return domain_lengths

def create_protein_lengths_from_fasta(fasta_dir: str, output_file: str):
    """
    Create protein lengths CSV from FASTA files.
    
    Expects files named like: 1a0a_A.fasta
    """
    protein_lengths = {}
    
    print(f"Scanning FASTA files in {fasta_dir}")
    
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fasta'):
            # Parse PDB and chain from filename
            base = filename.replace('.fasta', '')
            if '_' in base:
                parts = base.split('_')
                pdb_id = parts[0]
                chain_id = parts[1] if len(parts) > 1 else 'A'
                
                # Read sequence length
                filepath = os.path.join(fasta_dir, filename)
                with open(filepath, 'r') as f:
                    sequence = ''
                    for line in f:
                        if not line.startswith('>'):
                            sequence += line.strip()
                
                if sequence:
                    protein_lengths[(pdb_id, chain_id)] = len(sequence)
    
    print(f"Found {len(protein_lengths)} protein lengths")
    
    # Write to CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pdb_id', 'chain_id', 'length'])
        for (pdb_id, chain_id), length in sorted(protein_lengths.items()):
            writer.writerow([pdb_id, chain_id, length])
    
    print(f"Wrote protein lengths to {output_file}")
    return protein_lengths

def create_reference_from_database(db_connection_string: str):
    """
    Example of creating reference files from a database.
    
    This is a template - adapt to your database schema.
    """
    print("Database example - adapt this to your schema:")
    
    print("""
    -- Query for domain lengths
    SELECT domain_id, length 
    FROM ecod_domains 
    WHERE length IS NOT NULL;
    
    -- Query for protein lengths  
    SELECT pdb_id, chain_id, length
    FROM pdb_chains
    WHERE length IS NOT NULL;
    """)

def validate_reference_files(domain_file: str, protein_file: str):
    """
    Validate reference files have the expected format.
    """
    print("\nValidating reference files...")
    
    # Check domain lengths
    if os.path.exists(domain_file):
        with open(domain_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header != ['domain_id', 'length']:
                print(f"WARNING: Unexpected header in {domain_file}: {header}")
            
            # Check a few rows
            valid_rows = 0
            for i, row in enumerate(reader):
                if i >= 5:  # Just check first 5
                    break
                if len(row) == 2 and row[1].isdigit():
                    valid_rows += 1
                    print(f"  ✓ {row[0]}: {row[1]} residues")
            
            if valid_rows > 0:
                print(f"Domain lengths file looks valid")
    else:
        print(f"Domain lengths file not found: {domain_file}")
    
    # Check protein lengths
    if os.path.exists(protein_file):
        with open(protein_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header != ['pdb_id', 'chain_id', 'length']:
                print(f"WARNING: Unexpected header in {protein_file}: {header}")
            
            # Check a few rows
            valid_rows = 0
            for i, row in enumerate(reader):
                if i >= 5:  # Just check first 5
                    break
                if len(row) == 3 and row[2].isdigit():
                    valid_rows += 1
                    print(f"  ✓ {row[0]}_{row[1]}: {row[2]} residues")
            
            if valid_rows > 0:
                print(f"Protein lengths file looks valid")
    else:
        print(f"Protein lengths file not found: {protein_file}")

def main():
    """Example usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Prepare reference length files for mini_pyecod')
    parser.add_argument('--ecod-domains', help='ECOD domains file (tab-separated)')
    parser.add_argument('--fasta-dir', help='Directory containing FASTA files')
    parser.add_argument('--validate', action='store_true', help='Validate existing files')
    parser.add_argument('--domain-output', default='domain_lengths.csv', 
                        help='Output file for domain lengths')
    parser.add_argument('--protein-output', default='protein_lengths.csv',
                        help='Output file for protein lengths')
    
    args = parser.parse_args()
    
    if args.validate:
        validate_reference_files(args.domain_output, args.protein_output)
    else:
        if args.ecod_domains:
            create_domain_lengths_from_ecod(args.ecod_domains, args.domain_output)
        
        if args.fasta_dir:
            create_protein_lengths_from_fasta(args.fasta_dir, args.protein_output)
        
        if not args.ecod_domains and not args.fasta_dir:
            print("Example usage:")
            print("  # Create from ECOD domains file")
            print("  python prepare_reference_lengths.py --ecod-domains ecod.latest.domains.txt")
            print()
            print("  # Create from FASTA directory")
            print("  python prepare_reference_lengths.py --fasta-dir /data/fasta/")
            print()
            print("  # Validate existing files")
            print("  python prepare_reference_lengths.py --validate")
            print()
            create_reference_from_database("example")

if __name__ == "__main__":
    main()
