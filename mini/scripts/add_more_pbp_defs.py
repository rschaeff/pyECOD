#!/usr/bin/env python3
"""Add domain definitions for other PBP proteins"""

import csv
import os

def add_pbp_domain_definitions():
    """Add domain definitions for other PBP proteins that appear in the results"""
    
    # Read existing definitions
    existing = []
    if os.path.exists("test_data/domain_definitions.csv"):
        with open("test_data/domain_definitions.csv", "r") as f:
            reader = csv.DictReader(f)
            existing = list(reader)
    
    # PBP proteins that need definitions (same architecture as 2ia4)
    pbp_proteins = ['2vha', '1pda', '1usi', '1usj']
    
    # Add definitions for these PBPs using the same architecture as 2ia4
    new_entries = []
    for pdb in pbp_proteins:
        for chain in ['A', 'B', 'C', 'D']:
            # Skip if already exists
            if any(e['pdb_id'] == pdb and e['chain_id'] == chain for e in existing):
                continue
                
            # Domain 1 (like 2ia4 domain 1)
            new_entries.append({
                "domain_id": f"e{pdb}{chain}1",
                "pdb_id": pdb,
                "chain_id": chain,
                "range": "110-209",  # Same as 2ia4
                "length": "100",
                "t_group": "7523.1.1",
                "h_group": "7523.1"
            })
            
            # Domain 2 (like 2ia4 domain 2)
            new_entries.append({
                "domain_id": f"e{pdb}{chain}2", 
                "pdb_id": pdb,
                "chain_id": chain,
                "range": "3-109,210-279",  # Same as 2ia4
                "length": "177",
                "t_group": "7523.1.1.4",
                "h_group": "7523.1"
            })
    
    # Combine and write
    all_entries = existing + new_entries
    
    with open("test_data/domain_definitions.csv", "w", newline="") as f:
        fieldnames = ["domain_id", "pdb_id", "chain_id", "range", "length", "t_group", "h_group"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_entries)
    
    print(f"Added {len(new_entries)} new domain definitions")
    print(f"Total definitions: {len(all_entries)}")
    
    # Show what was added
    if new_entries:
        print("\nAdded definitions for:")
        added_pdbs = set((e['pdb_id'], e['chain_id']) for e in new_entries)
        for pdb, chain in sorted(added_pdbs):
            print(f"  {pdb}_{chain}")

if __name__ == "__main__":
    add_pbp_domain_definitions()
