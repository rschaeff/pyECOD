"""
Utility script to manage the reference blacklist
"""

import csv
import sys
from pathlib import Path
from datetime import date

def add_to_blacklist(pdb_id: str, chain_id: str, reason: str, 
                    blacklist_path: str = "test_data/reference_blacklist.csv"):
    """Add an entry to the blacklist"""
    
    # Check if file exists, create with header if not
    if not Path(blacklist_path).exists():
        with open(blacklist_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['pdb_id', 'chain_id', 'reason', 'date_added', 'added_by'])
    
    # Add new entry
    with open(blacklist_path, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            pdb_id.lower(),
            chain_id,
            reason,
            str(date.today()),
            "manual"
        ])
    
    print(f"Added {pdb_id}_{chain_id} to blacklist: {reason}")

def list_blacklist(blacklist_path: str = "test_data/reference_blacklist.csv"):
    """List all blacklisted entries"""
    
    if not Path(blacklist_path).exists():
        print("No blacklist file found")
        return
    
    print("Blacklisted reference chains:")
    print("=" * 60)
    
    with open(blacklist_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            print(f"{row['pdb_id']}_{row['chain_id']}")
            print(f"  Reason: {row['reason']}")
            print(f"  Added: {row['date_added']} by {row['added_by']}")
            print()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python manage_blacklist.py list")
        print("  python manage_blacklist.py add PDB_ID CHAIN_ID 'Reason for blacklisting'")
        sys.exit(1)
    
    cmd = sys.argv[1]
    
    if cmd == "list":
        list_blacklist()
    elif cmd == "add" and len(sys.argv) >= 5:
        pdb_id, chain_id, reason = sys.argv[2], sys.argv[3], sys.argv[4]
        add_to_blacklist(pdb_id, chain_id, reason)
    else:
        print("Invalid command or arguments")
