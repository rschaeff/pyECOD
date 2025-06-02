#!/usr/bin/env python3
"""Debug domain definitions loading"""

import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.decomposer import load_domain_definitions

def debug_domain_definitions():
    """Check what domain definitions we have"""
    
    domain_defs_csv = "test_data/domain_definitions.csv"
    
    if not os.path.exists(domain_defs_csv):
        print(f"Domain definitions file not found: {domain_defs_csv}")
        return
    
    domain_defs = load_domain_definitions(domain_defs_csv)
    
    print(f"\nLoaded definitions for {len(domain_defs)} proteins:")
    for (pdb, chain), domains in sorted(domain_defs.items())[:10]:
        print(f"\n{pdb}_{chain}:")
        for d in domains:
            print(f"  {d.domain_id}: {d.range} ({d.length} residues)")
            if d.t_group:
                print(f"    T-group: {d.t_group}, H-group: {d.h_group}")
    
    # Check specifically for 2ia4
    if ('2ia4', 'A') in domain_defs:
        print(f"\n2ia4_A domain architecture:")
        for d in domain_defs[('2ia4', 'A')]:
            print(f"  {d.domain_id}: {d.range} ({d.length} residues)")
    else:
        print("\nWARNING: No domain definitions found for 2ia4_A")
        print("This is why decomposition isn't working!")

if __name__ == "__main__":
    debug_domain_definitions()
