#!/usr/bin/env python3
"""
Quick test script to debug the PyMOL comparison issue

Let's create a simplified version to test step by step.
"""

import sys
import os
import xml.etree.ElementTree as ET
from pathlib import Path

def test_file_access():
    """Test that we can access the required files"""
    
    protein_id = "8ovp_A"
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    
    print("Testing file access...")
    
    # Check new domains file
    new_domains_file = f"/tmp/{protein_id}_range_cache.domains.xml"
    print(f"New domains file: {new_domains_file}")
    print(f"  Exists: {os.path.exists(new_domains_file)}")
    
    # Check old domains file
    old_domains_file = os.path.join(batch_dir, "domains", f"{protein_id}.develop291.domains.xml")
    print(f"Old domains file: {old_domains_file}")
    print(f"  Exists: {os.path.exists(old_domains_file)}")
    
    # Check structure file
    pdb_id = protein_id.split('_')[0]
    structure_file = f"/usr2/pdb/data/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif.gz"
    print(f"Structure file: {structure_file}")
    print(f"  Exists: {os.path.exists(structure_file)}")
    
    return all([
        os.path.exists(new_domains_file),
        os.path.exists(old_domains_file),
        os.path.exists(structure_file)
    ])

def parse_domains_simple(xml_path: str, is_old: bool = False):
    """Simple domain parsing for testing"""
    
    print(f"\nParsing {'old' if is_old else 'new'} domains from: {xml_path}")
    
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        domains = []
        
        if is_old:
            # Old format: <domain start="5" end="245" range="5-245" source="hhsearch" source_id="e4q0cB4" ...>
            for domain_elem in root.findall(".//domain"):
                range_str = domain_elem.get("range", "")
                source_id = domain_elem.get("source_id", "")
                domains.append({
                    "range": range_str,
                    "family": source_id,
                    "source": domain_elem.get("source", "")
                })
        else:
            # New format: <domain id="d1" range="2-248,491-517" family="2ia4" ...>
            for domain_elem in root.findall(".//domain"):
                range_str = domain_elem.get("range", "")
                family = domain_elem.get("family", "")
                domains.append({
                    "range": range_str,
                    "family": family,
                    "source": domain_elem.get("source", "")
                })
        
        print(f"  Found {len(domains)} domains:")
        for i, domain in enumerate(domains, 1):
            print(f"    {i}. {domain['family']}: {domain['range']}")
        
        return domains
        
    except Exception as e:
        print(f"  Error: {e}")
        return []

def create_simple_pymol_script(protein_id: str, old_domains: list, new_domains: list):
    """Create a simple PyMOL script for testing"""
    
    pdb_id = protein_id.split('_')[0]
    chain_id = protein_id.split('_')[1]
    structure_path = f"/usr2/pdb/data/structures/divided/mmCIF/{pdb_id[1:3]}/{pdb_id}.cif.gz"
    
    script_path = f"/tmp/{protein_id}_simple_comparison.pml"
    
    colors = ["red", "blue", "green", "cyan", "magenta", "orange"]
    
    with open(script_path, 'w') as f:
        f.write(f"# Simple PyMOL comparison for {protein_id}\n\n")
        
        # Load structure
        f.write(f"load {structure_path}, {pdb_id}\n")
        f.write(f"create old_alg, {pdb_id}\n")
        f.write(f"create new_alg, {pdb_id}\n")
        f.write(f"delete {pdb_id}\n\n")
        
        # Basic setup
        f.write("hide everything\n")
        f.write("show cartoon\n")
        f.write("color gray90, all\n\n")
        
        # Color old domains
        f.write(f"# OLD algorithm ({len(old_domains)} domains)\n")
        for i, domain in enumerate(old_domains):
            color = colors[i % len(colors)]
            range_str = domain['range']
            if '-' in range_str:
                # Handle simple and discontinuous ranges
                for segment in range_str.split(','):
                    if '-' in segment:
                        f.write(f"color {color}, old_alg and chain {chain_id} and resi {segment}\n")
        
        f.write(f"\n# NEW algorithm ({len(new_domains)} domains)\n")
        for i, domain in enumerate(new_domains):
            color = colors[i % len(colors)]
            range_str = domain['range']
            if '-' in range_str:
                for segment in range_str.split(','):
                    if '-' in segment:
                        f.write(f"color {color}, new_alg and chain {chain_id} and resi {segment}\n")
        
        # Grid mode
        f.write(f"\n# Grid mode setup\n")
        f.write("set grid_mode, 1\n")
        f.write("set grid_slot, 1, old_alg\n")
        f.write("set grid_slot, 2, new_alg\n")
        f.write("zoom all\n")
        f.write("orient\n")
        
        # Save session
        f.write(f"\nsave /tmp/{protein_id}_simple.pse\n")
    
    return script_path

def main():
    """Main test function"""
    
    print("PYMOL COMPARISON DEBUG")
    print("=" * 40)
    
    # Test file access
    if not test_file_access():
        print("\n❌ Missing required files!")
        print("Make sure you've run:")
        print("  python run_range_cache_test.py --protein 8ovp_A")
        return
    
    # Parse domains
    protein_id = "8ovp_A"
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    
    new_domains_file = f"/tmp/{protein_id}_range_cache.domains.xml"
    old_domains_file = os.path.join(batch_dir, "domains", f"{protein_id}.develop291.domains.xml")
    
    old_domains = parse_domains_simple(old_domains_file, is_old=True)
    new_domains = parse_domains_simple(new_domains_file, is_old=False)
    
    if not old_domains or not new_domains:
        print("\n❌ Failed to parse domain files!")
        return
    
    # Create simple script
    script_path = create_simple_pymol_script(protein_id, old_domains, new_domains)
    
    print(f"\n✅ Simple PyMOL script created: {script_path}")
    print(f"Test with: pymol {script_path}")
    print(f"Session saved to: /tmp/{protein_id}_simple.pse")
    
    print(f"\nComparison:")
    print(f"  Old algorithm: {len(old_domains)} domains")
    print(f"  New algorithm: {len(new_domains)} domains")

if __name__ == "__main__":
    main()
