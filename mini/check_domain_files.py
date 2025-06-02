#!/usr/bin/env python3
"""
Check available domain files and their contents for PyMOL comparison

This script helps troubleshoot issues by showing:
1. Which domain files exist
2. What domains are in each file
3. File format compatibility
"""

import sys
import os
import xml.etree.ElementTree as ET
from pathlib import Path

def check_domain_file_content(file_path, file_type="unknown"):
    """Check the content of a domain file"""
    
    if not os.path.exists(file_path):
        return f"❌ File not found: {file_path}"
    
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Count domains
        domains = root.findall(".//domain")
        
        if not domains:
            return f"⚠️  No domains found in {file_type} file"
        
        # Extract domain info
        domain_info = []
        for domain in domains:
            domain_id = domain.get("id", domain.get("domain_id", "unknown"))
            range_str = domain.get("range", "unknown")
            family = domain.get("family", domain.get("source_id", domain.get("t_group", "unknown")))
            
            domain_info.append({
                'id': domain_id,
                'range': range_str, 
                'family': family
            })
        
        result = f"✅ {len(domains)} domains in {file_type} file:\n"
        for i, info in enumerate(domain_info, 1):
            result += f"     {i}. {info['family']}: {info['range']}\n"
        
        return result.rstrip()
        
    except ET.ParseError as e:
        return f"❌ XML parsing error in {file_type} file: {e}"
    except Exception as e:
        return f"❌ Error reading {file_type} file: {e}"

def check_all_domain_files(protein_id, batch_dir="/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"):
    """Check all possible domain files for a protein"""
    
    print(f"Checking domain files for {protein_id}")
    print("=" * 60)
    
    # Possible file locations
    files_to_check = {
        "New domains (mini)": f"/tmp/{protein_id}_mini.domains.xml",
        "New domains (range_cache)": f"/tmp/{protein_id}_range_cache.domains.xml",
        "Old domains": f"{batch_dir}/domains/{protein_id}.develop291.domains.xml",
        "Domain summary": f"{batch_dir}/domains/{protein_id}.develop291.domain_summary.xml"
    }
    
    found_files = []
    
    for file_type, file_path in files_to_check.items():
        print(f"\n{file_type}:")
        print(f"  Path: {file_path}")
        
        if "summary" in file_type.lower():
            # For domain summary, just check existence
            if os.path.exists(file_path):
                print(f"  ✅ File exists")
                found_files.append((file_type, file_path))
            else:
                print(f"  ❌ File not found")
        else:
            # For domain files, check content
            result = check_domain_file_content(file_path, file_type)
            print(f"  {result}")
            
            if result.startswith("✅"):
                found_files.append((file_type, file_path))
    
    # Summary
    print(f"\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    domain_files = [f for f in found_files if "domain" in f[0].lower() and "summary" not in f[0].lower()]
    
    if len(domain_files) >= 2:
        print("✅ Ready for PyMOL comparison!")
        print("\nAvailable domain files:")
        for file_type, file_path in domain_files:
            print(f"  - {file_type}: {file_path}")
        
        print(f"\nTo generate comparison:")
        print(f"  python mini/scripts/run_pymol_comparison.py {protein_id}")
        
    else:
        print("❌ Not ready for comparison")
        print(f"Found {len(domain_files)} domain files, need at least 2 (old + new)")
        
        if not any("new" in f[0].lower() for f in found_files):
            print(f"\nTo generate new domains:")
            print(f"  python mini/scripts/quick_test.py {protein_id} --reference-lengths test_data/domain_lengths.csv")
        
        if not any("old" in f[0].lower() for f in found_files):
            print(f"\nCheck that old domain files exist in batch directory:")
            print(f"  ls {batch_dir}/domains/{protein_id}*")

def check_pymol_requirements():
    """Check if PyMOL and structure files are available"""
    
    print(f"\n" + "=" * 60)
    print("PYMOL REQUIREMENTS CHECK")
    print("=" * 60)
    
    # Check PyMOL availability
    import shutil
    pymol_path = shutil.which("pymol")
    if pymol_path:
        print(f"✅ PyMOL found: {pymol_path}")
    else:
        print(f"⚠️  PyMOL not found in PATH")
        print(f"   Make sure PyMOL is installed and available")
    
    # Check PDB repository
    pdb_paths = [
        "/usr2/pdb/data",
        "/data/pdb", 
        "/pdb"
    ]
    
    pdb_found = False
    for pdb_path in pdb_paths:
        if os.path.exists(pdb_path):
            print(f"✅ PDB repository found: {pdb_path}")
            pdb_found = True
            break
    
    if not pdb_found:
        print(f"⚠️  PDB repository not found")
        print(f"   Checked: {', '.join(pdb_paths)}")
        print(f"   Specify custom path with --pdb-repo")

def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Check domain files for PyMOL comparison'
    )
    parser.add_argument('protein_id', nargs='?', default='8ovp_A',
                        help='Protein ID to check (default: 8ovp_A)')
    parser.add_argument('--batch-dir',
                        default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                        help='Batch directory')
    parser.add_argument('--check-pymol', action='store_true',
                        help='Also check PyMOL requirements')
    
    args = parser.parse_args()
    
    # Check domain files
    check_all_domain_files(args.protein_id, args.batch_dir)
    
    # Optionally check PyMOL requirements
    if args.check_pymol:
        check_pymol_requirements()

if __name__ == "__main__":
    main()
