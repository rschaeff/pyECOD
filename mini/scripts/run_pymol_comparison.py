#!/usr/bin/env python3
"""
Run PyMOL comparison for 8ovp_A using existing mini_pyecod components

This script:
1. Checks that required files exist (old domains, new domains, structure)
2. Uses existing PyMOL comparison components
3. Generates side-by-side comparison visualization
4. Provides helpful error messages and next steps

Usage:
    python run_pymol_comparison.py 8ovp_A
    python run_pymol_comparison.py 8ovp_A --batch-dir /custom/batch/path
"""

import sys
import os
from pathlib import Path

# Add mini directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def check_required_files(protein_id, batch_dir, output_dir="/tmp"):
    """Check that all required files exist for comparison"""
    
    print(f"Checking required files for {protein_id}...")
    print("=" * 60)
    
    # Expected file paths
    files_to_check = {
        'new_domains': f"{output_dir}/{protein_id}_mini.domains.xml",
        'new_domains_alt': f"{output_dir}/{protein_id}_range_cache.domains.xml", 
        'old_domains': f"{batch_dir}/domains/{protein_id}.develop291.domains.xml",
        'domain_summary': f"{batch_dir}/domains/{protein_id}.develop291.domain_summary.xml"
    }
    
    # Check which files exist
    found_files = {}
    missing_files = []
    
    for file_type, file_path in files_to_check.items():
        if os.path.exists(file_path):
            found_files[file_type] = file_path
            print(f"✅ {file_type}: {file_path}")
        else:
            missing_files.append((file_type, file_path))
            print(f"❌ {file_type}: {file_path}")
    
    # Determine which new domains file to use
    new_domains_file = None
    if 'new_domains' in found_files:
        new_domains_file = found_files['new_domains']
    elif 'new_domains_alt' in found_files:
        new_domains_file = found_files['new_domains_alt']
    
    # Check critical files
    critical_missing = []
    if not new_domains_file:
        critical_missing.append("new domain partitions (run quick_test.py first)")
    if 'old_domains' not in found_files:
        critical_missing.append("old domain partitions")
    
    if critical_missing:
        print(f"\n❌ CRITICAL FILES MISSING:")
        for missing in critical_missing:
            print(f"   - {missing}")
        print(f"\nTo fix:")
        if not new_domains_file:
            print(f"   1. Run: python mini/scripts/quick_test.py {protein_id} --reference-lengths test_data/domain_lengths.csv")
        if 'old_domains' not in found_files:
            print(f"   2. Check that old domains exist in batch directory")
        return False, None, None
    
    old_domains_file = found_files['old_domains']
    
    print(f"\n✅ Ready for comparison:")
    print(f"   OLD: {old_domains_file}")
    print(f"   NEW: {new_domains_file}")
    
    return True, old_domains_file, new_domains_file

def run_comparison(protein_id, old_domains_file, new_domains_file, pdb_repo_path="/usr2/pdb/data"):
    """Run the PyMOL comparison using existing components"""
    
    print(f"\nGenerating PyMOL comparison for {protein_id}...")
    print("=" * 60)
    
    try:
        # Import existing comparison components
        from mini.scripts.pymol_comparison import PDBRepository, DomainParser, PyMOLVisualizer
        
        # Initialize components
        pdb_repo = PDBRepository(pdb_repo_path)
        parser = DomainParser()
        visualizer = PyMOLVisualizer(pdb_repo)
        
        # Parse domains from both files
        print("Parsing domain assignments...")
        old_domains = parser.parse_old_domains(old_domains_file)
        new_domains = parser.parse_new_domains(new_domains_file)
        
        print(f"Old algorithm: {len(old_domains)} domains")
        for i, domain in enumerate(old_domains, 1):
            print(f"  {i}. {domain.family}: {domain.range_str}")
        
        print(f"New algorithm: {len(new_domains)} domains")
        for i, domain in enumerate(new_domains, 1):
            print(f"  {i}. {domain.family}: {domain.range_str}")
        
        # Create PyMOL script
        output_dir = "/tmp/pymol_comparison"
        script_path = visualizer.create_comparison_script(
            protein_id, old_domains, new_domains, output_dir
        )
        
        print(f"\n✅ PyMOL script created: {script_path}")
        print(f"   Session will be saved to: {output_dir}/{protein_id}_comparison.pse")
        
        return script_path
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("Make sure you're running from the correct directory")
        return None
    except Exception as e:
        print(f"❌ Error generating comparison: {e}")
        import traceback
        traceback.print_exc()
        return None

def display_next_steps(script_path, protein_id):
    """Display next steps for the user"""
    
    print(f"\n" + "=" * 60)
    print("NEXT STEPS")
    print("=" * 60)
    
    if script_path:
        print(f"1. Launch PyMOL with the comparison:")
        print(f"   pymol {script_path}")
        print(f"")
        print(f"2. In PyMOL, you'll see:")
        print(f"   - LEFT SIDE: Previous algorithm results")
        print(f"   - RIGHT SIDE: Mini PyECOD results")
        print(f"   - Different colors for each domain")
        print(f"")
        print(f"3. The session is automatically saved to:")
        print(f"   /tmp/pymol_comparison/{protein_id}_comparison.pse")
        print(f"")
        print(f"4. To reload later:")
        print(f"   pymol /tmp/pymol_comparison/{protein_id}_comparison.pse")
    else:
        print(f"❌ Could not generate PyMOL script")
        print(f"Check error messages above")

def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Generate PyMOL comparison between old and new domain partitioning',
        epilog='Uses existing mini_pyecod PyMOL comparison components'
    )
    
    parser.add_argument('protein_id', nargs='?', default='8ovp_A',
                        help='Protein ID to compare (default: 8ovp_A)')
    parser.add_argument('--batch-dir',
                        default='/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424',
                        help='Batch directory containing old domain files')
    parser.add_argument('--pdb-repo', default='/usr2/pdb/data',
                        help='PDB repository path')
    parser.add_argument('--output-dir', default='/tmp',
                        help='Directory containing new domain files')
    
    args = parser.parse_args()
    
    print("PyMOL Domain Partitioning Comparison")
    print("=" * 60)
    print(f"Protein: {args.protein_id}")
    print(f"Batch dir: {args.batch_dir}")
    print(f"Output dir: {args.output_dir}")
    print(f"PDB repo: {args.pdb_repo}")
    
    # Check required files
    files_ok, old_domains_file, new_domains_file = check_required_files(
        args.protein_id, args.batch_dir, args.output_dir
    )
    
    if not files_ok:
        sys.exit(1)
    
    # Run comparison
    script_path = run_comparison(
        args.protein_id, old_domains_file, new_domains_file, args.pdb_repo
    )
    
    # Display next steps
    display_next_steps(script_path, args.protein_id)

if __name__ == "__main__":
    main()
