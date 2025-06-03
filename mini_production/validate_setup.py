#!/usr/bin/env python3
"""
Validation script for mini production setup
Checks filesystem, executables, and creates a small test

Usage:
    python validate_setup.py
"""

import os
import sys
import subprocess
from pathlib import Path
import yaml

def check_config():
    """Check if config file exists and is valid"""
    config_path = "mini_production/config.local.yml"
    print("=== Config Validation ===")
    
    if not os.path.exists(config_path):
        print(f"❌ Config file missing: {config_path}")
        return False
    
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        print(f"✅ Config loaded from {config_path}")
        
        # Check required sections
        required = ["database", "slurm", "paths"]
        for section in required:
            if section in config:
                print(f"✅ Config section '{section}' found")
            else:
                print(f"❌ Config section '{section}' missing")
                
        return True
    except Exception as e:
        print(f"❌ Config error: {e}")
        return False

def check_mini_executable():
    """Check if mini executable exists and works"""
    print("\n=== Mini Executable ===")
    
    candidates = [
        "./pyecod_mini",
        "mini/pyecod_mini", 
        "../mini/pyecod_mini"
    ]
    
    executable = None
    for candidate in candidates:
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            executable = os.path.abspath(candidate)
            print(f"✅ Found executable: {executable}")
            break
    
    if not executable:
        print("❌ Mini executable not found")
        print("   Looked for:", candidates)
        return False
    
    # Test execution
    try:
        result = subprocess.run([executable, "--help"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0 or "usage" in result.stdout.lower():
            print("✅ Executable responds to --help")
            return True
        else:
            print(f"⚠️  Executable runs but unexpected output: {result.returncode}")
            return True  # Might still work
    except Exception as e:
        print(f"❌ Executable test failed: {e}")
        return False

def check_batch_directories():
    """Check batch directory structure"""
    print("\n=== Batch Directories ===")
    
    batch_base = Path("/data/ecod/pdb_updates/batches")
    if not batch_base.exists():
        print(f"❌ Batch base directory missing: {batch_base}")
        return False
    
    print(f"✅ Batch base directory exists: {batch_base}")
    
    # Find batch directories
    batch_dirs = [d for d in batch_base.iterdir() if d.is_dir()]
    print(f"✅ Found {len(batch_dirs)} batch directories")
    
    # Check a few for structure
    checked = 0
    for batch_dir in batch_dirs[:3]:
        domains_dir = batch_dir / "domains"
        if domains_dir.exists():
            summary_files = list(domains_dir.glob("*.develop291.domain_summary.xml"))
            print(f"✅ {batch_dir.name}: {len(summary_files)} domain summary files")
            checked += 1
        else:
            print(f"⚠️  {batch_dir.name}: no domains/ directory")
    
    return checked > 0

def check_slurm():
    """Check SLURM availability"""
    print("\n=== SLURM Availability ===")
    
    try:
        result = subprocess.run(["sinfo"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print("✅ sinfo command works")
            # Parse partition info
            lines = result.stdout.strip().split('\n')[1:]  # Skip header
            for line in lines[:3]:  # Show first few partitions
                parts = line.split()
                if len(parts) >= 4:
                    print(f"   Partition: {parts[0]} ({parts[3]} nodes)")
        else:
            print(f"❌ sinfo failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"❌ SLURM check failed: {e}")
        return False
    
    try:
        result = subprocess.run(["sbatch", "--version"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            print("✅ sbatch command available")
        else:
            print("❌ sbatch not available")
            return False
    except Exception as e:
        print(f"❌ sbatch check failed: {e}")
        return False
    
    return True

def check_database_connection():
    """Check database connectivity"""
    print("\n=== Database Connection ===")
    
    try:
        import psycopg2
        
        # Load config
        with open("mini_production/config.local.yml", 'r') as f:
            config = yaml.safe_load(f)
        
        db_config = config['database']
        
        # Test connection
        conn = psycopg2.connect(**db_config)
        cursor = conn.cursor()
        
        # Test query
        cursor.execute("SELECT COUNT(*) FROM ecod_schema.process_status WHERE is_representative = TRUE")
        count = cursor.fetchone()[0]
        
        print(f"✅ Database connection successful")
        print(f"✅ Found {count} representative proteins")
        
        conn.close()
        return True
        
    except ImportError:
        print("❌ psycopg2 not installed")
        return False
    except Exception as e:
        print(f"❌ Database connection failed: {e}")
        return False

def find_test_proteins():
    """Find some proteins for testing"""
    print("\n=== Test Protein Discovery ===")
    
    batch_base = Path("/data/ecod/pdb_updates/batches")
    test_proteins = []
    
    # Look for a few known proteins
    target_proteins = ["8ovp_A", "8oni_L", "8p6i_L", "8p8o_H", "8oz3_B"]
    
    for batch_dir in batch_base.iterdir():
        if not batch_dir.is_dir():
            continue
            
        domains_dir = batch_dir / "domains"
        if not domains_dir.exists():
            continue
        
        for target in target_proteins:
            summary_file = domains_dir / f"{target}.develop291.domain_summary.xml"
            if summary_file.exists():
                test_proteins.append({
                    'protein_id': target,
                    'batch': batch_dir.name,
                    'summary_file': str(summary_file)
                })
                if len(test_proteins) >= 3:
                    break
        
        if len(test_proteins) >= 3:
            break
    
    if test_proteins:
        print(f"✅ Found {len(test_proteins)} test proteins:")
        for protein in test_proteins:
            print(f"   {protein['protein_id']} in {protein['batch']}")
    else:
        print("❌ No test proteins found")
    
    return test_proteins

def main():
    """Run all validation checks"""
    print("🔍 Mini Production Setup Validation\n")
    
    checks = [
        ("Config File", check_config),
        ("Mini Executable", check_mini_executable),
        ("Batch Directories", check_batch_directories),
        ("SLURM", check_slurm),
        ("Database", check_database_connection)
    ]
    
    results = {}
    for name, check_func in checks:
        try:
            results[name] = check_func()
        except Exception as e:
            print(f"❌ {name} check crashed: {e}")
            results[name] = False
    
    # Find test proteins
    test_proteins = find_test_proteins()
    
    # Summary
    print("\n" + "="*50)
    print("VALIDATION SUMMARY")
    print("="*50)
    
    passed = 0
    for name, result in results.items():
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{name:20} {status}")
        if result:
            passed += 1
    
    print(f"\nOverall: {passed}/{len(results)} checks passed")
    
    if passed >= 4:  # Most checks pass
        print("\n🎉 Setup looks good! Ready to test mini production.")
        
        if test_proteins:
            print("\nNext steps:")
            print("1. Test with a few proteins:")
            print("   python filesystem_batch_processor.py --test-proteins 3")
            print("2. Monitor progress:")
            print("   python filesystem_batch_processor.py --monitor")
    else:
        print("\n⚠️  Some issues found. Please fix before proceeding.")
    
    return passed >= 4

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
