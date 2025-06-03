#!/usr/bin/env python3
"""
Quick start script for mini_pyecod production testing

This script runs the essential setup and validation steps.
"""

import subprocess
import sys
import os
from pathlib import Path

def run_command(cmd, description="", check=True):
    """Run a command with error handling"""
    print(f"\n=== {description} ===" if description else "")
    print(f"Running: {' '.join(cmd)}")
    print("-" * 50)
    
    try:
        result = subprocess.run(cmd, check=check, capture_output=False)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed with exit code {e.returncode}")
        return False
    except Exception as e:
        print(f"❌ Error running command: {e}")
        return False

def main():
    """Quick start workflow"""
    
    print("MINI PYECOD PRODUCTION TESTING - QUICK START")
    print("=" * 60)
    
    steps = [
        {
            "cmd": ["python", "apply_deterministic_fix.py"],
            "description": "Apply Deterministic Sorting Fix",
            "required": True
        },
        {
            "cmd": ["python", "setup_production_testing.py", "--test-proteins-only", "-v"],
            "description": "Setup Test Environment",
            "required": True
        },
        {
            "cmd": ["python", "-m", "pytest", "--version"],
            "description": "Check pytest Installation",
            "required": True
        },
        {
            "cmd": ["python", "run_tests.py", "unit", "-v"],
            "description": "Run Unit Tests",
            "required": False
        },
        {
            "cmd": ["python", "run_tests.py", "primary", "-v"],
            "description": "Run Primary Test Case",
            "required": False
        }
    ]
    
    success_count = 0
    
    for i, step in enumerate(steps, 1):
        print(f"\n[{i}/{len(steps)}] {step['description']}")
        
        if run_command(step["cmd"], step["description"]):
            success_count += 1
            print("✅ Success")
        else:
            print("❌ Failed")
            if step["required"]:
                print("This step is required - please fix the issue before continuing")
                return 1
    
    print("\n" + "=" * 60)
    print(f"QUICK START COMPLETE: {success_count}/{len(steps)} steps successful")
    print("=" * 60)
    
    if success_count >= 3:  # At least the required steps
        print("\n✅ Mini PyECOD testing environment is ready!")
        print("\nNext steps:")
        print("1. Run regression tests: python run_production_tests.py --primary -v")
        print("2. View results: cat regression_output/regression_results.json")
        print("3. Add more test proteins to expand coverage")
        return 0
    else:
        print("\n❌ Setup incomplete - please fix issues above")
        return 1

if __name__ == "__main__":
    exit(main())
