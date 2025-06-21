#!/usr/bin/env python3
"""
Debug test failures for pyecod_mini

This script runs individual tests and captures detailed error messages
to help identify the root causes of the remaining failures.
"""

import subprocess
import sys
from pathlib import Path

def run_single_test(test_path, capture_output=True):
    """Run a single test and capture detailed output"""
    cmd = ["python", "-m", "pytest", test_path, "-v", "-s", "--tb=long"]
    
    print(f"\n{'='*70}")
    print(f"🔍 DEBUGGING: {test_path}")
    print(f"{'='*70}")
    print(f"Command: {' '.join(cmd)}")
    print("-" * 70)
    
    try:
        result = subprocess.run(cmd, capture_output=capture_output, text=True, timeout=60)
        
        if result.stdout:
            print("STDOUT:")
            print(result.stdout)
        
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        
        print(f"\nExit code: {result.returncode}")
        return result.returncode == 0, result.stdout, result.stderr
        
    except subprocess.TimeoutExpired:
        print("❌ Test timed out")
        return False, "", "Timeout"
    except Exception as e:
        print(f"❌ Error running test: {e}")
        return False, "", str(e)

def check_imports():
    """Check if we can import the required modules"""
    print("🔍 CHECKING IMPORTS")
    print("="*50)
    
    imports_to_check = [
        ("mini.core.writer", "write_domain_partition"),
        ("mini.core.models", "Domain"),
        ("mini.core.models", "Evidence"), 
        ("mini.core.models", "PartitionMetadata"),
        ("mini.core.partitioner", "partition_domains"),
        ("mini.core.sequence_range", "SequenceRange"),
    ]
    
    for module_name, item_name in imports_to_check:
        try:
            module = __import__(module_name, fromlist=[item_name])
            item = getattr(module, item_name)
            print(f"✅ {module_name}.{item_name}")
            
            # For classes, check their signature
            if item_name in ["PartitionMetadata", "Domain", "Evidence"]:
                import inspect
                sig = inspect.signature(item.__init__)
                print(f"   Signature: {sig}")
                
        except ImportError as e:
            print(f"❌ {module_name}.{item_name} - ImportError: {e}")
        except AttributeError as e:
            print(f"❌ {module_name}.{item_name} - AttributeError: {e}")
        except Exception as e:
            print(f"⚠️  {module_name}.{item_name} - Other error: {e}")

def inspect_writer_api():
    """Inspect the current writer API to understand the signature"""
    print("\n🔍 INSPECTING WRITER API")
    print("="*50)
    
    try:
        from mini.core.writer import write_domain_partition
        import inspect
        
        sig = inspect.signature(write_domain_partition)
        print(f"write_domain_partition signature: {sig}")
        
        # Get the docstring
        doc = inspect.getdoc(write_domain_partition)
        if doc:
            print(f"Docstring: {doc[:200]}...")
            
    except Exception as e:
        print(f"❌ Could not inspect writer API: {e}")

def inspect_partitioner_api():
    """Inspect the current partitioner API"""
    print("\n🔍 INSPECTING PARTITIONER API")
    print("="*50)
    
    try:
        from mini.core.partitioner import partition_domains
        import inspect
        
        sig = inspect.signature(partition_domains)
        print(f"partition_domains signature: {sig}")
        
        # Get the docstring
        doc = inspect.getdoc(partition_domains)
        if doc:
            print(f"Docstring: {doc[:200]}...")
            
    except Exception as e:
        print(f"❌ Could not inspect partitioner API: {e}")

def main():
    """Main debugging function"""
    print("🐛 PyECOD Mini Test Failure Debugger")
    print("="*50)
    
    # Check if we're in the right directory
    if not Path("tests").exists():
        print("❌ Please run from the mini/ directory")
        return 1
    
    # Check imports first
    check_imports()
    
    # Inspect APIs
    inspect_writer_api()
    inspect_partitioner_api()
    
    # Test individual failing tests
    test_cases = [
        "tests/test_writer.py::TestDomainWriter::test_write_basic_domain_partition",
        "tests/test_core.py::TestResidueBlocking::test_basic_residue_blocking",
        "tests/test_parser.py::TestDomainSummaryParsing::test_parse_valid_domain_summary",
    ]
    
    results = {}
    for test_case in test_cases:
        success, stdout, stderr = run_single_test(test_case)
        results[test_case] = {
            'success': success,
            'stdout': stdout,
            'stderr': stderr
        }
    
    # Summary
    print("\n" + "="*70)
    print("🔍 DEBUG SUMMARY")
    print("="*70)
    
    for test_case, result in results.items():
        status = "✅" if result['success'] else "❌"
        print(f"{status} {test_case}")
        if not result['success']:
            # Extract key error lines
            if result['stdout']:
                lines = result['stdout'].split('\n')
                error_lines = [line for line in lines if any(keyword in line.lower() 
                              for keyword in ['error', 'failed', 'assert', 'exception'])]
                if error_lines:
                    print(f"   Key errors: {error_lines[:3]}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
