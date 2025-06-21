#!/usr/bin/env python3
"""
Targeted fixes for specific test failures

Based on the most common issues seen in pyecod_mini tests.
"""

import os
import sys
from pathlib import Path

def check_partitionmetadata_exists():
    """Check if PartitionMetadata class exists and get its signature"""
    try:
        from mini.core.models import PartitionMetadata
        import inspect
        sig = inspect.signature(PartitionMetadata.__init__)
        print(f"✅ PartitionMetadata exists with signature: {sig}")
        return True, str(sig)
    except ImportError:
        print("❌ PartitionMetadata not found in mini.core.models")
        return False, None
    except Exception as e:
        print(f"❌ Error with PartitionMetadata: {e}")
        return False, None

def create_minimal_writer_test():
    """Create a minimal writer test to isolate the API issue"""
    content = '''#!/usr/bin/env python3
"""
Minimal writer test to debug API issues
"""

import pytest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

# Test basic imports first
def test_imports():
    """Test that we can import everything we need"""
    try:
        from mini.core.writer import write_domain_partition
        print("✅ Can import write_domain_partition")
        
        from mini.core.models import Domain
        print("✅ Can import Domain")
        
        from mini.core.sequence_range import SequenceRange  
        print("✅ Can import SequenceRange")
        
        # Try to import PartitionMetadata
        try:
            from mini.core.models import PartitionMetadata
            print("✅ Can import PartitionMetadata")
            
            # Check signature
            import inspect
            sig = inspect.signature(PartitionMetadata.__init__)
            print(f"PartitionMetadata signature: {sig}")
            
        except ImportError:
            print("❌ Cannot import PartitionMetadata")
            
            # Maybe it's called something else?
            import mini.core.models as models
            classes = [name for name in dir(models) if name[0].isupper()]
            print(f"Available classes in models: {classes}")
            
        # Check write_domain_partition signature
        import inspect
        sig = inspect.signature(write_domain_partition)
        print(f"write_domain_partition signature: {sig}")
        
    except Exception as e:
        print(f"Import error: {e}")
        raise

def test_basic_domain_creation():
    """Test creating a basic domain"""
    from mini.core.models import Domain
    from mini.core.sequence_range import SequenceRange
    
    domain = Domain(
        id="test_d1",
        range=SequenceRange.parse("1-100"), 
        family="test_family",
        evidence_count=1,
        source="test",
        evidence_items=[]
    )
    
    assert domain.id == "test_d1"
    assert str(domain.range) == "1-100"
    print("✅ Domain creation works")

if __name__ == "__main__":
    test_imports()
    test_basic_domain_creation()
'''
    
    test_file = Path("test_minimal_writer.py")
    test_file.write_text(content)
    print(f"✅ Created {test_file}")
    return test_file

def create_alternative_writer_test():
    """Create writer test with different approaches based on actual API"""
    content = '''#!/usr/bin/env python3
"""
Alternative writer test with API detection
"""

import pytest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_write_domain_partition_api_detection(tmp_path):
    """Test write_domain_partition with automatic API detection"""
    from mini.core.writer import write_domain_partition
    from mini.core.models import Domain
    from mini.core.sequence_range import SequenceRange
    import inspect
    
    # Create a test domain
    domain = Domain(
        id="d1",
        range=SequenceRange.parse("10-100"),
        family="test_family", 
        evidence_count=1,
        source="test",
        evidence_items=[]
    )
    
    output_file = tmp_path / "test_output.xml"
    
    # Check the actual signature
    sig = inspect.signature(write_domain_partition)
    params = list(sig.parameters.keys())
    print(f"write_domain_partition parameters: {params}")
    
    # Try different API patterns
    if len(params) >= 3 and 'metadata' in params[1].lower():
        # New API with metadata object
        print("Detected new API with metadata")
        try:
            from mini.core.models import PartitionMetadata
            metadata = PartitionMetadata(
                pdb_id="1abc",
                chain_id="A", 
                algorithm_version="test",
                reference="test",
                is_classified=True
            )
            write_domain_partition([domain], metadata, str(output_file))
        except ImportError:
            # Maybe metadata class has different name or structure
            print("PartitionMetadata not found, trying alternative...")
            # Create a simple object with required attributes
            class SimpleMetadata:
                def __init__(self):
                    self.pdb_id = "1abc"
                    self.chain_id = "A"
                    self.algorithm_version = "test"
                    self.reference = "test"
                    self.is_classified = True
            
            metadata = SimpleMetadata()
            write_domain_partition([domain], metadata, str(output_file))
            
    elif len(params) >= 4:
        # Old API with separate parameters
        print("Detected old API with separate parameters")
        write_domain_partition([domain], "1abc", "A", str(output_file))
    else:
        # Unknown API
        print(f"Unknown API pattern with parameters: {params}")
        raise ValueError(f"Cannot determine API pattern from parameters: {params}")
    
    # Verify file was created
    assert output_file.exists()
    print("✅ File created successfully")
    
    # Try to parse the XML
    import xml.etree.ElementTree as ET
    tree = ET.parse(output_file)
    root = tree.getroot()
    print(f"Root element: {root.tag}")
    print(f"Root attributes: {root.attrib}")

if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        test_write_domain_partition_api_detection(Path(tmp_dir))
'''
    
    test_file = Path("test_alternative_writer.py")
    test_file.write_text(content)
    print(f"✅ Created {test_file}")
    return test_file

def create_minimal_core_test():
    """Create minimal core test to debug quality threshold issues"""
    content = '''#!/usr/bin/env python3
"""
Minimal core test to debug quality threshold issues
"""

import pytest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_minimal_evidence_processing():
    """Test minimal evidence processing"""
    from mini.core.models import Evidence
    from mini.core.partitioner import partition_domains
    from mini.core.sequence_range import SequenceRange
    import inspect
    
    # Check partition_domains signature
    sig = inspect.signature(partition_domains)
    print(f"partition_domains signature: {sig}")
    
    # Create high-quality evidence that should definitely pass
    evidence = Evidence(
        type="domain_blast",
        source_pdb="test_protein",
        query_range=SequenceRange.parse("10-100"),
        confidence=1.0,  # Maximum confidence
        evalue=1e-100,   # Excellent e-value
        reference_length=91,
        reference_coverage=1.0,  # Perfect coverage
        domain_id="test_protein_A"
    )
    
    print(f"Created evidence: {evidence}")
    print(f"Evidence type: {evidence.type}")
    print(f"Evidence confidence: {evidence.confidence}")
    print(f"Evidence reference_coverage: {evidence.reference_coverage}")
    
    # Try to run partition_domains with different parameter sets
    try:
        # Basic call
        domains = partition_domains([evidence], sequence_length=200)
        print(f"✅ Basic call succeeded: {len(domains)} domains")
        return domains
        
    except Exception as e:
        print(f"❌ Basic call failed: {e}")
        
        # Try with verbose to see what's happening
        try:
            domains = partition_domains([evidence], sequence_length=200, verbose=True)
            print(f"✅ Verbose call succeeded: {len(domains)} domains")
            return domains
        except Exception as e2:
            print(f"❌ Verbose call also failed: {e2}")
            
            # Try with additional parameters that might be needed
            try:
                domains = partition_domains(
                    [evidence], 
                    sequence_length=200, 
                    domain_definitions={},
                    verbose=True
                )
                print(f"✅ Call with domain_definitions succeeded: {len(domains)} domains")
                return domains
            except Exception as e3:
                print(f"❌ All calls failed. Last error: {e3}")
                raise

if __name__ == "__main__":
    test_minimal_evidence_processing()
'''
    
    test_file = Path("test_minimal_core.py")
    test_file.write_text(content)
    print(f"✅ Created {test_file}")
    return test_file

def main():
    """Create targeted test files for debugging"""
    print("🎯 Creating Targeted Debug Tests")
    print("="*50)
    
    # Check current state
    has_metadata, metadata_sig = check_partitionmetadata_exists()
    
    # Create debug test files
    minimal_writer = create_minimal_writer_test()
    alternative_writer = create_alternative_writer_test() 
    minimal_core = create_minimal_core_test()
    
    print("\n📋 NEXT STEPS:")
    print("="*50)
    print("1. Run minimal tests to isolate issues:")
    print(f"   python {minimal_writer}")
    print(f"   python {alternative_writer}")
    print(f"   python {minimal_core}")
    print()
    print("2. Run with pytest:")
    print(f"   python -m pytest {minimal_writer} -v -s")
    print(f"   python -m pytest {alternative_writer} -v -s") 
    print(f"   python -m pytest {minimal_core} -v -s")
    print()
    print("3. Based on the output, we can create targeted fixes")

if __name__ == "__main__":
    main()
