#!/usr/bin/env python3
"""
Evidence Model Migration and Testing Script

This script helps test the new consolidated Evidence model against existing
evidence creation code and validates compatibility.
"""

import os
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_evidence_compatibility():
    """Test the new Evidence model against various input scenarios"""
    
    print("Testing Evidence Model Compatibility...")
    print("=" * 50)
    
    # Test 1: Dictionary-based creation (current DomainEvidence.from_dict equivalent)
    print("\n1. Testing dictionary-based creation...")
    
    blast_dict = {
        "type": "domain_blast",
        "source_id": "e4hluA1",
        "domain_id": "e4hluA1",
        "query_range": "5-256",
        "hit_range": "1-252",
        "evalue": 1e-50,
        "hsp_count": 1,
        "t_group": "2004.1.1",
        "h_group": "2004.1",
        "pdb_id": "4hlu",
        "chain_id": "A"
    }
    
    try:
        from ecod.models.pipeline.evidence import Evidence
        evidence = Evidence.from_dict(blast_dict)
        print(f"✓ Created evidence from dict: {evidence.type}, confidence: {evidence.confidence:.4f}")
        
        # Test backward compatibility - to_dict
        result_dict = evidence.to_dict()
        print(f"✓ Converted back to dict: {len(result_dict)} fields")
        
        # Verify key fields are preserved
        assert result_dict["type"] == blast_dict["type"]
        assert result_dict["source_id"] == blast_dict["source_id"]
        assert result_dict["query_range"] == blast_dict["query_range"]
        assert "range" in result_dict  # Backward compatibility alias
        print("✓ Dictionary compatibility verified")
        
    except Exception as e:
        print(f"✗ Dictionary test failed: {e}")
        return False
    
    # Test 2: HHSearch evidence creation
    print("\n2. Testing HHSearch evidence creation...")
    
    class MockHHSearchHit:
        def __init__(self):
            self.domain_id = "e7p7q02"
            self.hit_id = "hit_1"
            self.range = "254-524"
            self.hit_range = "1-270"
            self.probability = 99.9
            self.evalue = 1e-45
            self.score = 180.5
    
    try:
        hhsearch_hit = MockHHSearchHit()
        evidence = Evidence.from_hhsearch_hit(hhsearch_hit)
        print(f"✓ Created HHSearch evidence: confidence: {evidence.confidence:.4f}")
        
        # Verify confidence calculation from probability
        expected_confidence = 99.9 / 100.0
        assert abs(evidence.confidence - expected_confidence) < 0.001
        print("✓ Probability to confidence conversion correct")
        
    except Exception as e:
        print(f"✗ HHSearch test failed: {e}")
        return False
    
    # Test 3: BLAST evidence creation
    print("\n3. Testing BLAST evidence creation...")
    
    class MockBlastHit:
        def __init__(self):
            self.domain_id = "e1abc01"
            self.hit_id = "hit_2"
            self.range = "10-150"
            self.hit_range = "5-145"
            self.evalue = 1e-20
            self.hsp_count = 2
            self.pdb_id = "1abc"
            self.chain_id = "A"
            self.discontinuous = False
    
    try:
        blast_hit = MockBlastHit()
        evidence = Evidence.from_blast_hit(blast_hit, "domain_blast")
        print(f"✓ Created BLAST evidence: confidence: {evidence.confidence:.4f}")
        
        # Verify e-value to confidence conversion
        expected_confidence = 1.0 / (1.0 + 1e-20)
        assert abs(evidence.confidence - expected_confidence) < 0.001
        print("✓ E-value to confidence conversion correct")
        
    except Exception as e:
        print(f"✗ BLAST test failed: {e}")
        return False
    
    # Test 4: XML serialization/deserialization
    print("\n4. Testing XML serialization...")
    
    try:
        # Create evidence and serialize to XML
        evidence = Evidence.from_dict(blast_dict)
        xml_element = evidence.to_xml()
        
        # Verify XML structure
        assert xml_element.tag == "evidence"
        assert xml_element.get("type") == "domain_blast"
        assert xml_element.get("source_id") == "e4hluA1"
        
        # Test query_range as child element
        query_range_elem = xml_element.find("query_range")
        assert query_range_elem is not None
        assert query_range_elem.text == "5-256"
        
        print("✓ XML serialization successful")
        
        # Test deserialization
        evidence_from_xml = Evidence.from_xml(xml_element)
        assert evidence_from_xml.type == evidence.type
        assert evidence_from_xml.source_id == evidence.source_id
        assert evidence_from_xml.query_range == evidence.query_range
        assert evidence_from_xml.confidence == evidence.confidence
        
        print("✓ XML deserialization successful")
        
    except Exception as e:
        print(f"✗ XML test failed: {e}")
        return False
    
    # Test 5: Legacy compatibility aliases
    print("\n5. Testing backward compatibility aliases...")
    
    try:
        from ecod.models.pipeline.evidence import DomainEvidence, BlastEvidence, HHSearchEvidence
        
        # Test that aliases work
        domain_ev = DomainEvidence.from_dict(blast_dict)
        blast_ev = BlastEvidence.from_dict(blast_dict)
        hhsearch_ev = HHSearchEvidence.from_dict(blast_dict)
        
        assert isinstance(domain_ev, Evidence)
        assert isinstance(blast_ev, Evidence)
        assert isinstance(hhsearch_ev, Evidence)
        
        print("✓ Backward compatibility aliases working")
        
    except Exception as e:
        print(f"✗ Compatibility test failed: {e}")
        return False
    
    # Test 6: Multiple evidence consolidation
    print("\n6. Testing evidence consolidation...")
    
    try:
        evidence_dicts = [
            {
                "type": "domain_blast",
                "source_id": "e1abc01",
                "query_range": "10-50",
                "confidence": 0.8
            },
            {
                "type": "domain_blast", 
                "source_id": "e1abc01",  # Same source
                "query_range": "10-50",  # Same range
                "confidence": 0.9        # Higher confidence
            },
            {
                "type": "hhsearch",
                "source_id": "e1def01",
                "query_range": "60-120",
                "confidence": 0.7
            }
        ]
        
        evidence_list = Evidence.create_from_multiple(evidence_dicts, consolidate=True)
        assert len(evidence_list) == 2  # Should consolidate first two
        
        # Find the consolidated domain_blast evidence
        domain_blast_ev = next(ev for ev in evidence_list if ev.type == "domain_blast")
        assert domain_blast_ev.confidence == 0.9  # Should keep higher confidence
        
        print("✓ Evidence consolidation working")
        
    except Exception as e:
        print(f"✗ Consolidation test failed: {e}")
        return False
    
    print("\n" + "=" * 50)
    print("✓ All Evidence model tests passed!")
    print("\nNext steps:")
    print("1. Create the directory: mkdir -p ecod/models/pipeline")
    print("2. Save the Evidence model to: ecod/models/pipeline/evidence.py")
    print("3. Update imports in existing code gradually")
    print("4. Test with real data from your domain analysis pipeline")
    
    return True

def create_migration_imports():
    """Create temporary import compatibility layer"""
    
    migration_code = '''
# Temporary migration imports - add to models/__init__.py

# New consolidated imports
from ecod.models.pipeline.evidence import Evidence, DomainEvidence as NewDomainEvidence

# Backward compatibility - gradually replace these
try:
    from ecod.models.evidence import DomainEvidence as LegacyDomainEvidence
    # Use legacy for now, but warn about deprecation
    import warnings
    warnings.warn("models.evidence.DomainEvidence is deprecated. Use models.pipeline.evidence.Evidence", 
                  DeprecationWarning, stacklevel=2)
    DomainEvidence = LegacyDomainEvidence
except ImportError:
    # New model only
    DomainEvidence = NewDomainEvidence

# Update these imports to point to new models
__all__ = [
    # Priority models (new)
    'Evidence', 'DomainEvidence',
    
    # Existing models (keep for now)
    'Protein', 'ProteinSequence', 'ProteinStructure',
    'PDBChain', 'ChainSequence', 'PDBEntry',
    'DomainSequence', 'DomainDSSPDetail',
    'DomainRange', 'DomainRangeSegment', 'DomainClassification',
    'Batch', 'ProcessStatus', 'ProcessFile',
    'Job', 'JobItem', 'ECODVersion', 'ReferenceResource',
    'BlastHit', 'HHSearchHit', 'DomainSummaryModel', 'PipelineResult',
    'ProteinResult', 'ProteinProcessingResult', 'DomainPartitionResult'
]
'''
    
    print("Migration import code:")
    print("=" * 50)
    print(migration_code)
    print("=" * 50)

if __name__ == "__main__":
    print("Evidence Model Migration Script")
    print("Testing consolidated Evidence model...")
    
    # Test the new model
    success = test_evidence_compatibility()
    
    if success:
        print("\n\nGenerating migration import code...")
        create_migration_imports()
    else:
        print("\n\n✗ Tests failed. Please fix issues before proceeding.")
        sys.exit(1)
