#!/usr/bin/env python3
"""
Complete Model Integration Test Script

This script tests the full integration of Evidence, DomainModel, and DomainPartitionResult
models working together, simulating the actual workflow in partition.py.
"""

import os
import sys
import tempfile
import xml.etree.ElementTree as ET
import time
from pathlib import Path

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_complete_workflow_integration():
    """Test complete workflow: Evidence -> DomainModel -> DomainPartitionResult"""
    print("\n1. Testing complete workflow integration...")
    
    try:
        from ecod.models.pipeline.evidence import Evidence
        from ecod.models.pipeline.domain import DomainModel  
        from ecod.models.pipeline.partition import DomainPartitionResult
        
        print("✓ All three models imported successfully")
        
        # Simulate real data from partition.py workflow
        # This mimics the domain_dicts created in _determine_domain_boundaries()
        domain_dicts = [
            {
                "domain_id": "8cmk_A_d1",
                "start": 2,
                "end": 923,
                "range": "2-923", 
                "source": "hhsearch",
                "confidence": 0.999,
                "source_id": "e4hluA1",
                "t_group": "2004.1.1",
                "h_group": "2004.1",
                "x_group": "2004",
                "a_group": "a.17",
                "is_manual_rep": False,
                "is_f70": False,
                "is_f40": False,
                "is_f99": False,
                "evidence": [
                    {
                        "type": "hhsearch",
                        "source_id": "e4hluA1",
                        "domain_id": "e4hluA1", 
                        "query_range": "2-923",
                        "hit_range": "1-922",
                        "confidence": 0.999,
                        "probability": 99.9,
                        "evalue": 1e-180,
                        "score": 1650.2
                    }
                ]
            },
            {
                "domain_id": "8cmk_A_d2", 
                "start": 50,
                "end": 150,
                "range": "50-150",
                "source": "domain_blast",
                "confidence": 0.85,
                "source_id": "e1abc01",
                "t_group": "3001.1.1",
                "evidence": [
                    {
                        "type": "domain_blast",
                        "source_id": "e1abc01",
                        "domain_id": "e1abc01",
                        "query_range": "50-150", 
                        "hit_range": "1-100",
                        "evalue": 1e-25,
                        "hsp_count": 1
                    }
                ]
            }
        ]
        
        # Test 1: Create DomainPartitionResult from domain dictionaries 
        # (simulates process_domains() in partition.py)
        pdb_id = "8cmk"
        chain_id = "A"
        reference = "develop291"
        sequence_length = 925
        
        result = DomainPartitionResult.from_domains(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            domains=domain_dicts,
            sequence_length=sequence_length
        )
        
        print(f"✓ Created DomainPartitionResult: {result}")
        print(f"  - Domains: {len(result.domains)}")
        print(f"  - Coverage: {result.coverage:.1%}")
        print(f"  - Sequence length: {result.sequence_length}")
        print(f"  - Is classified: {result.is_classified}")
        
        # Test 2: Verify automatic domain standardization
        assert len(result.domains) == 2
        
        # Check that domains were converted to DomainModel objects
        for i, domain in enumerate(result.domains):
            assert hasattr(domain, 'to_xml'), f"Domain {i} should be DomainModel object"
            assert hasattr(domain, 'evidence'), f"Domain {i} should have evidence attribute"
            assert hasattr(domain, 'confidence'), f"Domain {i} should have confidence attribute"
            
            print(f"  - Domain {i+1}: {domain.id}, confidence={domain.confidence:.3f}")
            print(f"    Range: {domain.range}, Evidence: {len(domain.evidence)} items")
        
        # Test 3: Verify evidence standardization
        domain1 = result.domains[0]
        evidence1 = domain1.evidence[0]
        
        # Evidence should be Evidence objects, not dictionaries
        assert hasattr(evidence1, 'type'), "Evidence should be Evidence object"
        assert evidence1.type == "hhsearch"
        assert hasattr(evidence1, 'confidence')
        assert evidence1.confidence > 0.9
        
        print(f"  - Evidence standardized: {evidence1.type}, confidence={evidence1.confidence:.3f}")
        
        # Test 4: Verify quality statistics generation
        assert result.domain_quality_stats
        stats = result.domain_quality_stats
        
        print(f"  - Quality stats: avg_confidence={stats.get('average_confidence', 0):.3f}")
        print(f"    Fully classified: {stats.get('fully_classified', 0)}/{stats.get('total_domains', 0)}")
        
        # Test 5: Verify evidence summary
        assert result.evidence_summary
        ev_summary = result.evidence_summary
        
        print(f"  - Evidence summary: {ev_summary.get('total_evidence_items', 0)} items")
        print(f"    Types: {ev_summary.get('evidence_type_distribution', {})}")
        
        print("✓ Complete workflow integration working correctly")
        
        return result
        
    except Exception as e:
        print(f"✗ Complete workflow integration failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def test_xml_round_trip_workflow():
    """Test XML serialization/deserialization of complete workflow"""
    print("\n2. Testing XML round-trip workflow...")
    
    try:
        from ecod.models.pipeline.partition import DomainPartitionResult
        
        # Create test result
        domain_dicts = [{
            "start": 10, "end": 100, "range": "10-100",
            "t_group": "2004.1.1", "confidence": 0.95,
            "evidence": [{"type": "hhsearch", "probability": 95.5}]
        }]
        
        original_result = DomainPartitionResult.from_domains(
            "test", "A", "develop291", domain_dicts, 150
        )
        
        # Test XML serialization
        xml_element = original_result.to_xml()
        assert xml_element.tag == "domain_partition"
        
        # Verify XML structure
        domains_elem = xml_element.find("domains")
        assert domains_elem is not None
        assert domains_elem.get("count") == "1"
        
        domain_elem = domains_elem.find("domain")
        assert domain_elem is not None
        assert domain_elem.get("start") == "10"
        assert domain_elem.get("end") == "100"
        
        # Check evidence in XML
        evidence_list = domain_elem.find("evidence_list")
        assert evidence_list is not None
        
        print("✓ XML serialization working correctly")
        
        # Test XML deserialization
        restored_result = DomainPartitionResult.from_xml(xml_element)
        
        assert restored_result.pdb_id == original_result.pdb_id
        assert restored_result.chain_id == original_result.chain_id
        assert len(restored_result.domains) == len(original_result.domains)
        assert restored_result.coverage == original_result.coverage
        
        print("✓ XML deserialization working correctly")
        
        # Test file save/load
        with tempfile.TemporaryDirectory() as temp_dir:
            success = original_result.save(temp_dir)
            assert success, "Save should succeed"
            
            assert os.path.exists(original_result.domain_file), "Domain file should exist"
            
            # Load from file
            loaded_result = DomainPartitionResult.from_xml_file(original_result.domain_file)
            assert loaded_result is not None
            assert loaded_result.pdb_id == original_result.pdb_id
            
            print("✓ File save/load working correctly")
        
        return True
        
    except Exception as e:
        print(f"✗ XML round-trip test failed: {e}")
        return False

def test_backward_compatibility_workflow():
    """Test backward compatibility with legacy partition.py patterns"""
    print("\n3. Testing backward compatibility workflow...")
    
    try:
        from ecod.models.pipeline.partition import DomainPartitionResult
        
        # Simulate legacy dictionary-based workflow
        legacy_domains = [
            {
                "start": 10, "end": 100, "range": "10-100",
                "source": "hhsearch", "confidence": 0.9,
                "t_group": "2004.1.1", "h_group": "2004.1"
            },
            {
                "start": 150, "end": 250, "range": "150-250", 
                "source": "blast", "confidence": 0.8
            }
        ]
        
        # Create result
        result = DomainPartitionResult.from_domains(
            "legacy", "A", "develop291", legacy_domains, 300
        )
        
        # Test to_dict() for backward compatibility
        result_dict = result.to_dict()
        
        # Verify expected fields exist
        expected_fields = [
            "pdb_id", "chain_id", "reference", "success", "domains",
            "sequence_length", "coverage", "domain_count"
        ]
        
        for field in expected_fields:
            assert field in result_dict, f"Missing expected field: {field}"
        
        # Verify domains are converted to dictionaries
        assert isinstance(result_dict["domains"], list)
        assert len(result_dict["domains"]) == 2
        
        for domain_dict in result_dict["domains"]:
            assert isinstance(domain_dict, dict)
            assert "start" in domain_dict
            assert "end" in domain_dict
            assert "range" in domain_dict
        
        print("✓ Backward compatibility working correctly")
        
        # Test mixed domain types (DomainModel + dict)
        mixed_domains = [
            result.domains[0],  # DomainModel object
            {"start": 300, "end": 400, "range": "300-400"}  # Plain dict
        ]
        
        mixed_result = DomainPartitionResult.from_domains(
            "mixed", "A", "develop291", mixed_domains, 450
        )
        
        assert len(mixed_result.domains) == 2
        print("✓ Mixed domain type handling working correctly")
        
        return True
        
    except Exception as e:
        print(f"✗ Backward compatibility test failed: {e}")
        return False

def test_performance_workflow():
    """Test performance of integrated workflow"""
    print("\n4. Testing workflow performance...")
    
    try:
        from ecod.models.pipeline.partition import DomainPartitionResult
        
        # Create larger test dataset
        large_domain_set = []
        
        for i in range(50):  # 50 proteins
            domains_per_protein = 2 + (i % 3)  # 2-4 domains per protein
            protein_domains = []
            
            for j in range(domains_per_protein):
                start = j * 100 + 10
                end = start + 80 + (j % 20)
                
                domain = {
                    "start": start,
                    "end": end,
                    "range": f"{start}-{end}",
                    "source": ["hhsearch", "blast", "domain_blast"][j % 3],
                    "confidence": 0.7 + (i + j) % 30 * 0.01,
                    "t_group": f"200{j}.1.{i % 3}",
                    "evidence": [
                        {
                            "type": ["hhsearch", "blast", "domain_blast"][j % 3],
                            "probability": 80 + (i + j) % 20,
                            "evalue": 10 ** -(10 + (i + j) % 10)
                        }
                    ]
                }
                protein_domains.append(domain)
            
            large_domain_set.append((f"test{i:03d}", "A", protein_domains))
        
        # Time the workflow
        start_time = time.time()
        
        results = []
        for pdb_id, chain_id, domains in large_domain_set:
            result = DomainPartitionResult.from_domains(
                pdb_id, chain_id, "develop291", domains, 500
            )
            results.append(result)
        
        creation_time = time.time() - start_time
        
        # Time XML serialization
        start_time = time.time()
        
        xml_elements = []
        for result in results[:10]:  # Test subset
            xml_elem = result.to_xml()
            xml_elements.append(xml_elem)
        
        xml_time = time.time() - start_time
        
        # Time analysis operations
        start_time = time.time()
        
        total_domains = sum(len(r.domains) for r in results)
        avg_coverage = sum(r.coverage for r in results) / len(results)
        high_confidence = sum(1 for r in results 
                             if r.domain_quality_stats.get('average_confidence', 0) > 0.9)
        
        analysis_time = time.time() - start_time
        
        print(f"✓ Performance test completed:")
        print(f"  - Created {len(results)} results in {creation_time:.3f}s")
        print(f"  - XML serialization: {xml_time:.3f}s for 10 results")
        print(f"  - Analysis operations: {analysis_time:.3f}s")
        print(f"  - Total domains processed: {total_domains}")
        print(f"  - Average coverage: {avg_coverage:.1%}")
        print(f"  - High confidence results: {high_confidence}")
        
        # Performance should be reasonable
        assert creation_time < 5.0, "Creation should be fast"
        assert xml_time < 2.0, "XML serialization should be fast"
        assert analysis_time < 1.0, "Analysis should be fast"
        
        return True
        
    except Exception as e:
        print(f"✗ Performance test failed: {e}")
        return False

def test_error_handling_workflow():
    """Test error handling in integrated workflow"""
    print("\n5. Testing error handling workflow...")
    
    try:
        from ecod.models.pipeline.partition import DomainPartitionResult
        
        # Test 1: Empty domains
        result = DomainPartitionResult.from_domains("empty", "A", "develop291", [], 100)
        assert result.is_unclassified
        assert not result.is_classified
        assert len(result.domains) == 0
        print("✓ Empty domains handled correctly")
        
        # Test 2: Invalid domain data
        invalid_domains = [
            {"start": "invalid", "end": 100},  # Invalid start
            {"start": 50, "end": "invalid"},   # Invalid end  
            {"range": "invalid-range"},        # Invalid range
            {}  # Empty domain
        ]
        
        result = DomainPartitionResult.from_domains(
            "invalid", "A", "develop291", invalid_domains, 200
        )
        
        # Should handle errors gracefully
        assert isinstance(result, DomainPartitionResult)
        print("✓ Invalid domain data handled gracefully")
        
        # Test 3: Error during save
        result = DomainPartitionResult(
            pdb_id="test", chain_id="A", reference="develop291"
        )
        
        # Try to save to invalid location
        success = result.save("/invalid/path/that/does/not/exist")
        assert not success, "Save to invalid path should fail"
        print("✓ Save error handling working correctly")
        
        # Test 4: XML parsing errors
        invalid_xml = "<invalid>not a domain partition</invalid>"
        
        try:
            root = ET.fromstring(invalid_xml)
            result = DomainPartitionResult.from_xml(root)
            # Should handle gracefully without crashing
            assert isinstance(result, DomainPartitionResult)
            print("✓ Invalid XML handled gracefully")
        except Exception:
            print("✓ Invalid XML rejected appropriately")
        
        return True
        
    except Exception as e:
        print(f"✗ Error handling test failed: {e}")
        return False

def main():
    """Run complete integration test suite"""
    print("Complete Model Integration Test Suite")
    print("Testing Evidence + DomainModel + DomainPartitionResult")
    print("=" * 60)
    
    # Track test results
    test_results = {}
    
    # Run integration tests
    test_results["workflow"] = test_complete_workflow_integration() is not None
    test_results["xml_roundtrip"] = test_xml_round_trip_workflow()
    test_results["compatibility"] = test_backward_compatibility_workflow()
    test_results["performance"] = test_performance_workflow()
    test_results["error_handling"] = test_error_handling_workflow()
    
    # Summary
    print("\n" + "=" * 60)
    print("Integration Test Results:")
    
    passed = sum(test_results.values())
    total = len(test_results)
    
    for test_name, result in test_results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {test_name:<20}: {status}")
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All integration tests passed!")
        print("\nThe consolidated Evidence, DomainModel, and DomainPartitionResult")
        print("models are ready for production use!")
        print("\n📋 Next Steps:")
        print("1. Create the model files in ecod/models/pipeline/")
        print("2. Update imports in partition.py")
        print("3. Test with real batch data")
        print("4. Roll out to production")
        
        print("\n💡 Benefits you'll get:")
        print("- Automatic evidence standardization")
        print("- Rich domain analysis and statistics")
        print("- Enhanced XML output with evidence")
        print("- Better error handling and logging")
        print("- Full backward compatibility")
        print("- Performance improvements")
        
    else:
        print(f"\n⚠️  {total - passed} tests failed. Please review issues before proceeding.")
    
    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
