#!/usr/bin/env python3
"""
Simple Coverage Debugging Script

A standalone script to debug the coverage calculation issue without complex test infrastructure.
Run this directly to identify where the 98% -> 29.6% discrepancy occurs.
"""

import sys
import os
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

def create_test_xml(temp_dir):
    """Create test XML that should produce ~98% coverage"""
    xml_content = '''<?xml version="1.0" encoding="utf-8"?>
<blast_summ_doc>
    <blast_summ pdb="3hhp" chain="A"/>
    <blast_run program="blastp">
        <hits>
            <hit domain_id="e3hhpA1" evalues="1e-50" hsp_count="1">
                <query_reg>5-180</query_reg>
                <hit_reg>1-175</hit_reg>
            </hit>
            <hit domain_id="e3hhpA2" evalues="1e-45" hsp_count="1">
                <query_reg>175-350</query_reg>
                <hit_reg>1-175</hit_reg>
            </hit>
            <hit domain_id="e3hhpA3" evalues="1e-40" hsp_count="1">
                <query_reg>345-506</query_reg>
                <hit_reg>1-160</hit_reg>
            </hit>
        </hits>
    </blast_run>
    <hh_run program="hhsearch">
        <hits>
            <hit domain_id="e3hhpA1" probability="95.0" evalue="1e-25" score="80.0">
                <query_reg>10-185</query_reg>
                <hit_reg>5-180</hit_reg>
            </hit>
            <hit domain_id="e3hhpA2" probability="92.0" evalue="1e-22" score="75.0">
                <query_reg>180-355</query_reg>
                <hit_reg>3-178</hit_reg>
            </hit>
            <hit domain_id="e3hhpA3" probability="88.0" evalue="1e-20" score="70.0">
                <query_reg>350-500</query_reg>
                <hit_reg>1-150</hit_reg>
            </hit>
        </hits>
    </hh_run>
</blast_summ_doc>'''
    
    xml_file = Path(temp_dir) / "3hhp_A.summary.xml"
    with open(xml_file, 'w') as f:
        f.write(xml_content)
    
    return str(xml_file)

def test_manual_coverage_calculation():
    """Test manual coverage calculation to verify our understanding"""
    print("🧮 MANUAL COVERAGE CALCULATION TEST")
    print("="*50)
    
    sequence_length = 506
    expected_domains = [
        (5, 185),    # 181 residues
        (175, 355),  # 181 residues  
        (345, 500)   # 156 residues
    ]
    
    print(f"Sequence length: {sequence_length}")
    print(f"Expected domains: {expected_domains}")
    
    # Calculate coverage using union method (correct)
    covered_positions = set()
    for start, end in expected_domains:
        covered_positions.update(range(start, end + 1))
    
    coverage = len(covered_positions) / sequence_length
    print(f"Union method: {len(covered_positions)}/{sequence_length} = {coverage:.6f} ({coverage*100:.1f}%)")
    
    # Calculate coverage using simple sum method (incorrect but shows overlap)
    simple_sum = sum(end - start + 1 for start, end in expected_domains)
    simple_coverage = simple_sum / sequence_length
    print(f"Simple sum method: {simple_sum}/{sequence_length} = {simple_coverage:.6f} ({simple_coverage*100:.1f}%)")
    
    # Show overlaps
    print(f"\nOverlap analysis:")
    domain1_positions = set(range(5, 186))    # 5-185
    domain2_positions = set(range(175, 356))  # 175-355
    domain3_positions = set(range(345, 501))  # 345-500
    
    overlap12 = domain1_positions.intersection(domain2_positions)
    overlap23 = domain2_positions.intersection(domain3_positions)
    
    print(f"  Domain 1 & 2 overlap: {len(overlap12)} residues ({min(overlap12)}-{max(overlap12)})")
    print(f"  Domain 2 & 3 overlap: {len(overlap23)} residues ({min(overlap23)}-{max(overlap23)})")
    
    return coverage

def test_xml_parsing():
    """Test XML parsing stage"""
    print("\n📄 XML PARSING TEST")
    print("="*50)
    
    with tempfile.TemporaryDirectory() as temp_dir:
        xml_file = create_test_xml(temp_dir)
        print(f"Created test XML: {xml_file}")
        
        try:
            from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
            from ecod.pipelines.domain_analysis.partition.models import PartitionOptions
            
            options = PartitionOptions()
            analyzer = EvidenceAnalyzer(options)
            
            print("Parsing XML...")
            result = analyzer.parse_domain_summary(xml_file)
            
            if 'error' in result:
                print(f"❌ Parse error: {result['error']}")
                return None
            
            print(f"✅ Parse successful")
            
            # Check BLAST hits
            blast_hits = result.get('domain_blast_hits', [])
            print(f"Domain BLAST hits: {len(blast_hits)}")
            for i, hit in enumerate(blast_hits):
                domain_id = hit.get('domain_id', 'unknown')
                query_range = hit.get('query_range', 'unknown')
                print(f"  Hit {i+1}: {domain_id} -> {query_range}")
            
            # Check HHSearch hits  
            hh_hits = result.get('hhsearch_hits', [])
            print(f"HHSearch hits: {len(hh_hits)}")
            for i, hit in enumerate(hh_hits):
                domain_id = hit.get('domain_id', 'unknown')
                query_range = hit.get('query_range', 'unknown')
                probability = hit.get('probability', 'unknown')
                print(f"  Hit {i+1}: {domain_id} -> {query_range} (prob: {probability})")
            
            return result
            
        except Exception as e:
            print(f"❌ XML parsing failed: {str(e)}")
            import traceback
            print(f"Traceback:\n{traceback.format_exc()}")
            return None

def test_evidence_extraction(summary_data):
    """Test evidence extraction stage"""
    print("\n🧬 EVIDENCE EXTRACTION TEST")
    print("="*50)
    
    if not summary_data:
        print("❌ No summary data to process")
        return None
        
    try:
        from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
        from ecod.pipelines.domain_analysis.partition.models import PartitionOptions
        
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        
        print("Extracting evidence...")
        evidence_list = analyzer.extract_evidence_with_classification(summary_data)
        
        print(f"✅ Evidence extracted: {len(evidence_list)} items")
        
        # Group by domain
        evidence_by_domain = {}
        for evidence in evidence_list:
            domain_id = evidence.domain_id or evidence.source_id or 'unknown'
            if domain_id not in evidence_by_domain:
                evidence_by_domain[domain_id] = []
            evidence_by_domain[domain_id].append(evidence)
        
        print(f"Evidence grouped by domain: {len(evidence_by_domain)} groups")
        for domain_id, domain_evidence in evidence_by_domain.items():
            print(f"  {domain_id}: {len(domain_evidence)} evidence items")
            for evidence in domain_evidence:
                print(f"    {evidence.type}: range={evidence.query_range}, conf={evidence.confidence:.3f}")
        
        return evidence_list
        
    except Exception as e:
        print(f"❌ Evidence extraction failed: {str(e)}")
        import traceback
        print(f"Traceback:\n{traceback.format_exc()}")
        return None

def test_domain_processing(evidence_list):
    """Test domain processing stage"""
    print("\n🎯 DOMAIN PROCESSING TEST")
    print("="*50)
    
    if not evidence_list:
        print("❌ No evidence to process")
        return None
        
    try:
        from ecod.pipelines.domain_analysis.partition.processor import PartitionProcessor
        from ecod.pipelines.domain_analysis.partition.models import PartitionOptions, PartitionContext
        from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
        
        options = PartitionOptions()
        analyzer = EvidenceAnalyzer(options)
        processor = PartitionProcessor(options, analyzer)
        
        # Create processing context
        context = PartitionContext(
            pdb_id="3hhp",
            chain_id="A",
            reference="develop291",
            output_dir=Path("/tmp"),
            sequence_length=506
        )
        
        print("Processing evidence...")
        result = processor.process_evidence(evidence_list, context)
        
        print(f"✅ Processing completed")
        print(f"Success: {result.success}")
        print(f"Error: {result.error}")
        print(f"Domains found: {len(result.domains)}")
        print(f"Coverage: {result.coverage:.6f} ({result.coverage*100:.1f}%)")
        print(f"Residues assigned: {result.residues_assigned}")
        print(f"Sequence length: {result.sequence_length}")
        
        # Analyze domains
        for i, domain in enumerate(result.domains):
            print(f"  Domain {i+1}:")
            if hasattr(domain, 'range'):
                print(f"    Range: {domain.range}")
                print(f"    Size: {domain.size}")
                print(f"    Confidence: {domain.confidence:.3f}")
                print(f"    Source: {domain.source}")
            elif isinstance(domain, dict):
                print(f"    Dict domain: {domain}")
        
        return result
        
    except Exception as e:
        print(f"❌ Domain processing failed: {str(e)}")
        import traceback
        print(f"Traceback:\n{traceback.format_exc()}")
        return None

def test_service_integration():
    """Test full service integration"""
    print("\n🏭 SERVICE INTEGRATION TEST")
    print("="*50)
    
    with tempfile.TemporaryDirectory() as temp_dir:
        xml_file = create_test_xml(temp_dir)
        
        try:
            from ecod.core.context import ApplicationContext
            from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService
            from unittest.mock import Mock
            
            # Create minimal context
            context = Mock(spec=ApplicationContext)
            context.config_manager = Mock()
            context.config_manager.config = {
                'database': {
                    'host': 'localhost',
                    'port': 5432,
                    'database': 'ecod_test',
                    'user': 'test',
                    'password': 'test'
                },
                'reference': {
                    'current_version': 'develop291'
                }
            }
            context.config_manager.get_db_config.return_value = context.config_manager.config['database']
            context.config_manager.get_path.return_value = temp_dir
            
            print("Creating service...")
            service = DomainPartitionService(context)
            
            print("Running partition_protein...")
            result = service.partition_protein(
                pdb_id="3hhp",
                chain_id="A",
                summary_path=xml_file,
                output_dir=temp_dir
            )
            
            print(f"✅ Service completed")
            print(f"Success: {result.success}")
            print(f"Error: {result.error}")
            print(f"Domains found: {len(result.domains)}")
            print(f"Coverage: {result.coverage:.6f} ({result.coverage*100:.1f}%)")
            print(f"Sequence length: {result.sequence_length}")
            
            # Compare with expected
            expected_coverage = 0.98  # From debug script
            actual_coverage = result.coverage
            discrepancy = abs(expected_coverage - actual_coverage)
            
            print(f"\n📊 COMPARISON:")
            print(f"Expected coverage: {expected_coverage:.1%}")
            print(f"Actual coverage: {actual_coverage:.1%}")
            print(f"Discrepancy: {discrepancy:.1%}")
            
            if discrepancy > 0.1:  # More than 10% off
                print(f"🚨 SIGNIFICANT DISCREPANCY!")
            else:
                print(f"✅ Coverage within reasonable range")
            
            return result
            
        except Exception as e:
            print(f"❌ Service integration failed: {str(e)}")
            import traceback
            print(f"Traceback:\n{traceback.format_exc()}")
            return None

def main():
    """Run comprehensive coverage debugging"""
    print("🔍 COMPREHENSIVE COVERAGE DEBUGGING")
    print("="*60)
    print("This script tests each stage of the pipeline to identify")
    print("where the coverage calculation discrepancy occurs.\n")
    
    # Stage 1: Manual calculation (baseline)
    expected_coverage = test_manual_coverage_calculation()
    
    # Stage 2: XML parsing
    summary_data = test_xml_parsing()
    
    # Stage 3: Evidence extraction
    evidence_list = test_evidence_extraction(summary_data)
    
    # Stage 4: Domain processing
    processor_result = test_domain_processing(evidence_list)
    
    # Stage 5: Full service integration
    service_result = test_service_integration()
    
    # Final analysis
    print(f"\n📋 FINAL ANALYSIS")
    print("="*50)
    print(f"Expected coverage (manual): {expected_coverage:.1%}")
    
    if processor_result:
        print(f"Processor result coverage: {processor_result.coverage:.1%}")
        
    if service_result:
        print(f"Service result coverage: {service_result.coverage:.1%}")
    
    # Identify where the discrepancy occurs
    print(f"\n💡 DIAGNOSIS:")
    
    if not summary_data:
        print("❌ Issue in XML parsing stage")
    elif not evidence_list:
        print("❌ Issue in evidence extraction stage")
    elif not processor_result or not processor_result.success:
        print("❌ Issue in domain processing stage")
    elif not service_result or not service_result.success:
        print("❌ Issue in service integration stage")
    else:
        # All stages completed - check where coverage diverges
        if service_result.coverage < 0.5:
            print("🚨 Major coverage loss in service integration")
        elif processor_result.coverage < 0.5:
            print("🚨 Major coverage loss in domain processing")
        else:
            print("✅ Pipeline appears to be working correctly")
    
    print(f"\nFor detailed investigation, check the output above")
    print(f"to see exactly where the coverage calculation fails.")

if __name__ == "__main__":
    main()
