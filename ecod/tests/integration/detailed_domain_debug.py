#!/usr/bin/env python3
"""
Detailed Domain Processing Debug

This script digs deeper into why domain candidates aren't being created
from the evidence that was successfully extracted.
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

def debug_domain_processing_detailed():
    """Debug domain processing in detail to see where candidates are lost"""
    print("🔬 DETAILED DOMAIN PROCESSING DEBUG")
    print("="*60)
    
    with tempfile.TemporaryDirectory() as temp_dir:
        xml_file = create_test_xml(temp_dir)
        
        try:
            from ecod.pipelines.domain_analysis.partition.analyzer import EvidenceAnalyzer
            from ecod.pipelines.domain_analysis.partition.processor import PartitionProcessor
            from ecod.pipelines.domain_analysis.partition.models import PartitionOptions, PartitionContext
            
            # Create analyzer and extract evidence
            options = PartitionOptions()
            analyzer = EvidenceAnalyzer(options)
            
            print("Step 1: Parse XML and extract evidence...")
            summary_data = analyzer.parse_domain_summary(xml_file)
            evidence_list = analyzer.extract_evidence_with_classification(summary_data)
            
            print(f"  Evidence extracted: {len(evidence_list)} items")
            for i, ev in enumerate(evidence_list):
                print(f"    {i+1}: {ev.type}, range={ev.query_range}, conf={ev.confidence:.3f}")
            
            # Create processor
            processor = PartitionProcessor(options, analyzer)
            
            print(f"\nStep 2: Test evidence grouping...")
            evidence_groups = analyzer.group_evidence_by_position(evidence_list, window_size=50)
            
            print(f"  Evidence groups created: {len(evidence_groups)}")
            for group_key, group_evidence in evidence_groups.items():
                print(f"    Group {group_key}: {len(group_evidence)} evidence items")
                for ev in group_evidence:
                    print(f"      {ev.type}: {ev.query_range}")
            
            print(f"\nStep 3: Test continuous evidence processing...")
            
            # Call the processor's continuous evidence processing directly
            try:
                candidates = processor._process_continuous_evidence(evidence_list, 506)
                print(f"  Domain candidates created: {len(candidates)}")
                
                for i, candidate in enumerate(candidates):
                    print(f"    Candidate {i+1}: {candidate.start}-{candidate.end}, conf={candidate.confidence:.3f}")
                    
            except Exception as e:
                print(f"  ❌ Error in continuous evidence processing: {str(e)}")
                import traceback
                print(f"  Traceback:\n{traceback.format_exc()}")
                
                # Try to debug the error
                print(f"\nStep 3b: Debug the error...")
                debug_continuous_evidence_processing(processor, evidence_list)
            
        except Exception as e:
            print(f"❌ Debug failed: {str(e)}")
            import traceback
            print(f"Traceback:\n{traceback.format_exc()}")

def debug_continuous_evidence_processing(processor, evidence_list):
    """Debug the continuous evidence processing method step by step"""
    print("🔍 DEBUGGING CONTINUOUS EVIDENCE PROCESSING")
    
    try:
        # Manually step through the _process_continuous_evidence logic
        sequence_length = 506
        
        print(f"  Input: {len(evidence_list)} evidence items, sequence_length={sequence_length}")
        
        # Step 1: Group evidence by position
        print("  Step 1: Grouping evidence by position...")
        evidence_groups = processor.analyzer.group_evidence_by_position(
            evidence_list,
            window_size=processor.options.merge_gap_tolerance * 2
        )
        
        print(f"    Groups created: {len(evidence_groups)}")
        for key, items in evidence_groups.items():
            print(f"      Group {key}: {len(items)} items")
        
        if not evidence_groups:
            print("    ❌ No evidence groups created!")
            return
        
        # Step 2: Process each group
        print("  Step 2: Processing each group...")
        candidates = []
        
        for position_key, evidence_items in evidence_groups.items():
            print(f"    Processing group {position_key} with {len(evidence_items)} items...")
            
            # Skip groups with insufficient evidence
            if len(evidence_items) < 1:
                print(f"      Skipped: insufficient evidence")
                continue
            
            # Create EvidenceGroup from list
            from ecod.pipelines.domain_analysis.partition.models import EvidenceGroup
            group = EvidenceGroup(evidence_items=evidence_items)
            
            # Get best evidence from group
            best_evidence = group.get_best_evidence()
            if not best_evidence:
                print(f"      Skipped: no best evidence found")
                continue
                
            print(f"      Best evidence: {best_evidence.type}, range={best_evidence.query_range}")
            
            # Check if evidence has coverage info
            print(f"      Evidence type: {type(best_evidence)}")
            print(f"      Has reference_coverage attr: {hasattr(best_evidence, 'reference_coverage')}")
            
            # Try to extract range
            if best_evidence.query_range:
                ranges = processor._parse_range(best_evidence.query_range)
                print(f"      Parsed range: {ranges}")
                
                if ranges[0] != 0 and ranges[1] != 0:  # Valid range
                    start, end = ranges
                    
                    # Create candidate
                    from ecod.pipelines.domain_analysis.partition.models import DomainCandidate
                    candidate = DomainCandidate(
                        start=start,
                        end=end,
                        evidence_group=group,
                        source=best_evidence.type,
                        confidence=best_evidence.confidence or 0.0
                    )
                    
                    # Check domain size
                    if candidate.size >= processor.options.min_domain_size:
                        if not processor.options.max_domain_size or candidate.size <= processor.options.max_domain_size:
                            candidates.append(candidate)
                            print(f"      ✅ Created candidate: {start}-{end}, size={candidate.size}")
                        else:
                            print(f"      Skipped: too large ({candidate.size} > {processor.options.max_domain_size})")
                    else:
                        print(f"      Skipped: too small ({candidate.size} < {processor.options.min_domain_size})")
                else:
                    print(f"      Skipped: invalid range {ranges}")
            else:
                print(f"      Skipped: no query range")
        
        print(f"  Final result: {len(candidates)} candidates created")
        
        # Try merging
        if candidates:
            print("  Step 3: Testing merge logic...")
            merged = processor._merge_nearby_candidates(candidates)
            print(f"    After merging: {len(merged)} candidates")
        
    except Exception as e:
        print(f"  ❌ Debug error: {str(e)}")
        import traceback
        print(f"  Traceback:\n{traceback.format_exc()}")

def test_evidence_group_creation():
    """Test EvidenceGroup creation directly"""
    print("\n🧪 TESTING EVIDENCE GROUP CREATION")
    print("="*50)
    
    try:
        from ecod.models.pipeline.evidence import Evidence
        from ecod.pipelines.domain_analysis.partition.models import EvidenceGroup
        
        # Create test evidence
        evidence1 = Evidence(
            type="hhsearch",
            source_id="e3hhpA1",
            domain_id="e3hhpA1", 
            query_range="10-185",
            probability=95.0,
            confidence=None
        )
        
        evidence2 = Evidence(
            type="hhsearch",
            source_id="e3hhpA2",
            domain_id="e3hhpA2",
            query_range="180-355", 
            probability=92.0,
            confidence=None
        )
        
        print(f"Evidence 1: {evidence1.type}, range={evidence1.query_range}, conf={evidence1.confidence}")
        print(f"Evidence 2: {evidence2.type}, range={evidence2.query_range}, conf={evidence2.confidence}")
        
        # Test EvidenceGroup creation
        group = EvidenceGroup(evidence_items=[evidence1, evidence2])
        print(f"Group created with {len(group.evidence_items)} items")
        
        # Test best evidence selection
        best = group.get_best_evidence()
        if best:
            print(f"Best evidence: {best.type}, range={best.query_range}, conf={best.confidence}")
        else:
            print("❌ No best evidence found")
            
    except Exception as e:
        print(f"❌ EvidenceGroup test failed: {str(e)}")
        import traceback
        print(f"Traceback:\n{traceback.format_exc()}")

def main():
    """Run detailed domain processing debug"""
    debug_domain_processing_detailed()
    test_evidence_group_creation()

if __name__ == "__main__":
    main()
