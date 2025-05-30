#!/usr/bin/env python3
"""
Debug script to trace range calculation issues in domain partitioning.
Focus on where domain ranges are being incorrectly calculated or reported.
"""

import logging
from typing import Dict, Any, List, Tuple
import xml.etree.ElementTree as ET

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def debug_range_pipeline():
    """Debug the complete range processing pipeline"""
    
    print("🔍 DEBUGGING DOMAIN RANGE CALCULATION PIPELINE")
    print("="*60)
    
    # 1. XML Parsing Stage
    print("\n📄 STAGE 1: XML PARSING")
    debug_xml_range_extraction()
    
    # 2. Evidence Creation Stage  
    print("\n🧬 STAGE 2: EVIDENCE CREATION")
    debug_evidence_range_processing()
    
    # 3. Domain Boundary Detection Stage
    print("\n🎯 STAGE 3: DOMAIN BOUNDARY DETECTION") 
    debug_boundary_determination()
    
    # 4. Coverage Calculation Stage
    print("\n📊 STAGE 4: COVERAGE CALCULATION")
    debug_coverage_calculation()
    
    # 5. Final Range Reporting
    print("\n📋 STAGE 5: FINAL RANGE REPORTING")
    debug_final_range_reporting()

def debug_xml_range_extraction():
    """Debug range extraction from XML"""
    
    # Simulate realistic domain summary XML
    sample_xml = """<?xml version="1.0"?>
    <blast_summ_doc>
        <blast_summ pdb="3hhp" chain="A"/>
        <blast_run program="blastp">
            <hits>
                <hit domain_id="e3hhpA1" evalues="1e-50">
                    <query_reg>5-180</query_reg>
                    <hit_reg>1-175</hit_reg>
                </hit>
                <hit domain_id="e3hhpA2" evalues="1e-45">
                    <query_reg>175-350</query_reg>
                    <hit_reg>1-175</hit_reg>
                </hit>
                <hit domain_id="e3hhpA3" evalues="1e-40">
                    <query_reg>345-506</query_reg>
                    <hit_reg>1-160</hit_reg>
                </hit>
            </hits>
        </blast_run>
        <hh_run program="hhsearch">
            <hits>
                <hit domain_id="e3hhpA1" probability="95.0">
                    <query_reg>10-185</query_reg>
                    <hit_reg>5-180</hit_reg>
                </hit>
                <hit domain_id="e3hhpA2" probability="92.0">
                    <query_reg>180-355</query_reg>
                    <hit_reg>3-178</hit_reg>
                </hit>
                <hit domain_id="e3hhpA3" probability="88.0">
                    <query_reg>350-500</query_reg>
                    <hit_reg>1-150</hit_reg>
                </hit>
            </hits>
        </hh_run>
    </blast_summ_doc>"""
    
    root = ET.fromstring(sample_xml)
    
    print("🔍 Testing range extraction methods:")
    
    # Test current range extraction logic
    def extract_range_current(element, *tags):
        """Current range extraction from analyzer.py"""
        for tag in tags:
            child = element.find(tag)
            if child is not None and child.text:
                return child.text.strip()
            attr_value = element.get(tag, "")
            if attr_value:
                return attr_value
        return ""
    
    def parse_range_current(range_str):
        """Current range parsing from analyzer.py"""
        if not range_str or '-' not in range_str:
            return []
        
        ranges = []
        for segment in range_str.split(','):
            segment = segment.strip()
            if '-' in segment:
                try:
                    start, end = segment.split('-')
                    ranges.append((int(start), int(end)))
                except ValueError:
                    pass
        return ranges
    
    # Extract ranges from sample data
    blast_hits = root.find(".//blast_run/hits")
    if blast_hits is not None:
        for hit in blast_hits.findall("hit"):
            domain_id = hit.get("domain_id", "unknown")
            query_range = extract_range_current(hit, "query_reg", "query_range")
            hit_range = extract_range_current(hit, "hit_reg", "hit_range")
            
            query_parsed = parse_range_current(query_range)
            hit_parsed = parse_range_current(hit_range)
            
            print(f"  BLAST {domain_id}:")
            print(f"    Raw query_range: '{query_range}'")
            print(f"    Parsed query: {query_parsed}")
            print(f"    Query span: {query_parsed[0][1] - query_parsed[0][0] + 1 if query_parsed else 0} residues")
    
    hh_hits = root.find(".//hh_run/hits") 
    if hh_hits is not None:
        for hit in hh_hits.findall("hit"):
            domain_id = hit.get("domain_id", "unknown")
            query_range = extract_range_current(hit, "query_reg", "query_range")
            
            query_parsed = parse_range_current(query_range)
            
            print(f"  HHSearch {domain_id}:")
            print(f"    Raw query_range: '{query_range}'")
            print(f"    Parsed query: {query_parsed}")
            print(f"    Query span: {query_parsed[0][1] - query_parsed[0][0] + 1 if query_parsed else 0} residues")

def debug_evidence_range_processing():
    """Debug evidence creation and range handling"""
    
    print("🔍 Testing evidence range processing:")
    
    # Simulate evidence creation
    mock_blast_hit = {
        'domain_id': 'e3hhpA1', 
        'query_range': '5-180',
        'hit_range': '1-175',
        'evalues': '1e-50'
    }
    
    mock_hh_hit = {
        'domain_id': 'e3hhpA1',
        'query_range': '10-185', 
        'hit_range': '5-180',
        'probability': '95.0'
    }
    
    # Check what ranges get stored in Evidence objects
    print(f"  BLAST evidence query_range: '{mock_blast_hit['query_range']}'")
    print(f"  HHSearch evidence query_range: '{mock_hh_hit['query_range']}'")
    
    # Check range parsing
    blast_parsed = parse_range_string(mock_blast_hit['query_range'])
    hh_parsed = parse_range_string(mock_hh_hit['query_range'])
    
    print(f"  BLAST parsed: {blast_parsed} -> span: {blast_parsed[0][1] - blast_parsed[0][0] + 1 if blast_parsed else 0}")
    print(f"  HHSearch parsed: {hh_parsed} -> span: {hh_parsed[0][1] - hh_parsed[0][0] + 1 if hh_parsed else 0}")

def parse_range_string(range_str):
    """Helper to parse range string"""
    if not range_str or '-' not in range_str:
        return []
    
    ranges = []
    for segment in range_str.split(','):
        segment = segment.strip()
        if '-' in segment:
            try:
                start, end = segment.split('-')
                ranges.append((int(start), int(end)))
            except ValueError:
                pass
    return ranges

def debug_boundary_determination():
    """Debug how domain boundaries are determined from evidence"""
    
    print("🔍 Testing domain boundary determination:")
    
    # Simulate evidence with overlapping ranges
    evidence_ranges = [
        ('e3hhpA1_blast', [(5, 180)]),     # 176 residues
        ('e3hhpA1_hh', [(10, 185)]),       # 176 residues  
        ('e3hhpA2_blast', [(175, 350)]),   # 176 residues
        ('e3hhpA2_hh', [(180, 355)]),      # 176 residues
        ('e3hhpA3_blast', [(345, 506)]),   # 162 residues
        ('e3hhpA3_hh', [(350, 500)]),      # 151 residues
    ]
    
    print("  Input evidence ranges:")
    total_evidence_coverage = 0
    for name, ranges in evidence_ranges:
        span = sum(end - start + 1 for start, end in ranges)
        total_evidence_coverage += span
        print(f"    {name}: {ranges} -> {span} residues")
    
    print(f"  Total evidence coverage: {total_evidence_coverage} residues")
    print(f"  Average per evidence: {total_evidence_coverage / len(evidence_ranges):.1f} residues")
    
    # Simulate boundary determination logic
    # Group by domain ID
    domain_groups = {}
    for name, ranges in evidence_ranges:
        domain_id = name.split('_')[0]  # e3hhpA1, e3hhpA2, etc.
        if domain_id not in domain_groups:
            domain_groups[domain_id] = []
        domain_groups[domain_id].extend(ranges)
    
    print("  Grouped by domain:")
    predicted_domains = []
    for domain_id, all_ranges in domain_groups.items():
        # Find consensus boundaries
        starts = [r[0] for r in all_ranges]
        ends = [r[1] for r in all_ranges]
        
        consensus_start = min(starts)  # Or could be median, mode, etc.
        consensus_end = max(ends)      # Or could be median, mode, etc.
        consensus_span = consensus_end - consensus_start + 1
        
        print(f"    {domain_id}: ranges {all_ranges}")
        print(f"      -> consensus: {consensus_start}-{consensus_end} ({consensus_span} residues)")
        predicted_domains.append((consensus_start, consensus_end))
    
    return predicted_domains

def debug_coverage_calculation():
    """Debug coverage calculation logic"""
    
    print("🔍 Testing coverage calculation:")
    
    # Use domains from previous step
    predicted_domains = [
        (5, 185),    # Domain 1: 181 residues 
        (175, 355),  # Domain 2: 181 residues
        (345, 500)   # Domain 3: 156 residues  
    ]
    
    sequence_length = 506
    
    print(f"  Sequence length: {sequence_length}")
    print(f"  Predicted domains: {predicted_domains}")
    
    # Method 1: Simple sum (wrong if overlaps exist)
    simple_sum = sum(end - start + 1 for start, end in predicted_domains)
    simple_coverage = simple_sum / sequence_length
    
    print(f"  Method 1 (simple sum): {simple_sum} residues = {simple_coverage:.1%} coverage")
    
    # Method 2: Union of ranges (correct)
    covered_positions = set()
    for start, end in predicted_domains:
        covered_positions.update(range(start, end + 1))
    
    union_coverage = len(covered_positions) / sequence_length
    
    print(f"  Method 2 (union): {len(covered_positions)} residues = {union_coverage:.1%} coverage")
    
    # Check overlaps
    print("  Overlap analysis:")
    for i, (start1, end1) in enumerate(predicted_domains):
        for j, (start2, end2) in enumerate(predicted_domains[i+1:], i+1):
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            if overlap_start <= overlap_end:
                overlap_size = overlap_end - overlap_start + 1
                print(f"    Domain {i+1} and {j+1} overlap: {overlap_start}-{overlap_end} ({overlap_size} residues)")

def debug_final_range_reporting():
    """Debug final range reporting in DomainModel"""
    
    print("🔍 Testing final range reporting:")
    
    # Simulate what should happen vs what might be happening
    expected_domains = [
        {'id': '3hhp_A_d1', 'start': 5, 'end': 185, 'range': '5-185'},
        {'id': '3hhp_A_d2', 'start': 175, 'end': 355, 'range': '175-355'}, 
        {'id': '3hhp_A_d3', 'start': 345, 'end': 500, 'range': '345-500'}
    ]
    
    print("  Expected domains:")
    expected_coverage = 0
    for domain in expected_domains:
        size = domain['end'] - domain['start'] + 1
        expected_coverage += size
        print(f"    {domain['id']}: {domain['range']} ({size} residues)")
    
    print(f"  Expected total coverage: {expected_coverage} residues")
    
    # Check for potential issues
    print("\n  Potential issues to check:")
    print("    1. Are domain start/end being set correctly in DomainModel?")
    print("    2. Is the range string being calculated correctly?") 
    print("    3. Are discontinuous domains being handled properly?")
    print("    4. Is coverage calculation using the right method?")
    print("    5. Are overlapping domains being resolved incorrectly?")
    
    # Specific code areas to check
    print("\n  Code areas to examine:")
    print("    - DomainModel.range property calculation") 
    print("    - PartitionProcessor._process_continuous_evidence()")
    print("    - DomainCandidate.to_domain_model()")
    print("    - DomainPartitionResult.coverage calculation")
    print("    - Range parsing in analyzer._parse_range_comprehensive()")

if __name__ == "__main__":
    debug_range_pipeline()
