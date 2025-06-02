#!/usr/bin/env python3
"""Find test cases from domain summary files that actually exist"""

import os
import xml.etree.ElementTree as ET
import random
from collections import defaultdict

def analyze_existing_files(batch_dir, sample_size=50):
    """Analyze a sample of existing domain summary files"""
    
    domains_dir = os.path.join(batch_dir, "domains")
    
    # Get all domain summary files
    all_files = [f for f in os.listdir(domains_dir) 
                 if f.endswith(".domain_summary.xml") and not f.endswith(".domain_summary.xml.domain_summary.xml")]
    
    print(f"Found {len(all_files)} domain summary files")
    
    # Sample files to analyze
    random.seed(42)  # For reproducibility
    sample_files = random.sample(all_files, min(sample_size, len(all_files)))
    
    # Categorize by characteristics
    single_domain = []
    two_domain = []
    three_domain = []
    many_domain = []
    discontinuous = []
    has_chain_blast = []
    
    print(f"\nAnalyzing {len(sample_files)} sample files...")
    
    for i, filename in enumerate(sample_files):
        if i % 10 == 0:
            print(f"  Processed {i}/{len(sample_files)}...")
            
        filepath = os.path.join(domains_dir, filename)
        protein = filename.replace(".develop291.domain_summary.xml", "")
        
        try:
            tree = ET.parse(filepath)
            root = tree.getroot()
            
            # Check for chain BLAST
            chain_blast_hits = root.findall(".//chain_blast_run/hits/hit")
            if chain_blast_hits:
                has_chain_blast.append(protein)
            
            # Check for domain BLAST to estimate domain count
            domain_blast_hits = root.findall(".//blast_run/hits/hit")
            unique_families = set()
            
            for hit in domain_blast_hits[:30]:  # Check first 30
                domain_id = hit.get("domain_id", "")
                if len(domain_id) > 5:
                    family = domain_id[1:5]
                    unique_families.add(family)
            
            # Check for discontinuous
            has_discontinuous = False
            for hit in chain_blast_hits[:5]:
                query_reg = hit.find("query_reg")
                if query_reg is not None and query_reg.text and ',' in query_reg.text:
                    has_discontinuous = True
                    discontinuous.append(protein)
                    break
            
            # Categorize by estimated domain count
            domain_count = len(unique_families)
            if domain_count <= 1:
                single_domain.append(protein)
            elif domain_count == 2:
                two_domain.append(protein)
            elif domain_count == 3:
                three_domain.append(protein)
            elif domain_count >= 4:
                many_domain.append((protein, domain_count))
                
        except Exception as e:
            continue
    
    # Report findings
    print("\n" + "="*60)
    print("TEST CASE RECOMMENDATIONS FROM EXISTING FILES:")
    print("="*60)
    
    # Single domain with chain BLAST
    single_with_blast = [p for p in single_domain if p in has_chain_blast]
    if single_with_blast:
        print(f"\n1. Single domain with chain BLAST ({len(single_with_blast)} found):")
        for p in single_with_blast[:3]:
            print(f"   python quick_test.py {p}")
    
    # Two domains with chain BLAST
    two_with_blast = [p for p in two_domain if p in has_chain_blast]
    if two_with_blast:
        print(f"\n2. Two domains with chain BLAST ({len(two_with_blast)} found):")
        for p in two_with_blast[:3]:
            print(f"   python quick_test.py {p}")
    
    # Three domains
    if three_domain:
        print(f"\n3. Three domain proteins ({len(three_domain)} found):")
        for p in three_domain[:3]:
            print(f"   python quick_test.py {p}")
    
    # Discontinuous
    if discontinuous:
        print(f"\n4. Discontinuous domains ({len(discontinuous)} found):")
        for p in discontinuous[:3]:
            print(f"   python quick_test.py {p}")
    
    # Many domains
    if many_domain:
        print(f"\n5. Many domain proteins ({len(many_domain)} found):")
        many_domain.sort(key=lambda x: x[1], reverse=True)
        for p, count in many_domain[:3]:
            print(f"   python quick_test.py {p}  # ~{count} domains")
    
    # Create a structured test list
    print("\n" + "="*60)
    print("RECOMMENDED TEST ORDER:")
    print("="*60)
    
    test_cases = []
    
    # Add one from each category
    if single_with_blast:
        test_cases.append((single_with_blast[0], "Single domain baseline"))
    if two_with_blast:
        test_cases.append((two_with_blast[0], "Two domain decomposition"))
    if three_domain:
        test_cases.append((three_domain[0], "Three domain complexity"))
    if discontinuous and discontinuous[0] != "8ovp_A":
        test_cases.append((discontinuous[0], "Discontinuous test"))
    if many_domain:
        test_cases.append((many_domain[0][0], f"Many domains ({many_domain[0][1]})"))
    
    for i, (case, desc) in enumerate(test_cases, 1):
        print(f"{i}. {case} - {desc}")
    
    return test_cases

def verify_case_files(batch_dir, protein):
    """Verify all required files exist for a test case"""
    
    domain_file = os.path.join(batch_dir, "domains", f"{protein}.develop291.domain_summary.xml")
    blast_file = os.path.join(batch_dir, "blast/chain", f"{protein}.develop291.xml")
    
    return os.path.exists(domain_file), os.path.exists(blast_file)

if __name__ == "__main__":
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    
    # Analyze more files for better coverage
    test_cases = analyze_existing_files(batch_dir, sample_size=100)
    
    print("\n" + "="*60)
    print("VERIFICATION:")
    print("="*60)
    
    # Verify the recommended cases
    for case, desc in test_cases[:5]:
        has_domain, has_blast = verify_case_files(batch_dir, case)
        status = "✓ Ready" if has_domain else "✗ Missing files"
        print(f"{case}: {status}")
        if has_domain and not has_blast:
            print(f"  Warning: No chain BLAST file")
