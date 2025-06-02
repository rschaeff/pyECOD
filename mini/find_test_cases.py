#!/usr/bin/env python3
"""Find interesting test cases in the batch for mini_pyecod testing"""

import os
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

def analyze_batch_for_test_cases(batch_dir: str):
    """Analyze domain summary files to find interesting test cases"""
    
    domains_dir = os.path.join(batch_dir, "domains")
    
    if not os.path.exists(domains_dir):
        print(f"Domains directory not found: {domains_dir}")
        return
    
    # Categories of interest
    single_domain = []
    multi_domain = []
    discontinuous = []
    many_domains = []  # 4+ domains
    mixed_evidence = []  # Has multiple types of evidence
    
    # Get all domain summary files
    summary_files = [f for f in os.listdir(domains_dir) if f.endswith(".domain_summary.xml")]
    print(f"Found {len(summary_files)} domain summary files")
    
    for i, filename in enumerate(summary_files[:100]):  # Analyze first 100
        if i % 20 == 0:
            print(f"  Processed {i}/{min(100, len(summary_files))} files...")
            
        filepath = os.path.join(domains_dir, filename)
        
        try:
            tree = ET.parse(filepath)
            root = tree.getroot()
            
            # Extract PDB and chain
            parts = filename.split('.')
            pdb_chain = parts[0]
            if '_' not in pdb_chain:
                continue
                
            pdb_id, chain_id = pdb_chain.split('_', 1)
            
            # Count evidence types
            evidence_types = set()
            
            # Chain BLAST
            chain_blast_hits = root.findall(".//chain_blast_run/hits/hit")
            if chain_blast_hits:
                evidence_types.add('chain_blast')
                
            # Domain BLAST  
            domain_blast_hits = root.findall(".//blast_run/hits/hit")
            if domain_blast_hits:
                evidence_types.add('domain_blast')
                
            # HHsearch
            hhsearch_hits = root.findall(".//hh_run/hits/hit")
            if hhsearch_hits:
                evidence_types.add('hhsearch')
            
            # Skip if no evidence
            if not evidence_types:
                continue
            
            # Check for discontinuous ranges
            has_discontinuous = False
            domain_count_estimate = 0
            
            # Check chain BLAST for discontinuous
            for hit in chain_blast_hits[:5]:  # Check first few
                query_reg = hit.find("query_reg")
                if query_reg is not None and query_reg.text and ',' in query_reg.text:
                    has_discontinuous = True
                    # Estimate domains from discontinuous segments
                    domain_count_estimate = max(domain_count_estimate, len(query_reg.text.split(',')))
            
            # Estimate domain count from domain BLAST
            unique_families = set()
            for hit in domain_blast_hits[:20]:  # Check first 20
                domain_id = hit.get("domain_id", "")
                if len(domain_id) > 5:
                    family = domain_id[1:5]  # Extract PDB
                    unique_families.add(family)
            
            domain_count_estimate = max(domain_count_estimate, len(unique_families))
            
            # Categorize
            entry = {
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'evidence_types': evidence_types,
                'has_discontinuous': has_discontinuous,
                'domain_count_estimate': domain_count_estimate
            }
            
            if domain_count_estimate == 1:
                single_domain.append(entry)
            elif domain_count_estimate == 2 or domain_count_estimate == 3:
                multi_domain.append(entry)
            elif domain_count_estimate >= 4:
                many_domains.append(entry)
                
            if has_discontinuous:
                discontinuous.append(entry)
                
            if len(evidence_types) >= 2:
                mixed_evidence.append(entry)
                
        except Exception as e:
            continue
    
    # Report findings
    print(f"\n{'='*60}")
    print("INTERESTING TEST CASES FOUND:")
    print(f"{'='*60}")
    
    print(f"\n1. Single domain proteins: {len(single_domain)}")
    for entry in single_domain[:3]:
        print(f"   {entry['pdb_id']}_{entry['chain_id']} - evidence: {', '.join(entry['evidence_types'])}")
    
    print(f"\n2. Multi-domain proteins (2-3 domains): {len(multi_domain)}")
    for entry in multi_domain[:5]:
        print(f"   {entry['pdb_id']}_{entry['chain_id']} - domains: ~{entry['domain_count_estimate']}, evidence: {', '.join(entry['evidence_types'])}")
    
    print(f"\n3. Many-domain proteins (4+ domains): {len(many_domains)}")
    for entry in many_domains[:3]:
        print(f"   {entry['pdb_id']}_{entry['chain_id']} - domains: ~{entry['domain_count_estimate']}, evidence: {', '.join(entry['evidence_types'])}")
    
    print(f"\n4. Proteins with discontinuous domains: {len(discontinuous)}")
    for entry in discontinuous[:5]:
        print(f"   {entry['pdb_id']}_{entry['chain_id']} - domains: ~{entry['domain_count_estimate']}")
    
    print(f"\n5. Proteins with mixed evidence types: {len(mixed_evidence)}")
    for entry in mixed_evidence[:5]:
        print(f"   {entry['pdb_id']}_{entry['chain_id']} - evidence: {', '.join(entry['evidence_types'])}")
    
    # Suggest specific test cases
    print(f"\n{'='*60}")
    print("SUGGESTED TEST CASES:")
    print(f"{'='*60}")
    
    suggestions = []
    
    # Get one from each category
    if single_domain:
        s = single_domain[0]
        suggestions.append(f"{s['pdb_id']}_{s['chain_id']} - Simple single domain")
        
    if multi_domain:
        # Find one with chain_blast
        for m in multi_domain:
            if 'chain_blast' in m['evidence_types']:
                suggestions.append(f"{m['pdb_id']}_{m['chain_id']} - Multi-domain with chain BLAST")
                break
                
    if discontinuous:
        # Find one that isn't 8ovp
        for d in discontinuous:
            if d['pdb_id'] != '8ovp':
                suggestions.append(f"{d['pdb_id']}_{d['chain_id']} - Discontinuous domain")
                break
                
    if many_domains:
        suggestions.append(f"{many_domains[0]['pdb_id']}_{many_domains[0]['chain_id']} - Many domains (challenging)")
    
    for s in suggestions:
        print(f"  python test_framework.py {s}")

if __name__ == "__main__":
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"
    analyze_batch_for_test_cases(batch_dir)
