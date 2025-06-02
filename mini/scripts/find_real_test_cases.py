#!/usr/bin/env python3
"""Find test cases from domain summary files that actually exist"""

import os
import xml.etree.ElementTree as ET
import random
from collections import defaultdict
import re

def extract_protein_id(filename):
    """Extract protein ID from filename, handling various formats"""
    # Remove all extensions
    base = filename
    while True:
        new_base = base.replace('.domain_summary.xml', '')
        new_base = new_base.replace('.develop291', '')
        if new_base == base:
            break
        base = new_base
    return base

def analyze_existing_files(batch_dir, sample_size=50):
    """Analyze a sample of existing domain summary files"""

    domains_dir = os.path.join(batch_dir, "domains")

    # Get all domain summary files
    all_files = []
    for f in os.listdir(domains_dir):
        if '.domain_summary.xml' in f and os.path.isfile(os.path.join(domains_dir, f)):
            all_files.append(f)

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
        protein = extract_protein_id(filename)

        try:
            tree = ET.parse(filepath)
            root = tree.getroot()

            # Get sequence length if available
            seq_length = 0
            blast_summ = root.find(".//blast_summ")
            if blast_summ is not None:
                seq_length = int(blast_summ.get("length", "0"))

            # Check for chain BLAST
            chain_blast_hits = root.findall(".//chain_blast_run/hits/hit")
            if chain_blast_hits:
                has_chain_blast.append(protein)

            # Check for domain BLAST to estimate domain count
            domain_blast_hits = root.findall(".//blast_run/hits/hit")
            hhsearch_hits = root.findall(".//hh_run/hits/hit")

            # Collect unique protein families from hits
            unique_families = set()

            # From domain BLAST
            for hit in domain_blast_hits[:50]:  # Check first 50
                domain_id = hit.get("domain_id", "")
                if len(domain_id) > 5:
                    family = domain_id[1:5]
                    unique_families.add(family)

            # From HHSearch
            for hit in hhsearch_hits[:50]:
                hit_id = hit.get("hit_id", "") or hit.get("domain_id", "")
                if len(hit_id) > 5:
                    family = hit_id[1:5]
                    unique_families.add(family)

            # Check for discontinuous domains
            has_discontinuous = False
            all_hits = chain_blast_hits + domain_blast_hits + hhsearch_hits
            for hit in all_hits[:20]:
                query_reg = hit.find("query_reg")
                if query_reg is not None and query_reg.text and ',' in query_reg.text:
                    has_discontinuous = True
                    discontinuous.append(protein)
                    break

            # Better domain count estimation based on coverage
            domain_count = estimate_domain_count(root)

            # Categorize by estimated domain count
            if domain_count <= 1:
                single_domain.append(protein)
            elif domain_count == 2:
                two_domain.append(protein)
            elif domain_count == 3:
                three_domain.append(protein)
            elif domain_count >= 4:
                many_domain.append((protein, domain_count, seq_length))

        except Exception as e:
            print(f"    Error processing {filename}: {e}")
            continue

    # Report findings
    print("\n" + "="*60)
    print("TEST CASE RECOMMENDATIONS FROM EXISTING FILES:")
    print("="*60)

    # Single domain with chain BLAST
    single_with_blast = [p for p in single_domain if p in has_chain_blast]
    if single_with_blast:
        print(f"\n1. Single domain with chain BLAST ({len(single_with_blast)} found):")
        for p in single_with_blast[:5]:
            print(f"   {p}")

    # Two domains with chain BLAST
    two_with_blast = [p for p in two_domain if p in has_chain_blast]
    if two_with_blast:
        print(f"\n2. Two domains with chain BLAST ({len(two_with_blast)} found):")
        for p in two_with_blast[:5]:
            print(f"   {p}")

    # Three domains
    three_with_blast = [p for p in three_domain if p in has_chain_blast]
    if three_with_blast:
        print(f"\n3. Three domains with chain BLAST ({len(three_with_blast)} found):")
        for p in three_with_blast[:5]:
            print(f"   {p}")

    # Discontinuous
    if discontinuous:
        print(f"\n4. Discontinuous domains ({len(discontinuous)} found):")
        for p in discontinuous[:5]:
            print(f"   {p}")

    # Many domains
    if many_domain:
        print(f"\n5. Many domain proteins ({len(many_domain)} found):")
        many_domain.sort(key=lambda x: x[1], reverse=True)
        for p, count, length in many_domain[:5]:
            print(f"   {p}  # ~{count} domains, {length} residues")

    # Create test suite
    test_suite = create_test_suite(single_with_blast, two_with_blast,
                                   three_with_blast, discontinuous, many_domain)

    return test_suite

def estimate_domain_count(root):
    """Better domain count estimation based on range coverage"""
    ranges = []

    # Collect all query ranges from high-confidence hits
    for hit in root.findall(".//blast_run/hits/hit")[:30]:
        evalue = float(hit.get("evalues", "999"))
        if evalue < 1e-5:
            query_reg = hit.find("query_reg")
            if query_reg is not None and query_reg.text:
                ranges.append(parse_range(query_reg.text))

    for hit in root.findall(".//hh_run/hits/hit")[:30]:
        prob = float(hit.get("probability", "0"))
        if prob > 70:
            query_reg = hit.find("query_reg")
            if query_reg is not None and query_reg.text:
                ranges.append(parse_range(query_reg.text))

    # Merge overlapping ranges to estimate domain count
    if not ranges:
        return 0

    merged = merge_ranges(ranges)
    return len(merged)

def parse_range(range_str):
    """Parse range string like '10-50' or '10-50,60-80'"""
    segments = []
    for segment in range_str.split(','):
        parts = segment.strip().split('-')
        if len(parts) == 2:
            try:
                start = int(parts[0])
                end = int(parts[1])
                segments.append((start, end))
            except:
                pass
    return segments

def merge_ranges(ranges):
    """Merge overlapping ranges"""
    if not ranges:
        return []

    # Flatten all segments
    all_segments = []
    for r in ranges:
        all_segments.extend(r)

    if not all_segments:
        return []

    # Sort by start position
    all_segments.sort()

    # Merge overlapping
    merged = [all_segments[0]]
    for start, end in all_segments[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 20:  # Allow 20 residue gap
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged

def create_test_suite(single, two, three, discontinuous, many):
    """Create a balanced test suite"""
    test_suite = []

    # Core test cases
    if single:
        test_suite.append({
            'protein': single[0],
            'category': 'single_domain',
            'description': 'Single domain protein'
        })

    if two:
        # Try to get different types of two-domain proteins
        test_suite.append({
            'protein': two[0],
            'category': 'two_domain',
            'description': 'Two domain protein'
        })
        if len(two) > 3:
            test_suite.append({
                'protein': two[3],
                'category': 'two_domain_alt',
                'description': 'Alternative two domain protein'
            })

    if three:
        test_suite.append({
            'protein': three[0],
            'category': 'three_domain',
            'description': 'Three domain protein'
        })

    # Complex cases
    if discontinuous:
        for d in discontinuous[:2]:
            if d not in [t['protein'] for t in test_suite]:
                test_suite.append({
                    'protein': d,
                    'category': 'discontinuous',
                    'description': 'Discontinuous domain'
                })
                break

    if many:
        # Get a moderate and a complex case
        for p, count, length in many:
            if 4 <= count <= 6:
                test_suite.append({
                    'protein': p,
                    'category': 'moderate_multi',
                    'description': f'{count} domain protein ({length} residues)'
                })
                break

        for p, count, length in many:
            if count >= 10:
                test_suite.append({
                    'protein': p,
                    'category': 'complex_multi',
                    'description': f'{count} domain protein ({length} residues)'
                })
                break

    return test_suite

def verify_test_suite(batch_dir, test_suite):
    """Verify all required files exist for test cases"""
    print("\n" + "="*60)
    print("TEST SUITE VERIFICATION:")
    print("="*60)

    verified_suite = []

    for i, test in enumerate(test_suite, 1):
        protein = test['protein']

        # Check required files
        domain_file = os.path.join(batch_dir, "domains", f"{protein}.develop291.domain_summary.xml")
        blast_chain_file = os.path.join(batch_dir, "blast/chain", f"{protein}.develop291.xml")

        files_exist = {
            'domain_summary': os.path.exists(domain_file),
            'chain_blast': os.path.exists(blast_chain_file)
        }

        if files_exist['domain_summary']:
            print(f"{i}. {protein} - {test['description']}")
            print(f"   Category: {test['category']}")
            print(f"   Files: {'✓' if files_exist['chain_blast'] else '⚠ No chain BLAST'}")
            verified_suite.append(test)
        else:
            print(f"{i}. {protein} - SKIPPED (missing files)")

    return verified_suite

def write_test_script(test_suite, output_file="run_test_suite.sh"):
    """Write a shell script to run all tests"""
    with open(output_file, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# pyECOD mini test suite\n\n")

        f.write("echo '============================================================'\n")
        f.write("echo 'Running pyECOD mini test suite'\n")
        f.write("echo '============================================================'\n\n")

        for test in test_suite:
            f.write(f"echo '\\n--- Testing {test['protein']} ({test['description']}) ---'\n")
            f.write(f"python quick_test.py {test['protein']}\n")
            f.write("echo ''\n")

    os.chmod(output_file, 0o755)
    print(f"\nTest script written to: {output_file}")

if __name__ == "__main__":
    batch_dir = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424"

    # Analyze more files for better coverage
    test_suite = analyze_existing_files(batch_dir, sample_size=200)

    # Verify and finalize test suite
    verified_suite = verify_test_suite(batch_dir, test_suite)

    # Write test script
    if verified_suite:
        write_test_script(verified_suite)

        print("\n" + "="*60)
        print("QUICK TEST COMMANDS:")
        print("="*60)
        for test in verified_suite[:5]:
            print(f"python quick_test.py {test['protein']}")
