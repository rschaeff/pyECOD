#!/usr/bin/env python3
"""Test what's happening with evidence reference lengths"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from mini.models import Evidence
from ecod.core.sequence_range import SequenceRange

print("Testing evidence reference length handling...")
print("=" * 60)

# Test 1: Create evidence manually and check attributes
print("\nTest 1: Manual evidence creation")
ev1 = Evidence(
    type="test",
    source_pdb="1abc",
    query_range=SequenceRange.parse("1-100"),
    reference_length=100
)
print(f"Created evidence with reference_length=100")
print(f"ev1.reference_length = {ev1.reference_length}")
print(f"hasattr(ev1, 'reference_length') = {hasattr(ev1, 'reference_length')}")

# Test 2: Check what parser produces
print("\n\nTest 2: Parser output")
xml_path = "/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424/domains/8ovp_A.develop291.domain_summary.xml"

# Load reference data
domain_lengths = load_reference_lengths("test_data/domain_lengths.csv") 
protein_lengths = load_protein_lengths("test_data/protein_lengths.csv")

print(f"Loaded {len(domain_lengths)} domain lengths")
print(f"Loaded {len(protein_lengths)} protein lengths")

# Parse with require_reference_lengths=False to see all evidence
evidence = parse_domain_summary(xml_path,
                               reference_lengths=domain_lengths,
                               protein_lengths=protein_lengths,
                               require_reference_lengths=False,
                               verbose=True)

print(f"\nTotal evidence: {len(evidence)}")

# Check reference lengths by type
by_type_with_ref = {}
by_type_total = {}

for ev in evidence:
    etype = ev.type
    by_type_total[etype] = by_type_total.get(etype, 0) + 1
    if ev.reference_length is not None:
        by_type_with_ref[etype] = by_type_with_ref.get(etype, 0) + 1

print("\nReference lengths by evidence type:")
for etype in sorted(by_type_total.keys()):
    with_ref = by_type_with_ref.get(etype, 0)
    total = by_type_total[etype]
    print(f"  {etype}: {with_ref}/{total} have reference lengths")

# Sample some evidence to see details
print("\nSample evidence details:")
for etype in ['chain_blast', 'domain_blast', 'hhsearch']:
    type_evidence = [e for e in evidence if e.type == etype][:2]
    if type_evidence:
        print(f"\n{etype}:")
        for ev in type_evidence:
            print(f"  {ev.source_pdb} / {ev.domain_id}")
            print(f"    reference_length: {ev.reference_length}")
            print(f"    hasattr check: {hasattr(ev, 'reference_length')}")

# Test 3: Check if reference_length survives a list
print("\n\nTest 3: Reference length in list")
test_list = [ev1]
print(f"ev1 in list: {test_list[0].reference_length}")

# Test 4: What the partitioner sees
print("\n\nTest 4: Simulating partitioner check")
# This is what the partitioner does
with_ref_length = sum(1 for e in evidence if e.reference_length is not None)
print(f"Partitioner check: {with_ref_length}/{len(evidence)} have reference lengths")

# Check if it's an attribute access issue
print("\nChecking attribute access methods:")
if evidence:
    ev = evidence[0]
    print(f"Direct access: {ev.reference_length}")
    print(f"getattr: {getattr(ev, 'reference_length', 'NOT FOUND')}")
    print(f"__dict__: {'reference_length' in ev.__dict__ if hasattr(ev, '__dict__') else 'No __dict__'}")
