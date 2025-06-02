#!/usr/bin/env python3
"""Check if test data files have the entries we need"""

import os
import csv

print("Checking test data files for 8ovp entries...")
print("=" * 60)

# Check domain_lengths.csv
domain_file = "test_data/domain_lengths.csv"
if os.path.exists(domain_file):
    print(f"\nChecking {domain_file}:")
    with open(domain_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        print(f"  Header: {header}")
        
        # Look for any 8ovp-related entries
        found_8ovp = []
        total = 0
        for row in reader:
            total += 1
            if len(row) >= 1 and '8ovp' in row[0].lower():
                found_8ovp.append(row)
        
        print(f"  Total entries: {total}")
        print(f"  8ovp-related entries: {len(found_8ovp)}")
        for entry in found_8ovp:
            print(f"    {entry}")
            
        # Check for expected domains
        print("\n  Checking for expected domain families:")
        expected_domains = ['2ia4', '1gfl', '1ytf', '3mfh', '6dgv']
        f.seek(0)
        next(f)  # Skip header
        for domain in expected_domains:
            f.seek(0)
            next(f)  # Skip header
            found = False
            for row in csv.reader(f):
                if len(row) >= 1 and domain in row[0]:
                    print(f"    ✓ Found {domain}: {row}")
                    found = True
                    break
            if not found:
                print(f"    ✗ Missing {domain}")
else:
    print(f"ERROR: {domain_file} not found!")

# Check protein_lengths.csv
protein_file = "test_data/protein_lengths.csv"
if os.path.exists(protein_file):
    print(f"\n\nChecking {protein_file}:")
    with open(protein_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        print(f"  Header: {header}")
        
        # Count entries
        f.seek(0)
        next(f)  # Skip header
        total = sum(1 for _ in reader)
        print(f"  Total entries: {total}")
        
        # Look for 8ovp
        f.seek(0)
        next(f)  # Skip header
        for row in csv.reader(f):
            if len(row) >= 2 and '8ovp' in str(row[0]).lower():
                print(f"  Found 8ovp: {row}")
else:
    print(f"\nWARNING: {protein_file} not found!")
    print("This file is REQUIRED for chain BLAST evidence!")
    print("\nTo create it, you need a CSV with format:")
    print("  pdb_id,chain_id,length")
    print("  8ovp,A,518")
    print("  OR")
    print("  protein_id,length")  
    print("  8ovp_A,518")

# Show what we need
print("\n" + "=" * 60)
print("What we need for 8ovp_A domain partitioning:")
print("=" * 60)
print("\n1. In protein_lengths.csv (for chain BLAST):")
print("   - Entry for 8ovp_A or (8ovp, A) with length ~518")
print("\n2. In domain_lengths.csv (for domain BLAST/HHSearch):")
print("   - Entries for domains that hit 8ovp")
print("   - e.g., e2ia4A1, e1ytfA1, e3mfhA1, etc.")
print("\n3. Parser must use:")
print("   - protein_lengths for chain BLAST evidence")
print("   - domain_lengths for domain BLAST and HHSearch evidence")

# Create example protein_lengths.csv if missing
if not os.path.exists(protein_file):
    print(f"\n\nCreating example {protein_file}...")
    os.makedirs("test_data", exist_ok=True)
    with open(protein_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pdb_id', 'chain_id', 'length'])
        writer.writerow(['8ovp', 'A', '518'])
        writer.writerow(['2ia4', 'A', '238'])
        writer.writerow(['1ytf', 'A', '97'])
        writer.writerow(['3mfh', 'A', '90'])
    print(f"Created {protein_file} with example entries")
