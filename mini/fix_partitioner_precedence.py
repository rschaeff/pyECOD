#!/usr/bin/env python3
"""Quick fix for partitioner precedence"""

import fileinput
import sys

# Fix the precedence in partitioner.py
filename = "./partitioner.py"

# Read the file and find the precedence_map
with open(filename, 'r') as f:
    content = f.read()

# Replace the precedence map
old_map = """        precedence_map = {
            'chain_blast': 1,           # Highest - contains architecture info
            'chain_blast_decomposed': 1.5,  # Still very high
            'domain_blast': 2,          # Good but less comprehensive
            'hhsearch': 3              # Lowest precedence
        }"""

new_map = """        precedence_map = {
            'chain_blast_decomposed': 0.5,  # Highest - most specific, decomposed hits
            'chain_blast': 1,               # High - but prefer decomposed when available
            'domain_blast': 2,              # Good but less comprehensive
            'hhsearch': 3                   # Lowest precedence
        }"""

if old_map in content:
    content = content.replace(old_map, new_map)
    with open(filename, 'w') as f:
        f.write(content)
    print(f"Updated precedence in {filename}")
    print("chain_blast_decomposed now has highest priority (0.5)")
else:
    print("Could not find the precedence map to update")
    print("Please manually edit mini/partitioner.py and change:")
    print("  'chain_blast_decomposed': 1.5  ->  'chain_blast_decomposed': 0.5")
