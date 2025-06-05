#!/usr/bin/env python3
"""
Parse ECOD ref_range_cache.txt file for mini_pyecod

This replaces the database dump approach with the authoritative range cache file.
The ref_range_cache contains the actual domain definitions used by ECOD version291.

Format: ID    domain_id    range    pdb_id    chain_id
Example: 001520984    e2ia4A1    A:110-209    2ia4    A
"""

import csv
import os
from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass
from pathlib import Path
from .ecod_domains_parser import parse_ecod_domains_file
from .sequence_range import SequenceRange

@dataclass
class RangeCacheEntry:
    """Entry from the range cache file"""
    cache_id: str
    domain_id: str
    range_spec: str
    pdb_id: str
    chain_id: str
    range: SequenceRange
    length: int

def parse_range_cache(cache_file: str, verbose: bool = False) -> Dict[str, RangeCacheEntry]:
    """
    Parse the ECOD ref_range_cache.txt file.
    
    Args:
        cache_file: Path to ecod.develop291.range_cache.txt
        verbose: Whether to print parsing information
        
    Returns:
        Dict mapping domain_id -> RangeCacheEntry
    """
    entries = {}
    invalid_count = 0
    
    print(f"Parsing range cache: {cache_file}")
    
    try:
        with open(cache_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) != 5:
                    if verbose:
                        print(f"Line {line_num}: Expected 5 columns, got {len(parts)}")
                    invalid_count += 1
                    continue
                    
                cache_id, domain_id, range_spec, pdb_id, chain_id = parts
                
                # Parse the range specification (e.g., "A:110-209" or "A:3-109,A:210-279")
                try:
                    # Remove chain prefix and parse range
                    range_parts = []
                    for segment in range_spec.split(','):
                        segment = segment.strip()
                        if ':' in segment:
                            # Format: "A:110-209"
                            chain_part, range_part = segment.split(':', 1)
                            range_parts.append(range_part)
                        else:
                            # Assume just the range part
                            range_parts.append(segment)
                    
                    # Join back for SequenceRange parsing
                    range_str = ','.join(range_parts)
                    sequence_range = SequenceRange.parse(range_str)
                    length = sequence_range.total_length
                    
                    entry = RangeCacheEntry(
                        cache_id=cache_id,
                        domain_id=domain_id,
                        range_spec=range_spec,
                        pdb_id=pdb_id.lower(),  # Normalize to lowercase
                        chain_id=chain_id,
                        range=sequence_range,
                        length=length
                    )
                    
                    entries[domain_id] = entry
                    
                except Exception as e:
                    if verbose:
                        print(f"Line {line_num}: Failed to parse range '{range_spec}': {e}")
                    invalid_count += 1
                    continue
    
    except FileNotFoundError:
        print(f"ERROR: Range cache file not found: {cache_file}")
        return {}
    except Exception as e:
        print(f"ERROR: Failed to parse range cache: {e}")
        return {}
    
    print(f"Loaded {len(entries)} domain entries from range cache")
    if invalid_count > 0:
        print(f"Skipped {invalid_count} invalid entries")
    
    return entries

def create_domain_definitions_from_cache_with_ecod(cache_file: str, domains_file: str, output_file: str, verbose: bool = False):
    """Create domain_definitions.csv with ECOD classifications"""

    # Parse range cache (has seqid ranges)
    entries = parse_range_cache(cache_file, verbose)

    # Parse ECOD domains file for classifications and pdb_ranges
    from mini.ecod_domains_parser import parse_ecod_domains_file
    ecod_domains = parse_ecod_domains_file(domains_file, verbose)

    print(f"Writing domain definitions with ECOD classifications to {output_file}")

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['domain_id', 'pdb_id', 'chain_id', 'range', 'length',
                        'pdb_range', 't_group', 'h_group', 'x_group'])

        for domain_id, entry in sorted(entries.items()):
            # Default values
            t_group = ""
            h_group = ""
            x_group = ""
            pdb_range = ""

            # Find matching ECOD domain
            key = (entry.pdb_id, entry.chain_id)
            if key in ecod_domains:
                ecod_class = ecod_domains[key]
                for ecod_domain in ecod_class.domains:
                    if ecod_domain.ecod_domain_id == domain_id:
                        t_group = ecod_domain.t_id
                        h_group = ecod_domain.h_name
                        x_group = ecod_domain.x_name
                        pdb_range = ecod_domain.pdb_range  # Preserve for visualization
                        break

            writer.writerow([
                domain_id,
                entry.pdb_id,
                entry.chain_id,
                str(entry.range),      # seqid range for computation
                entry.length,
                pdb_range,             # pdb range for visualization
                t_group,
                h_group,
                x_group
            ])

    print(f"Wrote {len(entries)} domain definitions with ECOD classifications")

def create_domain_definitions_from_cache(cache_file: str, output_file: str, verbose: bool = False):
    """
    Create domain_definitions.csv from range cache file for decomposition.
    
    Args:
        cache_file: Path to ecod.develop291.range_cache.txt
        output_file: Output CSV file path
        verbose: Whether to print detailed information
    """
    entries = parse_range_cache(cache_file, verbose)
    
    if not entries:
        print("No entries parsed from cache file!")
        return
    
    print(f"Writing domain definitions to {output_file}")
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['domain_id', 'pdb_id', 'chain_id', 'range', 'length'])
        
        for domain_id, entry in sorted(entries.items()):
            writer.writerow([
                domain_id,
                entry.pdb_id,
                entry.chain_id,
                str(entry.range),
                entry.length
            ])
    
    print(f"Wrote {len(entries)} domain definitions")

def extract_protein_lengths_from_cache(cache_file: str, output_file: str, verbose: bool = False):
    """
    Extract protein chain lengths by finding the maximum extent for each chain.
    
    This estimates the full protein length by looking at the maximum range end
    across all domains for each protein chain.
    
    Args:
        cache_file: Path to ecod.develop291.range_cache.txt
        output_file: Output CSV file path
        verbose: Whether to print detailed information
    """
    entries = parse_range_cache(cache_file, verbose)
    
    if not entries:
        print("No entries parsed from cache file!")
        return
    
    # Group by protein chain and find maximum extent
    chain_extents = {}
    
    for entry in entries.values():
        key = (entry.pdb_id, entry.chain_id)
        
        # Get maximum position from this domain
        max_pos = 0
        for segment in entry.range.segments:
            max_pos = max(max_pos, segment.end)
        
        if key not in chain_extents:
            chain_extents[key] = max_pos
        else:
            chain_extents[key] = max(chain_extents[key], max_pos)
    
    print(f"Writing protein lengths to {output_file}")
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pdb_id', 'chain_id', 'length'])
        
        for (pdb_id, chain_id), length in sorted(chain_extents.items()):
            writer.writerow([pdb_id, chain_id, length])
    
    print(f"Wrote {len(chain_extents)} protein chain lengths")

def validate_cache_data(cache_file: str, test_domains: List[str] = None):
    """
    Validate that the cache file contains expected test domains.
    
    Args:
        cache_file: Path to ecod.develop291.range_cache.txt
        test_domains: List of domain IDs to check for (default: 2ia4 domains)
    """
    if test_domains is None:
        test_domains = ['e2ia4A1', 'e2ia4A2', 'e2ia4B1', 'e2ia4B2']
    
    entries = parse_range_cache(cache_file)
    
    print(f"\nValidating cache data for test domains:")
    print("=" * 50)
    
    for domain_id in test_domains:
        if domain_id in entries:
            entry = entries[domain_id]
            print(f"✓ {domain_id}: {entry.range_spec} ({entry.length} residues)")
        else:
            print(f"✗ {domain_id}: NOT FOUND")
    
    # Check for 8ovp protein (should not be in cache since it's new)
    print(f"\nChecking for 8ovp domains (should be empty):")
    ovp_domains = [d for d in entries.keys() if '8ovp' in d.lower()]
    if ovp_domains:
        print(f"Found {len(ovp_domains)} 8ovp domains: {ovp_domains}")
    else:
        print("No 8ovp domains found (expected - it's a new structure)")

def main():
    """Main function with CLI interface"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Parse ECOD range cache file for mini_pyecod',
        epilog='Creates reference files from authoritative ECOD range cache'
    )
    
    parser.add_argument('--cache-file', 
                        default='/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt',
                        help='Path to ECOD range cache file')
    parser.add_argument('--output-dir', default='test_data',
                        help='Output directory for generated files')
    parser.add_argument('--validate', action='store_true',
                        help='Validate cache file for test domains')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.cache_file):
        print(f"ERROR: Cache file not found: {args.cache_file}")
        print("Please check the path or download the file from ECOD database")
        return
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.validate:
        validate_cache_data(args.cache_file)
        return
    
    # Generate all reference files
    print("Generating reference files from range cache...")
    print("=" * 60)
    
    # 1. Domain lengths (for domain BLAST/HHsearch coverage)
    domain_lengths_file = os.path.join(args.output_dir, 'domain_lengths.csv')
    create_domain_lengths_from_cache(args.cache_file, domain_lengths_file, args.verbose)
    
    # 2. Domain definitions (for chain BLAST decomposition)
    domain_defs_file = os.path.join(args.output_dir, 'domain_definitions.csv')
    create_domain_definitions_from_cache(args.cache_file, domain_defs_file, args.verbose)
    
    # 3. Protein lengths (for chain BLAST coverage)
    protein_lengths_file = os.path.join(args.output_dir, 'protein_lengths.csv')
    extract_protein_lengths_from_cache(args.cache_file, protein_lengths_file, args.verbose)
    
    print("\n" + "=" * 60)
    print("REFERENCE FILES CREATED:")
    print("=" * 60)
    print(f"Domain lengths:     {domain_lengths_file}")
    print(f"Domain definitions: {domain_defs_file}")
    print(f"Protein lengths:    {protein_lengths_file}")
    
    print("\nTo test with 8ovp_A:")
    print(f"python quick_test.py 8ovp_A --reference-lengths {domain_lengths_file}")

if __name__ == "__main__":
    main()
