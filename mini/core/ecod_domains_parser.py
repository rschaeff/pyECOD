#!/usr/bin/env python3
"""
ECOD domains.txt file parser for production-quality validation

This module parses the ECOD domains file to extract T-group assignments
and domain boundaries for regression testing against the current pyecod engine.
"""

import csv
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from collections import defaultdict

from .sequence_range import SequenceRange

@dataclass
class EcodDomain:
    """ECOD domain entry from domains.txt"""
    uid: str
    ecod_domain_id: str
    manual_rep: str
    t_id: str              # T-group ID (e.g., "1.1.1")
    pdb: str
    chain: str
    pdb_range: str         # Original PDB range for visualization/literature (e.g., "A:-1-75,A:100-115")
    seqid_range: str       # Sequential range for computation (e.g., "1-77,102-117")
    unp_acc: str
    arch_name: str
    x_name: str
    h_name: str
    t_name: str
    f_name: str
    asm_status: str
    ligand: str

    # Computed fields
    sequence_range: Optional[SequenceRange] = None  # Parsed from seqid_range
    domain_length: int = 0

@dataclass
class EcodClassification:
    """Complete ECOD classification for a protein chain"""
    pdb_id: str
    chain_id: str
    domains: List[EcodDomain]
    
    @property
    def t_groups(self) -> Set[str]:
        """Get all T-group IDs for this chain"""
        return {d.t_id for d in self.domains}
    
    @property
    def t_names(self) -> Set[str]:
        """Get all T-group names for this chain"""
        return {d.t_name for d in self.domains}
    
    @property  
    def h_names(self) -> Set[str]:
        """Get all H-group names for this chain"""
        return {d.h_name for d in self.domains}
    
    @property
    def total_domains(self) -> int:
        """Total number of domains"""
        return len(self.domains)
    
    @property
    def total_coverage(self) -> int:
        """Total residues covered by domains"""
        return sum(d.domain_length for d in self.domains)

def parse_ecod_domains_file(domains_file: str, verbose: bool = False) -> Dict[Tuple[str, str], EcodClassification]:
    """Parse ECOD domains.txt file"""
    classifications = {}
    invalid_count = 0

    print(f"Parsing ECOD domains file: {domains_file}")

    try:
        with open(domains_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 16:
                    if verbose:
                        print(f"Line {line_num}: Expected 16 columns, got {len(parts)}")
                    invalid_count += 1
                    continue

                try:
                    # Parse ECOD domain entry
                    domain = EcodDomain(
                        uid=parts[0],
                        ecod_domain_id=parts[1],
                        manual_rep=parts[2],
                        t_id=parts[3],
                        pdb=parts[4].lower(),
                        chain=parts[5],
                        pdb_range=parts[6],      # Keep as string for visualization
                        seqid_range=parts[7],    # Use for computation
                        unp_acc=parts[8],
                        arch_name=parts[9],
                        x_name=parts[10],
                        h_name=parts[11],
                        t_name=parts[12],
                        f_name=parts[13],
                        asm_status=parts[14],
                        ligand=parts[15] if len(parts) > 15 else ""
                    )

                    # Parse sequence range from seqid_range (clean sequential numbering)
                    if domain.seqid_range and domain.seqid_range != 'NULL':
                        try:
                            domain.sequence_range = SequenceRange.parse(domain.seqid_range)
                            domain.domain_length = domain.sequence_range.total_length
                        except Exception as e:
                            if verbose:
                                print(f"Line {line_num}: Failed to parse seqid_range '{domain.seqid_range}': {e}")
                            invalid_count += 1
                            continue
                    else:
                        # Skip domains without valid seqid_range
                        invalid_count += 1
                        continue

                    # Group by (pdb_id, chain_id)
                    key = (domain.pdb, domain.chain)
                    if key not in classifications:
                        classifications[key] = EcodClassification(
                            pdb_id=domain.pdb,
                            chain_id=domain.chain,
                            domains=[]
                        )

                    classifications[key].domains.append(domain)

                except Exception as e:
                    if verbose:
                        print(f"Line {line_num}: Failed to parse domain entry: {e}")
                    invalid_count += 1
                    continue

    except FileNotFoundError:
        print(f"ERROR: ECOD domains file not found: {domains_file}")
        return {}
    except Exception as e:
        print(f"ERROR: Failed to parse ECOD domains file: {e}")
        return {}

    # Sort domains within each classification by sequence range start
    for classification in classifications.values():
        classification.domains.sort(key=lambda d: d.sequence_range.segments[0].start if d.sequence_range else 0)
    
    print(f"Loaded ECOD classifications for {len(classifications)} protein chains")
    if invalid_count > 0:
        print(f"Skipped {invalid_count} invalid entries")
    
    return classifications

def create_ecod_classifications_file(domains_file: str, output_file: str, 
                                   test_proteins: List[str] = None, verbose: bool = False):
    """
    Create a subset ECOD classifications file for testing
    
    Args:
        domains_file: Full ECOD domains.txt file
        output_file: Output CSV file for test data
        test_proteins: List of protein IDs to include (e.g., ["8ovp_A", "1abc_B"])
        verbose: Whether to print detailed information
    """
    classifications = parse_ecod_domains_file(domains_file, verbose)
    
    if not classifications:
        print("No classifications loaded from domains file!")
        return
    
    print(f"Creating ECOD classifications file: {output_file}")
    
    # Filter to test proteins if specified
    if test_proteins:
        test_keys = set()
        for protein_id in test_proteins:
            if '_' in protein_id:
                pdb_id, chain_id = protein_id.split('_', 1)
                test_keys.add((pdb_id.lower(), chain_id))
            else:
                # Include all chains for this PDB
                pdb_id = protein_id.lower()
                for key in classifications.keys():
                    if key[0] == pdb_id:
                        test_keys.add(key)
        
        filtered_classifications = {k: v for k, v in classifications.items() if k in test_keys}
        print(f"Filtered to {len(filtered_classifications)} chains for test proteins: {test_proteins}")
        classifications = filtered_classifications
    
    # Write CSV file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'pdb_id', 'chain_id', 'ecod_domain_id', 't_id', 't_name', 'h_name', 
            'arch_name', 'pdb_range', 'domain_length', 'manual_rep'
        ])
        
        count = 0
        for classification in classifications.values():
            for domain in classification.domains:
                writer.writerow([
                    domain.pdb,
                    domain.chain,
                    domain.ecod_domain_id,
                    domain.t_id,
                    domain.t_name,
                    domain.h_name,
                    domain.arch_name,
                    domain.pdb_range,
                    domain.domain_length,
                    domain.manual_rep
                ])
                count += 1
        
    print(f"Wrote {count} ECOD domain entries to {output_file}")

def load_ecod_classifications(classifications_file: str) -> Dict[Tuple[str, str], EcodClassification]:
    """
    Load ECOD classifications from CSV file
    
    Args:
        classifications_file: Path to CSV file created by create_ecod_classifications_file
        
    Returns:
        Dict mapping (pdb_id, chain_id) -> EcodClassification
    """
    classifications = {}
    
    if not os.path.exists(classifications_file):
        print(f"WARNING: ECOD classifications file not found: {classifications_file}")
        return classifications
    
    try:
        with open(classifications_file, 'r') as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                pdb_id = row['pdb_id'].lower()
                chain_id = row['chain_id']
                key = (pdb_id, chain_id)
                
                # Create EcodDomain object
                try:
                    range_str = row['pdb_range']
                    if ':' in range_str:
                        _, range_part = range_str.split(':', 1)
                    else:
                        range_part = range_str
                    
                    sequence_range = SequenceRange.parse(range_part)
                    domain_length = int(row['domain_length'])
                    
                except Exception as e:
                    print(f"Warning: Failed to parse range for {row['ecod_domain_id']}: {e}")
                    continue
                
                domain = EcodDomain(
                    uid="",
                    ecod_domain_id=row['ecod_domain_id'],
                    manual_rep=row['manual_rep'],
                    t_id=row['t_id'],
                    pdb=pdb_id,
                    chain=chain_id,
                    pdb_range=row['pdb_range'],
                    seqid_range="",
                    unp_acc="",
                    arch_name=row['arch_name'],
                    x_name="",
                    h_name=row['h_name'],
                    t_name=row['t_name'],
                    f_name="",
                    asm_status="",
                    ligand="",
                    sequence_range=sequence_range,
                    domain_length=domain_length
                )
                
                # Group by protein chain
                if key not in classifications:
                    classifications[key] = EcodClassification(
                        pdb_id=pdb_id,
                        chain_id=chain_id,
                        domains=[]
                    )
                
                classifications[key].domains.append(domain)
        
        # Sort domains within each classification
        for classification in classifications.values():
            classification.domains.sort(key=lambda d: d.sequence_range.segments[0].start if d.sequence_range else 0)
        
        print(f"Loaded ECOD classifications for {len(classifications)} protein chains")
        
    except Exception as e:
        print(f"ERROR: Failed to load ECOD classifications: {e}")
    
    return classifications

def get_ecod_summary(classifications: Dict[Tuple[str, str], EcodClassification]) -> Dict:
    """Get summary statistics for ECOD classifications"""
    
    if not classifications:
        return {
            'total_chains': 0,
            'total_domains': 0,
            'unique_t_groups': 0,
            'unique_h_groups': 0
        }
    
    total_domains = sum(c.total_domains for c in classifications.values())
    all_t_groups = set()
    all_h_groups = set()
    
    for classification in classifications.values():
        all_t_groups.update(classification.t_groups)
        all_h_groups.update(classification.h_names)
    
    return {
        'total_chains': len(classifications),
        'total_domains': total_domains,
        'unique_t_groups': len(all_t_groups),
        'unique_h_groups': len(all_h_groups),
        'avg_domains_per_chain': total_domains / len(classifications) if classifications else 0
    }

def validate_ecod_data_for_protein(protein_id: str, classifications: Dict[Tuple[str, str], EcodClassification]) -> bool:
    """
    Validate that we have ECOD data for a specific protein
    
    Args:
        protein_id: Protein ID (e.g., "8ovp_A")
        classifications: ECOD classifications dict
        
    Returns:
        True if ECOD data is available
    """
    if '_' in protein_id:
        pdb_id, chain_id = protein_id.split('_', 1)
        key = (pdb_id.lower(), chain_id)
        return key in classifications
    else:
        # Check if any chain for this PDB exists
        pdb_id = protein_id.lower()
        return any(key[0] == pdb_id for key in classifications.keys())

def main():
    """Command line interface for ECOD domains parsing"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Parse ECOD domains file for mini_pyecod testing',
        epilog='Creates ECOD classification data for production-quality validation'
    )
    
    parser.add_argument('--domains-file',
                        default='/data/ecod/database_versions/v291/ecod.develop291.domains.txt',
                        help='Path to ECOD domains.txt file')
    parser.add_argument('--output-dir', default='test_data',
                        help='Output directory for generated files')
    parser.add_argument('--test-proteins', nargs='+',
                        default=['8ovp_A', '1abc_A', '2ia4_A', '6dgv_A'],
                        help='Test proteins to include')
    parser.add_argument('--validate', action='store_true',
                        help='Validate domains file parsing')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.domains_file):
        print(f"ERROR: ECOD domains file not found: {args.domains_file}")
        print("Please check the path or download the file from ECOD database")
        return 1
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.validate:
        # Just validate parsing
        classifications = parse_ecod_domains_file(args.domains_file, args.verbose)
        summary = get_ecod_summary(classifications)
        
        print(f"\nECOD Domains File Validation:")
        print(f"=" * 50)
        print(f"Total protein chains: {summary['total_chains']:,}")
        print(f"Total domains: {summary['total_domains']:,}")
        print(f"Unique T-groups: {summary['unique_t_groups']:,}")
        print(f"Unique H-groups: {summary['unique_h_groups']:,}")
        print(f"Avg domains per chain: {summary['avg_domains_per_chain']:.2f}")
        
        # Check test proteins
        print(f"\nTest Protein Validation:")
        for protein_id in args.test_proteins:
            available = validate_ecod_data_for_protein(protein_id, classifications)
            status = "✓" if available else "✗"
            print(f"  {status} {protein_id}")
        
        return 0
    
    # Generate classification files
    print("Generating ECOD classification files...")
    print("=" * 60)
    
    # Create ECOD classifications for test proteins
    classifications_file = os.path.join(args.output_dir, 'ecod_classifications.csv')
    create_ecod_classifications_file(
        args.domains_file,
        classifications_file,
        args.test_proteins,
        args.verbose
    )
    
    print("\n" + "=" * 60)
    print("ECOD CLASSIFICATION FILES CREATED:")
    print("=" * 60)
    print(f"ECOD classifications: {classifications_file}")
    
    print("\nTo test ECOD integration:")
    print("python tests/test_ecod_regression.py")

if __name__ == "__main__":
    exit(main())
