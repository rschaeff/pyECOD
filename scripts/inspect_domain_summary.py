#!/usr/bin/env python3
"""
inspect_domain_summary.py - Inspect a specific domain summary file
"""

import os
import sys
import json
import logging
import argparse
from typing import Optional

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def inspect_summary_file(file_path: str):
    """Inspect a domain summary file"""
    logger = logging.getLogger("ecod.inspect")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    try:
        with open(file_path, 'r') as f:
            content = json.load(f)
        
        # Print summary information
        print(f"\nSummary file: {os.path.basename(file_path)}")
        print(f"PDB ID: {content.get('pdb_id')}")
        print(f"Chain ID: {content.get('chain_id')}")
        
        # Sequence information
        sequence = content.get('sequence', '')
        print(f"Sequence length: {len(sequence)}")
        if sequence:
            print(f"Sequence (first 30 aa): {sequence[:30]}...")
        
        # Domain information
        domains = content.get('domains', [])
        print(f"\nNumber of domains: {len(domains)}")
        
        if domains:
            print("\nDomains:")
            for i, domain in enumerate(domains[:5]):  # Show first 5 domains
                print(f"  Domain {i+1}:")
                print(f"    Range: {domain.get('range')}")
                print(f"    Score: {domain.get('score')}")
                print(f"    Evidence type: {domain.get('evidence_type')}")
                
                # Show domain classification if available
                if 'classification' in domain:
                    cls = domain['classification']
                    print(f"    Classification: {cls.get('t_group')}, {cls.get('h_group')}, {cls.get('x_group')}")
            
            if len(domains) > 5:
                print(f"  ... and {len(domains) - 5} more domains")
        
        # BLAST results
        blast_results = content.get('blast_results', [])
        print(f"\nNumber of BLAST results: {len(blast_results)}")
        
        if blast_results and len(blast_results) > 0:
            print("\nTop BLAST hits:")
            for i, hit in enumerate(blast_results[:3]):  # Show first 3 hits
                print(f"  Hit {i+1}:")
                print(f"    Query coverage: {hit.get('query_coverage', 'N/A')}")
                print(f"    Identity: {hit.get('identity', 'N/A')}")
                print(f"    E-value: {hit.get('evalue', 'N/A')}")
                print(f"    Hit ID: {hit.get('hit_id', 'N/A')}")
            
            if len(blast_results) > 3:
                print(f"  ... and {len(blast_results) - 3} more hits")
        
        # Check for errors
        errors = content.get('errors', [])
        if errors:
            print("\nErrors:")
            for error in errors:
                print(f"  {error}")
        
        return True
    
    except json.JSONDecodeError:
        logger.error(f"Invalid JSON in file: {file_path}")
        return False
    except Exception as e:
        logger.error(f"Error inspecting file: {e}")
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Inspect a domain summary file')
    parser.add_argument('file_path', type=str, help='Path to the domain summary file')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    success = inspect_summary_file(args.file_path)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())