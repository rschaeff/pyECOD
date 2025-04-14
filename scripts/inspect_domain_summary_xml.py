#!/usr/bin/env python3
"""
inspect_domain_summary_xml.py - Inspect a domain summary XML file
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from typing import Optional

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def inspect_summary_file(file_path: str):
    """Inspect a domain summary XML file"""
    logger = logging.getLogger("ecod.inspect")
    
    if not os.path.exists(file_path):
        logger.error(f"File not found: {file_path}")
        return False
    
    try:
        # Parse XML
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        print(f"\nSummary file: {os.path.basename(file_path)}")
        
        # Get basic information
        if root.tag == 'blast_summ_doc':
            blast_summ = root.find('blast_summ')
            if blast_summ is not None:
                pdb_id = blast_summ.get('pdb')
                chain_id = blast_summ.get('chain')
                print(f"PDB ID: {pdb_id}")
                print(f"Chain ID: {chain_id}")
        
        # Get BLAST hits
        hits_elem = root.find('.//hits')
        if hits_elem is not None:
            hit_elems = hits_elem.findall('hit')
            print(f"\nNumber of BLAST hits: {len(hit_elems)}")
            
            if hit_elems:
                print("\nTop BLAST hits:")
                for i, hit in enumerate(hit_elems[:3]):  # Show first 3 hits
                    pdb_id = hit.get('pdb_id')
                    chain_id = hit.get('chain_id')
                    evalues = hit.get('evalues')
                    
                    query_reg = hit.find('query_reg')
                    hit_reg = hit.find('hit_reg')
                    
                    print(f"  Hit {i+1}: {pdb_id}_{chain_id}")
                    print(f"    E-value: {evalues}")
                    if query_reg is not None and hit_reg is not None:
                        print(f"    Query region: {query_reg.text}")
                        print(f"    Hit region: {hit_reg.text}")
                
                if len(hit_elems) > 3:
                    print(f"  ... and {len(hit_elems) - 3} more hits")
        
        # Check for domains if present
        domains_elem = root.find('.//domains')
        if domains_elem is not None:
            domain_elems = domains_elem.findall('domain')
            print(f"\nNumber of domains: {len(domain_elems)}")
            
            if domain_elems:
                print("\nDomains:")
                for i, domain in enumerate(domain_elems[:5]):  # Show first 5 domains
                    domain_id = domain.get('id')
                    domain_range = domain.get('range')
                    
                    print(f"  Domain {i+1}: {domain_id}")
                    print(f"    Range: {domain_range}")
                    
                    # Extract classification if available
                    classification = domain.find('classification')
                    if classification is not None:
                        h_group = classification.get('h_group')
                        t_group = classification.get('t_group')
                        print(f"    Classification: {t_group}, {h_group}")
                
                if len(domain_elems) > 5:
                    print(f"  ... and {len(domain_elems) - 5} more domains")
        
        return True
    
    except ET.ParseError:
        logger.error(f"Invalid XML in file: {file_path}")
        return False
    except Exception as e:
        logger.error(f"Error inspecting file: {e}")
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Inspect a domain summary XML file')
    parser.add_argument('file_path', type=str, help='Path to the domain summary file')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    success = inspect_summary_file(args.file_path)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())