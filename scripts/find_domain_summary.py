#!/usr/bin/env python3
"""
find_domain_summary.py - Standalone script to find domain summary files
"""

import os
import sys
import logging
import argparse

def setup_logging(verbose=False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )

def find_domain_summary(pdb_id, chain_id, batch_path, ref_version="develop291", blast_only=False):
    """Find domain summary file for a protein chain
    
    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        batch_path: Base directory for batch files
        ref_version: Reference version
        blast_only: Whether to look for BLAST-only summary
        
    Returns:
        Path to domain summary file if found, None otherwise
    """
    logger = logging.getLogger("domain_finder")
    
    pdb_chain = f"{pdb_id}_{chain_id}"
    suffix = ".blast_only" if blast_only else ""
    
    logger.info(f"Looking for domain summary for {pdb_chain}")
    logger.info(f"Batch path: {batch_path}")
    logger.info(f"Reference version: {ref_version}")
    logger.info(f"Blast only: {blast_only}")
    
    # Define all possible locations for domain summary files
    domains_dir = os.path.join(batch_path, "domains")
    possible_paths = [
        # Standard paths
        os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domain_summary{suffix}.xml"),
        
        # Alternative standard paths (without reference version)
        os.path.join(domains_dir, f"{pdb_chain}.domain_summary{suffix}.xml"),
        
        # Legacy paths with different naming
        os.path.join(domains_dir, f"{pdb_chain}_domain_summary{suffix}.xml"),
        
        # Legacy paths in pdb_chain directory
        os.path.join(batch_path, pdb_chain, f"{pdb_chain}.{ref_version}.domain_summary{suffix}.xml"),
        os.path.join(batch_path, pdb_chain, f"{pdb_chain}.domain_summary{suffix}.xml"),
        
        # Legacy paths in ecod_dump directory
        os.path.join(batch_path, "ecod_dump", pdb_chain, f"{pdb_chain}.{ref_version}.domain_summary{suffix}.xml"),
        os.path.join(batch_path, "ecod_dump", pdb_chain, f"{pdb_chain}.domain_summary{suffix}.xml"),
        
        # Legacy paths with underscore
        os.path.join(domains_dir, f"{pdb_chain}_{ref_version}_domain_summary{suffix}.xml"),
        
        # Additional variants
        os.path.join(batch_path, f"{pdb_chain}.{ref_version}.domain_summary{suffix}.xml"),
        os.path.join(batch_path, f"{pdb_chain}_domain_summary{suffix}.xml")
    ]
    
    # Check each path
    for path in possible_paths:
        logger.debug(f"Checking path: {path}")
        if os.path.exists(path):
            file_size = os.path.getsize(path)
            logger.info(f"Found domain summary at: {path} (size: {file_size} bytes)")
            
            # Basic validation - ensure file is not empty
            if file_size == 0:
                logger.warning(f"File exists but is empty: {path}")
                continue
                
            return path
    
    logger.warning(f"No domain summary found for {pdb_chain}")
    return None

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Find domain summary files')
    parser.add_argument('--pdb-id', type=str, required=True,
                      help='PDB ID')
    parser.add_argument('--chain-id', type=str, required=True,
                      help='Chain ID')
    parser.add_argument('--batch-path', type=str, required=True,
                      help='Batch path')
    parser.add_argument('--ref-version', type=str, default='develop291',
                      help='Reference version')
    parser.add_argument('--blast-only', action='store_true',
                      help='Look for BLAST-only summary')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose)
    
    path = find_domain_summary(
        args.pdb_id, 
        args.chain_id,
        args.batch_path,
        args.ref_version,
        args.blast_only
    )
    
    if path:
        print(f"Found: {path}")
        return 0
    else:
        print(f"Not found: {args.pdb_id}_{args.chain_id}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
