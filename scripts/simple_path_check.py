#!/usr/bin/env python3
"""
simple_path_check.py - Very simple check for domain summary files
"""

import os
import sys
import logging
import argparse

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager

def setup_logging(verbose=False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )

def check_process_paths(process_id, config_path=None):
    """Simple check for domain summary paths for a process ID"""
    logger = logging.getLogger("simple_path_check")
    
    # Load configuration
    config_manager = ConfigManager(config_path)
    config = config_manager.config
    
    # Initialize database connection
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Get process details
    try:
        query = """
        SELECT 
            ps.id, p.pdb_id, p.chain_id, b.base_path, b.ref_version
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            ps.id = %s
        """
        
        row = db.execute_query(query, (process_id,))
        if not row or not row[0]:
            logger.error(f"Process {process_id} not found")
            return False
        
        # Extract data from the first row
        pdb_id = row[0][1]
        chain_id = row[0][2]
        batch_path = row[0][3]
        ref_version = row[0][4]
        
        logger.info(f"Checking paths for {pdb_id}_{chain_id} (Process {process_id})")
        logger.info(f"Batch path: {batch_path}")
        logger.info(f"Reference version: {ref_version}")
        
        # Simply check for domain summary files
        domains_dir = os.path.join(batch_path, "domains")
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Check domain summary files
        summary_paths = [
            os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domain_summary.xml"),
            os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml"),
            os.path.join(batch_path, pdb_chain, f"{pdb_chain}.domain_summary.xml"),
            os.path.join(batch_path, "ecod_dump", pdb_chain, f"{pdb_chain}.domain_summary.xml")
        ]
        
        logger.info("\nChecking domain summary files:")
        summary_found = False
        found_path = None
        
        for path in summary_paths:
            exists = os.path.exists(path)
            status = "✓" if exists else "✗"
            logger.info(f"  [{status}] {path}")
            
            if exists and not summary_found:
                summary_found = True
                found_path = path
        
        # Check BLAST files
        chain_blast_paths = [
            os.path.join(batch_path, "blast", "chain", f"{pdb_chain}.{ref_version}.xml"),
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.chainwise_blast.xml")
        ]
        
        logger.info("\nChecking chain BLAST files:")
        chain_blast_found = False
        
        for path in chain_blast_paths:
            exists = os.path.exists(path)
            status = "✓" if exists else "✗"
            logger.info(f"  [{status}] {path}")
            
            if exists:
                chain_blast_found = True
        
        domain_blast_paths = [
            os.path.join(batch_path, "blast", "domain", f"{pdb_chain}.{ref_version}.xml"),
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.blast.xml")
        ]
        
        logger.info("\nChecking domain BLAST files:")
        domain_blast_found = False
        
        for path in domain_blast_paths:
            exists = os.path.exists(path)
            status = "✓" if exists else "✗"
            logger.info(f"  [{status}] {path}")
            
            if exists:
                domain_blast_found = True
        
        # Summary
        logger.info("\nSummary:")
        if summary_found:
            logger.info(f"  [✓] Domain summary found at: {found_path}")
            logger.info("\nNext steps:")
            logger.info(f"  1. Run domain partition with process ID {process_id}")
        else:
            logger.error("  [✗] No domain summary file found!")
            if chain_blast_found and domain_blast_found:
                logger.info("\nSuggestion: Generate domain summary first")
            else:
                logger.info("\nSuggestion: Run BLAST analyses first")
        
        return True
    
    except Exception as e:
        logger.error(f"Error checking paths: {e}")
        return False

def check_file_paths(pdb_id, chain_id, batch_path, ref_version=None):
    """Simple check for domain summary paths for a PDB ID and chain"""
    logger = logging.getLogger("simple_path_check")
    
    if not ref_version:
        ref_version = "develop291"  # Default reference version
    
    logger.info(f"Checking paths for {pdb_id}_{chain_id}")
    logger.info(f"Batch path: {batch_path}")
    logger.info(f"Reference version: {ref_version}")
    
    # Simply check for domain summary files
    domains_dir = os.path.join(batch_path, "domains")
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    # Check domain summary files
    summary_paths = [
        os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domain_summary.xml"),
        os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml"),
        os.path.join(batch_path, pdb_chain, f"{pdb_chain}.domain_summary.xml"),
        os.path.join(batch_path, "ecod_dump", pdb_chain, f"{pdb_chain}.domain_summary.xml")
    ]
    
    logger.info("\nChecking domain summary files:")
    summary_found = False
    found_path = None
    
    for path in summary_paths:
        exists = os.path.exists(path)
        status = "✓" if exists else "✗"
        logger.info(f"  [{status}] {path}")
        
        if exists and not summary_found:
            summary_found = True
            found_path = path
    
    # Check BLAST files
    chain_blast_paths = [
        os.path.join(batch_path, "blast", "chain", f"{pdb_chain}.{ref_version}.xml"),
        os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.chainwise_blast.xml")
    ]
    
    logger.info("\nChecking chain BLAST files:")
    chain_blast_found = False
    
    for path in chain_blast_paths:
        exists = os.path.exists(path)
        status = "✓" if exists else "✗"
        logger.info(f"  [{status}] {path}")
        
        if exists:
            chain_blast_found = True
    
    domain_blast_paths = [
        os.path.join(batch_path, "blast", "domain", f"{pdb_chain}.{ref_version}.xml"),
        os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.blast.xml")
    ]
    
    logger.info("\nChecking domain BLAST files:")
    domain_blast_found = False
    
    for path in domain_blast_paths:
        exists = os.path.exists(path)
        status = "✓" if exists else "✗"
        logger.info(f"  [{status}] {path}")
        
        if exists:
            domain_blast_found = True
    
    # Summary
    logger.info("\nSummary:")
    if summary_found:
        logger.info(f"  [✓] Domain summary found at: {found_path}")
        logger.info("\nNext steps:")
        logger.info("  1. Run domain partition using this domain summary file")
    else:
        logger.error("  [✗] No domain summary file found!")
        if chain_blast_found and domain_blast_found:
            logger.info("\nSuggestion: Generate domain summary first")
        else:
            logger.info("\nSuggestion: Run BLAST analyses first")
    
    return True

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Simple check for domain summary paths')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--process-id', type=int,
                      help='Process ID to check')
    parser.add_argument('--pdb-id', type=str,
                      help='PDB ID to check')
    parser.add_argument('--chain-id', type=str,
                      help='Chain ID to check')
    parser.add_argument('--batch-path', type=str,
                      help='Batch path (required if using pdb-id and chain-id)')
    parser.add_argument('--ref-version', type=str, default='develop291',
                      help='Reference version')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose)
    
    logger = logging.getLogger("main")
    
    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)
    
    # Check if we're using process ID or direct PDB/chain IDs
    if args.process_id:
        return 0 if check_process_paths(args.process_id, args.config) else 1
    elif args.pdb_id and args.chain_id and args.batch_path:
        return 0 if check_file_paths(args.pdb_id, args.chain_id, args.batch_path, args.ref_version) else 1
    else:
        logger.error("Must provide either --process-id OR (--pdb-id, --chain-id, and --batch-path)")
        return 1

if __name__ == "__main__":
    sys.exit(main())
