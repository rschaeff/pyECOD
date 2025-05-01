#!/usr/bin/env python3
"""
domain_path_diagnosis.py - Diagnose domain summary path issues in pyECOD
"""

import os
import sys
import logging
import argparse
import psycopg2
from pathlib import Path

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.utils.path_utils import (
    get_standardized_paths, 
    get_all_evidence_paths,
    find_files_with_legacy_paths
)

def setup_logging(verbose=False, log_file=None):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )

def diagnose_process(process_id, config_path=None, fix_in_db=False):
    """Diagnose domain summary path issues for a specific process"""
    logger = logging.getLogger("path_diagnosis")
    
    # Load configuration
    config_manager = ConfigManager(config_path)
    config = config_manager.config
    
    # Initialize database connection
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Get process details
    query = """
    SELECT 
        ps.id as process_id, 
        p.pdb_id, 
        p.chain_id, 
        ps.batch_id,
        ps.relative_path,
        ps.current_stage,
        ps.status,
        b.base_path,
        b.ref_version
    FROM 
        ecod_schema.process_status ps
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    JOIN
        ecod_schema.batch b ON ps.batch_id = b.id
    WHERE 
        ps.id = %s
    """
    
    process = db.execute_dict_query(query, (process_id,))
    if not process:
        logger.error(f"Process {process_id} not found")
        return False
    
    process = process[0]
    pdb_id = process['pdb_id']
    chain_id = process['chain_id']
    batch_path = process['base_path']
    ref_version = process['ref_version']
    
    logger.info(f"Diagnosing paths for {pdb_id}_{chain_id} (Process {process_id})")
    logger.info(f"Batch path: {batch_path}")
    logger.info(f"Reference version: {ref_version}")
    
    # Check for domain summary file in database
    query = """
    SELECT pf.id, pf.file_type, pf.file_path, pf.file_exists
    FROM ecod_schema.process_file pf
    WHERE pf.process_id = %s AND pf.file_type LIKE '%summary'
    """
    
    file_records = db.execute_dict_query(query, (process_id,))
    
    logger.info(f"Database records for domain summary files:")
    if not file_records:
        logger.warning("No domain summary file records found in database")
    else:
        for record in file_records:
            file_path = record['file_path']
            full_path = os.path.join(batch_path, file_path) if not os.path.isabs(file_path) else file_path
            db_exists = record['file_exists']
            fs_exists = os.path.exists(full_path)
            
            status = "✓" if fs_exists else "✗"
            db_status = "✓" if db_exists else "✗"
            
            logger.info(f"  [{status}] {record['file_type']}: {full_path}")
            logger.info(f"     DB exists: [{db_status}], FS exists: [{status}]")
            
            if db_exists != fs_exists:
                logger.warning(f"     Database status mismatch! DB: {db_exists}, FS: {fs_exists}")
    
    # Get standardized paths
    standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)
    
    logger.info("\nStandard paths for domain summary files:")
    for file_type in ['domain_summary', 'blast_only_summary']:
        if file_type in standard_paths:
            path = standard_paths[file_type]
            exists = os.path.exists(path)
            status = "✓" if exists else "✗"
            
            logger.info(f"  [{status}] {file_type}: {path}")
    
    # Check for legacy paths
    legacy_files = find_files_with_legacy_paths(batch_path, pdb_id, chain_id, ref_version)
    
    logger.info("\nLegacy path search results:")
    for file_type in ['domain_summary', 'blast_only_summary']:
        if file_type in legacy_files:
            info = legacy_files[file_type]
            exists_at = info['exists_at']
            legacy_path = info['legacy_path']
            
            if exists_at:
                logger.info(f"  [✓] {file_type} found at: {exists_at}")
                if legacy_path and legacy_path != standard_paths[file_type]:
                    logger.info(f"     Using legacy path: {legacy_path}")
                    
                    # Update database if requested
                    if fix_in_db and file_records:
                        for record in file_records:
                            if record['file_type'] == file_type:
                                rel_path = os.path.relpath(legacy_path, batch_path)
                                db.update(
                                    "ecod_schema.process_file",
                                    {
                                        "file_path": rel_path,
                                        "file_exists": True,
                                        "file_size": os.path.getsize(legacy_path)
                                    },
                                    "id = %s",
                                    (record['id'],)
                                )
                                logger.info(f"     Updated database record to: {rel_path}")
            else:
                logger.info(f"  [✗] {file_type} not found via legacy path search")
    
    # Use comprehensive evidence path search
    logger.info("\nComprehensive evidence path search:")
    evidence_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, ref_version)
    
    for file_type in ['domain_summary', 'blast_only_summary']:
        if file_type in evidence_paths:
            info = evidence_paths[file_type]
            if info['exists_at']:
                logger.info(f"  [✓] {file_type} found at: {info['exists_at']}")
            else:
                logger.info(f"  [✗] {file_type} not found by comprehensive search")
    
    # Check related files
    logger.info("\nChecking related BLAST files that would be needed to generate a domain summary:")
    for file_type in ['chain_blast', 'domain_blast']:
        # Check standard paths
        standard_path = standard_paths.get(file_type)
        standard_exists = standard_path and os.path.exists(standard_path)
        
        # Check legacy paths
        legacy_info = legacy_files.get(file_type, {})
        legacy_path = legacy_info.get('exists_at')
        
        # Determine overall status
        if standard_exists:
            logger.info(f"  [✓] {file_type}: {standard_path}")
        elif legacy_path:
            logger.info(f"  [✓] {file_type} (legacy): {legacy_path}")
        else:
            logger.info(f"  [✗] {file_type}: Not found")
    
    # Check directory structure
    logger.info("\nDirectory structure check:")
    dirs_to_check = [
        os.path.join(batch_path, "domains"),
        os.path.join(batch_path, "blast"),
        os.path.join(batch_path, "blast", "chain"),
        os.path.join(batch_path, "blast", "domain")
    ]
    
    for directory in dirs_to_check:
        exists = os.path.isdir(directory)
        status = "✓" if exists else "✗"
        logger.info(f"  [{status}] {directory}")
    
    # Summarize findings
    logger.info("\nSummary:")
    
    # Determine if domain summary can be found
    domain_summary_found = False
    domain_summary_path = None
    
    for file_type in ['domain_summary', 'blast_only_summary']:
        if file_type in evidence_paths and evidence_paths[file_type]['exists_at']:
            domain_summary_found = True
            domain_summary_path = evidence_paths[file_type]['exists_at']
            logger.info(f"  [✓] {file_type} found at: {domain_summary_path}")
            break
    
    if not domain_summary_found:
        logger.error("  [✗] No domain summary file found! Domain partition will fail.")
        
        # Check if we have the necessary files to generate a domain summary
        chain_blast_found = (
            os.path.exists(standard_paths.get('chain_blast', '')) or
            (legacy_files.get('chain_blast', {}).get('exists_at') is not None)
        )
        
        domain_blast_found = (
            os.path.exists(standard_paths.get('domain_blast', '')) or
            (legacy_files.get('domain_blast', {}).get('exists_at') is not None)
        )
        
        if chain_blast_found and domain_blast_found:
            logger.info("  [✓] BLAST files found - domain summary can be generated")
        else:
            logger.error("  [✗] Missing BLAST files - domain summary cannot be generated")
            
            if not chain_blast_found:
                logger.error("     Missing chain BLAST results")
            if not domain_blast_found:
                logger.error("     Missing domain BLAST results")
    
    # Suggest next steps
    logger.info("\nRecommended actions:")
    
    if domain_summary_found:
        logger.info("  - Run domain partition using the EnhancedPartitionRunner script")
        logger.info(f"    python scripts/enhanced_run_domain_partition.py --config {config_path} --process-id {process_id}")
    elif chain_blast_found and domain_blast_found:
        logger.info("  - Generate domain summary first:")
        logger.info(f"    python scripts/run_domain_summary.py --config {config_path} --process-id {process_id}")
        logger.info("  - Then run domain partition")
    else:
        logger.info("  - Run missing BLAST analyses first")
        logger.info(f"    python scripts/run_blast.py --config {config_path} --process-id {process_id}")
        logger.info("  - Then generate domain summary")
        logger.info("  - Then run domain partition")
    
    return True

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Diagnose domain summary path issues')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--process-id', type=int, required=True,
                      help='Process ID to diagnose')
    parser.add_argument('--fix-in-db', action='store_true',
                      help='Fix file paths in database if legacy paths are found')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)
    
    success = diagnose_process(args.process_id, args.config, args.fix_in_db)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
