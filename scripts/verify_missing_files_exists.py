#!/usr/bin/env python3
"""
Verify Missing Files Exist on Disk

This script checks if domain partition files exist on disk for proteins 
that have no file records in the database, confirming the tracking bug theory.

Usage:
    python verify_missing_files_exist.py --config config.yml --batch-id 30 [options]
"""

import os
import sys
import logging
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import yaml
import glob
from pathlib import Path


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format_str)
    return logging.getLogger(__name__)


def parse_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_db_connection(config):
    db_config = config.get('database', {})
    try:
        conn = psycopg2.connect(
            host=db_config.get('host', 'dione'),
            port=db_config.get('port', 45000),
            dbname=db_config.get('name', 'ecod_protein'),
            user=db_config.get('user', 'ecod'),
            password=db_config.get('password', '')
        )
        return conn
    except psycopg2.Error as e:
        logging.error(f"Database connection error: {e}")
        raise


def get_missing_file_proteins(conn, batch_id):
    """Get proteins that should have partition files but no database records."""
    
    query = """
    SELECT 
        ep.source_id,
        ep.pdb_id,
        ep.chain_id,
        ps.batch_id,
        ps.status,
        ps.error_message,
        b.base_path,
        b.ref_version
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
    JOIN ecod_schema.batch b ON ps.batch_id = b.id
    LEFT JOIN pdb_analysis.partition_proteins pp ON (
        ep.source_id = (pp.pdb_id || '_' || pp.chain_id) 
        AND pp.batch_id = ps.batch_id
    )
    LEFT JOIN ecod_schema.process_file pf ON (
        ps.id = pf.process_id 
        AND pf.file_type = 'domain_partition'
    )
    WHERE ps.batch_id = %s
      AND ps.is_representative = true
      AND ps.status = 'error'
      AND pp.id IS NULL  -- Not in partition table
      AND pf.id IS NULL  -- Not in file tracking
    ORDER BY ep.source_id
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, (batch_id,))
        return cur.fetchall()


def check_file_exists(protein, base_data_dir="/data/ecod/pdb_updates"):
    """Check if partition files exist for a protein."""
    
    batch_id = protein['batch_id']
    pdb_id = protein['pdb_id'] 
    chain_id = protein['chain_id']
    ref_version = protein['ref_version'] or 'develop291'
    
    # Expected file paths
    batch_dir = Path(base_data_dir) / "batches" / str(batch_id)
    domains_dir = batch_dir / "domains"
    
    expected_files = [
        domains_dir / f"{pdb_id}_{chain_id}.{ref_version}.domains_v14.xml",
        batch_dir / "domain_partitions" / f"{pdb_id}_{chain_id}_{ref_version}_partition.xml"
    ]
    
    file_status = {}
    for file_path in expected_files:
        if file_path.exists():
            file_status[str(file_path)] = {
                'exists': True,
                'size': file_path.stat().st_size,
                'type': 'domain_summary' if 'domains_v14' in str(file_path) else 'domain_partition'
            }
        else:
            file_status[str(file_path)] = {
                'exists': False,
                'type': 'domain_summary' if 'domains_v14' in str(file_path) else 'domain_partition'
            }
    
    # Also check for any files matching the pattern
    domain_pattern = str(domains_dir / f"{pdb_id}_{chain_id}*.xml")
    partition_pattern = str(batch_dir / "domain_partitions" / f"{pdb_id}_{chain_id}*.xml")
    
    found_files = []
    found_files.extend(glob.glob(domain_pattern))
    found_files.extend(glob.glob(partition_pattern))
    
    return {
        'expected_files': file_status,
        'found_files': found_files,
        'has_any_files': len(found_files) > 0
    }


def main():
    parser = argparse.ArgumentParser(description='Verify missing files exist on disk')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True, help='Batch ID to check')
    parser.add_argument('--limit', type=int, default=50, help='Limit number of proteins to check')
    parser.add_argument('--base-data-dir', default='/data/ecod/pdb_updates', help='Base data directory')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    logger = setup_logging(args.verbose)

    try:
        config = parse_config(args.config)
        conn = get_db_connection(config)
    except Exception as e:
        logger.error(f"Setup error: {str(e)}")
        sys.exit(1)

    try:
        # Get proteins missing file records
        logger.info(f"Getting proteins with missing file records from batch {args.batch_id}...")
        missing_proteins = get_missing_file_proteins(conn, args.batch_id)
        
        if not missing_proteins:
            print(f"No proteins with missing file records found in batch {args.batch_id}")
            return
            
        print(f"\n{'='*80}")
        print(f"VERIFYING FILE EXISTENCE FOR BATCH {args.batch_id}")
        print(f"Found {len(missing_proteins)} proteins with missing database records")
        print(f"{'='*80}")
        
        # Check each protein
        verified_count = 0
        missing_count = 0
        files_exist_count = 0
        
        for i, protein in enumerate(missing_proteins[:args.limit]):
            source_id = protein['source_id']
            
            print(f"\n[{i+1}/{min(len(missing_proteins), args.limit)}] {source_id}")
            print(f"   Status: {protein['status']}")
            if protein['error_message']:
                print(f"   Error: {protein['error_message'][:100]}...")
            
            # Check files
            file_check = check_file_exists(protein, args.base_data_dir)
            
            if file_check['has_any_files']:
                files_exist_count += 1
                print(f"   🎯 FILES FOUND! (This confirms the bug)")
                
                for file_path in file_check['found_files']:
                    file_size = os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    file_type = 'domain_summary' if 'domains_v14' in file_path else 'domain_partition'
                    print(f"      ✅ {file_type}: {os.path.basename(file_path)} ({file_size} bytes)")
            else:
                missing_count += 1
                print(f"   ❌ No files found")
                
            verified_count += 1
        
        # Summary
        print(f"\n{'='*80}")
        print(f"VERIFICATION SUMMARY")
        print(f"{'='*80}")
        print(f"Proteins checked: {verified_count}")
        print(f"Files exist on disk: {files_exist_count} ({files_exist_count/verified_count*100:.1f}%)")
        print(f"No files found: {missing_count} ({missing_count/verified_count*100:.1f}%)")
        
        if files_exist_count > 0:
            print(f"\n🚨 BREAKTHROUGH CONFIRMED!")
            print(f"   {files_exist_count} proteins have files on disk but no database tracking")
            print(f"   This proves the file tracking bug theory!")
            print(f"\n🔧 RECOMMENDED ACTION:")
            print(f"   1. Create script to scan disk for missing partition files")
            print(f"   2. Import untracked files into database")
            print(f"   3. Add process_file records for tracking")
            print(f"   4. Update process status from 'error' to 'complete'")

    except Exception as e:
        logger.error(f"Error during verification: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
