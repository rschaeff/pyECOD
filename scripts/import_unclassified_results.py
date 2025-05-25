#!/usr/bin/env python3
"""
Import Unclassified Results - The Systematic Bug Fix

This script finds unclassified domain partition files that exist on disk but aren't
tracked in the database, then imports them properly and fixes the process status.

This addresses the systematic bug where unclassified results are treated as failures.

Usage:
    python import_unclassified_results.py --config config.yml --batch-id 13 [options]
"""

import os
import sys
import logging
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import yaml
import xml.etree.ElementTree as ET
import glob
from pathlib import Path
from datetime import datetime


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


def find_untracked_proteins(conn, batch_id, include_non_rep=False):
    """Find proteins that should have partition files but no database records."""

    # Build representative filter
    rep_filter = "" if include_non_rep else "AND ps.is_representative = true"

    query = f"""
    SELECT
        ep.id as protein_id,
        ps.id as process_id,
        ep.source_id,
        ep.pdb_id,
        ep.chain_id,
        ep.length,
        ps.batch_id,
        ps.status,
        ps.error_message,
        ps.is_representative,
        b.batch_name,
        b.base_path
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
    JOIN ecod_schema.batch b ON ps.batch_id = b.id
    LEFT JOIN pdb_analysis.partition_proteins pp ON (
        ep.source_id = (pp.pdb_id || '_' || pp.chain_id)
        AND pp.batch_id = ps.batch_id
    )
    WHERE ps.batch_id = %s
      AND ps.status = 'error'
      AND ps.current_stage = 'domain_partition_failed'
      AND pp.id IS NULL  -- Not in partition table
      {rep_filter}
    ORDER BY ep.source_id
    """

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, (batch_id,))
        return cur.fetchall()


def find_partition_file(protein, base_data_dir="/data/ecod/pdb_updates"):
    """Find partition file for a protein on disk."""

    batch_name = protein['batch_name']
    pdb_id = protein['pdb_id']
    chain_id = protein['chain_id']

    # Construct path based on batch name structure
    batch_dir = Path(base_data_dir) / "batches" / batch_name / "domains"

    # Look for partition file (domains.xml, not domain_summary.xml)
    partition_pattern = f"{pdb_id}_{chain_id}.develop291.domains.xml"
    partition_path = batch_dir / partition_pattern

    if partition_path.exists():
        return {
            'path': str(partition_path),
            'size': partition_path.stat().st_size,
            'exists': True
        }

    return {'exists': False}


def parse_partition_file(file_path):
    """Parse a domain partition file to extract data."""
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()

        # Extract basic info
        pdb_id = root.get('pdb_id')
        chain_id = root.get('chain_id')
        reference = root.get('reference', 'develop291')

        # Check classification status
        is_classified = root.get('is_classified', 'false').lower() == 'true'
        is_unclassified = root.get('is_unclassified', 'false').lower() == 'true'

        # If neither is explicitly set, infer from domain count
        if not is_classified and not is_unclassified:
            domains = root.findall('.//domain')
            is_classified = len(domains) > 0
            is_unclassified = len(domains) == 0

        # Extract metadata
        metadata = root.find('metadata')
        if metadata is not None:
            sequence_length = int(metadata.findtext('sequence_length', '0'))
            domain_count = int(metadata.findtext('domain_count', '0'))
            coverage = float(metadata.findtext('coverage', '0.0'))
            residues_assigned = int(metadata.findtext('residues_assigned', '0'))
        else:
            # Fall back to counting domains
            domains = root.findall('.//domain')
            domain_count = len(domains)
            sequence_length = 0
            coverage = 0.0
            residues_assigned = 0

        return {
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'reference': reference,
            'is_classified': is_classified,
            'is_unclassified': is_unclassified,
            'domain_count': domain_count,
            'coverage': coverage,
            'sequence_length': sequence_length,
            'residues_assigned': residues_assigned,
            'parse_success': True
        }

    except Exception as e:
        return {
            'parse_success': False,
            'error': str(e)
        }


def import_partition_to_database(conn, protein, file_info, parsed_data, dry_run=False):
    """Import a partition file to the database and fix process status."""

    if not parsed_data['parse_success']:
        return False, f"Parse failed: {parsed_data.get('error', 'Unknown error')}"

    try:
        with conn.cursor() as cur:
            # 1. Insert into pdb_analysis.partition_proteins
            insert_partition_query = """
            INSERT INTO pdb_analysis.partition_proteins (
                pdb_id, chain_id, batch_id,
                is_classified, is_unclassified,
                sequence_length, coverage,
                created_at, updated_at
            ) VALUES (
                %s, %s, %s, %s, %s, %s, %s, %s, %s
            ) RETURNING id
            """

            partition_data = (
                parsed_data['pdb_id'],
                parsed_data['chain_id'],
                protein['batch_id'],
                parsed_data['is_classified'],
                parsed_data['is_unclassified'],
                parsed_data['sequence_length'],
                parsed_data['coverage'],
                datetime.now(),
                datetime.now()
            )

            if not dry_run:
                cur.execute(insert_partition_query, partition_data)
                partition_id = cur.fetchone()[0]
            else:
                partition_id = "DRY_RUN"

            # 2. Add file record to process_file table
            insert_file_query = """
            INSERT INTO ecod_schema.process_file (
                process_id, file_type, file_path,
                file_exists, file_size, last_checked
            ) VALUES (
                %s, %s, %s, %s, %s, %s
            )
            """

            file_data = (
                protein['process_id'],
                'domain_partition',
                file_info['path'],
                True,
                file_info['size'],
                datetime.now()
            )

            if not dry_run:
                cur.execute(insert_file_query, file_data)

            # 3. Update process status to success
            if parsed_data['is_classified']:
                new_stage = 'classified'
                new_status = 'success'
                success_message = f"Classified result with {parsed_data['domain_count']} domains"
            else:
                new_stage = 'unclassified'
                new_status = 'success'
                success_message = f"Unclassified result (no domains found)"

            update_status_query = """
            UPDATE ecod_schema.process_status
            SET status = %s,
                current_stage = %s,
                error_message = %s,
                updated_at = %s
            WHERE id = %s
            """

            if not dry_run:
                cur.execute(update_status_query, (
                    new_status, new_stage, success_message, datetime.now(), protein['process_id']
                ))
                conn.commit()

            result_type = "classified" if parsed_data['is_classified'] else "unclassified"
            return True, f"Successfully imported {result_type} result (partition_id={partition_id})"

    except Exception as e:
        if not dry_run:
            conn.rollback()
        return False, f"Database error: {str(e)}"


def main():
    parser = argparse.ArgumentParser(description='Import unclassified results and fix systematic bug')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True, help='Batch ID to fix')
    parser.add_argument('--limit', type=int, help='Limit number of proteins to process')
    parser.add_argument('--base-data-dir', default='/data/ecod/pdb_updates', help='Base data directory')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without making changes')
    parser.add_argument('--include-non-rep', action='store_true', help='Include non-representative proteins')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    logger = setup_logging(args.verbose)

    try:
        config = parse_config(args.config)
        conn = get_db_connection(config)
    except Exception as e:
        logger.error(f"Setup error: {str(e)}")
        sys.exit(1)

    mode = "🧪 DRY RUN" if args.dry_run else "🚀 LIVE MODE"
    scope = "all proteins" if args.include_non_rep else "representative only"
    print(f"\n{'='*80}")
    print(f"IMPORTING UNCLASSIFIED RESULTS - {mode}")
    print(f"Batch ID: {args.batch_id}")
    print(f"Scope: {scope}")
    print(f"{'='*80}")

    try:
        # Find untracked proteins
        logger.info(f"Finding proteins with missing partition records...")
        untracked_proteins = find_untracked_proteins(conn, args.batch_id, args.include_non_rep)
        
        if not untracked_proteins:
            print(f"✅ No untracked proteins found in batch {args.batch_id}")
            print(f"All proteins appear to be properly tracked!")
            return
            
        print(f"🔍 Found {len(untracked_proteins)} proteins marked as errors")
        
        if args.limit:
            untracked_proteins = untracked_proteins[:args.limit]
            print(f"📊 Processing first {len(untracked_proteins)} proteins")
        
        # Process each protein
        processed = 0
        imported_classified = 0
        imported_unclassified = 0
        skipped_no_files = 0
        genuine_failures = 0
        errors = 0
        
        for i, protein in enumerate(untracked_proteins):
            source_id = protein['source_id']
            print(f"\n[{i+1}/{len(untracked_proteins)}] Processing {source_id} (length: {protein['length']})")
            
            # Find partition file
            file_info = find_partition_file(protein, args.base_data_dir)
            
            if not file_info['exists']:
                print(f"   ❌ No partition file found - genuine failure")
                genuine_failures += 1
                continue
            
            print(f"   📁 Found partition file ({file_info['size']} bytes)")
            
            # Parse the file
            parsed_data = parse_partition_file(file_info['path'])
            
            if not parsed_data['parse_success']:
                print(f"   ❌ Parse failed: {parsed_data.get('error', 'Unknown')}")
                errors += 1
                continue
            
            # Show what we found
            if parsed_data['is_classified']:
                result_type = f"✅ CLASSIFIED ({parsed_data['domain_count']} domains)"
                result_category = "classified"
            else:
                result_type = f"⚠️ UNCLASSIFIED (0 domains)"
                result_category = "unclassified"
            
            print(f"   🎯 Found: {result_type}")
            print(f"   📊 Coverage: {parsed_data['coverage']:.1%}, Length: {parsed_data['sequence_length']}")
            
            # Import to database
            success, message = import_partition_to_database(
                conn, protein, file_info, parsed_data, args.dry_run
            )
            
            if success:
                print(f"   ✅ {message}")
                if result_category == "classified":
                    imported_classified += 1
                else:
                    imported_unclassified += 1
            else:
                print(f"   ❌ Import failed: {message}")
                errors += 1
            
            processed += 1
        
        # Summary
        print(f"\n{'='*80}")
        print(f"PROCESSING SUMMARY - {mode}")
        print(f"{'='*80}")
        print(f"Proteins examined: {len(untracked_proteins)}")
        print(f"Proteins processed: {processed}")
        print(f"✅ Imported classified results: {imported_classified}")
        print(f"⚠️ Imported unclassified results: {imported_unclassified}")
        print(f"❌ Genuine failures (no files): {genuine_failures}")
        print(f"💥 Errors during processing: {errors}")
        
        total_imported = imported_classified + imported_unclassified
        total_fixed = total_imported
        
        if total_fixed > 0:
            print(f"\n🎉 SUCCESS IMPACT:")
            print(f"   • {total_fixed} proteins converted from 'error' to 'success'")
            print(f"   • {imported_unclassified} unclassified results properly recognized")
            print(f"   • Batch success rate improved significantly!")
            
            if args.dry_run:
                print(f"\n🔧 To apply these changes, run without --dry-run")
            else:
                print(f"\n🚀 Changes applied! Check the batch status to see improvement.")

    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
