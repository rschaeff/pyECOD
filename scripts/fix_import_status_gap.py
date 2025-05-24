#!/usr/bin/env python3
"""
Targeted Status Update for Successfully Imported Proteins

This script identifies proteins that have been successfully imported to pdb_analysis
but are not marked as imported in ecod_schema.process_status, and updates their status.

Usage:
    python fix_import_status_gap.py --config config.yml --batch-ids 57 58 60 61 63
"""

import os
import sys
import logging
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import yaml


def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format_str)
    return logging.getLogger(__name__)


def parse_config(config_path):
    """Parse configuration file."""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        logging.error(f"Error parsing config file: {str(e)}")
        raise


def get_db_connection(config):
    """Create database connection from config."""
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


def find_import_status_gaps(conn, batch_ids):
    """Find proteins that are imported but not marked as imported in process status."""
    
    # Query to find proteins with successful imports but outdated process status
    query = """
    SELECT 
        ps.id as process_id,
        ps.batch_id,
        p.pdb_id,
        p.chain_id,
        ps.current_stage,
        ps.status,
        pp.id as partition_id,
        pp.timestamp as import_timestamp
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    JOIN pdb_analysis.partition_proteins pp ON (
        pp.pdb_id = p.pdb_id 
        AND pp.chain_id = p.chain_id 
        AND pp.batch_id = ps.batch_id
    )
    WHERE ps.batch_id = ANY(%s)
      AND ps.is_representative = true
      AND ps.current_stage != 'pdb_analysis_imported'
      AND ps.status != 'error'
    ORDER BY ps.batch_id, p.pdb_id, p.chain_id
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, (batch_ids,))
        return cur.fetchall()


def update_import_status_batch(conn, gaps, dry_run=False):
    """Update process status for proteins with import gaps."""
    
    if not gaps:
        logging.info("No import status gaps found")
        return 0
    
    # Group by batch for reporting
    by_batch = {}
    for gap in gaps:
        batch_id = gap['batch_id']
        if batch_id not in by_batch:
            by_batch[batch_id] = []
        by_batch[batch_id].append(gap)
    
    logging.info(f"Found {len(gaps)} proteins with import status gaps across {len(by_batch)} batches:")
    for batch_id, batch_gaps in by_batch.items():
        logging.info(f"  Batch {batch_id}: {len(batch_gaps)} proteins")
    
    if dry_run:
        logging.info("DRY RUN - No changes will be made")
        return len(gaps)
    
    # Update process status for all gaps
    process_ids = [gap['process_id'] for gap in gaps]
    
    update_query = """
    UPDATE ecod_schema.process_status
    SET current_stage = 'pdb_analysis_imported',
        status = 'success',
        error_message = NULL,
        updated_at = CURRENT_TIMESTAMP
    WHERE id = ANY(%s)
    """
    
    with conn.cursor() as cur:
        cur.execute(update_query, (process_ids,))
        updated_count = cur.rowcount
    
    logging.info(f"Updated process status for {updated_count} proteins")
    return updated_count


def verify_batch_synchronization(conn, batch_ids):
    """Verify synchronization after updates."""
    
    query = """
    SELECT 
        ps.batch_id,
        COUNT(CASE WHEN ps.is_representative = true THEN 1 END) as rep_proteins,
        COUNT(pp.id) as imported_partitions,
        COUNT(CASE WHEN ps.is_representative = true 
                   AND ps.current_stage = 'pdb_analysis_imported' 
                   AND ps.status = 'success' THEN 1 END) as marked_imported
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    LEFT JOIN pdb_analysis.partition_proteins pp ON (
        pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id AND pp.batch_id = ps.batch_id
    )
    WHERE ps.batch_id = ANY(%s)
    GROUP BY ps.batch_id
    ORDER BY ps.batch_id
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, (batch_ids,))
        results = cur.fetchall()
    
    logging.info("Post-update synchronization status:")
    for result in results:
        batch_id = result['batch_id']
        rep_proteins = result['rep_proteins']
        imported = result['imported_partitions']
        marked_imported = result['marked_imported']
        
        if imported == marked_imported:
            status = "✅ SYNCHRONIZED"
        else:
            gap = imported - marked_imported
            status = f"⚠️ GAP: {gap}"
        
        logging.info(f"  Batch {batch_id}: {marked_imported}/{imported} marked imported ({status})")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Fix import status gaps for successfully imported proteins')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-ids', type=int, nargs='+', required=True,
                       help='Batch IDs to fix (e.g., 57 58 60 61 63)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be changed without making changes')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.verbose)

    # Parse config
    try:
        config = parse_config(args.config)
    except Exception as e:
        logger.error(f"Error parsing config file: {str(e)}")
        sys.exit(1)

    # Get database connection
    try:
        conn = get_db_connection(config)
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        sys.exit(1)

    try:
        conn.autocommit = False

        logger.info(f"Analyzing import status gaps for batches: {args.batch_ids}")

        # Find proteins with import status gaps
        gaps = find_import_status_gaps(conn, args.batch_ids)

        if not gaps:
            logger.info("No import status gaps found - all batches are synchronized!")
            return

        # Show sample of what will be updated
        logger.info("Sample proteins to be updated:")
        for i, gap in enumerate(gaps[:5]):
            logger.info(f"  {gap['pdb_id']}_{gap['chain_id']} (batch {gap['batch_id']}): "
                       f"{gap['current_stage']}/{gap['status']} -> pdb_analysis_imported/success")
        
        if len(gaps) > 5:
            logger.info(f"  ... and {len(gaps) - 5} more")

        # Update status
        updated_count = update_import_status_batch(conn, gaps, args.dry_run)

        if not args.dry_run and updated_count > 0:
            conn.commit()
            logger.info(f"Committed {updated_count} status updates")

            # Verify synchronization
            verify_batch_synchronization(conn, args.batch_ids)

        logger.info("Import status gap fix completed successfully")

    except Exception as e:
        logger.error(f"Error during import status fix: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        conn.rollback()
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
