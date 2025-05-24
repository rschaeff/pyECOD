#!/usr/bin/env python3
"""
Status Cleanup Script for Historical Batch Issues

This script fixes the synchronization issues between the process tracking system
and the pdb_analysis tables that occurred due to representative protein filtering
after batch creation.

Usage:
    python fix_status_synchronization.py --config config.yml [options]

Operations:
    1. Fix non-representative protein status (mark as filtered/skipped)
    2. Verify synchronization between ecod_schema and pdb_analysis
    3. Update batch completion status
    4. Generate status reports
"""

import os
import sys
import logging
import argparse
import datetime
import psycopg2
from psycopg2.extras import RealDictCursor
import yaml
from typing import Dict, List, Any, Optional, Tuple


def setup_logging(verbose=False, log_file=None):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(level=level, format=format_str, handlers=handlers)
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


def analyze_current_status(conn, batch_ids=None):
    """Analyze current status synchronization issues."""
    
    logger = logging.getLogger(__name__)
    
    # Build batch filter
    batch_filter = ""
    params = []
    if batch_ids:
        batch_filter = " AND b.id = ANY(%s)"
        params.append(batch_ids)

    # Query to check synchronization status
    query = f"""
    WITH batch_analysis AS (
        SELECT 
            b.id as batch_id,
            b.batch_name,
            COUNT(ps.id) as total_proteins,
            COUNT(CASE WHEN ps.is_representative = true THEN 1 END) as rep_proteins,
            COUNT(CASE WHEN ps.is_representative = false THEN 1 END) as non_rep_proteins,
            
            -- Process status tracking
            COUNT(CASE WHEN ps.is_representative = true 
                       AND ps.current_stage = 'domain_partition_complete' 
                       AND ps.status = 'success' THEN 1 END) as rep_marked_complete,
            COUNT(CASE WHEN ps.is_representative = false 
                       AND ps.current_stage = 'domain_partition_complete' 
                       AND ps.status = 'success' THEN 1 END) as non_rep_marked_complete,
            COUNT(CASE WHEN ps.is_representative = false 
                       AND ps.current_stage = 'filtered_non_representative' 
                       AND ps.status = 'skipped' THEN 1 END) as non_rep_filtered,
            
            -- Actual pdb_analysis data
            COUNT(pp.id) as imported_partitions
            
        FROM ecod_schema.batch b
        JOIN ecod_schema.process_status ps ON b.id = ps.batch_id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN pdb_analysis.partition_proteins pp ON (
            pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id AND pp.batch_id = ps.batch_id
        )
        WHERE 1=1 {batch_filter}
        GROUP BY b.id, b.batch_name
        ORDER BY b.id
    )
    SELECT 
        *,
        (rep_marked_complete - imported_partitions) as status_partition_gap,
        non_rep_marked_complete as problematic_non_rep,
        CASE 
            WHEN rep_marked_complete = imported_partitions 
                 AND non_rep_marked_complete = 0 
            THEN 'synchronized'
            WHEN rep_marked_complete > imported_partitions 
            THEN 'overstated_completion'
            WHEN rep_marked_complete < imported_partitions 
            THEN 'understated_completion'
            ELSE 'mixed_issues'
        END as sync_status
    FROM batch_analysis
    """

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, params)
        return cur.fetchall()


def fix_non_representative_status(conn, batch_ids=None, dry_run=False):
    """Fix status for non-representative proteins."""
    
    logger = logging.getLogger(__name__)
    
    # Build batch filter
    batch_filter = ""
    params = []
    if batch_ids:
        batch_filter = " AND batch_id = ANY(%s)"
        params.append(batch_ids)

    # Find non-representative proteins marked as complete
    find_query = f"""
    SELECT ps.id, ps.batch_id, p.pdb_id, p.chain_id, ps.current_stage, ps.status
    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    WHERE ps.is_representative = FALSE
      AND ps.current_stage IN ('domain_partition_complete', 'completed', 'classified')
      AND ps.status = 'success'
      {batch_filter}
    """

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(find_query, params)
        problematic_records = cur.fetchall()

    if not problematic_records:
        logger.info("No non-representative proteins need status correction")
        return 0

    logger.info(f"Found {len(problematic_records)} non-representative proteins marked as complete")

    # Show sample of what will be changed
    logger.info("Sample records to be updated:")
    for i, record in enumerate(problematic_records[:5]):
        logger.info(f"  {record['pdb_id']}_{record['chain_id']} (batch {record['batch_id']}): "
                   f"{record['current_stage']}/{record['status']}")
    
    if len(problematic_records) > 5:
        logger.info(f"  ... and {len(problematic_records) - 5} more")

    if dry_run:
        logger.info("DRY RUN - No changes made")
        return len(problematic_records)

    # Update their status
    process_ids = [record['id'] for record in problematic_records]
    
    update_query = """
    UPDATE ecod_schema.process_status 
    SET current_stage = 'filtered_non_representative',
        status = 'skipped',
        updated_at = CURRENT_TIMESTAMP
    WHERE id = ANY(%s)
    """

    with conn.cursor() as cur:
        cur.execute(update_query, (process_ids,))
        updated_count = cur.rowcount

    logger.info(f"Updated {updated_count} non-representative proteins to 'filtered' status")
    return updated_count


def update_batch_completion_status(conn, batch_ids=None, dry_run=False):
    """Update batch completion status based on current protein status."""
    
    logger = logging.getLogger(__name__)
    
    # Get batches to update
    if batch_ids:
        batch_filter = " WHERE id = ANY(%s)"
        params = [batch_ids]
    else:
        batch_filter = ""
        params = []

    batch_query = f"SELECT id, batch_name FROM ecod_schema.batch{batch_filter} ORDER BY id"
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(batch_query, params)
        batches = cur.fetchall()

    updated_batches = 0

    for batch in batches:
        batch_id = batch['id']
        
        # Calculate current completion status
        status_query = """
        SELECT 
            COUNT(CASE WHEN ps.is_representative = true THEN 1 END) as total_rep,
            COUNT(CASE WHEN ps.is_representative = true 
                       AND ps.current_stage = 'domain_partition_complete' 
                       AND ps.status = 'success' THEN 1 END) as complete_rep,
            COUNT(CASE WHEN ps.is_representative = true 
                       AND ps.status = 'error' THEN 1 END) as error_rep
        FROM ecod_schema.process_status ps
        WHERE ps.batch_id = %s
        """
        
        cur.execute(status_query, (batch_id,))
        result = cur.fetchone()
        
        total = result['total_rep']
        complete = result['complete_rep']
        errors = result['error_rep']
        
        if total == 0:
            continue

        # Determine appropriate batch status
        completion_rate = complete / total
        error_rate = errors / total
        
        if completion_rate >= 0.95:  # 95% or more complete
            if error_rate == 0:
                new_status = 'domain_partition_complete'
            elif error_rate < 0.1:  # Less than 10% errors
                new_status = 'domain_partition_complete_with_errors'
            else:
                new_status = 'domain_partition_partial'
        elif completion_rate >= 0.5:  # 50% or more complete
            new_status = 'domain_partition_partial'
        else:
            new_status = 'domain_partition_incomplete'

        # Get current batch status
        cur.execute("SELECT status FROM ecod_schema.batch WHERE id = %s", (batch_id,))
        current_status = cur.fetchone()['status']

        if current_status != new_status:
            logger.info(f"Batch {batch_id} ({batch['batch_name']}): "
                       f"updating status from '{current_status}' to '{new_status}' "
                       f"({complete}/{total} complete, {errors} errors)")

            if not dry_run:
                update_query = """
                UPDATE ecod_schema.batch 
                SET status = %s, 
                    completed_items = %s, 
                    updated_at = CURRENT_TIMESTAMP
                WHERE id = %s
                """
                cur.execute(update_query, (new_status, complete, batch_id))
                
            updated_batches += 1

    if dry_run:
        logger.info(f"DRY RUN - Would update {updated_batches} batch statuses")
    else:
        logger.info(f"Updated status for {updated_batches} batches")
        
    return updated_batches


def generate_status_report(conn, batch_ids=None, output_file=None):
    """Generate a comprehensive status report."""
    
    logger = logging.getLogger(__name__)
    
    # Analyze current status
    analysis = analyze_current_status(conn, batch_ids)
    
    # Prepare report
    report_lines = []
    report_lines.append("BATCH STATUS SYNCHRONIZATION REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"Generated: {datetime.datetime.now()}")
    report_lines.append("")
    
    # Summary
    total_batches = len(analysis)
    synchronized = sum(1 for batch in analysis if batch['sync_status'] == 'synchronized')
    problematic = total_batches - synchronized
    
    report_lines.append(f"SUMMARY:")
    report_lines.append(f"  Total batches analyzed: {total_batches}")
    report_lines.append(f"  Synchronized: {synchronized}")
    report_lines.append(f"  Problematic: {problematic}")
    report_lines.append("")
    
    # Detailed analysis
    report_lines.append("DETAILED ANALYSIS:")
    report_lines.append(f"{'Batch':>5} {'Name':20} {'Rep':>5} {'Complete':>8} {'Imported':>8} {'Gap':>5} {'NonRepIssue':>10} {'Status':15}")
    report_lines.append("-" * 80)
    
    for batch in analysis:
        report_lines.append(
            f"{batch['batch_id']:5d} "
            f"{(batch['batch_name'] or '')[:20]:20} "
            f"{batch['rep_proteins']:5d} "
            f"{batch['rep_marked_complete']:8d} "
            f"{batch['imported_partitions']:8d} "
            f"{batch['status_partition_gap']:5d} "
            f"{batch['problematic_non_rep']:10d} "
            f"{batch['sync_status']:15}"
        )
    
    report_lines.append("")
    
    # Issues summary
    report_lines.append("ISSUES FOUND:")
    overstated = [b for b in analysis if b['sync_status'] == 'overstated_completion']
    understated = [b for b in analysis if b['sync_status'] == 'understated_completion']
    mixed = [b for b in analysis if b['sync_status'] == 'mixed_issues']
    
    if overstated:
        report_lines.append(f"  Overstated completion (process tracking > actual imports): {len(overstated)} batches")
        for batch in overstated[:5]:
            report_lines.append(f"    Batch {batch['batch_id']}: {batch['status_partition_gap']} excess completions")
    
    if understated:
        report_lines.append(f"  Understated completion (process tracking < actual imports): {len(understated)} batches")
        for batch in understated[:5]:
            report_lines.append(f"    Batch {batch['batch_id']}: {abs(batch['status_partition_gap'])} missing completions")
    
    if mixed:
        report_lines.append(f"  Mixed issues: {len(mixed)} batches")
    
    # Non-representative protein issues
    non_rep_issues = sum(batch['problematic_non_rep'] for batch in analysis)
    if non_rep_issues > 0:
        report_lines.append(f"  Non-representative proteins incorrectly marked complete: {non_rep_issues}")
    
    report_text = "\n".join(report_lines)
    
    # Output report
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_text)
        logger.info(f"Status report written to {output_file}")
    else:
        print(report_text)
    
    return analysis


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Fix status synchronization issues in pyECOD batches')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-ids', type=int, nargs='+', 
                       help='Specific batch IDs to process (default: all)')
    parser.add_argument('--fix-non-rep', action='store_true',
                       help='Fix non-representative protein status')
    parser.add_argument('--update-batch-status', action='store_true',
                       help='Update batch completion status')
    parser.add_argument('--report', action='store_true',
                       help='Generate status report')
    parser.add_argument('--report-file', help='Output file for status report')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be changed without making changes')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--log-file', help='Path to log file')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.verbose, args.log_file)

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

        # Default to problematic batches if no specific action requested
        if not any([args.fix_non_rep, args.update_batch_status, args.report]):
            logger.info("No specific action requested. Running analysis and generating report.")
            args.report = True

        # Generate report
        if args.report:
            logger.info("Generating status synchronization report...")
            generate_status_report(conn, args.batch_ids, args.report_file)

        # Fix non-representative protein status
        if args.fix_non_rep:
            logger.info("Fixing non-representative protein status...")
            fixed_count = fix_non_representative_status(conn, args.batch_ids, args.dry_run)
            
            if not args.dry_run and fixed_count > 0:
                conn.commit()
                logger.info(f"Committed {fixed_count} status updates")

        # Update batch completion status
        if args.update_batch_status:
            logger.info("Updating batch completion status...")
            updated_count = update_batch_completion_status(conn, args.batch_ids, args.dry_run)
            
            if not args.dry_run and updated_count > 0:
                conn.commit()
                logger.info(f"Committed {updated_count} batch status updates")

        logger.info("Status cleanup completed successfully")

    except Exception as e:
        logger.error(f"Error during status cleanup: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        conn.rollback()
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
