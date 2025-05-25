#!/usr/bin/env python3
"""
Filesystem Reality Checker - FIXED VERSION

This script verifies the actual state of domain partition files on the filesystem
vs what the database tracking thinks is happening.

Updated based on actual filesystem structure:
- Files are directly in {batch_path}/domains/
- Pattern: {pdb_id}_{chain_id}.{ref_version}.domains.xml
- Summary files: {pdb_id}_{chain_id}.{ref_version}.domain_summary.xml
"""

import os
import sys
import logging
import argparse
import json
import yaml
import psycopg2
import glob
import xml.etree.ElementTree as ET
from psycopg2.extras import RealDictCursor
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Set
from dataclasses import dataclass, asdict
from pathlib import Path
import re
from collections import defaultdict, Counter


@dataclass
class FileSystemReality:
    """Reality check results for a batch."""
    batch_id: int
    batch_name: str
    batch_path: str

    # Database vs Filesystem comparison
    db_says_complete: int
    fs_actually_complete: int
    db_says_missing: int
    fs_actually_missing: int

    # File existence reality
    files_exist_not_in_db: int
    files_in_db_not_exist: int
    files_correctly_tracked: int

    # Partition analysis
    truly_unclassified: int
    legitimate_peptides: int
    files_missing_entirely: int
    sync_errors: int

    # Additional file analysis
    domain_partition_files: int
    domain_summary_files: int

    # Representative protein accuracy
    db_rep_count: int
    fs_rep_count: int
    rep_mismatch: int

    # Recommendations
    needs_import_run: bool
    needs_status_sync: bool
    needs_batch_stage_update: bool
    suggested_batch_status: str


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
            return yaml.safe_load(f)
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


def get_batch_filesystem_info(conn, batch_ids=None):
    """Get batch information needed for filesystem verification."""

    batch_filter = ""
    params = []
    if batch_ids:
        batch_filter = " AND b.id = ANY(%s)"
        params.append(batch_ids)

    query = f"""
    SELECT
        b.id as batch_id,
        b.batch_name,
        b.base_path,
        b.status as batch_status,
        b.ref_version,
        b.type as batch_type,
        COUNT(ps.id) as total_proteins,
        COUNT(CASE WHEN ps.is_representative = true THEN 1 END) as rep_proteins,
        COUNT(CASE WHEN ps.current_stage = 'domain_partition_complete'
                   AND ps.status = 'success' THEN 1 END) as db_complete_count,
        COUNT(CASE WHEN pf.file_exists = true
                   AND pf.file_type = 'domain_partition' THEN 1 END) as db_partition_files,
        COUNT(CASE WHEN pf.file_exists = true
                   AND pf.file_type = 'domain_summary' THEN 1 END) as db_summary_files
    FROM ecod_schema.batch b
    JOIN ecod_schema.process_status ps ON b.id = ps.batch_id
    LEFT JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE b.type IN ('pdb_hhsearch', 'domain_analysis')
    {batch_filter}
    GROUP BY b.id, b.batch_name, b.base_path, b.status, b.ref_version, b.type
    ORDER BY b.id
    """

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, params)
        return cur.fetchall()


def get_protein_details_for_batch(conn, batch_id):
    """Get detailed protein information for filesystem verification."""

    query = """
    SELECT
        ps.id as process_id,
        p.pdb_id,
        p.chain_id,
        p.source_id,
        ps.is_representative,
        ps.current_stage,
        ps.status,
        ps.relative_path,

        -- File records from database
        pf_partition.file_path as db_partition_path,
        pf_partition.file_exists as db_partition_exists,
        pf_summary.file_path as db_summary_path,
        pf_summary.file_exists as db_summary_exists,

        -- Actual partition data (if imported)
        pp.id as partition_record_id,
        pp.is_classified,
        pp.sequence_length,
        COUNT(pd.id) as domain_count

    FROM ecod_schema.process_status ps
    JOIN ecod_schema.protein p ON ps.protein_id = p.id
    LEFT JOIN ecod_schema.process_file pf_partition ON (
        ps.id = pf_partition.process_id AND pf_partition.file_type = 'domain_partition'
    )
    LEFT JOIN ecod_schema.process_file pf_summary ON (
        ps.id = pf_summary.process_id AND pf_summary.file_type = 'domain_summary'
    )
    LEFT JOIN pdb_analysis.partition_proteins pp ON (
        p.pdb_id = pp.pdb_id AND p.chain_id = pp.chain_id AND pp.batch_id = ps.batch_id
    )
    LEFT JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id

    WHERE ps.batch_id = %s
    GROUP BY ps.id, p.pdb_id, p.chain_id, p.source_id, ps.is_representative,
             ps.current_stage, ps.status, ps.relative_path,
             pf_partition.file_path, pf_partition.file_exists,
             pf_summary.file_path, pf_summary.file_exists,
             pp.id, pp.is_classified, pp.sequence_length
    ORDER BY p.pdb_id, p.chain_id
    """

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, (batch_id,))
        return cur.fetchall()


def find_actual_files(batch_path, ref_version):
    """Find all actual domain files on the filesystem.

    Based on actual structure:
    - Domain partition files: {batch_path}/domains/{pdb_id}_{chain_id}.{ref_version}.domains.xml
    - Domain summary files: {batch_path}/domains/{pdb_id}_{chain_id}.{ref_version}.domain_summary.xml
    """

    domains_dir = os.path.join(batch_path, 'domains')
    if not os.path.exists(domains_dir):
        return {}, {}

    # Find partition files (domains.xml)
    partition_pattern = os.path.join(domains_dir, f"*.{ref_version}.domains.xml")
    partition_files = {}

    for file_path in glob.glob(partition_pattern):
        basename = os.path.basename(file_path)
        # Extract protein ID from filename: {pdb_id}_{chain_id}.{ref_version}.domains.xml
        protein_part = basename.replace(f".{ref_version}.domains.xml", "")

        if '_' in protein_part:
            pdb_id, chain_id = protein_part.split('_', 1)  # Split on first underscore only
            protein_key = f"{pdb_id.lower()}_{chain_id}"
            partition_files[protein_key] = file_path

    # Find summary files (domain_summary.xml)
    summary_pattern = os.path.join(domains_dir, f"*.{ref_version}.domain_summary.xml")
    summary_files = {}

    for file_path in glob.glob(summary_pattern):
        basename = os.path.basename(file_path)
        # Extract protein ID from filename: {pdb_id}_{chain_id}.{ref_version}.domain_summary.xml
        protein_part = basename.replace(f".{ref_version}.domain_summary.xml", "")

        if '_' in protein_part:
            pdb_id, chain_id = protein_part.split('_', 1)  # Split on first underscore only
            protein_key = f"{pdb_id.lower()}_{chain_id}"
            summary_files[protein_key] = file_path

    return partition_files, summary_files


def analyze_partition_file_content(file_path):
    """Analyze the content of a partition file to understand classification status."""

    try:
        # Check file size first - very small files are often peptides/unclassified
        file_size = os.path.getsize(file_path)

        tree = ET.parse(file_path)
        root = tree.getroot()

        # Look for domain elements - could be in various locations
        domains = []

        # Try different common patterns for domain elements
        possible_paths = [
            './/domain',
            './/domains/domain',
            './/domain_partition/domains/domain',
            './/partition_domains/domain'
        ]

        for path in possible_paths:
            found_domains = root.findall(path)
            if found_domains:
                domains = found_domains
                break

        # Check for classification status indicators
        is_classified_attr = root.get('is_classified')
        is_classified = is_classified_attr == 'true' if is_classified_attr else None

        # Get sequence length if available
        sequence_length = 0

        # Try to find sequence length in various ways
        sequence_elem = root.find('.//sequence')
        if sequence_elem is not None:
            if sequence_elem.text:
                # Remove whitespace and count actual residues
                clean_seq = sequence_elem.text.replace('\n', '').replace(' ', '').replace('\t', '')
                sequence_length = len(clean_seq)
            elif sequence_elem.get('length'):
                try:
                    sequence_length = int(sequence_elem.get('length'))
                except ValueError:
                    pass

        # Also check for length attribute on root
        if sequence_length == 0:
            length_attr = root.get('sequence_length')
            if length_attr:
                try:
                    sequence_length = int(length_attr)
                except ValueError:
                    pass

        # Determine classification status based on content
        domain_count = len(domains)
        is_peptide = sequence_length > 0 and sequence_length < 30
        has_domains = domain_count > 0

        # Use explicit classification if available, otherwise infer
        if is_classified is not None:
            final_is_classified = is_classified and not is_peptide
        else:
            final_is_classified = has_domains and not is_peptide

        is_unclassified = not has_domains and not is_peptide

        return {
            'domain_count': domain_count,
            'sequence_length': sequence_length,
            'file_size': file_size,
            'is_peptide': is_peptide,
            'is_classified': final_is_classified,
            'is_unclassified': is_unclassified,
            'file_valid': True,
            'explicit_classification': is_classified
        }

    except Exception as e:
        return {
            'domain_count': 0,
            'sequence_length': 0,
            'file_size': 0,
            'is_peptide': False,
            'is_classified': False,
            'is_unclassified': False,
            'file_valid': False,
            'error': str(e)
        }


def verify_batch_reality(conn, batch_info, include_non_rep=False, logger=None):
    """Verify the reality of a single batch against filesystem."""

    if logger is None:
        logger = logging.getLogger(__name__)

    batch_id = batch_info['batch_id']
    batch_path = batch_info['base_path']
    ref_version = batch_info['ref_version']

    logger.info(f"🔍 Verifying batch {batch_id}: {batch_info['batch_name']}")
    logger.debug(f"   Batch path: {batch_path}")
    logger.debug(f"   Reference version: {ref_version}")

    # Get protein details from database
    protein_details = get_protein_details_for_batch(conn, batch_id)

    # Filter to representative proteins only if requested
    if not include_non_rep:
        protein_details = [p for p in protein_details if p['is_representative']]

    logger.info(f"   Checking {len(protein_details)} proteins")

    # Find actual files on filesystem
    actual_partition_files, actual_summary_files = find_actual_files(batch_path, ref_version)
    logger.info(f"   Found {len(actual_partition_files)} partition files, {len(actual_summary_files)} summary files")

    # Initialize counters
    db_says_complete = 0
    fs_actually_complete = 0
    db_says_missing = 0
    fs_actually_missing = 0

    files_exist_not_in_db = 0
    files_in_db_not_exist = 0
    files_correctly_tracked = 0

    truly_unclassified = 0
    legitimate_peptides = 0
    files_missing_entirely = 0
    sync_errors = 0

    # Track representative protein accuracy
    db_rep_count = len([p for p in protein_details if p['is_representative']])

    # Track what we find in files
    classified_count = 0
    peptide_count = 0
    unclassified_count = 0

    # Analyze each protein
    for protein in protein_details:
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        protein_key = f"{pdb_id.lower()}_{chain_id}"

        # Database status
        db_complete = (protein['current_stage'] == 'domain_partition_complete' and
                      protein['status'] == 'success')
        db_has_partition_file = protein['db_partition_exists']

        # Filesystem reality
        actual_partition_path = actual_partition_files.get(protein_key)
        actual_summary_path = actual_summary_files.get(protein_key)
        partition_file_exists = actual_partition_path is not None

        if db_complete:
            db_says_complete += 1

        # Check file tracking accuracy
        if db_has_partition_file and not partition_file_exists:
            files_in_db_not_exist += 1
            sync_errors += 1
        elif not db_has_partition_file and partition_file_exists:
            files_exist_not_in_db += 1
            sync_errors += 1
        elif db_has_partition_file and partition_file_exists:
            files_correctly_tracked += 1

        # Analyze file content if it exists
        if partition_file_exists:
            fs_actually_complete += 1

            content_analysis = analyze_partition_file_content(actual_partition_path)

            if content_analysis['is_peptide']:
                legitimate_peptides += 1
                peptide_count += 1
            elif content_analysis['is_classified']:
                classified_count += 1
            elif content_analysis['is_unclassified']:
                truly_unclassified += 1
                unclassified_count += 1
        else:
            files_missing_entirely += 1

    # Check for files that exist but aren't tracked for any protein in this batch
    tracked_proteins = {f"{p['pdb_id'].lower()}_{p['chain_id']}" for p in protein_details}
    untracked_partition_files = set(actual_partition_files.keys()) - tracked_proteins
    files_exist_not_in_db += len(untracked_partition_files)

    if untracked_partition_files:
        logger.debug(f"   Found {len(untracked_partition_files)} untracked partition files")

    # Determine suggestions
    needs_import_run = files_exist_not_in_db > 0 or sync_errors > len(protein_details) * 0.1
    needs_status_sync = abs(db_says_complete - fs_actually_complete) > len(protein_details) * 0.1
    needs_batch_stage_update = False

    # Suggest batch status based on reality
    completion_rate = (fs_actually_complete / len(protein_details)) * 100 if protein_details else 0

    if completion_rate >= 95:
        suggested_batch_status = "completed"
    elif completion_rate >= 80:
        suggested_batch_status = "domain_partition_complete"
    elif completion_rate >= 50:
        suggested_batch_status = "domain_partition_partial"
    else:
        suggested_batch_status = "domain_partition_incomplete"

    # Check if batch status update is needed
    if batch_info['batch_status'] != suggested_batch_status:
        needs_batch_stage_update = True

    logger.info(f"   📊 Results: {fs_actually_complete}/{len(protein_details)} complete "
                f"({completion_rate:.1f}%), {sync_errors} sync errors")
    logger.debug(f"   📈 Classification: {classified_count} classified, "
                 f"{peptide_count} peptides, {unclassified_count} unclassified")

    return FileSystemReality(
        batch_id=batch_id,
        batch_name=batch_info['batch_name'],
        batch_path=batch_path,
        db_says_complete=db_says_complete,
        fs_actually_complete=fs_actually_complete,
        db_says_missing=len(protein_details) - db_says_complete,
        fs_actually_missing=files_missing_entirely,
        files_exist_not_in_db=files_exist_not_in_db,
        files_in_db_not_exist=files_in_db_not_exist,
        files_correctly_tracked=files_correctly_tracked,
        truly_unclassified=truly_unclassified,
        legitimate_peptides=legitimate_peptides,
        files_missing_entirely=files_missing_entirely,
        sync_errors=sync_errors,
        domain_partition_files=len(actual_partition_files),
        domain_summary_files=len(actual_summary_files),
        db_rep_count=db_rep_count,
        fs_rep_count=db_rep_count,  # Assuming rep count is correct for now
        rep_mismatch=0,  # Would need additional analysis
        needs_import_run=needs_import_run,
        needs_status_sync=needs_status_sync,
        needs_batch_stage_update=needs_batch_stage_update,
        suggested_batch_status=suggested_batch_status
    )


def generate_reality_report(reality_checks: List[FileSystemReality], logger):
    """Generate a comprehensive reality check report."""

    if not reality_checks:
        logger.warning("No reality checks to report")
        return

    print("\n" + "="*100)
    print("🔍 FILESYSTEM REALITY CHECK REPORT")
    print("="*100)
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Batches analyzed: {len(reality_checks)}")
    print()

    # Summary statistics
    total_db_complete = sum(r.db_says_complete for r in reality_checks)
    total_fs_complete = sum(r.fs_actually_complete for r in reality_checks)
    total_sync_errors = sum(r.sync_errors for r in reality_checks)
    total_missing = sum(r.files_missing_entirely for r in reality_checks)
    total_untracked = sum(r.files_exist_not_in_db for r in reality_checks)
    total_partition_files = sum(r.domain_partition_files for r in reality_checks)
    total_summary_files = sum(r.domain_summary_files for r in reality_checks)

    print("📊 OVERALL SUMMARY:")
    print(f"   Database says complete: {total_db_complete}")
    print(f"   Filesystem actually complete: {total_fs_complete}")
    print(f"   Difference: {abs(total_db_complete - total_fs_complete)} ({abs(total_db_complete - total_fs_complete)/max(total_db_complete, 1)*100:.1f}%)")
    print(f"   Sync errors detected: {total_sync_errors}")
    print(f"   Files exist but not tracked: {total_untracked}")
    print(f"   Files truly missing: {total_missing}")
    print(f"   Total partition files found: {total_partition_files}")
    print(f"   Total summary files found: {total_summary_files}")
    print()

    # Issues by severity
    critical_batches = [r for r in reality_checks if r.sync_errors > 50 or r.needs_import_run]
    warning_batches = [r for r in reality_checks if r.needs_status_sync or r.needs_batch_stage_update]

    if critical_batches:
        print("🚨 CRITICAL ISSUES (Immediate Action Required):")
        for reality in critical_batches[:10]:
            print(f"   Batch {reality.batch_id}: {reality.batch_name}")
            if reality.sync_errors > 50:
                print(f"     • High sync errors: {reality.sync_errors}")
            if reality.files_exist_not_in_db > 0:
                print(f"     • Untracked files: {reality.files_exist_not_in_db}")
            if reality.needs_import_run:
                print(f"     • Needs import run")
            print()

    # Detailed batch analysis
    print("📋 DETAILED BATCH ANALYSIS:")
    print(f"{'ID':>5} {'Name':25} {'DB Complete':10} {'FS Complete':10} {'Sync Err':8} {'Missing':7} {'Sugg Status':15}")
    print("-" * 100)

    for reality in sorted(reality_checks, key=lambda x: x.sync_errors, reverse=True):
        db_complete = reality.db_says_complete
        fs_complete = reality.fs_actually_complete
        sync_errors = reality.sync_errors
        missing = reality.files_missing_entirely

        suggested = reality.suggested_batch_status[:15]

        print(f"{reality.batch_id:5d} {reality.batch_name[:25]:25} "
              f"{db_complete:10d} {fs_complete:10d} {sync_errors:8d} "
              f"{missing:7d} {suggested:15}")

    print()

    # Classification analysis
    total_peptides = sum(r.legitimate_peptides for r in reality_checks)
    total_unclassified = sum(r.truly_unclassified for r in reality_checks)

    if total_peptides > 0 or total_unclassified > 0:
        print("🧬 CLASSIFICATION ANALYSIS:")
        print(f"   Legitimate peptides (< 30 residues): {total_peptides}")
        print(f"   Truly unclassified (no homologs): {total_unclassified}")
        print(f"   Files missing entirely: {total_missing}")
        print()

    # Action recommendations
    print("🎯 RECOMMENDED ACTIONS:")

    # High priority actions
    high_priority = []
    if total_untracked > 100:
        high_priority.append("1. HIGH: Run import script to capture untracked partition files")
    if total_sync_errors > total_db_complete * 0.2:
        high_priority.append("2. HIGH: Fix status synchronization issues")

    # Medium priority actions
    medium_priority = []
    status_update_needed = sum(1 for r in reality_checks if r.needs_batch_stage_update)
    if status_update_needed > 5:
        medium_priority.append("3. MEDIUM: Update batch statuses based on filesystem reality")

    # Low priority actions
    low_priority = []
    if total_missing > 0:
        low_priority.append("4. LOW: Investigate truly missing files")

    for action in high_priority + medium_priority + low_priority:
        print(f"   {action}")

    print()

    # Success metrics for verified health
    print("✅ VERIFIED HEALTH CRITERIA:")
    print("   Based on filesystem reality, healthy pipeline should have:")
    accuracy = (total_fs_complete / max(total_db_complete, 1)) * 100 if total_db_complete > 0 else 0
    print(f"   • Database accuracy: {accuracy:.1f}% (current)")
    print(f"   • Target accuracy: ≥95%")
    sync_rate = (total_sync_errors/max(total_db_complete, 1)*100) if total_db_complete > 0 else 0
    print(f"   • Sync error rate: <5% (current: {sync_rate:.1f}%)")
    untracked_rate = (total_untracked/max(total_fs_complete, 1)*100) if total_fs_complete > 0 else 0
    print(f"   • Untracked files: <1% (current: {untracked_rate:.1f}%)")

    print("\n" + "="*100)


def export_reality_data(reality_checks: List[FileSystemReality], output_file: str, logger):
    """Export detailed reality check data to JSON."""

    export_data = {
        'metadata': {
            'generated_at': datetime.now().isoformat(),
            'total_batches': len(reality_checks),
            'analysis_type': 'filesystem_reality_check'
        },
        'summary': {
            'total_db_complete': sum(r.db_says_complete for r in reality_checks),
            'total_fs_complete': sum(r.fs_actually_complete for r in reality_checks),
            'total_sync_errors': sum(r.sync_errors for r in reality_checks),
            'total_partition_files': sum(r.domain_partition_files for r in reality_checks),
            'total_summary_files': sum(r.domain_summary_files for r in reality_checks),
            'batches_need_import': sum(1 for r in reality_checks if r.needs_import_run),
            'batches_need_status_sync': sum(1 for r in reality_checks if r.needs_status_sync),
            'batches_need_status_update': sum(1 for r in reality_checks if r.needs_batch_stage_update)
        },
        'batch_details': [asdict(reality) for reality in reality_checks]
    }

    with open(output_file, 'w') as f:
        json.dump(export_data, f, indent=2, default=str)

    logger.info(f"📄 Reality check data exported to: {output_file}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Verify filesystem reality vs database tracking')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--batch-ids', type=int, nargs='+',
                       help='Specific batch IDs to check (default: all)')
    parser.add_argument('--include-non-rep', action='store_true',
                       help='Include non-representative proteins in analysis')
    parser.add_argument('--export-data', help='Export detailed data to JSON file')
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
        # Get batch information
        logger.info("🔍 Loading batch information...")
        batch_info_list = get_batch_filesystem_info(conn, args.batch_ids)
        logger.info(f"Found {len(batch_info_list)} batches to verify")

        # Verify each batch against filesystem reality
        reality_checks = []
        for batch_info in batch_info_list:
            try:
                reality = verify_batch_reality(conn, batch_info, args.include_non_rep, logger)
                reality_checks.append(reality)
            except Exception as e:
                logger.error(f"Error verifying batch {batch_info['batch_id']}: {e}")
                if args.verbose:
                    import traceback
                    logger.error(traceback.format_exc())
                continue

        # Generate comprehensive report
        generate_reality_report(reality_checks, logger)

        # Export data if requested
        if args.export_data:
            export_reality_data(reality_checks, args.export_data, logger)

        # Summary
        logger.info(f"✅ Reality check complete: {len(reality_checks)} batches verified")
        
        critical_issues = sum(1 for r in reality_checks if r.needs_import_run or r.sync_errors > 50)
        if critical_issues > 0:
            logger.warning(f"⚠️  {critical_issues} batches have critical sync issues requiring immediate attention")
        else:
            logger.info("🎉 No critical sync issues detected")

    except Exception as e:
        logger.error(f"Reality check failed: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
