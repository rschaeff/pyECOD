#!/usr/bin/env python3
"""
Filesystem Reality Checker

This script verifies the actual state of domain partition files on the filesystem
vs what the database tracking thinks is happening. This is critical for establishing
accurate health baselines.

Usage:
    python filesystem_reality_checker.py --config config.yml [options]

Verification Areas:
1. File existence vs database records
2. Batch completion based on filesystem reality
3. Missing partition analysis (truly missing vs sync issues)
4. Peptide vs unclassified vs truly missing domains
5. Representative vs non-representative protein accuracy
6. Batch stage verification
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
                   AND pf.file_type = 'domain_partition' THEN 1 END) as db_partition_files
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
        p.source_id = (pp.pdb_id || '_' || pp.chain_id) AND pp.batch_id = ps.batch_id
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


def find_actual_partition_files(batch_path, ref_version):
    """Find all actual partition files on the filesystem."""
    
    # Common patterns for partition files
    partition_patterns = [
        f"{batch_path}/domains/*/{ref_version}/*_domains.xml",
        f"{batch_path}/domains/*/*_domains.xml",
        f"{batch_path}/partition/*/{ref_version}/*_partition.xml",
        f"{batch_path}/partition/*/*_partition.xml",
        f"{batch_path}/**/domains.xml",
        f"{batch_path}/**/*_domains.xml"
    ]
    
    found_files = {}
    
    for pattern in partition_patterns:
        files = glob.glob(pattern, recursive=True)
        for file_path in files:
            # Extract protein identifier from filename/path
            basename = os.path.basename(file_path)
            
            # Try different naming patterns
            if '_domains.xml' in basename:
                protein_id = basename.replace('_domains.xml', '')
            elif '_partition.xml' in basename:
                protein_id = basename.replace('_partition.xml', '')
            elif basename == 'domains.xml':
                # Look in parent directory for protein ID
                parent_dir = os.path.basename(os.path.dirname(file_path))
                protein_id = parent_dir
            else:
                continue
            
            # Normalize protein ID (handle various formats)
            if '_' in protein_id and len(protein_id.split('_')) == 2:
                pdb_id, chain_id = protein_id.split('_')
                normalized_id = f"{pdb_id.lower()}_{chain_id}"
                found_files[normalized_id] = file_path
    
    return found_files


def analyze_partition_file_content(file_path):
    """Analyze the content of a partition file to understand classification status."""
    
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Look for domain elements
        domains = root.findall('.//domain')
        
        # Get sequence length if available
        sequence_elem = root.find('.//sequence')
        sequence_length = 0
        if sequence_elem is not None and sequence_elem.text:
            sequence_length = len(sequence_elem.text.replace('\n', '').replace(' ', ''))
        
        # Check for length attribute in sequence tag
        if sequence_length == 0 and sequence_elem is not None:
            length_attr = sequence_elem.get('length')
            if length_attr:
                try:
                    sequence_length = int(length_attr)
                except ValueError:
                    pass
        
        # Determine classification status
        is_peptide = sequence_length > 0 and sequence_length < 30
        has_domains = len(domains) > 0
        is_classified = has_domains and not is_peptide
        is_unclassified = not has_domains and not is_peptide
        
        return {
            'domain_count': len(domains),
            'sequence_length': sequence_length,
            'is_peptide': is_peptide,
            'is_classified': is_classified,
            'is_unclassified': is_unclassified,
            'file_valid': True
        }
        
    except Exception as e:
        return {
            'domain_count': 0,
            'sequence_length': 0,
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
    
    # Get protein details from database
    protein_details = get_protein_details_for_batch(conn, batch_id)
    
    # Filter to representative proteins only if requested
    if not include_non_rep:
        protein_details = [p for p in protein_details if p['is_representative']]
    
    logger.info(f"   Checking {len(protein_details)} proteins")
    
    # Find actual partition files on filesystem
    actual_files = find_actual_partition_files(batch_path, ref_version)
    logger.info(f"   Found {len(actual_files)} actual partition files")
    
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
        actual_file_path = actual_files.get(protein_key)
        file_actually_exists = actual_file_path is not None
        
        if db_complete:
            db_says_complete += 1
        
        if db_has_partition_file and not file_actually_exists:
            files_in_db_not_exist += 1
            sync_errors += 1
        elif not db_has_partition_file and file_actually_exists:
            files_exist_not_in_db += 1
            sync_errors += 1
        elif db_has_partition_file and file_actually_exists:
            files_correctly_tracked += 1
        
        # Analyze file content if it exists
        if file_actually_exists:
            fs_actually_complete += 1
            
            content_analysis = analyze_partition_file_content(actual_file_path)
            
            if content_analysis['is_peptide']:
                legitimate_peptides += 1
            elif content_analysis['is_unclassified']:
                truly_unclassified += 1
        else:
            files_missing_entirely += 1
    
    # Check for files that exist but aren't tracked for any protein
    tracked_proteins = {f"{p['pdb_id'].lower()}_{p['chain_id']}" for p in protein_details}
    untracked_files = set(actual_files.keys()) - tracked_proteins
    files_exist_not_in_db += len(untracked_files)
    
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
    
    print("📊 OVERALL SUMMARY:")
    print(f"   Database says complete: {total_db_complete}")
    print(f"   Filesystem actually complete: {total_fs_complete}")
    print(f"   Difference: {abs(total_db_complete - total_fs_complete)} ({abs(total_db_complete - total_fs_complete)/max(total_db_complete, 1)*100:.1f}%)")
    print(f"   Sync errors detected: {total_sync_errors}")
    print(f"   Files exist but not tracked: {total_untracked}")
    print(f"   Files truly missing: {total_missing}")
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
    print(f"{'ID':>5} {'Name':20} {'DB Complete':10} {'FS Complete':10} {'Sync Err':8} {'Missing':7} {'Status':12} {'Suggested':12}")
    print("-" * 100)
    
    for reality in sorted(reality_checks, key=lambda x: x.sync_errors, reverse=True):
        db_complete = reality.db_says_complete
        fs_complete = reality.fs_actually_complete
        sync_errors = reality.sync_errors
        missing = reality.files_missing_entirely
        
        # Current batch status (first 12 chars)
        current_status = reality.batch_name.split('_')[-1] if '_' in reality.batch_name else "unknown"
        current_status = current_status[:12]
        
        suggested = reality.suggested_batch_status[:12]
        
        print(f"{reality.batch_id:5d} {reality.batch_name[:20]:20} "
              f"{db_complete:10d} {fs_complete:10d} {sync_errors:8d} "
              f"{missing:7d} {current_status:12} {suggested:12}")
    
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
    
    # Specific commands
    if critical_batches or warning_batches:
        print("🔧 SPECIFIC REPAIR COMMANDS:")
        
        if total_untracked > 0:
            print("   # Import untracked partition files")
            print("   python import_domain_partitions_with_repair.py --config config.yml --repair-mode --conflict-strategy update")
        
        if total_sync_errors > 0:
            print("   # Fix status synchronization")
            print("   python fix_status_synchronization.py --config config.yml --fix-non-rep --update-batch-status")
        
        if status_update_needed > 0:
            print("   # Update batch statuses")
            print("   python update_batch_statuses.py --config config.yml --based-on-filesystem")
        
        print()
    
    # Success metrics for verified health
    print("✅ VERIFIED HEALTH CRITERIA:")
    print("   Based on filesystem reality, healthy pipeline should have:")
    accuracy = (total_fs_complete / max(total_db_complete, 1)) * 100
    print(f"   • Database accuracy: {accuracy:.1f}% (current)")
    print(f"   • Target accuracy: ≥95%")
    print(f"   • Sync error rate: <5% (current: {total_sync_errors/max(total_db_complete, 1)*100:.1f}%)")
    print(f"   • Untracked files: <1% (current: {total_untracked/max(total_fs_complete, 1)*100:.1f}%)")
    
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
