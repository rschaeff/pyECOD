#!/usr/bin/env python3
"""
validate_batch_files.py - Validate and synchronize database records with filesystem

This script checks if files recorded in the database actually exist on the filesystem,
and updates the database records accordingly.
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from typing import Dict, Any, Optional, List, Tuple

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
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

def validate_file_exists(file_path: str, base_path: str, min_size: int = 10) -> bool:
    """
    Validate that a file exists and has a minimum size
    
    Args:
        file_path: Relative or absolute path to the file
        base_path: Base path to resolve relative paths
        min_size: Minimum file size in bytes
        
    Returns:
        True if file exists and has at least min_size bytes
    """
    # Handle absolute and relative paths
    if os.path.isabs(file_path):
        full_path = file_path
    else:
        full_path = os.path.join(base_path, file_path)
    
    # Check if file exists and has minimum size
    if os.path.exists(full_path) and os.path.getsize(full_path) >= min_size:
        return True
    return False

def check_xml_validity(file_path: str, base_path: str) -> bool:
    """
    Check if a file contains valid XML
    
    Args:
        file_path: Path to the XML file
        base_path: Base path to resolve relative paths
        
    Returns:
        True if file contains valid XML
    """
    # Handle absolute and relative paths
    if os.path.isabs(file_path):
        full_path = file_path
    else:
        full_path = os.path.join(base_path, file_path)
    
    # Check file existence first
    if not os.path.exists(full_path):
        return False
    
    # Try to parse XML
    try:
        tree = ET.parse(full_path)
        root = tree.getroot()
        return True
    except:
        return False

def get_batch_info(context, batch_id: int) -> Dict[str, Any]:
    """Get batch information from database"""
    query = """
    SELECT 
        id, batch_name, base_path, ref_version,
        total_items, completed_items, status
    FROM 
        ecod_schema.batch
    WHERE 
        id = %s
    """
    
    result = context.db.execute_query(query, (batch_id,))
    if not result:
        return {}
    
    return {
        'id': result[0][0],
        'name': result[0][1],
        'path': result[0][2],
        'reference': result[0][3],
        'total_items': result[0][4],
        'completed_items': result[0][5],
        'status': result[0][6]
    }

def validate_domain_summaries(context, batch_id: int, dry_run: bool = True, 
                            fix_errors: bool = False, limit: int = None,
                            xml_check: bool = False) -> Tuple[int, int, int]:
    """
    Validate domain summary files and update database if needed
    
    Args:
        context: Application context
        batch_id: Batch ID to validate
        dry_run: If True, don't make database changes
        fix_errors: If True, fix database errors for missing files
        limit: Maximum number of files to check
        xml_check: If True, validate XML content
        
    Returns:
        Tuple of (total files, valid files, fixed files)
    """
    logger = logging.getLogger("ecod.validate")
    
    # Get batch information
    batch_info = get_batch_info(context, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return 0, 0, 0
    
    base_path = batch_info['path']
    logger.info(f"Validating batch {batch_id} ({batch_info['name']}) with base path: {base_path}")
    
    # Query for domain summary files
    query = """
    SELECT 
        pf.id, pf.process_id, pf.file_path, pf.file_exists, pf.file_size,
        p.pdb_id, p.chain_id
    FROM 
        ecod_schema.process_file pf
    JOIN 
        ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
        AND pf.file_type = 'domain_summary'
    """
    
    if limit:
        query += f" LIMIT {limit}"
    
    results = context.db.execute_query(query, (batch_id,))
    
    # Track stats
    total_files = len(results)
    valid_files = 0
    fixed_files = 0
    invalid_files = []
    
    # Process each file
    for row in results:
        file_id = row[0]
        process_id = row[1]
        file_path = row[2]
        file_exists_db = row[3]
        file_size_db = row[4]
        pdb_id = row[5]
        chain_id = row[6]
        
        # Check if file actually exists with minimum size
        file_exists_fs = validate_file_exists(file_path, base_path, min_size=10)
        
        # Check XML validity if requested
        xml_valid = True
        if file_exists_fs and xml_check:
            xml_valid = check_xml_validity(file_path, base_path)
        
        # Determine if file is valid
        is_valid = file_exists_fs and xml_valid
        
        if is_valid:
            valid_files += 1
            
            # Update database if file exists but database says it doesn't
            if not file_exists_db and fix_errors and not dry_run:
                try:
                    # Get actual file size
                    full_path = os.path.join(base_path, file_path) if not os.path.isabs(file_path) else file_path
                    file_size = os.path.getsize(full_path)
                    
                    # Update database
                    update_query = """
                    UPDATE ecod_schema.process_file
                    SET file_exists = TRUE, file_size = %s, last_checked = NOW()
                    WHERE id = %s
                    """
                    
                    context.db.execute_query(update_query, (file_size, file_id))
                    fixed_files += 1
                    logger.info(f"Updated database record for {pdb_id}_{chain_id}: marked as exists")
                except Exception as e:
                    logger.error(f"Error updating database for {pdb_id}_{chain_id}: {str(e)}")
        else:
            # Track invalid files
            invalid_files.append({
                'file_id': file_id,
                'process_id': process_id,
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'file_path': file_path,
                'exists_db': file_exists_db,
                'exists_fs': file_exists_fs,
                'xml_valid': xml_valid if file_exists_fs else False
            })
            
            # Fix database if requested
            if file_exists_db and fix_errors and not dry_run:
                try:
                    # Update database to mark file as not existing
                    update_query = """
                    UPDATE ecod_schema.process_file
                    SET file_exists = FALSE, file_size = 0, last_checked = NOW()
                    WHERE id = %s
                    """
                    
                    context.db.execute_query(update_query, (file_id,))
                    fixed_files += 1
                    logger.info(f"Updated database record for {pdb_id}_{chain_id}: marked as not exists")
                    
                    # Also update process status
                    status_query = """
                    UPDATE ecod_schema.process_status
                    SET status = 'pending', current_stage = 'domain_summary'
                    WHERE id = %s
                    """
                    
                    context.db.execute_query(status_query, (process_id,))
                    logger.info(f"Reset process status for {pdb_id}_{chain_id}")
                except Exception as e:
                    logger.error(f"Error updating database for {pdb_id}_{chain_id}: {str(e)}")
    
    # Log summary
    logger.info(f"Validation summary for domain summaries:")
    logger.info(f"  Total files: {total_files}")
    logger.info(f"  Valid files: {valid_files} ({valid_files/total_files*100:.1f}%)")
    logger.info(f"  Invalid files: {len(invalid_files)} ({len(invalid_files)/total_files*100:.1f}%)")
    
    if fix_errors and not dry_run:
        logger.info(f"  Fixed database records: {fixed_files}")
    
    # Log batch status
    logger.info(f"Batch status: {batch_info['status']}")
    logger.info(f"Batch progress: {batch_info['completed_items']}/{batch_info['total_items']} completed")
    
    # Provide more details on invalid files
    if invalid_files:
        logger.info("First 10 invalid files:")
        for i, file in enumerate(invalid_files[:10]):
            issues = []
            if not file['exists_fs']:
                issues.append("missing on filesystem")
            elif not file['xml_valid']:
                issues.append("invalid XML")
            
            logger.info(f"  {i+1}. {file['pdb_id']}_{file['chain_id']}: {', '.join(issues)}")
    
    # Update batch status if fixed records
    if fixed_files > 0 and fix_errors and not dry_run:
        # Recalculate completed items
        count_query = """
        SELECT COUNT(*)
        FROM ecod_schema.process_status
        WHERE batch_id = %s AND status = 'success' AND current_stage = 'domain_partition_complete'
        """
        
        count_result = context.db.execute_query(count_query, (batch_id,))
        if count_result:
            completed_count = count_result[0][0]
            
            # Update batch status
            update_query = """
            UPDATE ecod_schema.batch
            SET completed_items = %s,
                status = CASE 
                    WHEN %s >= total_items THEN 'completed' 
                    ELSE 'processing' 
                END
            WHERE id = %s
            """
            
            context.db.execute_query(update_query, (completed_count, completed_count, batch_id))
            logger.info(f"Updated batch status: {completed_count}/{batch_info['total_items']} completed")
    
    return total_files, valid_files, fixed_files

def main():
    """Main function to validate batch files"""
    parser = argparse.ArgumentParser(description='Validate batch files and sync database')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to validate')
    parser.add_argument('--dry-run', action='store_true',
                      help='Check files but don\'t update database')
    parser.add_argument('--fix', action='store_true',
                      help='Fix database records for missing files')
    parser.add_argument('--xml-check', action='store_true',
                      help='Perform XML structure validation')
    parser.add_argument('--limit', type=int,
                      help='Limit number of files to check')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.validate")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Validate domain summary files
    total, valid, fixed = validate_domain_summaries(
        context, 
        args.batch_id, 
        args.dry_run, 
        args.fix, 
        args.limit,
        args.xml_check
    )
    
    # Summary
    if args.dry_run:
        logger.info("This was a dry run - no database changes were made")
        if args.fix:
            logger.info("Re-run with --fix without --dry-run to fix database records")
    
    if valid == total:
        logger.info("All files are valid and database is in sync")
        return 0
    else:
        logger.warning(f"Found {total - valid} invalid files")
        if not args.fix or args.dry_run:
            logger.info("Run with --fix to update database records")
        return 1

if __name__ == "__main__":
    sys.exit(main())