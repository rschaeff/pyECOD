#!/usr/bin/env python3
"""
validate_batch_files.py - Validate and synchronize database records with filesystem

This script checks if files recorded in the database actually exist on the filesystem,
and updates the database records accordingly.

Fixed version with enhanced debugging and explicit transaction control.
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

def fix_single_file_record(context, file_info, exists_on_fs, logger, dry_run=False):
    """
    Fix a single file record in the database
    
    Returns:
        True if record was fixed, False otherwise
    """
    file_id = file_info['file_id']
    process_id = file_info['process_id']
    pdb_id = file_info['pdb_id']
    chain_id = file_info['chain_id']
    exists_in_db = file_info['exists_db']
    file_path = file_info['file_path']
    base_path = file_info['base_path']
    
    # Determine if record needs fixing
    needs_update = (exists_in_db != exists_on_fs)
    
    # Enhanced debugging - print to console for immediate visibility
    print(f"DEBUG: Processing {pdb_id}_{chain_id}: exists_in_db={exists_in_db}, exists_on_fs={exists_on_fs}, needs_update={needs_update}, dry_run={dry_run}")
    logger.debug(f"Processing {pdb_id}_{chain_id}: exists_in_db={exists_in_db}, exists_on_fs={exists_on_fs}, needs_update={needs_update}, dry_run={dry_run}")
    
    if not needs_update:
        print(f"DEBUG: No update needed for {pdb_id}_{chain_id}")
        return False
        
    if dry_run:
        if exists_in_db and not exists_on_fs:
            print(f"DEBUG: Would update {pdb_id}_{chain_id} to not exist (dry run)")
            logger.debug(f"Would update {pdb_id}_{chain_id} to not exist (dry run)")
        elif not exists_in_db and exists_on_fs:
            print(f"DEBUG: Would update {pdb_id}_{chain_id} to exist (dry run)")
            logger.debug(f"Would update {pdb_id}_{chain_id} to exist (dry run)")
        return False
    
    # Log that we're going to make changes
    print(f"DEBUG: ATTEMPTING DATABASE UPDATE for {pdb_id}_{chain_id}")
    logger.debug(f"Attempting to fix record for {pdb_id}_{chain_id} (dry_run={dry_run})")
    
    try:
        if exists_in_db and not exists_on_fs:
            # Update database to mark file as not existing
            print(f"DEBUG: Updating {pdb_id}_{chain_id} to NOT exist")
            update_query = """
            UPDATE ecod_schema.process_file
            SET file_exists = FALSE, file_size = 0, last_checked = NOW()
            WHERE id = %s
            """
            
            # Use direct cursor for better control
            cursor = context.db.connection.cursor()
            cursor.execute(update_query, (file_id,))
            rows_affected = cursor.rowcount
            print(f"DEBUG: UPDATE query affected {rows_affected} rows for file_id {file_id}")
            logger.debug(f"UPDATE process_file query affected {rows_affected} rows for file_id {file_id}")
            
            # Explicitly commit
            context.db.connection.commit()
            print(f"DEBUG: Transaction committed for {pdb_id}_{chain_id}")
            cursor.close()
            
            if rows_affected > 0:
                print(f"DEBUG: Successfully updated {pdb_id}_{chain_id} as not existing")
                logger.info(f"Updated database record for {pdb_id}_{chain_id}: marked as not exists")
            else:
                print(f"DEBUG: No rows affected for {pdb_id}_{chain_id}")
                logger.warning(f"No rows affected when updating file record for {pdb_id}_{chain_id}")
            
            # Also update process status
            print(f"DEBUG: Updating process status for {pdb_id}_{chain_id}")
            status_query = """
            UPDATE ecod_schema.process_status
            SET status = 'pending', current_stage = 'domain_summary'
            WHERE id = %s
            """
            
            cursor = context.db.connection.cursor()
            cursor.execute(status_query, (process_id,))
            rows_affected = cursor.rowcount
            print(f"DEBUG: UPDATE process_status query affected {rows_affected} rows for process_id {process_id}")
            logger.debug(f"UPDATE process_status query affected {rows_affected} rows for process_id {process_id}")
            
            # Explicitly commit again
            context.db.connection.commit()
            cursor.close()
            
            if rows_affected > 0:
                print(f"DEBUG: Successfully reset process status for {pdb_id}_{chain_id}")
                logger.info(f"Reset process status for {pdb_id}_{chain_id}")
            else:
                print(f"DEBUG: No rows affected for process status of {pdb_id}_{chain_id}")
                logger.warning(f"No rows affected when updating process status for {pdb_id}_{chain_id}")
            
            return rows_affected > 0
            
        elif not exists_in_db and exists_on_fs:
            # Update database to mark file as existing
            print(f"DEBUG: Updating {pdb_id}_{chain_id} to EXIST")
            full_path = os.path.join(base_path, file_path) if not os.path.isabs(file_path) else file_path
            file_size = os.path.getsize(full_path)
            
            update_query = """
            UPDATE ecod_schema.process_file
            SET file_exists = TRUE, file_size = %s, last_checked = NOW()
            WHERE id = %s
            """
            
            # Use direct cursor for better control
            cursor = context.db.connection.cursor()
            cursor.execute(update_query, (file_size, file_id))
            rows_affected = cursor.rowcount
            print(f"DEBUG: UPDATE query affected {rows_affected} rows for file_id {file_id}")
            
            # Explicitly commit
            context.db.connection.commit()
            cursor.close()
            
            if rows_affected > 0:
                print(f"DEBUG: Successfully updated {pdb_id}_{chain_id} as existing")
                logger.info(f"Updated database record for {pdb_id}_{chain_id}: marked as exists")
            else:
                print(f"DEBUG: No rows affected for {pdb_id}_{chain_id}")
                logger.warning(f"No rows affected when updating file record for {pdb_id}_{chain_id}")
            
            return rows_affected > 0
    except Exception as e:
        print(f"DEBUG: ERROR fixing record for {pdb_id}_{chain_id}: {str(e)}")
        logger.error(f"Error fixing database record for {pdb_id}_{chain_id}: {str(e)}", exc_info=True)
        return False
        
    return False

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
    
    # Add debug logging for function parameters
    print(f"DEBUG: validate_domain_summaries called with: batch_id={batch_id}, dry_run={dry_run}, fix_errors={fix_errors}")
    logger.debug(f"validate_domain_summaries called with: batch_id={batch_id}, dry_run={dry_run}, fix_errors={fix_errors}")
    
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
    ORDER BY pf.id  -- Order by ID for consistent results
    """
    
    if limit:
        query += f" LIMIT {limit}"
    
    print(f"DEBUG: Executing query: {query}")
    print(f"DEBUG: With parameters: ({batch_id},)")
    
    results = context.db.execute_query(query, (batch_id,))
    print(f"DEBUG: Query returned {len(results)} rows")
    
    # CRITICAL CHECK - Print first few rows to verify data structure
    if results:
        print(f"DEBUG: First row sample: {results[0]}")
        if len(results) > 1:
            print(f"DEBUG: Second row sample: {results[1]}")
    else:
        print("DEBUG: NO RESULTS RETURNED FROM QUERY")
        logger.error("Query returned no results - database may be empty or query incorrect")
        return 0, 0, 0
    
    # Track stats
    total_files = len(results)
    valid_files = 0
    fixed_files = 0
    invalid_files = []
    
    # Add additional debugging
    print(f"DEBUG: Processing {total_files} file records with fix_errors={fix_errors} and dry_run={dry_run}")
    logger.debug(f"Processing {total_files} file records with fix_errors={fix_errors} and dry_run={dry_run}")
    
    # Process counter to track loop progress
    processed_count = 0
    
    # Process each file
    for row in results:
        processed_count += 1
        
        # Print progress for every 100 files
        if processed_count % 100 == 0 or processed_count == 1:
            print(f"DEBUG: Processing file {processed_count} of {total_files}")
        
        file_id = row[0]
        process_id = row[1]
        file_path = row[2]
        file_exists_db = row[3]
        file_size_db = row[4]
        pdb_id = row[5]
        chain_id = row[6]
        
        # First item - print detailed debug info
        if processed_count == 1:
            print(f"DEBUG: Processing first item: file_id={file_id}, process_id={process_id}, pdb_id={pdb_id}, chain_id={chain_id}")
            print(f"DEBUG: File path: {file_path}")
            print(f"DEBUG: File exists in DB: {file_exists_db}")
        
        # Check if file actually exists with minimum size
        file_exists_fs = validate_file_exists(file_path, base_path, min_size=10)
        
        # For first file, print detailed path info
        if processed_count == 1:
            full_path = os.path.join(base_path, file_path) if not os.path.isabs(file_path) else file_path
            print(f"DEBUG: Checking file exists at: {full_path}")
            print(f"DEBUG: File exists on filesystem: {file_exists_fs}")
            print(f"DEBUG: Base path exists: {os.path.exists(base_path)}")
            if not os.path.exists(base_path):
                logger.error(f"Base path does not exist: {base_path}")
        
        # Check XML validity if requested
        xml_valid = True
        if file_exists_fs and xml_check:
            xml_valid = check_xml_validity(file_path, base_path)
        
        # Determine if file is valid
        is_valid = file_exists_fs and xml_valid
        
        if is_valid:
            valid_files += 1
            
            # Update database if file exists but database says it doesn't
            if fix_errors and not file_exists_db:
                print(f"DEBUG: Valid file {pdb_id}_{chain_id} needs DB update (exists on FS but not in DB)")
                file_info = {
                    'file_id': file_id,
                    'process_id': process_id,
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'file_path': file_path,
                    'exists_db': file_exists_db,
                    'base_path': base_path
                }
                
                if fix_single_file_record(context, file_info, True, logger, dry_run):
                    fixed_files += 1
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
                'xml_valid': xml_valid if file_exists_fs else False,
                'base_path': base_path
            })
            
            # First few invalid files - print detailed debug info
            if len(invalid_files) <= 5:
                print(f"DEBUG: Invalid file {pdb_id}_{chain_id}: exists_in_db={file_exists_db}, exists_on_fs={file_exists_fs}")
                
            # Fix database if requested
            if fix_errors and file_exists_db:
                print(f"DEBUG: Invalid file {pdb_id}_{chain_id} needs DB update (exists in DB but not on FS)")
                if fix_single_file_record(context, invalid_files[-1], False, logger, dry_run):
                    fixed_files += 1
                    
                    # Detailed debugging for first few fixed files
                    if fixed_files <= 5:
                        print(f"DEBUG: Successfully fixed record for {pdb_id}_{chain_id}")
    
    print(f"DEBUG: Finished processing all {processed_count} files")
    print(f"DEBUG: Results: valid={valid_files}, invalid={len(invalid_files)}, fixed={fixed_files}")

    
    # Log summary
    logger.info(f"Validation summary for domain summaries:")
    logger.info(f"  Total files: {total_files}")
    logger.info(f"  Valid files: {valid_files} ({valid_files/total_files*100:.1f}%)")
    logger.info(f"  Invalid files: {len(invalid_files)} ({len(invalid_files)/total_files*100:.1f}%)")
    
    if fix_errors:
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
        try:
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
                
                rows_affected = context.db.execute_query_with_rowcount(update_query, (completed_count, completed_count, batch_id))
                logger.info(f"Updated batch status: {rows_affected} rows affected, {completed_count}/{batch_info['total_items']} completed")
        except Exception as e:
            logger.error(f"Error updating batch status: {str(e)}", exc_info=True)
    
    return total_files, valid_files, fixed_files

def main():
    """Main function to validate batch files"""
    parser = argparse.ArgumentParser(description='Validate batch files and sync database')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to validate')
    parser.add_argument('--dry-run', action='store_true',
                      help='Check files but don\'t update database (default is to fix)')
    parser.add_argument('--xml-check', action='store_true',
                      help='Perform XML structure validation')
    parser.add_argument('--limit', type=int,
                      help='Limit number of files to check')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--debug-db', action='store_true',
                      help='Show SQL queries and detailed database info')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.validate")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Add execute_query_with_rowcount method to database manager if not exists
    if not hasattr(context.db, 'execute_query_with_rowcount'):
        def execute_query_with_rowcount(sql, params=None):
            """Execute query and return number of affected rows"""
            cursor = context.db.connection.cursor()
            cursor.execute(sql, params)
            rowcount = cursor.rowcount
            context.db.connection.commit()
            cursor.close()
            return rowcount
        
        # Add method to database manager
        context.db.execute_query_with_rowcount = execute_query_with_rowcount
    
    # Enable database debugging if requested
    if args.debug_db:
        logger.info("Enabling detailed database logging")
        
        # Monkey patch to log queries
        original_execute = context.db.execute_query
        
        def debug_execute(sql, params=None):
            logger.debug(f"Executing SQL: {sql}")
            logger.debug(f"Parameters: {params}")
            result = original_execute(sql, params)
            logger.debug(f"Result: {result}")
            return result
        
        context.db.execute_query = debug_execute
    
    # Simplified logic: dry_run comes directly from args
    # By default, we fix records unless dry_run is specified
    dry_run = args.dry_run
    
    # Always attempt to fix records unless in dry_run mode
    fix_errors = not dry_run
    
    logger.debug(f"Running with dry_run={dry_run}, fix_errors={fix_errors}")
    
    total, valid, fixed = validate_domain_summaries(
        context, 
        args.batch_id, 
        dry_run,
        fix_errors,
        args.limit,
        args.xml_check
    )
    
    # Summary
    if dry_run:
        logger.info("This was a dry run - no database changes were made")
        logger.info("Run without --dry-run to fix database records")
    else:
        logger.info(f"Database records updated: {fixed} files fixed")
    
    if valid == total:
        logger.info("All files are valid and database is in sync")
        return 0
    else:
        logger.warning(f"Found {total - valid} invalid files")
        if dry_run:
            logger.info("Run without --dry-run to fix database records")
        return 1

if __name__ == "__main__":
    sys.exit(main())