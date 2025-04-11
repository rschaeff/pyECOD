#!/usr/bin/env python3
"""
fix_chain_blast_paths.py - Fix path duplication issues in chain blast file records
"""

import os
import sys
import logging
import argparse
from typing import Dict, Any, Optional, List

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

def fix_chain_blast_paths(context: Any, batch_id: int, dry_run: bool = False) -> int:
    """
    Fix path duplication issues in chain blast file records
    
    Args:
        context: Application context
        batch_id: Batch ID to fix
        dry_run: Whether to show changes without applying them
        
    Returns:
        Number of records fixed
    """
    logger = logging.getLogger("ecod.fix_paths")
    
    # Get batch information
    batch_query = "SELECT batch_name, base_path FROM ecod_schema.batch WHERE id = %s"
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch ID {batch_id} not found")
        return 0
    
    batch_name = batch_result[0][0]
    batch_base_path = batch_result[0][1]
    
    logger.info(f"Fixing chain blast paths for batch: {batch_name} (ID: {batch_id})")
    
    # Find all chain_blast_result file records for this batch
    file_query = """
    SELECT 
        pf.id, pf.file_path, p.pdb_id, p.chain_id
    FROM 
        ecod_schema.process_file pf
    JOIN 
        ecod_schema.process_status ps ON pf.process_id = ps.id
    JOIN 
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s AND pf.file_type = 'chain_blast_result'
    """
    
    file_results = context.db.execute_query(file_query, (batch_id,))
    
    if not file_results:
        logger.warning(f"No chain_blast_result files found for batch {batch_id}")
        return 0
    
    fixed_count = 0
    
    # Process each file record
    for file_result in file_results:
        file_id = file_result[0]
        file_path = file_result[1]
        pdb_id = file_result[2]
        chain_id = file_result[3]
        
        # Check if path has duplication issue
        # Look for patterns like /base/path/../../../../base/path/actual/path
        if "/../" in file_path:
            logger.info(f"Found path with duplication issue: {file_path}")
            
            # Find the correct batch subdirectory for this file
            # Check all possible batch_N directories
            pdb_chain = f"{pdb_id}_{chain_id}"
            possible_batch_dirs = ["batch_0", "batch_1", "batch_2"]
            
            # Construct the expected correct path
            corrected_path = None
            
            for batch_dir in possible_batch_dirs:
                # Try to find the file in this batch directory
                test_path = f"blast/chain/{batch_dir}/{pdb_chain}.chainwise_blast.xml"
                full_test_path = os.path.join(batch_base_path, test_path)
                
                if os.path.exists(full_test_path):
                    corrected_path = test_path
                    break
            
            # If not found with standard name, try alternative file naming patterns
            if not corrected_path:
                for batch_dir in possible_batch_dirs:
                    # Check for alternative file names
                    alt_names = [
                        f"{pdb_chain}.blast",
                        f"{pdb_chain}_blast.txt",
                        f"{pdb_chain}.xml"
                    ]
                    
                    for alt_name in alt_names:
                        test_path = f"blast/chain/{batch_dir}/{alt_name}"
                        full_test_path = os.path.join(batch_base_path, test_path)
                        
                        if os.path.exists(full_test_path):
                            corrected_path = test_path
                            break
                    
                    if corrected_path:
                        break
            
            # If we found a corrected path, update the record
            if corrected_path:
                logger.info(f"Found correct path: {corrected_path}")
                
                if not dry_run:
                    # Update the record
                    update_query = """
                    UPDATE ecod_schema.process_file
                    SET file_path = %s, file_exists = TRUE, last_checked = NOW()
                    WHERE id = %s
                    """
                    
                    context.db.execute_query(update_query, (corrected_path, file_id))
                    logger.info(f"Updated file path for ID {file_id}")
                else:
                    logger.info(f"Would update file path for ID {file_id} to {corrected_path}")
                
                fixed_count += 1
            else:
                # If we couldn't find the file, set exists to FALSE
                logger.warning(f"Could not find correct path for {pdb_chain}")
                
                if not dry_run:
                    # Update the record to mark as not existing
                    update_query = """
                    UPDATE ecod_schema.process_file
                    SET file_exists = FALSE, last_checked = NOW()
                    WHERE id = %s
                    """
                    
                    context.db.execute_query(update_query, (file_id,))
                    logger.info(f"Updated file existence flag for ID {file_id} to FALSE")
                else:
                    logger.info(f"Would update file existence flag for ID {file_id} to FALSE")
    
    if not dry_run:
        logger.info(f"Fixed {fixed_count} file paths out of {len(file_results)} records")
    else:
        logger.info(f"Would fix {fixed_count} file paths out of {len(file_results)} records")
    
    return fixed_count

def main():
    """Main function to fix chain blast paths"""
    parser = argparse.ArgumentParser(description='Fix path duplication issues in chain blast file records')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to fix')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show changes without applying them')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Fix chain blast paths
    fixed_count = fix_chain_blast_paths(context, args.batch_id, args.dry_run)
    
    return 0 if fixed_count > 0 else 1

if __name__ == "__main__":
    sys.exit(main())