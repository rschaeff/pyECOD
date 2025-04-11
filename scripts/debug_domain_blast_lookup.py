#!/usr/bin/env python3
"""
debug_domain_blast_lookup.py - Debug and fix domain_blast_result file lookup issues
"""

import os
import sys
import logging
import argparse
from typing import Dict, Any, Optional

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

def main():
    """Debug and fix domain_blast_result file lookup"""
    parser = argparse.ArgumentParser(description='Debug domain_blast_result file lookup')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to debug')
    parser.add_argument('--protein-id', type=int, required=True,
                      help='Protein ID to debug')
    parser.add_argument('--apply-fix', action='store_true',
                      help='Apply fix for domain_blast_result file path')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.debug_blast_lookup")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get protein and batch information
    query = """
    SELECT 
        p.id, p.source_id, p.pdb_id, p.chain_id,
        b.id as batch_id, b.batch_name, b.base_path, b.ref_version,
        ps.id as process_id
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    WHERE 
        p.id = %s AND b.id = %s
    """
    
    result = context.db.execute_query(query, (args.protein_id, args.batch_id))
    
    if not result:
        logger.error(f"Protein {args.protein_id} not found in batch {args.batch_id}")
        return 1
    
    protein_info = result[0]
    pdb_id = protein_info[2]
    chain_id = protein_info[3]
    batch_path = protein_info[6]
    process_id = protein_info[8]
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    logger.info(f"Debugging domain_blast_result lookup for protein: {pdb_chain} (ID: {args.protein_id})")
    
    # Get the domain_blast_result file record
    file_query = """
    SELECT 
        id, file_type, file_path, file_exists
    FROM 
        ecod_schema.process_file
    WHERE 
        process_id = %s AND file_type = 'domain_blast_result'
    """
    
    file_result = context.db.execute_query(file_query, (process_id,))
    
    if file_result:
        file_id = file_result[0][0]
        file_path = file_result[0][2]
        file_exists_flag = file_result[0][3]
        
        logger.info(f"Found domain_blast_result record in database:")
        logger.info(f"  ID: {file_id}")
        logger.info(f"  Path: {file_path}")
        logger.info(f"  Exists flag: {file_exists_flag}")
        
        full_path = os.path.join(batch_path, file_path)
        logger.info(f"  Full path: {full_path}")
        logger.info(f"  File exists: {os.path.exists(full_path)}")
        
        # Check if file exists at expected location
        if not os.path.exists(full_path):
            logger.warning(f"File doesn't exist at expected location, searching for alternatives")
            
            # Search for domain blast files
            possible_batch_dirs = ["batch_0", "batch_1", "batch_2"]
            found_path = None
            
            for batch_dir in possible_batch_dirs:
                possible_patterns = [
                    f"{pdb_chain}.domain_blast.xml",
                    f"{pdb_chain}.blast",
                    f"{pdb_chain}_blast.txt"
                ]
                
                for pattern in possible_patterns:
                    test_path = f"blast/domain/{batch_dir}/{pattern}"
                    full_test_path = os.path.join(batch_path, test_path)
                    
                    if os.path.exists(full_test_path):
                        found_path = test_path
                        logger.info(f"Found alternative path: {full_test_path}")
                        break
                
                if found_path:
                    break
            
            if found_path and args.apply_fix:
                # Update the record
                update_query = """
                UPDATE ecod_schema.process_file
                SET file_path = %s, file_exists = TRUE, last_checked = NOW()
                WHERE id = %s
                """
                
                context.db.execute_query(update_query, (found_path, file_id))
                logger.info(f"Updated domain_blast_result path to: {found_path}")
            elif found_path:
                logger.info(f"Would update domain_blast_result path to: {found_path}")
            else:
                logger.error(f"Couldn't find domain_blast_result file for {pdb_chain}")
    else:
        logger.warning(f"No domain_blast_result record found in database")
        
        # Search for domain blast files
        possible_batch_dirs = ["batch_0", "batch_1", "batch_2"]
        found_path = None
        
        for batch_dir in possible_batch_dirs:
            possible_patterns = [
                f"{pdb_chain}.domain_blast.xml",
                f"{pdb_chain}.blast",
                f"{pdb_chain}_blast.txt"
            ]
            
            for pattern in possible_patterns:
                test_path = f"blast/domain/{batch_dir}/{pattern}"
                full_test_path = os.path.join(batch_path, test_path)
                
                if os.path.exists(full_test_path):
                    found_path = test_path
                    logger.info(f"Found domain blast file: {full_test_path}")
                    break
            
            if found_path:
                break
        
        if found_path and args.apply_fix:
            # Create a new record
            insert_query = """
            INSERT INTO ecod_schema.process_file
            (process_id, file_type, file_path, file_exists, file_size, last_checked)
            VALUES (%s, %s, %s, TRUE, %s, NOW())
            """
            
            full_path = os.path.join(batch_path, found_path)
            file_size = os.path.getsize(full_path)
            
            context.db.execute_query(insert_query, (process_id, 'domain_blast_result', found_path, file_size))
            logger.info(f"Created domain_blast_result record with path: {found_path}")
        elif found_path:
            logger.info(f"Would create domain_blast_result record with path: {found_path}")
        else:
            logger.error(f"Couldn't find domain_blast_result file for {pdb_chain}")
    
    # Verify what happens when calling simplified_file_path_resolution
    logger.info("Testing file lookup through simplified_file_path_resolution...")
    
    # Need to import the DomainSummary class to test its method
    try:
        from ecod.pipelines.domain_analysis.summary import DomainSummary
        
        # Create an instance to test
        domain_summary = DomainSummary(args.config)
        
        # Check what happens when trying to resolve the domain_blast_result
        result_path = domain_summary.simplified_file_path_resolution(
            pdb_id, chain_id, 'domain_blast_result', batch_path
        )
        
        if result_path:
            logger.info(f"simplified_file_path_resolution found: {result_path}")
            logger.info(f"File exists: {os.path.exists(result_path)}")
        else:
            logger.error(f"simplified_file_path_resolution failed to find domain_blast_result")
            
            # Let's check the code in the function and what might be going wrong
            logger.info("Possible issues:")
            logger.info("1. File type mismatch - check if using 'domain_blast_result' vs 'domain_blast_results'")
            logger.info("2. Path joining issues - check if base path is correct")
            logger.info("3. Database lookup issues - check SQL query and parameters")
            
            # Debug by checking file type options
            types_to_check = ['domain_blast_result', 'domain_blast_results', 'blast_result']
            
            for file_type in types_to_check:
                logger.info(f"Checking with file_type = '{file_type}'...")
                
                # Query database directly
                check_query = """
                SELECT file_path
                FROM ecod_schema.process_file pf
                JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
                JOIN ecod_schema.protein p ON ps.protein_id = p.id
                WHERE p.pdb_id = %s AND p.chain_id = %s
                AND pf.file_type = %s
                AND pf.file_exists = TRUE
                ORDER BY pf.id DESC
                LIMIT 1
                """
                
                check_result = context.db.execute_query(check_query, (pdb_id, chain_id, file_type))
                
                if check_result:
                    logger.info(f"  Found record with file_type = '{file_type}': {check_result[0][0]}")
                else:
                    logger.info(f"  No record found with file_type = '{file_type}'")
    except ImportError as e:
        logger.error(f"Could not import DomainSummary: {e}")
    except Exception as e:
        logger.error(f"Error testing simplified_file_path_resolution: {e}", exc_info=True)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())