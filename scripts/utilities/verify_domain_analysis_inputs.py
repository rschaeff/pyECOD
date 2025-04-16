#!/usr/bin/env python3
"""
verify_domain_analysis_inputs.py - Verify input files for domain analysis

This script checks that BLAST result files exist and are correctly indexed in the database
for a given batch or specific protein.
"""

import os
import sys
import logging
import argparse
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

def verify_protein_files(context: Any, protein_id: int, batch_id: int, base_path: str) -> Dict[str, Any]:
    """
    Verify files for a specific protein
    
    Args:
        context: Application context
        protein_id: Protein ID
        batch_id: Batch ID
        base_path: Base path for ECOD data
        
    Returns:
        Dictionary with verification results
    """
    logger = logging.getLogger("ecod.verify_inputs")
    
    # Get protein and process information
    query = """
    SELECT 
        p.id, p.source_id, p.pdb_id, p.chain_id,
        ps.id as process_id,
        b.batch_name, b.base_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    WHERE 
        p.id = %s AND ps.batch_id = %s
    """
    
    results = context.db.execute_query(query, (protein_id, batch_id))
    
    if not results:
        logger.error(f"Protein {protein_id} not found in batch {batch_id}")
        return {"error": f"Protein {protein_id} not found in batch {batch_id}"}
    
    protein_info = results[0]
    process_id = protein_info[4]  # process_id
    source_id = protein_info[1]  # source_id
    pdb_id = protein_info[2]  # pdb_id
    chain_id = protein_info[3]  # chain_id
    batch_name = protein_info[5]  # batch_name
    batch_path = protein_info[6] if base_path is None else os.path.join(base_path, 'batches', batch_name)
    
    logger.info(f"Verifying files for protein: {source_id} (ID: {protein_id})")
    logger.info(f"Using batch path: {batch_path}")
    
    # Get process files from database
    file_query = """
    SELECT 
        id, file_type, file_path, file_exists
    FROM 
        ecod_schema.process_file
    WHERE 
        process_id = %s
    """
    
    file_results = context.db.execute_query(file_query, (process_id,))
    
    # Organize files by type
    db_files = {}
    for file_result in file_results:
        file_id = file_result[0]
        file_type = file_result[1]
        file_path = file_result[2]
        file_exists_flag = file_result[3]
        
        db_files[file_type] = {
            "id": file_id,
            "path": file_path,
            "exists_flag": file_exists_flag,
            "abs_path": os.path.join(batch_path, file_path) if not os.path.isabs(file_path) else file_path
        }
    
    # Check for required file types
    required_types = ['chain_blast_result', 'domain_blast_result', 'fasta']
    verification = {
        "protein_id": protein_id,
        "source_id": source_id,
        "process_id": process_id,
        "files": {}
    }
    
    # Check each file type
    for file_type in required_types:
        if file_type in db_files:
            file_info = db_files[file_type]
            abs_path = file_info["abs_path"]
            
            # Check if file exists on filesystem
            exists_on_disk = os.path.exists(abs_path)
            
            verification["files"][file_type] = {
                "in_db": True,
                "db_path": file_info["path"],
                "exists_on_disk": exists_on_disk,
                "exists_flag_correct": file_info["exists_flag"] == exists_on_disk,
                "abs_path": abs_path
            }
            
            # Log results
            status = "✓" if exists_on_disk else "✗"
            logger.info(f"{status} {file_type}: {abs_path}")
            
            if file_info["exists_flag"] != exists_on_disk:
                flag_status = "✗ DB flag incorrect"
                logger.warning(f"{flag_status} for {file_type}: DB says {file_info['exists_flag']}, actual is {exists_on_disk}")
        else:
            verification["files"][file_type] = {
                "in_db": False,
                "exists_on_disk": False,
                "exists_flag_correct": False
            }
            
            logger.warning(f"✗ {file_type}: Not indexed in database")
    
    # Look for files on filesystem (focusing on possible sub-batch directories)
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    # Function to search recursively for files matching the protein
    def find_files_for_protein(directory, filename_pattern):
        found_files = []
        
        if not os.path.exists(directory):
            return found_files
            
        for root, dirs, files in os.walk(directory):
            for file in files:
                if filename_pattern in file:
                    found_files.append(os.path.join(root, file))
                    
        return found_files
    
    # Search for FASTA files
    fasta_dir = os.path.join(batch_path, 'fastas')
    fasta_files = find_files_for_protein(fasta_dir, pdb_chain)
    
    # Search for chain BLAST files
    chain_blast_dir = os.path.join(batch_path, 'blast', 'chain')
    chain_blast_files = find_files_for_protein(chain_blast_dir, pdb_chain)
    
    # Search for domain BLAST files
    domain_blast_dir = os.path.join(batch_path, 'blast', 'domain')
    domain_blast_files = find_files_for_protein(domain_blast_dir, pdb_chain)
    
    # Add filesystem search results
    verification["filesystem"] = {
        "fasta_files": fasta_files,
        "chain_blast_files": chain_blast_files,
        "domain_blast_files": domain_blast_files
    }
    
    # Log filesystem search results
    if fasta_files and 'fasta' not in db_files:
        logger.info(f"Found FASTA files on filesystem not in DB: {fasta_files}")
    
    if chain_blast_files and 'chain_blast_result' not in db_files:
        logger.info(f"Found chain BLAST files on filesystem not in DB: {chain_blast_files}")
    
    if domain_blast_files and 'domain_blast_result' not in db_files:
        logger.info(f"Found domain BLAST files on filesystem not in DB: {domain_blast_files}")
    
    # Summarize results
    missing_files = [t for t in required_types if t not in db_files or not verification["files"][t]["exists_on_disk"]]
    
    if missing_files:
        verification["status"] = "incomplete"
        verification["missing_files"] = missing_files
        logger.warning(f"Missing required files: {', '.join(missing_files)}")
    else:
        verification["status"] = "ready"
        logger.info("All required files found and correctly indexed")
    
    return verification

def suggest_fixes(context: Any, verification: Dict[str, Any], apply_fixes: bool = False) -> List[str]:
    """
    Suggest fixes for issues found during verification
    
    Args:
        context: Application context
        verification: Verification results
        apply_fixes: Whether to apply suggested fixes
        
    Returns:
        List of suggested fixes
    """
    logger = logging.getLogger("ecod.verify_inputs")
    
    if verification.get("error"):
        return [f"Error: {verification['error']}"]
    
    fixes = []
    
    # Check for missing files in DB but found on filesystem
    if "filesystem" in verification:
        fs = verification["filesystem"]
        
        # Check FASTA files
        if fs["fasta_files"] and ("fasta" not in verification["files"] or not verification["files"]["fasta"]["in_db"]):
            fix = f"Add FASTA file record: {fs['fasta_files'][0]}"
            fixes.append(fix)
            
            if apply_fixes:
                try:
                    # Insert FASTA file record
                    file_path = os.path.relpath(fs["fasta_files"][0], verification.get("batch_path", ""))
                    exists = True
                    file_size = os.path.getsize(fs["fasta_files"][0])
                    
                    query = """
                    INSERT INTO ecod_schema.process_file
                    (process_id, file_type, file_path, file_exists, file_size, last_checked)
                    VALUES (%s, %s, %s, %s, %s, NOW())
                    """
                    
                    context.db.execute_query(query, (verification["process_id"], "fasta", file_path, exists, file_size))
                    logger.info(f"Applied fix: Added FASTA file record {file_path}")
                except Exception as e:
                    logger.error(f"Error applying fix: {e}")
        
        # Check chain BLAST files
        if fs["chain_blast_files"] and ("chain_blast_result" not in verification["files"] or not verification["files"]["chain_blast_result"]["in_db"]):
            fix = f"Add chain BLAST file record: {fs['chain_blast_files'][0]}"
            fixes.append(fix)
            
            if apply_fixes:
                try:
                    # Insert chain BLAST file record
                    file_path = os.path.relpath(fs["chain_blast_files"][0], verification.get("batch_path", ""))
                    exists = True
                    file_size = os.path.getsize(fs["chain_blast_files"][0])
                    
                    query = """
                    INSERT INTO ecod_schema.process_file
                    (process_id, file_type, file_path, file_exists, file_size, last_checked)
                    VALUES (%s, %s, %s, %s, %s, NOW())
                    """
                    
                    context.db.execute_query(query, (verification["process_id"], "chain_blast_result", file_path, exists, file_size))
                    logger.info(f"Applied fix: Added chain BLAST file record {file_path}")
                except Exception as e:
                    logger.error(f"Error applying fix: {e}")
        
        # Check domain BLAST files
        if fs["domain_blast_files"] and ("domain_blast_result" not in verification["files"] or not verification["files"]["domain_blast_result"]["in_db"]):
            fix = f"Add domain BLAST file record: {fs['domain_blast_files'][0]}"
            fixes.append(fix)
            
            if apply_fixes:
                try:
                    # Insert domain BLAST file record
                    file_path = os.path.relpath(fs["domain_blast_files"][0], verification.get("batch_path", ""))
                    exists = True
                    file_size = os.path.getsize(fs["domain_blast_files"][0])
                    
                    query = """
                    INSERT INTO ecod_schema.process_file
                    (process_id, file_type, file_path, file_exists, file_size, last_checked)
                    VALUES (%s, %s, %s, %s, %s, NOW())
                    """
                    
                    context.db.execute_query(query, (verification["process_id"], "domain_blast_result", file_path, exists, file_size))
                    logger.info(f"Applied fix: Added domain BLAST file record {file_path}")
                except Exception as e:
                    logger.error(f"Error applying fix: {e}")
    
    # Check for incorrect exists_flag
    for file_type, file_info in verification.get("files", {}).items():
        if file_info["in_db"] and not file_info["exists_flag_correct"]:
            exists_on_disk = file_info["exists_on_disk"]
            file_id = file_info.get("id")
            
            if file_id:
                fix = f"Update '{file_type}' exists_flag to {exists_on_disk}"
                fixes.append(fix)
                
                if apply_fixes:
                    try:
                        query = """
                        UPDATE ecod_schema.process_file
                        SET file_exists = %s, last_checked = NOW()
                        WHERE id = %s
                        """
                        
                        context.db.execute_query(query, (exists_on_disk, file_id))
                        logger.info(f"Applied fix: Updated exists_flag for {file_type}")
                    except Exception as e:
                        logger.error(f"Error applying fix: {e}")
    
    return fixes

def main():
    """Main function to verify domain analysis inputs"""
    parser = argparse.ArgumentParser(description='Verify input files for domain analysis')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to verify')
    parser.add_argument('--protein-id', type=int,
                      help='Specific protein ID to verify (omit to verify all in batch)')
    parser.add_argument('--limit', type=int, default=10,
                      help='Maximum number of proteins to verify (default: 10)')
    parser.add_argument('--base-path', type=str,
                      help='Base path for ECOD data structure (default: from batch record)')
    parser.add_argument('--apply-fixes', action='store_true',
                      help='Apply suggested fixes')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.verify_inputs")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Verify batch exists
    batch_query = "SELECT batch_name, base_path FROM ecod_schema.batch WHERE id = %s"
    batch_result = context.db.execute_query(batch_query, (args.batch_id,))
    
    if not batch_result:
        logger.error(f"Batch ID {args.batch_id} not found")
        return 1
    
    batch_name = batch_result[0][0]
    batch_path = batch_result[0][1]
    
    logger.info(f"Verifying inputs for batch: {batch_name} (ID: {args.batch_id})")
    
    # Get proteins to verify
    if args.protein_id:
        # Verify specific protein
        verification = verify_protein_files(context, args.protein_id, args.batch_id, args.base_path)
        
        # Suggest fixes
        fixes = suggest_fixes(context, verification, args.apply_fixes)
        
        if fixes:
            logger.info(f"Suggested fixes: {len(fixes)}")
            for fix in fixes:
                logger.info(f"  - {fix}")
        
        # Print summary
        if verification.get("status") == "ready":
            logger.info(f"Protein {args.protein_id} is ready for domain analysis")
        else:
            logger.warning(f"Protein {args.protein_id} is not ready for domain analysis: {verification.get('status')}")
    else:
        # Get proteins from the batch
        protein_query = """
        SELECT 
            p.id
        FROM 
            ecod_schema.protein p
        JOIN 
            ecod_schema.process_status ps ON p.id = ps.protein_id
        WHERE 
            ps.batch_id = %s
        ORDER BY 
            p.id
        LIMIT %s
        """
        
        protein_results = context.db.execute_query(protein_query, (args.batch_id, args.limit))
        
        if not protein_results:
            logger.warning(f"No proteins found in batch {args.batch_id}")
            return 0
        
        # Track results
        ready_count = 0
        incomplete_count = 0
        error_count = 0
        fixes_applied = 0
        
        # Verify each protein
        for protein_result in protein_results:
            protein_id = protein_result[0]
            
            verification = verify_protein_files(context, protein_id, args.batch_id, args.base_path)
            
            # Suggest and apply fixes
            fixes = suggest_fixes(context, verification, args.apply_fixes)
            
            if fixes:
                if args.apply_fixes:
                    fixes_applied += len(fixes)
                logger.info(f"Protein {protein_id}: {len(fixes)} fixes suggested")
            
            # Track status
            if verification.get("error"):
                error_count += 1
            elif verification.get("status") == "ready":
                ready_count += 1
            else:
                incomplete_count += 1
        
        # Print summary
        logger.info(f"Verification complete: {ready_count} ready, {incomplete_count} incomplete, {error_count} errors")
        
        if args.apply_fixes:
            logger.info(f"Applied {fixes_applied} fixes")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())