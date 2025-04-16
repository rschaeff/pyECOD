#!/usr/bin/env python3
"""
regenerate_chain_blast.py - Regenerate missing chain BLAST files for proteins

This script identifies proteins with missing chain BLAST files but valid domain BLAST files,
regenerates the chain BLAST files, and updates the database.
"""

import os
import sys
import time
import logging
import argparse
import subprocess
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Optional
from pathlib import Path

# Add parent directory to path for imports
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
    
    return logging.getLogger("ecod.regenerate_blast")

def find_proteins_with_missing_chain_blast(context, batch_id: int, limit: int = 10) -> List[Dict[str, Any]]:
    """
    Find proteins with missing chain BLAST files but valid domain BLAST files
    
    Args:
        context: Application context
        batch_id: Batch ID to search in
        limit: Maximum number of proteins to return
        
    Returns:
        List of protein dictionaries
    """
    logger = logging.getLogger("ecod.regenerate_blast")
    
    # First get batch information
    batch_query = """
    SELECT 
        id, batch_name, base_path, ref_version
    FROM 
        ecod_schema.batch
    WHERE 
        id = %s
    """
    
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return []
    
    batch_info = {
        'id': batch_result[0][0],
        'name': batch_result[0][1],
        'base_path': batch_result[0][2],
        'reference': batch_result[0][3]
    }
    
    logger.info(f"Analyzing batch {batch_id} ({batch_info['name']})")
    
    # Find proteins with domain BLAST but no chain BLAST
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, p.length,
        ps.id as process_id, ps.current_stage, ps.status
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.process_file pf_domain ON ps.id = pf_domain.process_id AND pf_domain.file_type = 'domain_blast_result'
    LEFT JOIN 
        ecod_schema.process_file pf_chain ON ps.id = pf_chain.process_id AND pf_chain.file_type = 'blast_result'
    WHERE 
        ps.batch_id = %s
        AND pf_domain.file_exists = TRUE
        AND (pf_chain.id IS NULL OR pf_chain.file_exists = FALSE)
        AND (ps.status = 'error' OR ps.current_stage LIKE 'domain_%_failed')
    ORDER BY 
        ps.updated_at DESC
    LIMIT %s
    """
    
    candidates = context.db.execute_query(query, (batch_id, limit))
    
    if not candidates:
        logger.info(f"No proteins found with missing chain BLAST in batch {batch_id}")
        return []
    
    logger.info(f"Found {len(candidates)} proteins with missing chain BLAST")
    
    # Process candidates
    results = []
    
    for candidate in candidates:
        protein_id = candidate[0]
        pdb_id = candidate[1]
        chain_id = candidate[2]
        length = candidate[3]
        process_id = candidate[4]
        current_stage = candidate[5]
        status = candidate[6]
        
        # Get FASTA file path
        fasta_query = """
        SELECT file_path 
        FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'fasta' AND file_exists = TRUE
        """
        
        fasta_result = context.db.execute_query(fasta_query, (process_id,))
        fasta_path = None
        
        if fasta_result:
            rel_path = fasta_result[0][0]
            fasta_path = os.path.join(batch_info['base_path'], rel_path)
            fasta_path = os.path.normpath(fasta_path)
            
            if not os.path.exists(fasta_path):
                logger.warning(f"FASTA file not found at {fasta_path}")
                fasta_path = None
        
        # Get domain BLAST file path
        domain_query = """
        SELECT file_path 
        FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'domain_blast_result' AND file_exists = TRUE
        """
        
        domain_result = context.db.execute_query(domain_query, (process_id,))
        domain_blast_path = None
        
        if domain_result:
            rel_path = domain_result[0][0]
            domain_blast_path = os.path.join(batch_info['base_path'], rel_path)
            domain_blast_path = os.path.normpath(domain_blast_path)
            
            if not os.path.exists(domain_blast_path):
                logger.warning(f"Domain BLAST file not found at {domain_blast_path}")
                domain_blast_path = None
        
        # Check if we have what we need
        if not fasta_path:
            logger.warning(f"Skipping {pdb_id}_{chain_id} - FASTA file not found")
            continue
        
        # Add to results
        results.append({
            'protein_id': protein_id,
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'length': length,
            'process_id': process_id,
            'current_stage': current_stage,
            'status': status,
            'fasta_path': fasta_path,
            'domain_blast_path': domain_blast_path,
            'batch_path': batch_info['base_path'],
            'reference': batch_info['reference']
        })
    
    return results

def run_chain_blast(protein: Dict[str, Any], context, config: Dict[str, Any]) -> Optional[str]:
    """
    Run chain BLAST for a protein
    
    Args:
        protein: Protein information dictionary
        context: Application context
        config: Configuration dictionary
        
    Returns:
        Path to generated chain BLAST file, or None if failed
    """
    logger = logging.getLogger("ecod.regenerate_blast")
    logger.info(f"Running chain BLAST for {protein['pdb_id']}_{protein['chain_id']}")
    
    # Get BLAST configuration
    blast_bin = config.get('tools', {}).get('blast_path', '/usr/bin')
    if not os.path.exists(blast_bin):
        logger.error(f"BLAST binary directory not found: {blast_bin}")
        return None
    
    blast_db = config.get('ecod', {}).get('pfam_blast_db', '/usr/data/blast/pfam')
    if not os.path.exists(f"{blast_db}.phr"):
        logger.error(f"BLAST database not found: {blast_db}")
        return None
    
    # Set up output directory
    pdb_chain = f"{protein['pdb_id']}_{protein['chain_id']}"
    chain_dir = os.path.join(protein['batch_path'], pdb_chain)
    os.makedirs(chain_dir, exist_ok=True)
    
    # Set up output file path
    output_file = os.path.join(chain_dir, f"{pdb_chain}.chainwise_blast.xml")
    
    # Check if file already exists
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        logger.info(f"Chain BLAST file already exists at {output_file}")
        return output_file
    
    # Run BLAST
    blast_cmd = [
        os.path.join(blast_bin, "blastp"),
        "-query", protein['fasta_path'],
        "-db", blast_db,
        "-outfmt", "5",  # XML format
        "-out", output_file,
        "-evalue", "1e-3",
        "-num_threads", "1",
        "-max_target_seqs", "10"
    ]
    
    try:
        logger.debug(f"Running command: {' '.join(blast_cmd)}")
        
        # Run blastp and capture output
        result = subprocess.run(
            blast_cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logger.info(f"Chain BLAST completed successfully: {output_file}")
            return output_file
        else:
            logger.error(f"Chain BLAST output file not found or empty")
            logger.error(f"BLAST stdout: {result.stdout}")
            logger.error(f"BLAST stderr: {result.stderr}")
            return None
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running BLAST: {e}")
        logger.error(f"BLAST stdout: {e.stdout}")
        logger.error(f"BLAST stderr: {e.stderr}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return None

def update_database(protein: Dict[str, Any], blast_file_path: str, context) -> bool:
    """
    Update database with new chain BLAST file
    
    Args:
        protein: Protein information dictionary
        blast_file_path: Path to the chain BLAST file
        context: Application context
        
    Returns:
        True if database was updated successfully, False otherwise
    """
    logger = logging.getLogger("ecod.regenerate_blast")
    
    try:
        # Get relative path
        rel_path = os.path.relpath(blast_file_path, protein['batch_path'])
        
        # Check if chain blast file already exists in database
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'blast_result'
        """
        
        result = context.db.execute_query(check_query, (protein['process_id'],))
        
        if result:
            # Update existing record
            file_id = result[0][0]
            update_query = """
            UPDATE ecod_schema.process_file
            SET file_path = %s, file_exists = TRUE, file_size = %s, last_checked = NOW()
            WHERE id = %s
            """
            
            file_size = os.path.getsize(blast_file_path)
            context.db.execute_query(update_query, (rel_path, file_size, file_id))
            logger.info(f"Updated existing blast_result record in database")
        else:
            # Create new record
            insert_query = """
            INSERT INTO ecod_schema.process_file
            (process_id, file_type, file_path, file_exists, file_size, last_checked)
            VALUES (%s, %s, %s, %s, %s, NOW())
            """
            
            file_size = os.path.getsize(blast_file_path)
            context.db.execute_query(insert_query, (protein['process_id'], 'blast_result', rel_path, True, file_size))
            logger.info(f"Added new blast_result record to database")
        
        # Update process status
        status_query = """
        UPDATE ecod_schema.process_status
        SET current_stage = 'blast_search', status = 'success', error_message = NULL, updated_at = NOW()
        WHERE id = %s
        """
        
        context.db.execute_query(status_query, (protein['process_id'],))
        logger.info(f"Updated process status to blast_search:success")
        
        return True
    except Exception as e:
        logger.error(f"Error updating database: {e}")
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Regenerate missing chain BLAST files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, default=31,
                      help='Batch ID to process (default: 31)')
    parser.add_argument('--limit', type=int, default=10,
                      help='Maximum number of proteins to process (default: 10)')
    parser.add_argument('--output-dir', type=str,
                      help='Override output directory')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--protein-id', type=int,
                      help='Process specific protein ID')
    
    args = parser.parse_args()
    logger = setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    config = context.config
    
    # Find proteins to process
    if args.protein_id:
        # Get specific protein information
        query = """
        SELECT 
            p.id, p.pdb_id, p.chain_id, p.length,
            ps.id as process_id, ps.current_stage, ps.status,
            b.id as batch_id, b.base_path, b.ref_version
        FROM 
            ecod_schema.protein p
        JOIN 
            ecod_schema.process_status ps ON p.id = ps.protein_id
        JOIN 
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            p.id = %s
        """
        
        result = context.db.execute_query(query, (args.protein_id,))
        
        if not result:
            logger.error(f"Protein {args.protein_id} not found")
            return 1
        
        protein_info = {
            'protein_id': result[0][0],
            'pdb_id': result[0][1],
            'chain_id': result[0][2],
            'length': result[0][3],
            'process_id': result[0][4],
            'current_stage': result[0][5],
            'status': result[0][6],
            'batch_id': result[0][7],
            'batch_path': result[0][8],
            'reference': result[0][9]
        }
        
        # Get FASTA file path
        fasta_query = """
        SELECT file_path 
        FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'fasta' AND file_exists = TRUE
        """
        
        fasta_result = context.db.execute_query(fasta_query, (protein_info['process_id'],))
        
        if fasta_result:
            rel_path = fasta_result[0][0]
            fasta_path = os.path.join(protein_info['batch_path'], rel_path)
            fasta_path = os.path.normpath(fasta_path)
            
            if os.path.exists(fasta_path):
                protein_info['fasta_path'] = fasta_path
            else:
                logger.error(f"FASTA file not found at {fasta_path}")
                return 1
        else:
            logger.error(f"No FASTA file found for protein {args.protein_id}")
            return 1
        
        proteins = [protein_info]
    else:
        # Find proteins with missing chain BLAST
        proteins = find_proteins_with_missing_chain_blast(context, args.batch_id, args.limit)
    
    if not proteins:
        logger.error("No proteins found to process")
        return 1
    
    logger.info(f"Found {len(proteins)} proteins to process")
    
    # Process proteins
    success_count = 0
    
    for protein in proteins:
        # Run chain BLAST
        blast_file = run_chain_blast(protein, context, config)
        
        if not blast_file:
            logger.error(f"Failed to generate chain BLAST for {protein['pdb_id']}_{protein['chain_id']}")
            continue
        
        # Update database
        if update_database(protein, blast_file, context):
            logger.info(f"Successfully updated database for {protein['pdb_id']}_{protein['chain_id']}")
            success_count += 1
        else:
            logger.error(f"Failed to update database for {protein['pdb_id']}_{protein['chain_id']}")
    
    logger.info(f"Processed {len(proteins)} proteins, {success_count} successful")
    
    # Generate domain summary commands
    if success_count > 0:
        print("\nRun these commands to generate domain summaries:")
        
        for protein in proteins:
            if os.path.exists(os.path.join(protein['batch_path'], f"{protein['pdb_id']}_{protein['chain_id']}", f"{protein['pdb_id']}_{protein['chain_id']}.chainwise_blast.xml")):
                cmd = f"python scripts/generate_domain_summary_v2.py --config config/config.yml --batch-id {protein.get('batch_id', args.batch_id)} --protein-id {protein['protein_id']} --blast-only -v"
                print(cmd)
    
    return 0 if success_count > 0 else 1

if __name__ == "__main__":
    sys.exit(main())