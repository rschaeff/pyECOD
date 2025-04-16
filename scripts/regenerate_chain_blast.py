#!/usr/bin/env python3
"""
regenerate_chain_blast.py - Regenerate missing chain BLAST files for proteins with domain BLAST

This script leverages the existing BlastPipeline to regenerate chain BLAST files
for proteins that have domain BLAST results but are missing chain BLAST files.
"""

import os
import sys
import logging
import argparse
import subprocess
from typing import Dict, List, Any, Optional
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.pipelines.blast_pipeline import BlastPipeline
from ecod.jobs import JobManager
from ecod.exceptions import PipelineError, JobSubmissionError

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
    
    return logging.getLogger("ecod.regenerate_chain_blast")

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
    logger = logging.getLogger("ecod.regenerate_chain_blast")
    
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
        p.id, p.pdb_id, p.chain_id, p.length, p.source_id,
        ps.id as process_id, ps.current_stage, ps.status, ps.error_message
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
    ORDER BY 
        ps.updated_at DESC
    LIMIT %s
    """
    
    candidates = context.db.execute_dict_query(query, (batch_id, limit))
    
    if not candidates:
        logger.info(f"No proteins found with missing chain BLAST in batch {batch_id}")
        return []
    
    logger.info(f"Found {len(candidates)} proteins with missing chain BLAST")
    
    # Get protein sequences for candidates
    proteins = []
    for candidate in candidates:
        protein_id = candidate['id']
        # Get sequence from protein_sequence table
        seq_query = """
        SELECT sequence FROM ecod_schema.protein_sequence WHERE protein_id = %s
        """
        seq_result = context.db.execute_query(seq_query, (protein_id,))
        
        if seq_result and seq_result[0][0]:
            candidate['sequence'] = seq_result[0][0]
            proteins.append(candidate)
        else:
            logger.warning(f"No sequence found for protein {protein_id}")
    
    return proteins

def generate_chain_blast_files(context, pipeline, proteins, batch_id):
    """Generate chain BLAST files for the specified proteins"""
    logger = logging.getLogger("ecod.regenerate_chain_blast")
    
    # Get batch path
    batch_query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return []
    
    batch_path = batch_result[0][0]
    
    # Create query_fastas directory if it doesn't exist
    fasta_dir = os.path.join(batch_path, "query_fastas")
    os.makedirs(fasta_dir, exist_ok=True)
    
    # Create chain blast directory if it doesn't exist
    chain_blast_dir = os.path.join(batch_path, "chain_blast_results")
    os.makedirs(chain_blast_dir, exist_ok=True)
    
    # Generate FASTA files
    fasta_paths = []
    for protein in proteins:
        protein_id = protein['id']
        process_id = protein['process_id']
        source_id = protein['source_id']
        sequence = protein['sequence']
        
        # Create FASTA file
        fasta_path = os.path.join(fasta_dir, f"{source_id}.fa")
        with open(fasta_path, 'w') as f:
            f.write(f">{source_id}\n{sequence}\n")
        
        # Register FASTA file if not already registered
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'fasta'
        """
        
        result = context.db.execute_query(check_query, (process_id,))
        
        if not result:
            # Register new FASTA file
            context.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "fasta",
                    "file_path": f"query_fastas/{source_id}.fa",
                    "file_exists": True,
                    "file_size": os.path.getsize(fasta_path)
                }
            )
            logger.info(f"Registered FASTA file for {source_id}")
        
        fasta_paths.append((process_id, fasta_path))
    
    if not fasta_paths:
        logger.warning("No FASTA files generated")
        return
    
    # Prepare for chain BLAST
    logger.info(f"Submitting chain BLAST for {len(fasta_paths)} proteins")
    
    try:
        # Get chain BLAST database
        chain_db = context.config.get('reference', {}).get('chain_db')
        if not chain_db:
            logger.error("Chain BLAST database not configured")
            return
        
        # Define batch size
        batch_size = context.config.get('pipeline', {}).get('batch_size', 10)
        
        # Use pipeline to run chain BLAST
        job_ids = pipeline.run_chain_blast(batch_id, batch_size)
        
        if job_ids:
            logger.info(f"Submitted {len(job_ids)} chain BLAST jobs")
            
            # Check job status
            logger.info("Checking job status...")
            pipeline.check_job_status(batch_id)
            
            logger.info("Chain BLAST jobs submitted successfully")
        else:
            logger.warning("No chain BLAST jobs were submitted")
    
    except Exception as e:
        logger.error(f"Error running chain BLAST: {e}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Regenerate missing chain BLAST files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, default=31,
                      help='Batch ID to process (default: 31)')
    parser.add_argument('--limit', type=int, default=10,
                      help='Maximum number of proteins to process (default: 10)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('--verbose', '-v', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--run-blast', action='store_true',
                      help='Run BLAST immediately (not just prepare)')
    parser.add_argument('--protein-id', type=int,
                      help='Process specific protein ID')
    
    args = parser.parse_args()
    logger = setup_logging(args.verbose, args.log_file)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Initialize job manager
    job_manager = JobManager(args.config)
    
    # Initialize BlastPipeline
    pipeline = BlastPipeline(context.db_manager, job_manager, context.config)
    
    # Process specific protein if requested
    if args.protein_id:
        query = """
        SELECT 
            p.id, p.pdb_id, p.chain_id, p.length, p.source_id,
            ps.id as process_id, ps.current_stage, ps.status,
            ps_seq.sequence
        FROM 
            ecod_schema.protein p
        JOIN 
            ecod_schema.process_status ps ON p.id = ps.protein_id
        JOIN 
            ecod_schema.protein_sequence ps_seq ON p.id = ps_seq.protein_id
        WHERE 
            p.id = %s
        """
        
        result = context.db.execute_dict_query(query, (args.protein_id,))
        
        if not result:
            logger.error(f"Protein {args.protein_id} not found")
            return 1
        
        proteins = result
    else:
        # Find proteins with missing chain BLAST
        proteins = find_proteins_with_missing_chain_blast(context, args.batch_id, args.limit)
    
    if not proteins:
        logger.error("No proteins found to process")
        return 1
    
    logger.info(f"Found {len(proteins)} proteins to process")
    
    # Print protein information
    print("\nProteins with missing chain BLAST files:")
    print(f"{'Protein ID':<10} {'PDB Chain':<12} {'Length':<8} {'Stage':<20} {'Status':<10}")
    print("-" * 70)
    
    for protein in proteins:
        print(f"{protein['id']:<10} {protein['pdb_id']}_{protein['chain_id']:<7} {protein['length']:<8} {protein['current_stage']:<20} {protein['status']:<10}")
    
    # Generate FASTA files and set up for BLAST
    generate_chain_blast_files(context, pipeline, proteins, args.batch_id)
    
    if args.run_blast:
        logger.info("Checking job status...")
        pipeline.check_job_status(args.batch_id)
        
        # Generate domain summary commands
        print("\nRun these commands to generate domain summaries:")
        for protein in proteins:
            cmd = f"python scripts/generate_domain_summary_v2.py --config config/config.yml --batch-id {args.batch_id} --protein-id {protein['id']} --blast-only -v"
            print(cmd)
    else:
        print("\nBLAST setup complete. Run the following to check job status:")
        print(f"python scripts/check_jobs.py --batch-id {args.batch_id}")
        print("\nOnce jobs complete, run these commands to generate domain summaries:")
        for protein in proteins:
            cmd = f"python scripts/generate_domain_summary_v2.py --config config/config.yml --batch-id {args.batch_id} --protein-id {protein['id']} --blast-only -v"
            print(cmd)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())