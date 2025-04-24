#!/usr/bin/env python3
"""
create_alt_rep_batches.py - Create and manage batches for novel alternative representatives
"""

import os
import sys
import logging
import argparse
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager

def setup_logging(verbose: bool = False, log_file: str = None):
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

def get_novel_alt_reps(db, batch_size: int = 4000, limit: int = None) -> List[List[Dict[str, Any]]]:
    """Get novel alternative representatives from the classification table and group them into batches"""
    logger = logging.getLogger("ecod.alt_rep_batches")
    
    limit_clause = f"LIMIT {limit}" if limit else ""
    
    query = f"""
    SELECT 
        rc.alt_rep_id,
        a.protein_id as alt_protein_id,
        a.pdb_id,
        a.chain_id,
        a.sequence_md5,
        a.sequence_length
    FROM 
        pdb_analysis.rep_classification rc
    JOIN 
        pdb_analysis.alt_representative_proteins a ON rc.alt_rep_id = a.id
    WHERE 
        rc.classification = 'novel'
    ORDER BY 
        a.id
    {limit_clause}
    """
    
    result = db.execute_dict_query(query)
    
    logger.info(f"Found {len(result)} novel alternative representatives")
    
    # Group into batches
    batches = []
    current_batch = []
    
    for protein in result:
        current_batch.append(protein)
        
        if len(current_batch) >= batch_size:
            batches.append(current_batch)
            current_batch = []
    
    # Add the last batch if it has any proteins
    if current_batch:
        batches.append(current_batch)
    
    logger.info(f"Grouped into {len(batches)} batches")
    return batches

def get_sequence_by_md5(db, md5_hash: str) -> Optional[str]:
    """Get protein sequence by MD5 hash"""
    query = """
    SELECT sequence 
    FROM pdb_analysis.alt_representative_proteins a
    WHERE sequence_md5 = %s
    LIMIT 1
    """
    
    result = db.execute_query(query, (md5_hash,))
    if result and result[0]:
        return result[0][0]
    
    # If not found in alt_representatives, try in protein_sequence
    query = """
    SELECT sequence 
    FROM ecod_schema.protein_sequence
    WHERE md5_hash = %s
    LIMIT 1
    """
    
    result = db.execute_query(query, (md5_hash,))
    if result and result[0]:
        return result[0][0]
    
    return None

def create_batch(db, proteins: List[Dict[str, Any]], batch_num: int, 
                batch_type: str = "alt_rep", ref_version: str = "develop291"):
    """Create a new batch for processing"""
    logger = logging.getLogger("ecod.alt_rep_batches")
    
    # Generate batch name with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    batch_name = f"alt_rep_batch_{batch_num:03d}_{timestamp}"
    
    # Set base path following existing pattern
    base_path = f"/data/ecod/pdb_updates/batches/{batch_name}"
    
    # Insert batch record
    batch_id = db.insert(
        "ecod_schema.batch",
        {
            "batch_name": batch_name,
            "base_path": base_path,
            "type": batch_type,
            "ref_version": ref_version,
            "total_items": len(proteins),
            "status": "created"
        },
        "id"
    )
    
    # Create directory structure
    os.makedirs(base_path, exist_ok=True)
    os.makedirs(os.path.join(base_path, "query_fastas"), exist_ok=True)
    os.makedirs(os.path.join(base_path, "results"), exist_ok=True)
    
    logger.info(f"Created batch {batch_name} with ID {batch_id} containing {len(proteins)} proteins")
    return batch_id, base_path

def ensure_protein_exists(db, protein_data: Dict[str, Any]) -> int:
    """Ensure protein exists in the main ECOD schema"""
    # Check if protein exists by source_id
    source_id = f"{protein_data['pdb_id']}_{protein_data['chain_id']}"
    
    query = """
    SELECT id FROM ecod_schema.protein
    WHERE pdb_id = %s AND chain_id = %s
    """
    
    result = db.execute_query(query, (protein_data['pdb_id'], protein_data['chain_id']))
    
    if result:
        return result[0][0]
    
    # Get sequence from MD5 hash
    sequence = get_sequence_by_md5(db, protein_data['sequence_md5'])
    
    # Insert new protein
    protein_id = db.insert(
        "ecod_schema.protein",
        {
            "pdb_id": protein_data['pdb_id'],
            "chain_id": protein_data['chain_id'],
            "source_id": source_id,
            "length": protein_data['sequence_length']
        },
        "id"
    )
    
    # Insert sequence
    if sequence:
        db.insert(
            "ecod_schema.protein_sequence",
            {
                "protein_id": protein_id,
                "sequence": sequence,
                "md5_hash": protein_data['sequence_md5']
            }
        )
    
    return protein_id

def register_proteins_in_batch(db, batch_id: int, base_path: str, proteins: List[Dict[str, Any]]):
    """Register proteins in a batch and create initial files"""
    logger = logging.getLogger("ecod.alt_rep_batches")
    fasta_dir = os.path.join(base_path, "query_fastas")
    
    for protein in proteins:
        # Determine relative path for this protein
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        rel_path = f"{pdb_id}_{chain_id}"
        
        # First create or get protein entry in ECOD schema
        protein_id = ensure_protein_exists(db, protein)
        
        # Get sequence for FASTA generation
        sequence = get_sequence_by_md5(db, protein['sequence_md5'])
        if not sequence:
            logger.warning(f"Could not find sequence for {pdb_id}_{chain_id} with MD5 {protein['sequence_md5']}")
            continue
        
        # Register in process_status
        process_id = db.insert(
            "ecod_schema.process_status",
            {
                "protein_id": protein_id,
                "batch_id": batch_id,
                "current_stage": "fasta_generated",
                "status": "pending",
                "relative_path": rel_path,
                "is_representative": True  # These are all representatives
            },
            "id"
        )
        
        # Generate FASTA file
        source_id = f"{pdb_id}_{chain_id}"
        fasta_path = os.path.join(fasta_dir, f"{source_id}.fa")
        with open(fasta_path, 'w') as f:
            f.write(f">{source_id}\n{sequence}\n")
        
        # Register FASTA file
        db.insert(
            "ecod_schema.process_file",
            {
                "process_id": process_id,
                "file_type": "fasta",
                "file_path": f"query_fastas/{source_id}.fa",
                "file_exists": True,
                "file_size": os.path.getsize(fasta_path)
            }
        )
        
        logger.info(f"Registered protein {source_id} in batch {batch_id}")

def main():
    parser = argparse.ArgumentParser(description='Create batches for novel alternative representatives')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-size', type=int, default=4000,
                      help='Number of proteins per batch')
    parser.add_argument('--batch-type', type=str, default='alt_rep',
                      help='Type label for the batches')
    parser.add_argument('--ref-version', type=str, default='develop291',
                      help='Reference version')
    parser.add_argument('--dry-run', action='store_true',
                      help='Show what would be done without making changes')
    parser.add_argument('--start-from', type=int, default=1,
                      help='Start batch numbering from this value')
    parser.add_argument('--limit', type=int,
                      help='Limit number of proteins to process (for testing)')
    parser.add_argument('--log-file', type=str,
                      help='Path to log file')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.alt_rep_batches")
    
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Get novel alternative representatives
    batches = get_novel_alt_reps(db, args.batch_size, args.limit)
    
    if args.dry_run:
        logger.info(f"Dry run: would create {len(batches)} batches")
        for i, batch in enumerate(batches):
            batch_num = i + args.start_from
            logger.info(f"  Batch {batch_num}: {len(batch)} proteins")
        return 0
    
    # Process batches
    for i, proteins in enumerate(batches):
        batch_num = i + args.start_from
        batch_id, base_path = create_batch(
            db, 
            proteins, 
            batch_num, 
            args.batch_type, 
            args.ref_version
        )
        
        register_proteins_in_batch(db, batch_id, base_path, proteins)
        
        logger.info(f"Completed setup for batch {batch_num}/{len(batches) + args.start_from - 1}")
    
    logger.info(f"Successfully created {len(batches)} batches for novel alternative representatives")
    return 0

if __name__ == "__main__":
    sys.exit(main())