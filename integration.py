#!/usr/bin/env python3
"""
integration.py - Example of integrating domain summary and partition into the pipeline
"""

import os
import sys
import argparse
import logging
import subprocess
from typing import Dict, Any, List, Optional, Tuple

# Import your domain modules
from domain_summary import DomainSummary
from domain_partition import DomainPartition

# Import pipeline components
# In a real implementation, you would import your existing pipeline modules
from ecod.config import ConfigManager
from ecod.db import DBManager

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

def run_domain_analysis_pipeline(batch_id: int, config_path: str = None, concise: bool = False):
    """Run the domain analysis pipeline (summary and partition)"""
    logger = logging.getLogger("ecod.integration")
    
    # Initialize config
    config_manager = ConfigManager(config_path)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    # Initialize domain analysis components
    domain_summary = DomainSummary(config_path)
    domain_partition = DomainPartition(config_path)
    
    # Get batch information
    query = """
    SELECT 
        b.id, b.batch_name, b.base_path, b.type, b.ref_version
    FROM 
        ecod_schema.batch b
    WHERE 
        b.id = %s
    """
    batch_rows = db.execute_dict_query(query, (batch_id,))
    
    if not batch_rows:
        logger.error(f"Batch {batch_id} not found")
        return False
    
    batch_info = batch_rows[0]
    batch_path = batch_info["base_path"]
    reference = batch_info["ref_version"]
    
    # Get proteins in this batch
    query = """
    SELECT 
        ps.id, p.id AS protein_id, p.pdb_id, p.chain_id, p.source_id,
        ps.current_stage, ps.status, ps.relative_path
    FROM 
        ecod_schema.process_status ps
    JOIN
        ecod_schema.protein p ON ps.protein_id = p.id
    WHERE 
        ps.batch_id = %s
        AND ps.status = 'success'
        AND EXISTS (
            SELECT 1 FROM ecod_schema.process_file pf
            WHERE pf.process_id = ps.id
            AND pf.file_type IN ('chain_blast_result', 'domain_blast_result')
            AND pf.file_exists = TRUE
        )
    """
    protein_rows = db.execute_dict_query(query, (batch_id,))
    
    if not protein_rows:
        logger.error(f"No proteins ready for domain analysis in batch {batch_id}")
        return False
    
    logger.info(f"Found {len(protein_rows)} proteins for domain analysis")
    
    # Process each protein
    for protein in protein_rows:
        pdb_id = protein["pdb_id"]
        chain_id = protein["chain_id"]
        source_id = protein["source_id"]
        process_id = protein["id"]
        
        logger.info(f"Processing {pdb_id}_{chain_id}")
        
        # Step 1: Create domain summary
        try:
            summary_file = domain_summary.create_summary(
                pdb_id,
                chain_id,
                reference,
                batch_path,
                concise
            )
            
            if not summary_file:
                logger.error(f"Failed to create domain summary for {pdb_id}_{chain_id}")
                continue
                
            # Register summary file
            db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "domain_summary",
                    "file_path": summary_file,
                    "file_exists": True
                }
            )
            
            # Update process status
            db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": "domain_summary_complete",
                    "status": "success"
                },
                "id = %s",
                (process_id,)
            )
            
            logger.info(f"Created domain summary for {pdb_id}_{chain_id}")
            
        except Exception as e:
            logger.error(f"Error creating domain summary for {pdb_id}_{chain_id}: {e}")
            continue
        
        # Step 2: Perform domain partition
        try:
            partition_file = domain_partition.partition_domains(
                pdb_id,
                chain_id,
                batch_path,
                "struct_seqid",  # Default input mode
                reference,
                concise
            )
            
            if not partition_file:
                logger.error(f"Failed to partition domains for {pdb_id}_{chain_id}")
                continue
                
            # Register partition file
            db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "domain_partition",
                    "file_path": partition_file,
                    "file_exists": True
                }
            )
            
            # Update process status
            db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": "domain_partition_complete",
                    "status": "success"
                },
                "id = %s",
                (process_id,)
            )
            
            logger.info(f"Created domain partition for {pdb_id}_{chain_id}")
            
        except Exception as e:
            logger.error(f"Error partitioning domains for {pdb_id}_{chain_id}: {e}")
            continue
    
    # Update batch status
    db.update(
        "ecod_schema.batch",
        {
            "completed_items": len(protein_rows)
        },
        "id = %s",
        (batch_id,)
    )
    
    logger.info(f"Domain analysis pipeline completed for batch {batch_id}")
    return True

def main():
    parser = argparse.ArgumentParser(description='Domain Analysis Pipeline Integration')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--concise', action='store_true',
                      help='Use concise summary/partition')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    success = run_domain_analysis_pipeline(args.batch_id, args.config, args.concise)
    
    if success:
        print(f"Domain analysis pipeline completed successfully for batch {args.batch_id}")
    else:
        print(f"Domain analysis pipeline failed for batch {args.batch_id}")
        sys.exit(1)

if __name__ == "__main__":
    main()