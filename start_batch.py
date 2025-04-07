#!/usr/bin/env python3
"""
start_batch.py - Script to start a batch processing in the PyECOD pipeline
"""

import os
import sys
import argparse
import logging
from pathlib import Path
import datetime
from typing import List, Dict, Any, Optional

# Add parent directory to path if needed
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from ecod.core.config import ConfigManager
from ecod.core.db_manager import DBManager
from ecod.core.models import Protein, Batch

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

class BatchProcessor:
    """Class to handle batch creation and processing"""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.db_config = self.config_manager.get_db_config()
        self.config = self.config_manager.config
        self.db = DBManager(self.db_config)
        self.logger = logging.getLogger("ecod.batch_processor")
        
    def get_unprocessed_proteins(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get proteins that haven't been processed yet"""
        query = """
        SELECT 
            p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.sequence
        FROM 
            ecod_schema.protein p
        JOIN
            ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        LEFT JOIN
            ecod_schema.process_status ps_status ON p.id = ps_status.protein_id
        WHERE 
            ps_status.id IS NULL
            AND ps.sequence IS NOT NULL
        ORDER BY 
            p.id
        LIMIT %s
        """
        rows = self.db.execute_dict_query(query, (limit,))
        return rows
    
    def create_batch(self, proteins: List[Dict[str, Any]], batch_type: str = "demo") -> int:
        """Create a new processing batch"""
        # Generate batch name with timestamp
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        batch_name = f"{batch_type}_batch_{timestamp}"
        
        # Get base output directory from config
        base_dir = self.config.get('paths', {}).get('output_dir', './output')
        
        # Create batch directory
        batch_path = os.path.join(base_dir, batch_name)
        os.makedirs(batch_path, exist_ok=True)
        
        # Create subdirectories
        os.makedirs(os.path.join(batch_path, "query_fastas"), exist_ok=True)
        os.makedirs(os.path.join(batch_path, "results"), exist_ok=True)
        
        # Insert batch record
        batch_id = self.db.insert(
            "ecod_schema.batch",
            {
                "batch_name": batch_name,
                "base_path": batch_path,
                "type": batch_type,
                "ref_version": self.config.get('reference', {}).get('current_version', 'demo_version'),
                "total_items": len(proteins),
                "status": "created"
            },
            "id"
        )
        
        # Register proteins in this batch
        self.register_proteins_in_batch(batch_id, batch_path, proteins)
        
        self.logger.info(f"Created batch {batch_name} with ID {batch_id} containing {len(proteins)} proteins")
        return batch_id
    
    def register_proteins_in_batch(self, batch_id: int, batch_path: str, proteins: List[Dict[str, Any]]) -> None:
        """Register proteins in a batch and create initial files"""
        fasta_dir = os.path.join(batch_path, "query_fastas")
        
        for protein in proteins:
            # Determine relative path for this protein
            pdb_id = protein['pdb_id']
            chain_id = protein['chain_id']
            rel_path = f"{pdb_id}_{chain_id}"
            
            # Register in process_status
            process_id = self.db.insert(
                "ecod_schema.process_status",
                {
                    "protein_id": protein['id'],
                    "batch_id": batch_id,
                    "current_stage": "fasta_generated",
                    "status": "pending",
                    "relative_path": rel_path
                },
                "id"
            )
            
            # Generate FASTA file
            fasta_path = os.path.join(fasta_dir, f"{protein['source_id']}.fa")
            with open(fasta_path, 'w') as f:
                f.write(f">{protein['source_id']}\n{protein['sequence']}\n")
            
            # Register FASTA file
            self.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "fasta",
                    "file_path": f"query_fastas/{protein['source_id']}.fa",
                    "file_exists": True,
                    "file_size": os.path.getsize(fasta_path),
                    "last_checked": "CURRENT_TIMESTAMP"
                }
            )
            
            self.logger.info(f"Registered protein {protein['source_id']} in batch {batch_id}")
    
    def start_processing(self, batch_id: int) -> None:
        """Start processing a batch (simulate processing)"""
        # In a real scenario, this would submit jobs to a cluster
        # For demonstration, we'll just update the status of the first item
        
        query = """
        SELECT 
            ps.id, p.pdb_id, p.chain_id
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE 
            ps.batch_id = %s
            AND ps.status = 'pending'
        LIMIT 1
        """
        
        rows = self.db.execute_dict_query(query, (batch_id,))
        if not rows:
            self.logger.warning(f"No pending items found in batch {batch_id}")
            return
        
        process_id = rows[0]['id']
        pdb_id = rows[0]['pdb_id']
        chain_id = rows[0]['chain_id']
        
        # Update process status to simulate starting processing
        self.db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "processing",
                "status": "processing",
                "updated_at": "CURRENT_TIMESTAMP"
            },
            "id = %s",
            (process_id,)
        )
        
        self.logger.info(f"Started processing protein {pdb_id}_{chain_id} in batch {batch_id}")
        
        # Create a mock job record
        job_id = self.db.insert(
            "ecod_schema.job",
            {
                "batch_id": batch_id,
                "job_type": "demo_processing",
                "slurm_job_id": "demo_job_1",
                "status": "submitted",
                "items_count": 1
            },
            "id"
        )
        
        # Link job to the process
        self.db.insert(
            "ecod_schema.job_item",
            {
                "job_id": job_id,
                "process_id": process_id
            }
        )
        
        self.logger.info(f"Created demo job {job_id} for protein {pdb_id}_{chain_id}")

def main():
    parser = argparse.ArgumentParser(description='PyECOD Batch Processing Tool')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--output-dir', type=str, 
                      help='Output directory for batch data (overrides config)')
    parser.add_argument('--limit', type=int, default=5,
                      help='Maximum number of proteins to include in batch')
    parser.add_argument('--batch-type', type=str, default='demo',
                      help='Type of batch to create (demo, blast, hhsearch)')
    parser.add_argument('--start', action='store_true',
                      help='Start processing after creating batch')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    processor = BatchProcessor(args.config)
    
    # Override output directory if specified
    if args.output_dir:
        processor.config['paths'] = processor.config.get('paths', {})
        processor.config['paths']['output_dir'] = args.output_dir
    
    # Get unprocessed proteins
    proteins = processor.get_unprocessed_proteins(args.limit)
    
    if not proteins:
        print("No unprocessed proteins found. Please insert some sample proteins first.")
        sys.exit(1)
    
    print(f"Found {len(proteins)} unprocessed proteins")
    
    # Create batch
    batch_id = processor.create_batch(proteins, args.batch_type)
    print(f"Created batch with ID: {batch_id}")
    
    # Start processing if requested
    if args.start:
        processor.start_processing(batch_id)
        print(f"Started processing batch {batch_id}")

if __name__ == "__main__":
    main()