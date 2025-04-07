#!/usr/bin/env python3
"""
finish_batch.py - Continue processing an existing batch
"""

import os
import sys
import argparse
import logging
from typing import Dict, Any, List

# Add parent directory to path if needed
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from ecod.core.config import ConfigManager
from ecod.core.db_manager import DBManager

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

class BatchProcessor:
    """Class to handle batch processing"""
    
    def __init__(self, config_path: str = None):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.db_config = self.config_manager.get_db_config()
        self.db = DBManager(self.db_config)
        self.logger = logging.getLogger("ecod.batch_processor")
    
    def get_pending_batches(self):
        """Get batches that are not completed"""
        query = """
        SELECT 
            id, batch_name, type, total_items, completed_items, status, base_path
        FROM 
            ecod_schema.batch
        WHERE 
            status != 'completed'
        ORDER BY 
            id
        """
        return self.db.execute_dict_query(query)
    
    def get_pending_items(self, batch_id: int, limit: int = 10):
        """Get pending items in a batch"""
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
            AND ps.status = 'pending'
        ORDER BY 
            ps.id
        LIMIT %s
        """
        return self.db.execute_dict_query(query, (batch_id, limit))
    
    def create_file_record(self, process_id: int, batch_path: str, relative_path: str, source_id: str):
        """Create a FASTA file and record it"""
        # Get protein sequence
        query = """
        SELECT 
            ps.sequence
        FROM 
            ecod_schema.protein p
        JOIN
            ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        JOIN
            ecod_schema.process_status pstat ON p.id = pstat.protein_id
        WHERE 
            pstat.id = %s
        """
        rows = self.db.execute_dict_query(query, (process_id,))
        if not rows:
            self.logger.error(f"No sequence found for process ID {process_id}")
            return False
        
        sequence = rows[0]['sequence']
        
        # Create directories
        fasta_dir = os.path.join(batch_path, "query_fastas")
        os.makedirs(fasta_dir, exist_ok=True)
        
        # Create FASTA file
        fasta_path = os.path.join(fasta_dir, f"{source_id}.fa")
        with open(fasta_path, 'w') as f:
            f.write(f">{source_id}\n{sequence}\n")
        
        # Register FASTA file
        file_id = self.db.insert(
            "ecod_schema.process_file",
            {
                "process_id": process_id,
                "file_type": "fasta",
                "file_path": f"query_fastas/{source_id}.fa",
                "file_exists": True,
                "file_size": os.path.getsize(fasta_path)
            },
            "id"
        )
        
        self.logger.info(f"Created file record with ID {file_id} for process {process_id}")
        return True
    
    def start_processing(self, process_id: int, batch_id: int, pdb_id: str, chain_id: str):
        """Start processing a batch item"""
        # Update process status to simulate starting processing
        self.db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "processing",
                "status": "processing"
            },
            "id = %s",
            (process_id,)
        )
        
        self.logger.info(f"Started processing protein {pdb_id}_{chain_id}")
        
        # Create a mock job record
        job_id = self.db.insert(
            "ecod_schema.job",
            {
                "batch_id": batch_id,
                "job_type": "demo_processing",
                "slurm_job_id": f"demo_job_{process_id}",
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
        
        self.logger.info(f"Created job {job_id} for process {process_id}")
        return job_id
    
    def complete_processing(self, process_id: int, batch_id: int):
        """Complete processing for a batch item"""
        # Update process status
        self.db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "completed",
                "status": "success"
            },
            "id = %s",
            (process_id,)
        )
        
        # Update job status
        query = """
        UPDATE 
            ecod_schema.job
        SET 
            status = 'completed',
            completion_time = CURRENT_TIMESTAMP
        WHERE 
            id IN (
                SELECT 
                    job_id
                FROM 
                    ecod_schema.job_item
                WHERE 
                    process_id = %s
            )
        """
        self.db.execute_query(query, (process_id,))
        
        # Update batch progress
        query = """
        UPDATE 
            ecod_schema.batch
        SET 
            completed_items = (
                SELECT 
                    COUNT(*)
                FROM 
                    ecod_schema.process_status
                WHERE 
                    batch_id = %s
                    AND status IN ('success', 'completed')
            )
        WHERE 
            id = %s
        """
        self.db.execute_query(query, (batch_id, batch_id))
        
        self.logger.info(f"Completed processing for process {process_id}")
        return True

def main():
    parser = argparse.ArgumentParser(description='Continue PyECOD Batch Processing')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='ID of the batch to process (if not specified, first pending batch is used)')
    parser.add_argument('--limit', type=int, default=10,
                      help='Maximum number of items to process')
    parser.add_argument('--complete', action='store_true',
                      help='Mark processing as completed')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    processor = BatchProcessor(args.config)
    
    # Get pending batches
    batches = processor.get_pending_batches()
    
    if not batches:
        print("No pending batches found.")
        return
    
    # Determine which batch to process
    batch_id = args.batch_id if args.batch_id else batches[0]['id']
    batch = next((b for b in batches if b['id'] == batch_id), None)
    
    if not batch:
        print(f"Batch ID {batch_id} not found or not pending.")
        return
    
    print(f"Processing batch {batch_id}: {batch['batch_name']}")
    
    # Get pending items in this batch
    items = processor.get_pending_items(batch_id, args.limit)
    
    if not items:
        print(f"No pending items found in batch {batch_id}.")
        return
    
    print(f"Found {len(items)} pending items in batch {batch_id}")
    
    for item in items:
        process_id = item['id']
        pdb_id = item['pdb_id']
        chain_id = item['chain_id']
        source_id = item['source_id']
        relative_path = item['relative_path']
        
        print(f"Processing item {process_id}: {pdb_id}_{chain_id}")
        
        # Create file record if needed
        processor.create_file_record(
            process_id, 
            batch['base_path'], 
            relative_path, 
            source_id
        )
        
        # Start processing
        job_id = processor.start_processing(
            process_id, 
            batch_id, 
            pdb_id, 
            chain_id
        )
        
        # Complete processing if requested
        if args.complete:
            processor.complete_processing(process_id, batch_id)
            print(f"Completed processing for item {process_id}")
        else:
            print(f"Started processing for item {process_id}, job ID: {job_id}")
    
    print(f"Processed {len(items)} items in batch {batch_id}")

if __name__ == "__main__":
    main()