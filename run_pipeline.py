#!/usr/bin/env python3
"""
run_pipeline.py - Orchestrate the full PyECOD processing pipeline
"""

import os
import sys
import argparse
import logging
import time
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add parent directory to path if needed
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from ecod.core.config import ConfigManager
from ecod.core.db_manager import DBManager
from ecod.core.job_manager import JobManager
from slurm_job_manager import SlurmJobManager

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

class PipelineOrchestrator:
    """Orchestrate the complete ECOD pipeline workflow"""
    
    def __init__(self, config_path: str = None):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        self.db_config = self.config_manager.get_db_config()
        self.db = DBManager(self.db_config)
        self.job_manager = JobManager(self.config)
        self.slurm_manager = SlurmJobManager(config_path)
        self.logger = logging.getLogger("ecod.pipeline")
    
    def get_unprocessed_proteins(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get proteins that haven't been processed yet"""
        # First, let's check if we can find any proteins at all
        check_query = """
        SELECT 
            p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.sequence
        FROM 
            ecod_schema.protein p
        JOIN
            ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        ORDER BY 
            p.id
        LIMIT %s
        """
        check_rows = self.db.execute_dict_query(check_query, (limit,))
        
        if check_rows:
            self.logger.info(f"Found {len(check_rows)} proteins in the database")
            # If there are proteins, look for those without completed processing
            query = """
            SELECT 
                p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.sequence
            FROM 
                ecod_schema.protein p
            JOIN
                ecod_schema.protein_sequence ps ON p.id = ps.protein_id
            LEFT JOIN (
                SELECT DISTINCT protein_id 
                FROM ecod_schema.process_status 
                WHERE status IN ('success', 'completed')
            ) ps_done ON p.id = ps_done.protein_id
            WHERE 
                ps_done.protein_id IS NULL
                AND ps.sequence IS NOT NULL
            ORDER BY 
                p.id
            LIMIT %s
            """
            rows = self.db.execute_dict_query(query, (limit,))
            return rows
        else:
            self.logger.warning("No proteins found in the database")
            return []
    
    def create_batch(self, proteins: List[Dict[str, Any]], batch_type: str = "full") -> int:
        """Create a new processing batch"""
        from datetime import datetime
        
        # Generate batch name with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M")
        batch_name = f"{batch_type}_batch_{timestamp}"
        
        # Get base output directory from config
        base_dir = self.config.get('paths', {}).get('output_dir', './output')
        
        # Create batch directory
        batch_path = os.path.join(base_dir, batch_name)
        os.makedirs(batch_path, exist_ok=True)
        
        # Create subdirectories
        os.makedirs(os.path.join(batch_path, "query_fastas"), exist_ok=True)
        os.makedirs(os.path.join(batch_path, "chain_blast_results"), exist_ok=True)
        os.makedirs(os.path.join(batch_path, "domain_blast_results"), exist_ok=True)
        os.makedirs(os.path.join(batch_path, "ecod_dump"), exist_ok=True)
        
        # Insert batch record
        batch_id = self.db.insert(
            "ecod_schema.batch",
            {
                "batch_name": batch_name,
                "base_path": batch_path,
                "type": batch_type,
                "ref_version": self.config.get('reference', {}).get('current_version', 'develop291'),
                "total_items": len(proteins),
                "status": "created"
            },
            "id"
        )
        
        # Register proteins in this batch
        self._register_proteins_in_batch(batch_id, batch_path, proteins)
        
        self.logger.info(f"Created batch {batch_name} with ID {batch_id} containing {len(proteins)} proteins")
        return batch_id
    
    def _register_proteins_in_batch(self, batch_id: int, batch_path: str, proteins: List[Dict[str, Any]]) -> None:
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
                    "file_size": os.path.getsize(fasta_path)
                }
            )
            
            self.logger.info(f"Registered protein {protein['source_id']} in batch {batch_id}")
    
    def run_blast_pipeline(self, batch_id: int, batch_size: int = 10, check_interval: int = 300) -> None:
        """Run the BLAST pipeline for a batch"""
        import subprocess
        from time import sleep
        
        # First, run chain BLAST
        self.logger.info("Running chain BLAST")
        subprocess.run([
            sys.executable, "blast_runner.py",
            "--config", self.config_manager.config_path,
            "--batch-id", str(batch_id),
            "--chain",
            "--batch-size", str(batch_size),
            "--verbose"
        ])
        
        # Wait for jobs to complete
        self._wait_for_jobs("chain_blast", batch_id, check_interval)
        
        # Then, run domain BLAST
        self.logger.info("Running domain BLAST")
        subprocess.run([
            sys.executable, "blast_runner.py",
            "--config", self.config_manager.config_path,
            "--batch-id", str(batch_id),
            "--domain",
            "--batch-size", str(batch_size),
            "--verbose"
        ])
        
        # Wait for jobs to complete
        self._wait_for_jobs("domain_blast", batch_id, check_interval)
        
        # Parse results
        self._parse_blast_results(batch_id)
    
    def run_hhsearch_pipeline(self, batch_id: int, threads: int = 8, check_interval: int = 300) -> None:
        """Run the HHSearch pipeline for a batch"""
        import subprocess
        from time import sleep
        
        # First, run HHblits profile generation
        self.logger.info("Running HHblits profile generation")
        subprocess.run([
            sys.executable, "hhsearch_runner.py",
            "--config", self.config_manager.config_path,
            "--batch-id", str(batch_id),
            "--hhblits",
            "--threads", str(threads),
            "--verbose"
        ])
        
        # Wait for jobs to complete
        self._wait_for_jobs("hhblits", batch_id, check_interval)
        
        # Then, run HHsearch
        self.logger.info("Running HHsearch")
        subprocess.run([
            sys.executable, "hhsearch_runner.py",
            "--config", self.config_manager.config_path,
            "--batch-id", str(batch_id),
            "--hhsearch",
            "--threads", str(threads),
            "--verbose"
        ])
        
        # Wait for jobs to complete
        self._wait_for_jobs("hhsearch", batch_id, check_interval)
        
        # Parse results
        self._parse_hhsearch_results(batch_id)
    
    def _wait_for_jobs(self, job_type: str, batch_id: int, check_interval: int) -> None:
        """Wait for jobs of a specific type to complete"""
        self.logger.info(f"Waiting for {job_type} jobs to complete (checking every {check_interval} seconds)")
        
        while True:
            # Check job status
            query = """
            SELECT COUNT(*) 
            FROM ecod_schema.job 
            WHERE batch_id = %s 
              AND job_type = %s 
              AND status = 'submitted'
            """
            result = self.db.execute_query(query, (batch_id, job_type))
            pending_jobs = result[0][0] if result else 0
            
            # If no more pending jobs, break
            if pending_jobs == 0:
                self.logger.info(f"All {job_type} jobs completed")
                break
            
            self.logger.info(f"{pending_jobs} {job_type} jobs still running, checking again in {check_interval} seconds")
            
            # Update job status
            self.slurm_manager.check_all_jobs(batch_id)
            
            # Wait before checking again
            time.sleep(check_interval)
    
    def _parse_blast_results(self, batch_id: int) -> None:
        """Parse BLAST results for a batch"""
        # Get processes with BLAST results
        query = """
        SELECT 
            ps.id
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 
            ps.batch_id = %s
            AND ps.status = 'success'
            AND pf.file_type IN ('chain_blast_result', 'domain_blast_result')
            AND pf.file_exists = TRUE
            AND NOT EXISTS (
                SELECT 1 FROM ecod_schema.process_file 
                WHERE process_id = ps.id AND file_type = 'blast_summary'
            )
        """
        processes = self.db.execute_query(query, (batch_id,))
        
        self.logger.info(f"Parsing BLAST results for {len(processes)} processes")
        
        import subprocess
        for process in processes:
            process_id = process[0]
            # Parse results
            subprocess.run([
                sys.executable, "blast_runner.py",
                "--config", self.config_manager.config_path,
                "--parse",
                "--process-id", str(process_id),
                "--verbose"
            ])
    
    def _parse_hhsearch_results(self, batch_id: int) -> None:
        """Parse HHSearch results for a batch"""
        # Get processes with HHSearch results
        query = """
        SELECT 
            ps.id
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 
            ps.batch_id = %s
            AND ps.status = 'success'
            AND pf.file_type = 'hhr'
            AND pf.file_exists = TRUE
            AND NOT EXISTS (
                SELECT 1 FROM ecod_schema.process_file 
                WHERE process_id = ps.id AND file_type = 'hh_summ_xml'
            )
        """
        processes = self.db.execute_query(query, (batch_id,))
        
        self.logger.info(f"Parsing HHSearch results for {len(processes)} processes")
        
        import subprocess
        for process in processes:
            process_id = process[0]
            # Parse results
            subprocess.run([
                sys.executable, "hhsearch_runner.py",
                "--config", self.config_manager.config_path,
                "--parse",
                "--process-id", str(process_id),
                "--verbose"
            ])
    
    def update_batch_status(self, batch_id: int) -> None:
        """Update batch completion status"""
        # Get completed items
        query = """
        SELECT 
            COUNT(*) 
        FROM 
            ecod_schema.process_status
        WHERE 
            batch_id = %s
            AND status IN ('success', 'completed')
        """
        result = self.db.execute_query(query, (batch_id,))
        completed_items = result[0][0] if result else 0
        
        # Get total items
        query = """
        SELECT total_items FROM ecod_schema.batch WHERE id = %s
        """
        result = self.db.execute_query(query, (batch_id,))
        total_items = result[0][0] if result else 0
        
        # Update batch status
        status = "completed" if completed_items == total_items else "processing"
        
        self.db.update(
            "ecod_schema.batch",
            {
                "completed_items": completed_items,
                "status": status
            },
            "id = %s",
            (batch_id,)
        )
        
        self.logger.info(f"Updated batch {batch_id} status: {completed_items}/{total_items} completed")
    
    def run_full_pipeline(self, limit: int = 10, batch_size: int = 10, 
                        threads: int = 8, check_interval: int = 300) -> Optional[int]:
        """Run the full pipeline"""
        # Get unprocessed proteins
        proteins = self.get_unprocessed_proteins(limit)
        
        if not proteins:
            self.logger.warning("No unprocessed proteins found.")
            return None
        
        self.logger.info(f"Found {len(proteins)} unprocessed proteins")
        
        # Create batch
        batch_id = self.create_batch(proteins, "full")
        
        # Run BLAST pipeline
        self.run_blast_pipeline(batch_id, batch_size, check_interval)
        
        # Run HHSearch pipeline
        self.run_hhsearch_pipeline(batch_id, threads, check_interval)
        
        # Update batch status
        self.update_batch_status(batch_id)
        
        self.logger.info(f"Pipeline completed for batch {batch_id}")
        return batch_id

def main():
    parser = argparse.ArgumentParser(description='PyECOD Pipeline Orchestrator')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--limit', type=int, default=10,
                      help='Maximum number of proteins to process')
    parser.add_argument('--batch-size', type=int, default=10,
                      help='Number of items per BLAST job')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads for HHblits/HHsearch')
    parser.add_argument('--check-interval', type=int, default=300,
                      help='Interval in seconds to check job status')
    parser.add_argument('--batch-id', type=int,
                      help='Use existing batch instead of creating a new one')
    parser.add_argument('--blast-only', action='store_true',
                      help='Run only the BLAST pipeline')
    parser.add_argument('--hhsearch-only', action='store_true',
                      help='Run only the HHSearch pipeline')
    parser.add_argument('--update-status', action='store_true',
                      help='Update batch status')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    orchestrator = PipelineOrchestrator(args.config)
    
    if args.batch_id:
        # Use existing batch
        batch_id = args.batch_id
        
        if args.blast_only:
            # Run only BLAST pipeline
            orchestrator.run_blast_pipeline(batch_id, args.batch_size, args.check_interval)
        
        elif args.hhsearch_only:
            # Run only HHSearch pipeline
            orchestrator.run_hhsearch_pipeline(batch_id, args.threads, args.check_interval)
        
        else:
            # Run full pipeline for existing batch
            orchestrator.run_blast_pipeline(batch_id, args.batch_size, args.check_interval)
            orchestrator.run_hhsearch_pipeline(batch_id, args.threads, args.check_interval)
        
        # Update batch status
        orchestrator.update_batch_status(batch_id)
        
    elif args.update_status and args.batch_id:
        # Just update batch status
        orchestrator.update_batch_status(args.batch_id)
    
    else:
        # Create new batch and run full pipeline
        batch_id = orchestrator.run_full_pipeline(
            args.limit, 
            args.batch_size, 
            args.threads, 
            args.check_interval
        )
        
        if batch_id:
            print(f"Pipeline completed for batch {batch_id}")
        else:
            print("Failed to create batch and run pipeline")

if __name__ == "__main__":
    main()