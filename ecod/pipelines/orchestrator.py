# ecod/pipelines/orchestrator.py
import logging
from typing import Dict, Any, List, Optional
from pathlib import Path
from datetime import datetime

from ..core.db_manager import DBManager
from ..core.job_manager import JobManager
from ..core.config import ConfigManager
from .blast_pipeline import BlastPipeline
from .hhsearch_pipeline import HHSearchPipeline

class PipelineOrchestrator:
    def __init__(self, config_path: Optional[str] = None):
        # Initialize configuration
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        
        # Setup logging
        self.logger = logging.getLogger("ecod.orchestrator")
        
        # Initialize managers
        self.db = DBManager(self.config_manager.get_db_config())
        self.job_manager = JobManager(self.config)
        
        # Initialize pipeline components
        self.blast = BlastPipeline(self.db, self.job_manager, self.config)
        self.hhsearch = HHSearchPipeline(self.db, self.job_manager, self.config)
        
    def run_full_pipeline(self, batch_size: int = 500, max_chains: int = None) -> None:
        """Run the complete pipeline from unclassified chains to final classification"""
        # 1. Get unclassified chains
        self.logger.info("Getting unclassified chains")
        chains = self.blast.get_unclassified_chains(max_chains)
        if not chains:
            self.logger.info("No unclassified chains found")
            return
            
        # 2. Create batch for BLAST
        self.logger.info(f"Creating BLAST batch for {len(chains)} chains")
        blast_batch = self.blast.create_batch(chains, batch_size)
        
        # 3. Run BLAST searches
        self.logger.info("Running chain-wise BLAST")
        chain_job_ids = self.blast.run_chain_blast(blast_batch.id, batch_size)
        
        self.logger.info("Running domain-wise BLAST")
        domain_job_ids = self.blast.run_domain_blast(blast_batch.id, batch_size)
        
        # 4. Create batch for HHSearch
        self.logger.info("Creating HHSearch batch")
        hhsearch_batch = self.hhsearch.create_batch(chains, batch_size)
        
        # 5. Generate profiles
        self.logger.info("Generating HHblits profiles")
        profile_job_ids = self.hhsearch.generate_profiles(hhsearch_batch.id, 8, "16G")
        
        self.logger.info(f"Pipeline started successfully!")
        self.logger.info(f"- BLAST batch ID: {blast_batch.id}, submitted {len(chain_job_ids) + len(domain_job_ids)} jobs")
        self.logger.info(f"- HHSearch batch ID: {hhsearch_batch.id}, submitted {len(profile_job_ids)} jobs")
        
    # Update the check_status method in the orchestrator
    def check_status(self, batch_id: Optional[int] = None) -> None:
        """Check status of all running jobs"""
        self.logger.info("Checking job status")
        
        # Check BLAST jobs
        self.blast.check_job_status(batch_id)
        
        # Check HHSearch jobs
        self.hhsearch.check_status(batch_id)
        
        # Update batch completion status
        self._update_batch_status(batch_id)
        
        # Parse HHSearch results for completed jobs
        self._parse_completed_hhsearch_results(batch_id)
        
    def _parse_completed_hhsearch_results(self, batch_id: Optional[int] = None) -> None:
        """Parse and generate summaries for completed HHSearch results"""
        query = """
            SELECT 
                ps.id, p.pdb_id, p.chain_id, ps.relative_path, b.base_path
            FROM 
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.batch b ON ps.batch_id = b.id
            JOIN
                ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE 
                ps.current_stage = 'hhsearch_complete'
                AND ps.status = 'success'
                AND pf.file_type = 'hhr'
                AND pf.file_exists = TRUE
                AND NOT EXISTS (
                    SELECT 1 FROM ecod_schema.process_file pf2
                    WHERE pf2.process_id = ps.id AND pf2.file_type = 'hh_summ_xml'
                )
        """
        
        if batch_id:
            query += " AND ps.batch_id = %s"
            rows = self.db.execute_dict_query(query, (batch_id,))
        else:
            rows = self.db.execute_dict_query(query)
        
        if rows:
            self.logger.info(f"Parsing {len(rows)} completed HHSearch results")
        self.hhsearch.parse_results(rows)
        
    def _update_batch_status(self, batch_id: Optional[int] = None) -> None:
        """Update batch completion status"""
        query = """
            UPDATE ecod_schema.batch
            SET completed_items = (
                SELECT COUNT(*) 
                FROM ecod_schema.process_status
                WHERE batch_id = ecod_schema.batch.id 
                AND status IN ('success', 'skipped')
            ),
            status = CASE 
                WHEN (SELECT COUNT(*) FROM ecod_schema.process_status 
                      WHERE batch_id = ecod_schema.batch.id 
                      AND status NOT IN ('success', 'skipped')) = 0 
                THEN 'completed' 
                ELSE status 
            END,
            completed_at = CASE 
                WHEN (SELECT COUNT(*) FROM ecod_schema.process_status 
                      WHERE batch_id = ecod_schema.batch.id 
                      AND status NOT IN ('success', 'skipped')) = 0 
                THEN CURRENT_TIMESTAMP 
                ELSE completed_at 
            END
        """
        
        if batch_id:
            query += " WHERE id = %s"
            self.db.execute_query(query, (batch_id,))
        else:
            query += " WHERE status != 'completed'"
            self.db.execute_query(query)
            
        # Get summary
        if batch_id:
            batch_query = """
                SELECT 
                    b.id, b.batch_name, b.status, b.total_items, b.completed_items
                FROM 
                    ecod_schema.batch b
                WHERE 
                    b.id = %s
            """
            rows = self.db.execute_dict_query(batch_query, (batch_id,))
        else:
            batch_query = """
                SELECT 
                    b.id, b.batch_name, b.status, b.total_items, b.completed_items
                FROM 
                    ecod_schema.batch b
                WHERE 
                    b.status != 'completed'
                ORDER BY 
                    b.id
            """
            rows = self.db.execute_dict_query(batch_query)
            
        for row in rows:
            self.logger.info(f"Batch {row['id']} ({row['batch_name']}): {row['completed_items']}/{row['total_items']} completed, status: {row['status']}")