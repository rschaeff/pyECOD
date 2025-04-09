#!/usr/bin/env python3
"""
Domain analysis pipeline module for the ECOD pipeline
Orchestrates domain summary and partition processes
"""

import os
import logging
from typing import Dict, Any, List, Optional, Tuple

from ecod.config import ConfigManager
from ecod.db.manager import DBManager
from ecod.exceptions import PipelineError, FileOperationError

from .summary import DomainSummary
from .partition import DomainPartition


class DomainAnalysisPipeline:
    """Pipeline for domain analysis - orchestrates summary and partition processes"""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize pipeline with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.pipeline")
        
        # Initialize components
        self.summary = DomainSummary(config_path)
        self.partition = DomainPartition(config_path)
        
    def run_pipeline(self, batch_id: int, blast_only: bool = False, limit: int = None) -> bool:
        """Run the complete domain analysis pipeline for a batch
        
        Args:
            batch_id: Batch ID
            blast_only: Whether to use blast_only summaries (no HHSearch)
            limit: Maximum number of proteins to process
            
        Returns:
            True if pipeline completed successfully
        """
        self.logger.info(f"Starting domain analysis pipeline for batch {batch_id}")
        
        # Get batch information
        db_config = self.config_manager.get_db_config()
        db = DBManager(db_config)
        
        query = """
        SELECT 
            id, batch_name, base_path, ref_version
        FROM 
            ecod_schema.batch
        WHERE 
            id = %s
        """
        
        try:
            rows = db.execute_dict_query(query, (batch_id,))
            if not rows:
                self.logger.error(f"Batch {batch_id} not found")
                return False
            
            batch = rows[0]
            base_path = batch["base_path"]
            reference = batch["ref_version"]
            
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {e}")
            return False
        
        # Step 1: Generate domain summaries for the batch
        self.logger.info(f"Generating domain summaries for batch {batch_id}")
        summary_results = self._run_domain_summary(batch_id, base_path, reference, blast_only, limit)
        
        if not summary_results:
            self.logger.warning(f"No domain summaries were created for batch {batch_id}")
            return False
        
        self.logger.info(f"Created {len(summary_results)} domain summaries")
        
        # Step 2: Partition domains based on the summaries
        self.logger.info(f"Partitioning domains for batch {batch_id}")
        partition_results = self._run_domain_partition(batch_id, base_path, reference, blast_only, limit)
        
        if not partition_results:
            self.logger.warning(f"No domains were partitioned for batch {batch_id}")
            return False
        
        self.logger.info(f"Created {len(partition_results)} domain partitions")
        
        # Update batch status
        self._update_batch_status(batch_id, db)
        
        return True
    
    def _run_domain_summary(self, batch_id: int, base_path: str, reference: str, 
                          blast_only: bool = False, limit: int = None) -> List[str]:
        """Run domain summary creation for a batch
        
        Args:
            batch_id: Batch ID
            base_path: Base path for batch files
            reference: Reference version
            blast_only: Whether to use blast_only summaries
            limit: Maximum number of proteins to process
            
        Returns:
            List of created summary files
        """
        # Get database connection
        db_config = self.config_manager.get_db_config()
        db = DBManager(db_config)
        
        # Get proteins from the batch that have BLAST results
        query = """
        SELECT 
            ps.id, p.pdb_id, p.chain_id, pf1.file_exists as chain_blast_exists, 
            pf2.file_exists as domain_blast_exists,
            CASE WHEN pf3.id IS NOT NULL THEN TRUE ELSE FALSE END as summary_exists
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN
            ecod_schema.process_file pf1 ON ps.id = pf1.process_id AND pf1.file_type = 'chain_blast_result'
        LEFT JOIN
            ecod_schema.process_file pf2 ON ps.id = pf2.process_id AND pf2.file_type = 'domain_blast_result'
        LEFT JOIN
            ecod_schema.process_file pf3 ON ps.id = pf3.process_id AND pf3.file_type = 'domain_summary'
        WHERE 
            ps.batch_id = %s
            AND ps.status IN ('success', 'processing')
        """
        
        if limit:
            query += f" LIMIT {limit}"
            
        try:
            rows = db.execute_dict_query(query, (batch_id,))
        except Exception as e:
            self.logger.error(f"Error querying batch proteins for domain summary: {e}")
            return []
        
        if not rows:
            self.logger.warning(f"No proteins found for domain summary in batch {batch_id}")
            return []
        
        # Process each protein
        summary_files = []
        
        for row in rows:
            pdb_id = row["pdb_id"]
            chain_id = row["chain_id"]
            process_id = row["id"]
            summary_exists = row.get("summary_exists", false)

            # Skip if both blast files don't exist
            if not (row.get("chain_blast_exists", False) and row.get("domain_blast_exists", False)):
                missing_files_count += 1
                self.logger.warning(
                    f"Skipping {pdb_id}_{chain_id}: Missing BLAST files " 
                    f"(chain: {row.get('chain_blast_exists', False)}, "
                    f"domain: {row.get('domain_blast_exists', False)})"
                )
                continue

            try:
                summary_file = self.summary.create_summary(
                    pdb_id,
                    chain_id,
                    reference,
                    base_path,
                    blast_only
                )
                
                if summary_file:
                    summary_files.append(summary_file)
                    
                    # Update process status
                    db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "domain_summary_complete",
                            "status": "success"
                        },
                        "id = %s",
                        (row["id"],)
                    )
                    
                    # Register summary file
                    check_query = 
                    """
                    SELECT id FROM ecod_schema.process_file
                    WHERE process_id = %s AND file_type = 'domain_summary'
                    """

                    existing = db.execute_query(check_query, (process_id,))

                    if existing:
                        db.update(
                            "ecod_schema.process_file",
                            {
                                "file_path": os.path.relpath(summary_file, base_path),
                                "file_exists": True,
                                "file_size": os.path.getsize(summary_file) if os.path.exists(summary_file) else 0
                            },
                            "id = %s",
                            (existing[0][0],)
                        )
                    else:
                        db.insert(
                        "ecod_schema.process_file",
                            {
                                "process_id": row["id"],
                                "file_type": "domain_summary",
                                "file_path": os.path.relpath(summary_file, base_path),
                                "file_exists": True,
                                "file_size": os.path.getsize(summary_file) if os.path.exists(summary_file) else 0
                            }
                        )

                    
            except Exception as e:
                self.logger.error(f"Error creating domain summary for {pdb_id}_{chain_id}: {e}")
                
                # Update process status
                db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": "domain_summary_failed",
                        "status": "error",
                        "error_message": str(e)
                    },
                    "id = %s",
                    (process_id,)
                )
        self.logger.info(f"Created {len(summary_files)} domain summaries")
        self.logger.info(f"Skipped {skipped_count} existing summaries")
        self.logger.info(f"Skipped {missing_files_count} proteins with missing BLAST files")
        
        return summary_files
    
    def _run_domain_partition(self, batch_id: int, base_path: str, reference: str, 
                            blast_only: bool = False, limit: int = None) -> List[str]:
        """Run domain partition for a batch
        
        Args:
            batch_id: Batch ID
            base_path: Base path for batch files
            reference: Reference version
            blast_only: Whether to use blast_only summaries
            limit: Maximum number of proteins to process
            
        Returns:
            List of created domain partition files
        """
        # Use the partition module to process the batch
        return self.partition.process_batch(batch_id, base_path, reference, blast_only, limit)
    
    def _update_batch_status(self, batch_id: int, db: DBManager) -> None:
        """Update batch status based on completion
        
        Args:
            batch_id: Batch ID
            db: Database manager instance
        """
        # Count total and completed items
        count_query = """
        SELECT 
            COUNT(*) AS total,
            SUM(CASE WHEN status = 'success' AND current_stage = 'domain_partition_complete' THEN 1 ELSE 0 END) AS completed,
            SUM(CASE WHEN status = 'error' THEN 1 ELSE 0 END) AS failed
        FROM 
            ecod_schema.process_status
        WHERE 
            batch_id = %s
        """
        
        try:
            rows = db.execute_dict_query(count_query, (batch_id,))
            if rows:
                total = rows[0]["total"] or 0
                completed = rows[0]["completed"] or 0
                failed = rows[0]["failed"] or 0
                
                # Calculate progress
                progress = (completed + failed) / total if total > 0 else 0
                batch_status = "completed" if progress >= 1.0 else "processing"
                
                # Update batch
                db.update(
                    "ecod_schema.batch",
                    {
                        "completed_items": completed,
                        "status": batch_status
                    },
                    "id = %s",
                    (batch_id,)
                )
                
                self.logger.info(f"Updated batch {batch_id} status: {completed}/{total} completed, {failed} failed")
                
        except Exception as e:
            self.logger.error(f"Error updating batch status: {e}")
    
    def analyze_domain(self, pdb_id: str, chain_id: str, output_dir: str, 
                     reference: str, blast_only: bool = False) -> Optional[Dict[str, str]]:
        """Run domain analysis for a single protein chain
        
        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            output_dir: Output directory
            reference: Reference version
            blast_only: Whether to use blast_only summaries
            
        Returns:
            Dictionary with paths to created files
        """
        results = {}
        
        # Step 1: Create domain summary
        try:
            summary_file = self.summary.create_summary(
                pdb_id,
                chain_id,
                reference,
                output_dir,
                blast_only
            )
            
            if summary_file:
                results["summary"] = summary_file
                self.logger.info(f"Created domain summary for {pdb_id}_{chain_id}")
            else:
                self.logger.error(f"Failed to create domain summary for {pdb_id}_{chain_id}")
                return None
                
        except Exception as e:
            self.logger.error(f"Error creating domain summary for {pdb_id}_{chain_id}: {e}")
            return None
        
        # Step 2: Partition domains
        try:
            partition_file = self.partition.partition_domains(
                pdb_id,
                chain_id,
                output_dir,
                'struct_seqid',  # Default input mode
                reference,
                blast_only
            )
            
            if partition_file:
                results["partition"] = partition_file
                self.logger.info(f"Created domain partition for {pdb_id}_{chain_id}")
            else:
                self.logger.error(f"Failed to create domain partition for {pdb_id}_{chain_id}")
                
        except Exception as e:
            self.logger.error(f"Error creating domain partition for {pdb_id}_{chain_id}: {e}")
        
        return results if results else None
    
    def get_pipeline_status(self, batch_id: int) -> Dict[str, Any]:
        """Get status of domain analysis pipeline for a batch
        
        Args:
            batch_id: Batch ID
            
        Returns:
            Status information dictionary
        """
        db_config = self.config_manager.get_db_config()
        db = DBManager(db_config)
        
        # Get batch information
        batch_query = """
        SELECT 
            id, batch_name, base_path, type, ref_version,
            total_items, completed_items, status,
            created_at, completed_at
        FROM 
            ecod_schema.batch
        WHERE 
            id = %s
        """
        
        # Get stage breakdown
        stage_query = """
        SELECT 
            current_stage, status, COUNT(*) AS count
        FROM 
            ecod_schema.process_status
        WHERE 
            batch_id = %s
        GROUP BY 
            current_stage, status
        ORDER BY 
            COUNT(*) DESC
        """
        
        try:
            batch_rows = db.execute_dict_query(batch_query, (batch_id,))
            if not batch_rows:
                return {"error": f"Batch {batch_id} not found"}
            
            batch = batch_rows[0]
            
            # Get stage breakdown
            stage_rows = db.execute_dict_query(stage_query, (batch_id,))
            stages = {}
            
            for row in stage_rows:
                stage = row["current_stage"]
                status = row["status"]
                count = row["count"]
                
                if stage not in stages:
                    stages[stage] = {}
                
                stages[stage][status] = count
            
            # Calculate progress
            total = batch["total_items"] or 0
            completed = batch["completed_items"] or 0
            progress = completed / total if total > 0 else 0
            
            # Create result
            return {
                "batch_id": batch["id"],
                "batch_name": batch["batch_name"],
                "reference": batch["ref_version"],
                "status": batch["status"],
                "progress": progress,
                "total_items": total,
                "completed_items": completed,
                "created_at": batch["created_at"],
                "completed_at": batch["completed_at"],
                "stages": stages
            }
            
        except Exception as e:
            self.logger.error(f"Error getting pipeline status: {e}")
            return {"error": str(e)}