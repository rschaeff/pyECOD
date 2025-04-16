#!/usr/bin/env python3
"""
Domain analysis pipeline module for the ECOD pipeline
Orchestrates domain summary and partition processes
"""

import os
import logging
from typing import Dict, Any, List, Optional, Tuple

from ecod.exceptions import PipelineError, FileOperationError
from ecod.core.context import ApplicationContext
from ecod.db import DBManager #This is probably not best practice

from ecod.pipelines.domain_analysis.summary import DomainSummary
from ecod.pipelines.domain_analysis.partition import DomainPartition


class DomainAnalysisPipeline:
    """Pipeline for domain analysis - orchestrates summary and partition processes"""
    
    def __init__(self, context=None):
        """Initialize pipeline with configuration"""
        #self.config_manager = ConfigManager(config_path)
        #self.config = self.config_manager.config
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.pipeline")
    
        self.context = context or ApplicationContext()

        # Initialize components
        self.summary = DomainSummary(self.context)
        self.partition = DomainPartition(self.context)
        
    def run_pipeline(self, batch_id: int, blast_only: bool = False, limit: int = None,
                   partition_only: bool = False, process_ids: List[int] = None) -> bool:
        """Run the complete domain analysis pipeline for a batch
        
        Args:
            batch_id: Batch ID
            blast_only: Whether to use blast_only summaries (no HHSearch)
            limit: Maximum number of proteins to process
            partition_only: Whether to run only the partition step
            process_ids: Optional list of specific process IDs to process
            
        Returns:
            True if pipeline completed successfully
        """
        self.logger.info(f"Starting domain analysis pipeline for batch {batch_id}")

        if self.context.config_manager.config.get('force_overwrite', False):
            self.logger.info(f"Force overwrite set, all pipeline data will be regenerated...")
        
        # Get batch information
        db_config = self.context.config_manager.get_db_config()
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
        
        # If process_ids are specified, handle them directly
        if process_ids:
            return self.process_proteins(batch_id, process_ids, blast_only, partition_only)
        
        # If partition only, skip domain summary generation
        if partition_only:
            # Verify summaries are complete
            if not self._verify_summary_completion(batch_id, db):
                self.logger.error(f"Cannot run partition only: Domain summaries are incomplete")
                return False
                
            # Run partition step
            self.logger.info(f"Running partition only for batch {batch_id}")
            partition_results = self.partition.process_batch(batch_id, base_path, reference, blast_only, limit)
            
            if not partition_results:
                self.logger.warning(f"No domains were partitioned for batch {batch_id}")
                return False
            
            self.logger.info(f"Created {len(partition_results)} domain partitions")
            return True
        
        # Standard pipeline - run both summary and partition
        # Step 1: Generate domain summaries for the batch
        self.logger.info(f"Generating domain summaries for batch {batch_id}")
        summary_results = self._run_domain_summary(batch_id, base_path, reference, blast_only, limit)
        
        if not summary_results:
            self.logger.warning(f"No domain summaries were created for batch {batch_id}")
            return False
        
        self.logger.info(f"Created {len(summary_results)} domain summaries")
        
        # Step 2: Partition domains based on the summaries
        self.logger.info(f"Partitioning domains for batch {batch_id}")
        partition_results = self.partition.process_batch(batch_id, base_path, reference, blast_only, limit)
        
        if not partition_results:
            self.logger.warning(f"No domains were partitioned for batch {batch_id}")
            return False
        
        self.logger.info(f"Created {len(partition_results)} domain partitions")
        
        # Update batch status
        self._update_batch_status(batch_id, db)
        
        return True

    def process_proteins(self, batch_id: int, protein_ids: List[int], 
                       blast_only: bool = False, partition_only: bool = False) -> bool:
        """Process domain analysis for specific proteins
        
        Args:
            batch_id: Batch ID
            protein_ids: List of protein IDs to process
            blast_only: Whether to use only BLAST results (no HHSearch)
            partition_only: Whether to run only the partition step
            
        Returns:
            True if successful
        """
        # Get database connection
        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)
        
        # Get batch information
        batch_query = """
        SELECT base_path, ref_version FROM ecod_schema.batch WHERE id = %s
        """
        batch_info = db.execute_dict_query(batch_query, (batch_id,))
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return False
            
        base_path = batch_info[0]['base_path']
        reference = batch_info[0]['ref_version']
        
        # Get process IDs for the specified proteins
        # If protein_ids is actually process_ids, adjust the query
        if isinstance(protein_ids[0], int) and protein_ids[0] > 10000:  # Heuristic: process IDs tend to be large
            self.logger.info("Treating input as process IDs rather than protein IDs")
            process_query = """
            SELECT ps.id, p.id as protein_id, p.pdb_id, p.chain_id 
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE ps.id IN %s AND ps.batch_id = %s
            """
            process_rows = db.execute_dict_query(process_query, (tuple(protein_ids), batch_id))
            process_ids = protein_ids
        else:
            process_query = """
            SELECT ps.id, p.id as protein_id, p.pdb_id, p.chain_id 
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.id IN %s AND ps.batch_id = %s
            """
            process_rows = db.execute_dict_query(process_query, (tuple(protein_ids), batch_id))
            process_ids = [row['id'] for row in process_rows]
        
        if not process_rows:
            self.logger.warning(f"No matching processes found for specified IDs in batch {batch_id}")
            return False
            
        # Log what we're processing
        proteins_str = ", ".join([f"{row['pdb_id']}_{row['chain_id']}" for row in process_rows[:5]])
        if len(process_rows) > 5:
            proteins_str += f" and {len(process_rows) - 5} more"
            
        self.logger.info(f"Processing {len(process_rows)} proteins: {proteins_str}")
        
        # If partition only, run only partition step
        if partition_only:
            return self.partition.process_specific_ids(
                batch_id, process_ids, base_path, reference, blast_only
            )
            
        # Standard process - run both summary and partition
        success_count = 0
        
        for row in process_rows:
            pdb_id = row['pdb_id']
            chain_id = row['chain_id']
            process_id = row['id']
            
            try:
                # Step 1: Create domain summary
                summary_file = self.summary.create_summary(
                    pdb_id, chain_id, reference, base_path, blast_only
                )
                
                if not summary_file:
                    self.logger.error(f"Failed to create domain summary for {pdb_id}_{chain_id}")
                    continue
                    
                # Step 2: Partition domains
                domain_file = self.partition.partition_domains(
                    pdb_id, chain_id, base_path, 'struct_seqid', reference, blast_only
                )
                
                if domain_file:
                    success_count += 1
                    
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
                    
                    check_query = """
                    SELECT id FROM ecod_schema.process_file
                    WHERE process_id= %s AND file_type = %s
                    """

                    existing = db.execute_query(check_query, (process_id, file_type))

                    if existing:
                        db.update(
                            "ecod_schema.process_file",
                            {
                                "file_path": domain_rel_path,
                                "file_exists": True,
                                "file_size": os.path.getsize(domain_file) if os.path.exists(domain_file) else 0
                            })
                    else:     
                        # Register domain file
                        domain_rel_path = os.path.relpath(domain_file, base_path)
                        db.insert(
                            "ecod_schema.process_file",
                            {
                                "process_id": process_id,
                                "file_type": "domain_partition",
                                "file_path": domain_rel_path,
                                "file_exists": True,
                                "file_size": os.path.getsize(domain_file) if os.path.exists(domain_file) else 0
                            }
                        )
                    
            except Exception as e:
                self.logger.error(f"Error processing domains for {pdb_id}_{chain_id}: {e}")
        
        self.logger.info(f"Successfully processed {success_count}/{len(process_rows)} proteins")
        return success_count > 0
    
    def _verify_summary_completion(self, batch_id: int, context: ApplicationContext) -> bool:
        """Verify that domain summaries are complete for a batch
        
        Args:
            batch_id: Batch ID to check
            db: Database manager instance
            
        Returns:
            True if all proteins have summaries OR force_overwrite is True
        """
        # Check if force_overwrite is set - if so, proceed regardless of summary status
        if self.context.config_manager.config.get('force_overwrite', False):
            self.logger.info("Force overwrite enabled, proceeding with partition regardless of summary status")
            return True
            
        # Query to check summary completion
        query = """
        SELECT 
            COUNT(*) as total,
            SUM(CASE WHEN EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                AND pf.file_type = 'domain_summary'
                AND pf.file_exists = TRUE
            ) THEN 1 ELSE 0 END) as complete_count
        FROM 
            ecod_schema.process_status ps
        WHERE 
            ps.batch_id = %s
        """
        
        results = context.db_manager.db.execute_dict_query(query, (batch_id,))[0]
        total = results.get('total', 0)
        complete = results.get('complete_count', 0)
        is_complete = total > 0 and total == complete
        
        self.logger.info(f"Summary completion: {complete}/{total} ({is_complete})")
        
        if not is_complete:
            self.logger.warning(f"Domain summaries incomplete: {complete}/{total}")
            
        return is_complete


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
        db_config = self.context.config_manager.get_db_config()
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
        skipped_count = 0
        missing_files_count = 0;
        
        for row in rows:
            pdb_id = row["pdb_id"]
            chain_id = row["chain_id"]
            process_id = row["id"]
            summary_exists = row.get("summary_exists", False)

            # Skip if both blast files don't exist
            if not (row.get("chain_blast_exists", False) and row.get("domain_blast_exists", False)):
                missing_files_count += 1
                self.logger.warning(
                    f"Skipping {pdb_id}_{chain_id}: Missing BLAST files " 
                    f"(chain: {row.get('chain_blast_exists', False)}, "
                    f"domain: {row.get('domain_blast_exists', False)})"
                )
                skipped_count += 1
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
                    check_query = """
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
        db_config = self.context.config_manager.get_db_config()
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