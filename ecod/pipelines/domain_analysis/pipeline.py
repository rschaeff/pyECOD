#!/usr/bin/env python3
"""
Domain analysis pipeline module for the ECOD pipeline
Orchestrates domain summary and partition processes
"""

import os
import logging
from typing import Dict, Any, List, Optional, Tuple

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.exceptions import PipelineError, FileOperationError
from ecod.pipelines.domain_analysis.summary import DomainSummaryService
from ecod.pipelines.domain_analysis.partition import DomainPartition
from ecod.models import PipelineResult, DomainSummaryModel, ProteinResult, ProteinProcessingResult

class DomainAnalysisPipeline:
    """Pipeline for domain analysis - orchestrates summary and partition processes"""
    
    def __init__(self, context=None):
        """
        Initialize pipeline with application context
        
        Args:
            config_path: Optional path to configuration file
        """
        # Initialize application context
        self.context = context or ApplicationContext()
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.pipeline")
        
        # Initialize components with shared context
        self.summary = DomainSummary(self.context)
        self.partition = DomainPartition(self.context)
        
    def run_pipeline(self, batch_id: int, blast_only: bool = False, limit: int = None,
                   partition_only: bool = False, process_ids: List[int] = None,
                   reps_only: bool = False, reset_failed: bool = True
    ) -> PipelineResult:
        """
        Run the complete domain analysis pipeline for a batch

        Args:
            batch_id: Batch ID to process
            blast_only: Whether to use only BLAST results (no HHSearch)
            limit: Maximum number of proteins to process
            partition_only: Whether to run only the partition step
            process_ids: Specific process IDs to run
            reps_only: Whether to process only representative proteins
            reset_failed: Whether to reset failed processes before running

        Returns:
            PipelineResult object containing success status and detailed statistics
        """
        self.logger.info(f"Starting domain analysis pipeline for batch {batch_id}")

        # Initialize pipeline result object
        result = PipelineResult(batch_id=batch_id)

        # Reset failed processes if requested
        if reset_failed:
            self._reset_failed_processes(batch_id, result)

        # Get batch information
        batch_info = self._get_batch_info(db, batch_id)
        if not batch_info:
            result.set_error("Batch not found")
            return result

        # Set batch info in result
        result.set_batch_info(batch_info)

        # If specific process_ids are provided, handle them directly
        if process_ids:
            return self._process_specific_proteins(batch_id, process_ids, blast_only,
                                                  partition_only, reps_only, result)

        # If partition only, verify summaries are complete
        if partition_only:
            if not self._verify_summary_completion(batch_id, blast_only):
                result.set_error("Domain summaries incomplete")
                return result

            return self._run_partition_only(batch_id, batch_info, blast_only, limit, reps_only, result)

        # Run complete pipeline (both summary and partition)
        process_rows = self._get_proteins_for_processing(batch_id, limit, reps_only)
        if not process_rows:
            result.set_error("No proteins found for processing")
            return result

        # Process each protein
        return self._process_all_proteins(process_rows, batch_info, blast_only, result)

    def _reset_failed_processes(self, batch_id: int, result: PipelineResult) -> None:
        """Reset failed processes and update result statistics"""
        self.logger.info("Resetting failed processes...")

        # Reset domain partition failures
        reset_count = self.reset_failed_processes(batch_id, 'domain_partition_failed')
        result.reset_stats["domain_partition_failed"] = reset_count

        # Reset domain summary failures
        reset_count = self.reset_failed_processes(batch_id, 'domain_summary_failed')
        result.reset_stats["domain_summary_failed"] = reset_count

    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information from database"""
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
            return rows[0] if rows else None
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {e}")
            return None

    def _process_specific_proteins(self, batch_id: int, process_ids: List[int], blast_only: bool,
                                 partition_only: bool, reps_only: bool,
                                 result: PipelineResult
    ) -> PipelineResult:
        """Process specific proteins by process ID"""
        self.logger.info(f"Processing {len(process_ids)} specific proteins")

        # Delegate to dedicated method
        processing_result = self.process_proteins(batch_id, process_ids, blast_only, partition_only, reps_only)

        # Update result from processing output
        if isinstance(processing_result, dict):
            result.success = processing_result.get("success", False)
            result.processing_stats = processing_result.get("processing_stats", {})
            if not result.success:
                result.set_error(processing_result.get("error", "Unknown error"))
        else:
            # Handle boolean return for backwards compatibility
            result.success = bool(processing_result)

        return result

    def _run_partition_only(self, batch_id: int, batch_info: Dict[str, Any],
                           blast_only: bool, limit: int, reps_only: bool,
                           result: PipelineResult
    ) -> PipelineResult:
        """Run only the partition step of the pipeline"""
        self.logger.info(f"Running partition only for batch {batch_id}")

        # Use the dedicated partition method from DomainPartition
        partition_result = self.partition.process_batch(
            batch_id,
            batch_info['base_path'],
            batch_info['ref_version'],
            blast_only,
            limit,
            reps_only
        )

        # Process the result
        if isinstance(partition_result, dict):
            result.success = partition_result.get("success", False)
            result.partition_stats = partition_result.get("stats", {})
            if not result.success:
                result.set_error(partition_result.get("error", "Unknown partition error"))
        else:
            # Handle legacy boolean/list return
            result.success = bool(partition_result)
            if isinstance(partition_result, list):
                result.partition_stats["files_created"] = len(partition_result)

            if not result.success:
                result.set_error("No domains were partitioned")

        return result

    def _get_proteins_for_processing(self, batch_id: int, limit: int = None,
                                   reps_only: bool = False
    ) -> List[ProteinResult]:
        """
        Get proteins that need processing from the database

        Args:
            batch_id: Batch ID to query
            limit: Maximum number of proteins to return
            reps_only: Whether to return only representative proteins

        Returns:
            List of ProteinResult objects with initial metadata
        """
        from ecod.models.pipeline import ProteinResult

        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)

        # Query needs to include protein_id for ProteinResult model
        query = """
        SELECT
            ps.id as process_id,
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.relative_path
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE
            ps.batch_id = %s
            AND ps.status IN ('success', 'processing')
        """

        if reps_only:
            query += " AND ps.is_representative = TRUE"

        if limit:
            query += f" LIMIT {limit}"

        try:
            rows = db.execute_dict_query(query, (batch_id,))

            # Convert database rows to ProteinResult objects
            result = []
            for row in rows:
                protein = ProteinResult(
                    protein_id=row['protein_id'],
                    pdb_id=row['pdb_id'],
                    chain_id=row['chain_id'],
                    process_id=row['process_id'],
                    # Initialize with default values for processing fields
                    success=False,
                    summary_success=False,
                    partition_success=False
                )
                result.append(protein)

            self.logger.info(f"Found {len(result)} proteins for processing in batch {batch_id}")
            return result

        except Exception as e:
            self.logger.error(f"Error querying proteins for batch {batch_id}: {e}")
            return []

    def _process_all_proteins(self, process_rows: List[Dict[str, Any]],
                            batch_info: Dict[str, Any], blast_only: bool,
                            result: PipelineResult
    ) -> PipelineResult:
        """Process all proteins in the batch"""
        base_path = batch_info['base_path']
        reference = batch_info['ref_version']

        # Initialize statistics
        summary_stats = {
            "total_proteins": len(process_rows),
            "summary_success": 0,
            "chain_blast_processed": 0,
            "domain_blast_processed": 0,
            "hhsearch_processed": 0,
            "total_hhsearch_hits": 0,
            "total_chain_blast_hits": 0,
            "total_domain_blast_hits": 0,
            "skipped_count": 0,
            "missing_blast_count": 0,
            "errors": []
        }

        partition_stats = {
            "total_proteins": len(process_rows),
            "partition_success": 0,
            "domains_created": 0,
            "errors": []
        }

        # Get database connection
        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)

        # Process each protein
        summary_models = []

        for row in process_rows:
            pdb_id = row['pdb_id']
            chain_id = row['chain_id']
            process_id = row['id']

            summary_model = self._process_single_protein(
                process_id, pdb_id, chain_id, base_path, reference,
                blast_only, db, summary_stats
            )

            if summary_model:
                summary_models.append(summary_model)

        # Update result with summary statistics
        result.summary_stats = summary_stats

        # Check if any summaries were created
        if not summary_models:
            result.set_error("No domain summaries were created")
            return result

        # Run domain partition based on summaries
        partition_result = self._run_domain_partition(
            base_path, reference, blast_only, partition_stats, db
        )

        result.partition_stats = partition_stats
        result.success = True

        return result

    def _process_single_protein(self, process_id: int, pdb_id: str, chain_id: str,
                              base_path: str, reference: str, blast_only: bool,
                              db: DBManager, summary_stats: Dict[str, Any]
    ) -> Optional[DomainSummaryModel]:
        """Process a single protein chain"""
        try:
            # Create domain summary using the model-based approach
            summary_model = self.summary.create_summary(
                pdb_id, chain_id, reference, base_path, blast_only
            )

            # Update statistics from model
            stats = summary_model.get_stats()
            for key, value in stats.items():
                if key in summary_stats:
                    if isinstance(value, bool) and value:
                        summary_stats[key] += 1
                    elif isinstance(value, (int, float)):
                        summary_stats[key] += value

            if summary_model.skipped:
                summary_stats['skipped_count'] += 1
                return summary_model

            if not summary_model.output_file_path:
                self.logger.error(f"Missing output file path for {pdb_id}_{chain_id}")
                summary_stats['errors'].append(f"Missing output path for {pdb_id}_{chain_id}")
                return None

            # Register the file
            summary_file = summary_model.output_file_path
            self._register_summary_file(process_id, os.path.relpath(summary_file, base_path), db)

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

            summary_stats['summary_success'] += 1
            return summary_model

        except Exception as e:
            self.logger.error(f"Error creating domain summary for {pdb_id}_{chain_id}: {e}")
            summary_stats['errors'].append(f"Error for {pdb_id}_{chain_id}: {str(e)}")

            # Update process status for error
            db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": "domain_summary_failed",
                    "status": "error",
                    "error_message": str(e)[:500]  # Truncate to avoid DB field size issues
                },
                "id = %s",
                (process_id,)
            )

            return None

    def _register_summary_file(self, process_id: int, file_path: str, db: DBManager) -> bool:
        """Register domain summary file in database"""
        try:
            # Check if record already exists
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = 'domain_summary'
            """
            
            existing = db.execute_query(query, (process_id,))
            
            if existing:
                # Update existing record
                db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": file_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    },
                    "id = %s",
                    (existing[0][0],)
                )
            else:
                # Insert new record
                db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "domain_summary",
                        "file_path": file_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    }
                )
            return True
        except Exception as e:
            self.logger.warning(f"Error registering summary file: {e}")
            return False

    def process_proteins(self, batch_id: int, protein_ids: List[int], 
                        blast_only: bool = False, partition_only: bool = False,
                        reps_only: bool = False
    ) -> ProteinProcessingResult:
        """
        Process domain analysis for specific proteins

        Args:
            batch_id: Batch ID
            protein_ids: List of protein IDs to process
            blast_only: Whether to use only BLAST results (no HHSearch)
            partition_only: Whether to run only the partition step
            reps_only: Whether to process only representative proteins

        Returns:
            ProteinProcessingResult object with processing statistics and status
        """
        # Create result object
        result = ProteinProcessingResult(batch_id=batch_id, protein_ids=protein_ids)

        # Get database connection
        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)

        # Get batch information
        batch_info = self._get_batch_info(db, batch_id)
        if not batch_info:
            result.set_error("Batch not found")
            return result

        # Store batch info in result
        result.batch_info = batch_info

        # Get process information for specified proteins
        process_rows = self._get_process_info(db, batch_id, protein_ids, reps_only)
        if not process_rows:
            result.set_error("No matching processes found for specified protein IDs")
            return result

        # Log what we're processing
        self._log_processing_summary(process_rows)

        # If partition only, run only partition step
        if partition_only:
            return self._process_partition_only(
                db, batch_id, process_rows, batch_info, blast_only, result
            )

        # Standard process - run both summary and partition
        return self._process_full_pipeline(
            db, process_rows, batch_info, blast_only, result
        )

    def _get_batch_info(self, db: DBManager, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information from database"""
        batch_query = """
        SELECT base_path, ref_version FROM ecod_schema.batch WHERE id = %s
        """

        try:
            batch_info = db.execute_dict_query(batch_query, (batch_id,))
            if batch_info:
                return batch_info[0]
            self.logger.error(f"Batch {batch_id} not found")
            return None
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {e}")
            return None

    def _get_process_info(self, db: DBManager, batch_id: int,
                        protein_ids: List[int], reps_only: bool
    ) -> List[Dict[str, Any]]:
        """Get process information for specified proteins"""
        process_query = """
        SELECT ps.id, p.id as protein_id, p.pdb_id, p.chain_id
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.id IN %s AND ps.batch_id = %s
        """

        if reps_only:
            process_query += " AND ps.is_representative = TRUE"

        try:
            return db.execute_dict_query(process_query, (tuple(protein_ids), batch_id))
        except Exception as e:
            self.logger.error(f"Error querying process information: {e}")
            return []

    def _log_processing_summary(self, process_rows: List[Dict[str, Any]]) -> None:
        """Log summary of proteins being processed"""
        proteins_str = ", ".join([f"{row['pdb_id']}_{row['chain_id']}" for row in process_rows[:5]])
        if len(process_rows) > 5:
            proteins_str += f" and {len(process_rows) - 5} more"

        self.logger.info(f"Processing {len(process_rows)} proteins: {proteins_str}")

    def _process_partition_only(self, db: DBManager, batch_id: int,
                              process_rows: List[Dict[str, Any]],
                              batch_info: Dict[str, Any], blast_only: bool,
                              result: ProteinProcessingResult
    ) -> ProteinProcessingResult:
        """Process only the partition step for specified proteins"""
        self.logger.info("Running partition only")
        base_path = batch_info['base_path']
        reference = batch_info['ref_version']

        # Extract process IDs
        process_ids = [row['id'] for row in process_rows]

        # Run partition through the DomainPartition module
        partition_success = self.partition.process_specific_ids(
            batch_id, process_ids, base_path, reference, blast_only
        )

        # Update result
        result.success = bool(partition_success)
        if not result.success:
            result.set_error("Partition process failed")

        return result

    def _process_full_pipeline(self, db: DBManager, process_rows: List[Dict[str, Any]],
                             batch_info: Dict[str, Any], blast_only: bool,
                             result: ProteinProcessingResult
    ) -> ProteinProcessingResult:
        """Process both summary and partition steps for specified proteins"""
        base_path = batch_info['base_path']
        reference = batch_info['ref_version']

        # Initialize tracking
        successful_proteins = []

        # Process each protein
        for row in process_rows:
            pdb_id = row['pdb_id']
            chain_id = row['chain_id']
            process_id = row['id']
            protein_id = row['protein_id']

            # Track this protein in result
            protein_result = ProteinResult(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                process_id=process_id
            )
            result.protein_results.append(protein_result)

            try:
                # Step 1: Create domain summary with standardized naming
                summary_model = self._process_summary_step(
                    pdb_id, chain_id, reference, base_path, blast_only, protein_result
                )

                if not summary_model:
                    continue

                # Step 2: Partition domains with standardized naming
                domain_file = self._process_partition_step(
                    pdb_id, chain_id, base_path, reference, blast_only, protein_result
                )

                if domain_file:
                    successful_proteins.append(protein_id)
                    protein_result.success = True

                    # Update process status
                    self._update_process_status(db, process_id, "domain_partition_complete", "success")

                    # Register domain file
                    self._register_domain_file(db, process_id, domain_file, base_path)

            except Exception as e:
                error_message = f"Error processing domains for {pdb_id}_{chain_id}: {e}"
                self.logger.error(error_message)

                # Update protein result
                protein_result.success = False
                protein_result.error = str(e)

                # Update process status for error
                self._update_process_status(db, process_id, "domain_partition_failed", "error", error_message)

        # Update overall result
        result.success = len(successful_proteins) > 0
        result.success_count = len(successful_proteins)
        result.total_count = len(process_rows)

        self.logger.info(f"Successfully processed {result.success_count}/{result.total_count} proteins")
        return result

    def _process_summary_step(self, pdb_id: str, chain_id: str, reference: str,
                            base_path: str, blast_only: bool,
                            protein_result: ProteinResult
    ) -> Optional[DomainSummaryModel]:
        """Process the summary step for a single protein"""
        try:
            # Create domain summary
            summary_model = self.summary.create_summary(
                pdb_id, chain_id, reference, base_path, blast_only
            )

            # Store summary file path in result
            protein_result.summary_file = summary_model.output_file_path

            # Updated property names based on the new model
            protein_result.summary_success = (
                summary_model.processed and not summary_model.skipped
            )

            # Track statistics
            protein_result.chain_blast_count = len(summary_model.chain_blast_hits)
            protein_result.domain_blast_count = len(summary_model.domain_blast_hits)
            protein_result.hhsearch_count = len(summary_model.hhsearch_hits)

            return summary_model
        except Exception as e:
            error_message = f"Failed to create domain summary: {e}"
            self.logger.error(error_message)
            protein_result.summary_success = False
            protein_result.summary_error = error_message
            return None

    def _process_partition_step(self, pdb_id: str, chain_id: str, base_path: str,
                              reference: str, blast_only: bool,
                              protein_result: ProteinResult
    ) -> Optional[str]:
        """Process the partition step for a single protein"""
        try:
            # Partition domains
            domain_file = self.partition.partition_domains(
                pdb_id, chain_id, base_path, 'struct_seqid', reference, blast_only
            )

            # Store domain file path in result
            protein_result.domain_file = domain_file
            protein_result.partition_success = bool(domain_file)

            return domain_file
        except Exception as e:
            error_message = f"Failed to partition domains: {e}"
            self.logger.error(error_message)
            protein_result.partition_success = False
            protein_result.partition_error = error_message
            return None

    def _update_process_status(self, db: DBManager, process_id: int,
                             stage: str, status: str, error_message: str = None
    ) -> None:
        """Update process status in database"""
        update_data = {
            "current_stage": stage,
            "status": status
        }

        if error_message:
            update_data["error_message"] = error_message[:500]  # Truncate to avoid field size issues

        try:
            db.update(
                "ecod_schema.process_status",
                update_data,
                "id = %s",
                (process_id,)
            )
        except Exception as e:
            self.logger.error(f"Error updating process status: {e}")

    def _register_domain_file(self, db: DBManager, process_id: int,
                            domain_file: str, base_path: str
    ) -> None:
        """Register domain file in database"""
        try:
            domain_rel_path = os.path.relpath(domain_file, base_path)
            file_type = "domain_partition"

            # Check if domain file is already registered
            check_query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """

            existing = db.execute_query(check_query, (process_id, file_type))

            if existing:
                # Update existing record
                db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": domain_rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(domain_file) if os.path.exists(domain_file) else 0
                    },
                    "id = %s",
                    (existing[0][0],)
                )
            else:
                # Insert new record
                db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": domain_rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(domain_file) if os.path.exists(domain_file) else 0
                    }
                )
        except Exception as e:
            self.logger.error(f"Error registering domain file: {e}")

    def _verify_summary_completion(self, batch_id: int, blast_only: bool = False) -> bool:
        """Verify that domain summaries are complete for a batch"""
        # Check if force_overwrite is set
        if self.context.config_manager.config.get('force_overwrite', False):
            self.logger.info("Force overwrite enabled, proceeding regardless of summary status")
            return True
            
        # Add blast_only suffix to file_type for precise matching
        summary_type = "domain_summary"
        if blast_only:
            summary_type = "domain_summary_blast_only"
    

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
        
        results = self.context.db.execute_dict_query(query, (batch_id,))[0]
        total = results.get('total', 0)
        complete = results.get('complete_count', 0)
        is_complete = total > 0 and total == complete
        
        self.logger.info(f"Summary completion: {complete}/{total} ({is_complete})")
        
        if is_complete:
            # Verify physical files exist
            self.logger.info("Database indicates summaries complete, checking actual files...")
            # Get a few samples to check
            sample_query = """
            SELECT pf.file_path, b.base_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.batch b ON ps.batch_id = b.id
            WHERE ps.batch_id = %s AND pf.file_type = 'domain_summary' AND ps.is_representative = TRUE
            LIMIT 5
            """
            sample_files = self.context.db.execute_dict_query(sample_query, (batch_id,))
            
            # Check if sample files exist
            for sample in sample_files:
                file_path = os.path.join(sample['base_path'], sample['file_path'])
                if not os.path.exists(file_path):
                    self.logger.warning(f"File marked as existing does not exist: {file_path}")
                    return False  # Don't run partition-only if files don't actually exist
        
        return is_complete

    def _run_domain_summary(self, batch_id: int, base_path: str, reference: str, 
                          blast_only: bool = False, limit: int = None, reps_only: bool = None
    ) -> List[str]:
        """Run domain summary creation for a batch
        
        Args:
            batch_id: Batch ID
            base_path: Base path for batch files
            reference: Reference version
            blast_only: Whether to use blast_only summaries
            limit: Maximum number of proteins to process
            reps_only: Whether to process only representative proteins
            
        Returns:
            List of created summary files
        """
        # Get database connection
        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)
        
        self.logger.info(f"Running domain summary for batch {batch_id} (blast_only={blast_only})")
        
        # Get proteins from the batch that have BLAST results
        query = """
        SELECT 
            ps.id, p.pdb_id, p.chain_id, pf1.file_exists as chain_blast_exists, 
            pf2.file_exists as domain_blast_exists,
            CASE WHEN pf3.id IS NOT NULL THEN TRUE ELSE FALSE END as summary_exists,
            CASE WHEN pf4.id IS NOT NULL THEN TRUE ELSE FALSE END as hhsearch_exists
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
        LEFT JOIN
            ecod_schema.process_file pf4 ON ps.id = pf4.process_id AND 
            pf4.file_type IN ('hhsearch_result', 'hhr', 'hhsearch_xml')
        WHERE 
            ps.batch_id = %s
            AND ps.status IN ('success', 'processing')
        """

        if reps_only:
            query += f" AND ps.is_representative IS TRUE"
            self.logger.info("Filtering for representative proteins (processes) only")
        
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
        
        # Log summary of what we're about to process
        self.logger.info(f"Found {len(rows)} proteins to process for batch {batch_id}")
        
        # Process each protein
        summary_files = []
        skipped_count = 0
        missing_blast_count = 0
        missing_hhsearch_count = 0
        hhsearch_processed_count = 0
        
        # Prepare a progress logger that logs every 10% or every 10 proteins, whichever is smaller
        total_proteins = len(rows)
        log_interval = max(1, min(total_proteins // 10, 10))
        
        for idx, row in enumerate(rows):
            pdb_id = row["pdb_id"]
            chain_id = row["chain_id"]
            process_id = row["id"]
            summary_exists = row.get("summary_exists", False)
            hhsearch_exists = row.get("hhsearch_exists", False)

            # Skip if both blast files don't exist
            if not (row.get("chain_blast_exists", False) and row.get("domain_blast_exists", False)):
                missing_blast_count += 1
                self.logger.warning(
                    f"Skipping {pdb_id}_{chain_id}: Missing BLAST files " 
                    f"(chain: {row.get('chain_blast_exists', False)}, "
                    f"domain: {row.get('domain_blast_exists', False)})"
                )
                continue
                
            # Log HHSearch file availability when not in blast_only mode
            if not blast_only:
                if not hhsearch_exists:
                    missing_hhsearch_count += 1
                    self.logger.warning(f"No HHSearch files in database for {pdb_id}_{chain_id}")

            # Log progress periodically
            if idx % log_interval == 0:
                self.logger.info(f"Processing protein {idx+1}/{total_proteins} ({(idx+1)/total_proteins:.1%})")

            try:
                summary_result = self.summary.create_summary(
                    pdb_id,
                    chain_id,
                    reference,
                    base_path,
                    blast_only
                )
                
                if isinstance(summary_result, dict):
                    # Enhanced return value with stats
                    summary_file = summary_result.get('file_path', '')
                    stats = summary_result.get('stats', {})
                    
                    if stats.get('hhsearch_processed', False) and not blast_only:
                        hhsearch_processed_count += 1
                        self.logger.info(
                            f"Processed {pdb_id}_{chain_id} with HHSearch: "
                            f"{stats.get('hhsearch_hits', 0)} hits"
                        )
                else:
                    # Backward compatibility for string return
                    summary_file = summary_result
                    
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
                        (process_id,)
                    )
                    
                    # Register summary file in database
                    file_path = os.path.relpath(summary_file, base_path)
                    file_type = "domain_summary"
                    file_size = os.path.getsize(summary_file) if os.path.exists(summary_file) else 0
                    
                    check_query = """
                    SELECT id FROM ecod_schema.process_file
                    WHERE process_id = %s AND file_type = %s
                    """
                    existing = db.execute_query(check_query, (process_id, file_type))

                    if existing:
                        # Update existing record
                        db.update(
                            "ecod_schema.process_file",
                            {
                                "file_path": file_path,
                                "file_exists": True,
                                "file_size": file_size
                            },
                            "id = %s",
                            (existing[0][0],)
                        )
                    else:
                        # Insert new record
                        db.insert(
                            "ecod_schema.process_file",
                            {
                                "process_id": process_id,
                                "file_type": file_type,
                                "file_path": file_path,
                                "file_exists": True,
                                "file_size": file_size
                            }
                        )
                    
            except Exception as e:
                self.logger.error(f"Error creating domain summary for {pdb_id}_{chain_id}: {e}")
                
                # Update process status for error
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
        
        # Log comprehensive statistics
        self.logger.info(f"Domain summary statistics for batch {batch_id}:")
        self.logger.info(f"- Total proteins: {len(rows)}")
        self.logger.info(f"- Created {len(summary_files)} domain summaries")
        self.logger.info(f"- Skipped {skipped_count} existing summaries")
        self.logger.info(f"- Skipped {missing_blast_count} proteins with missing BLAST files")
        
        if not blast_only:
            self.logger.info(f"- Missing HHSearch files: {missing_hhsearch_count}")
            self.logger.info(f"- Successfully processed with HHSearch: {hhsearch_processed_count}")
        
        return summary_files
    
    def _run_domain_partition(self, batch_id: int, base_path: str, reference: str, 
                            blast_only: bool = False, limit: int = None
    ) -> List[str]:
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
                     reference: str, blast_only: bool = False
    ) -> Dict[str, Any]:
        """Run domain analysis for a single protein chain using model-based approach

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            output_dir: Output directory
            reference: Reference version
            blast_only: Whether to use blast_only summaries

        Returns:
            Dictionary with analysis results including models and file paths
        """
        from ecod.utils.path_utils import get_all_evidence_paths
        from ecod.models.pipeline import DomainSummaryModel, ProteinResult
        from ecod.models.domain import Domain, DomainRange

        # Initialize structured result
        result = {
            "status": "failed",
            "messages": [],
            "files": {},
            "models": {},
            "stats": {}
        }

        self.logger.info(f"Starting domain analysis for {pdb_id}_{chain_id}")

        # Get all standardized paths for files
        paths = get_all_evidence_paths(output_dir, pdb_id, chain_id, reference)

        # Step 1: Create domain summary
        try:
            # Call summary module with model-based approach
            summary_result = self.summary.create_summary(
                pdb_id, chain_id, reference, output_dir, blast_only
            )

            # Handle both dictionary and legacy string returns
            if isinstance(summary_result, dict):
                summary_file = summary_result.get('file_path', '')
                summary_model = summary_result.get('summary', None)
                stats = summary_result.get('stats', {})

                if summary_file:
                    result["files"]["summary"] = summary_file
                    result["messages"].append(f"Created domain summary for {pdb_id}_{chain_id}")

                    # Add the model if available
                    if summary_model and isinstance(summary_model, DomainSummaryModel):
                        result["models"]["summary"] = summary_model
                        # Update summary model with file path if not set
                        if not summary_model.output_file_path:
                            summary_model.output_file_path = summary_file

                    # Add detailed stats
                    result["stats"].update(stats)

            else:
                # Legacy string return (backward compatibility)
                summary_file = summary_result
                if summary_file:
                    result["files"]["summary"] = summary_file
                    result["messages"].append(f"Created domain summary for {pdb_id}_{chain_id}")
                else:
                    error_msg = f"Failed to create domain summary for {pdb_id}_{chain_id}"
                    self.logger.error(error_msg)
                    result["messages"].append(error_msg)
                    return result

        except Exception as e:
            error_msg = f"Error creating domain summary for {pdb_id}_{chain_id}: {str(e)}"
            self.logger.error(error_msg)
            self.logger.exception("Exception details:")
            result["messages"].append(error_msg)
            return result

        # Step 2: Partition domains
        partition_result = self.partition.process_protein_domains(
            pdb_id=pdb_id,
            chain_id=chain_id,
            domain_summary_path=summary_file,
            output_dir=output_dir,
            reference=reference
        )

        except Exception as e:
            error_msg = f"Error creating domain partition for {pdb_id}_{chain_id}: {str(e)}"
            self.logger.error(error_msg)
            self.logger.exception("Exception details:")
            result["messages"].append(error_msg)
            # Continue with partial success if we at least have the summary

        # Set final status
        if "summary" in result["files"] and "partition" in result["files"]:
            result["status"] = "success"
        elif "summary" in result["files"]:
            result["status"] = "partial"

        # Create a ProteinResult model for standardized result format
        protein_result = ProteinResult(
            protein_id=0,  # Not available in this context
            pdb_id=pdb_id,
            chain_id=chain_id,
            process_id=0,  # Not available in this context
            success=(result["status"] == "success"),
            summary_success=("summary" in result["files"]),
            partition_success=("partition" in result["files"]),
            error=None if result["status"] != "failed" else result["messages"][-1],
            summary_file=result["files"].get("summary"),
            domain_file=result["files"].get("partition"),
            chain_blast_count=result["stats"].get("chain_blast_hits", 0),
            domain_blast_count=result["stats"].get("domain_blast_hits", 0),
            hhsearch_count=result["stats"].get("hhsearch_hits", 0)
        )

        # Add the protein result model
        result["models"]["protein_result"] = protein_result

        self.logger.info(f"Domain analysis complete for {pdb_id}_{chain_id} with status: {result['status']}")
        return result
    
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

    def reset_failed_processes(self, batch_id: int, stage: str = 'domain_partition_failed') -> int:
        """
        Reset processes that failed at a specific stage
        
        Args:
            batch_id: Batch ID
            stage: Failed stage to reset (default: 'domain_partition_failed')
            
        Returns:
            Number of reset processes
        """
        self.logger.info(f"Resetting failed processes at stage '{stage}' for batch {batch_id}")
        
        # Get database connection
        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)
        
        # Reset status for failed processes
        update_query = """
        UPDATE ecod_schema.process_status
        SET status = 'processing', 
            current_stage = 'initialized',
            error_message = NULL
        WHERE batch_id = %s
          AND current_stage = %s
          AND status = 'error'
        RETURNING id
        """
        
        try:
            # Execute update query and get affected rows
            result = db.execute_query(update_query, (batch_id, stage))
            reset_count = len(result) if result else 0
            
            self.logger.info(f"Reset {reset_count} failed processes for batch {batch_id}")
            return reset_count
            
        except Exception as e:
            self.logger.error(f"Error resetting failed processes: {e}")
            return 0
