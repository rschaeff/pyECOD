#!/usr/bin/env python3
"""
BlastPipeline module for running BLAST searches on protein chains
with standardized path handling
"""
import os
import logging
from typing import List, Dict, Any, Optional, Tuple

from ecod.core.base_pipeline import BasePipeline
from ecod.utils.path_utils import get_standardized_paths, get_file_db_path
from ecod.exceptions import PipelineError, ConfigurationError
from ecod.models import Protein, Batch

class BlastPipeline(BasePipeline):
    """Pipeline for running BLAST searches on protein chains"""

    def __init__(self, context=None):
        """
        Initialize with application context

        Args:
            context: Application context
        """
        super().__init__(context, logger_name="ecod.blast_pipeline")

        # Validate configuration after loading
        self._validate_config()

    def _load_configuration(self) -> None:
        """Load BLAST-specific configuration"""
        # Get tool paths from config
        tools = self.config.get('tools', {})
        self.blast_path = tools.get('blast_path', 'blastp')

        # Get reference paths
        ref = self.config.get('reference', {})
        self.chain_db = ref.get('chain_db')
        self.domain_db = ref.get('domain_db')

        # Default reference version
        self.default_ref_version = self.config.get('reference', {}).get('current_version', 'develop291')

    def _validate_config(self) -> None:
        """Validate configuration for the BLAST pipeline"""
        if not self.blast_path:
            error_msg = "BLAST path not configured"
            self.logger.error(error_msg)
            raise ConfigurationError(error_msg)

        if not self.chain_db:
            self.logger.warning("Chain BLAST database not configured")

        if not self.domain_db:
            self.logger.warning("Domain BLAST database not configured")

    def _get_batch_path(self, batch_id: int) -> Optional[str]:
        """Get the base path for a batch"""
        query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
        rows = self.db.execute_query(query, (batch_id,))
        if rows:
            return rows[0][0]
        return None

    def _get_batch_version(self, batch_id: int) -> str:
        """Get the reference version for a batch"""
        query = "SELECT ref_version FROM ecod_schema.batch WHERE id = %s"
        rows = self.db.execute_query(query, (batch_id,))
        if rows:
            return rows[0][0]
        return self.default_ref_version

    def get_unclassified_chains(self, limit: Optional[int] = None) -> List[Protein]:
        """Get unclassified protein chains from the database"""
        query = """
            SELECT
                p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.sequence
            FROM
                ecod_schema.protein p
            JOIN
                ecod_schema.protein_sequence ps ON p.id = ps.protein_id
            LEFT JOIN
                pdb_analysis.domain d ON p.id = d.protein_id
            LEFT JOIN
                ecod_schema.process_status ps_stat ON p.id = ps_stat.protein_id
            WHERE
                d.id IS NULL
                AND ps.sequence IS NOT NULL
                AND p.length > 30
                AND ps_stat.id IS NULL
            ORDER BY
                p.length DESC
        """

        if limit:
            query += f" LIMIT {limit}"

        rows = self.db.execute_dict_query(query)
        chains = [
            Protein(
                id=row['id'],
                pdb_id=row['pdb_id'],
                chain_id=row['chain_id'],
                source_id=row['source_id'],
                length=row['length'],
                sequence=row['sequence']
            )
            for row in rows
        ]

        self.logger.info(f"Retrieved {len(chains)} unclassified chains")
        return chains

    def create_batch(self, chains: List[Protein], batch_size: int = 100) -> Batch:
        """Create a new batch for BLAST processing"""
        # Generate batch name with timestamp
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M")
        batch_name = f"blast_batch_{timestamp}"

        # Create batch directory
        base_dir = self.config.get('paths', {}).get('output_dir', '/data/ecod/weekly_updates/batches')
        batch_path = os.path.join(base_dir, batch_name)
        os.makedirs(batch_path, exist_ok=True)

        # Create batch record in database
        batch_id = self.db.insert(
            "ecod_schema.batch",
            {
                "batch_name": batch_name,
                "base_path": batch_path,
                "type": "blast",
                "ref_version": self.default_ref_version,
                "total_items": len(chains),
                "status": "created"
            },
            "id"
        )

        # Create batch object
        batch = Batch(
            id=batch_id,
            batch_name=batch_name,
            base_path=batch_path,
            type="blast",
            ref_version=self.default_ref_version,
            total_items=len(chains),
            status="created"
        )

        # Register chains in this batch
        self._register_chains_in_batch(batch, chains)

        self.logger.info(f"Created BLAST batch: {batch_name} with {len(chains)} chains")
        return batch

    def _register_chains_in_batch(self, batch: Batch, chains: List[Protein]) -> None:
        """Register chains in a batch and create directory structure"""
        # Process each chain using standardized paths
        for chain in chains:
            # Get standardized paths
            paths = get_standardized_paths(
                batch_path=batch.base_path,
                pdb_id=chain.pdb_id,
                chain_id=chain.chain_id,
                ref_version=batch.ref_version,
                create_dirs=True
            )

            # Determine the relative path for database
            fasta_rel_path = get_file_db_path(batch.base_path, paths['fasta'])
            rel_dir = os.path.dirname(fasta_rel_path)

            # Register in process_status
            process_id = self.db.insert(
                "ecod_schema.process_status",
                {
                    "protein_id": chain.id,
                    "batch_id": batch.id,
                    "current_stage": "fasta_generated",
                    "status": "pending",
                    "relative_path": rel_dir
                },
                "id"
            )

            # Generate FASTA file
            with open(paths['fasta'], 'w') as f:
                f.write(f">{chain.source_id}\n{chain.sequence}\n")

            # Register FASTA file
            self.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "fasta",
                    "file_path": fasta_rel_path,
                    "file_exists": True
                }
            )

    def generate_fasta_files(self, batch_id: int) -> List[Tuple[int, str]]:
        """Get FASTA files for a batch that don't already have results"""
        # Get batch information
        batch_path = self._get_batch_path(batch_id)
        if not batch_path:
            self.logger.error(f"Batch {batch_id} not found")
            return []

        ref_version = self._get_batch_version(batch_id)

        # Find proteins that need processing
        query = """
            SELECT
                ps.id as process_id,
                p.id as protein_id,
                p.pdb_id,
                p.chain_id,
                p.source_id,
                p_seq.sequence
            FROM
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.protein_sequence p_seq ON p.id = p_seq.protein_id
            LEFT JOIN (
                SELECT DISTINCT process_id
                FROM ecod_schema.process_file
                WHERE file_type IN ('chain_blast_result', 'domain_blast_result')
                AND file_exists = TRUE
            ) pf ON ps.id = pf.process_id
            WHERE
                ps.batch_id = %s
                AND pf.process_id IS NULL
        """

        rows = self.db.execute_dict_query(query, (batch_id,))

        fasta_paths = []
        for row in rows:
            # Get standardized paths
            paths = get_standardized_paths(
                batch_path=batch_path,
                pdb_id=row['pdb_id'],
                chain_id=row['chain_id'],
                ref_version=ref_version,
                create_dirs=True
            )

            # Check if FASTA already exists
            if not os.path.exists(paths['fasta']):
                with open(paths['fasta'], 'w') as f:
                    f.write(f">{row['source_id']}\n{row['sequence']}\n")

            fasta_paths.append((row['process_id'], paths['fasta']))

        self.logger.info(f"Generated/found {len(fasta_paths)} FASTA files for batch {batch_id}")
        return fasta_paths

    def run_chain_blast(self, batch_id: int, batch_size: int = 100) -> List[str]:
        """Run chain-wise BLAST with improved error handling"""
        try:
            # Get FASTA files for this batch
            fasta_paths = self.generate_fasta_files(batch_id)
            if not fasta_paths:
                self.logger.warning(f"No FASTA files found for batch {batch_id}")
                return []

            # Get batch information
            batch_path = self._get_batch_path(batch_id)
            if not batch_path:
                error_msg = f"Batch {batch_id} not found"
                self.logger.error(error_msg)
                raise PipelineError(error_msg)

            ref_version = self._get_batch_version(batch_id)

            # Check for chain BLAST database
            if not self.chain_db:
                error_msg = "Chain BLAST database not configured"
                self.logger.error(error_msg)
                raise ConfigurationError(error_msg)

            # Create job template with standardized paths
            job_items = []
            for process_id, fasta_path in fasta_paths:
                # Extract PDB and chain IDs from path
                pdb_id, chain_id = self._extract_ids_from_path(fasta_path)

                if not pdb_id or not chain_id:
                    self.logger.warning(f"Could not extract PDB/chain IDs from {fasta_path}")
                    continue

                # Get standardized paths
                paths = get_standardized_paths(
                    batch_path=batch_path,
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    ref_version=ref_version,
                    create_dirs=True
                )

                # Create command
                cmd = f"{self.blast_path} -query {fasta_path} -db {self.chain_db} " \
                      f"-outfmt 5 -num_alignments 5000 -evalue 0.002 " \
                      f"-out {paths['chain_blast']}"

                job_items.append((process_id, fasta_path, cmd))

            # Group into batches
            batches = [job_items[i:i+batch_size] for i in range(0, len(job_items), batch_size)]

            # Create and submit jobs
            job_ids = []
            for batch_num, batch_items in enumerate(batches):
                # Create job script
                job_name = f"chain_blast_{batch_id}_{batch_num}"
                script_dir = os.path.join(batch_path, "scripts")
                os.makedirs(script_dir, exist_ok=True)
                script_path = os.path.join(script_dir, f"{job_name}.sh")

                # Create commands
                commands = [item[2] for item in batch_items]

                # Create script
                self.job_manager.create_job_script(
                    commands,
                    job_name,
                    script_dir,
                    threads=4,
                    memory="4G",
                    time="24:00:00"
                )

                # Submit job
                slurm_job_id = self.job_manager.submit_job(script_path)
                if not slurm_job_id:
                    self.logger.warning(f"Failed to submit job for {job_name}")
                    continue

                # Record job in database
                job_db_id = self.db.insert(
                    "ecod_schema.job",
                    {
                        "batch_id": batch_id,
                        "job_type": "chain_blast",
                        "slurm_job_id": slurm_job_id,
                        "items_count": len(batch_items)
                    },
                    "id"
                )

                # Record items in this job
                for process_id, _, _ in batch_items:
                    self.db.insert(
                        "ecod_schema.job_item",
                        {
                            "job_id": job_db_id,
                            "process_id": process_id
                        }
                    )

                    # Update process status
                    self.db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "chain_blast_running",
                            "status": "processing"
                        },
                        "id = %s",
                        (process_id,)
                    )

                job_ids.append(slurm_job_id)
                self.logger.info(f"Submitted job {slurm_job_id} for {job_name}")

            self.logger.info(f"Submitted {len(job_ids)} chain BLAST jobs for batch {batch_id}")
            return job_ids

        except ConfigurationError:
            # Re-raise configuration errors
            raise
        except Exception as e:
            error_msg = f"Unexpected error in chain BLAST pipeline: {str(e)}"
            self.logger.error(error_msg, exc_info=True)
            raise PipelineError(error_msg) from e

    def run_domain_blast(self, batch_id: int, batch_size: int = 100) -> List[str]:
        """Run domain-wise BLAST on a batch of sequences"""
        try:
            # Get FASTA files for this batch
            fasta_paths = self.generate_fasta_files(batch_id)
            if not fasta_paths:
                self.logger.warning(f"No FASTA files found for batch {batch_id}")
                return []

            # Get batch information
            batch_path = self._get_batch_path(batch_id)
            if not batch_path:
                error_msg = f"Batch {batch_id} not found"
                self.logger.error(error_msg)
                raise PipelineError(error_msg)

            ref_version = self._get_batch_version(batch_id)

            # Get domain BLAST database path
            if not self.domain_db:
                error_msg = "Domain BLAST database not configured"
                self.logger.error(error_msg)
                raise ConfigurationError(error_msg)

            # Create job template with standardized paths
            job_items = []
            for process_id, fasta_path in fasta_paths:
                # Extract PDB and chain IDs from path
                pdb_id, chain_id = self._extract_ids_from_path(fasta_path)

                if not pdb_id or not chain_id:
                    self.logger.warning(f"Could not extract PDB/chain IDs from {fasta_path}")
                    continue

                # Get standardized paths
                paths = get_standardized_paths(
                    batch_path=batch_path,
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    ref_version=ref_version,
                    create_dirs=True
                )

                # Create command
                cmd = f"{self.blast_path} -query {fasta_path} -db {self.domain_db} " \
                      f"-outfmt 5 -num_alignments 5000 -evalue 0.002 " \
                      f"-out {paths['domain_blast']}"

                job_items.append((process_id, fasta_path, cmd))

            # Group into batches
            batches = [job_items[i:i+batch_size] for i in range(0, len(job_items), batch_size)]

            # Create and submit jobs
            job_ids = []
            for batch_num, batch_items in enumerate(batches):
                # Create job script
                job_name = f"domain_blast_{batch_id}_{batch_num}"
                script_dir = os.path.join(batch_path, "scripts")
                os.makedirs(script_dir, exist_ok=True)
                script_path = os.path.join(script_dir, f"{job_name}.sh")

                # Create commands
                commands = [item[2] for item in batch_items]

                # Create script
                self.job_manager.create_job_script(
                    commands,
                    job_name,
                    script_dir,
                    threads=4,
                    memory="4G",
                    time="24:00:00"
                )

                # Submit job
                slurm_job_id = self.job_manager.submit_job(script_path)
                if not slurm_job_id:
                    self.logger.warning(f"Failed to submit job for {job_name}")
                    continue

                # Record job in database
                job_db_id = self.db.insert(
                    "ecod_schema.job",
                    {
                        "batch_id": batch_id,
                        "job_type": "domain_blast",
                        "slurm_job_id": slurm_job_id,
                        "items_count": len(batch_items)
                    },
                    "id"
                )

                # Record items in this job
                for process_id, _, _ in batch_items:
                    self.db.insert(
                        "ecod_schema.job_item",
                        {
                            "job_id": job_db_id,
                            "process_id": process_id
                        }
                    )

                    # Update process status
                    self.db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "domain_blast_running",
                            "status": "processing"
                        },
                        "id = %s",
                        (process_id,)
                    )

                job_ids.append(slurm_job_id)
                self.logger.info(f"Submitted job {slurm_job_id} for {job_name}")

            self.logger.info(f"Submitted {len(job_ids)} domain BLAST jobs for batch {batch_id}")
            return job_ids

        except ConfigurationError:
            # Re-raise configuration errors
            raise
        except Exception as e:
            error_msg = f"Unexpected error in domain BLAST pipeline: {str(e)}"
            self.logger.error(error_msg, exc_info=True)
            raise PipelineError(error_msg) from e

    def check_job_status(self, batch_id: Optional[int] = None) -> None:
        """Check status of submitted BLAST jobs and update database"""
        # Query to get submitted jobs
        query = """
            SELECT
                j.id, j.slurm_job_id, j.job_type, j.batch_id
            FROM
                ecod_schema.job j
            WHERE
                j.status = 'submitted'
        """

        if batch_id:
            query += " AND j.batch_id = %s"
            rows = self.db.execute_dict_query(query, (batch_id,))
        else:
            rows = self.db.execute_dict_query(query)

        self.logger.info(f"Checking status of {len(rows)} submitted jobs")

        for row in rows:
            job_id = row['id']
            slurm_job_id = row['slurm_job_id']
            job_type = row['job_type']

            # Check job status with Slurm
            status = self.job_manager.check_job_status(slurm_job_id)

            if status == "COMPLETED":
                # Update job status
                self.db.update(
                    "ecod_schema.job",
                    {
                        "status": "completed"
                    },
                    "id = %s",
                    (job_id,)
                )

                # Update item statuses and check for output files
                self._update_job_items(job_id, job_type)

            elif status in ["FAILED", "TIMEOUT", "CANCELLED"]:
                self.db.update(
                    "ecod_schema.job",
                    {
                        "status": "failed"
                    },
                    "id = %s",
                    (job_id,)
                )

    def _update_job_items(self, job_id: int, job_type: str) -> None:
        """Update status of items in a completed job"""
        # Get items in this job with PDB/chain information
        query = """
            SELECT
                ji.id, ji.process_id, b.base_path, b.ref_version,
                p.pdb_id, p.chain_id
            FROM
                ecod_schema.job_item ji
            JOIN
                ecod_schema.process_status ps ON ji.process_id = ps.id
            JOIN
                ecod_schema.batch b ON ps.batch_id = b.id
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            WHERE
                ji.job_id = %s
        """

        rows = self.db.execute_dict_query(query, (job_id,))

        for row in rows:
            process_id = row['process_id']
            batch_path = row['base_path']
            ref_version = row['ref_version']
            pdb_id = row['pdb_id']
            chain_id = row['chain_id']

            # Get standardized paths
            paths = get_standardized_paths(
                batch_path=batch_path,
                pdb_id=pdb_id,
                chain_id=chain_id,
                ref_version=ref_version,
                create_dirs=False
            )

            # Determine expected output file based on job type
            if job_type == "chain_blast":
                output_file = paths['chain_blast']
                file_type = "chain_blast_result"
            elif job_type == "domain_blast":
                output_file = paths['domain_blast']
                file_type = "domain_blast_result"
            else:
                self.logger.warning(f"Unknown job type: {job_type}")
                continue

            # Check if output file exists
            file_exists = os.path.exists(output_file)
            file_size = os.path.getsize(output_file) if file_exists else 0

            # Get relative path for database
            relative_output = get_file_db_path(batch_path, output_file)

            # Register file
            self.register_process_file(
                process_id,
                file_type,
                relative_output,
                file_exists,
                file_size
            )

            # Determine next stage and status
            next_stage = "blast_completed" if file_exists and file_size > 0 else "blast_failed"
            next_status = "success" if file_exists and file_size > 0 else "error"
            error_message = None if file_exists and file_size > 0 else f"BLAST output file missing or empty: {output_file}"

            # Update process status
            update_data = {
                "current_stage": next_stage,
                "status": next_status
            }

            if error_message:
                update_data["error_message"] = error_message

            self.db.update(
                "ecod_schema.process_status",
                update_data,
                "id = %s",
                (process_id,)
            )

    def register_process_file(self, process_id: int, file_type: str,
                            file_path: str, file_exists: bool = True,
                            file_size: Optional[int] = None) -> None:
        """Register a process file with proper duplicate handling"""
        try:
            # Check if file_size is provided
            if file_size is None and file_exists:
                try:
                    if os.path.exists(file_path):
                        file_size = os.path.getsize(file_path)
                    else:
                        file_size = 0
                except Exception:
                    file_size = 0

            # Check if record already exists
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """

            existing = self.db.execute_query(query, (process_id, file_type))

            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": file_path,
                        "file_exists": file_exists,
                        "file_size": file_size or 0
                    },
                    "id = %s",
                    (existing[0][0],)
                )
                self.logger.debug(f"Updated existing {file_type} file record for process {process_id}")
            else:
                # Insert new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": file_path,
                        "file_exists": file_exists,
                        "file_size": file_size or 0
                    }
                )
                self.logger.debug(f"Created new {file_type} file record for process {process_id}")
        except Exception as e:
            self.logger.warning(f"Error registering {file_type} file for process {process_id}: {e}")

    def _extract_ids_from_path(self, file_path: str) -> Tuple[Optional[str], Optional[str]]:
        """Extract PDB ID and chain ID from file path"""
        import os
        import re

        # Extract basename
        basename = os.path.basename(file_path)

        # Try to match [pdbid]_[chainid].fa pattern
        match = re.match(r'([a-zA-Z0-9]{4})_([a-zA-Z0-9])\.fa', basename)
        if match:
            return match.group(1), match.group(2)

        # Alternative formats
        match = re.match(r'([a-zA-Z0-9]{4})_([a-zA-Z0-9])', basename)
        if match:
            return match.group(1), match.group(2)

        return None, None

    def _count_existing_results(self, batch_id: int, file_type: str) -> int:
        """Count the number of existing result files for a batch"""
        query = """
            SELECT COUNT(DISTINCT ps.id)
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE ps.batch_id = %s
            AND pf.file_type = %s
            AND pf.file_exists = TRUE
        """

        try:
            rows = self.db.execute_query(query, (batch_id, file_type))
            if rows and rows[0][0]:
                return rows[0][0]
        except Exception as e:
            self.logger.error(f"Error counting existing {file_type} files: {e}")

        return 0

    def gather_results(self, batch_id: int) -> Dict[str, Dict[str, str]]:
        """Gather and organize BLAST results for a batch"""
        query = """
            SELECT
                ps.id, p.source_id, pf.file_path, pf.file_type
            FROM
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE
                ps.batch_id = %s
                AND ps.status = 'success'
                AND pf.file_type IN ('chain_blast_result', 'domain_blast_result')
                AND pf.file_exists = TRUE
        """

        rows = self.db.execute_dict_query(query, (batch_id,))

        results = {
            'chain_results': {},
            'domain_results': {}
        }

        batch_path = self._get_batch_path(batch_id)
        if not batch_path:
            self.logger.error(f"Batch {batch_id} not found")
            return results

        for row in rows:
            source_id = row['source_id']
            file_path = row['file_path']
            file_type = row['file_type']

            # Resolve relative path to absolute path
            abs_path = os.path.join(batch_path, file_path)

            if file_type == 'chain_blast_result':
                results['chain_results'][source_id] = abs_path
            elif file_type == 'domain_blast_result':
                results['domain_results'][source_id] = abs_path

        self.logger.info(f"Found {len(results['chain_results'])} chain results and {len(results['domain_results'])} domain results")
        return results
