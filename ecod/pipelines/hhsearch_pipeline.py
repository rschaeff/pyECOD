#!/usr/bin/env python3
"""
HHSearch Pipeline for protein domain analysis with standardized path handling
"""
import os
import logging
import re
from typing import List, Dict, Any, Optional, Tuple, Set
from datetime import datetime
import xml.etree.ElementTree as ET

from ecod.core.base_pipeline import BasePipeline
from ecod.utils.path_utils import get_standardized_paths, get_file_db_path
from ecod.exceptions import PipelineError, ConfigurationError
from ecod.core.context import ApplicationContext

class HHSearchPipeline(BasePipeline):
    """Pipeline for running HHSearch for domain classification"""

    def __init__(self, context=None):
        """
        Initialize with application context

        Args:
            context: Application context with shared resources
        """
        super().__init__(context, logger_name="ecod.hhsearch_pipeline")

        # Validate configuration after loading
        self._validate_config()

    def _load_configuration(self) -> None:
        """Load HHSearch-specific configuration"""
        # Get tool paths from config
        tools = self.config.get('tools', {})
        self.hhblits_path = tools.get('hhblits_path', 'hhblits')
        self.hhsearch_path = tools.get('hhsearch_path', 'hhsearch')
        self.hhmake_path = tools.get('hhmake_path', 'hhmake')

        # Get reference paths
        ref = self.config.get('reference', {})
        self.uniclust_db = ref.get('uniclust_db')
        self.ecod_ref_db = ref.get('ecod_hh_db')

        # Default reference version
        self.default_ref_version = self.config.get('reference', {}).get('current_version', 'develop291')

    def _validate_config(self) -> None:
        """Validate configuration for the HHSearch pipeline"""
        if not self.hhblits_path:
            self.logger.warning("HHblits path not configured, using default")

        if not self.uniclust_db:
            error_msg = "Uniclust database not configured"
            self.logger.error(error_msg)
            raise ConfigurationError(error_msg)

        if not self.ecod_ref_db:
            error_msg = "ECOD HHSearch database not configured"
            self.logger.error(error_msg)
            raise ConfigurationError(error_msg)

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

    def create_batch(self, chains: List[Dict[str, Any]], batch_size: int = 500) -> int:
        """Create a new batch for HHsearch processing"""
        # Generate batch name with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M")
        batch_name = f"hhsearch_batch_{timestamp}"

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
                "type": "hhsearch",
                "ref_version": self.default_ref_version,
                "total_items": len(chains),
                "status": "created"
            },
            "id"
        )

        # Register chains in this batch
        self._register_chains_in_batch(batch_path, batch_id, chains)

        self.logger.info(f"Created HHsearch batch: {batch_name} with {len(chains)} chains")
        return batch_id

    def _register_chains_in_batch(self, batch_path: str, batch_id: int, chains: List[Dict[str, Any]]) -> None:
        """Register chains in a batch using standardized paths"""
        # Process each chain
        for chain in chains:
            pdb_id = chain['pdb_id']
            chain_id = chain['chain_id']
            protein_id = chain['id']
            ref_version = self.default_ref_version

            # Get standardized paths
            paths = get_standardized_paths(
                batch_path=batch_path,
                pdb_id=pdb_id,
                chain_id=chain_id,
                ref_version=ref_version,
                create_dirs=True
            )

            # Determine the relative path for database storage
            fasta_rel_path = get_file_db_path(batch_path, paths['fasta'])
            rel_dir = os.path.dirname(fasta_rel_path)

            # Register in process_status
            process_id = self.db.insert(
                "ecod_schema.process_status",
                {
                    "protein_id": protein_id,
                    "batch_id": batch_id,
                    "current_stage": "fasta_generated",
                    "status": "pending",
                    "relative_path": rel_dir
                },
                "id"
            )

            # Create FASTA file
            with open(paths['fasta'], 'w') as f:
                f.write(f">{pdb_id}_{chain_id}\n{chain['sequence']}\n")

            # Register FASTA file in database
            self.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "fasta",
                    "file_path": fasta_rel_path,
                    "file_exists": True
                }
            )

    def get_chains_for_profile_generation(self, batch_id: int) -> List[Dict[str, Any]]:
        """Get chains ready for profile generation"""
        query = """
            SELECT
                ps.id, p.id as protein_id, p.pdb_id, p.chain_id,
                p_seq.sequence,
                ps.relative_path, b.base_path
            FROM
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.protein_sequence p_seq ON p.id = p_seq.protein_id
            JOIN
                ecod_schema.batch b ON ps.batch_id = b.id
            WHERE
                ps.batch_id = %s
                AND ps.current_stage != 'created'
                AND ps.status != 'failed'
                AND NOT EXISTS (
                    SELECT 1 FROM ecod_schema.process_file
                    WHERE process_id = ps.id AND file_type = 'a3m'
                )
        """

        rows = self.db.execute_dict_query(query, (batch_id,))
        return rows

    def generate_profiles(self, batch_id: int, threads: int = 8, memory: str = "16G", force: bool = False) -> List[str]:
        """Generate HHblits profiles for a batch"""
        chains = self.get_chains_for_profile_generation(batch_id)

        if not chains:
            self.logger.warning(f"No chains ready for profile generation in batch {batch_id}")
            return []

        job_ids = []
        for chain in chains:
            job_id = self._submit_hhblits_job(
                chain['id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_id,
                threads,
                memory,
                force
            )

            if job_id:
                job_ids.append(job_id)

        self.logger.info(f"Submitted {len(job_ids)} HHblits profile generation jobs")
        return job_ids

    def _submit_hhblits_job(self, process_id: int, pdb_id: str, chain_id: str,
                           batch_id: int, threads: int, memory: str,
                           force: bool = False) -> Optional[str]:
        """Submit a job to generate HHblits profile for a chain"""
        # Get batch information
        batch_path = self._get_batch_path(batch_id)
        if not batch_path:
            self.logger.error(f"Batch {batch_id} not found")
            return None

        ref_version = self._get_batch_version(batch_id)

        # Get standardized paths
        paths = get_standardized_paths(
            batch_path=batch_path,
            pdb_id=pdb_id,
            chain_id=chain_id,
            ref_version=ref_version,
            create_dirs=True
        )

        # Create a scripts directory for job scripts
        scripts_dir = os.path.join(batch_path, "scripts")
        os.makedirs(scripts_dir, exist_ok=True)

        # Define expected input/output files
        fa_file = paths['fasta']
        a3m_file = paths['a3m']
        hhm_file = paths['hhm']

        # Check if both a3m and hhm files already exist and we're not forcing a rerun
        if (os.path.exists(a3m_file) and os.path.getsize(a3m_file) > 0 and
            os.path.exists(hhm_file) and os.path.getsize(hhm_file) > 0 and
            not force):
            self.logger.info(f"Skipping HHblits for {pdb_id}_{chain_id}: a3m and hhm profiles already exist")

            # Also check if the profile is registered in the database
            check_query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = 'a3m'
            """
            existing_profile = self.db.execute_query(check_query, (process_id,))

            if not existing_profile:
                # Register the existing a3m file
                a3m_rel_path = get_file_db_path(batch_path, a3m_file)
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "a3m",
                        "file_path": a3m_rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(a3m_file)
                    }
                )

                # Register the existing hhm file
                hhm_rel_path = get_file_db_path(batch_path, hhm_file)
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "hhm",
                        "file_path": hhm_rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(hhm_file)
                    }
                )

                # Update process status
                self.db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": "profile_complete",
                        "status": "success"
                    },
                    "id = %s",
                    (process_id,)
                )

            # Return None to indicate no job was submitted
            return None

        # Create FASTA file if it doesn't exist
        if not os.path.exists(fa_file):
            # Get sequence from database
            seq_query = """
            SELECT ps.sequence
            FROM ecod_schema.protein p
            JOIN ecod_schema.protein_sequence ps ON p.id = ps.protein_id
            JOIN ecod_schema.process_status proc ON p.id = proc.protein_id
            WHERE proc.id = %s
            """
            seq_result = self.db.execute_query(seq_query, (process_id,))

            if seq_result and seq_result[0][0]:
                # Create directory if it doesn't exist
                os.makedirs(os.path.dirname(fa_file), exist_ok=True)

                # Write FASTA file
                with open(fa_file, 'w') as f:
                    f.write(f">{pdb_id}_{chain_id}\n{seq_result[0][0]}\n")

                # Register in database
                fasta_rel_path = get_file_db_path(batch_path, fa_file)

                # Check if FASTA file is already registered
                check_query = """
                SELECT id FROM ecod_schema.process_file
                WHERE process_id = %s AND file_type = 'fasta'
                """
                existing_file = self.db.execute_query(check_query, (process_id,))

                if not existing_file:
                    # Register file in database if not already registered
                    self.db.insert(
                        "ecod_schema.process_file",
                        {
                            "process_id": process_id,
                            "file_type": "fasta",
                            "file_path": fasta_rel_path,
                            "file_exists": True,
                            "file_size": os.path.getsize(fa_file)
                        }
                    )
            else:
                self.logger.warning(f"No sequence found for process_id {process_id}")
                return None

        # Create job script with clear command separation
        job_name = f"hhblits_{pdb_id}_{chain_id}"
        job_script = os.path.join(scripts_dir, f"{job_name}.sh")

        commands = [
            f"module purge",
            f"module load hh-suite",
            f"{self.hhblits_path} -i {fa_file} -oa3m {a3m_file} -ohhm {hhm_file} "
            f"-d {self.uniclust_db} "
            f"-n 3 -maxfilt 20000 -diff inf -id 99 -cov 50 "
            f"-cpu {threads}"
        ]

        # Create job script
        self.job_manager.create_job_script(
            commands,
            job_name,
            scripts_dir,
            threads=threads,
            memory=memory,
            time="24:00:00"
        )

        # Submit job
        slurm_job_id = self.job_manager.submit_job(job_script)
        if not slurm_job_id:
            return None

        # Record job in database
        job_db_id = self.db.insert(
            "ecod_schema.job",
            {
                "batch_id": batch_id,
                "job_type": "hhblits",
                "slurm_job_id": slurm_job_id,
                "items_count": 1,
                "status": "submitted"
            },
            "id"
        )

        # Record chain in this job
        self.db.insert(
            "ecod_schema.job_item",
            {
                "job_id": job_db_id,
                "process_id": process_id,
                "status": "pending"
            }
        )

        # Update process status
        self.db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "profile_running",
                "status": "processing"
            },
            "id = %s",
            (process_id,)
        )

        self.logger.info(f"Submitted HHblits job {slurm_job_id} for {pdb_id}_{chain_id}")
        return slurm_job_id

    def run_hhsearch(self, batch_id: int, threads: int = 8, memory: str = "16G") -> List[str]:
        """Run HHsearch for chains with completed profiles"""
        # Get chains with completed profiles
        chains = self._get_chains_with_profiles(batch_id)
        if not chains:
            self.logger.warning(f"No chains with completed profiles found in batch {batch_id}")
            return []

        job_ids = []
        for chain in chains:
            job_id = self._submit_hhsearch_job(
                chain['id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_id,
                threads,
                memory
            )

            if job_id:
                job_ids.append(job_id)

        self.logger.info(f"Submitted {len(job_ids)} HHsearch jobs")
        return job_ids

    def _get_chains_with_profiles(self, batch_id: int) -> List[Dict[str, Any]]:
        """Get chains with completed profiles"""
        query = """
            SELECT
                ps.id, p.id as protein_id, p.pdb_id, p.chain_id,
                ps.relative_path, b.base_path
            FROM
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.batch b ON ps.batch_id = b.id
            JOIN
                ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE
                ps.batch_id = %s
                AND ps.current_stage = 'profile_complete'
                AND ps.status = 'success'
                AND pf.file_type = 'a3m'
                AND pf.file_exists = TRUE
                AND NOT EXISTS (
                    SELECT 1 FROM ecod_schema.process_file
                    WHERE process_id = ps.id AND file_type = 'hhr'
                )
        """

        rows = self.db.execute_dict_query(query, (batch_id,))
        return rows

    def _submit_hhsearch_job(self, process_id: int, pdb_id: str, chain_id: str,
                           batch_id: int, threads: int, memory: str) -> Optional[str]:
        """Submit a job to run HHsearch for a single chain"""
        # Get batch information
        batch_path = self._get_batch_path(batch_id)
        if not batch_path:
            self.logger.error(f"Batch {batch_id} not found")
            return None

        ref_version = self._get_batch_version(batch_id)

        # Get standardized paths
        paths = get_standardized_paths(
            batch_path=batch_path,
            pdb_id=pdb_id,
            chain_id=chain_id,
            ref_version=ref_version,
            create_dirs=True
        )

        # Create a scripts directory for job scripts
        scripts_dir = os.path.join(batch_path, "scripts")
        os.makedirs(scripts_dir, exist_ok=True)

        # Define input/output files
        a3m_file = paths['a3m']
        hhm_file = paths['hhm']
        hhr_file = paths['hhr']
        hh_xml_file = paths['hh_xml']

        # Check if A3M exists (required for HHSearch)
        if not os.path.exists(a3m_file):
            self.logger.warning(f"A3M file not found: {a3m_file}")
            return None

        # Create job script
        job_name = f"hhsearch_{pdb_id}_{chain_id}"
        job_script = os.path.join(scripts_dir, f"{job_name}.sh")

        commands = [
            f"module purge",
            f"module load hh-suite",
            # First ensure HHM file exists
            f"{self.hhmake_path} -i {a3m_file} -o {hhm_file}",
            # Then run HHsearch against ECOD database
            f"{self.hhsearch_path} "
            f"-i {hhm_file} "
            f"-d {self.ecod_ref_db} "
            f"-o {hhr_file} "
            f"-cpu {threads} "
            f"-v 2 -p 60.0 -z 100 -b 100 -ssm 2 -sc 1 -aliw 80 -glob"
        ]

        self.job_manager.create_job_script(
            commands,
            job_name,
            scripts_dir,
            threads=threads,
            memory=memory,
            time="24:00:00"
        )

        # Submit job
        slurm_job_id = self.job_manager.submit_job(job_script)
        if not slurm_job_id:
            return None

        # Record job in database
        job_db_id = self.db.insert(
            "ecod_schema.job",
            {
                "batch_id": batch_id,
                "job_type": "hhsearch",
                "slurm_job_id": slurm_job_id,
                "items_count": 1,
                "status": "submitted"
            },
            "id"
        )

        # Record chain in this job
        self.db.insert(
            "ecod_schema.job_item",
            {
                "job_id": job_db_id,
                "process_id": process_id,
                "status": "pending"
            }
        )

        # Update process status
        self.db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "hhsearch_running",
                "status": "processing"
            },
            "id = %s",
            (process_id,)
        )

        return slurm_job_id

    def check_status(self, batch_id: Optional[int] = None) -> None:
        """Check status of submitted HHSearch jobs and update database"""
        query = """
            SELECT
                j.id, j.slurm_job_id, j.job_type, j.batch_id
            FROM
                ecod_schema.job j
            WHERE
                j.status = 'submitted'
                AND j.job_type IN ('hhblits', 'hhsearch')
        """

        if batch_id:
            query += " AND j.batch_id = %s"
            rows = self.db.execute_dict_query(query, (batch_id,))
        else:
            rows = self.db.execute_dict_query(query)

        self.logger.info(f"Checking status of {len(rows)} HHSearch/HHblits jobs")

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

                # Update item statuses to mark as failed
                self._update_failed_job_items(job_id, job_type)

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

            # Determine expected output files based on job type
            if job_type == "hhblits":
                # For HHblits, check for a3m and hhm files
                a3m_file = paths['a3m']
                hhm_file = paths['hhm']

                a3m_exists = os.path.exists(a3m_file) and os.path.getsize(a3m_file) > 0
                hhm_exists = os.path.exists(hhm_file) and os.path.getsize(hhm_file) > 0

                if a3m_exists and hhm_exists:
                    # Register both files in database
                    a3m_rel_path = get_file_db_path(batch_path, a3m_file)
                    hhm_rel_path = get_file_db_path(batch_path, hhm_file)

                    self.register_file(
                        process_id,
                        "a3m",
                        a3m_rel_path,
                        True,
                        os.path.getsize(a3m_file)
                    )

                    self.register_file(
                        process_id,
                        "hhm",
                        hhm_rel_path,
                        True,
                        os.path.getsize(hhm_file)
                    )

                    # Update process status
                    self.db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "profile_complete",
                            "status": "success"
                        },
                        "id = %s",
                        (process_id,)
                    )
                else:
                    # File missing or empty
                    error_files = []
                    if not a3m_exists:
                        error_files.append(a3m_file)
                    if not hhm_exists:
                        error_files.append(hhm_file)

                    self.db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "profile_failed",
                            "status": "error",
                            "error_message": f"Output files missing or empty: {', '.join(error_files)}"
                        },
                        "id = %s",
                        (process_id,)
                    )

            elif job_type == "hhsearch":
                # For HHSearch, check for hhr file
                hhr_file = paths['hhr']

                hhr_exists = os.path.exists(hhr_file) and os.path.getsize(hhr_file) > 0

                if hhr_exists:
                    # Register HHR file in database
                    hhr_rel_path = get_file_db_path(batch_path, hhr_file)

                    self.register_file(
                        process_id,
                        "hhr",
                        hhr_rel_path,
                        True,
                        os.path.getsize(hhr_file)
                    )

                    # Process HHR file to generate XML summary
                    hh_xml_file = paths['hh_xml']
                    os.makedirs(os.path.dirname(hh_xml_file), exist_ok=True)

                    try:
                        # Parse HHR file
                        hits = self._parse_hhr_file(hhr_file)

                        # Generate XML summary
                        if hits:
                            self._generate_summary_xml(pdb_id, chain_id, hits, hh_xml_file)

                            # Register XML file
                            hh_xml_rel_path = get_file_db_path(batch_path, hh_xml_file)

                            self.register_file(
                                process_id,
                                "hh_xml",
                                hh_xml_rel_path,
                                True,
                                os.path.getsize(hh_xml_file)
                            )
                    except Exception as e:
                        self.logger.error(f"Error processing HHR file {hhr_file}: {str(e)}")

                    # Update process status
                    self.db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "hhsearch_complete",
                            "status": "success"
                        },
                        "id = %s",
                        (process_id,)
                    )
                else:
                    # File missing or empty
                    self.db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "hhsearch_failed",
                            "status": "error",
                            "error_message": f"Output file missing or empty: {hhr_file}"
                        },
                        "id = %s",
                        (process_id,)
                    )
            else:
                self.logger.warning(f"Unknown job type: {job_type}")

    def _update_failed_job_items(self, job_id: int, job_type: str) -> None:
        """Update status of items in a failed job"""
        # Get items in this job
        query = """
            SELECT ji.process_id
            FROM ecod_schema.job_item ji
            WHERE ji.job_id = %s
        """

        rows = self.db.execute_query(query, (job_id,))

        for row in rows:
            process_id = row[0]

            # Determine failure stage based on job type
            if job_type == "hhblits":
                failure_stage = "profile_failed"
            elif job_type == "hhsearch":
                failure_stage = "hhsearch_failed"
            else:
                failure_stage = "processing_failed"

            # Update process status
            self.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": failure_stage,
                    "status": "error",
                    "error_message": f"Job failed or was cancelled"
                },
                "id = %s",
                (process_id,)
            )

    def process_specific_proteins(self, batch_id: int, protein_ids: List[int],
                             threads: int = 8, memory: str = "16G") -> List[str]:
        """Run HHSearch for specific proteins

        Args:
            batch_id: Batch ID
            protein_ids: List of protein IDs to process
            threads: Number of threads for HHSearch
            memory: Memory allocation for HHSearch

        Returns:
            List of job IDs
        """
        # Get processes for specified proteins
        query = """
            SELECT
                ps.id, p.pdb_id, p.chain_id
            FROM
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            WHERE
                ps.batch_id = %s
                AND ps.protein_id IN %s
        """

        try:
            rows = self.db.execute_dict_query(query, (batch_id, tuple(protein_ids)))
        except Exception as e:
            self.logger.error(f"Error querying processes for proteins: {e}")
            return []

        # For each protein, submit profile generation followed by HHSearch
        job_ids = []

        for row in rows:
            process_id = row['id']
            pdb_id = row['pdb_id']
            chain_id = row['chain_id']

            # First generate profile
            profile_job_id = self._submit_hhblits_job(
                process_id,
                pdb_id,
                chain_id,
                batch_id,
                threads,
                memory,
                False
            )

            if profile_job_id:
                job_ids.append(profile_job_id)

        return job_ids

    def _parse_hhr_file(self, hhr_file: str) -> List[Dict[str, Any]]:
        """Parse HHsearch result file (HHR format)"""
        hits = []

        try:
            with open(hhr_file, 'r') as f:
                lines = f.readlines()

            # Process header lines (get query info)
            query_name = None
            query_length = None
            for i, line in enumerate(lines):
                if line.startswith("Query "):
                    parts = line.strip().split()
                    query_name = parts[1]
                    continue
                if line.startswith("Match_columns "):
                    parts = line.strip().split()
                    query_length = int(parts[1])
                    break

            # Find the beginning of the hit table
            table_start = None
            for i, line in enumerate(lines):
                if line.startswith(" No Hit"):
                    table_start = i + 1
                    break

            if not table_start:
                return hits

            # Process hits
            current_hit = None
            in_alignment = False
            query_ali = ""
            template_ali = ""

            i = table_start
            while i < len(lines):
                line = lines[i].strip()

                # New hit begins with a line like " 1 e4tm9c1 etc"
                if re.match(r"^\s*\d+\s+\S+", line) and not in_alignment:
                    # Store previous hit if exists
                    if current_hit and 'query_ali' in current_hit and 'template_ali' in current_hit:
                        hits.append(current_hit)

                    # Parse hit line
                    parts = line.split()
                    hit_num = int(parts[0])
                    hit_id = parts[1]

                    # Find probability, e-value, etc.
                    j = i
                    while j < len(lines) and not lines[j].startswith(">"):
                        j += 1

                    prob = None
                    e_value = None
                    score = None

                    for k in range(i, j):
                        if "Probab=" in lines[k]:
                            prob_match = re.search(r"Probab=(\d+\.\d+)", lines[k])
                            if prob_match:
                                prob = float(prob_match.group(1))

                        if "E-value=" in lines[k]:
                            eval_match = re.search(r"E-value=(\S+)", lines[k])
                            if eval_match:
                                try:
                                    e_value = float(eval_match.group(1))
                                except ValueError:
                                    e_value = 999.0

                        if "Score=" in lines[k]:
                            score_match = re.search(r"Score=(\d+\.\d+)", lines[k])
                            if score_match:
                                score = float(score_match.group(1))

                    # Create new hit
                    current_hit = {
                        'hit_num': hit_num,
                        'hit_id': hit_id,
                        'probability': prob,
                        'e_value': e_value,
                        'score': score
                    }

                    query_ali = ""
                    template_ali = ""
                    in_alignment = False

                    # Skip to alignment section
                    while i < len(lines) and not lines[i].startswith("Q "):
                        i += 1

                    in_alignment = True
                    continue

                # Process alignment
                if in_alignment:
                    if line.startswith("Q "):
                        parts = re.split(r'\s+', line.strip())
                        if len(parts) >= 4:
                            if 'query_start' not in current_hit:
                                try:
                                    current_hit['query_start'] = int(parts[2])
                                except ValueError:
                                    pass
                            query_ali += parts[3]

                    elif line.startswith("T "):
                        parts = re.split(r'\s+', line.strip())
                        if len(parts) >= 4:
                            if 'template_start' not in current_hit:
                                try:
                                    current_hit['template_start'] = int(parts[2])
                                except ValueError:
                                    pass
                            template_ali += parts[3]

                    # Blank line ends alignment section
                    elif line == "" and query_ali and template_ali:
                        current_hit['query_ali'] = query_ali
                        current_hit['template_ali'] = template_ali
                        in_alignment = False

                i += 1

            # Add the last hit
            if current_hit and 'query_ali' in current_hit and 'template_ali' in current_hit:
                hits.append(current_hit)

        except Exception as e:
            self.logger.error(f"Error parsing HHR file {hhr_file}: {str(e)}")

        return hits

    def _generate_summary_xml(self, pdb_id: str, chain_id: str, hits: List[Dict[str, Any]],
                             output_path: str
    ) -> None:
        """Generate XML summary of HHsearch results"""
        # Create XML document
        root = ET.Element("hh_summ_doc")

        # Add metadata
        metadata = ET.SubElement(root, "metadata")
        ET.SubElement(metadata, "pdb_id").text = pdb_id
        ET.SubElement(metadata, "chain_id").text = chain_id
        ET.SubElement(metadata, "reference").text = self.default_ref_version
        ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Add hit list
        hit_list = ET.SubElement(root, "hh_hit_list")

        for hit in hits:
            hit_node = ET.SubElement(hit_list, "hh_hit")

            # Basic hit attributes
            hit_node.set("hit_num", str(hit['hit_num']))
            hit_node.set("hh_prob", str(hit.get('probability', 0)))
            hit_node.set("hh_score", str(hit.get('score', 0)))
            hit_node.set("hh_evalue", str(hit.get('e_value', 999)))

            # Extract ECOD domain ID
            if 'hit_id' in hit and re.match(r"[dge]\d\w{3}\w+\d+", hit['hit_id']):
                hit_node.set("ecod_domain_id", hit['hit_id'])

            # Process alignment
            if 'query_ali' in hit and 'template_ali' in hit:
                query_ranges, coverage = self._map_alignment_to_ranges(
                    hit.get('query_ali', ''),
                    hit.get('template_ali', ''),
                    hit.get('query_start', 1),
                    hit.get('template_start', 1)
                )

                # Add query range
                query_range_node = ET.SubElement(hit_node, "query_range")
                query_range_node.text = query_ranges

                # Add coverage
                template_range_node = ET.SubElement(hit_node, "template_seqid_range")
                template_range_node.set("coverage", f"{coverage:.2f}")

                # Add alignment
                alignment_node = ET.SubElement(hit_node, "alignment")
                ET.SubElement(alignment_node, "query_ali").text = hit.get('query_ali', '')
                ET.SubElement(alignment_node, "template_ali").text = hit.get('template_ali', '')

        # Write XML to file
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding='utf-8', xml_declaration=True)

        self.logger.info(f"Generated summary XML: {output_path}")

    def _map_alignment_to_ranges(self, query_ali: str, template_ali: str,
                               query_start: int, template_start: int
    ) -> Tuple[str, float]:
        """Map alignment to sequence ranges and calculate coverage"""
        query_positions = []

        q_pos = query_start

        for q_char in query_ali:
            if q_char != '-':
                query_positions.append(q_pos)
                q_pos += 1
            else:
                q_pos += 1

        # Calculate coverage
        total_aligned = sum(1 for c in query_ali if c != '-')
        total_length = len(query_ali.replace('-', ''))
        coverage = total_aligned / total_length if total_length > 0 else 0.0

        # Convert to ranges
        ranges = self._positions_to_ranges(query_positions)

        return ranges, coverage

    def _positions_to_ranges(self, positions: List[int]) -> str:
        """Convert a list of positions to a range string"""
        if not positions:
            return ""

        positions = sorted(positions)
        ranges = []

        start = positions[0]
        prev = start

        for pos in positions[1:]:
            if pos > prev + 1:
                ranges.append(f"{start}-{prev}")
                start = pos
            prev = pos

        ranges.append(f"{start}-{prev}")
        return ",".join(ranges)

    def register_file(self, process_id: int, file_type: str,
                     file_path: str, file_exists: bool = True,
                     file_size: Optional[int] = None) -> None:
        """Register a process file with proper duplicate handling"""
        try:
            # Calculate file size if not provided
            if file_size is None and file_exists and os.path.exists(file_path):
                file_size = os.path.getsize(file_path)

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

    def get_chains_for_processing(self, batch_id: int, processing_stage: str = "profile") -> List[Dict[str, Any]]:
        """Get chains ready for specific processing stage

        Args:
            batch_id: Batch ID to process
            processing_stage: Which stage to find candidates for ('profile' or 'search')
                - 'profile': Chains with BLAST results but no HHblits profiles
                - 'search': Chains with HHblits profiles but no HHSearch results

        Returns:
            List of chains with all necessary information for processing
        """
        if processing_stage == "profile":
            # Get chains that need profile generation
            query = """
            SELECT
                p.id, p.pdb_id, p.chain_id, p.source_id, p.length,
                ps.id as process_id, ps.relative_path,
                seq.sequence
            FROM
                ecod_schema.protein p
            JOIN
                ecod_schema.process_status ps ON p.id = ps.protein_id
            JOIN
                ecod_schema.process_file pf ON ps.id = pf.process_id
            LEFT JOIN
                ecod_schema.protein_sequence seq ON p.id = seq.protein_id
            WHERE
                ps.batch_id = %s
                AND ps.status = 'success'
                AND pf.file_type IN ('chain_blast_result', 'domain_blast_result')
                AND pf.file_exists = TRUE
                AND NOT EXISTS (
                    SELECT 1 FROM ecod_schema.process_file
                    WHERE process_id = ps.id AND file_type = 'a3m'
                )
            GROUP BY
                p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.id, ps.relative_path, seq.sequence
            HAVING
                COUNT(DISTINCT pf.file_type) >= 1
            """
            stage_description = "profile generation (BLAST results but no HHblits profiles)"

        elif processing_stage == "search":
            # Get chains that need HHSearch (have profiles but no search results)
            query = """
            SELECT
                p.id, p.pdb_id, p.chain_id, p.source_id, p.length,
                ps.id as process_id, ps.relative_path,
                seq.sequence
            FROM
                ecod_schema.protein p
            JOIN
                ecod_schema.process_status ps ON p.id = ps.protein_id
            JOIN
                ecod_schema.process_file pf ON ps.id = pf.process_id
            LEFT JOIN
                ecod_schema.protein_sequence seq ON p.id = seq.protein_id
            WHERE
                ps.batch_id = %s
                AND ps.status = 'success'
                AND pf.file_type = 'a3m'
                AND pf.file_exists = TRUE
                AND NOT EXISTS (
                    SELECT 1 FROM ecod_schema.process_file
                    WHERE process_id = ps.id AND file_type = 'hhr'
                )
            GROUP BY
                p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.id, ps.relative_path, seq.sequence
            """
            stage_description = "HHSearch (have profiles but no search results)"

        else:
            self.logger.error(f"Invalid processing stage: {processing_stage}")
            return []

        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            self.logger.info(f"Found {len(rows)} chains ready for {stage_description}")

            # Filter out chains without sequences
            result = [row for row in rows if row.get('sequence')]

            if len(result) < len(rows):
                self.logger.warning(f"Filtered out {len(rows) - len(result)} chains without sequence data")

            return result
        except Exception as e:
            self.logger.error(f"Error querying chains for {stage_description}: {str(e)}")
            return []

    def adaptive_hhsearch_pipeline(self, batch_id: int, threads: int = 8,
                                  memory: str = "16G", adaptive: bool = True) -> Dict[str, int]:
        """Run adaptive HHSearch pipeline that processes chains based on their current state

        Args:
            batch_id: Batch ID to process
            threads: Number of threads for HHblits/HHsearch
            memory: Memory allocation for jobs
            adaptive: Whether to wait for profiles before running search

        Returns:
            Dictionary with counts of jobs submitted
        """
        # Get chains for different processing stages
        profile_candidates = self.get_chains_for_processing(batch_id, "profile")
        search_candidates = self.get_chains_for_processing(batch_id, "search")

        self.logger.info(f"Found {len(profile_candidates)} chains needing profiles and "
                       f"{len(search_candidates)} chains ready for search")

        profile_job_ids = []
        search_job_ids = []

        # Process profile generation first
        if profile_candidates:
            self.logger.info("Starting profile generation jobs")
            profile_job_ids = self.generate_profiles(batch_id, threads, memory)

            if adaptive and profile_job_ids:
                self.logger.info("Waiting for profile generation to complete before search")
                # Wait for profiles to complete
                self.job_manager.wait_for_jobs(profile_job_ids)

                # Get updated list of search candidates after profiles complete
                search_candidates = self.get_chains_for_processing(batch_id, "search")

        # Process search for all candidates with profiles
        if search_candidates:
            self.logger.info(f"Running HHSearch for {len(search_candidates)} chains")
            search_job_ids = self.run_hhsearch(batch_id, threads, memory)

        return {
            "profile_jobs": len(profile_job_ids),
            "search_jobs": len(search_job_ids)
        }
