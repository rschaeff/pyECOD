# ecod/pipelines/blast_pipeline.py
import os
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

from ecod.db import DBManager
from ecod.jobs import JobManager
from ecod.models import Protein, Batch, ProcessStatus, ProcessFile
from ecod.exceptions import PipelineError, JobSubmissionError, FileOperationError, ConfigurationError

class BlastPipeline:
    def __init__(self, db_manager: DBManager, job_manager: JobManager, config: Dict[str, Any]):
        self.db = db_manager
        self.job_manager = job_manager
        self.config = config
        self.logger = logging.getLogger("ecod.blast_pipeline")
        
        # Get tool paths from config
        tools = self.config.get('tools', {})
        self.blast_path = tools.get('blast_path', 'blastp')

        #Validate configuration
        self._validate_config()

    def _validate_config(self) -> None:
        """Validate configuration for the BLAST pipeline"""
        if not self.blast_path:
            error_msg = "BLAST path not configured"
            self.logger.error(error_msg)
            raise ConfigurationError(error_msg)
        
        # Check if BLAST database paths are configured
        ref = self.config.get('reference', {})
        self.chain_db = ref.get('chain_db')
        self.domain_db = ref.get('domain_db')
        
        if not self.chain_db:
            self.logger.warning("Chain BLAST database not configured")
        
        if not self.domain_db:
            self.logger.warning("Domain BLAST database not configured")
        
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
                ecod_schema.process_status ps ON p.id = ps.protein_id
            WHERE 
                d.id IS NULL
                AND ps.sequence IS NOT NULL
                AND p.length > 30
                AND ps.id IS NULL
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
                "ref_version": self.config.get('reference', {}).get('current_version', 'develop291'),
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
            ref_version=self.config.get('reference', {}).get('current_version', 'develop291'),
            total_items=len(chains),
            status="created"
        )
        
        # Register chains in this batch
        self._register_chains_in_batch(batch, chains)
        
        self.logger.info(f"Created BLAST batch: {batch_name} with {len(chains)} chains")
        return batch
        
    def _register_chains_in_batch(self, batch: Batch, chains: List[Protein]) -> None:
        """Register chains in a batch and create directory structure"""
        # Create query_fastas directory
        fasta_dir = os.path.join(batch.base_path, "query_fastas")
        os.makedirs(fasta_dir, exist_ok=True)
        
        # Create chain/domain results directories
        chain_dir = os.path.join(batch.base_path, "chain_blast_results")
        domain_dir = os.path.join(batch.base_path, "domain_blast_results")
        os.makedirs(chain_dir, exist_ok=True)
        os.makedirs(domain_dir, exist_ok=True)
        
        # Register each chain
        for chain in chains:
            # Determine the relative path for this chain
            rel_path = f"{chain.pdb_id}_{chain.chain_id}"
            
            # Register in process_status
            process_id = self.db.insert(
                "ecod_schema.process_status",
                {
                    "protein_id": chain.id,
                    "batch_id": batch.id,
                    "current_stage": "fasta_generated",
                    "status": "pending",
                    "relative_path": rel_path
                },
                "id"
            )
            
            # Generate FASTA file
            fasta_path = os.path.join(fasta_dir, f"{chain.source_id}.fa")
            with open(fasta_path, 'w') as f:
                f.write(f">{chain.source_id}\n{chain.sequence}\n")
            
            # Register FASTA file
            self.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "fasta",
                    "file_path": f"query_fastas/{chain.source_id}.fa",
                    "file_exists": True
                }
            )
            
    def generate_fasta_files(self, batch_id: int) -> List[Tuple[int, str]]:
        """Get FASTA files for a batch"""
        query = """
            SELECT 
                ps.id as process_id, 
                p.source_id,
                ps.relative_path,
                p_seq.sequence  
            FROM 
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.protein_sequence p_seq ON p.id = p_seq.protein_id
            WHERE 
                ps.batch_id = %s
        """
        
        rows = self.db.execute_dict_query(query, (batch_id,))
        
        batch_path = self._get_batch_path(batch_id)
        if not batch_path:
            self.logger.error(f"Batch {batch_id} not found")
            return []
            
        fasta_dir = os.path.join(batch_path, "query_fastas")
        os.makedirs(fasta_dir, exist_ok=True)
        
        fasta_paths = []
        for row in rows:
            fasta_path = os.path.join(fasta_dir, f"{row['source_id']}.fa")
            
            # Check if FASTA already exists
            if not os.path.exists(fasta_path):
                with open(fasta_path, 'w') as f:
                    f.write(f">{row['source_id']}\n{row['sequence']}\n")
                    
            fasta_paths.append((row['process_id'], fasta_path))
            
        self.logger.info(f"Generated/found {len(fasta_paths)} FASTA files for batch {batch_id}")
        return fasta_paths
        
    def _get_batch_path(self, batch_id: int) -> Optional[str]:
        """Get the base path for a batch"""
        query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
        rows = self.db.execute_query(query, (batch_id,))
        if rows:
            return rows[0][0]
        return None
  
    def run_domain_blast(self, batch_id: int, batch_size: int = 100) -> List[str]:
        """Run domain-wise BLAST on a batch of sequences"""
        # Get FASTA files for this batch
        fasta_paths = self.generate_fasta_files(batch_id)
        if not fasta_paths:
            self.logger.warning(f"No FASTA files found for batch {batch_id}")
            return []
            
        # Get batch information
        batch_path = self._get_batch_path(batch_id)
        if not batch_path:
            self.logger.error(f"Batch {batch_id} not found")
            return []
            
        # Create domain blast directory
        domain_blast_dir = os.path.join(batch_path, "domain_blast_results")
        os.makedirs(domain_blast_dir, exist_ok=True)
        
        # Get domain BLAST database path
        domain_db = self.config.get('reference', {}).get('domain_db')
        if not domain_db:
            self.logger.error("Domain BLAST database not configured")
            return []
        
    # Create job template
    def create_blast_command(process_id, fasta_path):
        basename = os.path.basename(fasta_path).replace('.fa', '')
        output_file = os.path.join(domain_blast_dir, f"{basename}.domain_blast.xml")
        return (f"{self.blast_path} -query {fasta_path} -db {domain_db} "
               f"-outfmt 5 -num_alignments 5000 -evalue 0.002 "
               f"-out {output_file}")
    
        job_template = {
            'name': f"domain_blast_{batch_id}",
            'output_dir': domain_blast_dir,
            'command_template': create_blast_command,
            'threads': 4,
            'memory': '4G',
            'time': '24:00:00'
        }
        
        # Create batch jobs
        jobs = self.job_manager.create_batch_jobs(fasta_paths, batch_size, job_template)
        
        # Submit jobs and record in database
        job_ids = []
        for job in jobs:
            # Submit job
            slurm_job_id = self.job_manager.submit_job(job['script_path'])
            if not slurm_job_id:
                continue
                
            # Record job in database
            job_db_id = self.db.insert(
                "ecod_schema.job",
                {
                    "batch_id": batch_id,
                    "job_type": "domain_blast",
                    "slurm_job_id": slurm_job_id,
                    "items_count": len(job['items'])
                },
                "id"
            )
            
            # Record items in this job
            for process_id, _ in job['items']:
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
            
        self.logger.info(f"Submitted {len(job_ids)} domain BLAST jobs for batch {batch_id}")
        return job_ids if job_ids else []
        
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

    def run_chain_blast(self, batch_id: int, batch_size: int = 100) -> List[str]:
        """Run chain-wise BLAST with improved error handling"""
        try:
            # Get FASTA files for this batch
            fasta_paths = self.generate_fasta_files(batch_id)
            if not fasta_paths:
                error_msg = f"No FASTA files found for batch {batch_id}"
                self.logger.warning(error_msg)
                return []

            # NEW: Check how many already have BLAST results
            existing_results = self._count_existing_results(batch_id, "chain_blast_result")
            self.logger.info(f"Found {existing_results} existing chain BLAST results for batch {batch_id}")
     
            # NEW: Skip if all files already have results (unless force flag is set)
            if existing_results >= len(fasta_paths) and not self.config.get('force_overwrite', False):
            self.logger.info(f"All {len(fasta_paths)} proteins already have chain BLAST results, skipping")
            return []

            # Get batch information
            batch_path = self._get_batch_path(batch_id)
            if not batch_path:
                error_msg = f"Batch {batch_id} not found"
                self.logger.error(error_msg)
                raise PipelineError(error_msg)
                
            # Create chain blast directory
            chain_blast_dir = os.path.join(batch_path, "chain_blast_results")
            os.makedirs(chain_blast_dir, exist_ok=True)
            
            # Check for chain BLAST database
            if not self.chain_db:
                error_msg = "Chain BLAST database not configured"
                self.logger.error(error_msg)
                raise ConfigurationError(error_msg)
            
            # Create job template
            def create_blast_command(process_id, fasta_path):
                basename = os.path.basename(fasta_path).replace('.fa', '')
                output_file = os.path.join(chain_blast_dir, f"{basename}.chainwise_blast.xml")
                return (f"{self.blast_path} -query {fasta_path} -db {self.chain_db} "
                      f"-outfmt 5 -num_alignments 5000 -evalue 0.002 "
                      f"-out {output_file}")
            
            job_template = {
                'name': f"chain_blast_{batch_id}",
                'output_dir': chain_blast_dir,
                'command_template': create_blast_command,
                'threads': 4,
                'memory': '4G',
                'time': '24:00:00'
            }
            
            # Create batch jobs
            try:
                jobs = self.job_manager.create_batch_jobs(fasta_paths, batch_size, job_template)
            except Exception as e:
                error_msg = f"Error creating BLAST jobs: {str(e)}"
                self.logger.error(error_msg)
                raise PipelineError(error_msg) from e
            
            # Submit jobs and record in database
            job_ids = []
            for job in jobs:
                try:
                    # Submit job
                    slurm_job_id = self.job_manager.submit_job(job['script_path'])
                    if not slurm_job_id:
                        self.logger.warning(f"Failed to submit job for {job['name']}")
                        continue
                        
                    # Record job in database
                    job_db_id = self.db.insert(
                        "ecod_schema.job",
                        {
                            "batch_id": batch_id,
                            "job_type": "chain_blast",
                            "slurm_job_id": slurm_job_id,
                            "items_count": len(job['items'])
                        },
                        "id"
                    )
                    
                    # Record items in this job
                    for process_id, _ in job['items']:
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
                    self.logger.info(f"Submitted job {slurm_job_id} for {job['name']}")
                    
                except JobSubmissionError as e:
                    # Log the error but continue with other jobs
                    self.logger.error(f"Job submission error for {job['name']}: {str(e)}")
                    
                    # Update process status for affected items
                    for process_id, _ in job['items']:
                        self.db.update(
                            "ecod_schema.process_status",
                            {
                                "current_stage": "chain_blast_failed",
                                "status": "error",
                                "error_message": f"Job submission failed: {str(e)}"
                            },
                            "id = %s",
                            (process_id,)
                        )
                
            self.logger.info(f"Submitted {len(job_ids)} chain BLAST jobs for batch {batch_id}")
            return job_ids
            
        except ConfigurationError:
            # Re-raise configuration errors
            raise
        except Exception as e:
            error_msg = f"Unexpected error in chain BLAST pipeline: {str(e)}"
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
        # Get items in this job
        query = """
            SELECT 
                ji.id, ji.process_id, ps.relative_path, b.base_path
            FROM 
                ecod_schema.job_item ji
            JOIN
                ecod_schema.process_status ps ON ji.process_id = ps.id
            JOIN
                ecod_schema.batch b ON ps.batch_id = b.id
            WHERE 
                ji.job_id = %s
        """
        
        rows = self.db.execute_dict_query(query, (job_id,))
        
        for row in rows:
            process_id = row['process_id']
            relative_path = row['relative_path']
            base_path = row['base_path']
            
            # Determine expected output file
            file_type = f"{job_type}_result"
            if job_type == "chain_blast":
                output_dir = os.path.join(base_path, "chain_blast_results")
                basename = relative_path.split('/')[-1]
                output_file = f"{output_dir}/{basename}.chainwise_blast.xml"
                relative_output = f"chain_blast_results/{basename}.chainwise_blast.xml"
            elif job_type == "domain_blast":
                output_dir = os.path.join(base_path, "domain_blast_results")
                basename = relative_path.split('/')[-1]
                output_file = f"{output_dir}/{basename}.domain_blast.xml"
                relative_output = f"domain_blast_results/{basename}.domain_blast.xml"
            else:
                self.logger.warning(f"Unknown job type: {job_type}")
                continue
            
            # Check if output file exists
            file_exists = os.path.exists(output_file)
            file_size = os.path.getsize(output_file) if file_exists else 0
            
            # Update file registry
            file_id = self.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": file_type,
                    "file_path": relative_output,
                    "file_exists": file_exists,
                    "file_size": file_size
                },
                "id"
            )
            
            # Determine next stage
            next_stage = "blast_completed" if file_exists and file_size > 0 else "blast_failed"
            next_status = "success" if file_exists and file_size > 0 else "error"
            
            # Update process status
            self.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": next_stage,
                    "status": next_status
                },
                "id = %s",
                (process_id,)
            )
            
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
        
        for row in rows:
            source_id = row['source_id']
            file_path = row['file_path']
            file_type = row['file_type']
            
            if file_type == 'chain_blast_result':
                results['chain_results'][source_id] = file_path
            elif file_type == 'domain_blast_result':
                results['domain_results'][source_id] = file_path
        
        self.logger.info(f"Found {len(results['chain_results'])} chain results and {len(results['domain_results'])} domain results")
        return results

    def parse_and_store_blast_results(self, batch_id: int) -> Dict[str, int]:
        """
        Parse BLAST result files and store in database
        
        Args:
            batch_id: Batch ID
            
        Returns:
            Dictionary with counts of hits stored by type
        """
        # Get batch information
        batch_query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
        batch_rows = self.db.execute_query(batch_query, (batch_id,))
        
        if not batch_rows:
            self.logger.error(f"Batch {batch_id} not found")
            return {"domain_hits": 0, "chain_hits": 0}
        
        batch_path = batch_rows[0][0]
        
        # Get processes with BLAST result files
        query = """
            SELECT 
                ps.id, ps.protein_id, p.pdb_id, p.chain_id, pf.file_path, pf.file_type
            FROM 
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.process_file pf ON ps.id = pf.process_id
            WHERE 
                ps.batch_id = %s
                AND pf.file_type IN ('domain_blast_result', 'chain_blast_result')
                AND pf.file_exists = TRUE
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
        except Exception as e:
            self.logger.error(f"Error querying BLAST files: {e}")
            return {"domain_hits": 0, "chain_hits": 0}
        
        self.logger.info(f"Found {len(rows)} BLAST result files")
        
        # Organize files by protein and type
        blast_files = {}
        for row in rows:
            protein_id = row['protein_id']
            file_type = row['file_type']
            file_path = os.path.join(batch_path, row['file_path'])
            
            if protein_id not in blast_files:
                blast_files[protein_id] = {'protein_id': protein_id, 'pdb_id': row['pdb_id'], 'chain_id': row['chain_id']}
                
            blast_files[protein_id][file_type] = file_path
        
        # Process each protein's BLAST files
        domain_hit_count = 0
        chain_hit_count = 0
        
        for protein_id, files in blast_files.items():
            # Process domain BLAST results
            if 'domain_blast_result' in files:
                file_path = files['domain_blast_result']
                domain_hits = self._parse_domain_blast_file(file_path, protein_id, files['pdb_id'], files['chain_id'])
                
                if domain_hits:
                    domain_hit_count += self._store_domain_blast_hits(protein_id, domain_hits)
                    
            # Process chain BLAST results
            if 'chain_blast_result' in files:
                file_path = files['chain_blast_result']
                chain_hits = self._parse_chain_blast_file(file_path, protein_id, files['pdb_id'], files['chain_id'])
                
                if chain_hits:
                    chain_hit_count += self._store_chain_blast_hits(protein_id, chain_hits)
        
        self.logger.info(f"Stored {domain_hit_count} domain hits and {chain_hit_count} chain hits")
        return {"domain_hits": domain_hit_count, "chain_hits": chain_hit_count}

    def _parse_domain_blast_file(self, file_path: str, protein_id: int, 
                               pdb_id: str, chain_id: str) -> List[Dict[str, Any]]:
        """
        Parse domain BLAST XML file with the format observed in examples
        
        Args:
            file_path: Path to XML file
            protein_id: Protein ID
            pdb_id: PDB ID
            chain_id: Chain ID
            
        Returns:
            List of hit dictionaries
        """
        try:
            import xml.etree.ElementTree as ET
            tree = ET.parse(file_path)
            root = tree.getroot()
            
            hits = []
            
            # Process hits
            for iteration in root.findall(".//Iteration"):
                # Get query information
                query_len = int(iteration.findtext(".//Iteration_query-len", "0"))
                
                for hit in iteration.findall(".//Hit"):
                    hit_num = hit.findtext("Hit_num", "")
                    hit_def = hit.findtext("Hit_def", "")
                    hit_len = int(hit.findtext("Hit_len", "0"))
                    
                    # Parse hit definition based on the format in the example
                    # Format: e8b7oAAA1 AAA:4-183 002982980
                    domain_id = "unknown"
                    hit_pdb_id = "unknown"
                    hit_chain_id = "unknown"
                    
                    # Try to extract domain ID and other info
                    import re
                    # Pattern for the example: e8b7oAAA1 AAA:4-183 002982980
                    match = re.search(r'([edg]\d\w{3}\w+\d*)\s+([A-Z]+):(\d+)-(\d+)', hit_def)
                    if match:
                        domain_id = match.group(1)  # e8b7oAAA1
                        hit_chain_id = match.group(2)  # AAA
                        # Extract PDB ID from domain ID (first 5 characters)
                        if len(domain_id) >= 5:
                            hit_pdb_id = domain_id[1:5]  # 8b7o
                    
                    # Process HSPs for this hit
                    hsps = []
                    evalues = []
                    query_regions = []
                    hit_regions = []
                    query_seqs = []
                    hit_seqs = []
                    
                    for hsp in hit.findall(".//Hsp"):
                        hsp_evalue = float(hsp.findtext("Hsp_evalue", "999"))
                        
                        # Get alignment coordinates
                        hsp_query_from = int(hsp.findtext("Hsp_query-from", "0"))
                        hsp_query_to = int(hsp.findtext("Hsp_query-to", "0"))
                        hsp_hit_from = int(hsp.findtext("Hsp_hit-from", "0"))
                        hsp_hit_to = int(hsp.findtext("Hsp_hit-to", "0"))
                        
                        # Get aligned sequences
                        hsp_qseq = hsp.findtext("Hsp_qseq", "")
                        hsp_hseq = hsp.findtext("Hsp_hseq", "")
                        
                        # Calculate identity percentage
                        hsp_identity = int(hsp.findtext("Hsp_identity", "0"))
                        hsp_align_len = int(hsp.findtext("Hsp_align-len", "0"))
                        identity_pct = hsp_identity / hsp_align_len if hsp_align_len > 0 else 0
                        
                        # Store HSP details
                        hsps.append({
                            'evalue': hsp_evalue,
                            'query_from': hsp_query_from,
                            'query_to': hsp_query_to,
                            'hit_from': hsp_hit_from,
                            'hit_to': hsp_hit_to,
                            'qseq': hsp_qseq,
                            'hseq': hsp_hseq,
                            'identity': identity_pct
                        })
                        
                        evalues.append(str(hsp_evalue))
                        query_regions.append(f"{hsp_query_from}-{hsp_query_to}")
                        hit_regions.append(f"{hsp_hit_from}-{hsp_hit_to}")
                        query_seqs.append(hsp_qseq)
                        hit_seqs.append(hsp_hseq)
                    
                    # Skip hits with no HSPs
                    if not hsps:
                        continue
                    
                    # Find best identity percentage from all HSPs
                    best_identity = max(hsp['identity'] for hsp in hsps) if hsps else 0
                    
                    # Add hit
                    hits.append({
                        'domain_id': domain_id,
                        'pdb_id': hit_pdb_id,
                        'chain_id': hit_chain_id,
                        'evalues': ",".join(evalues),
                        'query_regions': ",".join(query_regions),
                        'hit_regions': ",".join(hit_regions),
                        'query_seqs': ",".join(query_seqs),
                        'hit_seqs': ",".join(hit_seqs),
                        'identity': best_identity
                    })
            
            return hits
        except Exception as e:
            self.logger.error(f"Error parsing domain BLAST file {file_path}: {e}")
            return []

    def _parse_chain_blast_file(self, file_path: str, protein_id: int, 
                              pdb_id: str, chain_id: str) -> List[Dict[str, Any]]:
        """
        Parse chain BLAST XML file with the format observed in examples
        
        Args:
            file_path: Path to XML file
            protein_id: Protein ID
            pdb_id: PDB ID
            chain_id: Chain ID
            
        Returns:
            List of hit dictionaries
        """
        try:
            import xml.etree.ElementTree as ET
            tree = ET.parse(file_path)
            root = tree.getroot()
            
            hits = []
            
            # Process hits
            for iteration in root.findall(".//Iteration"):
                # Get query information
                query_len = int(iteration.findtext(".//Iteration_query-len", "0"))
                
                for hit in iteration.findall(".//Hit"):
                    hit_num = hit.findtext("Hit_num", "")
                    hit_def = hit.findtext("Hit_def", "")
                    hit_len = int(hit.findtext("Hit_len", "0"))
                    
                    # Parse hit definition based on the format in the example
                    # Format: "1b3b F"
                    hit_pdb_id = "unknown"
                    hit_chain_id = "unknown"
                    
                    # Try to extract PDB ID and chain ID
                    import re
                    # Pattern for the example: "1b3b F"
                    match = re.search(r'(\d\w{3})\s+([A-Za-z0-9])', hit_def)
                    if match:
                        hit_pdb_id = match.group(1)  # 1b3b
                        hit_chain_id = match.group(2)  # F
                    
                    # Process HSPs for this hit
                    hsps = []
                    evalues = []
                    query_regions = []
                    hit_regions = []
                    query_seqs = []
                    hit_seqs = []
                    
                    for hsp in hit.findall(".//Hsp"):
                        hsp_evalue = float(hsp.findtext("Hsp_evalue", "999"))
                        
                        # Get alignment coordinates
                        hsp_query_from = int(hsp.findtext("Hsp_query-from", "0"))
                        hsp_query_to = int(hsp.findtext("Hsp_query-to", "0"))
                        hsp_hit_from = int(hsp.findtext("Hsp_hit-from", "0"))
                        hsp_hit_to = int(hsp.findtext("Hsp_hit-to", "0"))
                        
                        # Get aligned sequences
                        hsp_qseq = hsp.findtext("Hsp_qseq", "")
                        hsp_hseq = hsp.findtext("Hsp_hseq", "")
                        
                        # Calculate identity percentage
                        hsp_identity = int(hsp.findtext("Hsp_identity", "0"))
                        hsp_align_len = int(hsp.findtext("Hsp_align-len", "0"))
                        identity_pct = hsp_identity / hsp_align_len if hsp_align_len > 0 else 0
                        
                        # Store HSP details
                        hsps.append({
                            'evalue': hsp_evalue,
                            'query_from': hsp_query_from,
                            'query_to': hsp_query_to,
                            'hit_from': hsp_hit_from,
                            'hit_to': hsp_hit_to,
                            'qseq': hsp_qseq,
                            'hseq': hsp_hseq,
                            'identity': identity_pct
                        })
                        
                        evalues.append(str(hsp_evalue))
                        query_regions.append(f"{hsp_query_from}-{hsp_query_to}")
                        hit_regions.append(f"{hsp_hit_from}-{hsp_hit_to}")
                        query_seqs.append(hsp_qseq)
                        hit_seqs.append(hsp_hseq)
                    
                    # Skip hits with no HSPs
                    if not hsps:
                        continue
                    
                    # Find best identity percentage from all HSPs
                    best_identity = max(hsp['identity'] for hsp in hsps) if hsps else 0
                    
                    # Add hit
                    hits.append({
                        'pdb_id': hit_pdb_id,
                        'chain_id': hit_chain_id,
                        'evalues': ",".join(evalues),
                        'query_regions': ",".join(query_regions),
                        'hit_regions': ",".join(hit_regions),
                        'query_seqs': ",".join(query_seqs),
                        'hit_seqs': ",".join(hit_seqs),
                        'identity': best_identity
                    })
            
            return hits
        except Exception as e:
            self.logger.error(f"Error parsing chain BLAST file {file_path}: {e}")
            return []

    def _store_domain_blast_hits(self, protein_id: int, hits: List[Dict[str, Any]]) -> int:
        """
        Store domain BLAST hits in database
        
        Args:
            protein_id: Protein ID
            hits: List of hit dictionaries
            
        Returns:
            Number of hits stored
        """
        if not hits:
            return 0
        
        stored_count = 0
        
        for hit in hits:
            try:
                # Check if the minimum e-value is reasonable
                evalues = [float(e) for e in hit['evalues'].split(',') if e.strip()]
                min_evalue = min(evalues) if evalues else 10.0
                
                # Insert hit record
                self.db.insert(
                    "ecod_schema.domain_blast_hit",
                    {
                        "protein_id": protein_id,
                        "domain_id": hit['domain_id'],
                        "pdb_id": hit['pdb_id'],
                        "chain_id": hit['chain_id'],
                        "evalue": min_evalue,
                        "query_regions": hit['query_regions'],
                        "hit_regions": hit['hit_regions'],
                        "query_seqs": hit['query_seqs'],
                        "hit_seqs": hit['hit_seqs'],
                        "identity": hit['identity']
                    }
                )
                
                stored_count += 1
            except Exception as e:
                self.logger.error(f"Error storing domain BLAST hit for protein {protein_id}: {e}")
        
        return stored_count

    def _store_chain_blast_hits(self, protein_id: int, hits: List[Dict[str, Any]]) -> int:
        """
        Store chain BLAST hits in database
        
        Args:
            protein_id: Protein ID
            hits: List of hit dictionaries
            
        Returns:
            Number of hits stored
        """
        if not hits:
            return 0
        
        stored_count = 0
        
        for hit in hits:
            try:
                # Check if the minimum e-value is reasonable
                evalues = [float(e) for e in hit['evalues'].split(',') if e.strip()]
                min_evalue = min(evalues) if evalues else 10.0
                
                # Insert hit record
                self.db.insert(
                    "ecod_schema.chain_blast_hit",
                    {
                        "protein_id": protein_id,
                        "pdb_id": hit['pdb_id'],
                        "chain_id": hit['chain_id'],
                        "evalue": min_evalue,
                        "query_regions": hit['query_regions'],
                        "hit_regions": hit['hit_regions'],
                        "query_seqs": hit['query_seqs'],
                        "hit_seqs": hit['hit_seqs'],
                        "identity": hit['identity']
                    }
                )
                
                stored_count += 1
            except Exception as e:
                self.logger.error(f"Error storing chain BLAST hit for protein {protein_id}: {e}")
        
        return stored_count
        