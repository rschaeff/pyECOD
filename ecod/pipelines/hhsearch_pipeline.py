# ecod/pipelines/hhsearch_pipeline.py
import os
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import re
import xml.etree.ElementTree as ET
from datetime import datetime

from ..core.db_manager import DBManager
from ..core.job_manager import JobManager
from ..core.models import Batch, ProcessStatus

class HHSearchPipeline:
    def __init__(self, db_manager: DBManager, job_manager: JobManager, config: Dict[str, Any]):
        self.db = db_manager
        self.job_manager = job_manager
        self.config = config
        self.logger = logging.getLogger("ecod.hhsearch_pipeline")
        
        # Get tool paths from config
        tools = self.config.get('tools', {})
        self.hhblits_path = tools.get('hhblits_path', 'hhblits')
        self.hhsearch_path = tools.get('hhsearch_path', 'hhsearch')
        self.hhmake_path = tools.get('hhmake_path', 'hhmake')
        
        # Get reference paths
        ref = self.config.get('reference', {})
        self.uniclust_db = ref.get('uniclust_db')
        self.ecod_ref_db = ref.get('ecod_hh_db')
        
    def create_batch(self, chains: List[Dict[str, Any]], batch_size: int = 500) -> int:
        """Create a new batch for HHsearch processing"""
        # Generate batch name with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M")
        batch_name = f"hhsearch_batch_{timestamp}"
        
        # Create batch directory
        base_dir = self.config.get('paths', {}).get('output_dir', '/data/ecod/weekly_updates/batches')
        batch_path = os.path.join(base_dir, batch_name)
        dump_path = os.path.join(batch_path, "ecod_dump")
        os.makedirs(batch_path, exist_ok=True)
        os.makedirs(dump_path, exist_ok=True)
        
        # Create batch record in database
        batch_id = self.db.insert(
            "ecod_schema.batch",
            {
                "batch_name": batch_name,
                "base_path": batch_path,
                "type": "hhsearch",
                "ref_version": self.config.get('reference', {}).get('current_version', 'develop291'),
                "total_items": len(chains),
                "status": "created"
            },
            "id"
        )
        
        # Register chains in this batch
        for chain in chains:
            pdb_id = chain['pdb_id']
            chain_id = chain['chain_id']
            protein_id = chain['id']
            
            # Create relative path
            rel_path = f"{pdb_id}_{chain_id}"
            
            # Create chain directory
            chain_dir = os.path.join(dump_path, rel_path)
            os.makedirs(chain_dir, exist_ok=True)
            
            # Register in process_status
            process_id = self.db.insert(
                "ecod_schema.process_status",
                {
                    "protein_id": protein_id,
                    "batch_id": batch_id,
                    "current_stage": "fasta_generated",
                    "status": "pending",
                    "relative_path": rel_path
                },
                "id"
            )
            
            # Create FASTA file
            fasta_path = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.fa")
            with open(fasta_path, 'w') as f:
                f.write(f">{pdb_id}_{chain_id}\n{chain['sequence']}\n")
            
            # Register FASTA file
            self.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "fasta",
                    "file_path": f"ecod_dump/{rel_path}/{pdb_id}_{chain_id}.fa",
                    "file_exists": True
                }
            )
            
        self.logger.info(f"Created HHsearch batch: {batch_name} with {len(chains)} chains")
        return batch_id
        
    def get_chains_for_profile_generation(self, batch_id: int) -> List[Dict[str, Any]]:
        """Get chains ready for profile generation"""
        query = """
            SELECT 
                ps.id, p.id as protein_id, p.pdb_id, p.chain_id, p.sequence,
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
                AND ps.current_stage = 'fasta_generated'
                AND ps.status = 'pending'
        """
        
        rows = self.db.execute_dict_query(query, (batch_id,))
        return rows
        
    def generate_profiles(self, batch_id: int, threads: int = 8, memory: str = "16G") -> List[str]:
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
                chain['relative_path'],
                chain['base_path'],
                threads,
                memory
            )
            
            if job_id:
                job_ids.append(job_id)
                
        self.logger.info(f"Submitted {len(job_ids)} HHblits profile generation jobs")
        return job_ids
        
    def _submit_hhblits_job(self, process_id: int, pdb_id: str, chain_id: str, 
                           rel_path: str, base_path: str, threads: int, memory: str) -> Optional[str]:
        """Submit a job to generate HHblits profile for a chain"""
        chain_dir = os.path.join(base_path, "ecod_dump", rel_path)
        fa_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.fa")
        a3m_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.a3m")
        
        # Check if FASTA exists
        if not os.path.exists(fa_file):
            self.logger.warning(f"FASTA file not found: {fa_file}")
            return None
            
        # Create job script
        job_name = f"hhblits_{pdb_id}_{chain_id}"
        job_script = os.path.join(chain_dir, f"{job_name}.sh")
        
        commands = [
            f"module purge",
            f"module load hh-suite",
            f"{self.hhblits_path} -i {fa_file} -oa3m {a3m_file} "
            f"-d {self.uniclust_db} "
            f"-n 3 -maxfilt 20000 -diff inf -id 99 -cov 50 "
            f"-cpu {threads}"
        ]
        
        self.job_manager.create_job_script(
            commands,
            job_name,
            chain_dir,
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
                "items_count": 1
            },
            "id"
        )
        
        # Record chain in this job
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
                "current_stage": "profile_running",
                "status": "processing",
                "updated_at": "CURRENT_TIMESTAMP"
            },
            "id = %s",
            (process_id,)
        )
        
        return slurm_job_id
        
    # Additional methods for running HHsearch, parsing results, etc.