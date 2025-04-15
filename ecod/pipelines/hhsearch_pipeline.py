# ecod/pipelines/hhsearch_pipeline.py
import os
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import re
import xml.etree.ElementTree as ET
from datetime import datetime

from ecod.db import DBManager
from ecod.jobs import JobManager
from ecod.models import Batch, ProcessStatus

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
            check_query = """
            SELECT id FROM ecod_schema.process_status
            WHERE protein_id = %s and file_type = %s
            """

            existing = db.execute_query(check_query, (protein_id, file_type))

            process_id = None;
            if existing:
                process_id = self.db.update(
                    "ecod_schema.process_status",
                    {
                        "protein_id": protein_id,
                        "batch_id": batch_id,
                        "current_status": "fasta_generated",
                        "status": "pending",
                        "relative_path": rel_path
                    }, "id = %s"),
                (existing[0][0],)
            else:
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
                ps.id 
            FROM 
                ecod_schema.process_status ps
            WHERE 
                ps.batch_id = %s
                AND ps.protein_id IN %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id, tuple(protein_ids)))
            process_ids = [row['id'] for row in rows]
        except Exception as e:
            self.logger.error(f"Error querying processes for proteins: {e}")
            return []
        
        # Submit HHSearch jobs for these processes
        job_ids = []
        for process_id in process_ids:
            job_id = self._submit_hhsearch_job(process_id, threads, memory)
            if job_id:
                job_ids.append(job_id)
        
        return job_ids


    def _submit_hhblits_job(self, process_id: int, pdb_id: str, chain_id: str, 
                           rel_path: str, base_path: str, threads: int, memory: str
    ) -> Optional[str]:
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
                "status": "processing"
            },
            "id = %s",
            (process_id,)
        )
        
        return slurm_job_id
        
    # Additional methods for running HHsearch, parsing results, etc.
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
                chain['relative_path'],
                chain['base_path'],
                threads,
                memory,
                batch_id
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
            """
            
            rows = self.db.execute_dict_query(query, (batch_id,))
            return rows

        def _submit_hhsearch_job(self, process_id: int, pdb_id: str, chain_id: str, 
                                rel_path: str, base_path: str, threads: int, memory: str,
                                batch_id: int) -> Optional[str]:
            """Submit a job to run HHsearch for a single chain"""
            chain_dir = os.path.join(base_path, "ecod_dump", rel_path)
            a3m_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.a3m")
            hhm_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.hhm")
            result_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.{self.config.get('reference', {}).get('current_version', 'develop291')}.hhr")
            
            # Check if A3M exists
            if not os.path.exists(a3m_file):
                self.logger.warning(f"A3M file not found: {a3m_file}")
                return None
            
            # Create job script
            job_name = f"hhsearch_{pdb_id}_{chain_id}"
            job_script = os.path.join(chain_dir, f"{job_name}.sh")
            
            commands = [
                f"module purge",
                f"module load hh-suite",
                f"{self.hhmake_path} -i {a3m_file} -o {hhm_file}",
                f"{self.hhsearch_path} "
                f"-i {hhm_file} "
                f"-d {self.ecod_ref_db} "
                f"-o {result_file} "
                f"-cpu {threads} "
                f"-v 2 -p 60.0 -z 100 -b 100 -ssm 2 -sc 1 -aliw 80 -glob"
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
                    "job_type": "hhsearch",
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
        
        self.logger.info(f"Checking status of {len(rows)} HHSearch jobs")
        
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
                ji.id, ji.process_id, ps.relative_path, b.base_path, p.pdb_id, p.chain_id
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
        ref_version = self.config.get('reference', {}).get('current_version', 'develop291')
        
        for row in rows:
            process_id = row['process_id']
            relative_path = row['relative_path']
            base_path = row['base_path']
            pdb_id = row['pdb_id']
            chain_id = row['chain_id']
            
            # Determine expected output file and next stage
            if job_type == "hhblits":
                chain_dir = os.path.join(base_path, "ecod_dump", relative_path)
                output_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.a3m")
                relative_output = f"ecod_dump/{relative_path}/{pdb_id}_{chain_id}.a3m"
                file_type = "a3m"
                next_stage = "profile_complete"
                next_status = "success"
            elif job_type == "hhsearch":
                chain_dir = os.path.join(base_path, "ecod_dump", relative_path)
                output_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.{ref_version}.hhr")
                relative_output = f"ecod_dump/{relative_path}/{pdb_id}_{chain_id}.{ref_version}.hhr"
                file_type = "hhr"
                next_stage = "hhsearch_complete"
                next_status = "success"
            else:
                self.logger.warning(f"Unknown job type: {job_type}")
                continue
            
            # Check if output file exists
            file_exists = os.path.exists(output_file)
            file_size = os.path.getsize(output_file) if file_exists else 0
            
            if file_exists and file_size > 0:
                # Update file registry
                file_id = self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": relative_output,
                        "file_exists": True,
                        "file_size": file_size
                    },
                    "id"
                )
                
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
            else:
                # File missing or empty
                self.db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": "failed",
                        "status": "error",
                        "error_message": f"Output file missing or empty: {output_file}"
                    },
                    "id = %s",
                    (process_id,)
                )

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
        ET.SubElement(metadata, "reference").text = self.config.get('reference', {}).get('current_version', 'develop291')
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