#!/usr/bin/env python3
"""
hhsearch_runner.py - Execute HHblits and HHsearch for batches of proteins
"""

import os
import sys
import argparse
import logging
import re
import xml.etree.ElementTree as ET
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

# Add parent directory to path if needed
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from ecod.core.config import ConfigManager
from ecod.core.db_manager import DBManager
from ecod.core.job_manager import JobManager

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

class HHSearchRunner:
    """Class to run HHblits and HHsearch for batch processing"""
    
    def __init__(self, config_path: str = None):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        self.db_config = self.config_manager.get_db_config()
        self.db = DBManager(self.db_config)
        self.job_manager = JobManager(self.config)
        self.logger = logging.getLogger("ecod.hhsearch_runner")
        
        # Get HHsuite paths from config
        tools = self.config.get('tools', {})
        self.hhblits_path = tools.get('hhblits_path', 'hhblits')
        self.hhsearch_path = tools.get('hhsearch_path', 'hhsearch')
        self.hhmake_path = tools.get('hhmake_path', 'hhmake')
        
        # Get reference DBs from config
        reference = self.config.get('reference', {})
        self.uniclust_db = reference.get('uniclust_db')
        self.ecod_ref_db = reference.get('ecod_hh_db')
        self.current_version = reference.get('current_version', 'develop291')
        
        if not self.uniclust_db or not self.ecod_ref_db:
            self.logger.warning("HHsearch databases not configured properly")
    
    def get_pending_batches(self):
        """Get batches that are not completed"""
        query = """
        SELECT 
            id, batch_name, type, total_items, completed_items, status, base_path
        FROM 
            ecod_schema.batch
        WHERE 
            status != 'completed'
            AND type IN ('hhsearch', 'demo')
        ORDER BY 
            id
        """
        return self.db.execute_dict_query(query)
    
    def get_hhblits_ready_items(self, batch_id: int, limit: int = 10):
        """Get items ready for HHblits in a batch"""
        query = """
        SELECT 
            ps.id, p.id AS protein_id, p.pdb_id, p.chain_id, p.source_id,
            ps.current_stage, ps.status, ps.relative_path, b.base_path,
            b.id AS batch_id
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
            AND ps.status IN ('pending', 'success')
            AND ps.current_stage IN ('fasta_generated', 'completed')
            AND pf.file_type = 'fasta'
            AND pf.file_exists = TRUE
            AND NOT EXISTS (
                SELECT 1 FROM ecod_schema.process_file 
                WHERE process_id = ps.id AND file_type = 'a3m'
            )
        ORDER BY 
            ps.id
        LIMIT %s
        """
        return self.db.execute_dict_query(query, (batch_id, limit))
    
    def get_hhsearch_ready_items(self, batch_id: int, limit: int = 10):
        """Get items ready for HHsearch in a batch"""
        query = """
        SELECT 
            ps.id, p.id AS protein_id, p.pdb_id, p.chain_id, p.source_id,
            ps.current_stage, ps.status, ps.relative_path, b.base_path,
            b.id AS batch_id
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
            AND ps.status = 'success'
            AND ps.current_stage = 'profile_complete'
            AND pf.file_type = 'a3m'
            AND pf.file_exists = TRUE
            AND NOT EXISTS (
                SELECT 1 FROM ecod_schema.process_file 
                WHERE process_id = ps.id AND file_type = 'hhr'
            )
        ORDER BY 
            ps.id
        LIMIT %s
        """
        return self.db.execute_dict_query(query, (batch_id, limit))
    
    def prepare_hhblits_jobs(self, items: List[Dict[str, Any]], threads: int = 8,
                          batch_size: int = 1) -> List[Dict]:
        """Prepare HHblits jobs for a list of items"""
        # HHblits is CPU-intensive, so we typically run one per job
        jobs = []
        
        for i, item in enumerate(items):
            process_id = item['id']
            pdb_id = item['pdb_id']
            chain_id = item['chain_id']
            source_id = item['source_id']
            relative_path = item['relative_path']
            base_path = item['base_path']
            batch_id = item['batch_id']
            
            # Create ecod_dump directory structure
            dump_dir = os.path.join(base_path, "ecod_dump")
            chain_dir = os.path.join(dump_dir, relative_path)
            os.makedirs(chain_dir, exist_ok=True)
            
            # Define file paths
            fasta_path = os.path.join(base_path, "query_fastas", f"{source_id}.fa")
            if not os.path.exists(fasta_path):
                self.logger.warning(f"FASTA file not found: {fasta_path}")
                continue
                
            # Copy FASTA to chain directory if needed
            chain_fasta = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.fa")
            if not os.path.exists(chain_fasta):
                with open(fasta_path, 'r') as src, open(chain_fasta, 'w') as dst:
                    dst.write(src.read())
            
            # Define output paths
            a3m_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.a3m")
            
            # Create job name
            job_name = f"hhblits_{pdb_id}_{chain_id}"
            
            # Create commands
            commands = [
                "module purge",
                "module load hh-suite",
                f"{self.hhblits_path} -i {chain_fasta} -oa3m {a3m_file} "
                f"-d {self.uniclust_db} "
                f"-n 3 -maxfilt 20000 -diff inf -id 99 -cov 50 "
                f"-cpu {threads}"
            ]
            
            # Create job script
            script_path = self.job_manager.create_job_script(
                commands,
                job_name,
                chain_dir,
                threads=threads,
                memory="16G",
                time="24:00:00"
            )
            
            jobs.append({
                'name': job_name,
                'script_path': script_path,
                'items': [item],
                'output_dir': chain_dir,
                'output_file': a3m_file,
                'batch_id': batch_id
            })
        
        return jobs
    
    def prepare_hhsearch_jobs(self, items: List[Dict[str, Any]], threads: int = 8,
                           batch_size: int = 1) -> List[Dict]:
        """Prepare HHsearch jobs for a list of items"""
        # HHsearch is also CPU-intensive, so we typically run one per job
        jobs = []
        
        for i, item in enumerate(items):
            process_id = item['id']
            pdb_id = item['pdb_id']
            chain_id = item['chain_id']
            source_id = item['source_id']
            relative_path = item['relative_path']
            base_path = item['base_path']
            batch_id = item['batch_id']
            
            # Get chain directory
            dump_dir = os.path.join(base_path, "ecod_dump")
            chain_dir = os.path.join(dump_dir, relative_path)
            
            # Define file paths
            a3m_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.a3m")
            hhm_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.hhm")
            hhr_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.{self.current_version}.hhr")
            
            # Verify A3M exists
            if not os.path.exists(a3m_file):
                self.logger.warning(f"A3M file not found: {a3m_file}")
                continue
            
            # Create job name
            job_name = f"hhsearch_{pdb_id}_{chain_id}"
            
            # Create commands
            commands = [
                "module purge",
                "module load hh-suite",
                f"{self.hhmake_path} -i {a3m_file} -o {hhm_file}",
                f"{self.hhsearch_path} "
                f"-i {hhm_file} "
                f"-d {self.ecod_ref_db} "
                f"-o {hhr_file} "
                f"-cpu {threads} "
                f"-v 2 -p 60.0 -z 100 -b 100 -ssm 2 -sc 1 -aliw 80 -glob"
            ]
            
            # Create job script
            script_path = self.job_manager.create_job_script(
                commands,
                job_name,
                chain_dir,
                threads=threads,
                memory="16G",
                time="24:00:00"
            )
            
            jobs.append({
                'name': job_name,
                'script_path': script_path,
                'items': [item],
                'output_dir': chain_dir,
                'output_file': hhr_file,
                'batch_id': batch_id
            })
        
        return jobs
    
    def submit_jobs(self, jobs: List[Dict], job_type: str) -> List[str]:
        """Submit jobs to the cluster"""
        job_ids = []
        
        for job in jobs:
            # Submit job
            slurm_job_id = self.job_manager.submit_job(job['script_path'])
            if not slurm_job_id:
                self.logger.error(f"Failed to submit job {job['name']}")
                continue
            
            self.logger.info(f"Submitted job {job['name']} with ID {slurm_job_id}")
            
            # Record job in database
            job_db_id = self.db.insert(
                "ecod_schema.job",
                {
                    "batch_id": job['batch_id'],
                    "job_type": job_type,
                    "slurm_job_id": slurm_job_id,
                    "status": "submitted",
                    "items_count": len(job['items'])
                },
                "id"
            )
            
            # Record items in this job
            for item in job['items']:
                process_id = item['id']
                
                self.db.insert(
                    "ecod_schema.job_item",
                    {
                        "job_id": job_db_id,
                        "process_id": process_id
                    }
                )
                
                # Update process status
                stage = "profile_running" if job_type == "hhblits" else "hhsearch_running"
                self.db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": stage,
                        "status": "processing"
                    },
                    "id = %s",
                    (process_id,)
                )
            
            job_ids.append(slurm_job_id)
        
        return job_ids
    
    def check_running_jobs(self):
        """Check status of running HHsuite jobs"""
        query = """
        SELECT 
            j.id, j.slurm_job_id, j.job_type, j.batch_id, j.status
        FROM 
            ecod_schema.job j
        WHERE 
            j.status = 'submitted'
            AND j.job_type IN ('hhblits', 'hhsearch')
        ORDER BY
            j.id
        """
        
        jobs = self.db.execute_dict_query(query)
        
        for job in jobs:
            job_id = job['id']
            slurm_job_id = job['slurm_job_id']
            job_type = job['job_type']
            
            # Check job status
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
                
                # Update process items
                self._update_completed_job_items(job_id, job_type)
                
            elif status in ["FAILED", "TIMEOUT", "CANCELLED"]:
                # Update job status
                self.db.update(
                    "ecod_schema.job",
                    {
                        "status": "failed"
                    },
                    "id = %s",
                    (job_id,)
                )
                
                # Update process items to failed
                query = """
                UPDATE ecod_schema.process_status
                SET status = 'error', 
                    current_stage = %s,
                    error_message = %s
                WHERE id IN (
                    SELECT process_id
                    FROM ecod_schema.job_item
                    WHERE job_id = %s
                )
                """
                
                error_stage = "profile_failed" if job_type == "hhblits" else "hhsearch_failed"
                error_msg = f"{job_type} job failed: {status}"
                
                self.db.execute_query(query, (error_stage, error_msg, job_id))
    
    def _update_completed_job_items(self, job_id: int, job_type: str):
        """Update items for a completed job"""
        # Get items in this job
        query = """
        SELECT 
            ji.id, ji.process_id, ps.relative_path, b.base_path,
            p.pdb_id, p.chain_id, p.source_id
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
        
        items = self.db.execute_dict_query(query, (job_id,))
        
        for item in items:
            process_id = item['process_id']
            relative_path = item['relative_path']
            base_path = item['base_path']
            pdb_id = item['pdb_id']
            chain_id = item['chain_id']
            
            # Define expected output file
            chain_dir = os.path.join(base_path, "ecod_dump", relative_path)
            
            if job_type == "hhblits":
                output_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.a3m")
                file_type = "a3m"
                next_stage = "profile_complete"
                relative_path = f"ecod_dump/{relative_path}/{pdb_id}_{chain_id}.a3m"
            else:  # hhsearch
                output_file = os.path.join(chain_dir, f"{pdb_id}_{chain_id}.{self.current_version}.hhr")
                file_type = "hhr"
                next_stage = "hhsearch_complete"
                relative_path = f"ecod_dump/{relative_path}/{pdb_id}_{chain_id}.{self.current_version}.hhr"
            
            # Check if file exists
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                # Register file
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": relative_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(output_file)
                    }
                )
                
                # Update process status
                self.db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": next_stage,
                        "status": "success"
                    },
                    "id = %s",
                    (process_id,)
                )
            else:
                # Update process status to error
                error_stage = "profile_failed" if job_type == "hhblits" else "hhsearch_failed"
                self.db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": error_stage,
                        "status": "error",
                        "error_message": f"Output file not found or empty: {output_file}"
                    },
                    "id = %s",
                    (process_id,)
                )
    
    def parse_hhr_file(self, hhr_file: str) -> List[Dict[str, Any]]:
        """Parse HHSearch result file (HHR format)"""
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
    
    def generate_summary_xml(self, process_id: int):
        """Generate XML summary of HHsearch results for a process"""
        # Get HHR file
        query = """
        SELECT 
            pf.file_path, b.base_path, p.pdb_id, p.chain_id
        FROM 
            ecod_schema.process_file pf
        JOIN
            ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE 
            pf.process_id = %s
            AND pf.file_type = 'hhr'
        """
        
        files = self.db.execute_dict_query(query, (process_id,))
        
        if not files:
            self.logger.warning(f"No HHR file found for process {process_id}")
            return False
        
        file_info = files[0]
        file_path = os.path.join(file_info['base_path'], file_info['file_path'])
        pdb_id = file_info['pdb_id']
        chain_id = file_info['chain_id']
        
        if not os.path.exists(file_path):
            self.logger.warning(f"HHR file not found: {file_path}")
            return False
        
        # Parse HHR file
        hits = self.parse_hhr_file(file_path)
        
        if not hits:
            self.logger.warning(f"No hits found in HHR file: {file_path}")
            return False
        
        # Create summary directory
        summary_dir = os.path.join(os.path.dirname(os.path.dirname(file_path)), "hh_summaries")
        os.makedirs(summary_dir, exist_ok=True)
        
        # Create XML structure
        root = ET.Element("hh_summ_doc")
        
        # Add metadata
        metadata = ET.SubElement(root, "metadata")
        ET.SubElement(metadata, "pdb_id").text = pdb_id
        ET.SubElement(metadata, "chain_id").text = chain_id
        ET.SubElement(metadata, "reference").text = self.current_version
        ET.SubElement(metadata, "creation_date").text = "CURRENT_TIMESTAMP"
        
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
                query_ranges = self._map_alignment_to_ranges(hit)
                
                # Add query range
                query_range_node = ET.SubElement(hit_node, "query_range")
                query_range_node.text = query_ranges
                
                # Add alignment
                alignment_node = ET.SubElement(hit_node, "alignment")
                ET.SubElement(alignment_node, "query_ali").text = hit.get('query_ali', '')
                ET.SubElement(alignment_node, "template_ali").text = hit.get('template_ali', '')
        
        # Write XML to file
        summary_file = os.path.join(summary_dir, f"{pdb_id}_{chain_id}_hh_summary.xml")
        
        tree = ET.ElementTree(root)
        tree.write(summary_file, encoding='utf-8', xml_declaration=True)
        
        # Register file
        self.db.insert(
            "ecod_schema.process_file",
            {
                "process_id": process_id,
                "file_type": "hh_summ_xml",
                "file_path": f"hh_summaries/{pdb_id}_{chain_id}_hh_summary.xml",
                "file_exists": True,
                "file_size": os.path.getsize(summary_file)
            }
        )
        
        self.logger.info(f"Generated HHsearch summary for process {process_id}")
        return True
    
    def _map_alignment_to_ranges(self, hit: Dict[str, Any]) -> str:
        """Map alignment to sequence ranges"""
        query_ali = hit.get('query_ali', '')
        query_start = hit.get('query_start', 1)
        
        # Calculate sequence positions
        positions = []
        pos = query_start
        
        for char in query_ali:
            if char != '-':
                positions.append(pos)
                pos += 1
            else:
                pos += 1
        
        # Convert to ranges
        if not positions:
            return ""
        
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

def main():
    parser = argparse.ArgumentParser(description='HHSearch Runner for PyECOD')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Specific batch ID to process')
    parser.add_argument('--hhblits', action='store_true',
                      help='Run HHblits profile generation')
    parser.add_argument('--hhsearch', action='store_true',
                      help='Run HHsearch against ECOD database')
    parser.add_argument('--check', action='store_true',
                      help='Check status of running jobs')
    parser.add_argument('--parse', action='store_true',
                      help='Parse HHsearch results')
    parser.add_argument('--process-id', type=int,
                      help='Specific process ID to parse results for')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of CPU threads to use')
    parser.add_argument('--limit', type=int, default=10,
                      help='Maximum number of items to process')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    runner = HHSearchRunner(args.config)
    
    if args.check:
        # Check status of running jobs
        runner.check_running_jobs()
        return
    
    if args.parse and args.process_id:
        # Parse results for specific process
        runner.generate_summary_xml(args.process_id)
        return
    
    # Get pending batches
    if args.batch_id:
        batches = runner.get_pending_batches()
        batch = next((b for b in batches if b['id'] == args.batch_id), None)
        if not batch:
            print(f"Batch ID {args.batch_id} not found or not pending")
            return
        batches = [batch]
    else:
        batches = runner.get_pending_batches()
    
    if not batches:
        print("No pending batches found")
        return
    
    # Process batches
    for batch in batches:
        batch_id = batch['id']
        print(f"Processing batch {batch_id}: {batch['batch_name']}")
        
        # Run HHblits
        if args.hhblits:
            # Get items ready for HHblits
            items = runner.get_hhblits_ready_items(batch_id, args.limit)
            
            if not items:
                print(f"No items ready for HHblits in batch {batch_id}")
            else:
                print(f"Found {len(items)} items ready for HHblits in batch {batch_id}")
                
                # Prepare and submit jobs
                hhblits_jobs = runner.prepare_hhblits_jobs(items, args.threads, 1)
                job_ids = runner.submit_jobs(hhblits_jobs, "hhblits")
                print(f"Submitted {len(job_ids)} HHblits jobs")
        
        # Run HHsearch
        if args.hhsearch:
            # Get items ready for HHsearch
            items = runner.get_hhsearch_ready_items(batch_id, args.limit)
            
            if not items:
                print(f"No items ready for HHsearch in batch {batch_id}")
            else:
                print(f"Found {len(items)} items ready for HHsearch in batch {batch_id}")
                
                # Prepare and submit jobs
                hhsearch_jobs = runner.prepare_hhsearch_jobs(items, args.threads, 1)
                job_ids = runner.submit_jobs(hhsearch_jobs, "hhsearch")
                print(f"Submitted {len(job_ids)} HHsearch jobs")
        
        # If neither specified, show help
        if not args.hhblits and not args.hhsearch and not args.check and not args.parse:
            print("Please specify --hhblits, --hhsearch, --check, or --parse")
            parser.print_help()

if __name__ == "__main__":
    main()