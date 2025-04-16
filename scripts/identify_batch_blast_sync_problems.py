#!/usr/bin/env python3
"""
Script to identify and fix batch proteins with missing chain BLAST results 
but existing domain BLAST results, and analyze regeneration issues.
"""

import os
import sys
import logging
import yaml
from typing import Dict, List, Tuple, Set, Optional
import argparse
from datetime import datetime

# Add parent directory to Python path to access ecod packages
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, parent_dir)

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.exceptions import PipelineError, ConfigurationError
from ecod.pipelines.blast_pipeline import BlastPipeline

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"blast_sync_check_{datetime.now().strftime('%Y%m%d_%H%M')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("blast_sync_checker")

class BlastSyncChecker:
    """
    Class to check and fix synchronization issues between database and filesystem
    for BLAST results in ECOD pipeline.
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize with database connection and configuration"""
        try:
            # Load configuration
            self.config = self._load_config(config_path)
            
            # Set up database connection
            from ecod.db import DBManager
            self.db = DBManager(self.config.get('database', {}))
            self.job = 
            
            # Initialize blast pipeline with our config and db
            from ecod.pipelines.blast_pipeline import BlastPipeline
            self.blast_pipeline = BlastPipeline(self)
            
            self.logger = logger
            logger.info("BlastSyncChecker initialized successfully")
            
        except Exception as e:
            logger.error(f"Error initializing BlastSyncChecker: {e}")
            raise
    
    def _load_config(self, config_path: Optional[str] = None) -> Dict:
        """Load configuration from YAML file"""
        # If no config path provided, look in standard locations
        if not config_path:
            # Try standard config locations
            locations = [
                os.path.join(parent_dir, 'config', 'config.yml'),
                os.path.join(parent_dir, 'config.yml'),
                '/config/config.yml',
                os.path.expanduser('~/.ecod/config.yml')
            ]
            
            for loc in locations:
                if os.path.exists(loc):
                    config_path = loc
                    logger.info(f"Using configuration from {loc}")
                    break
        
        if not config_path or not os.path.exists(config_path):
            logger.error("No valid configuration file found")
            raise FileNotFoundError("No valid configuration file found")
        
        # Load YAML config
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
            
            # Ensure database section exists
            if 'database' not in config:
                logger.warning("No database section in config, using defaults")
                config['database'] = {
                    'host': 'localhost',
                    'port': 5432,
                    'user': 'postgres',
                    'password': '',
                    'database': 'ecod'
                }
            
            return config
            
        except Exception as e:
            logger.error(f"Error loading configuration from {config_path}: {e}")
            raise
    
    # This method allows us to be used as a context
    def is_force_overwrite(self):
        """Check if force overwrite is enabled in config"""
        return self.config.get('force_overwrite', False)
    
    def check_batch(self, batch_id: int) -> Dict:
        """
        Check a specific batch for proteins with missing chain BLAST results
        but existing domain BLAST results.
        
        Args:
            batch_id: Batch ID to check
            
        Returns:
            Dictionary with results summary
        """
        self.logger.info(f"Checking batch {batch_id}")
        
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return {"error": "Batch not found"}
        
        base_path = batch_info.get("base_path")
        self.logger.info(f"Batch path: {base_path}")
        
        # Identify proteins with domain BLAST but no chain BLAST
        mismatched_proteins = self._find_proteins_with_domain_blast_only(batch_id)
        
        self.logger.info(f"Found {len(mismatched_proteins)} proteins with domain BLAST but no chain BLAST")
        
        # Check database vs filesystem discrepancies
        db_only, fs_only, both_mismatch = self._check_db_fs_sync(batch_id, mismatched_proteins, base_path)
        
        self.logger.info(f"Found {len(db_only)} proteins with records in DB only")
        self.logger.info(f"Found {len(fs_only)} proteins with files on FS only")
        self.logger.info(f"Found {len(both_mismatch)} proteins with mismatches between DB and FS")
        
        # Attempt to regenerate missing chain BLAST results
        regeneration_results = {"regenerated": 0, "failed": 0, "failures": []}
        if mismatched_proteins:
            regeneration_results = self._regenerate_chain_blast(batch_id, mismatched_proteins)
            regenerated_count = regeneration_results.get("regenerated", 0)
            failed_count = regeneration_results.get("failed", 0)
            self.logger.info(f"Regenerated chain BLAST for {regenerated_count} proteins")
            self.logger.info(f"Failed to regenerate chain BLAST for {failed_count} proteins")
            
            # Analyze failures
            if failed_count > 0:
                failure_analysis = self._analyze_regeneration_failures(
                    batch_id, regeneration_results.get("failures", [])
                )
                self.logger.info("Failure analysis:")
                for reason, count in failure_analysis.items():
                    self.logger.info(f"  {reason}: {count}")
        
        # Return summary
        return {
            "batch_id": batch_id,
            "base_path": base_path,
            "mismatched_total": len(mismatched_proteins),
            "db_only": len(db_only),
            "fs_only": len(fs_only),
            "both_mismatch": len(both_mismatch),
            "regenerated": regeneration_results.get("regenerated", 0),
            "failed": regeneration_results.get("failed", 0),
            "failure_analysis": failure_analysis if "failure_analysis" in locals() else {},
            "mismatched_proteins": mismatched_proteins,
            "db_only_proteins": db_only,
            "fs_only_proteins": fs_only,
            "both_mismatch_proteins": both_mismatch
        }
    
    def _get_batch_info(self, batch_id: int) -> Dict:
        """Get batch information from database"""
        query = """
        SELECT id, batch_name, base_path, ref_version, total_items, status
        FROM ecod_schema.batch
        WHERE id = %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            if rows:
                return rows[0]
        except Exception as e:
            self.logger.error(f"Error getting batch info: {e}")
        
        return {}
    
    def _find_proteins_with_domain_blast_only(self, batch_id: int) -> List[Dict]:
        """
        Find proteins that have domain BLAST results but no chain BLAST results
        in the database.
        
        Args:
            batch_id: Batch ID to check
            
        Returns:
            List of protein dictionaries with domain BLAST but no chain BLAST
        """
        query = """
        WITH batch_proteins AS (
            SELECT 
                ps.id as process_id, 
                p.id as protein_id, 
                p.pdb_id, 
                p.chain_id, 
                ps.relative_path,
                p.source_id
            FROM 
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            WHERE 
                ps.batch_id = %s
        ),
        domain_blast_files AS (
            SELECT 
                bp.process_id,
                bp.protein_id,
                bp.pdb_id,
                bp.chain_id,
                bp.relative_path,
                bp.source_id,
                pf.file_path as domain_blast_path
            FROM 
                batch_proteins bp
            JOIN
                ecod_schema.process_file pf ON bp.process_id = pf.process_id
            WHERE 
                pf.file_type = 'domain_blast_result'
                AND pf.file_exists = TRUE
        ),
        chain_blast_files AS (
            SELECT 
                bp.process_id,
                pf.file_path as chain_blast_path
            FROM 
                batch_proteins bp
            JOIN
                ecod_schema.process_file pf ON bp.process_id = pf.process_id
            WHERE 
                pf.file_type = 'chain_blast_result'
                AND pf.file_exists = TRUE
        )
        SELECT 
            d.*
        FROM 
            domain_blast_files d
        LEFT JOIN
            chain_blast_files c ON d.process_id = c.process_id
        WHERE 
            c.process_id IS NULL
        """
        
        try:
            return self.db.execute_dict_query(query, (batch_id,))
        except Exception as e:
            self.logger.error(f"Error querying proteins with domain BLAST only: {e}")
            return []
    
    def _check_db_fs_sync(self, batch_id: int, proteins: List[Dict], base_path: str) -> Tuple[List[Dict], List[Dict], List[Dict]]:
        """
        Check for inconsistencies between database records and actual files
        on the filesystem.
        
        Args:
            batch_id: Batch ID to check
            proteins: List of protein dictionaries to check
            base_path: Base path of the batch
            
        Returns:
            Tuple of (db_only, fs_only, both_mismatch) proteins
        """
        db_only = []
        fs_only = []
        both_mismatch = []
        
        for protein in proteins:
            process_id = protein["process_id"]
            pdb_id = protein["pdb_id"]
            chain_id = protein["chain_id"]
            relative_path = protein["relative_path"]
            
            # Check domain BLAST file on filesystem
            domain_blast_path_db = protein.get("domain_blast_path", "")
            domain_blast_path_fs = os.path.join(base_path, domain_blast_path_db)
            domain_blast_exists_fs = os.path.exists(domain_blast_path_fs) if domain_blast_path_db else False
            
            # Check for chain BLAST file on filesystem (even though not in DB)
            chain_blast_paths = self._find_potential_chain_blast_paths(base_path, pdb_id, chain_id, relative_path)
            chain_blast_exists_fs = any(os.path.exists(path) for path in chain_blast_paths)
            
            if domain_blast_exists_fs and not chain_blast_exists_fs:
                # Domain BLAST exists in both DB and FS, but no chain BLAST anywhere
                db_only.append(protein)
            elif domain_blast_exists_fs and chain_blast_exists_fs:
                # Domain BLAST exists in both DB and FS, and chain BLAST exists on FS but not in DB
                fs_only.append(protein)
            elif not domain_blast_exists_fs:
                # Domain BLAST exists in DB but not on FS
                both_mismatch.append(protein)
            
        return db_only, fs_only, both_mismatch
    
    def _find_potential_chain_blast_paths(self, base_path: str, pdb_id: str, chain_id: str, relative_path: str) -> List[str]:
        """Find potential paths where chain BLAST files might exist"""
        # Common patterns for chain BLAST files
        patterns = [
            os.path.join(base_path, "chain_blast_results", f"{pdb_id}_{chain_id}.chainwise_blast.xml"),
            os.path.join(base_path, "chain_blast_results", f"{pdb_id.lower()}_{chain_id}.chainwise_blast.xml"),
            os.path.join(base_path, relative_path, f"{pdb_id}_{chain_id}.chainwise_blast.xml"),
            os.path.join(base_path, relative_path, "chain_blast_results", f"{pdb_id}_{chain_id}.chainwise_blast.xml")
        ]
        
        # Also check under 'ecod_dump' structure which is commonly used
        dump_path = os.path.join(base_path, "ecod_dump")
        if os.path.exists(dump_path):
            patterns.extend([
                os.path.join(dump_path, f"{pdb_id}_{chain_id}", f"{pdb_id}_{chain_id}.chainwise_blast.xml"),
                os.path.join(dump_path, f"{pdb_id}_{chain_id}", "chain_blast_results", f"{pdb_id}_{chain_id}.chainwise_blast.xml")
            ])
        
        # Check if any directories exist that might contain chain BLAST results
        chain_results_dir = os.path.join(base_path, "chain_blast_results")
        if os.path.exists(chain_results_dir):
            # Look for any files matching pattern for this protein
            for filename in os.listdir(chain_results_dir):
                if f"{pdb_id}_{chain_id}" in filename or f"{pdb_id.lower()}_{chain_id}" in filename:
                    patterns.append(os.path.join(chain_results_dir, filename))
        
        return patterns
    
    def _regenerate_chain_blast(self, batch_id: int, proteins: List[Dict]) -> Dict:
        """
        Attempt to regenerate chain BLAST results for proteins missing them.
        
        Args:
            batch_id: Batch ID
            proteins: List of protein dictionaries to regenerate chain BLAST for
            
        Returns:
            Dictionary with regeneration results
        """
        regenerated = []
        failures = []
        
        # Get process IDs for these proteins
        process_ids = [protein["process_id"] for protein in proteins]
        
        try:
            # Get FASTA files for these proteins
            fasta_files = self._get_fasta_files(batch_id, process_ids)
            
            # Check if we have all required FASTA files
            for protein in proteins:
                process_id = protein["process_id"]
                if process_id not in fasta_files or not fasta_files[process_id]:
                    self.logger.warning(f"No FASTA file found for protein {protein['pdb_id']}_{protein['chain_id']} (process_id: {process_id})")
                    failures.append({
                        "protein": protein,
                        "reason": "no_fasta_file"
                    })
                    continue
                
                # Attempt to run chain BLAST for this protein
                try:
                    success = self._run_single_chain_blast(batch_id, process_id, protein, fasta_files[process_id])
                    if success:
                        regenerated.append(protein)
                    else:
                        failures.append({
                            "protein": protein,
                            "reason": "blast_failed"
                        })
                except Exception as e:
                    self.logger.error(f"Error running chain BLAST for {protein['pdb_id']}_{protein['chain_id']}: {e}")
                    failures.append({
                        "protein": protein,
                        "reason": "exception",
                        "error": str(e)
                    })
            
            return {
                "regenerated": len(regenerated),
                "failed": len(failures),
                "regenerated_proteins": regenerated,
                "failures": failures
            }
            
        except Exception as e:
            self.logger.error(f"Error in chain BLAST regeneration: {e}")
            return {
                "regenerated": 0,
                "failed": len(proteins),
                "error": str(e),
                "failures": [{"protein": p, "reason": "general_error"} for p in proteins]
            }
    
    def _get_fasta_files(self, batch_id: int, process_ids: List[int]) -> Dict[int, str]:
        """
        Get FASTA files for specified processes.
        
        Args:
            batch_id: Batch ID
            process_ids: List of process IDs
            
        Returns:
            Dictionary mapping process_id to FASTA file path
        """
        query = """
        SELECT 
            pf.process_id, 
            pf.file_path,
            b.base_path
        FROM 
            ecod_schema.process_file pf
        JOIN
            ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            ps.batch_id = %s
            AND pf.process_id IN %s
            AND pf.file_type = 'fasta'
            AND pf.file_exists = TRUE
        """
        
        fasta_files = {}
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id, tuple(process_ids)))
            
            for row in rows:
                process_id = row["process_id"]
                file_path = row["file_path"]
                base_path = row["base_path"]
                
                full_path = os.path.join(base_path, file_path)
                if os.path.exists(full_path):
                    fasta_files[process_id] = full_path
                else:
                    self.logger.warning(f"FASTA file {full_path} does not exist")
            
            return fasta_files
            
        except Exception as e:
            self.logger.error(f"Error getting FASTA files: {e}")
            return {}
    
    def _run_single_chain_blast(self, batch_id: int, process_id: int, protein: Dict, fasta_path: str) -> bool:
        """
        Run chain BLAST for a single protein.
        
        Args:
            batch_id: Batch ID
            process_id: Process ID
            protein: Protein dictionary
            fasta_path: Path to FASTA file
            
        Returns:
            True if successful
        """
        # Get batch path
        batch_info = self._get_batch_info(batch_id)
        base_path = batch_info.get("base_path")
        if not base_path:
            self.logger.error(f"Failed to get base path for batch {batch_id}")
            return False
        
        # Get chain BLAST database
        chain_db = self.config.get('reference', {}).get('chain_db')
        if not chain_db:
            self.logger.error("Chain BLAST database not configured")
            return False
        
        # Get BLAST path
        blast_path = self.config.get('tools', {}).get('blast_path', 'blastp')
        
        # Create chain blast directory
        chain_blast_dir = os.path.join(base_path, "chain_blast_results")
        os.makedirs(chain_blast_dir, exist_ok=True)
        
        # Define output file
        pdb_id = protein["pdb_id"]
        chain_id = protein["chain_id"]
        source_id = protein.get("source_id", f"{pdb_id}_{chain_id}")
        output_file = os.path.join(chain_blast_dir, f"{source_id}.chainwise_blast.xml")
        
        # Check if output already exists
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            self.logger.info(f"Chain BLAST result already exists for {pdb_id}_{chain_id}: {output_file}")
            
            # Register the file in database
            relative_output = os.path.relpath(output_file, base_path)
            self._register_chain_blast_file(process_id, relative_output)
            return True
        
        # Build BLAST command
        cmd = f"{blast_path} -query {fasta_path} -db {chain_db} " \
              f"-outfmt 5 -num_alignments 5000 -evalue 0.002 " \
              f"-out {output_file}"
        
        self.logger.info(f"Running chain BLAST for {pdb_id}_{chain_id}: {cmd}")
        
        # Execute BLAST command
        import subprocess
        try:
            result = subprocess.run(cmd, shell=True, check=True, 
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Check if output file exists and is not empty
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                # Register file in database
                relative_output = os.path.relpath(output_file, base_path)
                self._register_chain_blast_file(process_id, relative_output)
                
                self.logger.info(f"Successfully generated chain BLAST for {pdb_id}_{chain_id}")
                return True
            else:
                self.logger.error(f"Chain BLAST completed but output file is empty or missing: {output_file}")
                return False
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error running chain BLAST for {pdb_id}_{chain_id}: {e}")
            self.logger.error(f"STDOUT: {e.stdout.decode() if e.stdout else ''}")
            self.logger.error(f"STDERR: {e.stderr.decode() if e.stderr else ''}")
            return False
        except Exception as e:
            self.logger.error(f"Unexpected error running chain BLAST for {pdb_id}_{chain_id}: {e}")
            return False
    
    def _register_chain_blast_file(self, process_id: int, file_path: str) -> None:
        """Register chain BLAST file in database"""
        # Check if record already exists
        query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'chain_blast_result'
        """
        
        try:
            existing = self.db.execute_query(query, (process_id,))
            
            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": file_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    },
                    "id = %s",
                    (existing[0][0],)
                )
                self.logger.debug(f"Updated existing chain_blast_result file record for process {process_id}")
            else:
                # Insert new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "chain_blast_result",
                        "file_path": file_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    }
                )
                self.logger.debug(f"Created new chain_blast_result file record for process {process_id}")
                
            # Update process status
            self.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": "chain_blast_complete",
                    "status": "success"
                },
                "id = %s",
                (process_id,)
            )
            
        except Exception as e:
            self.logger.error(f"Error registering chain BLAST file for process {process_id}: {e}")
    
    def _analyze_regeneration_failures(self, batch_id: int, failures: List[Dict]) -> Dict[str, int]:
        """
        Analyze why chain BLAST regeneration failed for some proteins.
        
        Args:
            batch_id: Batch ID
            failures: List of failure dictionaries
            
        Returns:
            Dictionary mapping failure reasons to counts
        """
        reason_counts = {}
        detailed_analysis = {}
        
        for failure in failures:
            reason = failure.get("reason", "unknown")
            
            if reason not in reason_counts:
                reason_counts[reason] = 0
                detailed_analysis[reason] = []
                
            reason_counts[reason] += 1
            detailed_analysis[reason].append(failure.get("protein", {}))
            
            # Perform additional analysis for specific failure reasons
            if reason == "blast_failed":
                protein = failure.get("protein", {})
                pdb_id = protein.get("pdb_id", "")
                chain_id = protein.get("chain_id", "")
                
                # Check for short sequences
                process_id = protein.get("process_id")
                seq_length = self._get_protein_sequence_length(process_id)
                
                if seq_length is not None:
                    if seq_length < 30:  # Arbitrary threshold for short sequences
                        if "short_sequence" not in reason_counts:
                            reason_counts["short_sequence"] = 0
                            detailed_analysis["short_sequence"] = []
                        reason_counts["short_sequence"] += 1
                        detailed_analysis["short_sequence"].append(protein)
        
        # Save detailed analysis to file
        self._save_failure_analysis(batch_id, detailed_analysis)
        
        return reason_counts
    
    def _get_protein_sequence_length(self, process_id: int) -> Optional[int]:
        """Get sequence length for a protein"""
        query = """
        SELECT p.length
        FROM ecod_schema.protein p
        JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
        WHERE ps.id = %s
        """
        
        try:
            rows = self.db.execute_query(query, (process_id,))
            if rows:
                return rows[0][0]
        except Exception as e:
            self.logger.error(f"Error getting protein sequence length: {e}")
            
        return None
    
    def _save_failure_analysis(self, batch_id: int, analysis: Dict) -> None:
        """Save detailed failure analysis to file"""
        filename = f"batch_{batch_id}_failure_analysis_{datetime.now().strftime('%Y%m%d_%H%M')}.txt"
        
        try:
            with open(filename, 'w') as f:
                f.write(f"Failure Analysis for Batch {batch_id}\n")
                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                for reason, proteins in analysis.items():
                    f.write(f"=== {reason.upper()} ({len(proteins)} proteins) ===\n")
                    for protein in proteins:
                        pdb_id = protein.get("pdb_id", "")
                        chain_id = protein.get("chain_id", "")
                        process_id = protein.get("process_id", "")
                        f.write(f"  - {pdb_id}_{chain_id} (Process ID: {process_id})\n")
                    f.write("\n")
                    
            self.logger.info(f"Saved failure analysis to {filename}")
            
        except Exception as e:
            self.logger.error(f"Error saving failure analysis: {e}")
    
    def fix_db_file_sync(self, batch_id: int) -> Dict:
        """
        Fix synchronization issues between database and filesystem.
        
        Args:
            batch_id: Batch ID to fix
            
        Returns:
            Dictionary with results summary
        """
        self.logger.info(f"Fixing DB-file sync issues for batch {batch_id}")
        
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return {"error": "Batch not found"}
        
        base_path = batch_info.get("base_path")
        
        # First check for discrepancies
        check_results = self.check_batch(batch_id)
        
        # Process FS-only cases (files exist but not in DB)
        fs_only_count = len(check_results.get("fs_only_proteins", []))
        if fs_only_count > 0:
            self.logger.info(f"Registering {fs_only_count} chain BLAST files that exist on filesystem but not in DB")
            fs_only_proteins = check_results.get("fs_only_proteins", [])
            
            for protein in fs_only_proteins:
                process_id = protein["process_id"]
                pdb_id = protein["pdb_id"]
                chain_id = protein["chain_id"]
                
                # Find the chain BLAST file
                chain_blast_paths = self._find_potential_chain_blast_paths(base_path, pdb_id, chain_id, protein["relative_path"])
                for path in chain_blast_paths:
                    if os.path.exists(path):
                        relative_path = os.path.relpath(path, base_path)
                        self._register_chain_blast_file(process_id, relative_path)
                        self.logger.info(f"Registered existing chain BLAST file for {pdb_id}_{chain_id}")
                        break
        
        # Process DB-only cases (missing files that are registered in DB)
        db_only_count = len(check_results.get("both_mismatch_proteins", []))
        if db_only_count > 0:
            self.logger.info(f"Updating DB records for {db_only_count} missing files")
            db_only_proteins = check_results.get("both_mismatch_proteins", [])
            
            for protein in db_only_proteins:
                process_id = protein["process_id"]
                
                # Update process_file record
                self._update_missing_file_record(process_id, "domain_blast_result")
                
                # If we find DB record for chain_blast_result, update that too
                self._update_missing_file_record(process_id, "chain_blast_result")
        
        # Regenerate missing chain BLAST files
        mismatched_count = check_results.get("mismatched_total", 0)
        if mismatched_count > 0:
            self.logger.info(f"Attempting to regenerate {mismatched_count} missing chain BLAST files")
            # We already attempted regeneration in check_batch
            regenerated = check_results.get("regenerated", 0)
            self.logger.info(f"Regenerated {regenerated} chain BLAST files")
        
        # Return updated summary
        return {
            "batch_id": batch_id,
            "fs_only_fixed": fs_only_count,
            "db_only_fixed": db_only_count,
            "regenerated": check_results.get("regenerated", 0),
            "failed": check_results.get("failed", 0),
            "total_fixed": fs_only_count + db_only_count + check_results.get("regenerated", 0)
        }
    
    def _update_missing_file_record(self, process_id: int, file_type: str) -> None:
        """Update process_file record for missing file"""
        query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = %s
        """
        
        try:
            existing = self.db.execute_query(query, (process_id, file_type))
            
            if existing:
                # Update existing record to mark file as not existing
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_exists": False,
                        "file_size": 0
                    },
                    "id = %s",
                    (existing[0][0],)
                )
                self.logger.debug(f"Updated {file_type} file record for process {process_id} to mark as missing")
        except Exception as e:
            self.logger.error(f"Error updating missing file record for process {process_id}: {e}")
    
    def run_batch_check(self, batch_ids: List[int]) -> None:
        """
        Run checks for multiple batches.
        
        Args:
            batch_ids: List of batch IDs to check
        """
        overall_results = {
            "batches_checked": len(batch_ids),
            "total_mismatched": 0,
            "total_regenerated": 0,
            "total_failed": 0,
            "batch_results": {}
        }
        
        for batch_id in batch_ids:
            try:
                self.logger.info(f"Processing batch {batch_id}")
                results = self.check_batch(batch_id)
                
                overall_results["total_mismatched"] += results.get("mismatched_total", 0)
                overall_results["total_regenerated"] += results.get("regenerated", 0)
                overall_results["total_failed"] += results.get("failed", 0)
                overall_results["batch_results"][batch_id] = results
                
                self.logger.info(f"Completed batch {batch_id}")
                
            except Exception as e:
                self.logger.error(f"Error processing batch {batch_id}: {e}")
                overall_results["batch_results"][batch_id] = {"error": str(e)}
        
        # Save overall results
        self._save_overall_results(overall_results)
        
        self.logger.info("Batch check completed")
        self.logger.info(f"Total mismatched proteins: {overall_results['total_mismatched']}")
        self.logger.info(f"Total regenerated: {overall_results['total_regenerated']}")
        self.logger.info(f"Total failed: {overall_results['total_failed']}")
    
    def _save_overall_results(self, results: Dict) -> None:
        """Save overall check results to file"""
        filename = f"batch_check_results_{datetime.now().strftime('%Y%m%d_%H%M')}.txt"
        
        try:
            with open(filename, 'w') as f:
                f.write(f"ECOD BLAST Sync Check Results\n")
                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write(f"Batches checked: {results['batches_checked']}\n")
                f.write(f"Total mismatched proteins: {results['total_mismatched']}\n")
                f.write(f"Total regenerated: {results['total_regenerated']}\n")
                f.write(f"Total failed: {results['total_failed']}\n\n")
                
                f.write("Batch Details:\n")
                f.write("=============\n\n")
                
                for batch_id, batch_results in results["batch_results"].items():
                    f.write(f"Batch {batch_id}:\n")
                    
                    if "error" in batch_results:
                        f.write(f"  ERROR: {batch_results['error']}\n")
                        continue
                    
                    f.write(f"  Mismatched proteins: {batch_results.get('mismatched_total', 0)}\n")
                    
                    if "db_only" in batch_results:
                        f.write(f"  DB only (no files): {batch_results['db_only']}\n")
                    
                    if "fs_only" in batch_results:
                        f.write(f"  FS only (not in DB): {batch_results['fs_only']}\n")
                    
                    if "both_mismatch" in batch_results:
                        f.write(f"  Both mismatch: {batch_results['both_mismatch']}\n")
                    
                    f.write(f"  Regenerated: {batch_results.get('regenerated', 0)}\n")
                    f.write(f"  Failed: {batch_results.get('failed', 0)}\n")
                    
                    if "failure_analysis" in batch_results and batch_results["failure_analysis"]:
                        f.write(f"  Failure reasons:\n")
                        for reason, count in batch_results["failure_analysis"].items():
                            f.write(f"    - {reason}: {count}\n")
                    
                    f.write("\n")
                    
            self.logger.info(f"Saved overall results to {filename}")
            
        except Exception as e:
            self.logger.error(f"Error saving overall results: {e}")

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Check and fix BLAST synchronization issues in ECOD pipeline.')
    parser.add_argument('--batch', type=int, help='Batch ID to check', required=False)
    parser.add_argument('--batch-list', type=str, help='File with list of batch IDs to check', required=False)
    parser.add_argument('--config', type=str, help='Path to configuration file', required=False)
    parser.add_argument('--fix', action='store_true', help='Attempt to fix synchronization issues')
    
    args = parser.parse_args()
    
    if not args.batch and not args.batch_list:
        parser.error("At least one of --batch or --batch-list is required")
    
    try:
        # Initialize checker
        checker = BlastSyncChecker(args.config)
        
        if args.batch:
            # Process single batch
            if args.fix:
                results = checker.fix_db_file_sync(args.batch)
                logger.info(f"Fixed {results.get('total_fixed', 0)} issues in batch {args.batch}")
            else:
                results = checker.check_batch(args.batch)
                logger.info(f"Found {results.get('mismatched_total', 0)} issues in batch {args.batch}")
        
        elif args.batch_list:
            # Process multiple batches
            batch_ids = []
            try:
                with open(args.batch_list, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and line.isdigit():
                            batch_ids.append(int(line))
            except Exception as e:
                logger.error(f"Error reading batch list file: {e}")
                return 1
            
            if not batch_ids:
                logger.error("No valid batch IDs found in batch list file")
                return 1
            
            if args.fix:
                for batch_id in batch_ids:
                    results = checker.fix_db_file_sync(batch_id)
                    logger.info(f"Fixed {results.get('total_fixed', 0)} issues in batch {batch_id}")
            else:
                checker.run_batch_check(batch_ids)
        
        logger.info("Operation completed successfully")
        return 0
        
    except Exception as e:
        logger.error(f"Unhandled error: {e}", exc_info=True)
        return 1

if __name__ == "__main__":
    exit(main())