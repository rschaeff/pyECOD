#!/usr/bin/env python3
"""
Script to identify and fix batch proteins with missing chain BLAST results 
but existing domain BLAST results, and analyze regeneration issues.
"""

import os
import sys
import logging
from typing import Dict, List, Tuple, Set, Optional
import argparse
from datetime import datetime
import subprocess
import yaml

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

# Add current directory to Python path to ensure ecod modules can be imported
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)
    logger.info(f"Added {parent_dir} to Python path")

# Import ecod modules
try:
    from ecod.core.context import ApplicationContext
    from ecod.pipelines.blast_pipeline import BlastPipeline
    from ecod.exceptions import PipelineError, ConfigurationError
    from ecod.db import DBManager
    logger.info("Successfully imported ecod modules")
except ImportError as e:
    logger.error(f"Error importing ecod modules: {e}")
    logger.error("Make sure you're running this script from the pyecod project directory")
    sys.exit(1)


def load_configuration(config_path=None):
    """
    Load configuration from YAML file, merging with local configuration for secrets
    
    Args:
        config_path: Path to main configuration file
        
    Returns:
        Merged configuration dictionary
    """
    # If no config path provided, look in standard locations
    if not config_path:
        # Try standard config locations
        locations = [
            os.path.join(parent_dir, 'config', 'config.yml'),
            os.path.join(parent_dir, 'config', 'config.yaml'),
            os.path.join(parent_dir, 'config.yml'),
            os.path.join(parent_dir, 'config.yaml'),
            '/config/config.yml',
            '/config/config.yaml',
            os.path.expanduser('~/.ecod/config.yml')
        ]
        
        for loc in locations:
            if os.path.exists(loc):
                config_path = loc
                logger.info(f"Using main configuration from {loc}")
                break
    
    if not config_path or not os.path.exists(config_path):
        logger.error("No valid main configuration file found")
        return None
    
    # Load main config
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Find and load local config with secrets
        config_dir = os.path.dirname(config_path)
        local_config_paths = [
            os.path.join(config_dir, 'config.local.yml'),
            os.path.join(config_dir, 'config.local.yaml'),
            os.path.join(parent_dir, 'config', 'config.local.yml'),
            os.path.join(parent_dir, 'config', 'config.local.yaml')
        ]
        
        # Try to load local config with secrets
        local_config = {}
        for local_path in local_config_paths:
            if os.path.exists(local_path):
                logger.info(f"Loading local configuration from {local_path}")
                with open(local_path, 'r') as f:
                    local_config = yaml.safe_load(f)
                break
        
        # Deep merge configs
        merged_config = deep_merge(config, local_config)
        
        # Ensure database section exists
        if 'database' not in merged_config:
            logger.warning("No database section in config, using defaults")
            merged_config['database'] = {
                'host': 'localhost',
                'port': 5432,
                'user': 'postgres',
                'password': '',
                'database': 'ecod'
            }
        
        return merged_config
        
    except Exception as e:
        logger.error(f"Error loading configuration from {config_path}: {e}")
        return None


def deep_merge(dict1, dict2):
    """
    Deep merge two dictionaries
    
    Args:
        dict1: First dictionary
        dict2: Second dictionary (values override dict1)
        
    Returns:
        Merged dictionary
    """
    if not isinstance(dict1, dict) or not isinstance(dict2, dict):
        return dict2
    
    merged = dict1.copy()
    
    for key, value in dict2.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
            
    return merged


def find_proteins_with_domain_blast_only(db: DBManager, batch_id: int) -> List[Dict]:
    """
    Find proteins that have domain BLAST results but no chain BLAST results
    in the database.
    
    Args:
        db: Database manager
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
        return db.execute_dict_query(query, (batch_id,))
    except Exception as e:
        logger.error(f"Error querying proteins with domain BLAST only: {e}")
        return []


def get_batch_info(db: DBManager, batch_id: int) -> Dict:
    """Get batch information from database"""
    query = """
    SELECT id, batch_name, base_path, ref_version, total_items, status
    FROM ecod_schema.batch
    WHERE id = %s
    """
    
    try:
        rows = db.execute_dict_query(query, (batch_id,))
        if rows:
            return rows[0]
    except Exception as e:
        logger.error(f"Error getting batch info: {e}")
    
    return {}


def find_potential_chain_blast_paths(base_path: str, pdb_id: str, chain_id: str, relative_path: str) -> List[str]:
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


def check_db_fs_sync(batch_id: int, proteins: List[Dict], base_path: str) -> Tuple[List[Dict], List[Dict], List[Dict]]:
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
        chain_blast_paths = find_potential_chain_blast_paths(base_path, pdb_id, chain_id, relative_path)
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


def get_fasta_files(db: DBManager, batch_id: int, process_ids: List[int]) -> Dict[int, str]:
    """
    Get FASTA files for specified processes.
    
    Args:
        db: Database manager
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
        rows = db.execute_dict_query(query, (batch_id, tuple(process_ids)))
        
        for row in rows:
            process_id = row["process_id"]
            file_path = row["file_path"]
            base_path = row["base_path"]
            
            full_path = os.path.join(base_path, file_path)
            if os.path.exists(full_path):
                fasta_files[process_id] = full_path
            else:
                logger.warning(f"FASTA file {full_path} does not exist")
        
        return fasta_files
        
    except Exception as e:
        logger.error(f"Error getting FASTA files: {e}")
        return {}


def register_chain_blast_file(db: DBManager, process_id: int, file_path: str) -> None:
    """Register chain BLAST file in database"""
    # Check if record already exists
    query = """
    SELECT id FROM ecod_schema.process_file
    WHERE process_id = %s AND file_type = 'chain_blast_result'
    """
    
    try:
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
            logger.debug(f"Updated existing chain_blast_result file record for process {process_id}")
        else:
            # Insert new record
            db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "chain_blast_result",
                    "file_path": file_path,
                    "file_exists": True,
                    "file_size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
                }
            )
            logger.debug(f"Created new chain_blast_result file record for process {process_id}")
            
        # Update process status
        db.update(
            "ecod_schema.process_status",
            {
                "current_stage": "chain_blast_complete",
                "status": "success"
            },
            "id = %s",
            (process_id,)
        )
        
    except Exception as e:
        logger.error(f"Error registering chain BLAST file for process {process_id}: {e}")


def run_single_chain_blast(config: Dict, batch_id: int, process_id: int, protein: Dict, fasta_path: str, db: DBManager) -> bool:
    """
    Run chain BLAST for a single protein.
    
    Args:
        config: Configuration dictionary
        batch_id: Batch ID
        process_id: Process ID
        protein: Protein dictionary
        fasta_path: Path to FASTA file
        db: Database manager
        
    Returns:
        True if successful
    """
    # Get batch path
    batch_info = get_batch_info(db, batch_id)
    base_path = batch_info.get("base_path")
    if not base_path:
        logger.error(f"Failed to get base path for batch {batch_id}")
        return False
    
    # Get chain BLAST database
    chain_db = config.get('reference', {}).get('chain_db')
    if not chain_db:
        logger.error("Chain BLAST database not configured")
        return False
    
    # Get BLAST path
    blast_path = config.get('tools', {}).get('blast_path', 'blastp')
    if not blast_path:
        logger.warning("BLAST path not configured, using default 'blastp'")
        blast_path = "blastp"
    
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
        logger.info(f"Chain BLAST result already exists for {pdb_id}_{chain_id}: {output_file}")
        
        # Register the file in database
        relative_output = os.path.relpath(output_file, base_path)
        register_chain_blast_file(db, process_id, relative_output)
        return True
    
    # Build BLAST command
    cmd = f"{blast_path} -query {fasta_path} -db {chain_db} " \
          f"-outfmt 5 -num_alignments 5000 -evalue 0.002 " \
          f"-out {output_file}"
    
    logger.info(f"Running chain BLAST for {pdb_id}_{chain_id}: {cmd}")
    
    # Execute BLAST command
    try:
        result = subprocess.run(cmd, shell=True, check=True, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Check if output file exists and is not empty
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            # Register file in database
            relative_output = os.path.relpath(output_file, base_path)
            register_chain_blast_file(db, process_id, relative_output)
            
            logger.info(f"Successfully generated chain BLAST for {pdb_id}_{chain_id}")
            return True
        else:
            logger.error(f"Chain BLAST completed but output file is empty or missing: {output_file}")
            return False
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running chain BLAST for {pdb_id}_{chain_id}: {e}")
        logger.error(f"STDOUT: {e.stdout.decode() if e.stdout else ''}")
        logger.error(f"STDERR: {e.stderr.decode() if e.stderr else ''}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error running chain BLAST for {pdb_id}_{chain_id}: {e}")
        return False


def regenerate_chain_blast(config: Dict, db: DBManager, batch_id: int, proteins: List[Dict]) -> Dict:
    """
    Attempt to regenerate chain BLAST results for proteins missing them.
    
    Args:
        config: Configuration dictionary
        db: Database manager
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
        fasta_files = get_fasta_files(db, batch_id, process_ids)
        
        # Check if we have all required FASTA files
        for protein in proteins:
            process_id = protein["process_id"]
            if process_id not in fasta_files or not fasta_files[process_id]:
                logger.warning(f"No FASTA file found for protein {protein['pdb_id']}_{protein['chain_id']} (process_id: {process_id})")
                failures.append({
                    "protein": protein,
                    "reason": "no_fasta_file"
                })
                continue
            
            # Attempt to run chain BLAST for this protein
            try:
                success = run_single_chain_blast(config, batch_id, process_id, protein, fasta_files[process_id], db)
                if success:
                    regenerated.append(protein)
                else:
                    failures.append({
                        "protein": protein,
                        "reason": "blast_failed"
                    })
            except Exception as e:
                logger.error(f"Error running chain BLAST for {protein['pdb_id']}_{protein['chain_id']}: {e}")
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
        logger.error(f"Error in chain BLAST regeneration: {e}")
        return {
            "regenerated": 0,
            "failed": len(proteins),
            "error": str(e),
            "failures": [{"protein": p, "reason": "general_error"} for p in proteins]
        }


def get_protein_sequence_length(db: DBManager, process_id: int) -> Optional[int]:
    """Get sequence length for a protein"""
    query = """
    SELECT p.length
    FROM ecod_schema.protein p
    JOIN ecod_schema.process_status ps ON p.id = ps.protein_id
    WHERE ps.id = %s
    """
    
    try:
        rows = db.execute_query(query, (process_id,))
        if rows:
            return rows[0][0]
    except Exception as e:
        logger.error(f"Error getting protein sequence length: {e}")
        
    return None


def analyze_regeneration_failures(db: DBManager, batch_id: int, failures: List[Dict]) -> Dict[str, int]:
    """
    Analyze why chain BLAST regeneration failed for some proteins.
    
    Args:
        db: Database manager
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
            seq_length = get_protein_sequence_length(db, process_id)
            
            if seq_length is not None:
                if seq_length < 30:  # Arbitrary threshold for short sequences
                    if "short_sequence" not in reason_counts:
                        reason_counts["short_sequence"] = 0
                        detailed_analysis["short_sequence"] = []
                    reason_counts["short_sequence"] += 1
                    detailed_analysis["short_sequence"].append(protein)
    
    # Save detailed analysis to file
    save_failure_analysis(batch_id, detailed_analysis)
    
    return reason_counts


def save_failure_analysis(batch_id: int, analysis: Dict) -> None:
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
                
        logger.info(f"Saved failure analysis to {filename}")
        
    except Exception as e:
        logger.error(f"Error saving failure analysis: {e}")


def update_missing_file_record(db: DBManager, process_id: int, file_type: str) -> None:
    """Update process_file record for missing file"""
    query = """
    SELECT id FROM ecod_schema.process_file
    WHERE process_id = %s AND file_type = %s
    """
    
    try:
        existing = db.execute_query(query, (process_id, file_type))
        
        if existing:
            # Update existing record to mark file as not existing
            db.update(
                "ecod_schema.process_file",
                {
                    "file_exists": False,
                    "file_size": 0
                },
                "id = %s",
                (existing[0][0],)
            )
            logger.debug(f"Updated {file_type} file record for process {process_id} to mark as missing")
    except Exception as e:
        logger.error(f"Error updating missing file record for process {process_id}: {e}")


def save_overall_results(results: Dict) -> None:
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
                
        logger.info(f"Saved overall results to {filename}")
        
    except Exception as e:
        logger.error(f"Error saving overall results: {e}")


def check_batch(config: Dict, db: DBManager, batch_id: int) -> Dict:
    """
    Check a specific batch for proteins with missing chain BLAST results
    but existing domain BLAST results.
    
    Args:
        config: Configuration dictionary
        db: Database manager
        batch_id: Batch ID to check
        
    Returns:
        Dictionary with results summary
    """
    logger.info(f"Checking batch {batch_id}")
    
    # Get batch information
    batch_info = get_batch_info(db, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return {"error": "Batch not found"}
    
    base_path = batch_info.get("base_path")
    logger.info(f"Batch path: {base_path}")
    
    # Identify proteins with domain BLAST but no chain BLAST
    mismatched_proteins = find_proteins_with_domain_blast_only(db, batch_id)
    
    logger.info(f"Found {len(mismatched_proteins)} proteins with domain BLAST but no chain BLAST")
    
    # Check database vs filesystem discrepancies
    db_only, fs_only, both_mismatch = check_db_fs_sync(batch_id, mismatched_proteins, base_path)
    
    logger.info(f"Found {len(db_only)} proteins with records in DB only")
    logger.info(f"Found {len(fs_only)} proteins with files on FS only")
    logger.info(f"Found {len(both_mismatch)} proteins with mismatches between DB and FS")
    
    # Attempt to regenerate missing chain BLAST results
    regeneration_results = {"regenerated": 0, "failed": 0, "failures": []}
    if mismatched_proteins:
        regeneration_results = regenerate_chain_blast(config, db, batch_id, mismatched_proteins)
        regenerated_count = regeneration_results.get("regenerated", 0)
        failed_count = regeneration_results.get("failed", 0)
        logger.info(f"Regenerated chain BLAST for {regenerated_count} proteins")
        logger.info(f"Failed to regenerate chain BLAST for {failed_count} proteins")
        
        # Analyze failures
        if failed_count > 0:
            failure_analysis = analyze_regeneration_failures(db, batch_id, regeneration_results.get("failures", []))
            logger.info("Failure analysis:")
            for reason, count in failure_analysis.items():
                logger.info(f"  {reason}: {count}")
    
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


def fix_db_file_sync(config: Dict, db: DBManager, batch_id: int) -> Dict:
    """
    Fix synchronization issues between database and filesystem.
    
    Args:
        config: Configuration dictionary
        db: Database manager
        batch_id: Batch ID to fix
        
    Returns:
        Dictionary with results summary
    """
    logger.info(f"Fixing DB-file sync issues for batch {batch_id}")
    
    # Get batch information
    batch_info = get_batch_info(db, batch_id)
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return {"error": "Batch not found"}
    
    base_path = batch_info.get("base_path")
    
    # First check for discrepancies
    check_results = check_batch(config, db, batch_id)
    
    # Process FS-only cases (files exist but not in DB)
    fs_only_count = len(check_results.get("fs_only_proteins", []))
    if fs_only_count > 0:
        logger.info(f"Registering {fs_only_count} chain BLAST files that exist on filesystem but not in DB")
        fs_only_proteins = check_results.get("fs_only_proteins", [])
        
        for protein in fs_only_proteins:
            process_id = protein["process_id"]
            pdb_id = protein["pdb_id"]
            chain_id = protein["chain_id"]
            
            # Find the chain BLAST file
            chain_blast_paths = find_potential_chain_blast_paths(base_path, pdb_id, chain_id, protein["relative_path"])
            for path in chain_blast_paths:
                if os.path.exists(path):
                    relative_path = os.path.relpath(path, base_path)
                    register_chain_blast_file(db, process_id, relative_path)
                    logger.info(f"Registered existing chain BLAST file for {pdb_id}_{chain_id}")
                    break
    
    # Process DB-only cases (missing files that are registered in DB)
    db_only_count = len(check_results.get("both_mismatch_proteins", []))
    if db_only_count > 0:
        logger.info(f"Updating DB records for {db_only_count} missing files")
        db_only_proteins = check_results.get("both_mismatch_proteins", [])
        
        for protein in db_only_proteins:
            process_id = protein["process_id"]
            
            # Update process_file record
            update_missing_file_record(db, process_id, "domain_blast_result")
            
            # If we find DB record for chain_blast_result, update that too
            update_missing_file_record(db, process_id, "chain_blast_result")
    
    # Regenerate missing chain BLAST files
    mismatched_count = check_results.get("mismatched_total", 0)
    if mismatched_count > 0:
        logger.info(f"Attempting to regenerate {mismatched_count} missing chain BLAST files")
        # We already attempted regeneration in check_batch
        regenerated = check_results.get("regenerated", 0)
        logger.info(f"Regenerated {regenerated} chain BLAST files")
    
    # Return updated summary
    return {
        "batch_id": batch_id,
        "fs_only_fixed": fs_only_count,
        "db_only_fixed": db_only_count,
        "regenerated": check_results.get("regenerated", 0),
        "failed": check_results.get("failed", 0),
        "total_fixed": fs_only_count + db_only_count + check_results.get("regenerated", 0)
    }


def run_batch_check(config: Dict, db: DBManager, batch_ids: List[int]) -> None:
    """
    Run checks for multiple batches.
    
    Args:
        config: Configuration dictionary
        db: Database manager
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
            logger.info(f"Processing batch {batch_id}")
            results = check_batch(config, db, batch_id)
            
            overall_results["total_mismatched"] += results.get("mismatched_total", 0)
            overall_results["total_regenerated"] += results.get("regenerated", 0)
            overall_results["total_failed"] += results.get("failed", 0)
            overall_results["batch_results"][batch_id] = results
            
            logger.info(f"Completed batch {batch_id}")
            
        except Exception as e:
            logger.error(f"Error processing batch {batch_id}: {e}")
            overall_results["batch_results"][batch_id] = {"error": str(e)}
    
    # Save overall results
    save_overall_results(overall_results)
    
    logger.info("Batch check completed")
    logger.info(f"Total mismatched proteins: {overall_results['total_mismatched']}")
    logger.info(f"Total regenerated: {overall_results['total_regenerated']}")
    logger.info(f"Total failed: {overall_results['total_failed']}")


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
        # Load configuration
        config = load_configuration(args.config)
        if not config:
            logger.error("Failed to load configuration")
            return 1
        
        # Initialize database connection
        db = DBManager(config.get('database', {}))
        logger.info("Database connection established")
        
        if args.batch:
            # Process single batch
            if args.fix:
                results = fix_db_file_sync(config, db, args.batch)
                logger.info(f"Fixed {results.get('total_fixed', 0)} issues in batch {args.batch}")
            else:
                results = check_batch(config, db, args.batch)
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
                    results = fix_db_file_sync(config, db, batch_id)
                    logger.info(f"Fixed {results.get('total_fixed', 0)} issues in batch {batch_id}")
            else:
                run_batch_check(config, db, batch_ids)
        
        logger.info("Operation completed successfully")
        return 0
        
    except Exception as e:
        logger.error(f"Unhandled error: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    exit(main())