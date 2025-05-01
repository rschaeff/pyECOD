#!/usr/bin/env python3
"""
Enhanced run_domain_partition.py - Run domain partition with improved path handling
"""

import os
import sys
import logging
import argparse
from pathlib import Path

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.pipelines.domain_analysis.pipeline import DomainAnalysisPipeline
from ecod.utils.path_utils import get_standardized_paths, get_all_evidence_paths

def setup_logging(verbose=False, log_file=None):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO

    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )

class EnhancedPartitionRunner:
    """Enhanced runner for domain partition with improved path handling"""

    def __init__(self, config_path=None):
        """Initialize with configuration"""
        self.logger = logging.getLogger("ecod.partition_runner")

        # Load configuration
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config

        # Initialize database connection
        db_config = self.config_manager.get_db_config()
        self.db = DBManager(db_config)

        # Needed for DomainAnalysisPipeline
        self.db_manager = self.db

        # Initialize domain analysis pipeline
        self.domain_pipeline = DomainAnalysisPipeline(self)

    def get_db(self):
        """Access method for database - mimics ApplicationContext interface"""
        return self.db

    def is_force_overwrite(self):
        """Check if force overwrite is enabled - mimics ApplicationContext interface"""
        return self.config.get('pipeline', {}).get('force_overwrite', False)

    def get_process_details(self, process_id):
        """Get detailed information about a specific process"""
        query = """
        SELECT
            ps.id as process_id,
            p.pdb_id,
            p.chain_id,
            ps.batch_id,
            ps.relative_path,
            ps.current_stage,
            ps.status,
            b.base_path,
            b.ref_version
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE
            ps.id = %s
        """

        results = self.db.execute_dict_query(query, (process_id,))
        if not results:
            return None
        return results[0]

    def verify_domain_summary_exists(self, pdb_id, chain_id, batch_path, ref_version, blast_only=False):
        """Check if domain summary file exists (standard or legacy paths)"""
        self.logger.info(f"Checking for domain summary: {pdb_id}_{chain_id}")

        # Use path_utils to check all standard and legacy paths
        summary_type = 'blast_only_summary' if blast_only else 'domain_summary'
        evidence_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, ref_version)

        if summary_type in evidence_paths and evidence_paths[summary_type]['exists_at']:
            path = evidence_paths[summary_type]['exists_at']
            self.logger.info(f"Found domain summary at: {path}")
            return path

        # Query database as fallback
        query = """
        SELECT pf.file_path, pf.file_exists
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        AND pf.file_type = %s
        LIMIT 1
        """

        results = self.db.execute_dict_query(query, (pdb_id, chain_id, summary_type))
        if results and results[0]['file_exists']:
            db_path = results[0]['file_path']
            full_path = os.path.join(batch_path, db_path) if not os.path.isabs(db_path) else db_path
            if os.path.exists(full_path):
                self.logger.info(f"Found domain summary via DB: {full_path}")
                return full_path

        self.logger.warning(f"No domain summary found for {pdb_id}_{chain_id}")
        return None

    def check_paths_for_process(self, process_id):
        """Check paths for a specific process with detailed diagnostics"""
        process = self.get_process_details(process_id)
        if not process:
            self.logger.error(f"Process {process_id} not found")
            return False

        pdb_id = process['pdb_id']
        chain_id = process['chain_id']
        batch_path = process['base_path']
        ref_version = process['ref_version']

        self.logger.info(f"Checking paths for {pdb_id}_{chain_id}")
        self.logger.info(f"Batch path: {batch_path}")
        self.logger.info(f"Reference version: {ref_version}")

        # Use path_utils to get all standard paths
        standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)
        self.logger.info("Standard paths:")
        for file_type, path in standard_paths.items():
            exists = os.path.exists(path)
            status = "EXISTS" if exists else "MISSING"
            self.logger.info(f"  {file_type:20s}: {status:8s} {path}")

        # Check for domain summary specifically
        summary_path = self.verify_domain_summary_exists(pdb_id, chain_id, batch_path, ref_version)
        if summary_path:
            self.logger.info(f"Domain summary confirmed at: {summary_path}")
        else:
            self.logger.error(f"MISSING DOMAIN SUMMARY for {pdb_id}_{chain_id}")

            # Check for related files to help diagnose the issue
            self.logger.info("Checking related files:")

            # Check FASTA
            fasta_path = standard_paths.get('fasta')
            if fasta_path and os.path.exists(fasta_path):
                self.logger.info(f"FASTA file exists: {fasta_path}")
            else:
                self.logger.warning(f"FASTA file missing: {fasta_path}")

            # Check BLAST results
            chain_blast = standard_paths.get('chain_blast')
            if chain_blast and os.path.exists(chain_blast):
                self.logger.info(f"Chain BLAST exists: {chain_blast}")
            else:
                self.logger.warning(f"Chain BLAST missing: {chain_blast}")

            domain_blast = standard_paths.get('domain_blast')
            if domain_blast and os.path.exists(domain_blast):
                self.logger.info(f"Domain BLAST exists: {domain_blast}")
            else:
                self.logger.warning(f"Domain BLAST missing: {domain_blast}")

        return summary_path is not None

    def run_domain_partition_for_process(self, process_id, force=False):
        """Run domain partition for a specific process with enhanced path handling"""
        process = self.get_process_details(process_id)
        if not process:
            self.logger.error(f"Process {process_id} not found")
            return False

        pdb_id = process['pdb_id']
        chain_id = process['chain_id']
        batch_path = process['base_path']
        ref_version = process['ref_version']

        # Verify domain summary exists before attempting partition
        summary_path = self.verify_domain_summary_exists(pdb_id, chain_id, batch_path, ref_version)
        if not summary_path:
            self.logger.error(f"Cannot run domain partition - domain summary missing for {pdb_id}_{chain_id}")

            # Update process status to indicate the issue
            self.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": "domain_partition_failed",
                    "status": "error",
                    "error_message": "Domain summary file not found"
                },
                "id = %s",
                (process_id,)
            )
            return False

        self.logger.info(f"Running domain partition for {pdb_id}_{chain_id}")

        # Use the domain pipeline to process this protein
        try:
            success = self.domain_pipeline.process_proteins(
                batch_id=process['batch_id'],
                protein_ids=[process_id],
                blast_only=False,  # Use full pipeline results
                partition_only=True  # Only run the partition step
            )

            # Additional verification to make sure the domain partition file was created
            paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version)
            domain_fn = paths['domain_partition']

            if os.path.exists(domain_fn):
                self.logger.info(f"Successfully created domain partition file: {domain_fn}")

                # Register the file in the database for tracking
                self._register_domain_file(process_id, os.path.relpath(domain_fn, batch_path))
                return True
            else:
                self.logger.error(f"Domain partition process completed but file not found: {domain_fn}")
                return False

        except Exception as e:
            self.logger.error(f"Error running domain partition for {pdb_id}_{chain_id}: {e}", exc_info=True)

            # Update process status
            self.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": "domain_partition_failed",
                    "status": "error",
                    "error_message": str(e)
                },
                "id = %s",
                (process_id,)
            )
            return False

    def _register_domain_file(self, process_id, file_path):
        """Register domain partition file in database with proper duplicate handling"""
        try:
            # Check if record already exists
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = 'domain_partition'
            """

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
            else:
                # Insert new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "domain_partition",
                        "file_path": file_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    }
                )
        except Exception as e:
            self.logger.warning(f"Error registering domain file: {e}")

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Run ECOD Domain Partition with Enhanced Path Handling')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int,
                      help='Batch ID to process')
    parser.add_argument('--process-id', type=int,
                      help='Specific process ID to check or process')
    parser.add_argument('--check-paths', action='store_true',
                      help='Check paths for process without running partition')
    parser.add_argument('--limit', type=int,
                      help='Limit the number of proteins to process')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed proteins')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)

    logger = logging.getLogger("main")

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    runner = EnhancedPartitionRunner(args.config)

    # Process a specific process ID
    if args.process_id:
        logger.info(f"Working with specific process ID: {args.process_id}")

        if args.check_paths:
            # Just check paths without running partition
            success = runner.check_paths_for_process(args.process_id)
            return 0 if success else 1
        else:
            # Run domain partition for specific process
            success = runner.run_domain_partition_for_process(args.process_id, args.force)
            return 0 if success else 1

    # Batch processing is not implemented in this example but would go here
    if args.batch_id:
        logger.error("Batch processing not implemented in this example")
        return 1

    logger.error("Please specify --process-id or --batch-id")
    return 1

if __name__ == "__main__":
    sys.exit(main())
