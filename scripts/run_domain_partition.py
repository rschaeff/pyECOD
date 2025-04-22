#!/usr/bin/env python3
"""
run_domain_partition.py - Run domain partition for representative processes in batch 31
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

class PartitionRunner:
    """Runner for domain partition"""
    
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
    
    def get_representative_process_ids(self, batch_id):
        """Get process IDs for representative proteins with full domain summaries"""
        query = """
        SELECT 
            ps.id as process_id, p.pdb_id, p.chain_id
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 
            ps.batch_id = %s
            AND ps.is_representative = TRUE
            AND pf.file_type = 'domain_summary'
            AND pf.file_exists = TRUE
            AND ps.current_stage = 'domain_summary_complete'
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            process_ids = [row['process_id'] for row in rows]
            
            self.logger.info(f"Found {len(process_ids)} representative processes with domain summaries")
            return process_ids
        except Exception as e:
            self.logger.error(f"Error finding representative processes: {str(e)}")
            return []

    def run_domain_partition(self, batch_id, limit=None, force=False):
        """Run domain partition for representative processes with full pipeline results"""
        self.logger.info(f"Running domain partition for batch {batch_id}")
        
        # Get batch info
        batch_query = """
        SELECT base_path, ref_version 
        FROM ecod_schema.batch 
        WHERE id = %s
        """
        
        batch_info = self.db.execute_dict_query(batch_query, (batch_id,))
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return False
        
        base_path = batch_info[0]['base_path']
        ref_version = batch_info[0]['ref_version']
        
        # Get representative process IDs
        process_ids = self.get_representative_process_ids(batch_id)
        
        if not process_ids:
            self.logger.error("No representative processes found with domain summaries")
            return False
        
        if limit and limit < len(process_ids):
            process_ids = process_ids[:limit]
            self.logger.info(f"Limited to {limit} processes")
        
        # Run partition specifically for these processes
        self.logger.info(f"Running domain partition for {len(process_ids)} representative processes")
        success = self.domain_pipeline.process_proteins(
            batch_id=batch_id,
            protein_ids=process_ids,
            blast_only=False,  # Use full pipeline results
            partition_only=True  # Only run the partition step
        )
        
        return success

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Run ECOD Domain Partition for Representative Processes')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, default=31,
                      help='Batch ID to process (default: 31)')
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
    logger.info(f"Starting domain partition for representative processes in batch {args.batch_id}")
    
    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)
    
    runner = PartitionRunner(args.config)
    success = runner.run_domain_partition(args.batch_id, args.limit, args.force)
    
    if success:
        logger.info(f"Successfully ran domain partition for batch {args.batch_id}")
        return 0
    else:
        logger.error(f"Failed to run domain partition for batch {args.batch_id}")
        return 1

if __name__ == "__main__":
    sys.exit(main())q