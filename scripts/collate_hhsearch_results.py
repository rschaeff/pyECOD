#!/usr/bin/env python3
"""
collate_hhsearch_results.py - Collate HHSearch and BLAST results for batch 31
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
from ecod.pipelines.domain_analysis.summary import DomainSummary

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

class CollationRunner:
    """Runner for collating HHSearch results"""
    
    def __init__(self, config_path=None):
        """Initialize with configuration"""
        self.logger = logging.getLogger("ecod.collation_runner")
        
        # Load configuration
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        
        # Initialize database connection
        db_config = self.config_manager.get_db_config()
        self.db = DBManager(db_config)
        
        # Initialize domain summary component with self as context
        self.domain_summary = DomainSummary(self)
    
    def get_db(self):
        """Access method for database - mimics ApplicationContext interface"""
        return self.db
        
    def is_force_overwrite(self):
        """Check if force overwrite is enabled - mimics ApplicationContext interface"""
        return self.config.get('pipeline', {}).get('force_overwrite', False)
    
    def collate_batch_results(self, batch_id, force=False, limit=None):
        """Collate HHSearch results with BLAST results for representative processes in a batch"""
        self.logger.info(f"Collating HHSearch and BLAST results for batch {batch_id}")
        
        # Get batch info
        batch_query = """
        SELECT id, batch_name, base_path, ref_version 
        FROM ecod_schema.batch 
        WHERE id = %s
        """
        
        batch_info = self.db.execute_dict_query(batch_query, (batch_id,))
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return False
        
        base_path = batch_info[0]['base_path']
        ref_version = batch_info[0]['ref_version']
        batch_name = batch_info[0]['batch_name']
        
        self.logger.info(f"Processing batch: {batch_name} with reference {ref_version}")
        
        # Get representative proteins with HHSearch results
        protein_query = """
        SELECT 
            p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'hhr'
        WHERE 
            ps.batch_id = %s
            AND ps.is_representative = TRUE
            AND pf.file_exists = TRUE
        ORDER BY
            p.pdb_id, p.chain_id
        """
        
        proteins = self.db.execute_dict_query(protein_query, (batch_id,))
        self.logger.info(f"Found {len(proteins)} representative proteins with HHSearch results")
        
        if limit and limit < len(proteins):
            proteins = proteins[:limit]
            self.logger.info(f"Limited to {limit} proteins")
        
        success_count = 0
        for protein in proteins:
            pdb_id = protein['pdb_id']
            chain_id = protein['chain_id']
            process_id = protein['process_id']
            
            self.logger.info(f"Processing {pdb_id}_{chain_id}")
            
            # Check if a full pipeline summary already exists
            if os.path.exists(regular_summary_path) and not force:
                # Verify this is a full pipeline summary (contains HHSearch evidence)
                is_full_summary = self._check_for_hhsearch_evidence(regular_summary_path)
                
                if is_full_summary:
                    self.logger.info(f"Full pipeline domain summary already exists for {pdb_id}_{chain_id}, skipping")
                    success_count += 1
                    continue
                else:
                    self.logger.info(f"Found blast-only summary for {pdb_id}_{chain_id}, replacing with full pipeline version")
                
                # Create domain summary with HHSearch evidence
                try:
                    # Call DomainSummary's create_summary with blast_only=False
                    summary_file = self.domain_summary.create_summary(
                        pdb_id, 
                        chain_id, 
                        ref_version, 
                        base_path, 
                        blast_only=False  # Use HHSearch results
                    )
            
            if summary_file:
                # Register summary file in database
                self._register_summary(process_id, summary_file, base_path)
                success_count += 1
                self.logger.info(f"Successfully created full domain summary for {pdb_id}_{chain_id}")
            else:
                self.logger.warning(f"Failed to create domain summary for {pdb_id}_{chain_id}")
        
            except Exception as e:
                self.logger.error(f"Error processing {pdb_id}_{chain_id}: {str(e)}")
        
        self.logger.info(f"Successfully collated results for {success_count}/{len(proteins)} proteins")
        return success_count > 0
        
    def _register_summary(self, process_id, summary_file, base_path):
        """Register domain summary in database"""
        try:
            relative_path = os.path.relpath(summary_file, base_path)
            file_size = os.path.getsize(summary_file)
            
            # Check if summary already registered
            check_query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = 'domain_summary'
            """
            
            existing = self.db.execute_query(check_query, (process_id,))
            
            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": relative_path,
                        "file_exists": True,
                        "file_size": file_size
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
                        "file_type": "domain_summary",
                        "file_path": relative_path,
                        "file_exists": True,
                        "file_size": file_size
                    }
                )
            
            # Update process status
            self.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": "domain_summary_complete",
                    "status": "success"
                },
                "id = %s",
                (process_id,)
            )
            return True
        except Exception as e:
            self.logger.error(f"Error registering summary: {str(e)}")
            return False

    def _check_for_hhsearch_evidence(self, summary_file):
        """Check if a domain summary contains HHSearch evidence"""
        try:
            import xml.etree.ElementTree as ET
            tree = ET.parse(summary_file)
            root = tree.getroot()
            
            # Look for HHSearch evidence section
            hhsearch_elem = root.find(".//hhsearch_evidence")
            if hhsearch_elem is None:
                return False
                
            # Check if it has any hit elements
            hits = hhsearch_elem.findall(".//hh_hit")
            return len(hits) > 0
        except Exception as e:
            self.logger.error(f"Error checking summary for HHSearch evidence: {str(e)}")
            return False

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Collate ECOD HHSearch Results with BLAST')
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
                      help='Force reprocessing of already processed results')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    logger = logging.getLogger("main")
    logger.info(f"Starting collation of BLAST and HHSearch results for batch {args.batch_id}")
    
    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)
    
    runner = CollationRunner(args.config)
    
    success = runner.collate_batch_results(args.batch_id, args.force, args.limit)
    
    if success:
        logger.info(f"Successfully collated results for batch {args.batch_id}")
        return 0
    else:
        logger.error(f"Failed to collate results for batch {args.batch_id}")
        return 1

if __name__ == "__main__":
    sys.exit(main())