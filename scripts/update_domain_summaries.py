#!/usr/bin/env python3
"""
update_domain_summaries.py - Update domain summaries with HHSearch evidence

This script:
1. Finds chains with HHSearch results but outdated domain summaries
2. Updates domain summaries to include HHSearch evidence
3. Enables better domain boundary determination with the enhanced evidence

Usage:
    python update_domain_summaries.py --batch-id <batch_id> [--force] [--chains <pdb_id>_<chain_id> ...]
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.summary import DomainSummary
from ecod.pipelines.domain_analysis.partition import DomainPartition


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


class DomainSummaryUpdater:
    """Update domain summaries with HHSearch evidence"""
    
    def __init__(self, context=None):
        """Initialize with application context"""
        self.context = context or ApplicationContext()
        self.logger = logging.getLogger("ecod.domain_summary_updater")
        
        # Initialize domain analysis components
        self.summary_processor = DomainSummary(self.context)
        self.domain_partition = DomainPartition(self.context)
    
    def update_batch_summaries(self, batch_id: int, force: bool = False, 
                             specific_chains: Optional[List[str]] = None) -> Tuple[int, int]:
        """Update domain summaries for a batch
        
        Args:
            batch_id: Batch ID to process
            force: Force regeneration of files even if they exist
            specific_chains: List of specific chain IDs to process (pdbid_chainid)
            
        Returns:
            Tuple of (updated_summaries, updated_partitions)
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return (0, 0)
        
        # Find chains that need updating
        chains_to_update = self._find_chains_for_update(batch_id, batch_info, specific_chains)
        if not chains_to_update:
            self.logger.warning(f"No chains need domain summary updates in batch {batch_id}")
            return (0, 0)
        
        # Update summaries and partitions
        updated_summaries = 0
        updated_partitions = 0
        
        for chain in chains_to_update:
            # Update domain summary
            summary_updated = self._update_domain_summary(
                chain['pdb_id'],
                chain['chain_id'],
                batch_info['ref_version'],
                batch_info['base_path'],
                force
            )
            
            if summary_updated:
                updated_summaries += 1
                
                # Update domain partition
                partition_updated = self._update_domain_partition(
                    chain['pdb_id'],
                    chain['chain_id'],
                    batch_info['ref_version'],
                    batch_info['base_path'],
                    force
                )
                
                if partition_updated:
                    updated_partitions += 1
        
        self.logger.info(f"Updated {updated_summaries} domain summaries and {updated_partitions} domain partitions")
        return (updated_summaries, updated_partitions)
    
    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information
        
        Args:
            batch_id: Batch ID
            
        Returns:
            Dictionary with batch information
        """
        query = """
        SELECT 
            id, 
            batch_name, 
            base_path, 
            ref_version
        FROM 
            ecod_schema.batch
        WHERE 
            id = %s
        """
        
        try:
            rows = self.context.db.execute_dict_query(query, (batch_id,))
            if rows:
                return rows[0]
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {e}")
            
        return None
    
    def _find_chains_for_update(self, batch_id: int, batch_info: Dict[str, Any], 
                             specific_chains: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """Find chains that need domain summary updates
        
        Args:
            batch_id: Batch ID
            batch_info: Batch information dictionary
            specific_chains: List of specific chain IDs to process (pdbid_chainid)
            
        Returns:
            List of dictionaries with chain information
        """
        chains = []
        
        if specific_chains:
            self.logger.info(f"Looking for domain summaries to update for {len(specific_chains)} specific chains")
            
            for chain_id in specific_chains:
                try:
                    pdb_id, chain_letter = chain_id.split('_')
                    
                    # Get process ID and check if it has HHSearch results
                    query = """
                    SELECT 
                        ps.id as process_id,
                        p.id as protein_id,
                        p.pdb_id,
                        p.chain_id,
                        ps.relative_path,
                        ps.current_stage
                    FROM 
                        ecod_schema.process_status ps
                    JOIN
                        ecod_schema.protein p ON ps.protein_id = p.id
                    LEFT JOIN
                        ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'hhsearch_xml'
                    WHERE 
                        ps.batch_id = %s
                        AND p.pdb_id = %s
                        AND p.chain_id = %s
                        AND pf.file_exists = TRUE
                    LIMIT 1
                    """
                    
                    rows = self.context.db.execute_dict_query(query, (batch_id, pdb_id, chain_letter))
                    if rows:
                        chains.append(rows[0])
                except ValueError:
                    self.logger.warning(f"Invalid chain ID format: {chain_id}, expected pdbid_chainid")
            
            self.logger.info(f"Found {len(chains)} chains that need domain summary updates")
        else:
            # Find chains with HHSearch results but without updated domain summaries
            # Corrected query to properly use the schema
            query = """
            SELECT 
                ps.id as process_id,
                p.id as protein_id,
                p.pdb_id,
                p.chain_id,
                ps.relative_path,
                ps.current_stage
            FROM 
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            JOIN
                ecod_schema.process_file pf1 ON ps.id = pf1.process_id AND pf1.file_type = 'hhsearch_xml'
            LEFT JOIN
                ecod_schema.process_file pf2 ON ps.id = pf2.process_id AND pf2.file_type = 'domain_summary'
            WHERE 
                ps.batch_id = %s
                AND pf1.file_exists = TRUE
                AND (
                    ps.current_stage = 'hhsearch_complete'
                    OR (ps.current_stage = 'domain_summary_complete' AND ps.updated_at < pf1.last_checked)
                    OR pf2.id IS NULL
                    OR pf2.file_exists = FALSE
                )
            ORDER BY ps.id
            """
            
            rows = self.context.db.execute_dict_query(query, (batch_id,))
            self.logger.info(f"Found {len(rows)} chains with HHSearch results that need domain summary updates")
            
            chains = rows
        
        return chains
    
    def _update_domain_summary(self, pdb_id: str, chain_id: str, ref_version: str, 
                            base_path: str, force: bool = False) -> bool:
        """Update domain summary with HHSearch evidence
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            ref_version: Reference version
            base_path: Base path
            force: Force regeneration even if file exists
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Set force regeneration in context if needed
            if force:
                self.context.config_manager.set_force_overwrite(True)
            
            # Use the DomainSummary class to regenerate the summary
            # This will automatically incorporate HHSearch evidence if available
            summary_path = self.summary_processor.create_summary(
                pdb_id,
                chain_id,
                ref_version,
                base_path,
                False  # Not blast_only
            )
            
            # Reset force setting
            if force:
                self.context.config_manager.set_force_overwrite(False)
            
            if summary_path and os.path.exists(summary_path):
                self.logger.info(f"Updated domain summary for {pdb_id}_{chain_id}")
                return True
            else:
                self.logger.error(f"Failed to update domain summary for {pdb_id}_{chain_id}")
                return False
        except Exception as e:
            self.logger.error(f"Error updating domain summary for {pdb_id}_{chain_id}: {str(e)}")
            return False
    
     def _update_domain_partition(self, pdb_id: str, chain_id: str, ref_version: str, 
                              base_path: str, force: bool = False) -> bool:
        """Update domain partition with improved boundaries
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            ref_version: Reference version
            base_path: Base path
            force: Force regeneration even if file exists
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Set force regeneration in context if needed
            if force:
                self.context.config_manager.set_force_overwrite(True)
            
            # Use the DomainPartition class to regenerate the partition
            # This will use the updated domain summary with HHSearch evidence
            partition_path = self.domain_partition.partition_domains(
                pdb_id,
                chain_id,
                base_path,
                'struct_seqid',  # Default input mode
                ref_version,
                False  # Not blast_only
            )
            
            # Reset force setting
            if force:
                self.context.config_manager.set_force_overwrite(False)
            
            if partition_path and os.path.exists(partition_path):
                self.logger.info(f"Updated domain partition for {pdb_id}_{chain_id}")
                return True
            else:
                self.logger.error(f"Failed to update domain partition for {pdb_id}_{chain_id}")
                return False
        except Exception as e:
            self.logger.error(f"Error updating domain partition for {pdb_id}_{chain_id}: {str(e)}")
            return False


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Update domain summaries with HHSearch evidence')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--chains', nargs='+', default=None,
                      help='Specific chains to process (format: pdbid_chainid)')
    parser.add_argument('--force', action='store_true',
                      help='Force regeneration of domain summaries and partitions')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    logger = logging.getLogger("ecod.update_domain_summaries")
    logger.info(f"Starting domain summary updates for batch {args.batch_id}")
    
    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Create updater
    updater = DomainSummaryUpdater(context)
    
    try:
        # Update domain summaries
        updated_summaries, updated_partitions = updater.update_batch_summaries(
            args.batch_id, 
            args.force, 
            args.chains
        )
        
        logger.info(f"Successfully updated {updated_summaries} domain summaries and {updated_partitions} domain partitions")
        
        if updated_summaries == 0:
            return 1
        
        return 0
    except Exception as e:
        logger.error(f"Error updating domain summaries: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())