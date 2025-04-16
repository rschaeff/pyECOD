#!/usr/bin/env python3
"""
remove_invalid_summaries.py - Remove domain summaries that fail validation
"""

import os
import sys
import logging
import argparse
from typing import List, Optional
import xml.etree.ElementTree as ET

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

class InvalidSummaryRemover:
    """Finds and removes invalid domain summaries"""
    
    def __init__(self, context, batch_id, batch_path):
        """Initialize with context and batch info"""
        self.context = context
        self.batch_id = batch_id
        self.batch_path = batch_path
        self.logger = logging.getLogger("ecod.remove_invalid_summaries")
        self.removed_count = 0
        
    def check_and_remove_file(self, file_path: str, pdb_id: str, chain_id: str, dry_run: bool = True) -> bool:
        """Check if a domain summary is invalid and remove it if so"""
        try:
            # Check if file exists
            if not os.path.exists(file_path):
                return False
                
            # Try to parse the XML
            tree = ET.parse(file_path)
            root = tree.getroot()
            
            # Look for common validation issues
            blast_summ = root.find(".//blast_summ")
            if blast_summ is None:
                self.logger.info(f"Removing invalid summary - missing blast_summ element: {file_path}")
                return self._remove_summary(file_path, pdb_id, chain_id, dry_run)
                
            # Check for missing required run elements
            chain_blast_run = root.find(".//chain_blast_run")
            domain_blast_run = root.find(".//blast_run")
            
            # Only remove if both elements are missing and we don't have the _no_hits flags
            no_chain_hits = blast_summ.get("chain_blast_no_hits") == "true"
            no_domain_hits = blast_summ.get("domain_blast_no_hits") == "true"
            
            if chain_blast_run is None and not no_chain_hits:
                self.logger.info(f"Removing invalid summary - missing chain_blast_run element: {file_path}")
                return self._remove_summary(file_path, pdb_id, chain_id, dry_run)
                
            if domain_blast_run is None and not no_domain_hits:
                self.logger.info(f"Removing invalid summary - missing blast_run element: {file_path}")
                return self._remove_summary(file_path, pdb_id, chain_id, dry_run)
                
            # File appears valid
            return False
            
        except Exception as e:
            self.logger.info(f"Removing invalid summary - parsing error: {file_path}, {str(e)}")
            return self._remove_summary(file_path, pdb_id, chain_id, dry_run)
            
    def _remove_summary(self, file_path: str, pdb_id: str, chain_id: str, dry_run: bool) -> bool:
        """Remove a summary file and update the database"""
        if dry_run:
            self.logger.info(f"[DRY RUN] Would remove: {file_path}")
            self.removed_count += 1
            return True
            
        try:
            # Get relative path
            rel_path = os.path.relpath(file_path, self.batch_path) if self.batch_path in file_path else file_path
            
            # Update database record first
            process_query = """
            SELECT pf.id
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE ps.batch_id = %s AND p.pdb_id = %s AND p.chain_id = %s
            AND pf.file_type = 'domain_summary'
            """
            
            result = self.context.db.execute_query(process_query, (self.batch_id, pdb_id, chain_id))
            
            if result and len(result) > 0:
                file_id = result[0][0]
                
                # Use the update method to set file_exists = False
                self.context.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_exists": False
                    },
                    "id = %s",
                    (file_id,)
                )
                
                self.logger.debug(f"Updated database record for {pdb_id}_{chain_id}")
            else:
                self.logger.warning(f"No database record found for {pdb_id}_{chain_id}")
            
            # Now remove the file
            os.remove(file_path)
            self.removed_count += 1
            self.logger.info(f"Removed file: {file_path}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error removing {file_path}: {str(e)}")
            return False

def remove_invalid_summaries(batch_id: int, config_path: str, 
                           dry_run: bool = True, filter_chain: Optional[str] = None):
    """Remove invalid domain summaries for a batch"""
    context = ApplicationContext(config_path)
    logger = logging.getLogger("ecod.remove_invalid_summaries")
    
    # Get batch path
    batch_query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return 1
        
    batch_path = batch_result[0][0]
    
    # Get all domain summaries for the batch
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, pf.file_path
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id AND pf.file_type = 'domain_summary'
    WHERE 
        ps.batch_id = %s AND pf.file_exists = TRUE
    ORDER BY p.pdb_id, p.chain_id
    """
    
    proteins = context.db.execute_query(query, (batch_id,))
    
    if not proteins:
        logger.error(f"No domain summaries found for batch {batch_id}")
        return 1
    
    logger.info(f"Found {len(proteins)} proteins with domain summaries")
    
    # Initialize validator
    remover = InvalidSummaryRemover(context, batch_id, batch_path)
    
    # Process each summary
    for i, protein in enumerate(proteins):
        pdb_id = protein[1]
        chain_id = protein[2]
        summary_path = protein[3]
        
        # Apply filter if specified
        if filter_chain and f"{pdb_id}_{chain_id}" != filter_chain:
            continue
        
        # Prepare absolute path if needed
        summary_path = os.path.join(batch_path, summary_path) if not os.path.isabs(summary_path) else summary_path
        
        # Check and remove if invalid
        remover.check_and_remove_file(summary_path, pdb_id, chain_id, dry_run)
        
        # Progress update
        if (i + 1) % 100 == 0 or (i + 1) == len(proteins):
            logger.info(f"Progress: {i+1}/{len(proteins)} processed")
    
    logger.info(f"Summary: {remover.removed_count} invalid files would be removed")
    
    if dry_run:
        logger.info("DRY RUN - No files were actually removed")
        logger.info("Run with --apply to actually remove files")
    
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Remove invalid domain summaries')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to validate')
    parser.add_argument('--filter', type=str,
                      help='Filter to specific protein (format: pdb_chain, e.g., 8gh6_P)')
    parser.add_argument('--apply', action='store_true',
                      help='Actually remove files (default is dry run)')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    sys.exit(remove_invalid_summaries(
        args.batch_id, 
        args.config, 
        dry_run=not args.apply, 
        filter_chain=args.filter
    ))