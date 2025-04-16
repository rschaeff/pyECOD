#!/usr/bin/env python3
"""
Script to fix the 63 proteins in batch 31 that are stuck at the domain_summary stage.
This script addresses various issues identified by the diagnostic tool and updates
the database to move these proteins to the domain_partition_complete stage.
"""

import os
import sys
import argparse
import logging
import xml.etree.ElementTree as ET
from datetime import datetime
import psycopg2
from psycopg2.extras import DictCursor
import yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("fix_batch31.log")
    ]
)
logger = logging.getLogger('fix_batch31')

class BatchFixer:
    """Class to handle fixing proteins stuck at domain_summary stage"""
    
    def __init__(self, config):
        """Initialize the batch fixer"""
        self.config = config
        self.conn = None
        self.cursor = None
        
        # Get database configuration
        self.db_config = self.config.get('database', {})
        
        # Get data paths from configuration
        data_paths = self.config.get('paths', {})
        self.data_dir = data_paths.get('data_dir', '/data/ecod')
        
        self.batch_dir = os.path.join(
            self.data_dir,
            data_paths.get('updates_dir', 'pdb_updates'),
            'batches',
            str(self.config['batch_id'])
        )
        
        self.reference_version = self.config.get('reference_version', 'develop291')
        logger.info(f"Using batch directory: {self.batch_dir}")
        logger.info(f"Using reference version: {self.reference_version}")
    
    def connect_db(self):
        """Connect to the database"""
        try:
            self.conn = psycopg2.connect(
                host=self.db_config.get('host', 'localhost'),
                port=self.db_config.get('port', 5432),
                dbname=self.db_config.get('dbname', 'ecod'),
                user=self.db_config.get('user', 'ecod'),
                password=self.db_config.get('password', '')
            )
            self.cursor = self.conn.cursor(cursor_factory=DictCursor)
            logger.info("Connected to database")
            return True
        except Exception as e:
            logger.error(f"Database connection error: {e}")
            return False
    
    def close_db(self):
        """Close the database connection"""
        if self.cursor:
            self.cursor.close()
        if self.conn:
            self.conn.close()
        logger.info("Database connection closed")
    
    def commit(self):
        """Commit database changes"""
        if self.conn:
            self.conn.commit()
    
    def rollback(self):
        """Rollback database changes"""
        if self.conn:
            self.conn.rollback()
    
    def get_stuck_proteins(self):
        """Get proteins stuck at domain_summary stage"""
        query = """
        SELECT p.id, p.pdb_id, p.chain_id, p.length, ps.status, ps.current_stage, ps.id as process_id
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.batch_id = %s
        AND ps.current_stage = 'domain_summary'
        AND ps.status = 'success'
        ORDER BY p.pdb_id, p.chain_id
        """
        self.cursor.execute(query, (self.config['batch_id'],))
        return self.cursor.fetchall()
    
    def get_file_path(self, process_id, file_type):
        """Get file path from database"""
        query = """
        SELECT file_path
        FROM ecod_schema.process_file
        WHERE process_id = %s
        AND file_type = %s
        """
        self.cursor.execute(query, (process_id, file_type))
        result = self.cursor.fetchone()
        return result['file_path'] if result else None
    
    def _parse_range(self, range_str):
        """Parse range string to extract start and end positions"""
        if not range_str:
            return None
        
        # Handle formats like "2-923" or "1:100"
        if '-' in range_str:
            parts = range_str.split('-')
            start = int(parts[0])
            end = int(parts[1])
            return (start, end)
        elif ':' in range_str:
            parts = range_str.split(':')
            start = int(parts[0])
            end = int(parts[1])
            return (start, end)
        else:
            # Handle single position
            try:
                pos = int(range_str)
                return (pos, pos)
            except ValueError:
                logger.warning(f"Could not parse range: {range_str}")
                return None
    
    def extract_domains_from_summary(self, file_path, min_domain_length=14):
        """
        Extract domain information from domain summary file
        
        Args:
            file_path: Path to domain summary XML file
            min_domain_length: Minimum domain length to consider valid
            
        Returns:
            List of domain ranges [(start1, end1), (start2, end2), ...]
        """
        if not os.path.exists(file_path):
            logger.error(f"Domain summary file not found: {file_path}")
            return []
        
        try:
            tree = ET.parse(file_path)
            root = tree.getroot()
            
            domains = []
            
            # Check for chain-level BLAST hits
            chain_blast_runs = root.findall('.//chain_blast_run')
            for chain_run in chain_blast_runs:
                hits = chain_run.findall('.//hit')
                for hit in hits:
                    query_regions = hit.findall('.//query_reg')
                    for qr in query_regions:
                        range_str = qr.get('range') or qr.text
                        if range_str:
                            domain_range = self._parse_range(range_str)
                            if domain_range:
                                start, end = domain_range
                                if end - start + 1 >= min_domain_length:
                                    domains.append(domain_range)
            
            # Check for domain-level BLAST hits
            domain_blast_runs = root.findall('.//blast_run')
            for domain_run in domain_blast_runs:
                hits = domain_run.findall('.//hit')
                for hit in hits:
                    query_regions = hit.findall('.//query_reg')
                    for qr in query_regions:
                        range_str = qr.get('range') or qr.text
                        if range_str:
                            domain_range = self._parse_range(range_str)
                            if domain_range:
                                start, end = domain_range
                                if end - start + 1 >= min_domain_length:
                                    domains.append(domain_range)
            
            # If no domains found, check for protein length
            if not domains:
                # Look for sequence length info
                seq_length_elem = root.find('.//sequence_length')
                if seq_length_elem is not None and seq_length_elem.text:
                    seq_length = int(seq_length_elem.text)
                    if seq_length >= min_domain_length:
                        domains.append((1, seq_length))
                else:
                    # Look for protein info
                    protein_elem = root.find('.//protein')
                    if protein_elem is not None:
                        length_attr = protein_elem.get('length')
                        if length_attr and int(length_attr) >= min_domain_length: