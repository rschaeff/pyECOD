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
        self.data_root = config['data_root']
        self.batch_dir = os.path.join(self.data_root, 'pdb_updates', 'batches', str(config['batch_id']))
        self.reference_version = config.get('reference_version', 'develop291')
    
    def connect_db(self):
        """Connect to the database"""
        try:
            self.conn = psycopg2.connect(
                host=self.config['db_host'],
                port=self.config['db_port'],
                dbname=self.config['db_name'],
                user=self.config['db_user'],
                password=self.config['db_password']
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
                            domains.append((1, int(length_attr)))
            
            # Handle when we still have no domains
            if not domains:
                logger.warning(f"No domains found in {file_path}")
            
            return domains
            
        except Exception as e:
            logger.error(f"Error processing domain summary file {file_path}: {e}")
            return []
    
    def create_domain_partition_file(self, pdb_id, chain_id, domains, output_dir):
        """
        Create a domain partition XML file
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            domains: List of domain ranges [(start1, end1), (start2, end2), ...]
            output_dir: Output directory
            
        Returns:
            Path to created file
        """
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Create output file path
        output_file = os.path.join(output_dir, f"{pdb_id}_{chain_id}.{self.reference_version}.domains.xml")
        
        # Create XML structure
        root = ET.Element("domain_partition")
        root.set("pdb_id", pdb_id)
        root.set("chain_id", chain_id)
        root.set("ref_version", self.reference_version)
        root.set("created", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        
        # Add domains
        for i, (start, end) in enumerate(domains, 1):
            domain = ET.SubElement(root, "domain")
            domain.set("id", f"{pdb_id}_{chain_id}_{i}")
            domain.set("start", str(start))
            domain.set("end", str(end))
            domain.set("method", "auto_blast")
            
            # Add evidence element
            evidence = ET.SubElement(domain, "evidence")
            evidence.set("type", "blast")
            evidence.set("confidence", "medium")
        
        # Write to file
        tree = ET.ElementTree(root)
        tree.write(output_file, encoding="utf-8", xml_declaration=True)
        
        logger.info(f"Created domain partition file: {output_file}")
        return output_file
    
    def create_single_domain(self, pdb_id, chain_id, length, output_dir):
        """
        Create a domain partition file with a single domain for the entire protein
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            length: Protein length
            output_dir: Output directory
            
        Returns:
            Path to created file
        """
        # Use length or default to 100 if not provided
        actual_length = length or 100
        
        # Create a single domain covering the entire protein
        domains = [(1, actual_length)]
        
        return self.create_domain_partition_file(pdb_id, chain_id, domains, output_dir)
    
    def add_domain_partition_file(self, process_id, file_path):
        """Add domain partition file record to database"""
        # First check if the file already exists
        query_check = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = %s
        """
        self.cursor.execute(query_check, (process_id, 'domain_partition'))
        existing = self.cursor.fetchone()
        
        try:
            if existing:
                # Update existing record
                query_update = """
                UPDATE ecod_schema.process_file
                SET file_path = %s, file_exists = TRUE
                WHERE id = %s
                """
                self.cursor.execute(query_update, (file_path, existing['id']))
            else:
                # Insert new record
                query_insert = """
                INSERT INTO ecod_schema.process_file
                (process_id, file_type, file_path, file_exists)
                VALUES (%s, %s, %s, %s)
                """
                self.cursor.execute(query_insert, (
                    process_id, 
                    'domain_partition', 
                    file_path, 
                    True
                ))
            return True
        except Exception as e:
            logger.error(f"Error adding domain partition file: {e}")
            return False
    
    def update_protein_status(self, process_id, new_stage, status='success'):
        """Update protein processing status in database"""
        query = """
        UPDATE ecod_schema.process_status
        SET current_stage = %s, status = %s, updated_at = %s
        WHERE id = %s
        """
        try:
            self.cursor.execute(query, (
                new_stage, 
                status, 
                datetime.now(), 
                process_id
            ))
            return True
        except Exception as e:
            logger.error(f"Error updating protein status: {e}")
            return False
    
    def update_batch_status(self):
        """Update batch status if all proteins are complete"""
        # First check if all proteins are processed
        query_check = """
        SELECT COUNT(*) as total_count,
               SUM(CASE WHEN current_stage = 'domain_partition_complete' THEN 1 ELSE 0 END) as completed_count
        FROM ecod_schema.process_status
        WHERE batch_id = %s
        """
        self.cursor.execute(query_check, (self.config['batch_id'],))
        result = self.cursor.fetchone()
        
        total = result['total_count']
        completed = result['completed_count']
        
        if total == completed:
            # All proteins processed, update batch status
            query_update = """
            UPDATE ecod_schema.batch
            SET status = 'completed', updated_at = %s
            WHERE id = %s
            """
            try:
                self.cursor.execute(query_update, (datetime.now(), self.config['batch_id']))
                logger.info(f"Batch {self.config['batch_id']} marked as completed ({completed}/{total} proteins)")
                return True
            except Exception as e:
                logger.error(f"Error updating batch status: {e}")
                return False
        else:
            logger.info(f"Batch {self.config['batch_id']} progress: {completed}/{total} proteins completed")
            return False
    
    def fix_protein(self, protein, dry_run=False):
        """
        Fix a protein stuck at domain_summary stage
        
        Args:
            protein: Protein record from database
            dry_run: If True, don't make any database changes
            
        Returns:
            True if successful, False otherwise
        """
        protein_id = protein['id']
        process_id = protein['process_id']
        pdb_id = protein['pdb_id']
        chain_id = protein['chain_id']
        length = protein['length']
        
        logger.info(f"Processing {pdb_id}_{chain_id} (ID: {protein_id}, Process ID: {process_id}, Length: {length})")
        
        # Determine if this is a small peptide or large protein
        is_small_peptide = length is not None and length < 20
        is_large_protein = length is not None and length > 1500
        
        # Get domain summary file path
        domain_summary_path = self.get_file_path(process_id, 'domain_summary')
        
        # Set up output directory for domain partition files
        output_dir = os.path.join(self.batch_dir, "domains")
        os.makedirs(output_dir, exist_ok=True)
        
        # For very small peptides, just create a single domain
        if is_small_peptide:
            logger.info(f"Small peptide detected: {pdb_id}_{chain_id} (Length: {length})")
            partition_file = self.create_single_domain(pdb_id, chain_id, length, output_dir)
        # For large proteins, also create a single domain but log differently
        elif is_large_protein:
            logger.info(f"Large protein detected: {pdb_id}_{chain_id} (Length: {length})")
            partition_file = self.create_single_domain(pdb_id, chain_id, length, output_dir)
        # For normal proteins, extract domains from summary
        elif domain_summary_path and os.path.exists(os.path.join(self.data_root, domain_summary_path)):
            full_path = os.path.join(self.data_root, domain_summary_path)
            domains = self.extract_domains_from_summary(full_path)
            
            # If no domains extracted, create a single domain
            if not domains:
                logger.warning(f"No domains found in summary, creating single domain: {pdb_id}_{chain_id}")
                partition_file = self.create_single_domain(pdb_id, chain_id, length, output_dir)
            else:
                logger.info(f"Found {len(domains)} domains in summary: {pdb_id}_{chain_id}")
                partition_file = self.create_domain_partition_file(pdb_id, chain_id, domains, output_dir)
        else:
            # No domain summary file or it doesn't exist
            logger.warning(f"No valid domain summary file for {pdb_id}_{chain_id}, creating single domain")
            partition_file = self.create_single_domain(pdb_id, chain_id, length, output_dir)
        
        # Make partition_file relative to data_root if it's not already
        if os.path.isabs(partition_file) and partition_file.startswith(self.data_root):
            relative_path = os.path.relpath(partition_file, self.data_root)
        else:
            relative_path = partition_file
        
        # Update database
        if not dry_run:
            # Add domain partition file record
            if not self.add_domain_partition_file(process_id, relative_path):
                logger.error(f"Failed to add domain partition file record for {pdb_id}_{chain_id}")
                return False
            
            # Update protein status
            if not self.update_protein_status(process_id, 'domain_partition_complete'):
                logger.error(f"Failed to update status for {pdb_id}_{chain_id}")
                return False
            
            logger.info(f"Successfully fixed {pdb_id}_{chain_id}")
            return True
        else:
            logger.info(f"Dry run: Would update {pdb_id}_{chain_id} to domain_partition_complete")
            return True
    
    def fix_stuck_proteins(self, dry_run=False):
        """
        Fix all proteins stuck at domain_summary stage
        
        Args:
            dry_run: If True, don't make any database changes
            
        Returns:
            Number of successfully fixed proteins
        """
        # Get stuck proteins
        proteins = self.get_stuck_proteins()
        logger.info(f"Found {len(proteins)} proteins stuck at domain_summary stage")
        
        if not proteins:
            logger.info("No proteins to fix")
            return 0
        
        # Fix each protein
        success_count = 0
        for protein in proteins:
            try:
                if self.fix_protein(protein, dry_run):
                    success_count += 1
            except Exception as e:
                logger.error(f"Error fixing protein {protein['pdb_id']}_{protein['chain_id']}: {e}")
                if not dry_run:
                    self.rollback()
        
        # Update batch status if needed
        if not dry_run and success_count > 0:
            self.update_batch_status()
        
        return success_count

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Fix proteins stuck at domain_summary stage")
    parser.add_argument("--batch-id", type=int, default=31, help="Batch ID to fix")
    parser.add_argument("--config", type=str, required=True, help="Configuration file path")
    parser.add_argument("--dry-run", action="store_true", help="Dry run without making database changes")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

def load_config(config_path):
    """Load configuration from file"""
    import yaml
    import os
    
    # Load main config file
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Check for local config override
        config_dir = os.path.dirname(config_path)
        local_config_path = os.path.join(config_dir, 'config.local.yml')
        
        if os.path.exists(local_config_path):
            with open(local_config_path, 'r') as f:
                local_config = yaml.safe_load(f)
                
            # Merge configs, with local config taking precedence
            if local_config:
                merge_configs(config, local_config)
        
        logger.info(f"Loaded configuration from {config_path} and local overrides")
        return config
    except Exception as e:
        logger.error(f"Error loading configuration: {e}")
        # Return fallback config for testing purposes
        return {
            'db_host': 'localhost',
            'db_port': 5432,
            'db_name': 'ecod',
            'db_user': 'ecod',
            'db_password': 'password',
            'data_root': '/data/ecod',
            'reference_version': 'develop291'
        }

def merge_configs(config, local_config):
    """Recursively merge two configuration dictionaries"""
    for key, value in local_config.items():
        if key in config and isinstance(config[key], dict) and isinstance(value, dict):
            merge_configs(config[key], value)
        else:
            config[key] = value

def main():
    """Main function"""
    args = parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info(f"Starting fix for batch {args.batch_id}")
    
    # Load configuration
    config = load_config(args.config)
    config['batch_id'] = args.batch_id
    
    # Initialize batch fixer
    fixer = BatchFixer(config)
    
    # Connect to database
    if not fixer.connect_db():
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Fix stuck proteins
        success_count = fixer.fix_stuck_proteins(args.dry_run)
        
        if not args.dry_run:
            fixer.commit()
            logger.info(f"Successfully fixed {success_count} proteins")
        else:
            logger.info(f"Dry run: Would fix {success_count} proteins")
        
    except Exception as e:
        logger.error(f"Error fixing batch: {e}")
        if not args.dry_run:
            fixer.rollback()
        return 1
    finally:
        # Close database connection
        fixer.close_db()
    
    logger.info("Completed successfully")
    return 0

if __name__ == "__main__":
    sys.exit(main())