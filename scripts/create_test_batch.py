#!/usr/bin/env python3
"""
create_test_batch.py - Create a test batch using existing protein data from pdb_analysis

This script creates a new batch in ecod_schema and uses proteins from pdb_analysis
to quickly set up a testing environment without needing full schema migration.
"""

import os
import sys
import argparse
import logging
import datetime
import random
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import application modules
from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.exceptions import ECODError
from ecod.error_handlers import handle_exceptions

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )
    return logging.getLogger("ecod.test_creator")

class TestBatchCreator:
    """Class to handle creating a test batch for PyECOD using existing data"""
    
    def __init__(self, config_path: str):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.db_config = self.config_manager.get_db_config()
        self.config = self.config_manager.config
        self.db = DBManager(self.db_config)
        self.logger = logging.getLogger("ecod.test_creator")
    
    def get_proteins_from_pdb_analysis(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get sample proteins from pdb_analysis schema
        
        Args:
            limit: Number of proteins to retrieve
            
        Returns:
            List of protein data
        """
        query = f"""
        SELECT 
            p.pdb_id, pc.chain_id, 
            CONCAT(p.pdb_id, '_', pc.chain_id) AS source_id,
            pc.sequence_length AS length,
            cs.sequence
        FROM 
            pdb_analysis.pdb_entries p
        JOIN
            pdb_analysis.pdb_chains pc ON p.id = pc.pdb_entry_id
        LEFT JOIN
            pdb_analysis.chain_sequences cs ON pc.id = cs.chain_id
        WHERE 
            cs.sequence IS NOT NULL
        ORDER BY 
            RANDOM()
        LIMIT {limit}
        """
        
        try:
            proteins = self.db.execute_dict_query(query)
            self.logger.info(f"Found {len(proteins)} random proteins in pdb_analysis schema")
            return proteins
        except Exception as e:
            self.logger.error(f"Error retrieving proteins: {str(e)}")
            return []
    
    def create_protein_records(self, protein_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Create protein records in ecod_schema
        
        Args:
            protein_data: List of protein data from pdb_analysis
            
        Returns:
            List of created protein data with IDs
        """
        created_proteins = []
        
        for protein in protein_data:
            pdb_id = protein['pdb_id']
            chain_id = protein['chain_id']
            source_id = f"{pdb_id}_{chain_id}"
            length = protein.get('length', 0)
            sequence = protein.get('sequence', '')
            
            # Check if protein already exists
            check_query = """
            SELECT id FROM ecod_schema.protein WHERE pdb_id = %s AND chain_id = %s
            """
            existing = self.db.execute_query(check_query, (pdb_id, chain_id))
            
            if existing:
                # Protein already exists, get the ID
                protein_id = existing[0][0]
                self.logger.info(f"Found existing protein {source_id} with ID {protein_id}")
            else:
                # Insert protein record
                protein_id = self.db.insert(
                    "ecod_schema.protein",
                    {
                        "pdb_id": pdb_id,
                        "chain_id": chain_id,
                        "source_id": source_id,
                        "length": length
                    },
                    "id"
                )
                
                # Insert sequence record if available
                if sequence:
                    self.db.insert(
                        "ecod_schema.protein_sequence",
                        {
                            "protein_id": protein_id,
                            "sequence": sequence
                        }
                    )
                
                self.logger.info(f"Created protein {source_id} with ID {protein_id}")
            
            # Add to created proteins list
            created_proteins.append({
                "id": protein_id,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "source_id": source_id,
                "length": length,
                "sequence": sequence
            })
            
        return created_proteins
    
    def create_batch(self, proteins: List[Dict[str, Any]], batch_type: str = "test") -> int:
        """Create a new batch containing the specified proteins
        
        Args:
            proteins: List of protein data
            batch_type: Type of batch to create
            
        Returns:
            Batch ID
        """
        # Generate batch name with timestamp
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        batch_name = f"{batch_type}_batch_{timestamp}"
        
        # Get base output directory from config
        base_dir = self.config.get('paths', {}).get('output_dir', './output')
        
        # Create batch directory
        batch_path = os.path.join(base_dir, batch_name)
        os.makedirs(batch_path, exist_ok=True)
        
        # Create subdirectories
        os.makedirs(os.path.join(batch_path, "query_fastas"), exist_ok=True)
        os.makedirs(os.path.join(batch_path, "chain_blast_results"), exist_ok=True)
        os.makedirs(os.path.join(batch_path, "domain_blast_results"), exist_ok=True)
        
        # Insert batch record
        batch_id = self.db.insert(
            "ecod_schema.batch",
            {
                "batch_name": batch_name,
                "base_path": batch_path,
                "type": batch_type,
                "ref_version": self.config.get('reference', {}).get('current_version', 'develop281'),
                "total_items": len(proteins),
                "status": "created"
            },
            "id"
        )
        
        # Register proteins in this batch
        self._register_proteins_in_batch(batch_id, batch_path, proteins)
        
        self.logger.info(f"Created batch {batch_name} with ID {batch_id} containing {len(proteins)} proteins")
        return batch_id
    
    def _register_proteins_in_batch(self, batch_id: int, batch_path: str, proteins: List[Dict[str, Any]]) -> None:
        """Register proteins in a batch and create initial files
        
        Args:
            batch_id: Batch ID
            batch_path: Path to batch directory
            proteins: List of protein data
        """
        fasta_dir = os.path.join(batch_path, "query_fastas")
        
        for protein in proteins:
            # Determine relative path for this protein
            pdb_id = protein['pdb_id']
            chain_id = protein['chain_id']
            source_id = protein['source_id']
            rel_path = f"{pdb_id}_{chain_id}"
            
            # Register in process_status
            process_id = self.db.insert(
                "ecod_schema.process_status",
                {
                    "protein_id": protein['id'],
                    "batch_id": batch_id,
                    "current_stage": "fasta_generated",
                    "status": "pending",
                    "relative_path": rel_path
                },
                "id"
            )
            
            # Generate FASTA file
            fasta_path = os.path.join(fasta_dir, f"{source_id}.fa")
            with open(fasta_path, 'w') as f:
                f.write(f">{source_id}\n{protein['sequence']}\n")
            
            # Register FASTA file
            self.db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": "fasta",
                    "file_path": f"query_fastas/{source_id}.fa",
                    "file_exists": True,
                    "file_size": os.path.getsize(fasta_path)
                }
            )
            
            self.logger.info(f"Registered protein {source_id} in batch {batch_id}")
    
    def link_existing_blast_results(self, source_id: str, process_id: int, 
                                   existing_batch_path: str, dest_batch_path: str) -> bool:
        """Link existing BLAST results for a protein from an existing batch
        
        Args:
            source_id: Protein source ID (e.g., "1k4t_D")
            process_id: Process ID in the new batch
            existing_batch_path: Path to existing batch
            dest_batch_path: Path to new batch
            
        Returns:
            True if successful
        """
        success = False
        
        # Find chainwise blast result
        chain_pattern = f"{source_id}.chainwise_blast.xml"
        chain_blast_path = None
        
        # Look in batch_N subdirectories
        for batch_dir in os.listdir(os.path.join(existing_batch_path, "chain_blast_results")):
            potential_path = os.path.join(existing_batch_path, "chain_blast_results", batch_dir, chain_pattern)
            if os.path.exists(potential_path):
                chain_blast_path = potential_path
                break
        
        # Find domain blast result
        domain_pattern = f"{source_id}.domain_blast.xml"
        domain_blast_path = None
        
        # Look in batch_N subdirectories
        for batch_dir in os.listdir(os.path.join(existing_batch_path, "domain_blast_results")):
            potential_path = os.path.join(existing_batch_path, "domain_blast_results", batch_dir, domain_pattern)
            if os.path.exists(potential_path):
                domain_blast_path = potential_path
                break
        
        # Create destination directories if they don't exist
        chain_result_dir = os.path.join(dest_batch_path, "chain_blast_results")
        domain_result_dir = os.path.join(dest_batch_path, "domain_blast_results")
        
        os.makedirs(chain_result_dir, exist_ok=True)
        os.makedirs(domain_result_dir, exist_ok=True)
        
        # Copy/link chain blast result
        if chain_blast_path:
            dest_chain_path = os.path.join(chain_result_dir, f"{source_id}.blast_chain")
            try:
                # Copy file content
                with open(chain_blast_path, 'rb') as src, open(dest_chain_path, 'wb') as dst:
                    dst.write(src.read())
                
                # Register file
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "chain_blast_result",
                        "file_path": f"chain_blast_results/{source_id}.blast_chain",
                        "file_exists": True,
                        "file_size": os.path.getsize(dest_chain_path)
                    }
                )
                
                self.logger.info(f"Copied chain BLAST result for {source_id}")
                success = True
            except Exception as e:
                self.logger.error(f"Error copying chain BLAST result for {source_id}: {str(e)}")
        
        # Copy/link domain blast result
        if domain_blast_path:
            dest_domain_path = os.path.join(domain_result_dir, f"{source_id}.blast_domain")
            try:
                # Copy file content
                with open(domain_blast_path, 'rb') as src, open(dest_domain_path, 'wb') as dst:
                    dst.write(src.read())
                
                # Register file
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "domain_blast_result",
                        "file_path": f"domain_blast_results/{source_id}.blast_domain",
                        "file_exists": True,
                        "file_size": os.path.getsize(dest_domain_path)
                    }
                )
                
                self.logger.info(f"Copied domain BLAST result for {source_id}")
                success = True
            except Exception as e:
                self.logger.error(f"Error copying domain BLAST result for {source_id}: {str(e)}")
        
        # Update process status if successful
        if success:
            self.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": "blast_complete",
                    "status": "success"
                },
                "id = %s",
                (process_id,)
            )
        
        return success
    
    def copy_existing_results(self, batch_id: int, source_batch_path: str) -> int:
        """Copy existing BLAST results from a source batch to the new batch
        
        Args:
            batch_id: New batch ID
            source_batch_path: Path to source batch with existing results
            
        Returns:
            Number of proteins with copied results
        """
        # Get batch info
        batch_query = """
        SELECT id, base_path FROM ecod_schema.batch WHERE id = %s
        """
        batch_result = self.db.execute_dict_query(batch_query, (batch_id,))
        
        if not batch_result:
            self.logger.error(f"Batch {batch_id} not found")
            return 0
        
        dest_batch_path = batch_result[0]['base_path']
        
        # Get processes in this batch
        process_query = """
        SELECT 
            ps.id as process_id, 
            p.source_id
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        WHERE 
            ps.batch_id = %s
        """
        processes = self.db.execute_dict_query(process_query, (batch_id,))
        
        if not processes:
            self.logger.error(f"No processes found for batch {batch_id}")
            return 0
        
        # Copy results for each protein
        count = 0
        for process in processes:
            source_id = process['source_id']
            process_id = process['process_id']
            
            if self.link_existing_blast_results(source_id, process_id, source_batch_path, dest_batch_path):
                count += 1
        
        # Update batch status with completed items
        self.db.update(
            "ecod_schema.batch",
            {
                "completed_items": count,
                "status": "processing" if count > 0 else "created"
            },
            "id = %s",
            (batch_id,)
        )
        
        return count

@handle_exceptions(exit_on_error=True)
def main():
    parser = argparse.ArgumentParser(description='Create a test batch for PyECOD using existing data')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--output-dir', type=str, 
                      help='Output directory for batch data (overrides config)')
    parser.add_argument('--protein-count', type=int, default=5,
                      help='Number of sample proteins to create')
    parser.add_argument('--batch-type', type=str, default='test',
                      help='Type of batch to create')
    parser.add_argument('--source-batch-path', type=str,
                      help='Path to existing batch with BLAST results to copy')
    parser.add_argument('--protein-id', type=str,
                      help='Specific protein ID (e.g. 1k4t_D) to include')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    logger = setup_logging(args.verbose)
    
    # Initialize creator
    creator = TestBatchCreator(args.config)
    
    # Override output directory if specified
    if args.output_dir:
        creator.config['paths'] = creator.config.get('paths', {})
        creator.config['paths']['output_dir'] = args.output_dir
    
    proteins = []
    
    # If specific protein ID requested
    if args.protein_id:
        # Parse PDB ID and chain ID
        parts = args.protein_id.split('_')
        if len(parts) != 2:
            logger.error(f"Invalid protein ID format: {args.protein_id}. Expected format: PDBID_CHAIN")
            return 1
        
        pdb_id, chain_id = parts
        
        # Query for this specific protein
        query = f"""
        SELECT 
            p.pdb_id, pc.chain_id, 
            CONCAT(p.pdb_id, '_', pc.chain_id) AS source_id,
            pc.sequence_length AS length,
            cs.sequence
        FROM 
            pdb_analysis.pdb_entries p
        JOIN
            pdb_analysis.pdb_chains pc ON p.id = pc.pdb_entry_id
        LEFT JOIN
            pdb_analysis.chain_sequences cs ON pc.id = cs.chain_id
        WHERE 
            p.pdb_id = '{pdb_id}' AND pc.chain_id = '{chain_id}'
            AND cs.sequence IS NOT NULL
        """
        
        specific_proteins = creator.db.execute_dict_query(query)
        if specific_proteins:
            proteins.extend(specific_proteins)
            logger.info(f"Found specific protein: {args.protein_id}")
        else:
            logger.warning(f"Protein {args.protein_id} not found or has no sequence")
    
    # Get random proteins if needed to reach requested count
    if len(proteins) < args.protein_count:
        more_needed = args.protein_count - len(proteins)
        random_proteins = creator.get_proteins_from_pdb_analysis(more_needed)
        proteins.extend(random_proteins)
    
    if not proteins:
        logger.error("No proteins found to create batch")
        print("Error: No proteins found to create batch")
        return 1
    
    # Create protein records
    created_proteins = creator.create_protein_records(proteins)
    print(f"Created/found {len(created_proteins)} protein records")
    
    # Create batch
    batch_id = creator.create_batch(created_proteins, args.batch_type)
    print(f"Created batch with ID: {batch_id}")
    
    # Copy existing BLAST results if source batch specified
    if args.source_batch_path:
        if os.path.exists(args.source_batch_path):
            count = creator.copy_existing_results(batch_id, args.source_batch_path)
            print(f"Copied BLAST results for {count} proteins from {args.source_batch_path}")
        else:
            logger.error(f"Source batch path not found: {args.source_batch_path}")
            print(f"Warning: Source batch path not found: {args.source_batch_path}")
    
    print("Test batch creation complete!")
    print(f"You can now run the PyECOD pipeline with: --batch-id {batch_id}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())