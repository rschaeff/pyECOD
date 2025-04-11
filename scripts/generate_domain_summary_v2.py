#!/usr/bin/env python3
"""
generate_domain_summary.py - Generate domain summary for a single protein
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from typing import Dict, Any, Optional

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.pipelines.domain_analysis.summary import DomainSummary

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
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

def fix_domain_summary(summary_processor):
    """
    Fix issues in the DomainSummary class
    """
    # Get the logger
    logger = logging.getLogger("ecod.fix_domain_summary")
    
    # Define a new create_summary method to fix the issues
    def fixed_create_summary(self, pdb_id, chain_id, reference, job_dump_dir, blast_only=False):
        """Fixed version of create_summary method"""
        logger.info(f"Creating summary for {pdb_id}_{chain_id} (reference: {reference})")
        
        # Define paths 
        pdb_chain = f"{pdb_id}_{chain_id}"
        chain_dir = os.path.join(job_dump_dir, pdb_chain)
        
        # Create output directory if it doesn't exist
        os.makedirs(chain_dir, exist_ok=True)
        
        # Define output file name
        summary_xml_file = (f"{pdb_chain}.{reference}.blast_summ.blast_only.xml" 
                          if blast_only else f"{pdb_chain}.{reference}.blast_summ.xml")
        full_output_path = os.path.join(chain_dir, summary_xml_file)
        
        if os.path.exists(full_output_path) and not self.config.get('force_overwrite', False):
            logger.info(f"Output file {full_output_path} already exists, skipping...")
            return full_output_path
        
        # Check if this is a peptide by reading the FASTA file
        fasta_path = os.path.join(chain_dir, f"{pdb_chain}.fa")
        sequence = None
        if os.path.exists(fasta_path):
            try:
                with open(fasta_path, 'r') as f:
                    lines = f.readlines()
                    if len(lines) > 1:
                        sequence = ''.join(lines[1:]).strip()
            except Exception as e:
                logger.warning(f"Error reading FASTA file: {e}")
        
        # Check for database-recorded FASTA file if not found
        if not sequence:
            logger.info("FASTA file not found in expected location, checking database")
            db_config = self.config_manager.get_db_config()
            from ecod.db.manager import DBManager
            db = DBManager(db_config)
            
            query = """
            SELECT pf.file_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            AND pf.file_type = 'fasta'
            AND pf.file_exists = TRUE
            ORDER BY pf.id DESC
            LIMIT 1
            """
            
            try:
                rows = db.execute_query(query, (pdb_id, chain_id))
                if rows:
                    db_fasta_path = rows[0][0]
                    full_fasta_path = os.path.join(job_dump_dir, db_fasta_path)
                    if os.path.exists(full_fasta_path):
                        logger.info(f"Found FASTA file in database: {full_fasta_path}")
                        with open(full_fasta_path, 'r') as f:
                            lines = f.readlines()
                            if len(lines) > 1:
                                sequence = ''.join(lines[1:]).strip()
            except Exception as e:
                logger.warning(f"Error querying database for FASTA file: {e}")
        
        # If sequence is found and is a peptide, create a special summary
        if sequence and len(sequence) < 30:
            logger.info(f"Sequence for {pdb_id}_{chain_id} is a peptide with length {len(sequence)}")
            
            # Create a special summary for peptides
            peptide_summary = ET.Element("blast_summ_doc")
            blast_summ = ET.SubElement(peptide_summary, "blast_summ")
            blast_summ.set("pdb", pdb_id)
            blast_summ.set("chain", chain_id)
            blast_summ.set("is_peptide", "true")
            blast_summ.set("sequence_length", str(len(sequence)))
            
            # Write output file
            os.makedirs(os.path.dirname(full_output_path), exist_ok=True)
            tree = ET.ElementTree(peptide_summary)
            tree.write(full_output_path, encoding='utf-8', xml_declaration=True)
            
            logger.info(f"Created peptide summary: {full_output_path}")
            return full_output_path
        
        # Create XML document root for regular proteins
        root = ET.Element("blast_summ_doc")
        blast_summ = ET.SubElement(root, "blast_summ")
        blast_summ.set("pdb", pdb_id)
        blast_summ.set("chain", chain_id)
        
        # Process self-comparison results
        self_comp_path = os.path.join(chain_dir, f"{pdb_chain}.self_comp.xml")
        if not os.path.exists(self_comp_path):
            logger.warning(f"No self comparison results for {pdb_id} {chain_id}")
            blast_summ.set("no_selfcomp", "true")
        else:
            # Use the original self-comparison processing
            self._process_self_comparison(self_comp_path, blast_summ)
        
        # Find and process chain blast results
        chain_blast_path = None
        
        # First try to find in database with 'blast_result' type (most common)
        db_config = self.config_manager.get_db_config()
        from ecod.db.manager import DBManager
        db = DBManager(db_config)
        
        query = """
        SELECT pf.file_path
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        AND pf.file_type = 'blast_result'
        AND pf.file_exists = TRUE
        ORDER BY pf.id DESC
        LIMIT 1
        """
        
        try:
            rows = db.execute_query(query, (pdb_id, chain_id))
            if rows:
                db_file_path = rows[0][0]
                full_path = os.path.join(job_dump_dir, db_file_path)
                full_path = os.path.normpath(full_path)
                
                if os.path.exists(full_path):
                    logger.info(f"Found chain blast file from database: {full_path}")
                    chain_blast_path = full_path
        except Exception as e:
            logger.error(f"Error querying database for chain blast file: {e}")
        
        # If not found, try some common patterns
        if not chain_blast_path:
            logger.info("Chain blast file not found in database, checking common locations")
            common_patterns = [
                os.path.join(job_dump_dir, "blast", "chain", f"{pdb_chain}.chainwise_blast.xml"),
                os.path.join(job_dump_dir, "blast", "chain", "batch_0", f"{pdb_chain}.chainwise_blast.xml"),
                os.path.join(job_dump_dir, "blast", "chain", "batch_1", f"{pdb_chain}.chainwise_blast.xml"),
                os.path.join(job_dump_dir, "blast", "chain", "batch_2", f"{pdb_chain}.chainwise_blast.xml"),
                os.path.join(job_dump_dir, f"{pdb_chain}", f"{pdb_chain}.chainwise_blast.xml")
            ]
            
            for pattern in common_patterns:
                if os.path.exists(pattern):
                    logger.info(f"Found chain blast file at common location: {pattern}")
                    chain_blast_path = pattern
                    break
        
        # Process chain blast if found
        if chain_blast_path:
            try:
                # Read file content
                with open(chain_blast_path, 'r') as f:
                    chain_blast_content = f.read()
                
                # Process chain blast
                if chain_blast_content:
                    self._process_chain_blast(chain_blast_path, blast_summ)
            except Exception as e:
                logger.error(f"Error processing chain blast file: {e}")
                blast_summ.set("chain_blast_error", "true")
        else:
            logger.error(f"No chain blast result file for {pdb_id} {chain_id}")
            blast_summ.set("no_chain_blast", "true")
        
        # Find and process domain blast results
        domain_blast_path = None
        
        # Try to find in database with 'domain_blast_result' type
        query = """
        SELECT pf.file_path
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        AND pf.file_type = 'domain_blast_result'
        AND pf.file_exists = TRUE
        ORDER BY pf.id DESC
        LIMIT 1
        """
        
        try:
            rows = db.execute_query(query, (pdb_id, chain_id))
            if rows:
                db_file_path = rows[0][0]
                full_path = os.path.join(job_dump_dir, db_file_path)
                full_path = os.path.normpath(full_path)
                
                if os.path.exists(full_path):
                    logger.info(f"Found domain blast file from database: {full_path}")
                    domain_blast_path = full_path
        except Exception as e:
            logger.error(f"Error querying database for domain blast file: {e}")
        
        # If not found, try some common patterns
        if not domain_blast_path:
            logger.info("Domain blast file not found in database, checking common locations")
            common_patterns = [
                os.path.join(job_dump_dir, "blast", "domain", f"{pdb_chain}.domain_blast.xml"),
                os.path.join(job_dump_dir, "blast", "domain", "batch_0", f"{pdb_chain}.domain_blast.xml"),
                os.path.join(job_dump_dir, "blast", "domain", "batch_1", f"{pdb_chain}.domain_blast.xml"),
                os.path.join(job_dump_dir, "blast", "domain", "batch_2", f"{pdb_chain}.domain_blast.xml"),
                os.path.join(job_dump_dir, f"{pdb_chain}", f"{pdb_chain}.domain_blast.xml")
            ]
            
            for pattern in common_patterns:
                if os.path.exists(pattern):
                    logger.info(f"Found domain blast file at common location: {pattern}")
                    domain_blast_path = pattern
                    break
        
        # Process domain blast if found
        if domain_blast_path:
            try:
                # Process domain blast
                self._process_blast(domain_blast_path, blast_summ)
            except Exception as e:
                logger.error(f"Error processing domain blast file: {e}")
                blast_summ.set("domain_blast_error", "true")
        else:
            logger.error(f"No domain blast result file for {pdb_id} {chain_id}")
            blast_summ.set("no_domain_blast", "true")
        
        # Process HHSearch results (skip if blast_only mode)
        if not blast_only:
            hhsearch_path = os.path.join(chain_dir, f"{pdb_chain}.{reference}.hh_summ.xml")
            if not os.path.exists(hhsearch_path):
                logger.warning(f"No hhsearch result file for {reference} {pdb_id} {chain_id}")
                blast_summ.set("no_hhsearch", "true")
            else:
                self._process_hhsearch(hhsearch_path, blast_summ)
        
        # Write output file
        os.makedirs(os.path.dirname(full_output_path), exist_ok=True)
        tree = ET.ElementTree(root)
        tree.write(full_output_path, encoding='utf-8', xml_declaration=True)
        
        logger.info(f"Created domain summary: {full_output_path}")
        return full_output_path
    
    # Replace the method
    summary_processor.create_summary = lambda pdb_id, chain_id, reference, job_dump_dir, blast_only=False: fixed_create_summary(summary_processor, pdb_id, chain_id, reference, job_dump_dir, blast_only)
    
    logger.info("Applied fixes to domain summary generation")
    return summary_processor

def main():
    """Main function to generate domain summary for a single protein"""
    parser = argparse.ArgumentParser(description='Generate domain summary for a single protein')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID containing the protein')
    parser.add_argument('--protein-id', type=int, required=True,
                      help='Protein ID to process')
    parser.add_argument('--blast-only', action='store_true',
                      help='Generate blast-only summary (no HHSearch)')
    parser.add_argument('--output-dir', type=str,
                      help='Override output directory (default: from batch path)')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.generate_summary")
    
    # Initialize application context
    context = ApplicationContext(args.config)
    
    # Get protein and batch information
    query = """
    SELECT 
        p.id, p.source_id, p.pdb_id, p.chain_id,
        b.id as batch_id, b.batch_name, b.base_path, b.ref_version
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.batch b ON ps.batch_id = b.id
    WHERE 
        p.id = %s AND b.id = %s
    """
    
    result = context.db.execute_query(query, (args.protein_id, args.batch_id))
    
    if not result:
        logger.error(f"Protein {args.protein_id} not found in batch {args.batch_id}")
        return 1
    
    protein_info = result[0]
    pdb_id = protein_info[2]
    chain_id = protein_info[3]
    batch_path = protein_info[6]
    reference = protein_info[7]
    
    logger.info(f"Generating domain summary for protein: {pdb_id}_{chain_id} (ID: {args.protein_id})")
    
    # Initialize domain summary processor
    summary_processor = DomainSummary(args.config)
    
    # Apply comprehensive fix
    summary_processor = fix_domain_summary(summary_processor)
    
    # Determine output directory
    output_dir = args.output_dir if args.output_dir else batch_path
    
    # Generate domain summary
    try:
        summary_file = summary_processor.create_summary(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference,
            job_dump_dir=output_dir,
            blast_only=args.blast_only
        )
        
        if summary_file:
            logger.info(f"Successfully generated domain summary: {summary_file}")
            
            # Register the summary file in the database
            process_query = """
            SELECT id FROM ecod_schema.process_status
            WHERE protein_id = %s AND batch_id = %s
            """
            process_result = context.db.execute_query(process_query, (args.protein_id, args.batch_id))
            
            if process_result:
                process_id = process_result[0][0]
                
                # Check if summary file record already exists
                check_query = """
                SELECT id FROM ecod_schema.process_file
                WHERE process_id = %s AND file_type = 'domain_summary'
                """
                
                existing = context.db.execute_query(check_query, (process_id,))
                
                if existing:
                    # Update existing record
                    update_query = """
                    UPDATE ecod_schema.process_file
                    SET file_path = %s, file_exists = TRUE, file_size = %s, last_checked = NOW()
                    WHERE id = %s
                    """
                    
                    rel_path = os.path.relpath(summary_file, batch_path)
                    file_size = os.path.getsize(summary_file)
                    
                    context.db.execute_query(update_query, (rel_path, file_size, existing[0][0]))
                    logger.info(f"Updated domain_summary file record in database")
                else:
                    # Create new record
                    insert_query = """
                    INSERT INTO ecod_schema.process_file
                    (process_id, file_type, file_path, file_exists, file_size, last_checked)
                    VALUES (%s, %s, %s, %s, %s, NOW())
                    """
                    
                    rel_path = os.path.relpath(summary_file, batch_path)
                    file_size = os.path.getsize(summary_file)
                    
                    context.db.execute_query(insert_query, (process_id, 'domain_summary', rel_path, True, file_size))
                    logger.info(f"Added domain_summary file record to database")
                
                # Update process status
                status_query = """
                UPDATE ecod_schema.process_status
                SET current_stage = 'domain_summary', status = 'success', updated_at = NOW()
                WHERE id = %s
                """
                
                context.db.execute_query(status_query, (process_id,))
                logger.info(f"Updated process status to domain_summary:success")
            
            return 0
        else:
            logger.error(f"Failed to generate domain summary")
            return 1
            
    except Exception as e:
        logger.error(f"Error generating domain summary: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())