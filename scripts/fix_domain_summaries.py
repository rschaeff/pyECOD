#!/usr/bin/env python3
"""
fix_domain_summaries.py - Repair incorrectly formatted domain summary files
"""

import os
import sys
import logging
import argparse
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext

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

def get_blast_file_paths(context, protein_id: int, batch_path: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Get the chain and domain BLAST file paths for a protein from the database
    
    Args:
        context: Application context
        protein_id: Protein ID
        batch_path: Base path for the batch
        
    Returns:
        Tuple of (chain_blast_path, domain_blast_path)
    """
    logger = logging.getLogger("ecod.fix_summaries")
    
    query = """
    SELECT 
        pf.file_type, pf.file_path
    FROM 
        ecod_schema.process_status ps
    JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE 
        ps.protein_id = %s
        AND pf.file_type IN ('chain_blast_result', 'domain_blast_result')
        AND pf.file_exists = TRUE
    """
    
    try:
        results = context.db.execute_query(query, (protein_id,))
        
        chain_blast_path = None
        domain_blast_path = None
        
        for file_type, file_path in results:
            # Convert to absolute path if needed
            abs_path = os.path.join(batch_path, file_path) if not os.path.isabs(file_path) else file_path
            
            if file_type == 'chain_blast_result':
                chain_blast_path = abs_path
            elif file_type == 'domain_blast_result':
                domain_blast_path = abs_path
        
        return chain_blast_path, domain_blast_path
        
    except Exception as e:
        logger.error(f"Error getting BLAST file paths: {e}")
        return None, None

def fix_summary_file(file_info: Tuple, batch_path: str, context) -> bool:
    """
    Fix an incorrectly formatted domain summary file
    
    Args:
        file_info: Tuple of (protein_id, pdb_id, chain_id, file_path, file_id)
        batch_path: Base path for the batch
        context: Application context
        
    Returns:
        True if fixed successfully, False otherwise
    """
    logger = logging.getLogger("ecod.fix_summaries")
    
    protein_id, pdb_id, chain_id, file_path, file_id = file_info
    
    # Get absolute path if needed
    abs_path = os.path.join(batch_path, file_path) if not os.path.isabs(file_path) else file_path
    
    if not os.path.exists(abs_path):
        logger.error(f"File not found: {abs_path}")
        return False
    
    try:
        # Check if file already has correct structure
        try:
            tree = ET.parse(abs_path)
            root = tree.getroot()
            
            if root.tag == "blast_summ_doc":
                logger.debug(f"File already has correct structure: {abs_path}")
                
                # Still update file location in database to use new structure
                domains_dir = os.path.join(batch_path, "domains")
                os.makedirs(domains_dir, exist_ok=True)
                
                pdb_chain = f"{pdb_id}_{chain_id}"
                new_path = os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml")
                
                # Copy to new location
                try:
                    os.makedirs(os.path.dirname(new_path), exist_ok=True)
                    import shutil
                    shutil.copy2(abs_path, new_path)
                    logger.info(f"Copied file to proper location: {new_path}")
                    
                    # Update database
                    rel_path = os.path.relpath(new_path, batch_path)
                    update_query = """
                    UPDATE ecod_schema.process_file
                    SET file_path = %s
                    WHERE id = %s
                    """
                    context.db.execute_query(update_query, (rel_path, file_id))
                    
                    return True
                except Exception as e:
                    logger.error(f"Error copying file to new location: {e}")
                    return False
                
        except Exception:
            logger.debug(f"Error parsing existing file: {abs_path} - will fix structure")
        
        # Get BLAST file paths from database
        chain_blast_file, domain_blast_file = get_blast_file_paths(context, protein_id, batch_path)
        
        # Create a new domain summary XML
        root = ET.Element("blast_summ_doc")
        blast_summ = ET.SubElement(root, "blast_summ")
        blast_summ.set("pdb", pdb_id)
        blast_summ.set("chain", chain_id)
        
        # Process chain blast if available
        if not chain_blast_file or not os.path.exists(chain_blast_file):
            blast_summ.set("no_chain_blast", "true")
            logger.warning(f"No chain blast file found for {pdb_id}_{chain_id}")
        else:
            try:
                chain_tree = ET.parse(chain_blast_file)
                chain_root = chain_tree.getroot()
                
                # Check if file has hits
                hits = chain_root.findall(".//Hit")
                if not hits:
                    blast_summ.set("chain_blast_no_hits", "true")
                    logger.debug(f"Chain blast file has no hits: {chain_blast_file}")
                else:
                    # Process the chain BLAST properly
                    chain_blast_run = ET.SubElement(blast_summ, "chain_blast_run")
                    
                    # Add program and version
                    program = chain_root.findtext(".//BlastOutput_program", "")
                    chain_blast_run.set("program", program)
                    
                    version = chain_root.findtext(".//BlastOutput_version", "")
                    chain_blast_run.set("version", version)
                    
                    # Add database
                    db = chain_root.findtext(".//BlastOutput_db", "")
                    db_node = ET.SubElement(chain_blast_run, "blast_db")
                    db_node.text = db
                    
                    # Add query info
                    query = chain_root.findtext(".//BlastOutput_query-def", "")
                    query_node = ET.SubElement(chain_blast_run, "blast_query")
                    query_node.text = query
                    
                    query_len = chain_root.findtext(".//BlastOutput_query-len", "")
                    query_len_node = ET.SubElement(chain_blast_run, "query_len")
                    query_len_node.text = query_len
                    
                    # Add hits section
                    hits_node = ET.SubElement(chain_blast_run, "hits")
                    
                    # Add a limited number of hits (max 20)
                    for hit_idx, hit in enumerate(hits[:20]):
                        hit_num = hit.findtext("Hit_num", "")
                        hit_def = hit.findtext("Hit_def", "")
                        
                        # Parse PDB ID and chain from hit definition
                        hit_pdb = "NA"
                        hit_chain = "NA"
                        if " " in hit_def:
                            parts = hit_def.split()
                            if len(parts) >= 2:
                                hit_pdb = parts[0]
                                hit_chain = parts[1]
                        
                        # Create hit element with basic attributes
                        hit_elem = ET.SubElement(hits_node, "hit")
                        hit_elem.set("num", hit_num)
                        hit_elem.set("pdb_id", hit_pdb)
                        hit_elem.set("chain_id", hit_chain)
                        
                        # Add simplified hit info for HSPs
                        hsps = hit.findall(".//Hsp")
                        if hsps:
                            # Get regions from first HSP
                            hsp = hsps[0]
                            hsp_query_from = hsp.findtext("Hsp_query-from", "")
                            hsp_query_to = hsp.findtext("Hsp_query-to", "")
                            hsp_hit_from = hsp.findtext("Hsp_hit-from", "")
                            hsp_hit_to = hsp.findtext("Hsp_hit-to", "")
                            
                            # Add regions
                            query_reg = ET.SubElement(hit_elem, "query_reg")
                            query_reg.text = f"{hsp_query_from}-{hsp_query_to}"
                            
                            hit_reg = ET.SubElement(hit_elem, "hit_reg")
                            hit_reg.text = f"{hsp_hit_from}-{hsp_hit_to}"
            except Exception as e:
                logger.error(f"Error processing chain blast file: {e}")
                blast_summ.set("chain_blast_error", "true")
        
        # Process domain blast if available
        if not domain_blast_file or not os.path.exists(domain_blast_file):
            blast_summ.set("no_domain_blast", "true")
            logger.warning(f"No domain blast file found for {pdb_id}_{chain_id}")
        else:
            try:
                domain_tree = ET.parse(domain_blast_file)
                domain_root = domain_tree.getroot()
                
                # Check if file has hits
                hits = domain_root.findall(".//Hit")
                if not hits:
                    blast_summ.set("domain_blast_no_hits", "true")
                    logger.debug(f"Domain blast file has no hits: {domain_blast_file}")
                else:
                    # Process the domain BLAST properly
                    blast_run = ET.SubElement(blast_summ, "blast_run")
                    
                    # Add program and version
                    program = domain_root.findtext(".//BlastOutput_program", "")
                    blast_run.set("program", program)
                    
                    version = domain_root.findtext(".//BlastOutput_version", "")
                    blast_run.set("version", version)
                    
                    # Add database
                    db = domain_root.findtext(".//BlastOutput_db", "")
                    db_node = ET.SubElement(blast_run, "blast_db")
                    db_node.text = db
                    
                    # Add query info
                    query = domain_root.findtext(".//BlastOutput_query-def", "")
                    query_node = ET.SubElement(blast_run, "blast_query")
                    query_node.text = query
                    
                    query_len = domain_root.findtext(".//BlastOutput_query-len", "")
                    query_len_node = ET.SubElement(blast_run, "query_len")
                    query_len_node.text = query_len
                    
                    # Add hits section
                    hits_node = ET.SubElement(blast_run, "hits")
                    
                    # Add a limited number of hits (max 20)
                    for hit_idx, hit in enumerate(hits[:20]):
                        hit_num = hit.findtext("Hit_num", "")
                        hit_def = hit.findtext("Hit_def", "")
                        
                        # Try to parse domain ID from hit definition
                        import re
                        domain_id = "NA"
                        hit_pdb = "NA"
                        hit_chain = "NA"
                        
                        domain_match = re.search(r"((d|g|e)(\d\w{3})\w+\d+)\s+(\w+):", hit_def)
                        if domain_match:
                            domain_id = domain_match.group(1)
                            hit_pdb = domain_match.group(3)
                            hit_chain = domain_match.group(4)
                        
                        # Create hit element
                        hit_elem = ET.SubElement(hits_node, "hit")
                        hit_elem.set("num", hit_num)
                        hit_elem.set("domain_id", domain_id)
                        hit_elem.set("pdb_id", hit_pdb)
                        hit_elem.set("chain_id", hit_chain)
                        
                        # Add simplified hit info for HSPs
                        hsps = hit.findall(".//Hsp")
                        if hsps:
                            # Get e-values
                            evalues = []
                            for hsp in hsps[:3]:  # Limit to first 3 HSPs
                                evalues.append(hsp.findtext("Hsp_evalue", ""))
                            
                            hit_elem.set("evalues", ",".join(evalues))
                            
                            # Get regions from first HSP
                            hsp = hsps[0]
                            hsp_query_from = hsp.findtext("Hsp_query-from", "")
                            hsp_query_to = hsp.findtext("Hsp_query-to", "")
                            hsp_hit_from = hsp.findtext("Hsp_hit-from", "")
                            hsp_hit_to = hsp.findtext("Hsp_hit-to", "")
                            
                            # Add regions
                            query_reg = ET.SubElement(hit_elem, "query_reg")
                            query_reg.text = f"{hsp_query_from}-{hsp_query_to}"
                            
                            hit_reg = ET.SubElement(hit_elem, "hit_reg")
                            hit_reg.text = f"{hsp_hit_from}-{hsp_hit_to}"
            except Exception as e:
                logger.error(f"Error processing domain blast file: {e}")
                blast_summ.set("domain_blast_error", "true")
        
        # Add self-comparison note if missing
        blast_summ.set("no_selfcomp", "true")
        
        # Create backup of original file
        backup_path = abs_path + ".bak"
        try:
            os.rename(abs_path, backup_path)
            logger.debug(f"Created backup: {backup_path}")
        except Exception as e:
            logger.warning(f"Could not create backup: {e}")
        
        # Prepare new file path in proper directory structure
        domains_dir = os.path.join(batch_path, "domains")
        os.makedirs(domains_dir, exist_ok=True)
        
        pdb_chain = f"{pdb_id}_{chain_id}"
        new_path = os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml")
        
        # Write the fixed file directly to the new location
        tree = ET.ElementTree(root)
        tree.write(new_path, encoding='utf-8', xml_declaration=True)
        logger.info(f"Created fixed domain summary: {new_path}")
        
        # Update database record to point to the new location
        rel_path = os.path.relpath(new_path, batch_path)
        update_query = """
        UPDATE ecod_schema.process_file
        SET file_path = %s
        WHERE id = %s
        """
        context.db.execute_query(update_query, (rel_path, file_id))
        logger.debug(f"Updated database record to point to: {rel_path}")
        
        return True
        
    except Exception as e:
        logger.error(f"Error fixing summary file for {pdb_id}_{chain_id}: {str(e)}", exc_info=True)
        return False

def fix_batch_summaries(batch_id: int, config_path: str, 
                      force: bool = False, 
                      limit: Optional[int] = None, 
                      dry_run: bool = False):
    """Fix all domain summaries in a batch"""
    context = ApplicationContext(config_path)
    logger = logging.getLogger("ecod.fix_summaries")
    
    # Get batch path
    batch_query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
    batch_result = context.db.execute_query(batch_query, (batch_id,))
    
    if not batch_result:
        logger.error(f"Batch {batch_id} not found")
        return 1
        
    batch_path = batch_result[0][0]
    
    # Create domains directory in the batch path
    domains_dir = os.path.join(batch_path, "domains")
    if not dry_run:
        os.makedirs(domains_dir, exist_ok=True)
    
    # Query for all domain summary files in the batch
    query = """
    SELECT 
        p.id, p.pdb_id, p.chain_id, pf.file_path, pf.id as file_id
    FROM 
        ecod_schema.protein p
    JOIN 
        ecod_schema.process_status ps ON p.id = ps.protein_id
    JOIN 
        ecod_schema.process_file pf ON ps.id = pf.process_id
    WHERE 
        ps.batch_id = %s 
        AND pf.file_type = 'domain_summary'
        AND pf.file_exists = TRUE
    ORDER BY p.pdb_id, p.chain_id
    """
    
    if limit:
        query += f" LIMIT {limit}"
    
    summaries = context.db.execute_query(query, (batch_id,))
    
    if not summaries:
        logger.error(f"No domain summaries found for batch {batch_id}")
        return 1
    
    logger.info(f"Found {len(summaries)} domain summaries to fix")
    
    if dry_run:
        logger.info("Dry run - no files will be modified")
        return 0
    
    # Process all summaries
    fixed_count = 0
    error_count = 0
    
    for i, summary in enumerate(summaries):
        # Fix the file
        if fix_summary_file(summary, batch_path, context):
            fixed_count += 1
        else:
            error_count += 1
        
        # Progress update
        if (i + 1) % 100 == 0 or (i + 1) == len(summaries):
            logger.info(f"Processed {i+1}/{len(summaries)} files. Fixed: {fixed_count}, Errors: {error_count}")
    
    logger.info(f"Summary fix complete: {fixed_count} fixed, {error_count} errors")
    
    # If requested, clean up the old directory structure
    if fixed_count > 0 and not dry_run:
        logger.info("Would you like to clean up the old directory structure? (y/n)")
        response = input().strip().lower()
        
        if response == 'y':
            logger.info("Cleaning up old directory structure...")
            for protein_id, pdb_id, chain_id, file_path, file_id in summaries:
                old_dir = os.path.join(batch_path, pdb_id, chain_id)
                if os.path.exists(old_dir):
                    try:
                        import shutil
                        shutil.rmtree(old_dir)
                        logger.info(f"Removed old directory: {old_dir}")
                    except Exception as e:
                        logger.error(f"Error removing directory {old_dir}: {e}")
    
    return 0

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Fix incorrectly formatted domain summary files')
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to fix')
    parser.add_argument('--limit', type=int,
                      help='Maximum number of files to fix')
    parser.add_argument('--force', action='store_true',
                      help='Fix files even if they appear to have correct structure')
    parser.add_argument('--dry-run', action='store_true',
                      help='Don\'t modify any files, just report')
    parser.add_argument('--cleanup', action='store_true',
                      help='Clean up old directory structure after fixing')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.fix_summaries")
    
    # Fix batch summaries
    return fix_batch_summaries(
        args.batch_id,
        args.config,
        args.force,
        args.limit,
        args.dry_run
    )

if __name__ == "__main__":
    sys.exit(main())