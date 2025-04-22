#!/usr/bin/env python3
"""
process_hhsearch_batch.py - Register and process HHSearch results for a batch
"""

import os
import sys
import argparse
import logging
import glob
from datetime import datetime

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.core.context import ApplicationContext
from ecod.pipelines.hhsearch import HHRToXMLConverter, HHSearchProcessor
from ecod.utils.hhsearch_utils import HHRParser

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

def calculate_range(alignment, start_pos):
    """Calculate range from alignment"""
    ranges = []
    current_range_start = None
    current_pos = start_pos
    
    for char in alignment:
        if char != '-':  # Not a gap
            if current_range_start is None:
                current_range_start = current_pos
            current_pos += 1
        else:  # Gap
            if current_range_start is not None:
                ranges.append(f"{current_range_start}-{current_pos-1}")
                current_range_start = None
    
    # Add the last range if exists
    if current_range_start is not None:
        ranges.append(f"{current_range_start}-{current_pos-1}")
    
    return ",".join(ranges)

def collate_evidence(chain_blast_xml, domain_blast_xml, hhsearch_xml, pdb_id, chain_id, ref_version):
    """Collate evidence from different sources to create domain summary"""
    logger = logging.getLogger("evidence_collator")
    
    try:
        # Parse XML files if they exist
        chain_blast_data = None
        domain_blast_data = None
        hhsearch_data = None
        
        if chain_blast_xml and os.path.exists(chain_blast_xml):
            try:
                chain_blast_data = ET.parse(chain_blast_xml)
            except Exception as e:
                logger.warning(f"Error parsing chain BLAST XML: {str(e)}")
        
        if domain_blast_xml and os.path.exists(domain_blast_xml):
            try:
                domain_blast_data = ET.parse(domain_blast_xml)
            except Exception as e:
                logger.warning(f"Error parsing domain BLAST XML: {str(e)}")
        
        if hhsearch_xml and os.path.exists(hhsearch_xml):
            try:
                hhsearch_data = ET.parse(hhsearch_xml)
            except Exception as e:
                logger.warning(f"Error parsing HHSearch XML: {str(e)}")
        
        # Create root element
        root = ET.Element("domain_summ_doc")
        
        # Add metadata
        metadata = ET.SubElement(root, "metadata")
        ET.SubElement(metadata, "pdb_id").text = pdb_id
        ET.SubElement(metadata, "chain_id").text = chain_id
        ET.SubElement(metadata, "reference").text = ref_version
        ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Add chain blast evidence if available
        chain_blast_elem = ET.SubElement(root, "chain_blast_evidence")
        if chain_blast_data is not None:
            blast_doc = chain_blast_data.find("blast_summ_doc")
            if blast_doc is not None:
                for child in blast_doc:
                    chain_blast_elem.append(ET.fromstring(ET.tostring(child)))
        
        # Add domain blast evidence if available
        domain_blast_elem = ET.SubElement(root, "domain_blast_evidence")
        if domain_blast_data is not None:
            blast_doc = domain_blast_data.find("blast_summ_doc")
            if blast_doc is not None:
                for child in blast_doc:
                    domain_blast_elem.append(ET.fromstring(ET.tostring(child)))
        
        # Add HHSearch evidence if available
        hhsearch_elem = ET.SubElement(root, "hhsearch_evidence")
        if hhsearch_data is not None:
            hh_doc = hhsearch_data.find("hh_summ_doc")
            if hh_doc is not None:
                for child in hh_doc:
                    hhsearch_elem.append(ET.fromstring(ET.tostring(child)))
        
        # Add stub for domain suggestions (placeholder)
        ET.SubElement(root, "domain_suggestions")
        
        # Convert to string
        rough_string = ET.tostring(root, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")
        
        return pretty_xml
        
    except Exception as e:
        logger.error(f"Error collating evidence: {str(e)}")
        return None

def process_batch(batch_id, limit=None, force=False, config_path=None):
    """Process HHSearch results for a batch"""
    logger = logging.getLogger("batch_processor")
    
    # Initialize application context
    context = ApplicationContext(config_path)
    
    # Create parser and converter
    parser = HHRParser(logger)
    converter = HHRToXMLConverter(logger)
    
    # Get batch info
    batch_query = """
    SELECT id, batch_name, base_path, ref_version 
    FROM ecod_schema.batch 
    WHERE id = %s
    """
    
    batch_info = context.db.execute_dict_query(batch_query, (batch_id,))
    if not batch_info:
        logger.error(f"Batch {batch_id} not found")
        return 0
    
    base_path = batch_info[0]['base_path']
    ref_version = batch_info[0]['ref_version']
    batch_name = batch_info[0]['batch_name']
    
    logger.info(f"Processing batch {batch_id} ({batch_name})")
    
    # Find all HHSearch HHR files - use the specific pattern
    hhsearch_dir = os.path.join(base_path, "hhsearch")
    hhr_pattern = os.path.join(hhsearch_dir, f"*.hhsearch.{ref_version}.hhr")
    
    hhr_files = glob.glob(hhr_pattern)
    logger.info(f"Found {len(hhr_files)} HHR files on disk")
    
    if not hhr_files:
        logger.warning(f"No HHR files found in {hhsearch_dir}")
        return 0
    
    if limit:
        hhr_files = hhr_files[:limit]
        logger.info(f"Limited processing to {limit} files")
    
    # Create domains directory if it doesn't exist
    domains_dir = os.path.join(base_path, "domains")
    os.makedirs(domains_dir, exist_ok=True)
    
    # Process each HHR file
    processed_count = 0
    skipped_count = 0
    for hhr_file in hhr_files:
        # Extract PDB and chain ID from filename
        filename = os.path.basename(hhr_file)
        logger.debug(f"Processing file: {filename}")    
  
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain
        pdb_parts = pdb_chain.split('_')
        
        if len(pdb_parts) != 2:
            logger.warning(f"Invalid filename format: {filename}")
            continue
        
        pdb_id = pdb_parts[0]
        chain_id = pdb_parts[1]
        
        # Define output paths
        hh_xml_path = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml")
        chain_blast_path = os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.chainwise_blast.xml")
        domain_blast_path = os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.blast.xml")
        domain_summary_path = os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domains.xml")
        
        # Skip if domain summary already exists and not forcing
        if os.path.exists(domain_summary_path) and not force:
            logger.debug(f"Domain summary already exists for {pdb_chain}, skipping")
            skipped_count += 1
            continue
        
        try:
            # Parse HHR file
            logger.debug(f"Parsing HHR file: {hhr_file}")
            hhr_data = parse_hhr_file(hhr_file)
            
            if not hhr_data:
                logger.warning(f"Failed to parse HHR file: {hhr_file}")
                continue
            
            # Convert to XML
            logger.debug(f"Converting HHR data to XML for {pdb_chain}")
            xml_string = convert_hhr_to_xml(hhr_data, pdb_id, chain_id, ref_version)
            
            if not xml_string:
                logger.warning(f"Failed to convert HHR data to XML for {pdb_chain}")
                continue
            
            # Save HHSearch XML
            logger.debug(f"Saving HHSearch XML: {hh_xml_path}")
            try:
                with open(hh_xml_path, 'w', encoding='utf-8') as f:
                    f.write(xml_string)
            except Exception as e:
                logger.warning(f"Failed to save HHSearch XML: {hh_xml_path}, error: {str(e)}")
                continue
            
            # Collate evidence and create domain summary
            logger.debug(f"Collating evidence for {pdb_chain}")
            summary_xml = collate_evidence(
                chain_blast_path,
                domain_blast_path,
                hh_xml_path,
                pdb_id,
                chain_id,
                ref_version
            )
            
            if not summary_xml:
                logger.warning(f"Failed to collate evidence for {pdb_chain}")
                continue
            
            # Save domain summary
            logger.debug(f"Saving domain summary: {domain_summary_path}")
            try:
                with open(domain_summary_path, 'w', encoding='utf-8') as f:
                    f.write(summary_xml)
            except Exception as e:
                logger.warning(f"Failed to save domain summary: {domain_summary_path}, error: {str(e)}")
                continue
            
            # Lookup process_id for this chain
            chain_query = """
            SELECT ps.id as process_id
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE ps.batch_id = %s AND p.pdb_id = %s AND p.chain_id = %s
            """
            
            process_result = context.db.execute_dict_query(chain_query, (batch_id, pdb_id, chain_id))
            
            if process_result:
                process_id = process_result[0]['process_id']
                
                # Register HHSearch XML in database
                register_file(context.db, process_id, "hhsearch_xml", hh_xml_path)
                
                # Register domain summary in database
                register_file(context.db, process_id, "domain_summary", domain_summary_path)
                
                # Update process status
                update_process_status(context.db, process_id, "domain_summary_complete")
            
            processed_count += 1
            
            if processed_count % 10 == 0:
                logger.info(f"Processed {processed_count} chains so far")
        except Exception as e:
            logger.error(f"Error processing {pdb_chain}: {str(e)}", exc_info=True)
            error_count += 1
            continue 

    logger.info(f"Successfully processed {processed_count} chains in batch {batch_id}")
    return processed_count

def register_file(db, process_id, file_type, file_path):
    """Register file in database"""
    logger = logging.getLogger("file_registrar")
    
    try:
        # Check if file exists
        if not os.path.exists(file_path):
            logger.warning(f"File not found: {file_path}")
            return False
        
        file_size = os.path.getsize(file_path)
        
        # Check if already registered
        check_query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = %s
        """
        
        existing = db.execute_query(check_query, (process_id, file_type))
        
        if existing:
            # Update existing record
            db.update(
                "ecod_schema.process_file",
                {
                    "file_path": file_path,
                    "file_exists": True,
                    "file_size": file_size
                },
                "id = %s",
                (existing[0][0],)
            )
        else:
            # Insert new record
            db.insert(
                "ecod_schema.process_file",
                {
                    "process_id": process_id,
                    "file_type": file_type,
                    "file_path": file_path,
                    "file_exists": True,
                    "file_size": file_size
                }
            )
        
        return True
    except Exception as e:
        logger.error(f"Error registering file: {str(e)}")
        return False

def update_process_status(db, process_id, stage, error_message=None):
    """Update process status in database"""
    logger = logging.getLogger("status_updater")
    
    try:
        status = "error" if error_message else "success"
        
        db.update(
            "ecod_schema.process_status",
            {
                "current_stage": stage,
                "status": status,
                "error_message": error_message
            },
            "id = %s",
            (process_id,)
        )
        
        return True
    except Exception as e:
        logger.error(f"Error updating process status: {str(e)}")
        return False

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Process ECOD HHSearch Results')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, required=True,
                      help='Batch ID to process')
    parser.add_argument('--limit', type=int,
                      help='Limit the number of files to process')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--force', action='store_true',
                      help='Force reprocessing of already processed results')


    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    logger = logging.getLogger("main")
    logger.info(f"Starting HHSearch results processing for batch {args.batch_id}")
    
    processed_count = process_batch(args.batch_id, args.limit, args.force, config_path=args.config)
    
    if processed_count > 0:
        logger.info(f"Successfully processed HHSearch results for {processed_count} chains")
        return 0
    else:
        logger.warning("No chains were processed")
        return 1

if __name__ == "__main__":
    sys.exit(main())