#!/usr/bin/env python3
"""
fix_domain_summaries.py - Fix domain summaries to include HHSearch hits
"""

import os
import glob
import argparse
import logging
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime

def setup_logging(verbose=False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def fix_domain_summary(batch_path, pdb_chain, ref_version):
    """Fix domain summary to include HHSearch hits"""
    logger = logging.getLogger("summary_fixer")
    
    hhsearch_dir = os.path.join(batch_path, "hhsearch")
    domains_dir = os.path.join(batch_path, "domains")
    
    # Define file paths
    xml_path = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml")
    summary_path = os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domains.xml")
    
    if not os.path.exists(xml_path):
        logger.warning(f"HHSearch XML file not found: {xml_path}")
        return False
    
    if not os.path.exists(summary_path):
        logger.warning(f"Domain summary file not found: {summary_path}")
        return False
    
    try:
        # Parse XML files
        hh_tree = ET.parse(xml_path)
        hh_root = hh_tree.getroot()
        
        summary_tree = ET.parse(summary_path)
        summary_root = summary_tree.getroot()
        
        # Extract pdb_id and chain_id from summary
        pdb_id = summary_root.find(".//pdb_id").text if summary_root.find(".//pdb_id") is not None else "unknown"
        chain_id = summary_root.find(".//chain_id").text if summary_root.find(".//chain_id") is not None else "unknown"
        
        # Find hhsearch_evidence section in summary
        hhsearch_elem = summary_root.find(".//hhsearch_evidence")
        if hhsearch_elem is None:
            logger.warning(f"No hhsearch_evidence section in summary for {pdb_chain}")
            return False
        
        # Clear existing content
        for child in list(hhsearch_elem):
            hhsearch_elem.remove(child)
        
        # Copy hit_list from HHSearch XML
        hit_list = hh_root.find(".//hh_hit_list")
        if hit_list is None:
            logger.warning(f"No hit_list found in HHSearch XML for {pdb_chain}")
            return False
        
        # Create a new hit_list in the hhsearch_evidence section
        new_hit_list = ET.SubElement(hhsearch_elem, "hh_hit_list")
        
        # Copy all hits
        hit_count = 0
        for hit in hit_list.findall("hh_hit"):
            # Create a copy of the hit element
            hit_string = ET.tostring(hit)
            new_hit = ET.fromstring(hit_string)
            new_hit_list.append(new_hit)
            hit_count += 1
        
        logger.info(f"Copied {hit_count} hits to domain summary for {pdb_chain}")
        
        # Save updated summary
        rough_string = ET.tostring(summary_root, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")
        
        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write(pretty_xml)
        
        logger.debug(f"Updated domain summary saved: {summary_path}")
        return True
        
    except Exception as e:
        logger.error(f"Error fixing domain summary for {pdb_chain}: {str(e)}")
        return False

def fix_batch_summaries(batch_path, ref_version="develop291", limit=None):
    """Fix all domain summaries in a batch"""
    logger = logging.getLogger("batch_fixer")
    
    # Find all HHSearch XML files
    hhsearch_dir = os.path.join(batch_path, "hhsearch")
    xml_pattern = os.path.join(hhsearch_dir, f"*.{ref_version}.hhsearch.xml")
    
    xml_files = glob.glob(xml_pattern)
    logger.info(f"Found {len(xml_files)} HHSearch XML files")
    
    if not xml_files:
        logger.warning(f"No HHSearch XML files found in {hhsearch_dir}")
        return 0
    
    # Limit processing if requested
    if limit and limit < len(xml_files):
        xml_files = xml_files[:limit]
        logger.info(f"Limited processing to {limit} files")
    
    # Process each XML file
    fixed_count = 0
    for xml_file in xml_files:
        # Extract pdb_chain from filename
        filename = os.path.basename(xml_file)
        parts = filename.split('.')
        pdb_chain = parts[0]  # Format: pdbid_chain
        
        # Fix domain summary
        success = fix_domain_summary(batch_path, pdb_chain, ref_version)
        
        if success:
            fixed_count += 1
            
        if fixed_count % 50 == 0 and fixed_count > 0:
            logger.info(f"Fixed {fixed_count} domain summaries so far")
    
    logger.info(f"Successfully fixed {fixed_count} domain summaries")
    return fixed_count

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Fix domain summaries to include HHSearch hits')
    parser.add_argument('--batch-path', type=str, required=True,
                      help='Path to batch directory')
    parser.add_argument('--ref-version', type=str, default="develop291",
                      help='Reference version')
    parser.add_argument('--limit', type=int,
                      help='Limit the number of files to process')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose)
    
    logger = logging.getLogger("main")
    logger.info(f"Fixing domain summaries for batch at {args.batch_path}")
    
    fixed_count = fix_batch_summaries(args.batch_path, args.ref_version, args.limit)
    
    if fixed_count > 0:
        logger.info(f"Successfully fixed {fixed_count} domain summaries")
        return 0
    else:
        logger.warning("No domain summaries were fixed")
        return 1

if __name__ == "__main__":
    main()