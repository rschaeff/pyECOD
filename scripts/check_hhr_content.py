#!/usr/bin/env python3
"""
check_hhr_content.py - Check content of HHR and XML files
"""

import os
import glob
import argparse
import logging
import xml.etree.ElementTree as ET
import random

def setup_logging(verbose=False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def count_hits_in_hhr(hhr_file):
    """Count the number of hits in an HHR file"""
    try:
        with open(hhr_file, 'r') as f:
            content = f.read()
        
        lines = content.split('\n')
        hit_count = 0
        
        # Find the table header
        table_start = None
        for i, line in enumerate(lines):
            if line.startswith(' No Hit'):
                table_start = i + 1
                break
        
        if not table_start:
            return 0
        
        # Count hit entries
        for i in range(table_start, len(lines)):
            line = lines[i].strip()
            if line and line[0].isdigit() and not line.startswith('Q ') and not line.startswith('T '):
                hit_count += 1
        
        return hit_count
    except Exception as e:
        logging.warning(f"Error reading HHR file {hhr_file}: {str(e)}")
        return 0

def check_xml_content(xml_file):
    """Check the content of an XML file"""
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # Look for hh_hit elements in hh_hit_list
        hit_list = root.find(".//hh_hit_list")
        if hit_list is None:
            return 0, []
        
        hits = hit_list.findall("hh_hit")
        
        hit_info = []
        for hit in hits:
            hit_num = hit.get("hit_num", "unknown")
            hit_id = hit.get("hit_id", "unknown")
            hit_info.append(f"{hit_num}: {hit_id}")
        
        return len(hits), hit_info
    except Exception as e:
        logging.warning(f"Error parsing XML file {xml_file}: {str(e)}")
        return 0, []

def check_summary_content(summary_file):
    """Check the content of a domain summary file"""
    try:
        tree = ET.parse(summary_file)
        root = tree.getroot()
        
        # Look for hh_hit elements in hhsearch_evidence section
        hh_evidence = root.find(".//hhsearch_evidence")
        if hh_evidence is None:
            return 0
        
        hit_list = hh_evidence.find(".//hh_hit_list")
        if hit_list is None:
            return 0
        
        hits = hit_list.findall("hh_hit")
        return len(hits)
    except Exception as e:
        logging.warning(f"Error parsing summary file {summary_file}: {str(e)}")
        return 0

def examine_files(batch_path, ref_version="develop291", sample_size=5):
    """Examine HHR, XML, and summary files"""
    logger = logging.getLogger("file_examiner")
    
    # Find all HHR files
    hhsearch_dir = os.path.join(batch_path, "hhsearch")
    domains_dir = os.path.join(batch_path, "domains")
    
    hhr_pattern = os.path.join(hhsearch_dir, f"*.{ref_version}.hhr")
    hhr_files = glob.glob(hhr_pattern)
    
    if not hhr_files:
        logger.error(f"No HHR files found in {hhsearch_dir}")
        return
    
    logger.info(f"Found {len(hhr_files)} HHR files")
    
    # Select random samples
    if len(hhr_files) > sample_size:
        sample_files = random.sample(hhr_files, sample_size)
    else:
        sample_files = hhr_files
    
    logger.info(f"Examining {len(sample_files)} sample files")
    
    for hhr_file in sample_files:
        filename = os.path.basename(hhr_file)
        pdb_chain = filename.split('.')[0]  # Format: pdbid_chain
        
        xml_file = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml")
        summary_file = os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domains.xml")
        
        # Check HHR file
        hhr_hits = count_hits_in_hhr(hhr_file)
        logger.info(f"{pdb_chain} HHR file has {hhr_hits} hits")
        
        # Check XML file
        xml_hits, hit_info = check_xml_content(xml_file)
        logger.info(f"{pdb_chain} XML file has {xml_hits} hits")
        
        # Show top 5 hits from XML
        if hit_info:
            top_hits = hit_info[:5]
            logger.info(f"{pdb_chain} top 5 hits: {', '.join(top_hits)}")
        
        # Check summary file
        summary_hits = check_summary_content(summary_file)
        logger.info(f"{pdb_chain} summary file has {summary_hits} hits")
        
        # Check contents of XML and summary
        if os.path.exists(xml_file):
            with open(xml_file, 'r') as f:
                xml_content = f.read()
            logger.debug(f"{pdb_chain} XML size: {len(xml_content)} bytes")
            
            # Check if content looks valid
            if "hh_summ_doc" in xml_content and "hh_hit_list" in xml_content:
                logger.debug(f"{pdb_chain} XML has valid structure")
            else:
                logger.warning(f"{pdb_chain} XML is missing expected elements")
        
        if os.path.exists(summary_file):
            with open(summary_file, 'r') as f:
                summary_content = f.read()
            logger.debug(f"{pdb_chain} summary size: {len(summary_content)} bytes")
            
            # Check if content looks valid
            if "domain_summ_doc" in summary_content and "hhsearch_evidence" in summary_content:
                logger.debug(f"{pdb_chain} summary has valid structure")
            else:
                logger.warning(f"{pdb_chain} summary is missing expected elements")
        
        logger.info("------------------------------------")

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Check content of HHR and XML files')
    parser.add_argument('--batch-path', type=str, required=True,
                      help='Path to batch directory')
    parser.add_argument('--ref-version', type=str, default="develop291",
                      help='Reference version')
    parser.add_argument('--sample-size', type=int, default=5,
                      help='Number of files to examine')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose)
    
    logger = logging.getLogger("main")
    logger.info(f"Examining HHR and XML files for batch at {args.batch_path}")
    
    examine_files(args.batch_path, args.ref_version, args.sample_size)
    
    return 0

if __name__ == "__main__":
    main()