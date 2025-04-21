#!/usr/bin/env python3
"""
check_hhr_results.py - Compare HHR files with XML summaries
"""

import os
import glob
import argparse
import logging
import xml.etree.ElementTree as ET

def setup_logging(verbose=False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def count_hits_in_xml(xml_file):
    """Count the number of hits in an XML file"""
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # Look for HH hits in hhsearch_evidence section
        hits = root.findall(".//hhsearch_evidence/hh_hit_list/hh_hit")
        return len(hits)
    except Exception as e:
        logging.warning(f"Error parsing XML file {xml_file}: {str(e)}")
        return 0

def check_batch(batch_path, ref_version="develop291", sample_size=10):
    """Check HHR results against XML summaries"""
    logger = logging.getLogger("batch_checker")
    
    # Find all HHR files
    hhsearch_dir = os.path.join(batch_path, "hhsearch")
    domains_dir = os.path.join(batch_path, "domains")
    
    hhr_pattern = os.path.join(hhsearch_dir, f"*.{ref_version}.hhr")
    xml_pattern = os.path.join(hhsearch_dir, f"*.{ref_version}.hhsearch.xml")
    summary_pattern = os.path.join(domains_dir, f"*.{ref_version}.domains.xml")
    
    hhr_files = glob.glob(hhr_pattern)
    xml_files = glob.glob(xml_pattern)
    summary_files = glob.glob(summary_pattern)
    
    logger.info(f"Found {len(hhr_files)} HHR files")
    logger.info(f"Found {len(xml_files)} HHSearch XML files")
    logger.info(f"Found {len(summary_files)} domain summary files")
    
    # Check if counts match
    if len(hhr_files) != len(xml_files):
        logger.warning(f"Number of HHR files ({len(hhr_files)}) doesn't match number of XML files ({len(xml_files)})")
    else:
        logger.info(f"Number of HHR files matches number of XML files: {len(hhr_files)}")
    
    if len(hhr_files) != len(summary_files):
        logger.warning(f"Number of HHR files ({len(hhr_files)}) doesn't match number of summary files ({len(summary_files)})")
    else:
        logger.info(f"Number of HHR files matches number of summary files: {len(summary_files)}")
    
    # Sample some files to check hit counts
    if sample_size > 0:
        # Limit to sample size
        summary_sample = summary_files[:sample_size] if sample_size < len(summary_files) else summary_files
        
        logger.info(f"Checking hit counts in {len(summary_sample)} sample summary files")
        
        hit_counts = {}
        for summary_file in summary_sample:
            filename = os.path.basename(summary_file)
            pdb_chain = filename.split('.')[0]  # Format: pdbid_chain
            
            # Count hits in summary XML
            hit_count = count_hits_in_xml(summary_file)
            hit_counts[pdb_chain] = hit_count
            
            logger.debug(f"{pdb_chain}: {hit_count} hits")
        
        # Calculate average hit count
        if hit_counts:
            avg_hits = sum(hit_counts.values()) / len(hit_counts)
            logger.info(f"Average number of hits per summary: {avg_hits:.2f}")
            
            # Find chains with most and least hits
            max_chain = max(hit_counts.items(), key=lambda x: x[1])
            min_chain = min(hit_counts.items(), key=lambda x: x[1])
            
            logger.info(f"Chain with most hits: {max_chain[0]} ({max_chain[1]} hits)")
            logger.info(f"Chain with least hits: {min_chain[0]} ({min_chain[1]} hits)")
    
    # Check for missing files
    missing_xml = []
    missing_summary = []
    
    for hhr_file in hhr_files:
        filename = os.path.basename(hhr_file)
        pdb_chain = filename.split('.')[0]  # Format: pdbid_chain
        
        expected_xml = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml")
        expected_summary = os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domains.xml")
        
        if not os.path.exists(expected_xml):
            missing_xml.append(pdb_chain)
        
        if not os.path.exists(expected_summary):
            missing_summary.append(pdb_chain)
    
    if missing_xml:
        logger.warning(f"{len(missing_xml)} chains are missing XML files")
        if len(missing_xml) <= 5:
            logger.warning(f"Missing XML for: {', '.join(missing_xml)}")
    else:
        logger.info("All HHR files have corresponding XML files")
    
    if missing_summary:
        logger.warning(f"{len(missing_summary)} chains are missing summary files")
        if len(missing_summary) <= 5:
            logger.warning(f"Missing summaries for: {', '.join(missing_summary)}")
    else:
        logger.info("All HHR files have corresponding summary files")

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Check HHR results against XML summaries')
    parser.add_argument('--batch-path', type=str, required=True,
                      help='Path to batch directory')
    parser.add_argument('--ref-version', type=str, default="develop291",
                      help='Reference version')
    parser.add_argument('--sample-size', type=int, default=10,
                      help='Number of files to sample for hit count check')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose)
    
    logger = logging.getLogger("main")
    logger.info(f"Checking HHR results for batch at {args.batch_path}")
    
    check_batch(args.batch_path, args.ref_version, args.sample_size)
    
    return 0

if __name__ == "__main__":
    main()