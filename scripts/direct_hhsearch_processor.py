#!/usr/bin/env python3
"""
direct_hhsearch_processor.py - Process HHSearch results directly from filesystem
"""

import os
import sys
import argparse
import logging
import glob
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime

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

def parse_hhr_file(hhr_file):
    """Parse HHR file into structured data"""
    logger = logging.getLogger("hhr_parser")
    
    try:
        with open(hhr_file, 'r') as f:
            content = f.read()
        
        # Parse header
        header = {}
        lines = content.split('\n')
        
        for line in lines:
            if line.startswith('Query '):
                parts = line.strip().split()
                if len(parts) > 1:
                    header['query_id'] = parts[1]
            elif line.startswith('Match_columns '):
                parts = line.strip().split()
                if len(parts) > 1:
                    header['match_columns'] = int(parts[1])
            
            # Stop at the beginning of the hit table
            if line.startswith(' No Hit'):
                break
        
        # Parse hits
        hits = []
        hit_table_start = None
        
        for i, line in enumerate(lines):
            if line.startswith(' No Hit'):
                hit_table_start = i + 1
                break
        
        if not hit_table_start:
            return {'header': header, 'hits': []}
        
        # Process hits
        current_hit = None
        query_ali = ""
        template_ali = ""
        
        i = hit_table_start
        while i < len(lines):
            line = lines[i].strip()
            
            # New hit begins with a line like " 1 e4tm9c1 etc"
            if line and line[0].isdigit() and not line.startswith('Q ') and not line.startswith('T '):
                # Store previous hit if exists
                if current_hit and 'query_ali' in current_hit and 'template_ali' in current_hit:
                    hits.append(current_hit)
                
                # Parse hit line
                parts = line.split()
                if len(parts) >= 2:
                    hit_num = int(parts[0])
                    hit_id = parts[1]
                    
                    # Create new hit
                    current_hit = {
                        'hit_num': hit_num,
                        'hit_id': hit_id,
                        'probability': None,
                        'e_value': None,
                        'score': None,
                        'query_ali': "",
                        'template_ali': ""
                    }
                    
                    # Find probability, e-value, score
                    j = i + 1
                    while j < len(lines) and j < i + 5:
                        if 'Probab=' in lines[j]:
                            prob_parts = lines[j].split('Probab=')[1].split()[0]
                            try:
                                current_hit['probability'] = float(prob_parts)
                            except ValueError:
                                pass
                                
                        if 'E-value=' in lines[j]:
                            eval_parts = lines[j].split('E-value=')[1].split()[0]
                            try:
                                current_hit['e_value'] = float(eval_parts)
                            except ValueError:
                                pass
                                
                        if 'Score=' in lines[j]:
                            score_parts = lines[j].split('Score=')[1].split()[0]
                            try:
                                current_hit['score'] = float(score_parts)
                            except ValueError:
                                pass
                                
                        j += 1
            
            # Process alignment
            if line.startswith('Q '):
                parts = line.split()
                if len(parts) >= 4:
                    if 'query_start' not in current_hit:
                        try:
                            current_hit['query_start'] = int(parts[2])
                        except ValueError:
                            pass
                    current_hit['query_ali'] += parts[3]
                    
            elif line.startswith('T '):
                parts = line.split()
                if len(parts) >= 4:
                    if 'template_start' not in current_hit:
                        try:
                            current_hit['template_start'] = int(parts[2])
                        except ValueError:
                            pass
                    current_hit['template_ali'] += parts[3]
            
            i += 1
        
        # Add the last hit
        if current_hit and 'query_ali' in current_hit and 'template_ali' in current_hit:
            hits.append(current_hit)
            
        return {'header': header, 'hits': hits}
        
    except Exception as e:
        logger.error(f"Error parsing HHR file {hhr_file}: {str(e)}")
        return None

def convert_hhr_to_xml(hhr_data, pdb_id, chain_id, ref_version):
    """Convert HHR data to XML format"""
    logger = logging.getLogger("hhr_converter")
    
    try:
        # Create root element
        root = ET.Element("hh_summ_doc")
        
        # Add metadata
        metadata = ET.SubElement(root, "metadata")
        ET.SubElement(metadata, "pdb_id").text = pdb_id
        ET.SubElement(metadata, "chain_id").text = chain_id
        ET.SubElement(metadata, "reference").text = ref_version
        ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Add hits
        hits_elem = ET.SubElement(root, "hh_hit_list")
        
        for hit in hhr_data.get('hits', []):
            hit_elem = ET.SubElement(hits_elem, "hh_hit")
            hit_elem.set("hit_num", str(hit.get('hit_num', 0)))
            hit_elem.set("hit_id", str(hit.get('hit_id', '')))
            hit_elem.set("probability", str(hit.get('probability', 0)))
            hit_elem.set("e_value", str(hit.get('e_value', 0)))
            hit_elem.set("score", str(hit.get('score', 0)))
            
            # Calculate and add query range
            if 'query_ali' in hit and 'query_start' in hit:
                query_range = calculate_range(hit['query_ali'], hit['query_start'])
                query_range_elem = ET.SubElement(hit_elem, "query_range")
                query_range_elem.text = query_range
            
            # Add alignment details
            alignment_elem = ET.SubElement(hit_elem, "alignment")
            ET.SubElement(alignment_elem, "query_ali").text = hit.get('query_ali', '')
            ET.SubElement(alignment_elem, "template_ali").text = hit.get('template_ali', '')
        
        # Convert to string
        rough_string = ET.tostring(root, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")
        
        return pretty_xml
        
    except Exception as e:
        logger.error(f"Error converting HHR data to XML: {str(e)}")
        return None

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

def parse_xml_file(xml_file):
    """Parse XML file safely"""
    logger = logging.getLogger("xml_parser")
    
    try:
        if not os.path.exists(xml_file):
            logger.debug(f"XML file not found: {xml_file}")
            return None
            
        tree = ET.parse(xml_file)
        return tree
    except Exception as e:
        logger.warning(f"Error parsing XML file {xml_file}: {str(e)}")
        return None

def process_batch(batch_path, ref_version="develop291", limit=None, force=False):
    """Process HHSearch results for a batch"""
    logger = logging.getLogger("batch_processor")
    
    logger.info(f"Processing batch at {batch_path} with reference version {ref_version}")
    
    # Find all HHR files
    hhsearch_dir = os.path.join(batch_path, "hhsearch")
    hhr_pattern = os.path.join(hhsearch_dir, f"*.{ref_version}.hhr")
    
    hhr_files = glob.glob(hhr_pattern)
    logger.info(f"Found {len(hhr_files)} HHR files on disk")
    
    if not hhr_files:
        logger.warning(f"No HHR files found in {hhsearch_dir}")
        return 0
    
    if limit:
        hhr_files = hhr_files[:limit]
        logger.info(f"Limited processing to {limit} files")
    
    # Create domains directory if it doesn't exist
    domains_dir = os.path.join(batch_path, "domains")
    os.makedirs(domains_dir, exist_ok=True)
    
    # Process each HHR file
    processed_count = 0
    for hhr_file in hhr_files:
        # Extract PDB and chain ID from filename
        filename = os.path.basename(hhr_file)
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
            continue
        
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
        
        # Create simple domain summary
        root = ET.Element("domain_summ_doc")
        
        # Add metadata
        metadata = ET.SubElement(root, "metadata")
        ET.SubElement(metadata, "pdb_id").text = pdb_id
        ET.SubElement(metadata, "chain_id").text = chain_id
        ET.SubElement(metadata, "reference").text = ref_version
        ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Add evidence sections
        ET.SubElement(root, "chain_blast_evidence")
        ET.SubElement(root, "domain_blast_evidence")
        
        # Add HHSearch evidence
        hhsearch_elem = ET.SubElement(root, "hhsearch_evidence")
        hhsearch_data = parse_xml_file(hh_xml_path)
        
        if hhsearch_data is not None:
            hh_doc = hhsearch_data.find("hh_summ_doc")
            if hh_doc is not None:
                # Copy all elements from hh_doc to hhsearch_elem
                for child in hh_doc:
                    hhsearch_elem.append(ET.fromstring(ET.tostring(child)))
        
        # Add domain suggestions section (placeholder)
        ET.SubElement(root, "domain_suggestions")
        
        # Convert to string
        rough_string = ET.tostring(root, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")
        
        # Save domain summary
        logger.debug(f"Saving domain summary: {domain_summary_path}")
        try:
            with open(domain_summary_path, 'w', encoding='utf-8') as f:
                f.write(pretty_xml)
        except Exception as e:
            logger.warning(f"Failed to save domain summary: {domain_summary_path}, error: {str(e)}")
            continue
        
        processed_count += 1
        
        if processed_count % 10 == 0:
            logger.info(f"Processed {processed_count} chains so far")
    
    logger.info(f"Successfully processed {processed_count} chains")
    return processed_count

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description='Process ECOD HHSearch Results')
    parser.add_argument('--batch-path', type=str, required=True,
                      help='Path to batch directory')
    parser.add_argument('--ref-version', type=str, default="develop291",
                      help='Reference version')
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
    logger.info(f"Starting HHSearch results processing for batch at {args.batch_path}")
    
    processed_count = process_batch(
        args.batch_path, 
        args.ref_version, 
        args.limit, 
        args.force
    )
    
    if processed_count > 0:
        logger.info(f"Successfully processed HHSearch results for {processed_count} chains")
        return 0
    else:
        logger.warning("No chains were processed")
        return 1

if __name__ == "__main__":
    sys.exit(main())