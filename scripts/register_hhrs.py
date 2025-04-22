#!/usr/bin/env python3
"""
ECOD HHR Registration and XML Conversion Script

This script:
1. Finds HHR files that exist but aren't tracked in the database
2. Registers them in the database
3. Converts HHR files to XML format for domain analysis
4. Registers the XML files in the database

Usage:
    python register_hhrs.py --batch <batch_id> [--force] [--chains <pdb_id>_<chain_id> ...]
"""

import os
import sys
import logging
import argparse
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple, Set

# Import ECOD modules
from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.exceptions import PipelineError, ConfigurationError


def setup_logging(log_level=logging.INFO):
    """Configure logging"""
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=log_level, format=log_format)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Register and convert HHR files')
    
    parser.add_argument('--batch', type=int, required=True, 
                        help='Batch ID to process')
    
    parser.add_argument('--force', action='store_true', default=False,
                        help='Force regeneration of files even if they exist')
    
    parser.add_argument('--chains', nargs='+', default=None,
                        help='Specific chains to process (format: pdbid_chainid)')
    
    parser.add_argument('--config', type=str, default=None,
                        help='Path to configuration file')
    
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO', help='Logging level')
    
    return parser.parse_args()


class HHRParser:
    """Parse HHSearch result files (HHR format)"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def parse(self, hhr_file):
        """Parse HHR file and extract structured data
        
        Args:
            hhr_file: Path to HHR file
            
        Returns:
            Dictionary with parsed data
        """
        try:
            with open(hhr_file, 'r') as f:
                content = f.read()
                
            # Parse header information
            header = self._parse_header(content)
            
            # Parse hits
            hits = self._parse_hits(content)
            
            return {
                'header': header,
                'hits': hits
            }
        except Exception as e:
            self.logger.error(f"Error parsing HHR file {hhr_file}: {str(e)}")
            return None
    
    def _parse_header(self, content):
        """Parse header section of HHR file"""
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
            elif line.startswith('No_of_seqs '):
                parts = line.strip().split()
                if len(parts) > 1:
                    header['no_of_seqs'] = int(parts[1])
            
            # Stop at the beginning of the hit table
            if line.startswith(' No Hit'):
                break
                
        return header
    
    def _parse_hits(self, content):
        """Parse hits section of HHR file"""
        hits = []
        lines = content.split('\n')
        
        # Find the beginning of the hit table
        hit_table_start = None
        for i, line in enumerate(lines):
            if line.startswith(' No Hit'):
                hit_table_start = i + 1
                break
        
        if not hit_table_start:
            return hits
        
        # Process hits
        current_hit = None
        in_alignment = False
        query_ali = ""
        template_ali = ""
        
        i = hit_table_start
        while i < len(lines):
            line = lines[i].strip()
            
            # New hit begins with a line like " 1 e4tm9c1 etc"
            if line and re.match(r'^\s*\d+\s+\S+', line) and not in_alignment:
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
                    while j < len(lines) and not lines[j].startswith('>'):
                        if 'Probab=' in lines[j]:
                            prob_match = re.search(r'Probab=(\d+\.\d+)', lines[j])
                            if prob_match:
                                current_hit['probability'] = float(prob_match.group(1))
                        
                        if 'E-value=' in lines[j]:
                            eval_match = re.search(r'E-value=(\S+)', lines[j])
                            if eval_match:
                                try:
                                    current_hit['e_value'] = float(eval_match.group(1))
                                except ValueError:
                                    pass
                        
                        if 'Score=' in lines[j]:
                            score_match = re.search(r'Score=(\d+\.\d+)', lines[j])
                            if score_match:
                                current_hit['score'] = float(score_match.group(1))
                        
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
            
        return hits


class HHRToXMLConverter:
    """Convert HHR parsed data to XML format"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def convert(self, hhr_data, pdb_id, chain_id, ref_version):
        """Convert HHR data to XML
        
        Args:
            hhr_data: Parsed HHR data
            pdb_id: PDB ID
            chain_id: Chain ID
            ref_version: Reference version
            
        Returns:
            XML string
        """
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
                hit_elem.set("hit_id", hit.get('hit_id', ''))
                hit_elem.set("hh_prob", str(hit.get('probability', 0)))
                hit_elem.set("hh_evalue", str(hit.get('e_value', 0)))
                hit_elem.set("hh_score", str(hit.get('score', 0)))
                
                # Extract ECOD domain ID if it matches the pattern
                if re.match(r'[dge]\d\w{3}\w+\d+', hit.get('hit_id', '')):
                    hit_elem.set("ecod_domain_id", hit.get('hit_id', ''))
                
                # Extract and add query range
                if 'query_ali' in hit and 'query_start' in hit:
                    query_range = self._calculate_range(hit['query_ali'], hit['query_start'])
                    query_range_elem = ET.SubElement(hit_elem, "query_range")
                    query_range_elem.text = query_range
                
                # Add template range with coverage calculation
                if 'template_ali' in hit and 'template_start' in hit:
                    template_range = self._calculate_range(hit['template_ali'], hit['template_start'])
                    template_range_elem = ET.SubElement(hit_elem, "template_seqid_range")
                    template_range_elem.text = template_range
                    
                    # Calculate coverage
                    coverage = self._calculate_coverage(hit['query_ali'], hit['template_ali'])
                    template_range_elem.set("ungapped_coverage", str(coverage))
                    template_range_elem.set("coverage", str(coverage))
                
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
            self.logger.error(f"Error converting HHR data to XML: {str(e)}")
            return None
    
    def save(self, xml_string, output_path):
        """Save XML string to file
        
        Args:
            xml_string: XML string
            output_path: Output file path
            
        Returns:
            True if successful
        """
        try:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(xml_string)
            return True
        except Exception as e:
            self.logger.error(f"Error saving XML to {output_path}: {str(e)}")
            return False
    
    def _calculate_range(self, alignment, start_pos):
        """Calculate range from alignment
        
        Args:
            alignment: Alignment string (with gaps)
            start_pos: Starting position
            
        Returns:
            Range string in format "start-end,start-end,..."
        """
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
    
    def _calculate_coverage(self, query_ali, template_ali):
        """Calculate coverage between query and template alignments
        
        Args:
            query_ali: Query alignment string
            template_ali: Template alignment string
            
        Returns:
            Coverage as a float between 0 and 1
        """
        if not query_ali or not template_ali or len(query_ali) != len(template_ali):
            return 0.0
            
        # Count aligned (non-gap) positions
        aligned_positions = sum(1 for q, t in zip(query_ali, template_ali) if q != '-' and t != '-')
        total_template_positions = sum(1 for t in template_ali if t != '-')
        
        if total_template_positions == 0:
            return 0.0
            
        return aligned_positions / total_template_positions


class HHRRegistrar:
    """Class for finding and registering HHR files"""
    
    def __init__(self, context):
        self.context = context
        self.db = context.db
        self.logger = logging.getLogger("ecod.hhr_registrar")
        self.parser = HHRParser(self.logger)
        self.converter = HHRToXMLConverter(self.logger)
    
    def process_batch(self, batch_id, force=False, specific_chains=None):
        """Process a batch for HHR registration and conversion
        
        Args:
            batch_id: Batch ID to process
            force: Force regeneration of files even if they exist
            specific_chains: List of specific chain IDs to process
            
        Returns:
            Tuple of (registered_hhr_count, converted_xml_count)
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return (0, 0)
        
        # Find chains with HHR files
        chains = self._find_chains_with_hhr(batch_id, batch_info, specific_chains)
        if not chains:
            self.logger.warning(f"No chains with HHR files found in batch {batch_id}")
            return (0, 0)
        
        # Register and convert
        registered_count = 0
        converted_count = 0
        
        for chain in chains:
            hhr_registered = self._register_hhr(chain, batch_info, force)
            xml_converted = self._convert_hhr_to_xml(chain, batch_info, force)
            
            if hhr_registered:
                registered_count += 1
            
            if xml_converted:
                converted_count += 1
        
        self.logger.info(f"Processed {len(chains)} chains: registered {registered_count} HHR files, converted {converted_count} to XML")
        return (registered_count, converted_count)
    
    def _get_batch_info(self, batch_id):
        """Get batch information
        
        Args:
            batch_id: Batch ID
            
        Returns:
            Dictionary with batch information
        """
        query = """
        SELECT 
            id, 
            batch_name, 
            base_path, 
            ref_version
        FROM 
            ecod_schema.batch
        WHERE 
            id = %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            if rows:
                return rows[0]
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {e}")
            
        return None
    
    def _find_chains_with_hhr(self, batch_id, batch_info, specific_chains=None):
        """Find chains with HHR files
        
        Args:
            batch_id: Batch ID
            batch_info: Batch information dictionary
            specific_chains: List of specific chain IDs to process
            
        Returns:
            List of dictionaries with chain information
        """
        chains = []
        
        if specific_chains:
            self.logger.info(f"Looking for HHR files for {len(specific_chains)} specific chains")
            
            for chain_id in specific_chains:
                try:
                    pdb_id, chain_letter = chain_id.split('_')
                    
                    # Get process ID
                    query = """
                    SELECT 
                        ps.id as process_id,
                        p.id as protein_id,
                        p.pdb_id,
                        p.chain_id,
                        ps.relative_path
                    FROM 
                        ecod_schema.process_status ps
                    JOIN
                        ecod_schema.protein p ON ps.protein_id = p.id
                    WHERE 
                        ps.batch_id = %s
                        AND p.pdb_id = %s
                        AND p.chain_id = %s
                    LIMIT 1
                    """
                    
                    rows = self.db.execute_dict_query(query, (batch_id, pdb_id, chain_letter))
                    if rows:
                        chain = rows[0]
                        
                        # Find HHR file
                        hhr_path = self._find_hhr_file(chain, batch_info)
                        if hhr_path:
                            chain['hhr_path'] = hhr_path
                            chains.append(chain)
                except ValueError:
                    self.logger.warning(f"Invalid chain ID format: {chain_id}, expected pdbid_chainid")
            
            self.logger.info(f"Found HHR files for {len(chains)} out of {len(specific_chains)} specified chains")
        else:
            # Get all chains in the batch
            query = """
            SELECT 
                ps.id as process_id,
                p.id as protein_id,
                p.pdb_id,
                p.chain_id,
                ps.relative_path
            FROM 
                ecod_schema.process_status ps
            JOIN
                ecod_schema.protein p ON ps.protein_id = p.id
            WHERE 
                ps.batch_id = %s
                AND ps.status IN ('success', 'processing')
            """
            
            rows = self.db.execute_dict_query(query, (batch_id,))
            
            for chain in rows:
                hhr_path = self._find_hhr_file(chain, batch_info)
                if hhr_path:
                    chain['hhr_path'] = hhr_path
                    chains.append(chain)
            
            self.logger.info(f"Found HHR files for {len(chains)} out of {len(rows)} chains in batch {batch_id}")
        
        return chains
    
    def _find_hhr_file(self, chain, batch_info):
        """Find HHR file for a chain
        
        Args:
            chain: Chain information dictionary
            batch_info: Batch information dictionary
            
        Returns:
            Path to HHR file or None if not found
        """
        pdb_id = chain['pdb_id']
        chain_id = chain['chain_id']
        ref_version = batch_info['ref_version']
        base_path = batch_info['base_path']
        
        # Check standard locations
        potential_paths = [
            # New standard location (flat structure in hhsearch dir)
            os.path.join(base_path, "hhsearch", f"{pdb_id}_{chain_id}.{ref_version}.hhr"),
            
            # Old chain-specific directory structure
            os.path.join(base_path, chain['relative_path'], f"{pdb_id}_{chain_id}.{ref_version}.hhr"),
            
            # Other potential locations
            os.path.join(base_path, f"{pdb_id}_{chain_id}", f"{pdb_id}_{chain_id}.{ref_version}.hhr"),
            os.path.join(base_path, "ecod_dump", f"{pdb_id}_{chain_id}", f"{pdb_id}_{chain_id}.{ref_version}.hhr")
        ]
        
        for path in potential_paths:
            if os.path.exists(path) and os.path.getsize(path) > 0:
                return path
        
        return None
    
    def _register_hhr(self, chain, batch_info, force=False):
        """Register HHR file in database
        
        Args:
            chain: Chain information dictionary
            batch_info: Batch information dictionary
            force: Force registration even if already registered
            
        Returns:
            True if successful, False otherwise
        """
        process_id = chain['process_id']
        pdb_id = chain['pdb_id']
        chain_id = chain['chain_id']
        pdb_chain = f"{pdb_id}_{chain_id}"
        hhr_path = chain['hhr_path']
        
        # Check if already registered
        query = """
        SELECT id, file_path FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'hhr'
        """
        
        rows = self.db.execute_query(query, (process_id,))
        
        # Get relative path
        rel_path = os.path.relpath(hhr_path, batch_info['base_path'])
        
        if rows:
            # Already registered
            if not force:
                self.logger.info(f"HHR file already registered for {pdb_chain}, skipping")
                return True
            
            # Update registration
            try:
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(hhr_path)
                    },
                    "id = %s",
                    (rows[0][0],)
                )
                self.logger.info(f"Updated HHR registration for {pdb_chain}")
                return True
            except Exception as e:
                self.logger.error(f"Error updating HHR registration: {e}")
                return False
        else:
            # Register new file
            try:
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "hhr",
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(hhr_path)
                    }
                )
                
                # Update process status
                self.db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": "hhsearch_complete",
                        "status": "success"
                    },
                    "id = %s",
                    (process_id,)
                )
                
                self.logger.info(f"Registered HHR file for {pdb_chain}")
                return True
            except Exception as e:
                self.logger.error(f"Error registering HHR file: {e}")
                return False
    
    def _convert_hhr_to_xml(self, chain, batch_info, force=False):
        """Convert HHR file to XML and register
        
        Args:
            chain: Chain information dictionary
            batch_info: Batch information dictionary
            force: Force conversion even if XML already exists
            
        Returns:
            True if successful, False otherwise
        """
        process_id = chain['process_id']
        pdb_id = chain['pdb_id']
        chain_id = chain['chain_id']
        ref_version = batch_info['ref_version']
        base_path = batch_info['base_path']
        hhr_path = chain['hhr_path']
        
        # Define XML path
        hhsearch_dir = os.path.join(base_path, "hhsearch")
        os.makedirs(hhsearch_dir, exist_ok=True)
        xml_path = os.path.join(hhsearch_dir, f"{pdb_id}_{chain_id}.{ref_version}.hhsearch.xml")
        
        # Check if XML already exists
        if os.path.exists(xml_path) and os.path.getsize(xml_path) > 0 and not force:
            self.logger.info(f"XML file already exists for {pdb_id}_{chain_id}, checking registration")
            
            # Make sure it's registered
            return self._register_xml(process_id, xml_path, base_path)
        
        # Parse HHR file
        hhr_data = self.parser.parse(hhr_path)
        if not hhr_data:
            self.logger.error(f"Failed to parse HHR file: {hhr_path}")
            return False
        
        # Convert to XML
        xml_string = self.converter.convert(hhr_data, pdb_id, chain_id, ref_version)
        if not xml_string:
            self.logger.error(f"Failed to convert HHR to XML for {pdb_id}_{chain_id}")
            return False
        
        # Save XML
        if not self.converter.save(xml_string, xml_path):
            self.logger.error(f"Failed to save XML file: {xml_path}")
            return False
        
        # Register XML
        return self._register_xml(process_id, xml_path, base_path)
    
    def _register_xml(self, process_id, xml_path, base_path):
        """Register XML file in database
        
        Args:
            process_id: Process ID
            xml_path: Path to XML file
            base_path: Base path for relative path calculation
            
        Returns:
            True if successful, False otherwise
        """
        # Check if already registered
        query = """
        SELECT id FROM ecod_schema.process_file
        WHERE process_id = %s AND file_type = 'hhsearch_xml'
        """
        
        rows = self.db.execute_query(query, (process_id,))
        
        # Get relative path
        rel_path = os.path.relpath(xml_path, base_path)
        
        try:
            if rows:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(xml_path)
                    },
                    "id = %s",
                    (rows[0][0],)
                )
            else:
                # Create new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "hhsearch_xml",
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(xml_path)
                    }
                )
            
            return True
        except Exception as e:
            self.logger.error(f"Error registering XML file: {e}")
            return False


def main():
    """Main entry point"""
    args = parse_arguments()
    
    # Set up logging
    log_level = getattr(logging, args.log_level)
    setup_logging(log_level)
    
    logger = logging.getLogger("ecod.register_hhrs")
    logger.info("Starting HHR registration and conversion")
    
    # Initialize application context
    try:
        context = ApplicationContext(args.config)
    except ConfigurationError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    
    # Create registrar
    registrar = HHRRegistrar(context)
    
    try:
        # Process batch
        registered, converted = registrar.process_batch(
            args.batch, 
            args.force, 
            args.chains
        )
        
        logger.info(f"Successfully registered {registered} HHR files and converted {converted} to XML")
        
        if registered == 0 and converted == 0:
            return 1
        
        return 0
    except Exception as e:
        logger.error(f"Error processing batch: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())