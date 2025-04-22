#!/usr/bin/env python3
"""
HHSearch utility module for parsing and converting HHR files.
Consolidated from multiple implementations to avoid code duplication.
"""

import os
import re
import logging
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import Dict, Any, List, Optional, Tuple


class HHRParser:
    """Parse HHSearch result files (HHR format)"""
    
    def __init__(self, logger=None):
        """Initialize parser with logger"""
        self.logger = logger or logging.getLogger("ecod.utils.hhr_parser")
    
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
            match = re.match(r'^\s*(\d+)\s+(\S+)', line)
            if match and not in_alignment:
                # Store previous hit if exists
                if current_hit and 'query_ali' in current_hit and 'template_ali' in current_hit:
                    hits.append(current_hit)
                
                # Parse hit line
                hit_num = int(match.group(1))
                hit_id = match.group(2)
                
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
    
    def __init__(self, logger=None):
        """Initialize converter with logger"""
        self.logger = logger or logging.getLogger("ecod.utils.hhr_converter")
    
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
                
                # Store probability under both attribute names for backward compatibility
                probability = hit.get('probability', 0)
                hit_elem.set("probability", str(probability))
                hit_elem.set("hh_prob", str(probability))
                
                # Store e-value under both attribute names for backward compatibility
                e_value = hit.get('e_value', 0)
                hit_elem.set("e_value", str(e_value))
                hit_elem.set("hh_evalue", str(e_value))
                
                # Store score under both attribute names for backward compatibility
                score = hit.get('score', 0)
                hit_elem.set("score", str(score))
                hit_elem.set("hh_score", str(score))
                
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
            # Ensure directory exists
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