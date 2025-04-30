#!/usr/bin/env python3
"""
BlastResultParser - A parser for BLAST XML result files

This module provides utilities for parsing and extracting information
from BLAST XML output files used in the pyECOD pipeline.
"""

import os
import re
import logging
import xml.etree.ElementTree as ET
from typing import Dict, List, Any, Optional, Tuple, Union

class BlastResultParser:
    """
    Parser for BLAST XML result files.
    
    This class handles parsing of BLAST XML output files,
    extracting hit information, and providing statistics
    about search results.
    """
    
    def __init__(self, logger=None):
        """Initialize the parser with optional logger
        
        Args:
            logger: Logger instance for output
        """
        self.logger = logger or logging.getLogger("blast_parser")
    
    def parse(self, file_path: str) -> Optional[Dict[str, Any]]:
        """Parse a BLAST XML result file
        
        Args:
            file_path: Path to BLAST XML file
            
        Returns:
            Dictionary with parsed data or None if parsing failed
        """
        try:
            # Check if file exists
            if not os.path.exists(file_path):
                self.logger.error(f"File not found: {file_path}")
                return None
                
            # Determine file type from path
            blast_type = self._determine_blast_type(file_path)
            
            # Parse the XML file
            tree = ET.parse(file_path)
            root = tree.getroot()
            
            # Extract basic info
            result = {
                'program': self._safe_find_text(root, ".//BlastOutput_program"),
                'db': self._safe_find_text(root, ".//BlastOutput_db"),
                'query_id': self._safe_find_text(root, ".//BlastOutput_query-ID"),
                'query_def': self._safe_find_text(root, ".//BlastOutput_query-def"),
                'query_len': self._safe_find_int(root, ".//BlastOutput_query-len"),
                'hits': [],
                'type': blast_type,
                'hit_count': 0,
                'has_hits': False
            }
            
            # Extract hits
            result['hits'] = self._extract_hits(root, blast_type)
            result['hit_count'] = len(result['hits'])
            result['has_hits'] = result['hit_count'] > 0
            
            # Store best hit if available
            if result['hits']:
                result['best_hit'] = result['hits'][0]
            
            return result
            
        except ET.ParseError as e:
            self.logger.error(f"XML parsing error for {file_path}: {str(e)}")
            return None
        except Exception as e:
            self.logger.error(f"Error parsing BLAST file {file_path}: {str(e)}")
            self.logger.exception("Stack trace:")
            return None
    
    def _determine_blast_type(self, file_path: str) -> str:
        """Determine BLAST type from file path
        
        Args:
            file_path: Path to BLAST file
            
        Returns:
            BLAST type ('chain_blast' or 'domain_blast')
        """
        path = file_path.lower().replace('\\', '/')
        
        if 'chain_blast' in path or 'chainwise' in path or '/blast/chain/' in path:
            return 'chain_blast'
        elif 'domain_blast' in path or '/blast/domain/' in path:
            return 'domain_blast'
        
        # Default to chain_blast if not determinable
        return 'chain_blast'
    
    def _safe_find_text(self, element: ET.Element, xpath: str, default: str = "") -> str:
        """Safely find text using XPath
        
        Args:
            element: XML Element
            xpath: XPath to search
            default: Default value if not found
            
        Returns:
            Text content or default if not found
        """
        result = element.find(xpath)
        if result is not None and result.text:
            return result.text.strip()
        return default
    
    def _safe_find_int(self, element: ET.Element, xpath: str, default: int = 0) -> int:
        """Safely find integer using XPath
        
        Args:
            element: XML Element
            xpath: XPath to search
            default: Default value if not found
            
        Returns:
            Integer value or default if not found
        """
        result = element.find(xpath)
        if result is not None and result.text:
            try:
                return int(result.text.strip())
            except ValueError:
                pass
        return default
    
    def _safe_find_float(self, element: ET.Element, xpath: str, default: float = 0.0) -> float:
        """Safely find float using XPath
        
        Args:
            element: XML Element
            xpath: XPath to search
            default: Default value if not found
            
        Returns:
            Float value or default if not found
        """
        result = element.find(xpath)
        if result is not None and result.text:
            try:
                return float(result.text.strip())
            except ValueError:
                pass
        return default
    
    def _extract_hits(self, root: ET.Element, blast_type: str) -> List[Dict[str, Any]]:
        """Extract hits from BLAST XML
        
        Args:
            root: Root XML element
            blast_type: BLAST type ('chain_blast' or 'domain_blast')
            
        Returns:
            List of hit dictionaries
        """
        hits = []
        
        # Find all Hit elements
        hit_elements = root.findall(".//Hit")
        
        for hit_idx, hit_elem in enumerate(hit_elements):
            hit_id = self._safe_find_text(hit_elem, "./Hit_id")
            hit_def = self._safe_find_text(hit_elem, "./Hit_def")
            hit_len = self._safe_find_int(hit_elem, "./Hit_len")
            
            # Parse hit_def to get PDB ID and chain ID
            pdb_id, chain_id, domain_id = self._parse_hit_def(hit_def, blast_type)
            
            # Extract HSPs
            hsps = []
            for hsp_elem in hit_elem.findall("./Hit_hsps/Hsp"):
                hsp = {
                    'hsp_num': self._safe_find_int(hsp_elem, "./Hsp_num"),
                    'bit_score': self._safe_find_float(hsp_elem, "./Hsp_bit-score"),
                    'score': self._safe_find_int(hsp_elem, "./Hsp_score"),
                    'evalue': self._safe_find_float(hsp_elem, "./Hsp_evalue", 999.0),
                    'query_from': self._safe_find_int(hsp_elem, "./Hsp_query-from"),
                    'query_to': self._safe_find_int(hsp_elem, "./Hsp_query-to"),
                    'hit_from': self._safe_find_int(hsp_elem, "./Hsp_hit-from"),
                    'hit_to': self._safe_find_int(hsp_elem, "./Hsp_hit-to"),
                    'identity': self._safe_find_int(hsp_elem, "./Hsp_identity"),
                    'positive': self._safe_find_int(hsp_elem, "./Hsp_positive"),
                    'align_len': self._safe_find_int(hsp_elem, "./Hsp_align-len"),
                    'query_seq': self._safe_find_text(hsp_elem, "./Hsp_qseq"),
                    'hit_seq': self._safe_find_text(hsp_elem, "./Hsp_hseq"),
                    'midline': self._safe_find_text(hsp_elem, "./Hsp_midline")
                }
                
                # Calculate identity percentage
                if hsp['align_len'] > 0:
                    hsp['identity_percent'] = (hsp['identity'] / hsp['align_len']) * 100
                else:
                    hsp['identity_percent'] = 0.0
                
                hsps.append(hsp)
            
            # Sort HSPs by evalue
            hsps.sort(key=lambda x: x['evalue'])
            
            # Create hit object
            hit = {
                'hit_num': hit_idx + 1,
                'hit_id': hit_id,
                'hit_def': hit_def,
                'hit_len': hit_len,
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'domain_id': domain_id,
                'hsps': hsps,
                'hsp_count': len(hsps),
                'has_multiple_hsps': len(hsps) > 1,
                'best_evalue': min([hsp['evalue'] for hsp in hsps]) if hsps else 999.0
            }
            
            # Get query coverage from best HSP
            if hsps and 'query_from' in hsps[0] and 'query_to' in hsps[0]:
                hit['query_range'] = f"{hsps[0]['query_from']}-{hsps[0]['query_to']}"
            
            # Get hit coverage from best HSP
            if hsps and 'hit_from' in hsps[0] and 'hit_to' in hsps[0]:
                hit['hit_range'] = f"{hsps[0]['hit_from']}-{hsps[0]['hit_to']}"
            
            hits.append(hit)
        
        # Sort hits by evalue
        hits.sort(key=lambda x: x['best_evalue'])
        
        return hits
    
    def _parse_hit_def(self, hit_def: str, blast_type: str) -> Tuple[str, str, str]:
        """Parse hit definition to extract PDB ID, chain ID, and domain ID
        
        Args:
            hit_def: Hit definition string
            blast_type: BLAST type ('chain_blast' or 'domain_blast')
            
        Returns:
            Tuple of (pdb_id, chain_id, domain_id)
        """
        pdb_id = ""
        chain_id = ""
        domain_id = ""
        
        # Common PDB chain pattern: pdb|1ABC|A or 1ABC_A
        pdb_chain_match = re.search(r'(?:pdb\|)?([a-zA-Z0-9]{4})(?:\||_)([a-zA-Z0-9])', hit_def)
        if pdb_chain_match:
            pdb_id = pdb_chain_match.group(1)
            chain_id = pdb_chain_match.group(2)
        
        # Domain ID pattern for domain BLAST (e.g., e4abcA1)
        if blast_type == 'domain_blast':
            domain_match = re.search(r'([dge]\d[a-zA-Z0-9]{3}[a-zA-Z0-9]\d+)', hit_def)
            if domain_match:
                domain_id = domain_match.group(1)
        
        return pdb_id, chain_id, domain_id
    
    def get_summary(self, file_path: str) -> Dict[str, Any]:
        """Get a summary of BLAST results
        
        Args:
            file_path: Path to BLAST XML file
            
        Returns:
            Dictionary with summary information
        """
        result = self.parse(file_path)
        if not result:
            return {
                'success': False,
                'error': 'Failed to parse file',
                'hit_count': 0,
                'has_hits': False
            }
        
        # Calculate some additional statistics
        hit_count = len(result['hits'])
        has_hits = hit_count > 0
        lowest_evalue = min([hit['best_evalue'] for hit in result['hits']]) if result['hits'] else None
        
        summary = {
            'success': True,
            'query_len': result['query_len'],
            'hit_count': hit_count,
            'has_hits': has_hits,
            'type': result['type'],
            'db': result['db'],
            'lowest_evalue': lowest_evalue
        }
        
        # Add information about best hit
        if has_hits:
            best_hit = result['hits'][0]
            summary['best_hit'] = {
                'hit_id': best_hit['hit_id'],
                'pdb_id': best_hit['pdb_id'],
                'chain_id': best_hit['chain_id'],
                'domain_id': best_hit['domain_id'],
                'evalue': best_hit['best_evalue'],
                'query_range': best_hit.get('query_range', '')
            }
        
        return summary
    
    def check_file_valid(self, file_path: str) -> bool:
        """Check if BLAST XML file is valid
        
        Args:
            file_path: Path to BLAST XML file
            
        Returns:
            Boolean indicating if file is valid
        """
        try:
            # Check if file exists
            if not os.path.exists(file_path):
                self.logger.error(f"File not found: {file_path}")
                return False
            
            # Check if file is empty
            if os.path.getsize(file_path) == 0:
                self.logger.error(f"File is empty: {file_path}")
                return False
            
            # Try parsing XML
            tree = ET.parse(file_path)
            root = tree.getroot()
            
            # Check for BlastOutput element
            if root.tag != 'BlastOutput':
                self.logger.error(f"Not a valid BLAST XML file: {file_path}")
                return False
            
            # Check for program element
            program = root.find(".//BlastOutput_program")
            if program is None or not program.text:
                self.logger.error(f"Missing program information: {file_path}")
                return False
            
            return True
            
        except ET.ParseError as e:
            self.logger.error(f"XML parsing error for {file_path}: {str(e)}")
            return False
        except Exception as e:
            self.logger.error(f"Error checking BLAST file {file_path}: {str(e)}")
            return False
    
    def extract_hits_as_list(self, file_path: str) -> List[Dict[str, Any]]:
        """Extract hits as a simple list for use in domain analysis
        
        Args:
            file_path: Path to BLAST XML file
            
        Returns:
            List of simplified hit dictionaries
        """
        result = self.parse(file_path)
        if not result or not result['hits']:
            return []
        
        simplified_hits = []
        
        for hit in result['hits']:
            # Only use best HSP for each hit
            best_hsp = hit['hsps'][0] if hit['hsps'] else {}
            
            simplified_hit = {
                'hit_id': hit['hit_id'],
                'hit_num': hit['hit_num'],
                'pdb_id': hit['pdb_id'],
                'chain_id': hit['chain_id'],
                'domain_id': hit['domain_id'],
                'evalue': best_hsp.get('evalue', 999.0),
                'bit_score': best_hsp.get('bit_score', 0.0),
                'identity_percent': best_hsp.get('identity_percent', 0.0),
                'query_range': hit.get('query_range', ''),
                'hit_range': hit.get('hit_range', '')
            }
            
            simplified_hits.append(simplified_hit)
        
        return simplified_hits
