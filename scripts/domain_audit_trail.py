import os
import sys
import xml.etree.ElementTree as ET
import argparse
import logging
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import defaultdict
import psycopg2
from psycopg2.extras import RealDictCursor
import pandas as pd
import numpy as np
from Bio import SeqIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("DomainAuditTrail")

class DomainAuditTrail:
    """
    Generate audit trails for domain partitioning evidence in pyECOD.
    
    This class extracts evidence from BLAST summaries and visualizes the decision
    process for domain boundary determination.
    """
    
    def __init__(self, config_path, local_config_path=None):
        """
        Initialize the audit trail generator.
        
        Args:
            config_path: Path to the configuration file
            local_config_path: Path to local configuration file (optional)
        """
        self.config_path = config_path
        self.local_config_path = local_config_path
        self.config = self._load_config(config_path, local_config_path)
        self.conn = self._connect_db()
        self.schema = self.config.get('database', {}).get('schema', 'ecod_schema')
        
        # Load domain colors from config or use defaults
        self.domain_colors = self.config.get('visualization', {}).get('domain_colors', [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ])
    
    def _load_config(self, config_path, local_config_path=None):
        """
        Load configuration from a YAML file and merge with local config.
        
        This ensures secrets like database passwords are maintained from
        the local configuration.
        
        Args:
            config_path: Path to the main configuration file
            local_config_path: Path to local configuration file (optional)
        """
        import yaml
        import os
        
        logger.debug(f"Loading main configuration from {config_path}")
        # Load main configuration
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Check for local configuration
        if local_config_path is None:
            # Use default location if not specified
            local_config_path = os.path.join(os.path.dirname(config_path), 'config.local.yml')
        
        if os.path.exists(local_config_path):
            logger.debug(f"Merging with local configuration from {local_config_path}")
            with open(local_config_path, 'r') as f:
                local_config = yaml.safe_load(f)
                # Merge configurations with local taking precedence
                config = self._deep_merge(config, local_config)
        else:
            logger.debug(f"No local configuration found at {local_config_path}")
        
        return config
    
    def _deep_merge(self, source, destination):
        """
        Deep merge two dictionaries with destination taking precedence.
        """
        for key, value in source.items():
            if key in destination:
                if isinstance(value, dict) and isinstance(destination[key], dict):
                    destination[key] = self._deep_merge(value, destination[key])
            else:
                destination[key] = value
        return destination
    
    def _connect_db(self):
        """Connect to the PostgreSQL database."""
        db_config = self.config.get('database', {})
        conn = psycopg2.connect(
            host=db_config.get('host', 'dione'),
            port=db_config.get('port', 45000),
            database=db_config.get('database', 'ecod_protein'),
            user=db_config.get('user', 'ecod'),
            password=db_config.get('password', '')
        )
        
        # Set schema search path if provided
        schema = db_config.get('schema', 'ecod_schema')
        if schema:
            cursor = conn.cursor()
            try:
                # Set search path to include the schema
                cursor.execute(f"SET search_path TO {schema}, public")
                conn.commit()
                logger.debug(f"Search path set to {schema}, public")
                
                # Verify schema exists
                cursor.execute("SELECT current_schema()")
                current_schema = cursor.fetchone()[0]
                logger.debug(f"Current schema: {current_schema}")
                
                # Test table exists
                cursor.execute(f"SELECT EXISTS (SELECT 1 FROM information_schema.tables WHERE table_schema = '{schema}' AND table_name = 'protein')")
                table_exists = cursor.fetchone()[0]
                if not table_exists:
                    logger.warning(f"Table '{schema}.protein' does not exist!")
                    # Try to find the correct schema
                    cursor.execute("SELECT table_schema FROM information_schema.tables WHERE table_name = 'protein' LIMIT 1")
                    result = cursor.fetchone()
                    if result:
                        correct_schema = result[0]
                        logger.info(f"Found 'protein' table in schema '{correct_schema}' instead")
                        # Update schema for this session
                        cursor.execute(f"SET search_path TO {correct_schema}, public")
                        conn.commit()
                        # Store the correct schema for future queries
                        self.schema = correct_schema
                    else:
                        logger.error("Could not find 'protein' table in any schema!")
                else:
                    self.schema = schema
            finally:
                cursor.close()
        
        return conn

    def get_protein_info(self, pdb_id, chain_id):
        """Get protein information from the database."""
        cursor = self.conn.cursor(cursor_factory=RealDictCursor)
        try:
            query = f"""
                SELECT p.id, p.pdb_id, p.chain_id, p.source_id, p.length, 
                       ps.sequence, ps.md5_hash, b.ref_version
                FROM {self.schema}.protein p
                JOIN {self.schema}.protein_sequence ps ON p.id = ps.protein_id
                JOIN {self.schema}.process_status proc ON p.id = proc.protein_id
                JOIN {self.schema}.batch b ON proc.batch_id = b.id
                WHERE p.pdb_id = %s AND p.chain_id = %s
                LIMIT 1
            """
            cursor.execute(query, (pdb_id, chain_id))
            result = cursor.fetchone()
            if not result:
                logger.error(f"Protein {pdb_id}_{chain_id} not found in database")
                return None
            return result
        except Exception as e:
            logger.error(f"Database error in get_protein_info: {str(e)}")
            return None
        finally:
            cursor.close()

    def get_batch_info(self, batch_id):
        """Get batch information from the database."""
        cursor = self.conn.cursor(cursor_factory=RealDictCursor)
        try:
            query = f"""
                SELECT * FROM {self.schema}.batch WHERE id = %s
            """
            cursor.execute(query, (batch_id,))
            return cursor.fetchone()
        except Exception as e:
            logger.error(f"Database error in get_batch_info: {str(e)}")
            return None
        finally:
            cursor.close()

    def get_file_paths(self, protein_id, batch_id):
        """Get file paths for a protein from the database."""
        cursor = self.conn.cursor(cursor_factory=RealDictCursor)
        try:
            query = f"""
                SELECT pf.file_type, pf.file_path, pf.file_exists
                FROM {self.schema}.process_file pf
                JOIN {self.schema}.process_status ps ON pf.process_id = ps.id
                WHERE ps.protein_id = %s AND ps.batch_id = %s
            """
            cursor.execute(query, (protein_id, batch_id))
            files = {}
            for row in cursor.fetchall():
                files[row['file_type']] = {
                    'path': row['file_path'],
                    'exists': row['file_exists']
                }
            return files
        except Exception as e:
            logger.error(f"Database error in get_file_paths: {str(e)}")
            return {}
        finally:
            cursor.close()

    def parse_blast_summary(self, xml_path):
        """Parse BLAST summary XML file to extract hit information."""
        logger.info(f"Parsing BLAST summary from {xml_path}")
        if not os.path.exists(xml_path):
            logger.error(f"BLAST summary file not found: {xml_path}")
            return None
        
        try:
            # Parse the XML
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            # Extract hits from both chain-level and domain-level BLAST runs
            hits = []
            
            # For debugging
            logger.debug(f"Root tag: {root.tag}")
            
            # Look for chain level blast hits - in the format <chain_blast_run><hits><hit>
            chain_blast = root.find('.//chain_blast_run')
            if chain_blast is not None:
                hits_elem = chain_blast.find('hits')
                if hits_elem is not None:
                    logger.debug(f"Found chain_blast_run with {len(hits_elem)} hits")
                    for hit in hits_elem.findall('hit'):
                        # Extract basic hit info
                        hit_info = {
                            'type': 'chain',
                            'hit_id': hit.get('pdb_id', '') + '_' + hit.get('chain_id', ''),
                            'hit_desc': '',
                            'hit_len': int(hit.get('hsp_count', 0)),
                            'bit_score': 0.0,
                            'evalue': float(hit.get('evalues', 0.0)),
                            'query_coverage': 0,
                            'hsps': []
                        }
                        
                        # Extract query region directly from hit (not in HSPs)
                        query_reg = hit.find('query_reg')
                        if query_reg is not None and query_reg.text:
                            query_range = query_reg.text
                            if '-' in query_range:
                                query_from, query_to = map(int, query_range.split('-'))
                                # Create a single HSP for this hit
                                hsp_info = {
                                    'hsp_num': 1,
                                    'bit_score': hit_info['bit_score'],
                                    'evalue': hit_info['evalue'],
                                    'identity': 0.0,
                                    'align_len': 0,
                                    'query_from': query_from,
                                    'query_to': query_to,
                                    'hit_from': 0,
                                    'hit_to': 0
                                }
                                hit_info['hsps'].append(hsp_info)
                        
                        hits.append(hit_info)
            
            # Look for domain level blast hits - in the format <blast_run><hits><hit>
            domain_blast = root.find('.//blast_run')
            if domain_blast is not None:
                hits_elem = domain_blast.find('hits')
                if hits_elem is not None:
                    logger.debug(f"Found blast_run with {len(hits_elem)} hits")
                    for hit in hits_elem.findall('hit'):
                        # Extract basic hit info
                        hit_info = {
                            'type': 'domain',
                            'hit_id': hit.get('domain_id', ''),
                            'hit_desc': '',
                            'hit_len': int(hit.get('hsp_count', 0)),
                            'bit_score': 0.0,
                            'evalue': float(hit.get('evalues', 0.0)),
                            'query_coverage': 0,
                            'hsps': []
                        }
                        
                        # Extract query region directly from hit (not in HSPs)
                        query_reg = hit.find('query_reg')
                        if query_reg is not None and query_reg.text:
                            query_range = query_reg.text
                            if '-' in query_range:
                                query_from, query_to = map(int, query_range.split('-'))
                                # Create a single HSP for this hit
                                hsp_info = {
                                    'hsp_num': 1,
                                    'bit_score': hit_info['bit_score'],
                                    'evalue': hit_info['evalue'],
                                    'identity': 0.0,
                                    'align_len': 0,
                                    'query_from': query_from,
                                    'query_to': query_to,
                                    'hit_from': 0,
                                    'hit_to': 0
                                }
                                hit_info['hsps'].append(hsp_info)
                        
                        hits.append(hit_info)
            
            logger.debug(f"Extracted {len(hits)} hits with {sum(len(hit['hsps']) for hit in hits)} HSPs")
            
            # Debug output for the first few hits
            for i, hit in enumerate(hits[:5]):
                if hit['hsps']:
                    logger.debug(f"Hit {i+1}: {hit['hit_id']} ({len(hit['hsps'])} HSPs)")
                    for j, hsp in enumerate(hit['hsps']):
                        logger.debug(f"  HSP {j+1}: {hsp['query_from']}-{hsp['query_to']} (E-value: {hit['evalue']})")
            
            return hits
        except Exception as e:
            logger.error(f"Error parsing BLAST summary: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return None
    
    def _parse_blast_summary_fallback(self, xml_path):
        """Alternative parser for BLAST summary files with different structure."""
        logger.info("Using fallback parser for BLAST summary")
        try:
            # Try to extract HSPs directly from the file using regex
            import re
            hits = []
            
            with open(xml_path, 'r') as f:
                content = f.read()
            
            # Pattern for finding hits
            hit_pattern = r'<(hit|Hit)\s+([^>]*)>(.*?)</\1>'
            hit_matches = re.finditer(hit_pattern, content, re.DOTALL)
            
            for i, hit_match in enumerate(hit_matches):
                hit_attrs = hit_match.group(2)
                hit_content = hit_match.group(3)
                
                # Extract hit ID and other attributes
                hit_id_match = re.search(r'id="([^"]+)"', hit_attrs)
                hit_id = hit_id_match.group(1) if hit_id_match else f"hit_{i}"
                
                hit_info = {
                    'type': 'unknown',
                    'hit_id': hit_id,
                    'hit_desc': '',
                    'hit_len': 0,
                    'bit_score': 0.0,
                    'evalue': 0.0,
                    'query_coverage': 0,
                    'hsps': []
                }
                
                # Look for basic hit information
                for field in ['hit_desc', 'hit_len', 'bit_score', 'evalue']:
                    field_match = re.search(f'<{field}>(.*?)</{field}>', hit_content, re.DOTALL)
                    if field_match:
                        value = field_match.group(1).strip()
                        if field in ['hit_len']:
                            hit_info[field] = int(value) if value.isdigit() else 0
                        elif field in ['bit_score', 'evalue']:
                            try:
                                hit_info[field] = float(value)
                            except (ValueError, TypeError):
                                pass
                        else:
                            hit_info[field] = value
                
                # Extract HSPs - try to find any elements with query/hit ranges
                hsp_pattern = r'<(hsp|Hsp)\s*[^>]*>(.*?)</\1>'
                hsp_matches = re.finditer(hsp_pattern, hit_content, re.DOTALL)
                
                for hsp_match in hsp_matches:
                    hsp_content = hsp_match.group(2)
                    hsp_info = {
                        'hsp_num': len(hit_info['hsps']) + 1,
                        'bit_score': 0.0,
                        'evalue': 0.0,
                        'identity': 0.0,
                        'align_len': 0,
                        'query_from': 0,
                        'query_to': 0,
                        'hit_from': 0, 
                        'hit_to': 0
                    }
                    
                    # Extract basic HSP fields
                    for field in ['hsp_bit_score', 'hsp_evalue', 'hsp_identity', 'hsp_align_len']:
                        field_match = re.search(f'<{field}>(.*?)</{field}>', hsp_content, re.DOTALL)
                        if field_match:
                            value = field_match.group(1).strip()
                            if field in ['hsp_align_len']:
                                hsp_info[field[4:]] = int(value) if value.isdigit() else 0
                            elif field in ['hsp_bit_score', 'hsp_evalue', 'hsp_identity']:
                                try:
                                    hsp_info[field[4:]] = float(value)
                                except (ValueError, TypeError):
                                    pass
                    
                    # Look for query range - try multiple patterns
                    query_range_patterns = [
                        r'<query_reg\s+from="(\d+)"\s+to="(\d+)"',
                        r'<query_range>(\d+)-(\d+)</query_range>',
                        r'<hsp_query_from>(\d+)</hsp_query_from>.*?<hsp_query_to>(\d+)</hsp_query_to>'
                    ]
                    
                    for pattern in query_range_patterns:
                        query_match = re.search(pattern, hsp_content, re.DOTALL)
                        if query_match:
                            try:
                                hsp_info['query_from'] = int(query_match.group(1))
                                hsp_info['query_to'] = int(query_match.group(2))
                                break
                            except (ValueError, TypeError):
                                pass
                    
                    # Only add HSP if we have valid query range
                    if hsp_info['query_from'] > 0 and hsp_info['query_to'] > 0:
                        hit_info['hsps'].append(hsp_info)
                
                # Only add hit if it has HSPs
                if hit_info['hsps']:
                    hits.append(hit_info)
            
            logger.info(f"Fallback parser extracted {len(hits)} hits with a total of {sum(len(hit['hsps']) for hit in hits)} HSPs")
            return hits
        except Exception as e:
            logger.error(f"Fallback parser failed: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return []
    
    def _extract_hit_info(self, hit_element, hit_type):
        """Extract information from a hit element in the XML."""
        hit_info = {
            'type': hit_type,
            'hit_id': hit_element.get('id', ''),
            'hit_desc': hit_element.findtext('hit_desc', ''),
            'hit_len': int(hit_element.findtext('hit_len', 0)),
            'bit_score': float(hit_element.findtext('bit_score', 0)),
            'evalue': float(hit_element.findtext('evalue', 0)),
            'query_coverage': 0,
            'hsps': []
        }
        
        # Extract HSPs
        hsps_elem = hit_element.find('hsps')
        if hsps_elem is not None:
            for hsp in hsps_elem.findall('hsp'):
                hsp_info = {
                    'hsp_num': int(hsp.findtext('hsp_num', 0)),
                    'bit_score': float(hsp.findtext('hsp_bit_score', 0)),
                    'evalue': float(hsp.findtext('hsp_evalue', 0)),
                    'identity': float(hsp.findtext('hsp_identity', 0)),
                    'align_len': int(hsp.findtext('hsp_align_len', 0)),
                    'query_from': 0,
                    'query_to': 0,
                    'hit_from': 0, 
                    'hit_to': 0
                }
                
                # Extract query region - try different element names and formats
                query_from = None
                query_to = None
                
                # Try explicit query_reg element with from/to attributes
                query_reg = hsp.find('query_reg')
                if query_reg is not None:
                    query_from = query_reg.get('from')
                    query_to = query_reg.get('to')
                
                # Try hsp_query_from and hsp_query_to elements
                if query_from is None:
                    query_from = hsp.findtext('hsp_query_from')
                if query_to is None:
                    query_to = hsp.findtext('hsp_query_to')
                
                # Try query_range element as text (format: "start-end")
                if query_from is None or query_to is None:
                    query_range = hsp.findtext('query_range')
                    if query_range and '-' in query_range:
                        range_parts = query_range.split('-')
                        if len(range_parts) == 2:
                            query_from = range_parts[0]
                            query_to = range_parts[1]
                
                # Convert to integers
                try:
                    if query_from:
                        hsp_info['query_from'] = int(query_from)
                    if query_to:
                        hsp_info['query_to'] = int(query_to)
                except (ValueError, TypeError):
                    logger.warning(f"Could not parse query range: {query_from}-{query_to}")
                
                # Extract hit region
                hit_reg = hsp.find('hit_reg')
                if hit_reg is not None:
                    hsp_info['hit_from'] = int(hit_reg.get('from', 0))
                    hsp_info['hit_to'] = int(hit_reg.get('to', 0))
                else:
                    # Try hsp_hit_from and hsp_hit_to elements
                    hit_from = hsp.findtext('hsp_hit_from')
                    hit_to = hsp.findtext('hsp_hit_to')
                    if hit_from:
                        hsp_info['hit_from'] = int(hit_from)
                    if hit_to:
                        hsp_info['hit_to'] = int(hit_to)
                
                # Only add HSP if we have valid query range
                if hsp_info['query_from'] > 0 and hsp_info['query_to'] > 0:
                    hit_info['hsps'].append(hsp_info)
                else:
                    logger.warning(f"Skipping HSP without valid query range: {hsp_info}")
        
        # Calculate query coverage
        if hit_info['hsps']:
            # Create a coverage array
            max_pos = max(hsp['query_to'] for hsp in hit_info['hsps'])
            coverage = [0] * (max_pos + 1)
            
            # Mark covered positions
            for hsp in hit_info['hsps']:
                for i in range(hsp['query_from'], hsp['query_to'] + 1):
                    if i < len(coverage):
                        coverage[i] = 1
            
            # Calculate coverage percentage
            covered_positions = sum(coverage)
            total_positions = len(coverage) - 1  # Subtract 1 because we 0-indexed
            hit_info['query_coverage'] = (covered_positions / total_positions) * 100 if total_positions > 0 else 0
        
        return hit_info

    def parse_domain_partition(self, xml_path):
        """Parse domain partition XML file to extract defined domains."""
        logger.info(f"Parsing domain partition from {xml_path}")
        if not os.path.exists(xml_path):
            logger.error(f"Domain partition file not found: {xml_path}")
            return None
        
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            domains = []
            
            # Look for domains in domain_list/domain structure
            domain_list = root.find('.//domain_list')
            if domain_list is not None:
                for domain_elem in domain_list.findall('./domain'):
                    domain_info = {
                        'domain_id': domain_elem.get('domain_id', ''),
                        'range': domain_elem.get('range', ''),
                        'from': 0,
                        'to': 0,
                        'confidence': 1.0,  # Default confidence
                        'evidence': []
                    }
                    
                    # Parse range
                    range_parts = domain_info['range'].split('-')
                    if len(range_parts) == 2:
                        domain_info['from'] = int(range_parts[0])
                        domain_info['to'] = int(range_parts[1])
                    
                    # Extract evidence from match elements
                    evidence_elems = domain_elem.findall('./evidence/match')
                    for evidence_elem in evidence_elems:
                        evidence_info = {
                            'type': evidence_elem.get('type', ''),
                            'source': evidence_elem.get('domain_id', ''),
                            'score': 0,
                            'evalue': float(evidence_elem.get('evalue', 0)),
                            'range': ''
                        }
                        
                        # Extract query range if available
                        query_range_elem = evidence_elem.find('query_range')
                        if query_range_elem is not None and query_range_elem.text:
                            evidence_info['range'] = query_range_elem.text
                        
                        domain_info['evidence'].append(evidence_info)
                    
                    domains.append(domain_info)
            
            return domains
        except Exception as e:
            logger.error(f"Error parsing domain partition: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def visualize_domains(self, protein_info, blast_hits, domains, output_dir):
        """
        Generate visualization showing protein, BLAST hits, and domain partitions.
        
        Args:
            protein_info: Dictionary with protein information
            blast_hits: List of BLAST hit dictionaries
            domains: List of domain dictionaries
            output_dir: Directory to save output files
        """
        protein_length = protein_info['length']
        sequence = protein_info['sequence']
        pdb_chain = f"{protein_info['pdb_id']}_{protein_info['chain_id']}"
        
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Plot parameters
        plot_height = len(blast_hits) + len(domains) + 3  # +3 for protein, title, and spacing
        hit_height = 0.7
        domain_height = 0.8
        
        # Limit the number of blast hits to display (for better visualization)
        max_hits_to_display = 50
        if len(blast_hits) > max_hits_to_display:
            # Sort by evalue (best first) and take top hits
            sorted_hits = sorted(blast_hits, key=lambda x: x.get('evalue', float('inf')))
            displayed_hits = sorted_hits[:max_hits_to_display]
            logger.info(f"Limiting visualization to top {max_hits_to_display} hits (out of {len(blast_hits)})")
        else:
            displayed_hits = blast_hits
        
        # Adjust plot height based on displayed hits
        plot_height = len(displayed_hits) + len(domains) + 3
        
        # Draw protein
        plt.barh(plot_height - 1, protein_length, height=1.0, color='lightgray', alpha=0.5)
        plt.text(protein_length / 2, plot_height - 1, f"{pdb_chain} (1-{protein_length})", 
                 ha='center', va='center', fontweight='bold')
        
        # Draw BLAST hits
        for i, hit in enumerate(displayed_hits):
            y_pos = plot_height - 2 - i
            
            # Draw HSPs
            for j, hsp in enumerate(hit['hsps']):
                width = hsp['query_to'] - hsp['query_from'] + 1
                x_pos = hsp['query_from']
                
                # Use darker color for higher bit scores
                # Handle case where bit_score might be 0 or missing
                color_intensity = 0.7  # Default intensity if no bit_score
                if 'bit_score' in hsp and hsp['bit_score'] > 0:
                    # Find max bit score in this hit's HSPs
                    max_bit_score = max((h.get('bit_score', 0) for h in hit['hsps']), default=1.0)
                    if max_bit_score > 0:  # Avoid division by zero
                        color_intensity = min(1.0, hsp['bit_score'] / max_bit_score)
                
                # Choose a base color based on hit type
                base_color = '#1f77b4' if hit['type'] == 'chain' else '#ff7f0e'
                
                plt.barh(y_pos, width, height=hit_height, left=x_pos, 
                         color=base_color, alpha=0.5 + (color_intensity * 0.5))
                
                # Add HSP label if space permits
                if width > 30:
                    plt.text(x_pos + width/2, y_pos, f"{hsp['query_from']}-{hsp['query_to']}", 
                             ha='center', va='center', fontsize=8)
            
            # Add hit label
            hit_name = hit['hit_id'].split('|')[-1] if '|' in hit['hit_id'] else hit['hit_id']
            hit_evalue = f" (E: {hit.get('evalue'):.1e})" if hit.get('evalue', 0) > 0 else ""
            plt.text(0, y_pos, f"{hit_name}{hit_evalue}", 
                     ha='right', va='center', fontsize=9, fontweight='bold')
        
        # Draw domain partitions
        for i, domain in enumerate(domains):
            y_pos = len(domains) - i
            
            # Get domain range
            if 'from' in domain and 'to' in domain:
                domain_start, domain_end = domain['from'], domain['to']
            elif 'range' in domain and '-' in domain['range']:
                domain_start, domain_end = map(int, domain['range'].split('-'))
            else:
                logger.warning(f"Could not determine range for domain {domain.get('domain_id', 'unknown')}")
                continue
                
            width = domain_end - domain_start + 1
            
            plt.barh(y_pos, width, height=domain_height, left=domain_start, 
                     color=self.domain_colors[i % len(self.domain_colors)], alpha=0.7)
            
            # Add domain label
            domain_id = domain.get('domain_id', f"Domain {i+1}")
            plt.text(domain_start + width/2, y_pos, f"{domain_id} ({domain_start}-{domain_end})", 
                     ha='center', va='center', fontweight='bold')
        
        # Set plot parameters
        plt.xlabel('Residue Position')
        plt.ylabel('Domain / BLAST Hit')
        plt.title(f'Domain Partitioning Evidence for {pdb_chain}')
        plt.xlim(0, protein_length + 1)
        plt.ylim(0, plot_height)
        plt.tight_layout()
        
        # Save figure
        output_path = os.path.join(output_dir, f"{pdb_chain}_domain_audit.png")
        plt.savefig(output_path, dpi=300)
        plt.close()
        
        logger.info(f"Saved visualization to {output_path}")
        return output_path

    def generate_coverage_report(self, protein_info, blast_hits, domains, output_dir):
        """
        Generate a report on alignment coverage for domain hits.
        
        Args:
            protein_info: Dictionary with protein information
            blast_hits: List of BLAST hit dictionaries
            domains: List of domain dictionaries
            output_dir: Directory to save output files
        """
        pdb_chain = f"{protein_info['pdb_id']}_{protein_info['chain_id']}"
        protein_length = protein_info['length']
        
        # Initialize coverage arrays
        protein_coverage = [0] * (protein_length + 1)
        domain_coverage = {}
        
        # Map HSPs to residues
        for hit in blast_hits:
            for hsp in hit['hsps']:
                # Update protein coverage
                for i in range(hsp['query_from'], hsp['query_to'] + 1):
                    if i <= protein_length:
                        protein_coverage[i] = 1
        
        # Map HSPs to domains
        for hit in blast_hits:
            for hsp in hit['hsps']:
                # Track whether this HSP was added to any domain
                hsp_added = False
                
                # Check if HSP overlaps with any domain
                for domain_id, domain_data in domain_coverage.items():
                    domain_start, domain_end = map(int, domain_data['range'].split('-'))
                    
                    # Check if HSP overlaps with domain
                    if (hsp['query_from'] <= domain_end and 
                        hsp['query_to'] >= domain_start):
                        
                        # Calculate overlap
                        overlap_start = max(hsp['query_from'], domain_start)
                        overlap_end = min(hsp['query_to'], domain_end)
                        overlap_length = overlap_end - overlap_start + 1
                        domain_length = domain_data['length']
                        
                        # Add HSP to domain's supporting evidence
                        domain_data['hsps'].append({
                            'hit_id': hit['hit_id'],
                            'query_range': f"{hsp['query_from']}-{hsp['query_to']}",
                            'overlap_range': f"{overlap_start}-{overlap_end}",
                            'overlap_length': overlap_length,
                            'overlap_percent': (overlap_length / domain_length) * 100 if domain_length > 0 else 0,
                            'bit_score': hsp['bit_score'],
                            'evalue': hsp['evalue']
                        })
                        
                        hsp_added = True
                
                # Debug output for HSPs that don't match any domain
                if not hsp_added and logger.isEnabledFor(logging.DEBUG):
                    logger.debug(f"HSP {hsp['query_from']}-{hsp['query_to']} from hit {hit['hit_id']} doesn't match any domain")
        
        # Count covered residues per domain and remove duplicates
        for domain_id, domain_data in domain_coverage.items():
            # Get domain range
            domain_start, domain_end = map(int, domain_data['range'].split('-'))
            domain_range = set(range(domain_start, domain_end + 1))
            
            # Get covered positions from all HSPs
            covered_positions = set()
            for hsp in domain_data['hsps']:
                overlap_start, overlap_end = map(int, hsp['overlap_range'].split('-'))
                covered_positions.update(range(overlap_start, overlap_end + 1))
            
            # Keep only positions within domain range
            covered_positions = covered_positions.intersection(domain_range)
            
            # Update covered residues count
            domain_data['covered_residues'] = len(covered_positions)
            
            # Update coverage percentage
            domain_length = domain_data['length']
            if domain_length > 0:
                domain_data['coverage_percent'] = (len(covered_positions) / domain_length) * 100
            
            # Add debug info for low coverage domains
            if domain_data['coverage_percent'] < 10 and logger.isEnabledFor(logging.DEBUG):
                logger.debug(f"Domain {domain_id} ({domain_data['range']}) has low coverage: {domain_data['coverage_percent']:.2f}%")
                logger.debug(f"  Evidence count: {domain_data['evidence_count']}")
                logger.debug(f"  Supporting HSPs: {len(domain_data['hsps'])}")
                
                # Look at the evidence directly from the domain definition
                for domain in domains:
                    if domain['domain_id'] == domain_id:
                        logger.debug(f"  Domain definition evidence:")
                        for evidence in domain['evidence']:
                            logger.debug(f"    Type: {evidence['type']}, Source: {evidence['source']}, E-value: {evidence['evalue']}, Range: {evidence['range']}")
        
        # Generate report
        report = {
            'protein': {
                'pdb_chain': pdb_chain,
                'length': protein_length,
                'covered_residues': sum(protein_coverage),
                'coverage_percent': (sum(protein_coverage) / protein_length) * 100 if protein_length > 0 else 0
            },
            'domains': domain_coverage,
            'blast_hits': {
                'total': len(blast_hits),
                'chain_level': sum(1 for hit in blast_hits if hit['type'] == 'chain'),
                'domain_level': sum(1 for hit in blast_hits if hit['type'] == 'domain')
            },
            'short_alignments': [
                {
                    'hit_id': hit['hit_id'],
                    'hsp_query_range': f"{hsp['query_from']}-{hsp['query_to']}",
                    'length': hsp['query_to'] - hsp['query_from'] + 1,
                    'bit_score': hsp['bit_score'],
                    'evalue': hsp['evalue']
                }
                for hit in blast_hits
                for hsp in hit['hsps']
                if hsp['query_to'] - hsp['query_from'] + 1 < 30  # Short alignment threshold
            ]
        }
        
        # Save report
        output_path = os.path.join(output_dir, f"{pdb_chain}_coverage_report.json")
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Saved coverage report to {output_path}")
        
        # Also generate a text report
        text_report_path = os.path.join(output_dir, f"{pdb_chain}_coverage_report.txt")
        with open(text_report_path, 'w') as f:
            f.write(f"Domain Partition Audit Trail for {pdb_chain}\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Protein Length: {protein_length} residues\n")
            f.write(f"Overall Coverage: {report['protein']['coverage_percent']:.2f}% ({report['protein']['covered_residues']} residues)\n")
            f.write(f"Total BLAST Hits: {report['blast_hits']['total']} (Chain: {report['blast_hits']['chain_level']}, Domain: {report['blast_hits']['domain_level']})\n\n")
            
            f.write("Domain Summary:\n")
            f.write("-" * 80 + "\n")
            for domain_id, domain_data in report['domains'].items():
                f.write(f"Domain {domain_id} ({domain_data['range']}):\n")
                f.write(f"  Length: {domain_data['length']} residues\n")
                f.write(f"  Coverage: {domain_data['coverage_percent']:.2f}% ({domain_data['covered_residues']} residues)\n")
                f.write(f"  Evidence Count: {domain_data['evidence_count']}\n")
                f.write(f"  Supporting HSPs: {len(domain_data['hsps'])}\n")
                
                if domain_data['hsps']:
                    f.write("\n  Top 5 HSPs by bit score:\n")
                    sorted_hsps = sorted(domain_data['hsps'], key=lambda x: x['bit_score'], reverse=True)
                    
                    for i, hsp in enumerate(sorted_hsps[:5]):
                        f.write(f"    {i+1}. {hsp['hit_id'].split('|')[-1] if '|' in hsp['hit_id'] else hsp['hit_id']}\n")
                        f.write(f"       Query: {hsp['query_range']}, Overlap: {hsp['overlap_range']} ({hsp['overlap_percent']:.2f}%)\n")
                        f.write(f"       Bit Score: {hsp['bit_score']:.2f}, E-value: {hsp['evalue']}\n")
                
                f.write("\n")
            
            if report['short_alignments']:
                f.write("\nShort Alignments (< 30 residues):\n")
                f.write("-" * 80 + "\n")
                
                for i, aln in enumerate(sorted(report['short_alignments'], key=lambda x: x['length'])):
                    f.write(f"{i+1}. {aln['hit_id'].split('|')[-1] if '|' in aln['hit_id'] else aln['hit_id']}\n")
                    f.write(f"   Range: {aln['hsp_query_range']} (Length: {aln['length']})\n")
                    f.write(f"   Bit Score: {aln['bit_score']:.2f}, E-value: {aln['evalue']}\n")
        
        logger.info(f"Saved text report to {text_report_path}")
        return output_path, text_report_path
    
    def process_protein(self, pdb_id, chain_id, batch_id=None, output_dir=None, import_to_db=True):
        """
        Process a single protein and generate audit trail.
        
        Args:
            pdb_id: PDB ID of the protein
            chain_id: Chain ID of the protein
            batch_id: Batch ID (optional, will be determined from database if not provided)
            output_dir: Directory to save output files (optional)
            import_to_db: Whether to import results to database (default: True)
            
        Returns:
            Dictionary with paths to generated files, or False if processing failed
        """
        # Get protein information
        protein_info = self.get_protein_info(pdb_id, chain_id)
        if not protein_info:
            return False
        
        # Determine batch if not provided
        if batch_id is None:
            cursor = self.conn.cursor(cursor_factory=RealDictCursor)
            try:
                query = f"""
                    SELECT batch_id FROM {self.schema}.process_status 
                    WHERE protein_id = %s 
                    ORDER BY id DESC LIMIT 1
                """
                cursor.execute(query, (protein_info['id'],))
                result = cursor.fetchone()
                if result:
                    batch_id = result['batch_id']
                else:
                    logger.error(f"No batch found for protein {pdb_id}_{chain_id}")
                    return False
            except Exception as e:
                logger.error(f"Database error retrieving batch ID: {str(e)}")
                return False
            finally:
                cursor.close()
        
        # Get batch information
        batch_info = self.get_batch_info(batch_id)
        if not batch_info:
            logger.error(f"Batch {batch_id} not found")
            return False
        
        # Log batch information for debugging
        logger.debug(f"Batch information: ID={batch_id}, Name={batch_info.get('batch_name')}, Path={batch_info.get('base_path')}")
        
        # Determine output directory
        if output_dir is None:
            output_dir = os.path.join(batch_info['base_path'], "audit")
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Created output directory: {output_dir}")
        
        # Get file paths
        file_paths = self.get_file_paths(protein_info['id'], batch_id)
        logger.debug(f"Retrieved {len(file_paths)} file paths from database")
        for file_type, file_info in file_paths.items():
            logger.debug(f"  - {file_type}: {file_info['path']} (Exists: {file_info['exists']})")
        
        # Find BLAST summary file
        blast_summary_path = None
        for file_type in ['blast_summ', 'blast_summary', 'domain_summary']:
            if file_type in file_paths and file_paths[file_type]['exists']:
                file_path = file_paths[file_type]['path']
                # Check if path is relative or absolute
                if not os.path.isabs(file_path):
                    # Combine with batch path
                    file_path = os.path.join(batch_info['base_path'], file_path)
                blast_summary_path = file_path
                logger.debug(f"Using {file_type} path: {blast_summary_path}")
                break
        
        if not blast_summary_path:
            # Try to construct the path
            batch_path = batch_info['base_path']
            ref_version = batch_info['ref_version']
            blast_summary_path = os.path.join(
                batch_path, "domains", 
                f"{pdb_id}_{chain_id}.{ref_version}.domains.xml"
            )
            
            if not os.path.exists(blast_summary_path):
                # Try alternate path patterns
                alt_paths = [
                    os.path.join(batch_path, "domains", f"{pdb_id}_{chain_id}.domains.xml"),
                    os.path.join(batch_path, "domains", f"{pdb_id}_{chain_id}.{ref_version}.domains_v14.xml"),
                    os.path.join(batch_path, "domains", f"{pdb_id}_{chain_id}.{ref_version}.blast_summ.xml"),
                    os.path.join(batch_path, "blast", f"{pdb_id}_{chain_id}.{ref_version}.blast_summ.xml"),
                    # Additional patterns for blast_only files
                    os.path.join(batch_path, "domains", f"{pdb_id}_{chain_id}.{ref_version}.blast_summ.blast_only.xml"),
                    os.path.join(batch_path, "blast", f"{pdb_id}_{chain_id}.{ref_version}.blast_summ.blast_only.xml"),
                    # Try to find any file with the PDB ID and chain ID in the domains directory
                    os.path.join(batch_path, "domains", f"{pdb_id}_{chain_id}*.xml")
                ]
                
                # Check for explicit paths
                for path in alt_paths:
                    if '*' in path:
                        # Handle wildcard paths with glob
                        import glob
                        matching_files = glob.glob(path)
                        if matching_files:
                            # Use the first matching file
                            blast_summary_path = matching_files[0]
                            logger.info(f"Found BLAST summary file using wildcard: {blast_summary_path}")
                            break
                    elif os.path.exists(path):
                        blast_summary_path = path
                        logger.info(f"Found BLAST summary file at: {blast_summary_path}")
                        break
                
                # If still not found, list directory contents for debugging
                if not os.path.exists(blast_summary_path):
                    domains_dir = os.path.join(batch_path, "domains")
                    if os.path.exists(domains_dir):
                        import glob
                        files = glob.glob(os.path.join(domains_dir, f"{pdb_id}_{chain_id}*"))
                        if files:
                            logger.info(f"Found {len(files)} files for {pdb_id}_{chain_id} in domains directory:")
                            for file in files:
                                logger.info(f"  - {file}")
                            # Use the first file found
                            blast_summary_path = files[0]
                            logger.info(f"Using file: {blast_summary_path}")
                        else:
                            logger.warning(f"No files found for {pdb_id}_{chain_id} in domains directory")
                    else:
                        logger.warning(f"Domains directory not found: {domains_dir}")
        
        if not blast_summary_path or not os.path.exists(blast_summary_path):
            logger.error(f"BLAST summary file not found for {pdb_id}_{chain_id}")
            return False
        
        # Find domain partition file
        domain_file_path = None
        for file_type in ['domain_file', 'partition', 'domain_partition']:
            if file_type in file_paths and file_paths[file_type]['exists']:
                file_path = file_paths[file_type]['path']
                # Check if path is relative or absolute
                if not os.path.isabs(file_path):
                    # Combine with batch path
                    file_path = os.path.join(batch_info['base_path'], file_path)
                domain_file_path = file_path
                logger.debug(f"Using {file_type} path: {domain_file_path}")
                break
        
        if not domain_file_path:
            # Try to construct the path for domain partition file
            batch_path = batch_info['base_path']
            ref_version = batch_info['ref_version']
            
            # Define the specific domain partition file pattern
            domain_pattern = f"{pdb_id}_{chain_id}.{ref_version}.domains_v14.xml"
            domain_file_path = os.path.join(batch_path, "domains", domain_pattern)
            
            if not os.path.exists(domain_file_path):
                logger.warning(f"Domain partition file not found at: {domain_file_path}")
                logger.info(f"Will try to infer domains from BLAST hits for {pdb_id}_{chain_id}")
        
        # Parse BLAST summary
        if not blast_summary_path or not os.path.exists(blast_summary_path):
            logger.error(f"BLAST summary file not found for {pdb_id}_{chain_id}: {blast_summary_path}")
            return False
            
        logger.info(f"Parsing BLAST summary from {blast_summary_path}")
        blast_hits = self.parse_blast_summary(blast_summary_path)
        if not blast_hits:
            logger.error(f"Failed to parse BLAST summary for {pdb_id}_{chain_id}")
            return False
        
        # Parse domain partition
        domains = None
        if domain_file_path and os.path.exists(domain_file_path):
            logger.info(f"Parsing domain partition from {domain_file_path}")
            domains = self.parse_domain_partition(domain_file_path)
        
        
        # If domain partition doesn't exist or failed to parse, create a placeholder
        if not domains:
            logger.warning(f"Domain partition file not found or failed to parse for {pdb_id}_{chain_id}")
            
            # Try to infer domains from BLAST hits
            domains = self._infer_domains_from_blast(blast_hits, protein_info['length'])
            
            if not domains:
                logger.error(f"Could not infer domains for {pdb_id}_{chain_id}")
                return False
        
        # Visualize domains
        viz_path = self.visualize_domains(protein_info, blast_hits, domains, output_dir)
        
        # Generate coverage report
        report_path, text_report_path = self.generate_coverage_report(
            protein_info, blast_hits, domains, output_dir
        )
        
        result = {
            'visualization': viz_path,
            'json_report': report_path,
            'text_report': text_report_path
        }
        
        # Import results to database if requested
        if import_to_db:
            try:
                self._import_results_to_db(
                    protein_info['id'], 
                    batch_id,
                    viz_path,
                    report_path,
                    text_report_path
                )
                logger.info(f"Imported audit results to database for {pdb_id}_{chain_id}")
            except Exception as e:
                logger.error(f"Failed to import results to database: {str(e)}")
        
        logger.info(f"Processing complete for {pdb_id}_{chain_id}")
        return result
        
    def _import_results_to_db(self, protein_id, batch_id, viz_path, json_report_path, text_report_path):
        """Import audit results to database."""
        # Load the JSON report
        with open(json_report_path, 'r') as f:
            report_data = json.load(f)
        
        # Convert to JSONB format for PostgreSQL
        report_jsonb = json.dumps(report_data)
        
        # Check if the domain_audit table exists
        cursor = self.conn.cursor()
        try:
            try:
                cursor.execute(f"SELECT 1 FROM {self.schema}.domain_audit LIMIT 1")
                table_exists = True
            except psycopg2.errors.UndefinedTable:
                logger.warning(f"Table {self.schema}.domain_audit doesn't exist. Skipping database import.")
                table_exists = False
            
            if not table_exists:
                return None
            
            # Call the import function
            query = f"""
                SELECT {self.schema}.import_domain_audit_report(
                    %s, %s, %s, %s, %s, %s::jsonb
                )
            """
            cursor.execute(query, (
                protein_id,
                batch_id,
                json_report_path,
                text_report_path,
                viz_path,
                report_jsonb
            ))
            
            audit_id = cursor.fetchone()[0]
            self.conn.commit()
            return audit_id
        except Exception as e:
            self.conn.rollback()
            if "function" in str(e) and "does not exist" in str(e):
                logger.warning(f"Function {self.schema}.import_domain_audit_report doesn't exist. Skipping database import.")
                return None
            else:
                logger.error(f"Database error during import: {str(e)}")
                raise e
        finally:
            cursor.close()
    
    def _infer_domains_from_blast(self, blast_hits, protein_length):
        """Infer domains from BLAST hits when domain partition is not available."""
        if not blast_hits:
            return None
        
        # Collect HSP ranges
        ranges = []
        for hit in blast_hits:
            for hsp in hit['hsps']:
                ranges.append((hsp['query_from'], hsp['query_to']))
        
        if not ranges:
            return None
        
        # Sort ranges by start position
        ranges.sort()
        
        # Merge overlapping ranges
        merged_ranges = []
        current_start, current_end = ranges[0]
        
        for start, end in ranges[1:]:
            if start <= current_end + 30:  # Allow small gaps (up to 30 residues)
                current_end = max(current_end, end)
            else:
                merged_ranges.append((current_start, current_end))
                current_start, current_end = start, end
        
        merged_ranges.append((current_start, current_end))
        
        # Create domain objects
        domains = []
        for i, (start, end) in enumerate(merged_ranges):
            domains.append({
                'domain_id': str(i+1),
                'range': f"{start}-{end}",
                'from': start,
                'to': end,
                'confidence': 0.5,  # Placeholder confidence
                'evidence': []  # No explicit evidence
            })
        
        # If no domains were found, create a single domain spanning the entire protein
        if not domains:
            domains.append({
                'domain_id': '1',
                'range': f"1-{protein_length}",
                'from': 1,
                'to': protein_length,
                'confidence': 0.5,
                'evidence': []
            })
        
        return domains

    def process_batch(self, batch_id, limit=None, output_dir=None, import_to_db=True):
        """
        Process all proteins in a batch.
        
        Args:
            batch_id: Batch ID to process
            limit: Limit the number of proteins to process (optional)
            output_dir: Directory to save output files (optional)
            import_to_db: Whether to import results to database (default: True)
            
        Returns:
            List of processing results, or False if processing failed
        """
        # Get batch information
        batch_info = self.get_batch_info(batch_id)
        if not batch_info:
            logger.error(f"Batch {batch_id} not found")
            return False
        
        # Determine output directory
        if output_dir is None:
            output_dir = os.path.join(batch_info['base_path'], "audit")
            os.makedirs(output_dir, exist_ok=True)
        
        # Get proteins in batch
        cursor = self.conn.cursor(cursor_factory=RealDictCursor)
        try:
            # Check if we should skip already audited proteins
            if import_to_db:
                query = f"""
                    SELECT p.id, p.pdb_id, p.chain_id
                    FROM {self.schema}.protein p
                    JOIN {self.schema}.process_status ps ON p.id = ps.protein_id
                    LEFT JOIN {self.schema}.domain_audit da ON p.id = da.protein_id AND ps.batch_id = da.batch_id
                    WHERE ps.batch_id = %s AND da.id IS NULL  -- Only proteins not yet audited
                    ORDER BY p.pdb_id, p.chain_id
                """
            else:
                query = f"""
                    SELECT p.id, p.pdb_id, p.chain_id
                    FROM {self.schema}.protein p
                    JOIN {self.schema}.process_status ps ON p.id = ps.protein_id
                    WHERE ps.batch_id = %s
                    ORDER BY p.pdb_id, p.chain_id
                """
            
            if limit:
                query += f" LIMIT {limit}"
            
            cursor.execute(query, (batch_id,))
            proteins = cursor.fetchall()
            
            if not proteins:
                logger.info(f"No proteins found in batch {batch_id} that need processing")
                return False
            
            logger.info(f"Processing {len(proteins)} proteins in batch {batch_id}")
            
            # Process each protein
            results = []
            for protein in proteins:
                result = self.process_protein(
                    protein['pdb_id'], 
                    protein['chain_id'], 
                    batch_id, 
                    output_dir,
                    import_to_db
                )
                
                if result:
                    results.append({
                        'pdb_id': protein['pdb_id'],
                        'chain_id': protein['chain_id'],
                        'output': result
                    })
            
            # Generate batch summary
            summary_path = os.path.join(output_dir, f"batch_{batch_id}_summary.json")
            with open(summary_path, 'w') as f:
                json.dump({
                    'batch_id': batch_id,
                    'batch_name': batch_info['batch_name'],
                    'ref_version': batch_info['ref_version'],
                    'processed_proteins': len(results),
                    'proteins': results
                }, f, indent=2)
            
            logger.info(f"Processed {len(results)} proteins in batch {batch_id}")
            logger.info(f"Batch summary saved to {summary_path}")
            
            return results
        except Exception as e:
            logger.error(f"Database error in process_batch: {str(e)}")
            return False
        finally:
            cursor.close()

def main():
    """Main function to run the audit trail generator."""
    parser = argparse.ArgumentParser(description='Generate audit trails for domain partitioning in pyECOD')
    parser.add_argument('--config', required=True, help='Path to the configuration file')
    parser.add_argument('--local-config', help='Path to local configuration file (overrides default config.local.yaml)')
    parser.add_argument('--pdb-id', help='PDB ID')
    parser.add_argument('--chain-id', help='Chain ID')
    parser.add_argument('--batch-id', type=int, help='Batch ID')
    parser.add_argument('--output-dir', help='Output directory')
    parser.add_argument('--limit', type=int, help='Limit the number of proteins to process')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--no-db-import', action='store_true', help='Skip importing results to database')
    parser.add_argument('--suspicious-only', action='store_true', help='Only process proteins with suspicious domains')
    parser.add_argument('--min-coverage', type=float, default=50.0, help='Minimum domain coverage threshold (percent)')
    parser.add_argument('--max-short-length', type=int, default=30, help='Maximum length for short alignments')
    parser.add_argument('--schema', help='Database schema name (overrides config value)')
    args = parser.parse_args()
    
    # Set up logging
    if args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug logging enabled")
    
    # Create audit trail object
    audit_trail = DomainAuditTrail(args.config, local_config_path=args.local_config)
    
    # Override schema if provided
    if args.schema:
        audit_trail.schema = args.schema
        logger.info(f"Using schema from command line: {args.schema}")
    
    # Log the actual schema being used
    logger.info(f"Using database schema: {audit_trail.schema}")
    
    # Process suspicious proteins if requested
    if args.suspicious_only and args.batch_id:
        cursor = audit_trail.conn.cursor(cursor_factory=RealDictCursor)
        try:
            logger.info(f"Finding proteins with suspicious domains in batch {args.batch_id}...")
            
            # Find proteins with suspicious domains
            query = f"""
                SELECT p.pdb_id, p.chain_id
                FROM {audit_trail.schema}.protein p
                JOIN {audit_trail.schema}.process_status ps ON p.id = ps.protein_id
                LEFT JOIN {audit_trail.schema}.domain_audit da ON p.id = da.protein_id AND ps.batch_id = da.batch_id
                WHERE ps.batch_id = %s
                AND (
                    da.id IS NULL OR  -- Not yet audited
                    da.min_domain_coverage < %s OR  -- Low coverage
                    da.has_short_alignments = TRUE  -- Has short alignments
                )
                ORDER BY p.pdb_id, p.chain_id
            """
            
            try:
                cursor.execute(query, (args.batch_id, args.min_coverage))
                suspicious_proteins = cursor.fetchall()
            except psycopg2.errors.UndefinedTable:
                # Table domain_audit might not exist yet
                logger.info(f"The domain_audit table was not found in schema {audit_trail.schema}")
                logger.info("Checking for proteins that have not been audited yet")
                
                # Simpler query without domain_audit
                query = f"""
                    SELECT p.pdb_id, p.chain_id
                    FROM {audit_trail.schema}.protein p
                    JOIN {audit_trail.schema}.process_status ps ON p.id = ps.protein_id
                    WHERE ps.batch_id = %s
                    ORDER BY p.pdb_id, p.chain_id
                """
                cursor.execute(query, (args.batch_id,))
                suspicious_proteins = cursor.fetchall()
            
            if not suspicious_proteins:
                logger.info(f"No suspicious proteins found in batch {args.batch_id}")
                return
            
            logger.info(f"Found {len(suspicious_proteins)} proteins with suspicious domains")
            
            # Apply limit if specified
            if args.limit and args.limit < len(suspicious_proteins):
                suspicious_proteins = suspicious_proteins[:args.limit]
                logger.info(f"Processing {args.limit} proteins due to limit option")
            
            # Process each suspicious protein
            processed_count = 0
            for protein in suspicious_proteins:
                logger.info(f"Processing protein {protein['pdb_id']}_{protein['chain_id']}")
                result = audit_trail.process_protein(
                    protein['pdb_id'], 
                    protein['chain_id'], 
                    args.batch_id, 
                    args.output_dir,
                    not args.no_db_import
                )
                
                if result:
                    processed_count += 1
            
            print(f"Processed {processed_count} suspicious proteins from batch {args.batch_id}")
        except Exception as e:
            logger.error(f"Error processing suspicious proteins: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
        finally:
            cursor.close()
        
    elif args.pdb_id and args.chain_id:
        # Process a single protein
        try:
            result = audit_trail.process_protein(
                args.pdb_id, 
                args.chain_id, 
                args.batch_id, 
                args.output_dir,
                not args.no_db_import
            )
            
            if result:
                print(f"Generated audit trail for {args.pdb_id}_{args.chain_id}")
                print(f"Visualization: {result['visualization']}")
                print(f"JSON Report: {result['json_report']}")
                print(f"Text Report: {result['text_report']}")
            else:
                print(f"Failed to generate audit trail for {args.pdb_id}_{args.chain_id}")
        except Exception as e:
            logger.error(f"Error processing protein {args.pdb_id}_{args.chain_id}: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
    
    elif args.batch_id:
        # Process a batch
        try:
            results = audit_trail.process_batch(
                args.batch_id, 
                args.limit, 
                args.output_dir,
                not args.no_db_import
            )
            
            if results:
                print(f"Generated audit trails for {len(results)} proteins in batch {args.batch_id}")
            else:
                print(f"Failed to generate audit trails for batch {args.batch_id}")
        except Exception as e:
            logger.error(f"Error processing batch {args.batch_id}: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
    
    else:
        print("Error: Must provide either --pdb-id and --chain-id or --batch-id")

if __name__ == "__main__":
    main()