#!/usr/bin/env python3
"""
HHSearch result registration module for ECOD pipeline
Handles HHR file processing, XML conversion, and database registration
"""

import os
import logging
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
import re
from typing import Dict, Any, List, Optional, Tuple

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.exceptions import PipelineError, FileOperationError


class HHResultRegistrar:
    """Register HHSearch results and convert to XML for domain analysis"""
    
    def __init__(self, context=None):
        """Initialize with application context"""
        self.context = context or ApplicationContext()
        self.db = self.context.db
        self.config = self.context.config_manager.config
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.hhresult_registrar")
    
    def register_batch_results(self, batch_id: int, force_regenerate: bool = False) -> int:
        """
        Find, convert, and register HHSearch results for a batch
        
        Args:
            batch_id: Batch ID to process
            force_regenerate: Force regeneration of XML files
            
        Returns:
            Number of registered files
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 0
        
        # Get chains with HHR files
        chains = self._find_chains_with_hhr(batch_id)
        if not chains:
            self.logger.warning(f"No chains with HHR files found in batch {batch_id}")
            return 0
        
        # Register files
        registered_count = 0
        for chain in chains:
            if self._register_chain_results(
                chain['process_id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_info['ref_version'],
                batch_info['base_path'],
                force_regenerate
            ):
                registered_count += 1
                
        self.logger.info(f"Registered HHR results for {registered_count} out of {len(chains)} chains")
        return registered_count
    
    def register_specific_chains(self, batch_id: int, chain_ids: List[str], 
                               force_regenerate: bool = False) -> int:
        """
        Find, convert, and register HHSearch results for specific chains
        
        Args:
            batch_id: Batch ID
            chain_ids: List of chain IDs to process (format: "pdbid_chainid")
            force_regenerate: Force regeneration of XML files
            
        Returns:
            Number of registered files
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 0
        
        # Get specific chains
        chains = []
        for chain_id in chain_ids:
            try:
                pdb_id, chain_letter = chain_id.split('_')
                result = self._find_specific_chain(batch_id, pdb_id, chain_letter)
                if result:
                    chains.append(result)
            except ValueError:
                self.logger.warning(f"Invalid chain ID format: {chain_id}, expected pdbid_chainid")
        
        if not chains:
            self.logger.warning(f"No specified chains found in batch {batch_id}")
            return 0
        
        # Register files
        registered_count = 0
        for chain in chains:
            if self._register_chain_results(
                chain['process_id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_info['ref_version'],
                batch_info['base_path'],
                force_regenerate
            ):
                registered_count += 1
                
        self.logger.info(f"Registered HHR results for {registered_count} out of {len(chains)} chains")
        return registered_count
    
    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information"""
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
            self.logger.error(f"Error retrieving batch information: {str(e)}")
            
        return None
    
    def _find_chains_with_hhr(self, batch_id: int) -> List[Dict[str, Any]]:
        """Find chains with HHR files in a batch"""
        # First attempt to find chains with registered HHR files
        query = """
        SELECT 
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            pf.file_path as hhr_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 
            ps.batch_id = %s
            AND pf.file_type = 'hhr'
            AND pf.file_exists = TRUE
        """
        
        chains = self.db.execute_dict_query(query, (batch_id,))
        if chains:
            self.logger.info(f"Found {len(chains)} chains with registered HHR files")
            return chains
        
        # If no registered HHR files, search on filesystem using process entries
        self.logger.info("No registered HHR files found, searching filesystem...")
        
        query = """
        SELECT 
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            b.base_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            ps.batch_id = %s
            AND ps.status = 'success'
        """
        
        result = []
        chains = self.db.execute_dict_query(query, (batch_id,))
        
        if not chains:
            return []
            
        # Check if HHR files exist in filesystem
        ref_version = self._get_batch_ref_version(batch_id)
        for chain in chains:
            # Check standard locations
            potential_hhr_paths = [
                # New standard location (flat structure in hhsearch dir)
                os.path.join(chain['base_path'], "hhsearch", f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr"),
                
                # Old chain-specific directory structure
                os.path.join(chain['base_path'], chain['relative_path'], f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr"),
                
                # Other potential locations
                os.path.join(chain['base_path'], f"{chain['pdb_id']}_{chain['chain_id']}", f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr")
            ]
            
            for path in potential_hhr_paths:
                if os.path.exists(path) and os.path.getsize(path) > 0:
                    chain['hhr_path'] = path
                    result.append(chain)
                    break
        
        self.logger.info(f"Found {len(result)} chains with unregistered HHR files")
        return result
    
    def _find_specific_chain(self, batch_id: int, pdb_id: str, chain_id: str) -> Optional[Dict[str, Any]]:
        """Find a specific chain in a batch"""
        query = """
        SELECT 
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            b.base_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            ps.batch_id = %s
            AND p.pdb_id = %s
            AND p.chain_id = %s
        LIMIT 1
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id, pdb_id, chain_id))
            if rows:
                chain = rows[0]
                
                # Check for HHR file
                ref_version = self._get_batch_ref_version(batch_id)
                potential_hhr_paths = [
                    # New standard location (flat structure in hhsearch dir)
                    os.path.join(chain['base_path'], "hhsearch", f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr"),
                    
                    # Old chain-specific directory structure
                    os.path.join(chain['base_path'], chain['relative_path'], f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr"),
                    
                    # Other potential locations
                    os.path.join(chain['base_path'], f"{chain['pdb_id']}_{chain['chain_id']}", f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr")
                ]
                
                for path in potential_hhr_paths:
                    if os.path.exists(path) and os.path.getsize(path) > 0:
                        chain['hhr_path'] = path
                        return chain
                
                self.logger.warning(f"No HHR file found for {pdb_id}_{chain_id}")
                return None
        except Exception as e:
            self.logger.error(f"Error finding chain {pdb_id}_{chain_id}: {str(e)}")
            
        return None
    
    def _get_batch_ref_version(self, batch_id: int) -> str:
        """Get reference version for a batch"""
        query = "SELECT ref_version FROM ecod_schema.batch WHERE id = %s"
        
        try:
            rows = self.db.execute_query(query, (batch_id,))
            if rows:
                return rows[0][0]
        except Exception as e:
            self.logger.error(f"Error getting reference version: {str(e)}")
            
        # Default to configured reference version
        return self.config.get('reference', {}).get('current_version', 'develop291')
    
    def _register_chain_results(self, process_id: int, pdb_id: str, chain_id: str, 
                              ref_version: str, base_path: str, 
                              force_regenerate: bool = False) -> bool:
        """Register HHSearch results for a chain"""
        try:
            # Define standard paths
            pdb_chain = f"{pdb_id}_{chain_id}"
            hhsearch_dir = os.path.join(base_path, "hhsearch")
            os.makedirs(hhsearch_dir, exist_ok=True)
            
            hhr_file = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhr")
            xml_file = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml")
            
            # Check if XML file already exists (unless force regenerate)
            if os.path.exists(xml_file) and os.path.getsize(xml_file) > 0 and not force_regenerate:
                self.logger.info(f"XML file already exists for {pdb_chain}, using existing file")
                
                # Register XML file if needed
                self._check_and_register_file(process_id, 'hhsearch_xml', xml_file, base_path)
                
                # Update process status
                self._update_process_status(process_id, "hhsearch_complete")
                
                return True
            
            # Check if HHR file exists
            if not os.path.exists(hhr_file) or os.path.getsize(hhr_file) == 0:
                self.logger.warning(f"HHR file missing or empty: {hhr_file}")
                
                # Try to find HHR file in alternative locations
                found_hhr = self._find_and_move_hhr(pdb_id, chain_id, ref_version, base_path, hhr_file)
                if not found_hhr:
                    return False
            
            # Register HHR file
            self._check_and_register_file(process_id, 'hhr', hhr_file, base_path)
            
            # Convert HHR to XML
            self.logger.info(f"Converting {hhr_file} to XML...")
            parser = HHRParser(self.logger)
            hhr_data = parser.parse(hhr_file)
            
            if not hhr_data:
                self.logger.error(f"Failed to parse HHR file: {hhr_file}")
                return False
            
            # Convert to XML
            converter = HHRToXMLConverter(self.logger)
            xml_string = converter.convert(hhr_data, pdb_id, chain_id, ref_version)
            
            if not xml_string:
                self.logger.error(f"Failed to convert HHR data to XML for {pdb_chain}")
                return False
            
            # Save XML file
            with open(xml_file, 'w', encoding='utf-8') as f:
                f.write(xml_string)
            
            # Register XML file
            self._check_and_register_file(process_id, 'hhsearch_xml', xml_file, base_path)
            
            # Update process status
            self._update_process_status(process_id, "hhsearch_complete")
            
            self.logger.info(f"Successfully registered HHSearch results for {pdb_chain}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error registering results for {pdb_id}_{chain_id}: {str(e)}")
            return False
    
    def _find_and_move_hhr(self, pdb_id: str, chain_id: str, ref_version: str, 
                         base_path: str, target_path: str) -> bool:
        """Find HHR file in alternative locations and move to standard location"""
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Check alternative locations
        alternative_paths = [
            os.path.join(base_path, pdb_chain, f"{pdb_chain}.{ref_version}.hhr"),
            os.path.join(base_path, "ecod_dump", pdb_chain, f"{pdb_chain}.{ref_version}.hhr"),
            os.path.join(base_path, "ecod_dump", f"{pdb_chain}.{ref_version}.hhr")
        ]
        
        for path in alternative_paths:
            if os.path.exists(path) and os.path.getsize(path) > 0:
                self.logger.info(f"Found HHR file at alternative location: {path}")
                
                # Copy file to standard location
                try:
                    os.makedirs(os.path.dirname(target_path), exist_ok=True)
                    
                    # Use shutil.copy2 to preserve metadata
                    import shutil
                    shutil.copy2(path, target_path)
                    
                    self.logger.info(f"Copied HHR file to standard location: {target_path}")
                    return True
                except Exception as e:
                    self.logger.error(f"Error copying HHR file: {str(e)}")
                    return False
        
        self.logger.error(f"Could not find HHR file for {pdb_chain} in any location")
        return False
    
    def _check_and_register_file(self, process_id: int, file_type: str, 
                                file_path: str, base_path: str) -> bool:
        """Check if file is registered in database and register if not"""
        try:
            # Check if file already registered
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """
            
            existing = self.db.execute_query(query, (process_id, file_type))
            
            # Get relative path (for database storage)
            rel_path = os.path.relpath(file_path, base_path)
            
            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path)
                    },
                    "id = %s",
                    (existing[0][0],)
                )
                self.logger.info(f"Updated {file_type} record for process {process_id}")
            else:
                # Create new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path)
                    }
                )
                self.logger.info(f"Created new {file_type} record for process {process_id}")
            
            return True
        except Exception as e:
            self.logger.error(f"Error registering file {file_path}: {str(e)}")
            return False
    
    def _update_process_status(self, process_id: int, stage: str, error_message: str = None) -> bool:
        """Update process status in database"""
        try:
            status = "error" if error_message else "success"
            
            update_data = {
                "current_stage": stage,
                "status": status
            }
            
            if error_message:
                update_data["error_message"] = error_message
                
            self.db.update(
                "ecod_schema.process_status",
                update_data,
                "id = %s",
                (process_id,)
            )
            
            return True
        except Exception as e:
            self.logger.error(f"Error updating process status for {process_id}: {str(e)}")
            return False


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
        if current_hit and current_hit['query_ali'] and current_hit['template_ali']:
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


def register_hhsearch_results(context, batch_id, force_regenerate=False):
    """
    Find HHR files and register them in the database, converting to XML if needed
    
    Args:
        context: Application context
        batch_id: Batch ID to process
        force_regenerate: Force regeneration of XML files
        
    Returns:
        Number of registered files
    """
    registrar = HHResultRegistrar(context)
    return registrar.register_batch_results(batch_id, force_regenerate)