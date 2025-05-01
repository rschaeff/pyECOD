"""
HHSearch Processor Module

This module provides functionality for processing HHSearch results (HHR files)
and converting them to XML format for further processing in the ECOD pipeline.

Classes:
    HHRToXMLConverter: Convert HHR parsed data to XML format
    DomainEvidenceCollator: Collate evidence from different sources to create domain suggestions
    HHSearchProcessor: Process HHSearch results and integrate with BLAST evidence
"""
import os
import logging
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import Dict, Any, List, Optional, Tuple

from ecod.utils.hhsearch_utils import HHRParser
from ecod.utils.path_utils import get_standardized_paths, find_files_with_legacy_paths
from ecod.utils.xml_utils import ensure_dict, ensure_list_of_dicts
from ecod.utils.model_conversion import create_domain_summary, xml_to_hits
from ecod.models.pipeline import BlastHit, HHSearchHit



class HHRToXMLConverter:
    """Convert HHR parsed data to XML format"""
    
    def __init__(self, logger=None):
        """Initialize converter with logger"""
        self.logger = logger or logging.getLogger("hhr_to_xml_converter")

    def convert(self, hhr_data: Dict[str, Any], pdb_id: str, chain_id: str, ref_version: str, 
               min_probability: float = 20.0) -> Optional[str]:
        """
        Convert HHR data to XML
        
        Args:
            hhr_data: Parsed HHR data
            pdb_id: PDB ID
            chain_id: Chain ID
            ref_version: Reference version
            min_probability: Minimum probability threshold for including hits (default 20%)
            
        Returns:
            XML string
        """
        try:
            import xml.etree.ElementTree as ET
            from xml.dom import minidom
            from datetime import datetime
            
            self.logger.info(f"Converting HHR data to XML for {pdb_id}_{chain_id}")
            
            # Create root element
            root = ET.Element("hh_summ_doc")
            
            # Add metadata
            metadata = ET.SubElement(root, "metadata")
            ET.SubElement(metadata, "pdb_id").text = pdb_id
            ET.SubElement(metadata, "chain_id").text = chain_id
            ET.SubElement(metadata, "reference").text = ref_version
            ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            ET.SubElement(metadata, "min_probability").text = str(min_probability)
            
            # Add hits
            hits_elem = ET.SubElement(root, "hh_hit_list")
            
            # Filter hits by probability
            filtered_hits = [hit for hit in hhr_data.get('hits', [])
                            if float(hit.get('probability', 0)) >= min_probability]

            self.logger.info(f"Converting {len(filtered_hits)} hits (filtered from {len(hhr_data.get('hits', []))}) with probability >= {min_probability}%")
            
            for hit in filtered_hits:
                hit_elem = ET.SubElement(hits_elem, "hh_hit")
                
                # Add hit attributes
                hit_elem.set("hit_num", str(hit.get('hit_num', 0)))
                hit_elem.set("hit_id", hit.get('hit_id', ''))
                
                # Add probability
                probability = hit.get('probability', 0)
                hit_elem.set("probability", str(probability))
                
                # Add e-value
                e_value = hit.get('e_value', 0)
                hit_elem.set("e_value", str(e_value))
                
                # Add score
                score = hit.get('score', 0)
                hit_elem.set("score", str(score))
                
                # Add ECOD domain ID if hit_id matches pattern
                hit_id = hit.get('hit_id', '')
                if re.match(r'[dge]\d\w{3}\w+\d+', hit_id):
                    hit_elem.set("ecod_domain_id", hit_id)
                
                # Add query range with compact format
                if 'query_range' in hit:
                    query_range_elem = ET.SubElement(hit_elem, "query_range")
                    self._format_range_as_xml(hit['query_range'], query_range_elem)
                
                # Add template range with compact format
                if 'template_range' in hit:
                    template_range_elem = ET.SubElement(hit_elem, "template_seqid_range")
                    self._format_range_as_xml(hit['template_range'], template_range_elem)
                    
                    # Add coverage if available
                    if 'identity_percentage' in hit:
                        template_range_elem.set("identity", str(hit['identity_percentage']))
                    
                    if 'similarity_percentage' in hit:
                        template_range_elem.set("similarity", str(hit['similarity_percentage']))
                    
                    # Add coverage if calculated
                    if 'aligned_cols' in hit and 'length' in hit and hit['length'] > 0:
                        coverage = hit['aligned_cols'] / hit['length']
                        template_range_elem.set("coverage", str(coverage))
                
                # Add alignment details
                alignment_elem = ET.SubElement(hit_elem, "alignment")
                
                # Add alignment sequences if available
                if 'alignment_blocks' in hit and hit['alignment_blocks']:
                    block = hit['alignment_blocks'][0]  # Use first block
                    
                    query_ali_elem = ET.SubElement(alignment_elem, "query_ali")
                    query_ali_elem.text = block.get('query_seq', '')
                    
                    template_ali_elem = ET.SubElement(alignment_elem, "template_ali")
                    template_ali_elem.text = block.get('template_seq', '')
                    
                    if 'match_seq' in block and block['match_seq']:
                        match_elem = ET.SubElement(alignment_elem, "match_sequence")
                        match_elem.text = block['match_seq']
                    
                    if 'consensus_q' in block and block['consensus_q']:
                        consensus_q_elem = ET.SubElement(alignment_elem, "query_consensus")
                        consensus_q_elem.text = block['consensus_q']
                    
                    if 'consensus_t' in block and block['consensus_t']:
                        consensus_t_elem = ET.SubElement(alignment_elem, "template_consensus")
                        consensus_t_elem.text = block['consensus_t']
            
            # Convert to string
            rough_string = ET.tostring(root, 'utf-8')
            reparsed = minidom.parseString(rough_string)
            pretty_xml = reparsed.toprettyxml(indent="  ")
            
            self.logger.info(f"Successfully converted HHR data to XML, size: {len(pretty_xml)} characters")
            
            return pretty_xml
            
        except Exception as e:
            self.logger.error(f"Error converting HHR data to XML: {str(e)}")
            self.logger.exception("Full traceback:")
            return None

    # Improved method for formatting range information in XML
    def _format_range_as_xml(self, range_str: str, parent_elem: ET.Element) -> None:
        """
        Format range string as XML with a more comprehensive structure.

        Args:
            range_str: Range string in format "start-end,start-end,..."
            parent_elem: Parent XML element to add range elements to
        """
        if not range_str or '-' not in range_str:
            self.logger.warning(f"Invalid range string: {range_str}")
            return

        # Set the full range as an attribute for backwards compatibility
        parent_elem.text = range_str

        # Parse and create structured elements for each range segment
        range_parts = range_str.split(',')

        # Add a count of segments
        if len(range_parts) > 1:
            parent_elem.set("segments", str(len(range_parts)))

        # Parse the first segment to get start/end for common use
        if '-' in range_parts[0]:
            try:
                start, end = map(int, range_parts[0].split('-'))
                parent_elem.set("start", str(start))
                parent_elem.set("end", str(end))
            except ValueError:
                self.logger.warning(f"Could not parse range: {range_parts[0]}")
    
    def save(self, xml_string: str, output_path: str) -> bool:
        """
        Save XML string to file
        
        Args:
            xml_string: XML string
            output_path: Output file path
            
        Returns:
            True if successful
        """
        try:
            # Log the output path for debugging
            self.logger.info(f"Attempting to save XML to path: '{output_path}'")
            
            # Validate output path
            if not output_path:
                self.logger.error("Output path is empty or None")
                return False
            
            # Get directory path
            dir_path = os.path.dirname(output_path)
            self.logger.info(f"Directory path: '{dir_path}'")
            
            # Ensure directory exists
            if dir_path:
                os.makedirs(dir_path, exist_ok=True)
                self.logger.info(f"Created/verified directory: {dir_path}")
            else:
                self.logger.warning("No directory specified in output path, saving to current directory")
            
            # Write the file
            with open(output_path, 'w', encoding='utf-8') as f:
                self.logger.info(f"Writing {len(xml_string)} characters to file")
                f.write(xml_string)
            
            # Verify file was created
            if os.path.exists(output_path):
                file_size = os.path.getsize(output_path)
                self.logger.info(f"XML file successfully saved: {output_path} ({file_size} bytes)")
                return True
            else:
                self.logger.error(f"File was not created: {output_path}")
                return False
                
        except Exception as e:
            self.logger.error(f"Error saving XML to '{output_path}': {str(e)}")
            self.logger.exception("Full traceback:")
            return False

class DomainEvidenceCollator:
    """Collate evidence from different sources to create domain suggestions (THIS CLASS HAS BEEN DEPRECATED)"""
    
    def __init__(self, logger):
        self.logger = logger
        self.logger.warning("WARNING! DomainEvidenceCollator has been deprecated! Please move to model-based processing")
    
    def collate(self, chain_blast_dict, domain_blast_dict, hhsearch_dict, pdb_id, chain_id, ref_version):
        """Collate evidence from different sources

        Args:
            chain_blast_dict: Chain-wise BLAST data as dictionary
            domain_blast_dict: Domain-wise BLAST data as dictionary
            hhsearch_dict: HHSearch data as dictionary
            pdb_id: PDB ID
            chain_id: Chain ID
            ref_version: Reference version

        Returns:
            XML string with collated evidence and domain suggestions
        """
        try:
            # Create root element for output
            root = ET.Element("domain_summ_doc")

            # Add metadata
            metadata = ET.SubElement(root, "metadata")
            ET.SubElement(metadata, "pdb_id").text = pdb_id
            ET.SubElement(metadata, "chain_id").text = chain_id
            ET.SubElement(metadata, "reference").text = ref_version
            ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Process chain blast evidence
            chain_blast_elem = ET.SubElement(root, "chain_blast_evidence")
            if chain_blast_dict and 'blast_summ_doc' in chain_blast_dict:
                self._dict_to_xml(chain_blast_elem, chain_blast_dict['blast_summ_doc'])

            # Process domain blast evidence
            domain_blast_elem = ET.SubElement(root, "domain_blast_evidence")
            if domain_blast_dict and 'blast_summ_doc' in domain_blast_dict:
                self._dict_to_xml(domain_blast_elem, domain_blast_dict['blast_summ_doc'])

            # Process HHSearch evidence
            hhsearch_elem = ET.SubElement(root, "hhsearch_evidence")
            if hhsearch_dict and 'hh_summ_doc' in hhsearch_dict:
                self._dict_to_xml(hhsearch_elem, hhsearch_dict['hh_summ_doc'])

            # Generate domain suggestions
            suggestions_elem = ET.SubElement(root, "domain_suggestions")
            self._generate_domain_suggestions(
                suggestions_elem,
                chain_blast_dict,
                domain_blast_dict,
                hhsearch_dict
            )

            # Convert to string
            rough_string = ET.tostring(root, 'utf-8')
            reparsed = minidom.parseString(rough_string)
            pretty_xml = reparsed.toprettyxml(indent="  ")

            return pretty_xml

        except Exception as e:
            self.logger.error(f"Error collating evidence: {str(e)}")
            return None

    def _dict_to_xml(self, parent_elem, data_dict):
        """Convert dictionary to XML elements and append to parent"""
        if not data_dict:
            return

        # Process attributes
        for key, value in data_dict.items():
            if key.startswith('@'):
                # This is an attribute
                parent_elem.set(key[1:], str(value))
            elif key == '_text':
                # This is text content
                parent_elem.text = str(value)
            elif isinstance(value, list):
                # This is a list of child elements
                for item in value:
                    child_elem = ET.SubElement(parent_elem, key)
                    if isinstance(item, dict):
                        self._dict_to_xml(child_elem, item)
                    else:
                        child_elem.text = str(item)
            elif isinstance(value, dict):
                # This is a child element
                child_elem = ET.SubElement(parent_elem, key)
                self._dict_to_xml(child_elem, value)
            else:
                # This is a simple child element with text
                child_elem = ET.SubElement(parent_elem, key)
                child_elem.text = str(value)

    def _generate_domain_suggestions(self, parent, chain_blast_dict, domain_blast_dict, hhsearch_dict):
        """Generate domain suggestions based on evidence

        Args:
            parent: Parent element to append suggestions to
            chain_blast_dict: Chain-wise BLAST data as dictionary
            domain_blast_dict: Domain-wise BLAST data as dictionary
            hhsearch_dict: HHSearch data as dictionary
        """
        # Collect all domain boundaries from different sources
        boundaries = []

        # Extract from chain BLAST
        chain_boundaries = self._extract_boundaries_from_blast(chain_blast_dict)
        boundaries.extend(chain_boundaries)

        # Extract from domain BLAST
        domain_boundaries = self._extract_boundaries_from_blast(domain_blast_dict)
        boundaries.extend(domain_boundaries)

        # Extract from HHSearch
        hh_boundaries = self._extract_boundaries_from_hhsearch(hhsearch_dict)
        boundaries.extend(hh_boundaries)

        # Merge overlapping boundaries
        merged_boundaries = self._merge_boundaries(boundaries)

        # Create domain suggestions
        for i, boundary in enumerate(merged_boundaries):
            domain_elem = ET.SubElement(parent, "domain")
            domain_elem.set("id", f"domain_{i+1}")
            domain_elem.set("start", str(boundary[0]))
            domain_elem.set("end", str(boundary[1]))

            # Add evidence sources
            evidence_elem = ET.SubElement(domain_elem, "evidence")
            for source in boundary[2]:
                source_elem = ET.SubElement(evidence_elem, "source")
                source_elem.text = source

    def _extract_boundaries_from_blast(self, blast_dict):
        """Extract domain boundaries from BLAST data dictionary"""
        boundaries = []

        if not blast_dict:
            return boundaries

        # Find hit elements with query ranges - handle common structures
        hits = []

        # Check for chain blast run
        if 'blast_summ_doc' in blast_dict:
            doc = blast_dict['blast_summ_doc']

            # Check for chain blast pattern
            if 'chain_blast_run' in doc:
                chain_run = doc['chain_blast_run']
                if isinstance(chain_run, list) and len(chain_run) > 0:
                    chain_run = chain_run[0]  # Get first element if list

                if 'hits' in chain_run:
                    hits_container = chain_run['hits']
                    if isinstance(hits_container, list) and len(hits_container) > 0:
                        hits_container = hits_container[0]

                    if 'hit' in hits_container:
                        blast_hits = hits_container['hit']
                        if isinstance(blast_hits, list):
                            hits.extend(blast_hits)
                        else:
                            hits.append(blast_hits)

            # Check for domain blast pattern
            if 'blast_run' in doc:
                blast_run = doc['blast_run']
                if isinstance(blast_run, list) and len(blast_run) > 0:
                    blast_run = blast_run[0]  # Get first element if list

                if 'hits' in blast_run:
                    hits_container = blast_run['hits']
                    if isinstance(hits_container, list) and len(hits_container) > 0:
                        hits_container = hits_container[0]

                    if 'hit' in hits_container:
                        blast_hits = hits_container['hit']
                        if isinstance(blast_hits, list):
                            hits.extend(blast_hits)
                        else:
                            hits.append(blast_hits)

        # Process each hit
        for hit in hits:
            # Find query range
            query_range = None

            # Check for direct query_range field
            if 'query_range' in hit:
                query_range = hit['query_range']
                if isinstance(query_range, list) and len(query_range) > 0:
                    if isinstance(query_range[0], dict) and '_text' in query_range[0]:
                        query_range = query_range[0]['_text']

            # Check for query_reg field
            if not query_range and 'query_reg' in hit:
                query_reg = hit['query_reg']
                if isinstance(query_reg, list) and len(query_reg) > 0:
                    if isinstance(query_reg[0], dict) and '_text' in query_reg[0]:
                        query_range = query_reg[0]['_text']
                    else:
                        query_range = str(query_reg[0])
                else:
                    query_range = str(query_reg)

            if query_range:
                # Parse range in format "start-end,start-end,..."
                for range_str in str(query_range).split(','):
                    if '-' in range_str:
                        try:
                            start, end = map(int, range_str.split('-'))
                            boundaries.append((start, end, ["blast"]))
                        except ValueError:
                            self.logger.warning(f"Could not parse range: {range_str}")

        return boundaries

    def _extract_boundaries_from_hhsearch(self, hhsearch_dict):
        """Extract domain boundaries from HHSearch data dictionary"""
        boundaries = []

        if not hhsearch_dict:
            return boundaries

        # Find hit elements with query ranges - handle common structures
        hits = []

        # Check for hh_summ_doc structure
        if 'hh_summ_doc' in hhsearch_dict:
            doc = hhsearch_dict['hh_summ_doc']

            # Check for hh_hit_list
            if 'hh_hit_list' in doc:
                hit_list = doc['hh_hit_list']
                if isinstance(hit_list, list) and len(hit_list) > 0:
                    hit_list = hit_list[0]  # Get first element if list

                if 'hh_hit' in hit_list:
                    hh_hits = hit_list['hh_hit']
                    if isinstance(hh_hits, list):
                        hits.extend(hh_hits)
                    else:
                        hits.append(hh_hits)

        # Process each hit
        for hit in hits:
            # Find query range
            query_range = None

            # Check for query_range field
            if 'query_range' in hit:
                query_range = hit['query_range']
                if isinstance(query_range, list) and len(query_range) > 0:
                    if isinstance(query_range[0], dict) and '_text' in query_range[0]:
                        query_range = query_range[0]['_text']
                    else:
                        query_range = str(query_range[0])
                else:
                    query_range = str(query_range)

            if query_range:
                # Parse range in format "start-end,start-end,..."
                for range_str in str(query_range).split(','):
                    if '-' in range_str:
                        try:
                            start, end = map(int, range_str.split('-'))
                            boundaries.append((start, end, ["hhsearch"]))
                        except ValueError:
                            self.logger.warning(f"Could not parse range: {range_str}")

        return boundaries

    def _merge_boundaries(self, boundaries):
        """Merge overlapping domain boundaries"""
        if not boundaries:
            return []

        # Sort boundaries by start position
        sorted_boundaries = sorted(boundaries, key=lambda x: x[0])

        merged = []
        current = sorted_boundaries[0]

        for next_boundary in sorted_boundaries[1:]:
            # Check if boundaries overlap
            if next_boundary[0] <= current[1]:
                # Merge boundaries
                current = (
                    current[0],
                    max(current[1], next_boundary[1]),
                    list(set(current[2] + next_boundary[2]))
                )
            else:
                # No overlap, add current to merged list and update current
                merged.append(current)
                current = next_boundary

        # Add the last boundary
        merged.append(current)

        return merged

class HHSearchProcessor:
    """Process HHSearch results and integrate with BLAST evidence"""
    
    def __init__(self, context):
        """Initialize with application context"""
        self.context = context
        self.config = context.config_manager.config
        self.db = context.db
        self.logger = logging.getLogger("ecod.hhsearch_processor")
        
        self.parser = HHRParser(self.logger)
        self.converter = HHRToXMLConverter(self.logger)
        #self.collator = DomainEvidenceCollator(self.logger)
    
    def _get_batch_info(self, batch_id):
        """Get batch information
        
        Args:
            batch_id: Batch ID
            
        Returns:
            Batch information dictionary
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
        
        results = self.db.execute_dict_query(query, (batch_id,))
        if not results:
            return None
        
        return results[0]
    
    def _get_chains_with_hhsearch(self, batch_id):
        """Get chains with completed HHSearch results
        
        Args:
            batch_id: Batch ID
            
        Returns:
            List of chain dictionaries
        """
        # Check database for registered HHR files
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
            AND NOT EXISTS (
                SELECT 1 FROM ecod_schema.process_file 
                WHERE process_id = ps.id AND file_type = 'domain_summary'
            )
        LIMIT 10
        """
        
        self.logger.info(f"Looking for chains with HHR files in batch {batch_id}")
        chains = self.db.execute_dict_query(query, (batch_id,))
        
        if not chains:
            self.logger.warning("No chains found with HHR files in database. This is unexpected since files were registered.")
            self.logger.info("Checking with a simpler query...")
            
            # Try a simpler query to see if any HHR files are registered
            simple_query = """
            SELECT COUNT(*) 
            FROM ecod_schema.process_file
            WHERE file_type = 'hhr'
            AND process_id IN (
                SELECT id FROM ecod_schema.process_status
                WHERE batch_id = %s
            )
            """
            count = self.db.execute_query(simple_query, (batch_id,))
            self.logger.info(f"Found {count[0][0]} HHR files registered for this batch")
            
            # Check what file types are registered
            types_query = """
            SELECT file_type, COUNT(*) 
            FROM ecod_schema.process_file
            WHERE process_id IN (
                SELECT id FROM ecod_schema.process_status
                WHERE batch_id = %s
            )
            GROUP BY file_type
            """
            types = self.db.execute_dict_query(types_query, (batch_id,))
            for t in types:
                self.logger.info(f"File type '{t['file_type']}': {t['count']}")
        
        self.logger.info(f"Found {len(chains)} chains with HHR files")
        return chains
    
    def _get_file_paths(self, batch_info: Dict, pdb_id: str, chain_id: str,
                       ref_version: str
    ) -> Dict[str, str]:
        """Get standardized file paths for a protein chain"""
        self.logger.info(f"===== _GET_FILE_PATHS CALLED! =====")
        from ecod.utils.path_utils import get_standardized_paths, find_files_with_legacy_paths

        # Get standardized paths
        paths = get_standardized_paths(batch_info['base_path'], pdb_id, chain_id, ref_version)

        #Debug logging
        self.logger.info(f"========== STANDARDIZED PATHS ==========")
        self.logger.info(f"PDB ID: {pdb_id}, Chain ID: {chain_id}")
        self.logger.info(f"Chain BLAST standardized path: {paths['chain_blast']}")
        self.logger.info(f"Domain BLAST standardized path: {paths['domain_blast']}")

        # Check if files exist at standard paths, if not check legacy paths
        legacy_files = find_files_with_legacy_paths(batch_info['base_path'], pdb_id, chain_id, ref_version)

        #Debug logging
        self.logger.info(f"========== LEGACY FILES ==========")
        if 'chain_blast' in legacy_files:
            self.logger.info(f"Chain BLAST legacy path: {legacy_files['chain_blast'].get('exists_at', 'Not found')}")
        if 'domain_blast' in legacy_files:
            self.logger.info(f"Domain BLAST legacy path: {legacy_files['domain_blast'].get('exists_at', 'Not found')}")

        # For each file type, use 'exists_at' from legacy_files if the file doesn't exist at the standard path
        for file_type in paths:
            if file_type in legacy_files and not os.path.exists(paths[file_type]) and legacy_files[file_type]['exists_at']:
                paths[file_type] = legacy_files[file_type]['exists_at']
        
        return paths

    def _parse_xml(self, xml_path):
        """Parse XML file
        
        Args:
            xml_path: Path to XML file
            
        Returns:
            ElementTree object or None if failed
        """
        try:
            return ET.parse(xml_path)
        except Exception as e:
            self.logger.error(f"Error parsing XML file {xml_path}: {str(e)}")
            return None
    
    def _register_file(self, process_id, file_type, file_path, file_size):
        """Register file in database
        
        Args:
            process_id: Process ID
            file_type: File type
            file_path: File path
            file_size: File size
            
        Returns:
            True if successful
        """
        try:
            # Check if file already registered
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """
            
            existing = self.db.execute_query(query, (process_id, file_type))
            
            if existing:
                # Update existing record
                self.db.update(
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
                # Create new record
                self.db.insert(
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
            self.logger.error(f"Error registering file: {str(e)}")
            return False
    
    def _update_process_status(self, process_id, stage, error_message=None):
        """Update process status in database
        
        Args:
            process_id: Process ID
            stage: Current stage
            error_message: Optional error message
            
        Returns:
            True if successful
        """
        try:
            status = "error" if error_message else "success"
            
            self.db.update(
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
            self.logger.error(f"Error updating process status: {str(e)}")
            return False
    
    def process_batch(self, batch_id, force=False):
        """Process HHSearch results for a batch
        
        Args:
            batch_id: Batch ID to process
            force: Force reprocessing of already processed results
            
        Returns:
            Number of successfully processed chains
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 0
        
        # Get chains with completed HHSearch results
        chains = self._get_chains_with_hhsearch(batch_id)
        if not chains:
            self.logger.warning(f"No chains with completed HHSearch results found in batch {batch_id}")
            return 0
        
        self.logger.info(f"Found {len(chains)} chains with HHSearch results to process")
        
        processed_count = 0
        for chain in chains:
            success = self._process_chain(
                chain['pdb_id'], 
                chain['chain_id'], 
                chain['process_id'],
                batch_info,
                batch_info['ref_version']
            )
            
            if success:
                processed_count += 1
                
        self.logger.info(f"Successfully processed {processed_count} out of {len(chains)} chains")
        return processed_count
    
    def process_protein(self, protein_id, force=False):
        """Process HHSearch results for a specific protein
        
        Args:
            protein_id: Protein ID
            force: Force reprocessing of already processed results
            
        Returns:
            True if successful
        """
        # Get protein and process information
        query = """
        SELECT 
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            ps.batch_id,
            b.base_path,
            b.ref_version
        FROM 
            ecod_schema.protein p
        JOIN
            ecod_schema.process_status ps ON p.id = ps.protein_id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            p.id = %s
        """
        
        results = self.db.execute_dict_query(query, (protein_id,))
        if not results:
            self.logger.error(f"Protein {protein_id} not found")
            return False
        
        protein_info = results[0]
        
        # Check if HHR file exists
        batch_info = {
            'id': protein_info['batch_id'],
            'base_path': protein_info['base_path'],
            'ref_version': protein_info['ref_version']
        }
        
        paths = self._get_file_paths(
            batch_info, 
            protein_info['pdb_id'], 
            protein_info['chain_id'], 
            protein_info['ref_version']
        )
        
        if not os.path.exists(paths['hhr']):
            self.logger.error(f"HHR file not found: {paths['hhr']}")
            return False
        
        # Process the chain
        return self._process_chain(
            protein_info['pdb_id'],
            protein_info['chain_id'],
            protein_info['process_id'],
            batch_info,
            protein_info['ref_version'],
            force
        )
    
    def _process_chain(self, pdb_id, chain_id, process_id, batch_info, ref_version, force=False):
        """Process HHSearch results for a chain"""
        try:
            self.logger.info(f"===== _PROCESS_CHAIN CALLED! =====")

            # Get file paths
            #paths = self._get_file_paths(batch_info, pdb_id, chain_id, ref_version)
            paths = get_all_evidence_paths(batch_info['base_path'], pdb_id, chain_id, ref_version)

            # Debug logging
            self.logger.info(f"========== DEBUG PATHS ==========")
            self.logger.info(f"PDB ID: {pdb_id}, Chain ID: {chain_id}")
            self.logger.info(f"Batch path: {batch_info['base_path']}")
            self.logger.info(f"Chain BLAST path: {paths['chain_blast']}")
            self.logger.info(f"Chain BLAST exists: {os.path.exists(paths['chain_blast'])}")
            self.logger.info(f"Domain BLAST path: {paths['domain_blast']}")
            self.logger.info(f"Domain BLAST exists: {os.path.exists(paths['domain_blast'])}")
            self.logger.info(f"=================================")

            # Create domain summary model from evidence
            summary = create_domain_summary(pdb_id, chain_id, ref_version, paths)

            # Check if domain summary already exists
            domain_summary_path = paths['domain_summary']['standard_path']

            if os.path.exists(domain_summary_path) and not force:

                # Register domain summary in database
                self._register_file(
                    process_id,
                    "domain_summary",
                    paths['domain_summary'],
                    os.path.getsize(paths['domain_summary'])
                )

                # Update process status
                self._update_process_status(process_id, "domain_summary_complete")

                return True

            # Parse HHR file
            # Process HHR file if needed
            hhr_path = paths['hhr']['exists_at']
            if not paths['hh_xml']['exists_at'] and hhr_path:
                self.logger.error(f"Failed to parse HHR file: {paths['hhr']}")
                return False

            # Convert to XML
            self.logger.info(f"Converting HHR data to XML for {pdb_id}_{chain_id}")
            xml_string = self.converter.convert(hhr_data, pdb_id, chain_id, ref_version)
            if not xml_string:
                self.logger.error(f"Failed to convert HHR data to XML for {pdb_id}_{chain_id}")
                return False

            # Save HHSearch XML
            self.logger.info(f"Saving HHSearch XML: {paths['hh_xml']}")
            if not self.converter.save(xml_string, paths['hh_xml']):
                self.logger.error(f"Failed to save HHSearch XML: {paths['hh_xml']}")
                return False

            # Register HHSearch XML in database
            self._register_file(
                process_id,
                "hhsearch_xml",
                paths['hh_xml'],
                os.path.getsize(paths['hh_xml'])
            )

            # Parse XML files to dictionaries
            chain_blast_dict = None
            domain_blast_dict = None
            hhsearch_dict = None

            # Process Chain BLAST hits (if available) using models
            if os.path.exists(paths['chain_blast']):
                try:
                    self.logger.info(f"Processing chain BLAST: {paths['chain_blast']}")
                    summary.chain_blast_hits = xml_to_hits(
                        paths['chain_blast'], "chain_blast"
                    )

                    # Ensure hit_type is set correctly
                    for hit in summary.chain_blast_hits:
                        hit.hit_type = "chain_blast"

                except Exception as e:
                    self.logger.error(f"Error processing chain BLAST: {e}")
                    summary.errors["chain_blast_error"] = True
            else:
                self.logger.warning(f"Chain BLAST file not found: {paths['chain_blast']}")
                summary.errors["no_chain_blast"] = True

            # Process Domain BLAST hits (if available) using models
            if os.path.exists(paths['domain_blast']):
                try:
                    self.logger.info(f"Processing domain BLAST: {paths['domain_blast']}")
                    summary.domain_blast_hits = xml_to_hits(
                        paths['domain_blast'], "domain_blast"
                    )

                    # Ensure hit_type is set correctly
                    for hit in summary.domain_blast_hits:
                        hit.hit_type = "domain_blast"

                except Exception as e:
                    self.logger.error(f"Error processing domain BLAST: {e}")
                    summary.errors["domain_blast_error"] = True
            else:
                self.logger.warning(f"Domain BLAST file not found: {paths['domain_blast']}")
                summary.errors["no_domain_blast"] = True

            # Process HHSearch hits (if available) using models
            if os.path.exists(paths['hh_xml']):
                try:
                    self.logger.info(f"Processing HHSearch results: {paths['hh_xml']}")
                    summary.hhsearch_hits = xml_to_hits(
                        paths['hh_xml'], "hhsearch"
                    )
                except Exception as e:
                    self.logger.error(f"Error processing HHSearch results: {e}")
                    summary.errors["hhsearch_error"] = True
            else:
                self.logger.warning(f"HHSearch XML file not found: {paths['hh_xml']}")
                summary.errors["no_hhsearch"] = True

            # Generate XML and save domain summary
            self.logger.info(f"Generating domain summary for {pdb_id}_{chain_id}")
            root = summary.to_xml()
            
            if not summary_xml:
                self.logger.error(f"Failed to collate evidence for {pdb_id}_{chain_id}")
                return False
            
            # Save domain summary
            self.logger.info(f"Saving domain summary: {paths['domain_summary']}")
            root = summary.to_xml()
            save_xml(root, domain_summary_path)

            # Register in database
            self._register_file(
                process_id,
                "domain_summary",
                domain_summary_path,
                os.path.getsize(domain_summary_path)
            )
            
            # Update process status
            self._update_process_status(process_id, "domain_summary_complete")
            
            self.logger.info(f"Successfully processed HHSearch results for {pdb_id}_{chain_id}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error processing chain {pdb_id}_{chain_id}: {str(e)}")
            self._update_process_status(process_id, "error", str(e))
            return False
