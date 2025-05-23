#!/usr/bin/env python3
"""
Domain summary module for the ECOD pipeline
Processes and integrates BLAST and HHSearch results for domain analysis
"""

import os
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.exceptions import PipelineError, FileOperationError
from ecod.models import BlastHit, HHSearchHit, DomainSummaryModel
from ecod.utils.file import read_sequence_from_fasta, find_fasta_file

class DomainSummary:
    """Process and integrate BLAST and HHSearch results for a protein chain"""
    
    def __init__(self, context=None):
        """Initialize with configuration"""
        self.context = context
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.summary")
        
        # Set default thresholds
        self.hsp_evalue_threshold = 0.005
        self.hit_coverage_threshold = 0.7

    def simplified_file_path_resolution(self, pdb_id, chain_id, file_type, job_dump_dir):
        """Simplified method to locate files, handling absolute paths"""
        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)
        
        # Log what we're looking for
        self.logger.debug(f"Looking for {file_type} for {pdb_id}_{chain_id}")
        
        # Special handling for HHSearch file types
        if file_type in ['hhr', 'hhsearch_xml', 'hhsearch_result']:
            reference = self.context.config_manager.config.get('reference', {}).get('current_version', 'develop291')
            
            # Only look for standard naming pattern
            standard_pattern = f"{pdb_id}_{chain_id}.{reference}.hhsearch.xml"
            standard_path = os.path.join(job_dump_dir, "hhsearch", standard_pattern)
            
            if os.path.exists(standard_path) and os.path.getsize(standard_path) > 0:
                self.logger.info(f"Found standard HHSearch XML: {standard_path}")
                return [standard_path]
            else:
                self.logger.warning(f"No standard HHSearch XML found: {standard_path}")
                return []
        
        # Query database for file path
        query = """
        SELECT pf.file_path
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        AND pf.file_type = %s
        AND pf.file_exists = TRUE
        ORDER BY pf.id DESC
        LIMIT 1
        """
        
        try:
            rows = db.execute_query(query, (pdb_id, chain_id, file_type))
            if rows:
                # Get the path from database
                db_file_path = rows[0][0]
                self.logger.info(f"Found {file_type} file from database: {db_file_path}")
                
                # Use path as-is if it's absolute, otherwise combine with job_dump_dir
                if os.path.isabs(db_file_path):
                    full_path = db_file_path
                else:
                    full_path = os.path.join(job_dump_dir, db_file_path)
                
                # Normalize the path to resolve any '..' components
                full_path = os.path.normpath(full_path)
                
                if os.path.exists(full_path):
                    self.logger.info(f"File exists at path: {full_path}")
                    return [full_path]  # Return as a list for compatibility
                else:
                    self.logger.warning(f"File exists in database but not on filesystem: {full_path}")
        except Exception as e:
            self.logger.error(f"Error querying database for {file_type} file: {e}")
        
        # File not found or database query failed
        self.logger.error(f"No {file_type} file found for {pdb_id}_{chain_id}")
        return []  # Return empty list to maintain expected return type    

    # Add this helper method to check if a BLAST file has hits:
    def _check_blast_has_hits(self, blast_path):
        """Check if a BLAST XML file contains any hits"""
        try:
            tree = ET.parse(blast_path)
            root = tree.getroot()
            
            # Find all Hit elements
            hits = root.findall(".//Hit")
            if len(hits) == 0:
                self.logger.warning(f"BLAST file has no hits: {blast_path}")
                return False
            return True
        except Exception as e:
            self.logger.error(f"Error checking BLAST hits: {e}")
            return False

    def create_summary(self, pdb_id: str, chain_id: str, reference: str, 
                      job_dump_dir: str, blast_only: bool = False) -> DomainSummaryModel:
        """Create domain summary for a protein chain

        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            reference: Reference version
            job_dump_dir: Base path for files
            blast_only: Whether to use only BLAST results (no HHSearch)

        Returns:
            DomainSummaryModel containing summary data and file information
        """
        # Create a DomainSummary model object
        summary = DomainSummaryModel(
            pdb_id=pdb_id,
            chain_id=chain_id,
            reference=reference
        )

        # Define output path
        domains_dir = os.path.join(job_dump_dir, "domains")
        os.makedirs(domains_dir, exist_ok=True)
        suffix = ".blast_only" if blast_only else ""
        output_filename = f"{pdb_id}_{chain_id}.{reference}.domain_summary{suffix}.xml"
        output_path = os.path.join(domains_dir, output_filename)

        # Add file information to the model
        summary.output_file_path = output_path
        summary.processed = False
        summary.skipped = False

        # Check for existing file
        if os.path.exists(output_path) and not self.context.is_force_overwrite():
            self.logger.warning(f"Output file {output_path} already exists, skipping...")
            summary.skipped = True
            return summary

        # Get sequence from FASTA
        sequence = read_fasta_sequence(find_fasta_file(pdb_id, chain_id, job_dump_dir))
        summary.sequence_length = len(sequence) if sequence else 0
        summary.sequence = sequence

        # Special handling for peptides
        if sequence and len(sequence) < 30:
            summary.is_peptide = True
            # Write output file
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            root = summary.to_xml()
            tree = ET.ElementTree(root)
            tree.write(output_path, encoding='utf-8', xml_declaration=True)
            summary.processed = True
            return summary

        # Process self-comparison if available
        self_comp_path = self._find_self_comparison(pdb_id, chain_id, job_dump_dir)
        if not self_comp_path:
            summary.errors["no_selfcomp"] = True
        else:
            summary.self_comparison_hits = self._process_self_comparison_to_dict(self_comp_path)

        # Process chain BLAST results
        chain_blast_paths = self.simplified_file_path_resolution(
            pdb_id, chain_id, 'chain_blast_result', job_dump_dir
        )

        if not chain_blast_paths:
            summary.errors["no_chain_blast"] = True
        else:
            chain_blast_file = chain_blast_paths[0]
            try:
                if self._check_blast_has_hits(chain_blast_file):
                    chain_hits = self._process_chain_blast_to_dict(chain_blast_file)
                    summary.chain_blast_hits = chain_hits
                else:
                    summary.errors["chain_blast_no_hits"] = True
            except Exception as e:
                self.logger.error(f"Error processing chain BLAST: {e}")
                summary.errors["chain_blast_error"] = True

        # Process domain BLAST results
        domain_blast_paths = self.simplified_file_path_resolution(
            pdb_id, chain_id, 'domain_blast_result', job_dump_dir
        )

        if not domain_blast_paths:
            summary.errors["no_domain_blast"] = True
        else:
            domain_blast_file = domain_blast_paths[0]
            try:
                if self._check_blast_has_hits(domain_blast_file):
                    domain_hits = self._process_blast(domain_blast_file)
                    summary.domain_blast_hits = domain_hits
                else:
                    summary.errors["domain_blast_no_hits"] = True
            except Exception as e:
                self.logger.error(f"Error processing domain BLAST: {e}")
                summary.errors["domain_blast_error"] = True

        # Process HHSearch results (skip if blast_only mode)
        if not blast_only:
            standard_pattern = f"{pdb_id}_{chain_id}.{reference}.hhsearch.xml"
            hhsearch_path = os.path.join(job_dump_dir, "hhsearch", standard_pattern)

            if os.path.exists(hhsearch_path) and os.path.getsize(hhsearch_path) > 0:
                try:
                    hhsearch_hits = self._process_hhsearch_to_dict(hhsearch_path)
                    summary.hhsearch_hits = hhsearch_hits
                except Exception as e:
                    self.logger.error(f"Error processing HHSearch: {e}")
                    summary.errors["hhsearch_error"] = True
            else:
                summary.errors["no_hhsearch"] = True

        # Convert to XML and write output file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        root = summary.to_xml()
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding='utf-8', xml_declaration=True)
        summary.processed = True

        return summary

    def _read_fasta_sequence(self, fasta_path: str) -> Optional[str]:
        """Read sequence from a FASTA file"""
        if not os.path.exists(fasta_path):
            return None
        
        try:
            with open(fasta_path, 'r') as f:
                lines = f.readlines()
            
            # Skip header line
            sequence = ""
            for line in lines[1:]:
                sequence += line.strip()
            
            return sequence
        except Exception as e:
            self.logger.error(f"Error reading FASTA file {fasta_path}: {e}")
            return None

    def _process_self_comparison(self, self_comp_path: str, parent_node: ET.Element) -> None:
        """Process structural self-comparison results"""
        try:
            tree = ET.parse(self_comp_path)
            root = tree.getroot()
            
            # Create self-comparison section
            self_comp_run = ET.SubElement(parent_node, "self_comp_run")
            self_comp_run.set("programs", "dali")
            
            hits = ET.SubElement(self_comp_run, "hits")
            
            # Process structural repeats
            for repeat in root.findall(".//structural_repeat"):
                aligner = repeat.get("aligner", "")
                zscore = repeat.get("zscore", "")
                
                ref_range = repeat.findtext("ref_range", "")
                mob_range = repeat.findtext("mob_range", "")
                
                hit = ET.SubElement(hits, "hit")
                hit.set("aligner", aligner)
                hit.set("z_score", zscore)
                
                query_reg = ET.SubElement(hit, "query_reg")
                query_reg.text = ref_range
                
                hit_reg = ET.SubElement(hit, "hit_reg")
                hit_reg.text = mob_range
            
            # Process sequence repeats
            seq_comp_run = ET.SubElement(parent_node, "self_comp_run")
            seq_comp_run.set("programs", "hhrepid")
            
            repeat_set_list = ET.SubElement(seq_comp_run, "repeat_set_list")
            
            for repeat_set in root.findall(".//sequence_repeat_set"):
                aligner = repeat_set.get("aligner", "")
                type_val = repeat_set.get("type", "")
                
                # Collect all ranges in this set
                set_ranges = []
                for seq_repeat in repeat_set.findall("sequence_repeat"):
                    range_val = seq_repeat.findtext("range", "")
                    set_ranges.append(range_val)
                
                # Create repeat set element
                set_elem = ET.SubElement(repeat_set_list, "repeat_set")
                set_elem.set("aligner", aligner)
                set_elem.set("type", type_val)
                
                seqid_range = ET.SubElement(set_elem, "seqid_range")
                seqid_range.text = ",".join(set_ranges)
            
        except Exception as e:
            self.logger.error(f"Error processing self comparison: {e}")
    
    def _process_chain_blast(self, chain_blast_path: str, parent_node: ET.Element) -> None:
        """Process chain-wise BLAST results"""
        try:
            tree = ET.parse(chain_blast_path)
            root = tree.getroot()
            
            # Create BLAST run section
            blast_run = ET.SubElement(parent_node, "chain_blast_run")
            
            program = root.findtext(".//BlastOutput_program", "")
            blast_run.set("program", program)
            
            version = root.findtext(".//BlastOutput_version", "")
            blast_run.set("version", version)
            
            # Add database info
            db = root.findtext(".//BlastOutput_db", "")
            db_node = ET.SubElement(blast_run, "blast_db")
            db_node.text = db
            
            # Add query info
            query = root.findtext(".//BlastOutput_query-def", "")
            query_node = ET.SubElement(blast_run, "blast_query")
            query_node.text = query
            
            query_len = root.findtext(".//BlastOutput_query-len", "")
            query_len_node = ET.SubElement(blast_run, "query_len")
            query_len_node.text = query_len
            
            # Process hits
            hits_node = ET.SubElement(blast_run, "hits")
            
            for iteration in root.findall(".//Iteration"):
                for hit in iteration.findall(".//Hit"):
                    hit_num = hit.findtext("Hit_num", "")
                    hit_def = hit.findtext("Hit_def", "")
                    hit_len = int(hit.findtext("Hit_len", "0"))
                    
                    # Parse PDB ID and chain from hit definition
                    hit_pdb = "NA"
                    hit_chain = "NA"
                    if " " in hit_def:
                        parts = hit_def.split()
                        if len(parts) >= 2:
                            hit_pdb = parts[0]
                            hit_chain = parts[1]
                    
                    # Process HSPs
                    valid_hsps = self._process_hit_hsps(hit, int(query_len))
                    
                    if valid_hsps:
                        # Create hit element
                        hit_elem = ET.SubElement(hits_node, "hit")
                        hit_elem.set("num", hit_num)
                        hit_elem.set("pdb_id", hit_pdb)
                        hit_elem.set("chain_id", hit_chain)
                        hit_elem.set("hsp_count", str(len(valid_hsps["query_regions"])))
                        hit_elem.set("evalues", ",".join(map(str, valid_hsps["evalues"])))
                        
                        # Add query regions
                        query_reg = ET.SubElement(hit_elem, "query_reg")
                        query_reg.text = ",".join(valid_hsps["query_regions"])
                        
                        # Add hit regions
                        hit_reg = ET.SubElement(hit_elem, "hit_reg")
                        hit_reg.text = ",".join(valid_hsps["hit_regions"])
                        
                        # Add sequences
                        query_seq = ET.SubElement(hit_elem, "query_seq")
                        query_seq.text = ",".join(valid_hsps["query_seqs"])
                        
                        hit_seq = ET.SubElement(hit_elem, "hit_seq")
                        hit_seq.text = ",".join(valid_hsps["hit_seqs"])
        
        except Exception as e:
            self.logger.error(f"Error processing chain BLAST: {e}")

    def _stitch_hsps(self, hsps, domain_id=None, query_length=0):
        """
        Stitch together HSPs that represent parts of the same discontinuous domain.
        
        Args:
            hsps (list): List of HSP dictionaries with query/hit coordinates and scores
            domain_id (str): Domain ID for logging purposes
            query_length (int): Length of query sequence
            
        Returns:
            list: List of stitched HSPs representing complete domains
        """
        self.logger.debug(f"Stitching {len(hsps)} HSPs for domain {domain_id}")
        
        # If only one HSP, no stitching needed
        if len(hsps) <= 1:
            return hsps
            
        # Sort HSPs by query start position
        sorted_hsps = sorted(hsps, key=lambda x: x.get('query_from', 0))
        
        # Define maximum allowed gap between consecutive HSPs (adjust as needed)
        max_gap = 30  # residues
        
        # Initialize stitched HSPs list
        stitched_hsps = []
        current_group = [sorted_hsps[0]]
        
        # Iterate through remaining HSPs
        for i in range(1, len(sorted_hsps)):
            hsp = sorted_hsps[i]
            prev_hsp = current_group[-1]
            
            # Check if current HSP can be stitched to the current group
            query_gap = hsp.get('query_from', 0) - prev_hsp.get('query_to', 0)
            
            if query_gap <= max_gap:
                # Can be stitched - add to current group
                current_group.append(hsp)
            else:
                # Cannot be stitched - finalize current group and start new one
                stitched_hsp = self._merge_hsp_group(current_group, domain_id)
                stitched_hsps.append(stitched_hsp)
                current_group = [hsp]
        
        # Add the final group
        if current_group:
            stitched_hsp = self._merge_hsp_group(current_group, domain_id)
            stitched_hsps.append(stitched_hsp)
        
        self.logger.debug(f"Stitched {len(hsps)} HSPs into {len(stitched_hsps)} groups for domain {domain_id}")
        return stitched_hsps

    def _merge_hsp_group(self, hsp_group, domain_id=None):
        """
        Merge a group of HSPs into a single stitched HSP
        
        Args:
            hsp_group (list): List of HSP dictionaries to merge
            domain_id (str): Domain ID for reference
            
        Returns:
            dict: Merged HSP with combined attributes
        """
        if not hsp_group:
            return {}
        
        # If only one HSP in group, return as is
        if len(hsp_group) == 1:
            return hsp_group[0]
        
        # Create a new HSP that combines attributes from the group
        merged_hsp = {
            'domain_id': domain_id or hsp_group[0].get('domain_id', ''),
            'query_from': min(hsp.get('query_from', float('inf')) for hsp in hsp_group),
            'query_to': max(hsp.get('query_to', 0) for hsp in hsp_group),
            'hit_from': min(hsp.get('hit_from', float('inf')) for hsp in hsp_group),
            'hit_to': max(hsp.get('hit_to', 0) for hsp in hsp_group),
            'evalue': min(hsp.get('evalue', 999) for hsp in hsp_group),  # Take best e-value
            'discontinuous': True,  # Mark as discontinuous
            'segments': []  # Store individual segments
        }
        
        # Record information about each segment
        for hsp in hsp_group:
            segment = {
                'query_from': hsp.get('query_from', 0),
                'query_to': hsp.get('query_to', 0),
                'query_range': f"{hsp.get('query_from', 0)}-{hsp.get('query_to', 0)}",
                'hit_from': hsp.get('hit_from', 0),
                'hit_to': hsp.get('hit_to', 0),
                'hit_range': f"{hsp.get('hit_from', 0)}-{hsp.get('hit_to', 0)}",
                'evalue': hsp.get('evalue', 999)
            }
            merged_hsp['segments'].append(segment)
        
        # Add gap information for debugging/analysis
        gaps = []
        for i in range(len(hsp_group) - 1):
            current = hsp_group[i]
            next_hsp = hsp_group[i + 1]
            gap_size = next_hsp.get('query_from', 0) - current.get('query_to', 0)
            gaps.append(str(gap_size))
        
        merged_hsp['gap_info'] = ",".join(gaps)
        
        return merged_hsp

    def _process_blast(self, blast_path: str) -> List[BlastHit]:
        """Process domain BLAST results and return a list of BlastHit objects"""
        try:
            tree = ET.parse(blast_path)
            root = tree.getroot()
            
            # Get query info
            query_len_str = root.findtext(".//BlastOutput_query-len", "0")
            query_length = int(query_len_str) if query_len_str else 0
            
            # Group HSPs by domain for stitching
            domain_hsps = {}
            
            for iteration in root.findall(".//Iteration"):
                for hit in iteration.findall(".//Hit"):
                    hit_num = hit.findtext("Hit_num", "")
                    hit_def = hit.findtext("Hit_def", "")
                    hit_len = int(hit.findtext("Hit_len", "0"))
                    
                    # Parse domain ID, PDB ID, and chain from hit definition
                    hit_domain_id = "NA"
                    hit_pdb = "NA"
                    hit_chain = "NA"
                    
                    import re
                    # Pattern for the example: e8b7oAAA1 AAA:4-183 002982980
                    domain_match = re.search(r"((d|g|e)(\d\w{3})\w+\d*)\s+(\w+):", hit_def)
                    if domain_match:
                        hit_domain_id = domain_match.group(1)  # e8b7oAAA1
                        hit_pdb = domain_match.group(3)  # 8b7o
                        hit_chain = domain_match.group(4)  # AAA
                    
                    # Process HSPs for this hit
                    hit_hsps = []
                    
                    for hsp in hit.findall(".//Hsp"):
                        hsp_evalue = float(hsp.findtext("Hsp_evalue", "999"))
                        
                        # Get alignment coordinates
                        hsp_query_from = int(hsp.findtext("Hsp_query-from", "0"))
                        hsp_query_to = int(hsp.findtext("Hsp_query-to", "0"))
                        
                        hsp_hit_from = int(hsp.findtext("Hsp_hit-from", "0"))
                        hsp_hit_to = int(hsp.findtext("Hsp_hit-to", "0"))
                        
                        # Apply threshold filters
                        if hsp_evalue < self.hsp_evalue_threshold:
                            # Create HSP dictionary
                            hsp_dict = {
                                'domain_id': hit_domain_id,
                                'query_from': hsp_query_from,
                                'query_to': hsp_query_to,
                                'hit_from': hsp_hit_from,
                                'hit_to': hsp_hit_to,
                                'evalue': hsp_evalue
                            }
                            hit_hsps.append(hsp_dict)
                            
                            # Add to domain-specific group for stitching
                            if hit_domain_id not in domain_hsps:
                                domain_hsps[hit_domain_id] = []
                            domain_hsps[hit_domain_id].append(hsp_dict)
            
            # Perform HSP stitching for each domain
            stitched_domain_hsps = {}
            
            for domain_id, hsps in domain_hsps.items():
                if len(hsps) > 1:
                    self.logger.debug(f"Attempting to stitch {len(hsps)} HSPs for domain {domain_id}")
                    stitched_hsps = self._stitch_hsps(hsps, domain_id, query_length)
                    stitched_domain_hsps[domain_id] = stitched_hsps
                else:
                    stitched_domain_hsps[domain_id] = hsps
            
            # Create BlastHit objects from stitched HSPs
            blast_hits = []
            
            for domain_id, stitched_hsps in stitched_domain_hsps.items():
                for hsp in stitched_hsps:
                    # Extract domain components from ID
                    hit_match = re.search(r"([edg])(\d\w{3})(\w+)(\d+)", domain_id)
                    if hit_match:
                        hit_pdb = hit_match.group(2)
                        hit_chain = hit_match.group(3)
                    else:
                        hit_pdb = ""
                        hit_chain = ""
                    
                    # Format ranges correctly
                    if hsp.get('discontinuous', False) and 'segments' in hsp:
                        query_range = ",".join(seg['query_range'] for seg in hsp['segments'])
                        hit_range = ",".join(seg['hit_range'] for seg in hsp['segments'] if 'hit_range' in seg)
                    else:
                        query_range = f"{hsp.get('query_from', 0)}-{hsp.get('query_to', 0)}"
                        hit_range = f"{hsp.get('hit_from', 0)}-{hsp.get('hit_to', 0)}"
                    
                    # Create BlastHit object
                    blast_hit = BlastHit(
                        hit_id = str(len(blast_hits) + 1),
                        domain_id = domain_id,
                        pdb_id = hit_pdb,
                        chain_id = hit_chain,
                        evalue = hsp.get('evalue', 999.0),
                        hsp_count = len(hsp.get('segments', [])) if hsp.get('discontinuous', False) else 1,
                        hit_type = "domain_blast",
                        range = query_range,
                        hit_range = hit_range,
                        discontinuous = hsp.get('discontinuous', False)
                    )
                    
                    # Parse ranges
                    blast_hit.parse_ranges()
                    
                    blast_hits.append(blast_hit)
            
            return blast_hits
        
        except Exception as e:
            self.logger.error(f"Error processing domain BLAST: {e}")
            return []

    def _process_hhsearch_to_dict(self, hhsearch_path: str) -> List[HHSearchHit]:
        """Process HHSearch results and return a list of HHSearchHit objects"""
        try:
            tree = ET.parse(hhsearch_path)
            root = tree.getroot()
            
            # Verify root tag
            if root.tag != "hh_summ_doc":
                self.logger.error(f"Invalid XML format: expected root tag 'hh_summ_doc', found '{root.tag}'")
                return []
            
            # Find hits with the standard path
            hit_elements = root.findall(".//hh_hit_list/hh_hit")
            
            if not hit_elements:
                self.logger.warning(f"No hits found in HHSearch file: {hhsearch_path}")
                return []
            
            # Process each hit
            hits = []
            for hh_hit in hit_elements:
                # Extract attributes using the standard format
                domain_id = hh_hit.get("ecod_domain_id", "")
                hit_id = hh_hit.get("hit_id", "")
                probability = float(hh_hit.get("probability", "0"))
                e_value = float(hh_hit.get("e_value", "999"))
                score = float(hh_hit.get("score", "0"))
                
                # Create hit object
                hit = HHSearchHit(
                    hit_id=hit_id if hit_id else domain_id,
                    domain_id=domain_id,
                    probability=probability,
                    evalue=e_value,
                    score=score
                )
                
                # Process query range
                query_range_elem = hh_hit.find("query_range")
                if query_range_elem is not None:
                    range_value = query_range_elem.get("range", "")
                    hit.range = range_value
                    hit.parse_ranges()
                
                # Process template range
                template_range_elem = hh_hit.find("template_seqid_range")
                if template_range_elem is not None:
                    range_value = template_range_elem.get("range", "")
                    hit.hit_range = range_value
                
                hits.append(hit)
            
            return hits
        
        except Exception as e:
            self.logger.error(f"Error processing HHSearch results: {e}", exc_info=True)
            return []

    def _process_hit_hsps(self, hit_node: ET.Element, query_len: int) -> Optional[Dict[str, List[str]]]:
        """Process HSPs for a hit and return valid segments"""
        hit_len = int(hit_node.findtext("Hit_len", "0"))
        
        # Initialize tracking variables
        query_regions = []
        hit_regions = []
        query_seqs = []
        hit_seqs = []
        evalues = []
        
        # Track used and unused residues
        unused_hit_seqid = list(range(1, hit_len + 1))
        unused_query_seqid = list(range(1, query_len + 1))
        used_hit_seqid = []
        used_query_seqid = []
        
        hit_align_len = 0
        
        for hsp in hit_node.findall(".//Hsp"):
            hsp_evalue = float(hsp.findtext("Hsp_evalue", "999"))
            hsp_align_len = int(hsp.findtext("Hsp_align-len", "0"))
            
            hsp_query_from = int(hsp.findtext("Hsp_query-from", "0"))
            hsp_query_to = int(hsp.findtext("Hsp_query-to", "0"))
            
            hsp_hit_from = int(hsp.findtext("Hsp_hit-from", "0"))
            hsp_hit_to = int(hsp.findtext("Hsp_hit-to", "0"))
            
            # Get sequences
            hsp_qseq = hsp.findtext("Hsp_qseq", "")
            hsp_hseq = hsp.findtext("Hsp_hseq", "")
            
            # Check for overlap with previously used regions
            hsp_query_seqid = list(range(hsp_query_from, hsp_query_to + 1))
            hsp_hit_seqid = list(range(hsp_hit_from, hsp_hit_to + 1))
            
            used_query_residues = self._residue_coverage(hsp_query_seqid, used_query_seqid)
            used_hit_residues = self._residue_coverage(hsp_hit_seqid, used_hit_seqid)
            
            # Apply threshold filters
            if (hsp_evalue < self.hsp_evalue_threshold and
                used_hit_residues < 10 and used_query_residues < 5):
                
                # Add HSP to the hit
                query_regions.append(f"{hsp_query_from}-{hsp_query_to}")
                hit_regions.append(f"{hsp_hit_from}-{hsp_hit_to}")
                query_seqs.append(hsp_qseq)
                hit_seqs.append(hsp_hseq)
                evalues.append(str(hsp_evalue))
                
                hit_align_len += hsp_align_len
                
                # Update used/unused tracking
                for i in hsp_query_seqid:
                    if i not in used_query_seqid:
                        used_query_seqid.append(i)
                    if i in unused_query_seqid:
                        unused_query_seqid.remove(i)
                        
                for i in hsp_hit_seqid:
                    if i not in used_hit_seqid:
                        used_hit_seqid.append(i)
                    if i in unused_hit_seqid:
                        unused_hit_seqid.remove(i)
        
        # Check hit coverage
        hit_diff_tol = 50
        query_diff_tol = 50
        
        if (abs(hit_align_len - hit_len) < hit_diff_tol and 
            abs(query_len - hit_len) < query_diff_tol):
            
            # Return valid segments
            return {
                "query_regions": query_regions,
                "hit_regions": hit_regions,
                "query_seqs": query_seqs,
                "hit_seqs": hit_seqs,
                "evalues": evalues
            }
        
        return None
    
    def _residue_coverage(self, set1: List[int], set2: List[int]) -> int:
        """Count residues in set1 that are also in set2"""
        overlap = 0
        for residue in set1:
            if residue in set2:
                overlap += 1
        return overlap

    def _find_hhsearch_file(self, pdb_id: str, chain_id: str, reference: str, job_dump_dir: str) -> str:
        """
        Find HHSearch result file with strict standard pattern
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            reference: Reference version
            job_dump_dir: Base directory for output
            
        Returns:
            Path to HHSearch file or empty string if not found
        """
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Only accept the standard naming pattern
        standard_pattern = f"{pdb_chain}.{reference}.hhsearch.xml"
        standard_path = os.path.join(job_dump_dir, "hhsearch", standard_pattern)
        
        # Check if the standard pattern exists
        if os.path.exists(standard_path) and os.path.getsize(standard_path) > 0:
            self.logger.info(f"Found standard HHSearch XML: {standard_path}")
            return standard_path
        
        # If non-standard files exist, log warnings but don't use them
        non_standard_patterns = [
            f"{pdb_chain}.{reference}.hhr",
            f"{pdb_chain}.hhsearch.{reference}.xml",
            f"{pdb_chain}.{reference}.hh_summ.xml",
            f"{pdb_chain}.hhsearch.hhr",
            f"{pdb_chain}.hhr"
        ]
        
        search_locations = [
            os.path.join(job_dump_dir, "hhsearch"),
            os.path.join(job_dump_dir, pdb_chain),
            os.path.join(job_dump_dir, "ecod_dump", pdb_chain)
        ]
        
        # Check for non-standard files and log warnings
        for location in search_locations:
            for pattern in non_standard_patterns:
                path = os.path.join(location, pattern)
                if os.path.exists(path) and os.path.getsize(path) > 0:
                    self.logger.warning(f"Found non-standard HHSearch file: {path}")
                    self.logger.warning(f"Only accepting standard pattern: {standard_pattern}")
        
        # Return empty string if standard file not found
        self.logger.error(f"No standard HHSearch file found for {pdb_chain}")
        return ""

    def _parse_hhr(self, hhr_file: str) -> Dict[str, Any]:
        """
        Parse HHR file to extract structured data
        
        Args:
            hhr_file: Path to HHR file
            
        Returns:
            Dictionary with parsed data
        """
        try:
            with open(hhr_file, 'r') as f:
                content = f.read()
                
            # Parse header information
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
                    
            # Find the beginning of the hit table
            hit_table_start = None
            for i, line in enumerate(lines):
                if line.startswith(' No Hit'):
                    hit_table_start = i + 1
                    break
            
            if not hit_table_start:
                return {'header': header, 'hits': []}
            
            # Process hits
            hits = []
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
                
            return {
                'header': header,
                'hits': hits
            }
                
        except Exception as e:
            self.logger.error(f"Error parsing HHR file {hhr_file}: {e}")
            return {'header': {}, 'hits': []}

    def _calculate_range(self, alignment: str, start_pos: int) -> str:
        """
        Calculate range from alignment
        
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

    def _calculate_coverage(self, query_ali: str, template_ali: str) -> float:
        """
        Calculate coverage between query and template alignments
        
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

    def _process_self_comparison_to_dict(self, self_comp_path: str) -> List[Dict[str, Any]]:
        """Process self-comparison results and return a list of dictionaries"""
        try:
            tree = ET.parse(self_comp_path)
            root = tree.getroot()
            
            results = []
            
            # Process structural repeats
            for repeat in root.findall(".//structural_repeat"):
                aligner = repeat.get("aligner", "")
                zscore = float(repeat.get("zscore", "0"))
                
                ref_range = repeat.findtext("ref_range", "")
                mob_range = repeat.findtext("mob_range", "")
                
                results.append({
                    "type": "structural",
                    "aligner": aligner,
                    "z_score": zscore,
                    "query_range": ref_range,
                    "hit_range": mob_range
                })
            
            # Process sequence repeats
            for repeat_set in root.findall(".//sequence_repeat_set"):
                aligner = repeat_set.get("aligner", "")
                type_val = repeat_set.get("type", "")
                
                # Collect all ranges in this set
                set_ranges = []
                for seq_repeat in repeat_set.findall("sequence_repeat"):
                    range_val = seq_repeat.findtext("range", "")
                    set_ranges.append(range_val)
                
                results.append({
                    "type": "sequence",
                    "aligner": aligner,
                    "repeat_type": type_val,
                    "ranges": set_ranges
                })
            
            return results
            
        except Exception as e:
            self.logger.error(f"Error processing self comparison: {e}")
            return []

    def _find_self_comparison(self, pdb_id: str, chain_id: str, job_dump_dir: str) -> str:
        """Find self-comparison file for a protein chain
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            job_dump_dir: Base directory for job files
            
        Returns:
            Path to self-comparison file or empty string if not found
        """
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Check standard locations
        self_comp_path = os.path.join(job_dump_dir, "self_comparisons", f"{pdb_chain}.self_comp.xml")
        
        # If not in new location, check old
        if not os.path.exists(self_comp_path):
            alt_self_comp = os.path.join(job_dump_dir, pdb_chain, f"{pdb_chain}.self_comp.xml")
            if os.path.exists(alt_self_comp):
                self.logger.info(f"Found self-comparison at alternate location: {alt_self_comp}")
                return alt_self_comp
        else:
            self.logger.info(f"Found self-comparison file: {self_comp_path}")
            return self_comp_path
        
        self.logger.warning(f"No self-comparison file found for {pdb_id}_{chain_id}")
        return ""

    def _process_chain_blast_to_dict(self, chain_blast_path: str) -> List[BlastHit]:
        """Process chain-wise BLAST results and return a list of BlastHit objects"""
        try:
            tree = ET.parse(chain_blast_path)
            root = tree.getroot()
            
            # Get query info
            query_len_str = root.findtext(".//BlastOutput_query-len", "0")
            query_length = int(query_len_str) if query_len_str else 0
            
            # Process hits
            blast_hits = []
            
            for iteration in root.findall(".//Iteration"):
                for hit_elem in iteration.findall(".//Hit"):
                    hit_num = hit_elem.findtext("Hit_num", "")
                    hit_def = hit_elem.findtext("Hit_def", "")
                    hit_len = int(hit_elem.findtext("Hit_len", "0"))
                    
                    # Parse PDB ID and chain from hit definition
                    hit_pdb = "NA"
                    hit_chain = "NA"
                    if " " in hit_def:
                        parts = hit_def.split()
                        if len(parts) >= 2:
                            hit_pdb = parts[0]
                            hit_chain = parts[1]
                    
                    # Process HSPs
                    valid_hsps = self._process_hit_hsps(hit_elem, query_length)
                    
                    if valid_hsps:
                        # Create BlastHit object
                        blast_hit = BlastHit(
                            hit_id=hit_num,
                            pdb_id=hit_pdb,
                            chain_id=hit_chain,
                            hsp_count=len(valid_hsps["query_regions"]),
                            evalues=[float(e) for e in valid_hsps["evalues"]],
                            evalue=float(valid_hsps["evalues"][0]) if valid_hsps["evalues"] else 999.0,
                            hit_type="chain_blast",
                            range=",".join(valid_hsps["query_regions"]),
                            hit_range=",".join(valid_hsps["hit_regions"]),
                            query_seq=",".join(valid_hsps["query_seqs"]),
                            hit_seq=",".join(valid_hsps["hit_seqs"])
                        )
                        
                        # Parse ranges
                        blast_hit.parse_ranges()
                        
                        blast_hits.append(blast_hit)
            
            return blast_hits
        
        except Exception as e:
            self.logger.error(f"Error processing chain BLAST: {e}")
            return []
