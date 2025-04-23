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
                     job_dump_dir: str, blast_only: bool = False
    ) -> dict:
        """Create domain summary for a protein chain
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            reference: Reference version
            job_dump_dir: Base directory for output
            blast_only: Whether to use only BLAST results (no HHSearch)
            
        Returns:
            Dictionary with file path and statistics
        """
        # Define paths and check for existing files
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Initialize stats for return
        stats = {
            "chain_blast_processed": False,
            "domain_blast_processed": False,
            "hhsearch_processed": False,
            "hhsearch_hits": 0,
            "chain_blast_hits": 0,
            "domain_blast_hits": 0,
            "self_comp_processed": False
        }
        
        # Define output directory and filename once
        domains_dir = os.path.join(job_dump_dir, "domains")
        os.makedirs(domains_dir, exist_ok=True)
        
        # Define the output filename with proper context information
        suffix = ".blast_only" if blast_only else ""
        output_filename = f"{pdb_chain}.{reference}.domain_summary{suffix}.xml"
        output_path = os.path.join(domains_dir, output_filename)
        
        self.logger.info(f"Creating domain summary for {pdb_chain}, reference: {reference}, blast_only: {blast_only}")

        # Check for existing file
        if os.path.exists(output_path) and not self.context.is_force_overwrite():
            self.logger.warning(f"Output file {output_path} already exists, skipping...")
            return {"file_path": output_path, "stats": stats, "skipped": True}
        
        # Get protein sequence from FASTA file
        fasta_path = None
        
        # Try standard locations for FASTA file
        potential_fasta_paths = [
            os.path.join(job_dump_dir, "fastas", f"{pdb_chain}.fa"),
            os.path.join(job_dump_dir, "fastas", f"{pdb_chain}.fasta"),
            os.path.join(job_dump_dir, "fastas", "batch_0", f"{pdb_chain}.fa"),
            os.path.join(job_dump_dir, "fastas", "batch_0", f"{pdb_chain}.fasta"),
            os.path.join(job_dump_dir, "fastas", "batch_1", f"{pdb_chain}.fa"),
            os.path.join(job_dump_dir, "fastas", "batch_1", f"{pdb_chain}.fasta"),
        ]
        
        for path in potential_fasta_paths:
            if os.path.exists(path):
                fasta_path = path
                self.logger.info(f"Found FASTA file at: {fasta_path}")
                break
        
        # If not found, query database as last resort
        if not fasta_path:
            db_config = self.context.config_manager.get_db_config()
            db = DBManager(db_config)
            query = """
            SELECT pf.file_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            AND pf.file_type = 'fasta'
            AND pf.file_exists = TRUE
            LIMIT 1
            """
            try:
                rows = db.execute_query(query, (pdb_id, chain_id))
                if rows:
                    db_fasta_path = rows[0][0]
                    full_fasta_path = os.path.join(job_dump_dir, db_fasta_path)
                    if os.path.exists(full_fasta_path):
                        self.logger.info(f"Found FASTA file in database: {full_fasta_path}")
                        fasta_path = full_fasta_path
            except Exception as e:
                self.logger.warning(f"Error querying database for FASTA file: {e}")
        
        sequence = self._read_fasta_sequence(fasta_path)
        
        # Special handling for peptides
        if sequence and len(sequence) < 30:
            self.logger.warning(f"Sequence for {pdb_id}_{chain_id} is a peptide with length {len(sequence)}")
            
            # Create a special summary for peptides
            peptide_summary = ET.Element("blast_summ_doc")
            blast_summ = ET.SubElement(peptide_summary, "blast_summ")
            blast_summ.set("pdb", pdb_id)
            blast_summ.set("chain", chain_id)
            blast_summ.set("is_peptide", "true")
            blast_summ.set("sequence_length", str(len(sequence)))
            
            # Write to the defined output_path
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            tree = ET.ElementTree(peptide_summary)
            tree.write(output_path, encoding='utf-8', xml_declaration=True)
            
            self.logger.info(f"Created peptide summary: {output_path}")
            return {"file_path": output_path, "stats": stats, "is_peptide": True}
        
        # Create XML document root
        root = ET.Element("blast_summ_doc")
        
        # Create summary node
        blast_summ = ET.SubElement(root, "blast_summ")
        blast_summ.set("pdb", pdb_id)
        blast_summ.set("chain", chain_id)
        
        # Process self-comparison results
        self_comp_path = os.path.join(job_dump_dir, "self_comparisons", f"{pdb_chain}.self_comp.xml")
        
        # If not in new location, check old
        if not os.path.exists(self_comp_path):
            alt_self_comp = os.path.join(job_dump_dir, pdb_chain, f"{pdb_chain}.self_comp.xml")
            if os.path.exists(alt_self_comp):
                self_comp_path = alt_self_comp
        
        if not os.path.exists(self_comp_path):
            self.logger.warning(f"No self comparison results for {pdb_id} {chain_id}")
            blast_summ.set("no_selfcomp", "true")
        else:
            self._process_self_comparison(self_comp_path, blast_summ)
            stats["self_comp_processed"] = True
        
        # Find BLAST files using database
        chain_blast_paths = self.simplified_file_path_resolution(
            pdb_id, chain_id, 'chain_blast_result', job_dump_dir
        )

        if not chain_blast_paths:
            self.logger.error(f"No chain blast result file for {reference} {pdb_id} {chain_id}")
            blast_summ.set("no_chain_blast", "true")
        else:
            # Use the first file in the list
            chain_blast_file = chain_blast_paths[0]
            
            # Check if file has hits
            try:
                tree = ET.parse(chain_blast_file)
                root_elem = tree.getroot()
                hits = root_elem.findall(".//Hit")
                
                if not hits:
                    self.logger.warning(f"Chain BLAST file has no hits: {chain_blast_file}")
                    blast_summ.set("chain_blast_no_hits", "true")
                else:
                    self._process_chain_blast(chain_blast_file, blast_summ)
                    stats["chain_blast_processed"] = True
                    stats["chain_blast_hits"] = len(hits)
                    self.logger.info(f"Processed chain BLAST file with {len(hits)} hits")
            except Exception as e:
                self.logger.error(f"Error processing chain BLAST: {e}")
                blast_summ.set("chain_blast_error", "true")
        
        # Find domain BLAST files
        domain_blast_paths = self.simplified_file_path_resolution(
            pdb_id, chain_id, 'domain_blast_result', job_dump_dir
        )

        if not domain_blast_paths:
            self.logger.error(f"No domain blast result file for {reference} {pdb_id} {chain_id}")
            blast_summ.set("no_domain_blast", "true")
        else:
            # Use the first file in the list
            domain_blast_file = domain_blast_paths[0]
            
            # Check if file has hits
            try:
                tree = ET.parse(domain_blast_file)
                root_elem = tree.getroot()
                hits = root_elem.findall(".//Hit")
                
                if not hits:
                    self.logger.warning(f"Domain BLAST file has no hits: {domain_blast_file}")
                    blast_summ.set("domain_blast_no_hits", "true")
                else:
                    self._process_blast(domain_blast_file, blast_summ)
                    stats["domain_blast_processed"] = True
                    stats["domain_blast_hits"] = len(hits)
                    self.logger.info(f"Processed domain BLAST file with {len(hits)} hits")
            except Exception as e:
                self.logger.error(f"Error processing domain BLAST: {e}")
                blast_summ.set("domain_blast_error", "true")
        
        # Process HHSearch results (skip if blast_only mode)
        if not blast_only:
            # Try all potential HHSearch file patterns and locations
            hhsearch_found = False
            hhsearch_path = ""
            
            # Define all possible file patterns
            file_patterns = [
                f"{pdb_chain}.{reference}.hh_summ.xml",
                f"{pdb_chain}.{reference}.hhsearch.xml",
                f"{pdb_chain}.{reference}.hhr"
            ]
            
            # Define all possible locations
            locations = [
                os.path.join(job_dump_dir, "hhsearch"),
                os.path.join(job_dump_dir, pdb_chain),
                os.path.join(job_dump_dir, "ecod_dump", pdb_chain)
            ]
            
            # Try each combination
            for location in locations:
                for pattern in file_patterns:
                    path = os.path.join(location, pattern)
                    if os.path.exists(path) and os.path.getsize(path) > 0:
                        hhsearch_path = path
                        hhsearch_found = True
                        self.logger.info(f"Found HHSearch file: {path}")
                        break
                if hhsearch_found:
                    break
                    
            # If still not found, try database
            if not hhsearch_found:
                try:
                    db_config = self.context.config_manager.get_db_config()
                    db = DBManager(db_config)
                    
                    query = """
                    SELECT pf.file_path
                    FROM ecod_schema.process_file pf
                    JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
                    JOIN ecod_schema.protein p ON ps.protein_id = p.id
                    WHERE p.pdb_id = %s AND p.chain_id = %s
                    AND pf.file_type IN ('hhsearch_result', 'hhr', 'hhsearch_xml')
                    AND pf.file_exists = TRUE
                    LIMIT 1
                    """
                    
                    rows = db.execute_query(query, (pdb_id, chain_id))
                    if rows:
                        db_path = rows[0][0]
                        full_path = os.path.join(job_dump_dir, db_path) if not os.path.isabs(db_path) else db_path
                        if os.path.exists(full_path) and os.path.getsize(full_path) > 0:
                            hhsearch_path = full_path
                            hhsearch_found = True
                            self.logger.info(f"Found HHSearch file from database: {full_path}")
                except Exception as e:
                    self.logger.warning(f"Error querying database for HHSearch file: {e}")
            
            if not hhsearch_found:
                self.logger.warning(f"No hhsearch result file for {reference} {pdb_id} {chain_id}")
                blast_summ.set("no_hhsearch", "true")
            else:
                try:
                    # Process the file and count hits
                    if hhsearch_path.endswith('.xml'):
                        tree = ET.parse(hhsearch_path)
                        root_elem = tree.getroot()
                        hits = root_elem.findall(".//hh_hit")
                        hit_count = len(hits)
                        
                        self._process_hhsearch(hhsearch_path, blast_summ)
                        stats["hhsearch_processed"] = True
                        stats["hhsearch_hits"] = hit_count
                        self.logger.info(f"Processed HHSearch XML with {hit_count} hits")
                        
                    elif hhsearch_path.endswith('.hhr'):
                        # Parse HHR file to get hit count
                        hhr_data = self._parse_hhr(hhsearch_path)
                        hit_count = len(hhr_data.get('hits', []))
                        
                        self._process_hhr_directly(hhsearch_path, blast_summ)
                        stats["hhsearch_processed"] = True
                        stats["hhsearch_hits"] = hit_count
                        self.logger.info(f"Processed HHSearch HHR file with {hit_count} hits")
                except Exception as e:
                    self.logger.error(f"Error processing HHSearch file: {e}")
                    blast_summ.set("hhsearch_error", "true")
        
        # Write output file to new structure only
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        tree = ET.ElementTree(root)
        if os.path.exists(output_path) and self.context.is_force_overwrite():
            self.logger.info(f"Force overwrite enabled - regenerating {output_path}")

        tree.write(output_path, encoding='utf-8', xml_declaration=True)
        
        # Log comprehensive summary
        self.logger.info(f"Created domain summary: {output_path}")
        self.logger.info(f"Summary stats for {pdb_chain}: "
                       f"chain_hits={stats['chain_blast_hits']}, "
                       f"domain_hits={stats['domain_blast_hits']}, "
                       f"hhsearch_hits={stats['hhsearch_hits']}")
        
        return {"file_path": output_path, "stats": stats}

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

    def _process_blast(self, blast_path: str, parent_node: ET.Element) -> None:
        """Process domain BLAST results with HSP stitching for discontinuous domains"""
        try:
            tree = ET.parse(blast_path)
            root = tree.getroot()
            
            # Create BLAST run section
            blast_run = ET.SubElement(parent_node, "blast_run")
            
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
            
            # Process hits with HSP stitching
            hits_node = ET.SubElement(blast_run, "hits")
            
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
            query_length = int(query_len) if query_len else 0
            stitched_domain_hsps = {}
            
            for domain_id, hsps in domain_hsps.items():
                if len(hsps) > 1:
                    self.logger.debug(f"Attempting to stitch {len(hsps)} HSPs for domain {domain_id}")
                    stitched_hsps = self._stitch_hsps(hsps, domain_id, query_length)
                    stitched_domain_hsps[domain_id] = stitched_hsps
                else:
                    stitched_domain_hsps[domain_id] = hsps
            
            # Create hit elements with stitched HSPs
            for domain_id, stitched_hsps in stitched_domain_hsps.items():
                for hsp in stitched_hsps:
                    # Extract domain components from ID
                    hit_match = re.search(r"([edg])(\d\w{3})(\w+)(\d+)", domain_id)
                    if hit_match:
                        hit_pdb = hit_match.group(2)
                        hit_chain = hit_match.group(3)
                    
                    # Create hit element
                    hit_elem = ET.SubElement(hits_node, "hit")
                    hit_elem.set("domain_id", domain_id)
                    hit_elem.set("pdb_id", hit_pdb)
                    hit_elem.set("chain_id", hit_chain)
                    
                    # Format ranges correctly
                    if hsp.get('discontinuous', False) and 'segments' in hsp:
                        query_range = ",".join(seg['query_range'] for seg in hsp['segments'])
                        hit_range = ",".join(seg['hit_range'] for seg in hsp['segments'] if seg['hit_range'])
                    else:
                        query_range = f"{hsp.get('query_from', 0)}-{hsp.get('query_to', 0)}"
                        hit_range = f"{hsp.get('hit_from', 0)}-{hsp.get('hit_to', 0)}"
                    
                    # Add query regions
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = query_range
                    
                    # Add hit regions
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit_range
                    
                    # Add HSP count and evalue
                    hit_elem.set("hsp_count", "1" if not hsp.get('discontinuous', False) else str(len(hsp.get('segments', []))))
                    hit_elem.set("evalues", str(hsp.get('evalue', 999)))
                    
                    # Add discontinuous flag if needed
                    if hsp.get('discontinuous', False):
                        hit_elem.set("discontinuous", "true")
                        if 'gap_info' in hsp:
                            hit_elem.set("gap_info", hsp['gap_info'])
        
        except Exception as e:
            self.logger.error(f"Error processing domain BLAST: {e}")

    def _process_hhsearch(self, hhsearch_path: str, parent_node: ET.Element) -> None:
        """
        Process HHSearch results (compatible with multiple XML formats)
        
        Args:
            hhsearch_path: Path to HHSearch result file
            parent_node: Parent XML node to add HHSearch results to
        """
        try:
            tree = ET.parse(hhsearch_path)
            root = tree.getroot()
            
            # Create HHSearch run section
            hh_run = ET.SubElement(parent_node, "hh_run")
            hh_run.set("program", "hhsearch")
            
            hits_node = ET.SubElement(hh_run, "hits")
            
            # Determine XML format based on root tag
            is_new_format = root.tag == "hh_summ_doc"
            self.logger.debug(f"Processing HHSearch file: {hhsearch_path} (new format: {is_new_format})")
            
            # Adjust XPath based on document structure
            if is_new_format:
                # Try the newer structure with hh_hit_list container first
                hit_xpath = ".//hh_hit_list/hh_hit"
                hit_elements = root.findall(hit_xpath)
                
                # If no hits found, try the simpler direct path
                if not hit_elements:
                    hit_xpath = ".//hh_hit"
                    hit_elements = root.findall(hit_xpath)
            else:
                hit_xpath = ".//hh_hit"
                hit_elements = root.findall(hit_xpath)
            
            hit_num = 1
            for hh_hit in hit_elements:
                # Skip obsolete structures
                if hh_hit.get("structure_obsolete", "") == "true":
                    continue
                
                # Create hit element
                hit_elem = ET.SubElement(hits_node, "hit")
                
                # Extract attributes using appropriate keys based on format
                if is_new_format:
                    domain_id = hh_hit.get("ecod_domain_id", "")
                    prob = hh_hit.get("probability", "")
                    score = hh_hit.get("score", "")
                    evalue = hh_hit.get("e_value", "")
                else:
                    domain_id = hh_hit.get("ecod_domain_id", "")
                    prob = hh_hit.get("hh_prob", "")
                    score = hh_hit.get("hh_score", "")
                    evalue = hh_hit.get("hh_evalue", "")
                
                # Set attributes on output element
                hit_elem.set("domain_id", domain_id)
                hit_elem.set("num", str(hit_num))
                
                # Store both new and old naming for compatibility
                hit_elem.set("probability", prob)
                hit_elem.set("hh_prob", prob)
                
                hit_elem.set("score", score)
                hit_elem.set("hh_score", score)
                
                if evalue:
                    hit_elem.set("evalue", evalue)
                
                # Extract and add query region - handle both formats
                query_range_elem = hh_hit.find("query_range")
                if query_range_elem is not None:
                    # New format - get range from attribute
                    query_range = query_range_elem.get("range", "")
                    if not query_range:
                        # Fallback to constructing from start/end
                        start = query_range_elem.get("start", "")
                        end = query_range_elem.get("end", "")
                        if start and end:
                            query_range = f"{start}-{end}"
                else:
                    # Old format - get text content
                    query_range = hh_hit.findtext("query_reg", "")
                
                # Only add if we have range information
                if query_range:
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = query_range
                    self.logger.debug(f"Added query range: {query_range}")
                else:
                    self.logger.warning(f"No query range found for hit {hit_num}")
                
                # Extract and add hit region - handle both formats
                hit_range_elem = hh_hit.find("template_seqid_range")
                if hit_range_elem is not None:
                    # New format - get range from attribute
                    hit_range = hit_range_elem.get("range", "")
                    if not hit_range:
                        # Fallback to constructing from start/end
                        start = hit_range_elem.get("start", "")
                        end = hit_range_elem.get("end", "")
                        if start and end:
                            hit_range = f"{start}-{end}"
                    
                    # Get coverage information
                    hit_cover = hit_range_elem.get("identity", "") or hit_range_elem.get("ungapped_coverage", "")
                else:
                    # Old format - get text content
                    hit_range = hh_hit.findtext("hit_reg", "") or hh_hit.findtext("template_reg", "")
                    hit_cover = ""
                
                # Only add if we have range information
                if hit_range:
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = hit_range
                    self.logger.debug(f"Added hit range: {hit_range}")
                    
                    if hit_cover:
                        hit_elem.set("hit_cover", hit_cover)
                else:
                    self.logger.warning(f"No hit range found for hit {hit_num}")
                
                hit_num += 1
                
            self.logger.info(f"Processed {hit_num-1} HHSearch hits")
        
        except Exception as e:
            self.logger.error(f"Error processing HHSearch results: {e}", exc_info=True)

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

    def _find_hhsearch_file(self, pdb_id: str, chain_id: str, reference: str, 
                           job_dump_dir: str
    ) -> str:
        """
        Find HHSearch result file with better fallback options
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            reference: Reference version
            job_dump_dir: Base directory for output
            
        Returns:
            Path to HHSearch file or empty string if not found
        """
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Try standard locations
        potential_paths = [
            # New standard location (flat structure in hhsearch dir)
            os.path.join(job_dump_dir, "hhsearch", f"{pdb_chain}.{reference}.hhr"),
            os.path.join(job_dump_dir, "hhsearch", f"{pdb_chain}.{reference}.hhsearch.xml"),
            
            # Old chain-specific directory structure
            os.path.join(job_dump_dir, pdb_chain, f"{pdb_chain}.{reference}.hhr"),
            os.path.join(job_dump_dir, pdb_chain, f"{pdb_chain}.{reference}.hhsearch.xml"),
            
            # Other potential locations
            os.path.join(job_dump_dir, "ecod_dump", pdb_chain, f"{pdb_chain}.{reference}.hhr"),
            os.path.join(job_dump_dir, "ecod_dump", pdb_chain, f"{pdb_chain}.{reference}.hhsearch.xml")
        ]
        
        for path in potential_paths:
            if os.path.exists(path) and os.path.getsize(path) > 0:
                self.logger.info(f"Found HHSearch file at: {path}")
                return path
        
        # If not found in standard locations, query database
        try:
            db_config = self.context.config_manager.get_db_config()
            db = DBManager(db_config)
            
            # Try to find process ID first
            process_query = """
            SELECT ps.id
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            JOIN ecod_schema.batch b ON ps.batch_id = b.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            ORDER BY ps.id DESC
            LIMIT 1
            """
            
            process_rows = db.execute_query(process_query, (pdb_id, chain_id))
            if not process_rows:
                self.logger.warning(f"No process found for {pdb_id}_{chain_id}")
                return ""
            
            process_id = process_rows[0][0]
            
            # Look for HHSearch files
            file_query = """
            SELECT file_path
            FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type IN ('hhr', 'hhsearch_xml')
            ORDER BY file_type
            """
            
            file_rows = db.execute_query(file_query, (process_id,))
            if file_rows:
                file_path = os.path.join(job_dump_dir, file_rows[0][0])
                if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                    self.logger.info(f"Found HHSearch file from database: {file_path}")
                    return file_path
        except Exception as e:
            self.logger.error(f"Error querying database for HHSearch file: {e}")
        
        self.logger.warning(f"No HHSearch file found for {pdb_id}_{chain_id}")
        return ""

    def _process_hhr_directly(self, hhr_file: str, parent_node: ET.Element) -> None:
        """
        Process HHR file directly if XML is not available
        
        Args:
            hhr_file: Path to HHR file
            parent_node: Parent XML node to add HHSearch results to
        """
        try:
            # Parse HHR file
            hhr_data = self._parse_hhr(hhr_file)
            if not hhr_data:
                self.logger.error(f"Failed to parse HHR file: {hhr_file}")
                return
            
            # Create HHSearch run section
            hh_run = ET.SubElement(parent_node, "hh_run")
            hh_run.set("program", "hhsearch")
            
            hits_node = ET.SubElement(hh_run, "hits")
            
            # Process hits
            for hit_num, hit in enumerate(hhr_data.get('hits', []), 1):
                # Create hit element
                hit_elem = ET.SubElement(hits_node, "hit")
                
                # Add basic attributes
                hit_elem.set("num", str(hit_num))
                hit_id = hit.get('hit_id', '')
                hit_elem.set("hit_id", hit_id)
                
                # Extract domain ID if it matches pattern
                if re.match(r'[dge]\d\w{3}\w+\d+', hit_id):
                    hit_elem.set("domain_id", hit_id)
                
                # Add probability
                prob = hit.get('probability')
                if prob is not None:
                    hit_elem.set("probability", str(prob))
                    hit_elem.set("hh_prob", str(prob))  # For compatibility
                
                # Add e-value
                evalue = hit.get('e_value')
                if evalue is not None:
                    hit_elem.set("evalue", str(evalue))
                
                # Add score
                score = hit.get('score')
                if score is not None:
                    hit_elem.set("score", str(score))
                    hit_elem.set("hh_score", str(score))  # For compatibility
                
                # Process alignment and extract query range
                if 'query_ali' in hit and 'query_start' in hit:
                    query_range = self._calculate_range(hit['query_ali'], hit['query_start'])
                    query_reg = ET.SubElement(hit_elem, "query_reg")
                    query_reg.text = query_range
                
                # Process template range
                if 'template_ali' in hit and 'template_start' in hit:
                    template_range = self._calculate_range(hit['template_ali'], hit['template_start'])
                    hit_reg = ET.SubElement(hit_elem, "hit_reg")
                    hit_reg.text = template_range
                    
                    # Calculate coverage
                    coverage = self._calculate_coverage(hit['query_ali'], hit['template_ali'])
                    hit_elem.set("hit_cover", str(coverage))
            
            self.logger.info(f"Processed {len(hhr_data.get('hits', []))} hits from HHR file")
            
        except Exception as e:
            self.logger.error(f"Error processing HHR file: {e}")

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
