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

from ecod.config import ConfigManager
from ecod.db.manager import DBManager
from ecod.exceptions import PipelineError, FileOperationError


class DomainSummary:
    """Process and integrate BLAST and HHSearch results for a protein chain"""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.summary")
        
        # Set default thresholds
        self.hsp_evalue_threshold = 0.005
        self.hit_coverage_threshold = 0.7

    def simplified_file_path_resolution(self, pdb_id, chain_id, file_type, job_dump_dir):
        """Simplified method to locate files, handling absolute paths"""
        db_config = self.config_manager.get_db_config()
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
                     job_dump_dir: str, blast_only: bool = False) -> str:
        """Create domain summary for a protein chain"""
        # Define paths and check for existing files
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Use the new directory structure
        domains_dir = os.path.join(job_dump_dir, "domains")
        os.makedirs(domains_dir, exist_ok=True)
        
        # ALWAYS use the standardized output path
        output_path = os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml")

        # Check for existing file
        if os.path.exists(output_path) and not self.config.get('force_overwrite', False):
            self.logger.warning(f"Output file {output_path} already exists, skipping...")
            return output_path

        # Determine FASTA file path - try multiple potential locations
        fasta_path = None

        # Try standard locations first with sub-batch structure
        potential_fasta_paths = [
            os.path.join(job_dump_dir, "fastas", f"{pdb_chain}.fa"),  # Main fastas dir
            os.path.join(job_dump_dir, "fastas", "batch_0", f"{pdb_chain}.fa"),  # Common sub-batch
            os.path.join(job_dump_dir, "fastas", "batch_1", f"{pdb_chain}.fa"),  # Alternative sub-batch
            os.path.join(job_dump_dir, "fastas", "batch_2", f"{pdb_chain}.fa"),  # Another alternative
            os.path.join(job_dump_dir, pdb_chain, f"{pdb_chain}.fa")  # Legacy location
        ]

        # Try each path
        for path in potential_fasta_paths:
            if os.path.exists(path):
                fasta_path = path
                self.logger.info(f"Found FASTA file at: {fasta_path}")
                break

        # If still not found, query the database as a last resort
        if not fasta_path:
            db_config = self.config_manager.get_db_config()
            db = DBManager(db_config)
            query = """
            SELECT pf.file_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            AND pf.file_type = 'fasta'
            AND pf.file_exists = TRUE
            ORDER BY pf.id DESC
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
        
        # Define output file in the new structure
        output_path = os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml")

        # Check for existing file
        if os.path.exists(output_path) and not self.config.get('force_overwrite', False):
            self.logger.warning(f"Output file {output_path} already exists, skipping...")
            return output_path
        
        # Check for peptides
        if sequence and len(sequence) < 30:  # Typical cutoff for peptides
            self.logger.warning(f"Sequence for {pdb_id}_{chain_id} is a peptide with length {len(sequence)}")
            
            # Create a special summary for peptides
            peptide_summary = ET.Element("blast_summ_doc")
            blast_summ = ET.SubElement(peptide_summary, "blast_summ")
            blast_summ.set("pdb", pdb_id)
            blast_summ.set("chain", chain_id)
            blast_summ.set("is_peptide", "true")
            blast_summ.set("sequence_length", str(len(sequence)))
            
            # Write output file to new location
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            tree = ET.ElementTree(peptide_summary)
            tree.write(output_path, encoding='utf-8', xml_declaration=True)
            
            self.logger.info(f"Created peptide summary: {output_path}")
            return output_path
        
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
            except Exception as e:
                self.logger.error(f"Error processing domain BLAST: {e}")
                blast_summ.set("domain_blast_error", "true")
        
        # Process HHSearch results (skip if blast_only mode)
        if not blast_only:
            hhsearch_path = os.path.join(job_dump_dir, "hhsearch", f"{pdb_chain}.{reference}.hh_summ.xml")
            
            # Try old location if not found
            if not os.path.exists(hhsearch_path):
                alt_hhsearch = os.path.join(job_dump_dir, pdb_chain, f"{pdb_chain}.{reference}.hh_summ.xml")
                if os.path.exists(alt_hhsearch):
                    hhsearch_path = alt_hhsearch
                
            if not os.path.exists(hhsearch_path):
                self.logger.warning(f"No hhsearch result file for {reference} {pdb_id} {chain_id}")
                blast_summ.set("no_hhsearch", "true")
            else:
                self._process_hhsearch(hhsearch_path, blast_summ)
        
        # Write output file to new structure only
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding='utf-8', xml_declaration=True)
        
        self.logger.info(f"Created domain summary: {output_path}")
        return output_path

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
    
    def _process_blast(self, blast_path: str, parent_node: ET.Element) -> None:
        """Process domain BLAST results"""
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
            
            # Process hits
            hits_node = ET.SubElement(blast_run, "hits")
            
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
                    domain_match = re.search(r"((d|g|e)(\d\w{3})\w+\d+)\s+(\w+):", hit_def)
                    if domain_match:
                        hit_domain_id = domain_match.group(1)
                        hit_pdb = domain_match.group(3)
                        hit_chain = domain_match.group(4)
                    
                    # Process HSPs for this hit
                    hit_align_len = 0
                    query_segs = []
                    hit_segs = []
                    evals = []
                    hsp_count = 0
                    
                    # Track used and unused residues
                    unused_hit_seqid = list(range(1, hit_len + 1))
                    unused_query_seqid = list(range(1, int(query_len) + 1))
                    used_hit_seqid = []
                    used_query_seqid = []
                    
                    for hsp in hit.findall(".//Hsp"):
                        hsp_evalue = float(hsp.findtext("Hsp_evalue", "999"))
                        hsp_align_len = int(hsp.findtext("Hsp_align-len", "0"))
                        
                        hsp_query_from = int(hsp.findtext("Hsp_query-from", "0"))
                        hsp_query_to = int(hsp.findtext("Hsp_query-to", "0"))
                        
                        hsp_hit_from = int(hsp.findtext("Hsp_hit-from", "0"))
                        hsp_hit_to = int(hsp.findtext("Hsp_hit-to", "0"))
                        
                        # Check for overlap with previously used regions
                        hsp_query_seqid = list(range(hsp_query_from, hsp_query_to + 1))
                        hsp_hit_seqid = list(range(hsp_hit_from, hsp_hit_to + 1))
                        
                        used_query_residues = self._residue_coverage(hsp_query_seqid, used_query_seqid)
                        used_hit_residues = self._residue_coverage(hsp_hit_seqid, used_hit_seqid)
                        
                        # Apply threshold filters
                        if (hsp_evalue < self.hsp_evalue_threshold and
                            used_hit_residues < 10 and used_query_residues < 5):
                            
                            # Add HSP to the hit
                            query_segs.append(f"{hsp_query_from}-{hsp_query_to}")
                            hit_segs.append(f"{hsp_hit_from}-{hsp_hit_to}")
                            evals.append(str(hsp_evalue))
                            hsp_count += 1
                            
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
                    
                    # Check coverage threshold
                    if (hit_align_len / hit_len) > self.hit_coverage_threshold:
                        # Sort regions by start position
                        def sort_key(region):
                            return int(region.split('-')[0])
                        
                        query_segs.sort(key=sort_key)
                        hit_segs.sort(key=sort_key)
                        
                        # Create hit element
                        hit_elem = ET.SubElement(hits_node, "hit")
                        hit_elem.set("num", hit_num)
                        hit_elem.set("domain_id", hit_domain_id)
                        hit_elem.set("pdb_id", hit_pdb)
                        hit_elem.set("chain_id", hit_chain)
                        hit_elem.set("hsp_count", str(hsp_count))
                        hit_elem.set("evalues", ",".join(evals))
                        
                        # Add query regions
                        query_reg = ET.SubElement(hit_elem, "query_reg")
                        query_reg.text = ",".join(query_segs)
                        
                        # Add hit regions
                        hit_reg = ET.SubElement(hit_elem, "hit_reg")
                        hit_reg.text = ",".join(hit_segs)
        
        except Exception as e:
            self.logger.error(f"Error processing domain BLAST: {e}")
    
    def _process_hhsearch(self, hhsearch_path: str, parent_node: ET.Element) -> None:
        """Process HHSearch results"""
        try:
            tree = ET.parse(hhsearch_path)
            root = tree.getroot()
            
            # Create HHSearch run section
            hh_run = ET.SubElement(parent_node, "hh_run")
            hh_run.set("program", "hhsearch")
            hh_run.set("db", "hora_full")
            
            hits_node = ET.SubElement(hh_run, "hits")
            
            hit_num = 1
            for hh_hit in root.findall(".//hh_hit"):
                # Skip obsolete structures
                if hh_hit.get("structure_obsolete", "") == "true":
                    continue
                
                # Create hit element
                hit_elem = ET.SubElement(hits_node, "hit")
                
                domain_id = hh_hit.get("ecod_domain_id", "")
                hit_elem.set("domain_id", domain_id)
                hit_elem.set("num", str(hit_num))
                
                prob = hh_hit.get("hh_prob", "")
                hit_elem.set("hh_prob", prob)
                
                score = hh_hit.get("hh_score", "")
                hit_elem.set("hh_score", score)
                
                # Add query and hit regions
                query_range = hh_hit.findtext("query_range", "")
                query_reg = ET.SubElement(hit_elem, "query_reg")
                query_reg.text = query_range
                
                hit_range = hh_hit.findtext("template_seqid_range", "")
                hit_reg = ET.SubElement(hit_elem, "hit_reg")
                hit_reg.text = hit_range
                
                hit_cover = hh_hit.find("template_seqid_range").get("ungapped_coverage", "")
                hit_elem.set("hit_cover", hit_cover)
                
                hit_num += 1
        
        except Exception as e:
            self.logger.error(f"Error processing HHSearch results: {e}")

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