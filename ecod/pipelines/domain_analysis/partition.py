#!/usr/bin/env python3
"""
Domain partition module for the ECOD pipeline
Determines protein domain boundaries and classifications
"""

import os
import re
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Set

from ecod.config import ConfigManager
from ecod.db.manager import DBManager
from ecod.exceptions import PipelineError, FileOperationError
from ecod.core.context import ApplicationContext


class DomainPartition:
    """Determine domain boundaries and classifications from search results"""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.partition")
        
        # Set default thresholds
        self.new_coverage_threshold = 0.9
        self.old_coverage_threshold = 0.05
        self.dali_significance_threshold = 5.0
        self.hit_coverage_threshold = 0.7
        self.gap_tol = 20
        
        # Load reference data
        self.ref_range_cache = {}
        self.ref_domain_uid_lookup = {}
        self.ref_chain_domains = {}

        # Initialize classification caches
        self.domain_classification_cache = {}
        self.domain_id_classification_cache = {}

    def process_specific_ids(self, batch_id: int, process_ids: List[int], 
                        dump_dir: str, reference: str, blast_only: bool = False) -> bool:
        """Process domain partition for specific process IDs
        
        Args:
            batch_id: Batch ID
            process_ids: List of process IDs to process
            dump_dir: Base directory for output
            reference: Reference version
            blast_only: Whether to use only blast summaries (No HHsearch)
            
        Returns:
            True if all specified processes were processed successfully
        """
        # Get database connection
        db_config = self.config_manager.get_db_config()
        db = DBManager(db_config)
        
        # Get specific protein details by process IDs
        query = """
        SELECT 
            ps.id, p.pdb_id, p.chain_id, ps.relative_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            ps.id IN %s
            AND ps.batch_id = %s
        """
        
        try:
            rows = db.execute_dict_query(query, (tuple(process_ids), batch_id))
        except Exception as e:
            self.logger.error(f"Error querying protein details: {e}")
            return False
        
        if not rows:
            self.logger.warning(f"No proteins found with specified process IDs in batch {batch_id}")
            return False
        
        # Process each protein
        success_count = 0
        
        for row in rows:
            pdb_id = row["pdb_id"]
            chain_id = row["chain_id"]
            process_id = row["id"]
            
            try:
                # Get domain summary path
                summary_query = """
                SELECT file_path 
                FROM ecod_schema.process_file
                WHERE process_id = %s AND file_type = 'domain_summary'
                """
                summary_result = db.execute_query(summary_query, (process_id,))
                
                if not summary_result:
                    self.logger.warning(f"No domain summary found for {pdb_id}_{chain_id} (process_id: {process_id})")
                    continue
                    
                # Verify summary exists in filesystem
                summary_path = os.path.join(dump_dir, summary_result[0][0])
                if not os.path.exists(summary_path):
                    self.logger.warning(f"Domain summary file not found: {summary_path}")
                    continue
                
                # Run partition for this protein
                self.logger.info(f"Processing domains for {pdb_id}_{chain_id} (process_id: {process_id})")
                
                domain_file = self.partition_domains(
                    pdb_id,
                    chain_id,
                    dump_dir,
                    'struct_seqid',  # Default input mode
                    reference,
                    blast_only
                )
                
                if domain_file:
                    success_count += 1
                    
                    # Update process status
                    db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "domain_partition_complete",
                            "status": "success"
                        },
                        "id = %s",
                        (process_id,)
                    )
                    
                    # Register domain file
                    db.insert(
                        "ecod_schema.process_file",
                        {
                            "process_id": process_id,
                            "file_type": "domain_partition",
                            "file_path": os.path.relpath(domain_file, dump_dir),
                            "file_exists": True,
                            "file_size": os.path.getsize(domain_file) if os.path.exists(domain_file) else 0
                        }
                    )
                    
            except Exception as e:
                self.logger.error(f"Error processing domains for {pdb_id}_{chain_id}: {e}")
                
                # Update process status
                db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": "domain_partition_failed",
                        "status": "error",
                        "error_message": str(e)
                    },
                    "id = %s",
                    (process_id,)
                )
        
        self.logger.info(f"Processed domains for {success_count}/{len(rows)} proteins from specified process IDs")
        return success_count > 0

    def process_batch(self, batch_id: int, dump_dir: str, reference: str, blast_only: bool = False, limit: int = None) -> List[str]:
        """Process domain partition for a batch of proteins
        
        Args:
            batch_id: Batch ID
            dump_dir: Base directory for output
            reference: Reference version
            blast_only: Whether to use only blast summaries (No HHsearch)
            limit: Maximum number of proteins to process
            
        Returns:
            List of generated domain files
        """
        # Get database connection
        db_config = self.config_manager.get_db_config()
        db = DBManager(db_config)
        
        # Get proteins from the batch
        query = """
        SELECT 
            ps.id, p.pdb_id, p.chain_id, ps.relative_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 
            ps.batch_id = %s
            AND ps.status IN ('success', 'processing')
            AND pf.file_type = 'domain_summary'
            AND pf.file_exists = TRUE
        """
        
        if limit:
            query += f" LIMIT {limit}"
            
        try:
            rows = db.execute_dict_query(query, (batch_id,))
        except Exception as e:
            self.logger.error(f"Error querying batch proteins: {e}")
            return []
        
        if not rows:
            self.logger.warning(f"No proteins found for domain analysis in batch {batch_id}")
            return []
        
        # Process each protein
        domain_files = []
        
        for row in rows:
            pdb_id = row["pdb_id"]
            chain_id = row["chain_id"]
            
            try:
                domain_file = self.partition_domains(
                    pdb_id,
                    chain_id,
                    dump_dir,
                    'struct_seqid',  # Default input mode
                    reference,
                    blast_only
                )
                
                if domain_file:
                    domain_files.append(domain_file)
                    
                    # Update process status
                    db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "domain_partition_complete",
                            "status": "success"
                        },
                        "id = %s",
                        (row["id"],)
                    )
                    
                    # Register domain file
                    db.insert(
                        "ecod_schema.process_file",
                        {
                            "process_id": row["id"],
                            "file_type": "domain_partition",
                            "file_path": os.path.relpath(domain_file, dump_dir),
                            "file_exists": True,
                            "file_size": os.path.getsize(domain_file) if os.path.exists(domain_file) else 0
                        }
                    )
                    
            except Exception as e:
                self.logger.error(f"Error processing domains for {pdb_id}_{chain_id}: {e}")
                
                # Update process status
                db.update(
                    "ecod_schema.process_status",
                    {
                        "current_stage": "domain_partition_failed",
                        "status": "error",
                        "error_message": str(e)
                    },
                    "id = %s",
                    (row["id"],)
                )
        
        self.logger.info(f"Processed domains for {len(domain_files)} proteins from batch {batch_id}")
        return domain_files
        
    def _get_domain_classification(self, ecod_uid: int) -> Optional[Dict[str, Any]]:
        """Get domain classification from database with caching"""
        # Check cache first
        if ecod_uid in self.domain_classification_cache:
            return self.domain_classification_cache[ecod_uid]
            
        # Query database
        db_config = self.config_manager.get_db_config()
        db = DBManager(db_config)
        
        query = """
        SELECT 
            d.t_group, d.h_group, d.x_group, d.a_group,
            d.is_manual_rep, d.is_f70, d.is_f40, d.is_f99
        FROM 
            pdb_analysis.domain d
        WHERE 
            d.ecod_uid = %s
        """
        
        try:
            rows = db.execute_dict_query(query, (ecod_uid,))
            if rows:
                classification = {
                    "t_group": rows[0].get("t_group"),
                    "h_group": rows[0].get("h_group"),
                    "x_group": rows[0].get("x_group"),
                    "a_group": rows[0].get("a_group"),
                    "is_manual_rep": rows[0].get("is_manual_rep", False),
                    "is_f70": rows[0].get("is_f70", False),
                    "is_f40": rows[0].get("is_f40", False),
                    "is_f99": rows[0].get("is_f99", False)
                }
                
                # Cache the result
                self.domain_classification_cache[ecod_uid] = classification
                return classification
        except Exception as e:
            self.logger.error(f"Error getting domain classification for {ecod_uid}: {e}")
            
        return None

    def _get_domain_classification_by_id(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """Get domain classification from database by domain ID with caching"""
        # Check cache first
        if domain_id in self.domain_id_classification_cache:
            return self.domain_id_classification_cache[domain_id]
            
        # Query database
        db_config = self.config_manager.get_db_config()
        db = DBManager(db_config)
        
        query = """
        SELECT 
            d.t_group, d.h_group, d.x_group, d.a_group,
            d.is_manual_rep, d.is_f70, d.is_f40, d.is_f99
        FROM 
            pdb_analysis.domain d
        WHERE 
            d.domain_id = %s
        """
        
        try:
            rows = db.execute_dict_query(query, (domain_id,))
            if rows:
                classification = {
                    "t_group": rows[0].get("t_group"),
                    "h_group": rows[0].get("h_group"),
                    "x_group": rows[0].get("x_group"),
                    "a_group": rows[0].get("a_group"),
                    "is_manual_rep": False,  # New domains are not manual reps
                    "is_f70": False,
                    "is_f40": False,
                    "is_f99": False
                }
                
                # Cache the result
                self.domain_id_classification_cache[domain_id] = classification
                return classification
        except Exception as e:
            self.logger.error(f"Error getting classification for {domain_id}: {e}")
            
        return None
        
    def load_reference_data(self, reference: str) -> None:
        """Load reference domain classifications"""
        # In a real implementation, this would load data from a database
        # For this example, we'll simulate loading from a database
        self.logger.info(f"Loading reference data for {reference}")
        
        # Get database connection from config
        db_config = self.config_manager.get_db_config()
        db = DBManager(db_config)
        
        # Query to get reference domain ranges
        query = """
        SELECT 
            d.ecod_uid, d.domain_id, d.range, p.source_id
        FROM 
            pdb_analysis.domain d
        JOIN
            pdb_analysis.protein p ON d.protein_id = p.id
        """
        
        # Execute query
        try:
            result = db.execute_dict_query(query)
            
            # Process results
            for domain in result:
                uid = domain["ecod_uid"]
                domain_id = domain["domain_id"]
                range_str = domain["range"]
                source_id = domain["source_id"]
                
                # Store in reference cache
                if source_id not in self.ref_range_cache:
                    self.ref_range_cache[source_id] = []
                self.ref_range_cache[source_id].append({
                    "uid": uid,
                    "domain_id": domain_id,
                    "range": range_str
                })
                
                # Store in uid lookup
                self.ref_domain_uid_lookup[domain_id] = uid
            
            # Transform to chain-wise structure
            self._transform_reference_chain_wise()
                
            self.logger.info(f"Loaded {len(self.ref_range_cache)} chains with domains")
        except Exception as e:
            self.logger.error(f"Error loading reference data: {e}")
    
    def _transform_reference_chain_wise(self) -> None:
        """Transform reference data to chain-wise format"""
        for source_id, domains in self.ref_range_cache.items():
            if source_id not in self.ref_chain_domains:
                self.ref_chain_domains[source_id] = []
            
            # Sort domains by position
            domains_sorted = sorted(domains, key=lambda d: self._get_start_position(d["range"]))
            
            self.ref_chain_domains[source_id] = domains_sorted
    
    def _get_start_position(self, range_str: str) -> int:
        """Get the start position from a range string"""
        # Handle chain:position format (e.g., "A:1")
        if ":" in range_str:
            parts = range_str.split(":")
            # Extract just the position number
            range_str = parts[1]
        
        if "-" in range_str:
            parts = range_str.split("-")
            try:
                return int(parts[0])
            except ValueError:
                return 0
        elif "," in range_str:
            # Multi-segment range
            first_segment = range_str.split(",")[0]
            return self._get_start_position(first_segment)
        else:
            try:
                return int(range_str)
            except ValueError:
                return 0
    
    def partition_domains(self, pdb_id: str, chain_id: str, dump_dir: str, input_mode: str, reference: str, blast_only: bool = False) -> str:
        """Partition domains for a single protein chain"""
        # Load reference data if not already loaded
        if not self.ref_range_cache:
            self.load_reference_data(reference)
        
        # Define paths
        pdb_chain = f"{pdb_id}_{chain_id}"
        chain_dir = os.path.join(dump_dir, pdb_chain)  # Original chain dir
        domain_prefix = "domains_v12"  # Consistent with Perl script
        domain_fn = os.path.join(chain_dir, f"{domain_prefix}.{pdb_chain}.{reference}.xml")
        
        if os.path.exists(domain_fn) and not self.config.get('force_overwrite', False):
            self.logger.warning(f"Domain file {domain_fn} already exists, skipping...")
            return domain_fn
        
        # Define the domains directory
        domains_dir = os.path.join(dump_dir, "domains")
        
        # Look for domain summary file in the new standard location
        blast_summ_fn = os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml")
        
        # If not in standard location, try the legacy location
        if not os.path.exists(blast_summ_fn):
            alt_blast_summ = os.path.join(chain_dir, f"{pdb_chain}.{reference}.blast_summ{''.join(['.blast_only' if blast_only else ''])}.xml")
            if os.path.exists(alt_blast_summ):
                blast_summ_fn = alt_blast_summ
                self.logger.info(f"Found domain summary in legacy location: {blast_summ_fn}")
        
        if not os.path.exists(blast_summ_fn):
            self.logger.error(f"Blast summary file not found: {blast_summ_fn}")
            return None
        
        # Create domain document
        domains_doc = ET.Element("domain_doc")
        domains_doc.set("pdb", self._safe_str(pdb_id))
        domains_doc.set("chain", self._safe_str(chain_id))
        domains_doc.set("reference", self._safe_str(reference))
        
        # Process chain sequence - try multiple locations
        fasta_path = None
        potential_fasta_paths = [
            # New structure with batch subdirectories
            os.path.join(dump_dir, "fastas", f"{pdb_chain}.fa"),
            os.path.join(dump_dir, "fastas", "batch_0", f"{pdb_chain}.fasta"),
            os.path.join(dump_dir, "fastas", "batch_1", f"{pdb_chain}.fasta"),
            os.path.join(dump_dir, "fastas", "batch_2", f"{pdb_chain}.fasta"),
            os.path.join(dump_dir, "fastas", "batch_3", f"{pdb_chain}.fasta"),
            # Legacy location
            os.path.join(chain_dir, f"{pdb_chain}.fa")
        ]
        
        # Try each potential path
        for path in potential_fasta_paths:
            if os.path.exists(path):
                fasta_path = path
                self.logger.info(f"Found FASTA file at: {fasta_path}")
                break
        
        sequence = self._read_fasta_sequence(fasta_path)
        if not sequence:
            self.logger.error(f"Failed to read sequence from {fasta_path}")
            return None
        
        sequence_length = len(sequence)
        
        # Process blast summary
        blast_data = self._process_blast_summary(blast_summ_fn)
        
        # Determine domain boundaries
        domains = self._determine_domain_boundaries(blast_data, sequence_length, pdb_chain)
        
        # Assign classifications
        self._assign_domain_classifications(domains, blast_data, pdb_chain)
        
        # Create domain elements with enhanced debugging
        domain_list = ET.SubElement(domains_doc, "domain_list")

        for i, d in enumerate(domains):
            self.logger.debug(f"Processing domain {i+1} for XML: {d}")
            
            # Check for None values that would cause serialization issues
            for key, value in d.items():
                if value is None:
                    self.logger.warning(f"Domain {i+1} has None value for {key} - replacing with empty string")
                    d[key] = ""
            
            try:
                domain_elem = ET.SubElement(domain_list, "domain")
                domain_elem.set("pdb", str(pdb_id))
                domain_elem.set("chain", str(chain_id))
                
                # Check for required range attribute
                if "range" not in d or not d["range"]:
                    self.logger.error(f"Domain {i+1} missing required range attribute")
                    continue
                    
                domain_elem.set("range", str(d["range"]))
                
                # Add classification attributes if present
                for attr in ["t_group", "h_group", "x_group", "a_group"]:
                    if attr in d and d[attr]:
                        domain_elem.set(attr, str(d[attr]))
                
                # Add flags if present
                for flag in ["is_manual_rep", "is_f70", "is_f40", "is_f99"]:
                    if flag in d and d[flag]:
                        domain_elem.set(flag, "true")
                
                # Add range element
                range_elem = ET.SubElement(domain_elem, "range")
                range_elem.text = str(d["range"])
                
                # Add classification evidence
                if "evidence" in d and d["evidence"]:
                    evidence = ET.SubElement(domain_elem, "evidence")
                    for e_idx, e in enumerate(d["evidence"]):
                        if e is None:
                            self.logger.warning(f"Domain {i+1} has None evidence item at index {e_idx}")
                            continue
                        
                        self.logger.debug(f"Processing evidence item {e_idx+1}: {e}")
                        
                        # Check for None values in evidence
                        for e_key, e_value in e.items():
                            if e_value is None:
                                self.logger.warning(f"Evidence item {e_idx+1} has None value for {e_key}")
                                e[e_key] = ""
                        
                        try:
                            ev_item = ET.SubElement(evidence, "match")
                            ev_item.set("domain_id", str(e.get("domain_id", "")))
                            ev_item.set("type", str(e.get("type", "")))
                            
                            if "evalue" in e and e["evalue"] is not None:
                                ev_item.set("evalue", str(e["evalue"]))
                            if "probability" in e and e["probability"] is not None:
                                ev_item.set("probability", str(e["probability"]))
                                
                            query_range = ET.SubElement(ev_item, "query_range")
                            query_range.text = str(e.get("query_range", ""))
                            
                            hit_range = ET.SubElement(ev_item, "hit_range")
                            hit_range.text = str(e.get("hit_range", ""))
                        except Exception as e_err:
                            self.logger.error(f"Error creating evidence item {e_idx+1}: {e_err}")
                else:
                    self.logger.debug(f"Domain {i+1} has no evidence")
            
            except Exception as d_err:
                self.logger.error(f"Error creating domain {i+1}: {d_err}")

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

    def _process_blast_summary(self, blast_summ_fn: str) -> Dict[str, Any]:
        """Process BLAST summary file to extract hit information"""
        try:
            tree = ET.parse(blast_summ_fn)
            root = tree.getroot()
            
            result = {
                "chain_blast": [],
                "domain_blast": [],
                "hhsearch": [],
                "self_comparison": []
            }
            
            # Process chain BLAST hits
            for chain_hit in root.findall(".//chain_blast_run/hits/hit"):
                hit_data = {
                    "num": chain_hit.get("num", ""),
                    "pdb_id": chain_hit.get("pdb_id", ""),
                    "chain_id": chain_hit.get("chain_id", ""),
                    "evalues": chain_hit.get("evalues", "").split(","),
                    "query_regions": chain_hit.findtext("query_reg", "").split(","),
                    "hit_regions": chain_hit.findtext("hit_reg", "").split(",")
                }
                result["chain_blast"].append(hit_data)
            
            # Process domain BLAST hits
            for domain_hit in root.findall(".//blast_run/hits/hit"):
                hit_data = {
                    "num": domain_hit.get("num", ""),
                    "domain_id": domain_hit.get("domain_id", ""),
                    "pdb_id": domain_hit.get("pdb_id", ""),
                    "chain_id": domain_hit.get("chain_id", ""),
                    "evalues": domain_hit.get("evalues", "").split(","),
                    "query_regions": domain_hit.findtext("query_reg", "").split(","),
                    "hit_regions": domain_hit.findtext("hit_reg", "").split(",")
                }
                result["domain_blast"].append(hit_data)
            
            # Process HHSearch hits (if available)
            for hh_hit in root.findall(".//hh_run/hits/hit"):
                hit_data = {
                    "num": hh_hit.get("num", ""),
                    "domain_id": hh_hit.get("domain_id", ""),
                    "probability": float(hh_hit.get("hh_prob", "0")),
                    "score": float(hh_hit.get("hh_score", "0")),
                    "query_region": hh_hit.findtext("query_reg", ""),
                    "hit_region": hh_hit.findtext("hit_reg", "")
                }
                result["hhsearch"].append(hit_data)
            
            # Process self-comparison (if available)
            for self_hit in root.findall(".//self_comp_run[@programs='dali']/hits/hit"):
                hit_data = {
                    "aligner": self_hit.get("aligner", ""),
                    "z_score": float(self_hit.get("z_score", "0")),
                    "query_region": self_hit.findtext("query_reg", ""),
                    "hit_region": self_hit.findtext("hit_reg", "")
                }
                result["self_comparison"].append(hit_data)
            
            return result
            
        except Exception as e:
            self.logger.error(f"Error processing BLAST summary file {blast_summ_fn}: {e}")
            return {"chain_blast": [], "domain_blast": [], "hhsearch": [], "self_comparison": []}

    def _determine_domain_boundaries(self, blast_data: Dict[str, Any], sequence_length: int, pdb_chain: str) -> List[Dict[str, Any]]:
        """Determine domain boundaries from BLAST and HHSearch results"""
        domains = []
        
        # Check if we have reference domains for this chain
        if pdb_chain in self.ref_chain_domains:
            self.logger.info(f"Found reference domains for {pdb_chain}")
            # Use reference domains as starting point
            ref_domains = self.ref_chain_domains[pdb_chain]
            for ref_domain in ref_domains:
                domains.append({
                    "range": ref_domain["range"],
                    "reference": True,
                    "domain_id": ref_domain["domain_id"],
                    "uid": ref_domain["uid"]
                })
            return domains
        
        # Step 1: Use self-comparison to identify internal repeats
        repeats = self._identify_repeats(blast_data["self_comparison"], sequence_length)
        
        # Step 2: Use domain BLAST hits to identify domains
        blast_domains = self._identify_domains_from_blast(blast_data["domain_blast"], sequence_length)
        
        # Step 3: Use HHSearch hits for more sensitive domain detection
        hhsearch_domains = self._identify_domains_from_hhsearch(blast_data["hhsearch"], sequence_length)
        
        # Step 4: Merge and refine domain boundaries
        all_candidate_domains = repeats + blast_domains + hhsearch_domains
        
        # Sort by quality (probability/e-value/score)
        all_candidate_domains.sort(key=lambda d: d.get("quality", 0), reverse=True)
        
        # Resolve overlaps and gaps
        final_domains = self._resolve_domain_boundaries(all_candidate_domains, sequence_length)
        
        return final_domains

    def _identify_repeats(self, self_comparison: List[Dict[str, Any]], sequence_length: int) -> List[Dict[str, Any]]:
        """Identify internal repeats from self-comparison results"""
        repeats = []
        
        for hit in self_comparison:
            z_score = hit.get("z_score", 0)
            
            # Filter by significance
            if z_score < self.dali_significance_threshold:
                continue
            
            query_region = hit.get("query_region", "")
            hit_region = hit.get("hit_region", "")
            
            # Create repeat domains
            repeats.append({
                "range": query_region,
                "quality": z_score,
                "type": "repeat",
                "evidence": [{
                    "type": "self_comparison",
                    "query_range": query_region,
                    "hit_range": hit_region,
                    "z_score": z_score
                }]
            })
            
            repeats.append({
                "range": hit_region,
                "quality": z_score,
                "type": "repeat",
                "evidence": [{
                    "type": "self_comparison",
                    "query_range": hit_region,
                    "hit_range": query_region,
                    "z_score": z_score
                }]
            })
        
        return repeats

    # Add this helper function to the DomainPartition class:
    def _safe_str(self, value):
        """Convert value to string safely, handling None values"""
        if value is None:
            return ""
        return str(value)

    def _identify_domains_from_blast(self, domain_blast: List[Dict[str, Any]], sequence_length: int) -> List[Dict[str, Any]]:
        """Identify domains from BLAST results"""
        domains = []
        
        for hit in domain_blast:
            # Get domain information
            domain_id = hit.get("domain_id", "")
            if domain_id == "NA":
                continue
                
            # Process query regions
            query_regions = hit.get("query_regions", [])
            if not query_regions or query_regions[0] == "":
                continue
                
            # Combine regions into a single range
            combined_range = self._combine_ranges(query_regions)
            
            # Calculate quality based on e-values
            e_values = hit.get("evalues", [])
            if e_values and e_values[0] != "":
                try:
                    min_evalue = min(float(e) for e in e_values)
                    quality = 100.0 / (1.0 + min_evalue)  # Higher quality for lower e-values
                except ValueError:
                    quality = 0
            else:
                quality = 0
            
            # Create domain
            domains.append({
                "range": combined_range,
                "quality": quality,
                "type": "blast",
                "evidence": [{
                    "type": "blast",
                    "domain_id": domain_id,
                    "query_range": combined_range,
                    "hit_range": self._combine_ranges(hit.get("hit_regions", [])),
                    "evalue": min_evalue if 'min_evalue' in locals() else 999
                }]
            })
        
        return domains

    def _identify_domains_from_hhsearch(self, hhsearch: List[Dict[str, Any]], sequence_length: int) -> List[Dict[str, Any]]:
        """Identify domains from HHSearch results"""
        domains = []
        
        for hit in hhsearch:
            # Get domain information
            domain_id = hit.get("domain_id", "")
            if domain_id == "NA":
                continue
                
            # Get probability and region
            probability = hit.get("probability", 0)
            if probability < 50:  # Threshold for HHSearch
                continue
                
            query_region = hit.get("query_region", "")
            if not query_region:
                continue
                
            # Create domain
            domains.append({
                "range": query_region,
                "quality": probability,
                "type": "hhsearch",
                "evidence": [{
                    "type": "hhsearch",
                    "domain_id": domain_id,
                    "query_range": query_region,
                    "hit_range": hit.get("hit_region", ""),
                    "probability": probability
                }]
            })
        
        return domains

    def _resolve_domain_boundaries(self, candidate_domains: List[Dict[str, Any]], sequence_length: int) -> List[Dict[str, Any]]:
        """Resolve domain boundaries by handling overlaps and gaps"""
        if not candidate_domains:
            return [{"range": f"1-{sequence_length}", "type": "full_chain"}]
        
        # Sort domains by starting position
        candidate_domains.sort(key=lambda d: self._get_start_position(d["range"]))
        
        # Initialize with the first domain
        final_domains = [candidate_domains[0]]
        
        # Process remaining domains
        for domain in candidate_domains[1:]:
            # Check overlap with existing domains
            overlaps = False
            
            for existing in final_domains:
                overlap_percentage = self._calculate_overlap_percentage(
                    domain["range"], existing["range"], sequence_length
                )
                
                if overlap_percentage > self.old_coverage_threshold:
                    # If significant overlap, keep the higher quality one
                    if domain.get("quality", 0) > existing.get("quality", 0):
                        # Replace existing with this one
                        existing.update(domain)
                    overlaps = True
                    break
            
            if not overlaps:
                final_domains.append(domain)
        
        # Sort by position again
        final_domains.sort(key=lambda d: self._get_start_position(d["range"]))
        
        # Handle gaps between domains
        if len(final_domains) > 1:
            gap_domains = []
            
            for i in range(len(final_domains) - 1):
                current_end = self._get_end_position(final_domains[i]["range"])
                next_start = self._get_start_position(final_domains[i+1]["range"])
                
                gap_size = next_start - current_end - 1
                
                if gap_size > self.gap_tol:
                    # Add a domain for the gap
                    gap_domains.append({
                        "range": f"{current_end+1}-{next_start-1}",
                        "type": "gap",
                        "quality": 0
                    })
            
            # Check for gap at the beginning
            first_start = self._get_start_position(final_domains[0]["range"])
            if first_start > 1 + self.gap_tol:
                gap_domains.append({
                    "range": f"1-{first_start-1}",
                    "type": "gap",
                    "quality": 0
                })
            
            # Check for gap at the end
            last_end = self._get_end_position(final_domains[-1]["range"])
            if sequence_length - last_end > self.gap_tol:
                gap_domains.append({
                    "range": f"{last_end+1}-{sequence_length}",
                    "type": "gap",
                    "quality": 0
                })
            
            # Add gap domains
            final_domains.extend(gap_domains)
            
            # Sort again
            final_domains.sort(key=lambda d: self._get_start_position(d["range"]))
        
        return final_domains

    def _assign_domain_classifications(self, domains: List[Dict[str, Any]], blast_data: Dict[str, Any], pdb_chain: str) -> None:
        """Assign ECOD classifications to domains"""
        self.logger.info(f"Starting domain classification assignment for {pdb_chain}")
        self.logger.debug(f"Number of domains to assign: {len(domains)}")
        
        # Debug blast data contents
        self.logger.debug(f"BLAST data summary:")
        self.logger.debug(f"  Chain BLAST hits: {len(blast_data.get('chain_blast', []))}")
        self.logger.debug(f"  Domain BLAST hits: {len(blast_data.get('domain_blast', []))}")
        self.logger.debug(f"  HHSearch hits: {len(blast_data.get('hhsearch', []))}")
        
        # Show top domain blast hits for debugging
        for i, hit in enumerate(blast_data.get('domain_blast', [])[:3]):
            self.logger.debug(f"  Domain BLAST hit {i+1}:")
            self.logger.debug(f"    Domain ID: {hit.get('domain_id', 'NA')}")
            self.logger.debug(f"    Query regions: {hit.get('query_regions', [])}")
            self.logger.debug(f"    Hit regions: {hit.get('hit_regions', [])}")
            self.logger.debug(f"    E-values: {hit.get('evalues', [])}")
        
        # First, check for reference domains
        reference_classifications = {}
        
        # Check if we have direct reference for this chain
        if pdb_chain in self.ref_chain_domains:
            self.logger.info(f"Found reference domains for {pdb_chain}")
            for ref_domain in self.ref_chain_domains[pdb_chain]:
                domain_id = ref_domain["domain_id"]
                uid = ref_domain["uid"]
                
                # Get classifications from cache or database
                classification = self._get_domain_classification(uid)
                if classification:
                    reference_classifications[domain_id] = classification
                    self.logger.debug(f"Got classification for {domain_id}: {classification}")
        
        # Assign classifications to domains
        for i, domain in enumerate(domains):
            self.logger.debug(f"Assigning classification to domain {i+1}: {domain.get('range', 'unknown_range')}")
            
            # If domain has reference, use it directly
            if "reference" in domain and domain["reference"]:
                domain_id = domain.get("domain_id", "")
                if domain_id in reference_classifications:
                    domain.update(reference_classifications[domain_id])
                    self.logger.debug(f"Updated domain with reference classification: {domain}")
                continue
            
            # Check evidence for classification
            if "evidence" not in domain:
                self.logger.debug(f"No evidence found for domain {i+1}")
                continue
                
            # Debug evidence
            self.logger.debug(f"Evidence for domain {i+1}: {len(domain.get('evidence', []))} items")
            
            # Find the best evidence (highest probability/lowest e-value)
            best_evidence = None
            best_score = 0
            
            for j, evidence in enumerate(domain["evidence"]):
                self.logger.debug(f"Evidence item {j+1}: {evidence}")
                
                domain_id = evidence.get("domain_id", "")
                if not domain_id or domain_id == "NA":
                    self.logger.debug(f"  Skipping evidence {j+1} - no valid domain_id")
                    continue
                
                # Calculate score based on evidence type
                if evidence["type"] == "hhsearch":
                    score = evidence.get("probability", 0)
                    self.logger.debug(f"  HHSearch evidence - probability: {score}")
                elif evidence["type"] == "blast":
                    e_value = evidence.get("evalue", 999)
                    score = 100.0 / (1.0 + e_value) if e_value < 10 else 0
                    self.logger.debug(f"  BLAST evidence - evalue: {e_value}, score: {score}")
                else:
                    score = 0
                    self.logger.debug(f"  Unknown evidence type: {evidence['type']}")
                
                self.logger.debug(f"  Evidence score for {domain_id}: {score}")
                
                if score > best_score:
                    best_score = score
                    best_evidence = evidence
                    self.logger.debug(f"  New best evidence: {domain_id} with score {score}")
            
            if best_evidence:
                self.logger.debug(f"Best evidence found: {best_evidence}")
                domain_id = best_evidence.get("domain_id", "")
                
                # Get classifications for this domain from cache or database
                classification = self._get_domain_classification_by_id(domain_id)
                if classification:
                    self.logger.debug(f"Classification for {domain_id}: {classification}")
                    domain.update(classification)
                    self.logger.debug(f"Domain after update: {domain}")
                else:
                    self.logger.warning(f"No classification found for domain_id {domain_id}")
                    # Add empty classification to prevent None values
                    domain.update({
                        "t_group": "",
                        "h_group": "",
                        "x_group": "",
                        "a_group": ""
                    })
                    self.logger.debug(f"Added empty classification placeholders")
            else:
                self.logger.warning(f"No best evidence found for domain {i+1}")
                # Add empty classification to prevent None values
                domain.update({
                    "t_group": "",
                    "h_group": "",
                    "x_group": "",
                    "a_group": ""
                })
                self.logger.debug(f"Added empty classification placeholders")

        # After assignment, check for missing values
        for i, domain in enumerate(domains):
            self.logger.debug(f"Final domain {i+1}: {domain}")
            for key, value in domain.items():
                if value is None:
                    self.logger.warning(f"Domain {i+1} has None value for {key}")
                
    def _calculate_overlap_percentage(self, range1: str, range2: str, sequence_length: int) -> float:
        """Calculate the percentage of overlap between two ranges"""
        # Convert ranges to sets of positions
        positions1 = self._range_to_positions(range1)
        positions2 = self._range_to_positions(range2)
        
        # Calculate overlap
        overlap = len(positions1.intersection(positions2))
        
        # Calculate coverage relative to the smaller range
        min_length = min(len(positions1), len(positions2))
        if min_length == 0:
            return 0.0
            
        return overlap / min_length

    def _range_to_positions(self, range_str: str) -> Set[int]:
        """Convert a range string to a set of positions"""
        positions = set()
        
        if not range_str:
            return positions
            
        # Handle multi-segment ranges
        segments = range_str.split(",")
        
        for segment in segments:
            if "-" in segment:
                try:
                    start, end = map(int, segment.split("-"))
                    positions.update(range(start, end + 1))
                except ValueError:
                    continue
            else:
                try:
                    positions.add(int(segment))
                except ValueError:
                    continue
                    
        return positions

    def _combine_ranges(self, ranges: List[str]) -> str:
        """Combine multiple range strings into a single range"""
        if not ranges or ranges[0] == "":
            return ""
        
        # Extract all positions
        positions = []
        for r in ranges:
            if "-" in r:
                try:
                    start, end = map(int, r.split("-"))
                    positions.extend(range(start, end + 1))
                except ValueError:
                    continue
            elif r.isdigit():
                positions.append(int(r))
        
        if not positions:
            return ""
        
        # Sort and remove duplicates
        positions = sorted(set(positions))
        
        # Convert to ranges
        result = []
        start = positions[0]
        prev = start
        
        for pos in positions[1:]:
            if pos > prev + 1:
                result.append(f"{start}-{prev}")
                start = pos
            prev = pos
        
        result.append(f"{start}-{prev}")
        return ",".join(result)

    def _get_end_position(self, range_str: str) -> int:
        """Get the end position from a range string"""
        if "-" in range_str:
            parts = range_str.split("-")
            return int(parts[1])
        elif "," in range_str:
            # Multi-segment range - get the last segment
            segments = range_str.split(",")
            return self._get_end_position(segments[-1])
        else:
            try:
                return int(range_str)
            except ValueError:
                    return 0

    def partition_domains_assembly(self, pdb_id: str, chain_ids: List[str], dump_dir: str, input_mode: str, reference: str, blast_only: bool = False) -> str:
        """Partition domains for a multi-chain assembly"""
        # Load reference data if not already loaded
        if not self.ref_range_cache:
            self.load_reference_data(reference)
        
        # Define paths
        pdb_chains = f"{pdb_id}_{''.join(chain_ids)}"
        asm_dir = os.path.join(dump_dir, pdb_chains)
        domain_prefix = "domains_v12"
        domain_fn = os.path.join(asm_dir, f"{domain_prefix}.{pdb_chains}.{reference}.xml")
        
        if os.path.exists(domain_fn) and not self.config.get('force_overwrite', False):
            self.logger.warning(f"Domain file {domain_fn} already exists, skipping...")
            return domain_fn
        
        # Create domain document for assembly
        domains_doc = ET.Element("domain_doc")
        domains_doc.set("pdb", self._safe_str(pdb_id))
        domains_doc.set("chain", self._safe_str(chain_id))
        domains_doc.set("reference", self._safe_str(reference))

        # Process each chain separately
        all_domains = []
        
        for chain_id in chain_ids:
            pdb_chain = f"{pdb_id}_{chain_id}"
            chain_dir = os.path.join(dump_dir, pdb_chain)
            
            # Get blast summary for this chain
            blast_summ_fn = os.path.join(chain_dir, 
                                       f"{pdb_chain}.{reference}.blast_summ{''.join(['.blast_only' if blast_only else ''])}.xml")
            
            if not os.path.exists(blast_summ_fn):
                self.logger.warning(f"Blast summary file not found for {pdb_chain}: {blast_summ_fn}")
                continue
            
            # Get chain sequence
            fasta_path = os.path.join(chain_dir, f"{pdb_chain}.fa")
            sequence = self._read_fasta_sequence(fasta_path)
            sequence_length = len(sequence) if sequence else 0
            
            if not sequence:
                self.logger.warning(f"Failed to read sequence for {pdb_chain}")
                continue
            
            # Process blast summary
            blast_data = self._process_blast_summary(blast_summ_fn)
            
            # Determine domain boundaries
            domains = self._determine_domain_boundaries(blast_data, sequence_length, pdb_chain)
            
            # Assign classifications
            self._assign_domain_classifications(domains, blast_data, pdb_chain)
            
            # Store chain ID with each domain
            for domain in domains:
                domain["chain"] = chain_id
            
            all_domains.extend(domains)
        
        # Check for inter-chain domains
        assembly_domains = self._detect_assembly_domains(all_domains, pdb_id, chain_ids)
        
        # Add domains to the document
        domain_list = ET.SubElement(domains_doc, "domain_list")
        
        for d in assembly_domains:
            domain_elem = ET.SubElement(domain_list, "domain")
            domain_elem.set("pdb", self._safe_str(pdb_id))
            domain_elem.set("chain", self._safe_str(chain_id))
            domain_elem.set("range", self._safe_str(d["range"]))
            
            # Add classification attributes if present
            for attr in ["t_group", "h_group", "x_group", "a_group"]:
                if attr in d:
                    domain_elem.set(attr, d[attr])
            
            # Add flags if present
            for flag in ["is_manual_rep", "is_f70", "is_f40", "is_f99"]:
                if flag in d and d[flag]:
                    domain_elem.set(flag, "true")
            
            # For assembly domains
            if "assembly" in d and d["assembly"]:
                domain_elem.set("assembly", "true")
                domain_elem.set("chains", d.get("chains", ""))
            
            # Add range element
            range_elem = ET.SubElement(domain_elem, "range")
            range_elem.text = d["range"]
            
            # Add classification evidence
            if "evidence" in d:
                evidence = ET.SubElement(domain_elem, "evidence")
                for e in d["evidence"]:
                    ev_item = ET.SubElement(evidence, "match")
                    ev_item.set("domain_id", e.get("domain_id", ""))
                    ev_item.set("type", e.get("type", ""))
                    
                    if "evalue" in e:
                        ev_item.set("evalue", str(e["evalue"]))
                    if "probability" in e:
                        ev_item.set("probability", str(e["probability"]))
                        
                    query_range = ET.SubElement(ev_item, "query_range")
                    query_range.text = e.get("query_range", "")
                    
                    hit_range = ET.SubElement(ev_item, "hit_range")
                    hit_range.text = e.get("hit_range", "")
        
        # Write output file
        os.makedirs(os.path.dirname(domain_fn), exist_ok=True)
        tree = ET.ElementTree(domains_doc)
        tree.write(domain_fn, encoding='utf-8', xml_declaration=True)
        
        self.logger.info(f"Created assembly domain partition file: {domain_fn}")
        return domain_fn

    def _detect_assembly_domains(self, all_domains: List[Dict[str, Any]], 
                               pdb_id: str, chain_ids: List[str]
    ) -> List[Dict[str, Any]]:
        """Detect inter-chain domains in an assembly"""
        # This is a simplified implementation - in a real-world scenario,
        # we would need to analyze inter-chain contacts and structure
        
        # Start with individual chain domains
        assembly_domains = []
        
        # Add chain information to the domain range
        for domain in all_domains:
            # Make a copy to avoid modifying the original
            domain_copy = domain.copy()
            
            if "chain" in domain_copy:
                chain = domain_copy["chain"]
                original_range = domain_copy["range"]
                domain_copy["range"] = f"{chain}:{original_range}"
            
            assembly_domains.append(domain_copy)
        
        # In a full implementation, we would now:
        # 1. Look for domains from different chains that interact
        # 2. Analyze inter-chain contacts
        # 3. Determine if domains form multi-chain assembly domains
        
        # Example: detect domains from different chains with matching classifications
        if len(chain_ids) > 1:
            # Group domains by classification
            classification_groups = {}
            
            for domain in all_domains:
                if "t_group" in domain and "h_group" in domain:
                    key = f"{domain.get('t_group')}_{domain.get('h_group')}"
                    if key not in classification_groups:
                        classification_groups[key] = []
                    classification_groups[key].append(domain)
            
            # Look for matching classifications across chains
            for key, group in classification_groups.items():
                # Check if domains are from different chains
                chains = set(domain.get("chain", "") for domain in group)
                if len(chains) > 1:
                    # Create an assembly domain
                    chain_ranges = []
                    evidence = []
                    
                    for domain in group:
                        chain = domain.get("chain", "")
                        range_str = domain.get("range", "")
                        chain_ranges.append(f"{chain}:{range_str}")
                        
                        # Collect evidence
                        if "evidence" in domain:
                            evidence.extend(domain["evidence"])
                    
                    # Create combined domain
                    assembly_domain = {
                        "range": ",".join(chain_ranges),
                        "chains": ",".join(chains),
                        "assembly": True,
                        "t_group": group[0].get("t_group"),
                        "h_group": group[0].get("h_group"),
                        "x_group": group[0].get("x_group"),
                        "a_group": group[0].get("a_group"),
                        "evidence": evidence
                    }
                    
                    assembly_domains.append(assembly_domain)
        
        return assembly_domains