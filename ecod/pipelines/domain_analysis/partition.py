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
                        dump_dir: str, reference: str, blast_only: bool = False
    ) -> bool:
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

    def register_domain_file(self, process_id, file_path):
        """Register domain partition file in database with proper duplicate handling"""
        try:
            # Check if record already exists
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = 'domain_partition'
            """
            
            existing = self.db.execute_query(query, (process_id,))
            
            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": file_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    },
                    "id = %s",
                    (existing[0][0],)
                )
            else:
                # Insert new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "domain_partition",
                        "file_path": file_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    }
                )
        except Exception as e:
            self.logger.warning(f"Error registering domain file: {e}")
        
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

    def partition_domains(self, pdb_id: str, chain_id: str, dump_dir: str, input_mode: str, reference: str, blast_only: bool = False) -> str:
        """Partition domains for a single protein chain"""
        # Load reference data if not already loaded
        if not self.ref_range_cache:
            self.load_reference_data(reference)
        
        # Define paths
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Define the domains directory
        domains_dir = os.path.join(dump_dir, "domains")
        os.makedirs(domains_dir, exist_ok=True)
        
        # Set the domain output file path with the new naming convention
        domain_prefix = "domains_v14"
        domain_fn = os.path.join(domains_dir, f"{pdb_chain}.{reference}.{domain_prefix}.xml")
        
        if os.path.exists(domain_fn) and not self.config.get('force_overwrite', False):
            self.logger.warning(f"Domain file {domain_fn} already exists, skipping...")
            return domain_fn

        # Create domain document
        domains_doc = ET.Element("domain_doc")
        domains_doc.set("pdb", self._safe_str(pdb_id))
        domains_doc.set("chain", self._safe_str(chain_id))
        domains_doc.set("reference", self._safe_str(reference))
        
        # Look for domain summary file in the domains directory  
        blast_summ_fn = os.path.join(domains_dir, f"{pdb_chain}.domain_summary.xml")
        
        # Check for FASTA file
        fastas_dir = os.path.join(dump_dir, "fastas", "batch_0")  # Adjust batch folder as needed
        fasta_fn = os.path.join(fastas_dir, f"{pdb_chain}.fasta")
        
        # Log file paths for debugging
        self.logger.info(f"Found FASTA file at: {fasta_fn}")
        
        sequence = self._read_fasta_sequence(fasta_fn)
        if not sequence:
            self.logger.error(f"Failed to read sequence from {fasta_fn}")
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

            # Write output file
        os.makedirs(os.path.dirname(domain_fn), exist_ok=True)
        tree = ET.ElementTree(domains_doc)
        tree.write(domain_fn, encoding='utf-8', xml_declaration=True)
        
        self.logger.info(f"Created domain partition file: {domain_fn}")
        return domain_fn

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
        """
        Process domain summary XML file and extract all relevant data
        
        Args:
            blast_summ_fn: Path to domain summary XML file
            
        Returns:
            Dict with parsed data
        """
        logger = logging.getLogger("ecod.domain_partition")
        logger.info(f"Processing domain summary: {blast_summ_fn}")
        
        # Check file existence
        if not os.path.exists(blast_summ_fn):
            logger.error(f"Domain summary file not found: {blast_summ_fn}")
            return {"error": "File not found"}
        
        # Initialize result structure
        result = {
            "sequence": "",
            "sequence_length": 0,
            "chain_id": "",
            "blast_hits": [],
            "domain_blast_hits": []
        }
        
        try:
            # Parse XML
            tree = ET.parse(blast_summ_fn)
            root = tree.getroot()
            
            # Log root tag and structure for debugging
            logger.debug(f"XML root tag: {root.tag}, attributes: {root.attrib}")
            
            # Get blast_summ element
            blast_summ = root.find(".//blast_summ")
            if blast_summ is None:
                logger.error("No blast_summ element found in XML")
                return {"error": "Invalid XML format: No blast_summ element"}
                
            # Extract chain ID from attributes
            if blast_summ.get("pdb") and blast_summ.get("chain"):
                result["chain_id"] = f"{blast_summ.get('pdb')}_{blast_summ.get('chain')}"
                logger.debug(f"Found chain_id in attributes: {result['chain_id']}")
            
            # Extract query length
            query_len_elem = root.find(".//query_len")
            if query_len_elem is not None and query_len_elem.text:
                try:
                    result["sequence_length"] = int(query_len_elem.text.strip())
                    logger.debug(f"Found sequence length: {result['sequence_length']}")
                except ValueError:
                    logger.warning(f"Invalid sequence length: {query_len_elem.text}")
            
            # Process chain-level BLAST hits
            chain_blast_hits = []
            for hit_elem in root.findall(".//chain_blast_run/hits/hit"):
                hit = dict(hit_elem.attrib)
                
                # Extract query and hit regions
                query_reg_elem = hit_elem.find("query_reg")
                if query_reg_elem is not None and query_reg_elem.text:
                    hit["range"] = query_reg_elem.text.strip()
                    hit["range_parsed"] = self._parse_range(hit["range"])
                
                # Convert evalue
                if "evalues" in hit:
                    try:
                        evalue_str = hit["evalues"]
                        hit["evalue"] = float(evalue_str) if evalue_str != "0.0" else 1e-200
                    except (ValueError, TypeError):
                        hit["evalue"] = 1.0  # Default to a poor E-value
                
                # Add identity and coverage (not directly available, estimate from sequences)
                query_seq_elem = hit_elem.find("query_seq")
                hit_seq_elem = hit_elem.find("hit_seq")
                
                if query_seq_elem is not None and hit_seq_elem is not None:
                    query_seq = query_seq_elem.text.strip()
                    hit_seq = hit_seq_elem.text.strip()
                    
                    if query_seq and hit_seq and len(query_seq) == len(hit_seq):
                        # Calculate identity percentage
                        matches = sum(1 for q, h in zip(query_seq, hit_seq) if q == h)
                        identity = (matches / len(query_seq)) * 100
                        hit["identity"] = identity
                        
                        # Estimate coverage from region
                        if "range_parsed" in hit and hit["range_parsed"]:
                            total_covered = sum(end - start + 1 for start, end in hit["range_parsed"])
                            hit["query_coverage"] = (total_covered / result["sequence_length"]) * 100 if result["sequence_length"] > 0 else 0
                
                chain_blast_hits.append(hit)
            
            result["blast_hits"] = chain_blast_hits
            logger.debug(f"Found {len(chain_blast_hits)} chain BLAST hits")
            
            # Process domain-level BLAST hits
            domain_blast_hits = []
            for hit_elem in root.findall(".//blast_run/hits/hit"):
                hit = dict(hit_elem.attrib)
                
                # Extract query region
                query_reg_elem = hit_elem.find("query_reg")
                if query_reg_elem is not None and query_reg_elem.text:
                    hit["range"] = query_reg_elem.text.strip()
                    hit["range_parsed"] = self._parse_range(hit["range"])
                
                # Convert evalue
                if "evalues" in hit:
                    try:
                        evalue_str = hit["evalues"]
                        hit["evalue"] = float(evalue_str) if evalue_str != "0.0" else 1e-200
                    except (ValueError, TypeError):
                        hit["evalue"] = 1.0  # Default to a poor E-value
                
                # Add domain ID if present
                if "domain_id" in hit:
                    hit["target_id"] = hit["domain_id"]
                
                # Add identity and coverage (estimate from region)
                if "range_parsed" in hit and hit["range_parsed"]:
                    # Estimate identity (not directly available)
                    hit["identity"] = 90.0  # Use a reasonable default
                    
                    # Estimate coverage from region
                    total_covered = sum(end - start + 1 for start, end in hit["range_parsed"])
                    hit["query_coverage"] = (total_covered / result["sequence_length"]) * 100 if result["sequence_length"] > 0 else 0
                
                domain_blast_hits.append(hit)
            
            result["domain_blast_hits"] = domain_blast_hits
            logger.debug(f"Found {len(domain_blast_hits)} domain BLAST hits")
            
            return result
            
        except ET.ParseError as e:
            logger.error(f"XML parsing error: {str(e)}")
            return {"error": f"XML parsing error: {str(e)}"}
        except Exception as e:
            logger.error(f"Unexpected error processing domain summary: {str(e)}", exc_info=True)
            return {"error": f"Error processing domain summary: {str(e)}"}

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

    def _parse_range(self, range_str: str) -> List[Tuple[int, int]]:
        """
        Parse a range string like "1-100,150-200" into a list of tuples [(1,100), (150,200)]
        
        Args:
            range_str: Range string to parse
            
        Returns:
            List of (start, end) tuples
        """
        if not range_str or range_str.strip() == "":
            return []
        
        ranges = []
        parts = range_str.split(',')
        
        for part in parts:
            if '-' in part:
                try:
                    start, end = part.split('-')
                    start_num = int(re.sub(r'[^0-9]', '', start))
                    end_num = int(re.sub(r'[^0-9]', '', end))
                    ranges.append((start_num, end_num))
                except (ValueError, IndexError):
                    pass
        
        return ranges

    def _identify_domains_from_blast(self, blast_hits: List[Dict[str, Any]], sequence_length: int) -> List[Dict[str, Any]]:
        """
        Identify domain boundaries from BLAST hits
        
        Args:
            blast_hits: List of BLAST hits with parsed attributes
            sequence_length: Length of the protein sequence
            
        Returns:
            List of domain dictionaries with boundaries
        """
        logger = logging.getLogger("ecod.domain_partition")
        
        if not blast_hits or sequence_length <= 0:
            logger.warning("No BLAST hits or invalid sequence length")
            return []
        
        # Define significance thresholds (can be moved to config)
        thresholds = {
            "evalue": 1e-3,
            "identity": 30.0,
            "query_coverage": 50.0
        }
        
        # Filter significant hits
        significant_hits = []
        for hit in blast_hits:
            # Skip hits without required attributes
            if not all(key in hit for key in ["evalue", "identity", "query_coverage"]):
                continue
                
            # Skip hits without ranges
            if "range_parsed" not in hit or not hit["range_parsed"]:
                continue
                
            # Check significance
            if (hit["evalue"] <= thresholds["evalue"] and
                hit["identity"] >= thresholds["identity"] and
                hit["query_coverage"] >= thresholds["query_coverage"]):
                significant_hits.append(hit)
        
        logger.info(f"Found {len(significant_hits)}/{len(blast_hits)} significant BLAST hits")
        
        if not significant_hits:
            logger.warning("No significant BLAST hits found")
            return []
        
        # Analyze coverage with significant hits
        position_coverage = [0] * (sequence_length + 1)  # 1-indexed
        hit_regions = []
        
        for hit in significant_hits:
            hit_region = {
                "target_id": hit.get("target_id", "unknown"),
                "evalue": hit["evalue"],
                "identity": hit["identity"],
                "query_coverage": hit["query_coverage"],
                "ranges": []
            }
            
            for start, end in hit["range_parsed"]:
                # Track coverage
                for i in range(max(1, start), min(sequence_length + 1, end + 1)):
                    position_coverage[i-1] += 1
                
                hit_region["ranges"].append({"start": start, "end": end})
            
            hit_regions.append(hit_region)
        
        # Find contiguous regions
        regions = []
        region_start = None
        
        for i in range(sequence_length):
            if position_coverage[i] > 0:
                if region_start is None:
                    region_start = i + 1
            else:
                if region_start is not None:
                    regions.append({
                        "start": region_start,
                        "end": i,
                        "size": i - region_start + 1
                    })
                    region_start = None
        
        # Add final region if exists
        if region_start is not None:
            regions.append({
                "start": region_start,
                "end": sequence_length,
                "size": sequence_length - region_start + 1
            })
        
        # Merge small gaps between regions (< 30 residues)
        merged_regions = []
        if regions:
            merged_regions.append(regions[0])
            
            for region in regions[1:]:
                prev_region = merged_regions[-1]
                
                if region["start"] - prev_region["end"] <= 30:
                    # Merge regions
                    prev_region["end"] = region["end"]
                    prev_region["size"] = prev_region["end"] - prev_region["start"] + 1
                else:
                    # Add as separate region
                    merged_regions.append(region)
        
        logger.info(f"Identified {len(merged_regions)} domains from BLAST hits")
        
        # Log details for each identified domain
        for i, region in enumerate(merged_regions):
            logger.debug(f"Domain {i+1}: {region['start']}-{region['end']} (size: {region['size']})")
        
        return merged_regions
    
    def _identify_domains_from_hhsearch(self, hhsearch_hits: List[Dict[str, Any]], sequence_length: int) -> List[Dict[str, Any]]:
        """
        Identify domain boundaries from HHSearch hits
        
        Args:
            hhsearch_hits: List of HHSearch hits with parsed attributes
            sequence_length: Length of the protein sequence
            
        Returns:
            List of domain dictionaries with boundaries
        """
        logger = logging.getLogger("ecod.domain_partition")
        
        if not hhsearch_hits or sequence_length <= 0:
            logger.warning("No HHSearch hits or invalid sequence length")
            return []
        
        # Define significance thresholds (can be moved to config)
        thresholds = {
            "probability": 90.0,
            "evalue": 1e-3
        }
        
        # Filter significant hits
        significant_hits = []
        for hit in hhsearch_hits:
            # Skip hits without required attributes
            if not all(key in hit for key in ["probability", "evalue"]):
                continue
                
            # Skip hits without ranges
            if "range_parsed" not in hit or not hit["range_parsed"]:
                continue
                
            # Check significance (either condition)
            if (hit["probability"] >= thresholds["probability"] or
                hit["evalue"] <= thresholds["evalue"]):
                significant_hits.append(hit)
        
        logger.info(f"Found {len(significant_hits)}/{len(hhsearch_hits)} significant HHSearch hits")
        
        if not significant_hits:
            logger.warning("No significant HHSearch hits found")
            return []
        
        # Analyze coverage with significant hits
        position_coverage = [0] * (sequence_length + 1)  # 1-indexed
        hit_regions = []
        
        for hit in significant_hits:
            hit_region = {
                "target_id": hit.get("target_id", "unknown"),
                "probability": hit["probability"],
                "evalue": hit["evalue"],
                "ranges": []
            }
            
            for start, end in hit["range_parsed"]:
                # Track coverage
                for i in range(max(1, start), min(sequence_length + 1, end + 1)):
                    position_coverage[i-1] += 1
                
                hit_region["ranges"].append({"start": start, "end": end})
            
            hit_regions.append(hit_region)
        
        # Find contiguous regions
        regions = []
        region_start = None
        
        for i in range(sequence_length):
            if position_coverage[i] > 0:
                if region_start is None:
                    region_start = i + 1
            else:
                if region_start is not None:
                    regions.append({
                        "start": region_start,
                        "end": i,
                        "size": i - region_start + 1
                    })
                    region_start = None
        
        # Add final region if exists
        if region_start is not None:
            regions.append({
                "start": region_start,
                "end": sequence_length,
                "size": sequence_length - region_start + 1
            })
        
        # Merge small gaps between regions (< 30 residues)
        merged_regions = []
        if regions:
            merged_regions.append(regions[0])
            
            for region in regions[1:]:
                prev_region = merged_regions[-1]
                
                if region["start"] - prev_region["end"] <= 30:
                    # Merge regions
                    prev_region["end"] = region["end"]
                    prev_region["size"] = prev_region["end"] - prev_region["start"] + 1
                else:
                    # Add as separate region
                    merged_regions.append(region)
        
        logger.info(f"Identified {len(merged_regions)} domains from HHSearch hits")
        
        # Log details for each identified domain
        for i, region in enumerate(merged_regions):
            logger.debug(f"Domain {i+1}: {region['start']}-{region['end']} (size: {region['size']})")
        
        return merged_regions
    
    def _analyze_chainwise_hits_for_domains(self, chain_blast_hits):
        """
        Analyze chainwise BLAST hits to identify multiple domains from reference chains
        
        Args:
            chain_blast_hits: List of chain BLAST hits
            
        Returns:
            List of domain candidates with proper structure for _determine_domain_boundaries
        """

        # Add debugging to see what we're working with
        self.logger.debug(f"Chain BLAST hits count: {len(chain_blast_hits) if chain_blast_hits else 0}")

        if chain_blast_hits and len(chain_blast_hits) > 0:
            self.logger.debug(f"First chain hit keys: {chain_blast_hits[0].keys()}")

        domain_candidates = []
        
        for hit in chain_blast_hits:
            hit_pdb_id = hit.get("pdb_id", "")
            hit_chain_id = hit.get("chain_id", "")
            source_id = f"{hit_pdb_id}_{hit_chain_id}"

            self.logger.debug(f"Checking reference domains for chain: {source_id}")             
            
            # Get reference domains for this chain
            reference_domains = self._get_reference_chain_domains(source_id)
            self.logger.debug(f"Found {len(reference_domains)} reference domains for {source_id}")
            
            # Log query and hit regions
            if "query_regions" in hit and "hit_regions" in hit:
                self.logger.debug(f"  Query regions: {hit.get('query_regions', '')}")
                self.logger.debug(f"  Hit regions: {hit.get('hit_regions', '')}")
                
            self.logger.info(f"Found multi-domain reference chain: {source_id} with {len(reference_domains)} domains")
            
            # Parse query and hit regions from alignment
            query_regions = []
            hit_regions = []
            
            if "query_regions" in hit and "hit_regions" in hit:
                query_region_strs = hit.get("query_regions", "").split(",")
                hit_region_strs = hit.get("hit_regions", "").split(",")
                
                for q_range, h_range in zip(query_region_strs, hit_region_strs):
                    if not q_range or not h_range:
                        continue
                    
                    try:
                        q_start, q_end = map(int, q_range.split("-"))
                        h_start, h_end = map(int, h_range.split("-"))
                        
                        query_regions.append((q_start, q_end))
                        hit_regions.append((h_start, h_end))
                    except (ValueError, AttributeError):
                        continue
            
            # Skip if no alignment regions were parsed
            if not query_regions or not hit_regions:
                continue
            
            # Map reference domains to query coordinates
            for domain in reference_domains:
                domain_range = domain.get("range", "")
                domain_ranges = self._parse_range(domain_range)
                
                # Skip domains without a valid range
                if not domain_ranges:
                    continue
                    
                # Find corresponding query positions
                mapped_ranges = []
                
                for d_start, d_end in domain_ranges:
                    for i, (h_start, h_end) in enumerate(hit_regions):
                        # Check for overlap
                        if h_end < d_start or h_start > d_end:
                            continue
                            
                        # Calculate overlap region
                        overlap_start = max(d_start, h_start)
                        overlap_end = min(d_end, h_end)
                        
                        if overlap_end >= overlap_start:
                            # Map to query coordinates
                            q_start, q_end = query_regions[i]
                            
                            # Calculate proportion
                            h_len = h_end - h_start
                            if h_len == 0:
                                continue
                                
                            q_len = q_end - q_start
                            
                            # Map start position
                            start_ratio = (overlap_start - h_start) / h_len
                            q_mapped_start = int(q_start + start_ratio * q_len)
                            
                            # Map end position
                            end_ratio = (overlap_end - h_start) / h_len
                            q_mapped_end = int(q_start + end_ratio * q_len)
                            
                            mapped_ranges.append((q_mapped_start, q_mapped_end))
                
                # Only add domain if we could map it
                if mapped_ranges:
                    # Calculate domain start and end
                    # Use the earliest start and latest end from all mapped ranges
                    start = min(r[0] for r in mapped_ranges)
                    end = max(r[1] for r in mapped_ranges)
                    
                    # Create domain candidate with structure expected by _determine_domain_boundaries
                    domain_candidate = {
                        "start": start,
                        "end": end,
                        "size": end - start + 1,
                        "t_group": domain.get("t_group", ""),
                        "h_group": domain.get("h_group", ""),
                        "x_group": domain.get("x_group", ""),
                        "a_group": domain.get("a_group", ""),
                        "evidence": [{
                            "type": "chain_blast",
                            "domain_id": domain.get("domain_id", ""),
                            "query_range": f"{start}-{end}",
                            "hit_range": domain_range
                        }]
                    }
                    
                    domain_candidates.append(domain_candidate)
                    self.logger.info(f"Mapped domain from {source_id}: {start}-{end} ({domain_candidate['size']} residues)")
        
        self.logger.info(f"Found {len(domain_candidates)} domain candidates from chain BLAST hits")
        return domain_candidates

    def _determine_domain_boundaries(self, blast_data: Dict[str, Any], sequence_length: int, pdb_chain: str) -> List[Dict[str, Any]]:
        """
        Determine final domain boundaries using all available evidence
        
        Args:
            blast_data: Dictionary with processed domain summary data
            sequence_length: Length of the protein sequence
            pdb_chain: PDB chain identifier
            
        Returns:
            List of domain dictionaries
        """
        logger = logging.getLogger("ecod.domain_partition")
        logger.info(f"Determining domain boundaries for {pdb_chain} (length: {sequence_length})")
        
        # Skip if essential data is missing
        if "error" in blast_data or sequence_length <= 0:
            logger.warning(f"Invalid input data for {pdb_chain}")
            return [{
                "domain_num": 1,
                "range": f"1-{sequence_length}",
                "size": sequence_length,
                "confidence": "low",
                "source": "whole_chain",
                "reason": "Missing or invalid input data"
            }]
        
        # Get domains from different sources
        blast_domains = self._identify_domains_from_blast(blast_data.get("blast_hits", []), sequence_length)
        hhsearch_domains = self._identify_domains_from_hhsearch(blast_data.get("hhsearch_hits", []), sequence_length)
        # Add self-comparison and domain BLAST analysis if needed
        
        logger.debug(f"BLAST data keys: {blast_data.keys()}")
        if "chain_blast_hits" in blast_data:
            logger.debug(f"Chain BLAST hits count: {len(blast_data['chain_blast_hits'])}")
        
            # Print example of first hit structure
            if blast_data["chain_blast_hits"]:
                logger.debug(f"First chain hit keys: {blast_data['chain_blast_hits'][0].keys()}")

        chain_domain_candidates = self._analyze_chainwise_hits_for_domains(blast_data.get("blast_hits", []))

        # Log domain counts from each source
        logger.info(f"Domain candidates: {len(blast_domains)} from BLAST, {len(hhsearch_domains)} from HHSearch")
        
        # Consolidate domains from different sources
        all_domains = []
        all_domains.extend([{**d, "source": "blast"} for d in blast_domains])
        all_domains.extend([{**d, "source": "hhsearch"} for d in hhsearch_domains])
        all_domains.extend([{**d, "source": "chain_blast_domain"} for d in chain_domain_candidates])
        
        # Track coverage with all evidence
        position_evidence = [[] for _ in range(sequence_length + 1)]  # 1-indexed
        
        for domain in all_domains:
            source = domain["source"]
            for i in range(max(1, domain["start"]), min(sequence_length + 1, domain["end"] + 1)):
                if source not in position_evidence[i-1]:
                    position_evidence[i-1].append(source)
        
        # Find contiguous regions with consistent evidence
        final_regions = []
        region_start = None
        region_evidence = set()
        
        for i in range(sequence_length):
            evidence = set(position_evidence[i])
            
            if evidence:
                if region_start is None:
                    region_start = i + 1
                    region_evidence = evidence
                elif evidence != region_evidence:
                    # Evidence changed, finalize current region
                    final_regions.append({
                        "start": region_start,
                        "end": i,
                        "size": i - region_start + 1,
                        "evidence": list(region_evidence)
                    })
                    # Start new region
                    region_start = i + 1
                    region_evidence = evidence
            else:
                if region_start is not None:
                    final_regions.append({
                        "start": region_start,
                        "end": i,
                        "size": i - region_start + 1,
                        "evidence": list(region_evidence)
                    })
                    region_start = None
                    region_evidence = set()
        
        # Add final region if exists
        if region_start is not None:
            final_regions.append({
                "start": region_start,
                "end": sequence_length,
                "size": sequence_length - region_start + 1,
                "evidence": list(region_evidence)
            })
        
        # Merge small gaps between regions (< 30 residues)
        merged_regions = []
        if final_regions:
            merged_regions.append(final_regions[0])
            
            for region in final_regions[1:]:
                prev_region = merged_regions[-1]
                
                if region["start"] - prev_region["end"] <= 30:
                    # Merge regions, combining evidence
                    prev_region["end"] = region["end"]
                    prev_region["size"] = prev_region["end"] - prev_region["start"] + 1
                    prev_region["evidence"] = list(set(prev_region["evidence"] + region["evidence"]))
                else:
                    # Add as separate region
                    merged_regions.append(region)
        
        # If no regions found or all regions too small, use whole chain
        if not merged_regions or all(r["size"] < 30 for r in merged_regions):
            logger.warning(f"No significant domains found for {pdb_chain}, using whole chain")
            return [{
                "domain_num": 1,
                "range": f"1-{sequence_length}",
                "size": sequence_length,
                "confidence": "low",
                "source": "whole_chain",
                "reason": "No significant domains found"
            }]
        
        # Prepare final domain list
        domains = []
        for i, region in enumerate(merged_regions):
            # Calculate confidence based on evidence
            confidence = "high" if len(region["evidence"]) > 1 else "medium"
            if region["size"] < 50:
                confidence = "low"
            
            # Format range string
            range_str = f"{region['start']}-{region['end']}"
            
            # Record sources
            source = "+".join(region["evidence"])
            reason = f"Domain supported by {len(region['evidence'])} evidence types: {source}"
            
            domains.append({
                "domain_num": i + 1,
                "range": range_str,
                "size": region["size"],
                "confidence": confidence,
                "source": source,
                "reason": reason
            })
        
        logger.info(f"Final domain boundaries for {pdb_chain}: {len(domains)} domains")
        for domain in domains:
            logger.debug(f"Domain {domain['domain_num']}: {domain['range']} ({domain['size']} residues), "
                         f"confidence: {domain['confidence']}, source: {domain['source']}")
        
        return domains

    def _get_reference_chain_domains(self, source_id):
        """Get domain information for a reference chain"""
        self.logger.debug(f"Looking up reference domains for chain: {source_id}")
        
        if source_id in self.ref_chain_domains:
            domains = self.ref_chain_domains[source_id]
            self.logger.debug(f"Found {len(domains)} domains for {source_id} in reference cache")

            # Print details of first few domains for debugging
            for i, domain in enumerate(domains[:3]):  # Limit to first 3 domains
                self.logger.debug(f"  Domain {i+1}: {domain.get('domain_id', 'unknown')}, "
                              f"range: {domain.get('range', 'unknown')}, "
                              f"t_group: {domain.get('t_group', 'unknown')}, "
                              f"h_group: {domain.get('h_group', 'unknown')}")


            return domains
        
        # Try lowercase version
        lower_source_id = source_id.lower()
        if lower_source_id in self.ref_chain_domains:
            domains = self.ref_chain_domains[lower_source_id]
            self.logger.debug(f"Found {len(domains)} domains for {lower_source_id} (lowercase) in reference cache")
            return domains
        
        self.logger.debug(f"No reference domains found for {source_id}")
        return []

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