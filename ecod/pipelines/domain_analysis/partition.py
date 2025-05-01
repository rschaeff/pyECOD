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

from ecod.exceptions import PipelineError, FileOperationError
from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.models.pipeline import BlastHit, HHSearchHit
from ecod.models.pipeline.domain_analysis import DomainPartitionResult
from ecod.utils.xml_utils import ensure_dict, ensure_list_of_dicts
from ecod.utils.path_utils import get_standardized_paths, get_all_evidence_paths
from ecod.utils.file import find_fasta_file, read_sequence_from_fasta


class DomainPartition:
    """Determine domain boundaries and classifications from search results"""
    
    def __init__(self, context=None):
        """Initialize with configuration"""
        self.context = context or ApplicationContext()
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
        db_config = self.context.config_manager.get_db_config()
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
                    self.register_domain_file(process_id, os.path.relpath(domain_file, dump_dir), db)
                    
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

    def process_batch(self, batch_id: int, batch_path: str, reference: str,
                     blast_only: bool = False, limit: int = None,
                     reps_only: bool = False
    ) -> List[DomainPartitionResult]:
        """
        Process domains for a batch of proteins

        Args:
            batch_id: Batch ID
            batch_path: Batch base path
            reference: Reference version
            blast_only: Whether to use only BLAST results (no HHSearch)
            limit: Maximum number of proteins to process
            reps_only: Whether to process only representative proteins

        Returns:
            List of DomainPartitionResult models
        """
        logger = logging.getLogger(__name__)
        results = []

        # Get proteins to process (implementation depends on your database setup)
        proteins = self._get_proteins_to_process(batch_id, limit, reps_only)

        if reps_only:
            logger.info("Filtering for representative proteins (processes) only")

        # Process each protein
        for protein in proteins:
            pdb_id = protein["pdb_id"]
            chain_id = protein["chain_id"]

            # Find domain summary file
            domain_summary_path = self._find_domain_summary(
                batch_path, pdb_id, chain_id, reference, blast_only
            )

            if not domain_summary_path:
                logger.warning(f"No domain summary found for {pdb_id}_{chain_id}")
                result = DomainPartitionResult(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    reference=reference,
                    success=False,
                    error="Domain summary not found"
                )
                results.append(result)
                continue

            # Process domains
            result = self.process_domains(
                pdb_id, chain_id, domain_summary_path, batch_path, reference
            )

            # Store result
            results.append(result)

            try:
                # Update database status if needed
                self._update_process_status(protein["process_id"], result)
            except Exception as e:
                logger.error(f"Error updating process status: {str(e)}")

        logger.info(f"Processed domains for {len(proteins)} proteins from batch {batch_id}")
        return results

    def register_domain_file(self, process_id, file_path, db):
        """Register domain partition file in database with proper duplicate handling"""
        try:
            # Check if record already exists
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = 'domain_partition'
            """
            
            existing = db.execute_query(query, (process_id,))
            
            if existing:
                # Update existing record
                db.update(
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
                db.insert(
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
        db_config = self.context.config_manager.get_db_config()
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
        db_config = self.context.config_manager.get_db_config()
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
        db_config = self.context.config_manager.get_db_config()
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

    def partition_domains(self, pdb_id: str, chain_id: str, dump_dir: str, input_mode: str, 
                         reference: str, blast_only: bool = False
    ) -> Dict[str, Any]:
        """
        Partition domains for a protein chain based on domain summary

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            dump_dir: Base directory for I/O operations
            input_mode: Input mode ('struct_seqid', etc.)
            reference: Reference version
            blast_only: Whether to use BLAST-only summaries

        Returns:
            Dictionary with file path, domain models, and processing information
        """
        from ecod.models.domain import Domain, DomainRange, DomainRangeSegment

        # Initialize result structure
        result = {
            "file_path": None,
            "domains": [],
            "success": False,
            "stats": {
                "domain_count": 0,
                "coverage": 0.0,
                "discontinuous_domains": 0
            },
            "messages": []
        }

        self.logger.info(f"Partitioning domains for {pdb_id}_{chain_id}")

        # Load reference data if not already loaded
        if not self.ref_range_cache:
            self.load_reference_data(reference)



        # Define paths
        pdb_chain = f"{pdb_id}_{chain_id}"

        # Get standardized paths using path_utils
        paths = get_standardized_paths(dump_dir, pdb_id, chain_id, reference, create_dirs=True)

        # Set domain output file path based on standard paths
        file_type = 'blast_only_partition' if blast_only else 'domain_partition'
        domain_fn = paths[file_type]
        result["file_path"] = domain_fn

        force_overwrite = self.context.is_force_overwrite()

        # Check if file already exists
        if os.path.exists(domain_fn) and not force_overwrite:
            self.logger.info(f"Domain file {domain_fn} already exists, skipping...")
            result["success"] = True
            result["messages"].append(f"Domain file already exists: {domain_fn}")

            # ...keep existing code to load domains from file...

            return result

        # Get domain summary path
        summary_type = 'blast_only_summary' if blast_only else 'domain_summary'

        # Try to find domain summary using path_utils (checks standard and legacy paths)
        evidence_paths = get_all_evidence_paths(dump_dir, pdb_id, chain_id, reference)
        if summary_type in evidence_paths and evidence_paths[summary_type]['exists_at']:
            domain_summ_fn = evidence_paths[summary_type]['exists_at']
            self.logger.info(f"Found domain summary at: {domain_summ_fn}")
        else:
            # Fallback to traditional lookup if path_utils doesn't find it
            domain_summ_fn = self._find_domain_summary(pdb_id, chain_id, dump_dir, blast_only)

        if not os.path.exists(domain_summ_fn):
            self.logger.error(f"Domain summary file not found for {pdb_id}_{chain_id}")
            result["messages"].append(f"Domain summary file not found")
            return result

        if not os.path.exists(domain_summ_fn):
            self.logger.error(f"Domain summary file not found: {domain_summ_fn}")
            result["messages"].append(f"Domain summary file not found: {domain_summ_fn}")
            return result

        # Process the summary file
        blast_data = self._process_domain_summary(domain_summ_fn)

        # Check for errors in processing
        if "error" in blast_data:
            self.logger.error(f"Error processing domain summary: {blast_data['error']}")
            result["messages"].append(f"Error processing domain summary: {blast_data['error']}")
            return result

        # Read FASTA sequence
        fasta_path = find_fasta_file(pdb_id, chain_id, dump_dir)
        sequence = self._read_fasta_sequence(fasta_path)
        if not sequence:
            self.logger.error(f"Failed to read sequence from {fasta_path}")
            result["messages"].append(f"Failed to read sequence from {fasta_path}")
            return result

        sequence_length = len(sequence)

        # Determine domain boundaries using model data
        domains = self._determine_domain_boundaries(blast_data, sequence_length, pdb_chain)

        # Check if any domains were found
        if not domains or all(not d.get("evidence", []) for d in domains):
            # Create unclassified document
            domain_doc = self._create_unclassified_document(pdb_id, chain_id, reference, sequence_length)

            # Write output file
            os.makedirs(os.path.dirname(domain_fn), exist_ok=True)
            tree = ET.ElementTree(domain_doc)
            tree.write(domain_fn, encoding='utf-8', xml_declaration=True)

            self.logger.info(f"Created unclassified domain document for {pdb_chain}")
            result["success"] = True
            result["messages"].append(f"Created unclassified domain document: {domain_fn}")
            return result

        # Assign classifications
        self._assign_domain_classifications(domains, blast_data, pdb_chain)

        # Create XML document
        domain_doc, domain_stats = self._create_domain_document(pdb_id, chain_id, reference, domains, sequence_length)

        # Write output file
        os.makedirs(os.path.dirname(domain_fn), exist_ok=True)
        tree = ET.ElementTree(domain_doc)
        tree.write(domain_fn, encoding='utf-8', xml_declaration=True)

        self.logger.info(f"Created domain partition file: {domain_fn}")

        # Create Domain models for result
        domain_models = []
        for domain_dict in domains:
            # Extract domain attributes from dictionary
            domain_id = domain_dict.get("domain_id", "")
            domain_range = domain_dict.get("range", "")

            # Create domain model
            domain = Domain(
                domain_id=domain_id,
                range=domain_range,
                t_group=domain_dict.get("t_group", ""),
                h_group=domain_dict.get("h_group", ""),
                x_group=domain_dict.get("x_group", ""),
                a_group=domain_dict.get("a_group", ""),
                is_manual_rep=domain_dict.get("is_manual_rep", False),
                is_f70=domain_dict.get("is_f70", False),
                is_f40=domain_dict.get("is_f40", False),
                is_f99=domain_dict.get("is_f99", False)
            )

            domain_models.append(domain)

        # Update result with domain models and statistics
        result["domains"] = domain_models
        result["stats"].update(domain_stats)
        result["success"] = True
        result["messages"].append(f"Created domain partition with {len(domain_models)} domains")

        return result

    def process_domains(self, pdb_id: str, chain_id: str,
                       domain_summary_path: str,
                       output_dir: str,
                       reference: str = "develop291"
    ) -> DomainPartitionResult:
        """
        Process domains for a protein chain using domain summary information

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            domain_summary_path: Path to domain summary file
            output_dir: Output directory for domain file
            reference: Reference version (default: develop291)

        Returns:
            DomainPartitionResult model with processing results
        """
        logger = logging.getLogger(__name__)
        try:
            logger.info(f"Starting domain classification assignment for {pdb_id}_{chain_id}")

            # Create domain partition result model
            result = DomainPartitionResult(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference=reference
            )

            # Ensure domains directory exists
            domains_dir = os.path.join(output_dir, "domains")
            os.makedirs(domains_dir, exist_ok=True)

            # Set output file path
            domain_file = os.path.join(
                domains_dir,
                f"{pdb_id}_{chain_id}.{reference}.domains.xml"
            )
            result.domain_file = domain_file

            # Process domains using existing implementation
            domains = self._process_domains_internal(
                pdb_id, chain_id, domain_summary_path, domain_file, reference
            )

            # Handle result
            if not domains or isinstance(domains, str) and domains.startswith("ERROR:"):
                if isinstance(domains, str):
                    result.error = domains
                else:
                    result.error = "No domains found or processing failed"
                result.success = False
                result.is_unclassified = True
                return result

            # Add domains to result
            result.domains = domains if isinstance(domains, list) else []
            result.is_classified = len(result.domains) > 0
            result.is_unclassified = len(result.domains) == 0

            # Save to file
            if result.is_classified:
                if result.save():
                    logger.info(f"Created domain partition file: {domain_file}")
                else:
                    logger.error(f"Failed to save domain partition file: {domain_file}")
                    result.error = "Failed to save domain partition file"
                    result.success = False
            else:
                # Create unclassified domain document
                logger.info(f"Created unclassified domain document for {pdb_id}_{chain_id}")
                self._create_unclassified_document(pdb_id, chain_id, domain_file, reference)
                result.save()

            return result

        except Exception as e:
            logger.error(f"Error processing domains for {pdb_id}_{chain_id}: {str(e)}")
            result = DomainPartitionResult(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference=reference,
                success=False,
                error=str(e)
            )
            return result

    def _process_domains_internal(self, pdb_id, chain_id, domain_summary_path, domain_file, reference):
        """Internal implementation of domain processing - existing code"""
        # This is where your current domain processing logic would go
        # Just return the domains list instead of creating XML directly
        # ...

    def _create_unclassified_document(self, pdb_id, chain_id, domain_file, reference):
        """Create an unclassified domain document"""
        # Create simple XML for unclassified domain
        root = ET.Element("domain_partition")
        root.set("pdb_id", pdb_id)
        root.set("chain_id", chain_id)
        root.set("reference", reference)
        root.set("is_unclassified", "true")

        # Create domains element (empty)
        domains_elem = ET.SubElement(root, "domains")

        # Write to file
        tree = ET.ElementTree(root)
        tree.write(domain_file, encoding="utf-8", xml_declaration=True)

    def _create_domain_document(self, pdb_id: str, chain_id: str, reference: str, domains: List[Dict[str, Any]], sequence_length: int) -> Tuple[ET.Element, Dict[str, Any]]:
        """Create XML document from domain data and return statistics"""
        domains_doc = ET.Element("domain_doc")
        domains_doc.set("pdb", self._safe_str(pdb_id))
        domains_doc.set("chain", self._safe_str(chain_id))
        domains_doc.set("reference", self._safe_str(reference))

        # Create domain list
        domain_list = ET.SubElement(domains_doc, "domain_list")

        # Track statistics for chain coverage
        used_residues = set()
        discontinuous_count = 0

        # Process each domain
        for domain_dict in domains:
            self._add_domain_to_document(domain_dict, domain_list, used_residues, discontinuous_count)

        # Add statistics section
        statistics_elem = ET.SubElement(domains_doc, "statistics")

        # Calculate coverage
        unused_residues = sequence_length - len(used_residues)
        coverage = len(used_residues) / sequence_length if sequence_length > 0 else 0

        coverage_elem = ET.SubElement(statistics_elem, "coverage")
        coverage_elem.set("used_res", str(len(used_residues)))
        coverage_elem.set("unused_res", str(unused_residues))
        coverage_elem.set("total_res", str(sequence_length))
        coverage_elem.text = f"{coverage:.6f}"

        # Add discontinuous domains info
        disc_elem = ET.SubElement(statistics_elem, "discontinuous_domains")
        disc_elem.set("count", str(discontinuous_count))

        # Collect statistics
        stats = {
            "domain_count": len(domains),
            "coverage": coverage,
            "discontinuous_domains": discontinuous_count,
            "used_residues": len(used_residues),
            "unused_residues": unused_residues
        }

        return domains_doc, stats

    def _add_domain_to_document(self, domain: Dict[str, Any], domain_list: ET.Element, used_residues: Set[int], discontinuous_count: int) -> None:
        """Add a domain to the XML document"""
        try:
            # Create domain element
            domain_elem = ET.SubElement(domain_list, "domain")
            domain_elem.set("domain_id", domain.get("domain_id", ""))
            domain_elem.set("pdb", domain.get("pdb", ""))
            domain_elem.set("chain", domain.get("chain", ""))

            # Check for required range attribute
            if "range" not in domain or not domain["range"]:
                self.logger.error(f"Domain missing required range attribute")
                return

            domain_elem.set("range", domain["range"])

            # Add high_confidence flag if this domain was from the alternate method
            if domain.get("protected", False):
                domain_elem.set("high_confidence", "true")

            # Add classification attributes if present
            for attr in ["t_group", "h_group", "x_group", "a_group"]:
                if attr in domain and domain[attr]:
                    domain_elem.set(attr, str(domain[attr]))

            # Add flags if present
            for flag in ["is_manual_rep", "is_f70", "is_f40", "is_f99"]:
                if flag in domain and domain[flag]:
                    domain_elem.set(flag, "true")

            # Add range element
            range_elem = ET.SubElement(domain_elem, "range")
            range_elem.text = domain["range"]

            # Check if this is a discontinuous domain
            if "," in domain["range"]:
                discontinuous_count += 1

                # Add ungapped range with gap tolerance
                if "ungapped_range" in domain:
                    ungapped_elem = ET.SubElement(domain_elem, "ungapped_range")
                    ungapped_elem.text = domain["ungapped_range"]
                    ungapped_elem.set("gap_tolerance", str(domain.get("gap_tolerance", 20)))
                else:
                    # Create ungapped range from segments
                    segments = domain["range"].split(",")
                    if segments:
                        try:
                            first_start = int(segments[0].split("-")[0])
                            last_end = int(segments[-1].split("-")[-1])
                            ungapped_elem = ET.SubElement(domain_elem, "ungapped_range")
                            ungapped_elem.text = f"{first_start}-{last_end}"
                            ungapped_elem.set("gap_tolerance", "20")
                        except (ValueError, IndexError):
                            pass

            # Track used residues for coverage calculation
            domain_ranges = domain["range"].split(",")
            for segment in domain_ranges:
                if "-" in segment:
                    try:
                        start, end = map(int, segment.split("-"))
                        used_residues.update(range(start, end + 1))
                    except (ValueError, IndexError):
                        pass

            # Add classification evidence
            if "evidence" in domain and domain["evidence"]:
                evidence = self._summarize_domain_evidence(domain["evidence"], domain.get("range", ""))
                self._add_evidence_to_domain(evidence, domain_elem)

        except Exception as e:
            self.logger.error(f"Error adding domain to document: {e}")

    def _add_evidence_to_domain(self, evidence: List[Dict[str, Any]], domain_elem: ET.Element) -> None:
        """Add evidence items to a domain element"""
        if not evidence:
            return

        evidence_elem = ET.SubElement(domain_elem, "evidence")

        for e in evidence:
            try:
                match_elem = ET.SubElement(evidence_elem, "match")
                match_elem.set("domain_id", str(e.get("domain_id", "")))
                match_elem.set("type", str(e.get("type", "")))

                # Add T-group if available
                if "t_group" in e and e["t_group"]:
                    match_elem.set("t_group", str(e["t_group"]))

                if "evalue" in e and e["evalue"] is not None:
                    match_elem.set("evalue", str(e["evalue"]))
                if "probability" in e and e["probability"] is not None:
                    match_elem.set("probability", str(e["probability"]))

                if "segment_idx" in e:
                    match_elem.set("segment_idx", str(e.get("segment_idx", 0)))
                if "segment_range" in e:
                    match_elem.set("segment_range", str(e.get("segment_range", "")))

                query_range = ET.SubElement(match_elem, "query_range")
                query_range.text = str(e.get("query_range", ""))

                # Ensure hit range is populated
                hit_range = ET.SubElement(match_elem, "hit_range")
                hit_range_text = str(e.get("hit_range", ""))
                if not hit_range_text and "domain_id" in e:
                    # Try to derive hit range from reference domains
                    domain_id = e.get("domain_id", "")
                    hit_range_text = self._get_domain_range_by_id(domain_id)
                hit_range.text = hit_range_text
            except Exception as e_err:
                self.logger.error(f"Error adding evidence item: {e_err}")

    def _find_domain_summary(self, pdb_id: str, chain_id: str, dump_dir: str, blast_only: bool = False) -> str:
        """Find domain summary file using path utilities to check standard and legacy paths

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            dump_dir: Base directory for I/O operations
            blast_only: Whether to look for BLAST-only summary

        Returns:
            Path to domain summary file if found, empty string otherwise
        """
        from ecod.utils.path_utils import get_all_evidence_paths, resolve_file_path

        # Get current reference version (could also be passed as a parameter)
        reference = self.context.config_manager.config.get('reference', {}).get('current_version', 'develop291')

        # Use path_utils to get standard and legacy paths
        file_type = 'blast_only_summary' if blast_only else 'domain_summary'

        try:
            # Use path_utils to check all possible paths
            evidence_paths = get_all_evidence_paths(dump_dir, pdb_id, chain_id, reference)

            # Check if we found the domain summary
            if file_type in evidence_paths and evidence_paths[file_type]['exists_at']:
                summary_path = evidence_paths[file_type]['exists_at']
                self.logger.info(f"Found domain summary using path_utils: {summary_path}")
                return summary_path

            self.logger.debug(f"Domain summary not found via path_utils for {pdb_id}_{chain_id}")
        except Exception as e:
            self.logger.warning(f"Error using path_utils to find domain summary: {e}")

        # Fallback to database lookup if path_utils approach didn't work
        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)
        query = """
        SELECT pf.file_path
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        AND pf.file_type = %s
        AND pf.file_exists = TRUE
        LIMIT 1
        """

        # Set the file type for db query
        db_file_type = 'blast_only_summary' if blast_only else 'domain_summary'

        try:
            rows = db.execute_query(query, (pdb_id, chain_id, db_file_type))
            if rows:
                db_summ_path = rows[0][0]
                full_summ_path = resolve_file_path(dump_dir, db_summ_path)
                if os.path.exists(full_summ_path):
                    self.logger.info(f"Found domain summary in database: {full_summ_path}")
                    return full_summ_path
        except Exception as e:
            self.logger.warning(f"Error querying database for domain summary: {e}")

        # If all else fails, return empty string
        self.logger.warning(f"Domain summary not found for {pdb_id}_{chain_id}")
        return ""

    def _get_domain_range_by_id(self, domain_id: str) -> str:
        """Get the range string for a domain by its ID"""
        if not domain_id:
            return ""
            
        # Check if we have this domain in reference data
        for source_id, domains in self.ref_range_cache.items():
            for domain in domains:
                if domain.get("domain_id") == domain_id:
                    return domain.get("range", "")
        
        # If not found in cache, try database lookup
        try:
            db_config = self.context.config_manager.get_db_config()
            db = DBManager(db_config)
            query = """
            SELECT range FROM pdb_analysis.domain WHERE domain_id = %s LIMIT 1
            """
            rows = db.execute_query(query, (domain_id,))
            if rows and rows[0][0]:
                return rows[0][0]
        except Exception as e:
            self.logger.error(f"Error getting range for domain {domain_id}: {e}")
            
        return ""

    def _process_domain_summary(self, domain_summary_fn: str) -> Dict[str, Any]:
        """Process domain summary XML file and extract all relevant data"""
        logger = logging.getLogger("ecod.domain_partition")
        logger.info(f"Processing domain summary: {domain_summary_fn}")
        
        # Check file existence
        if not os.path.exists(domain_summary_fn):
            logger.error(f"Domain summary file not found: {domain_summary_fn}")
            return {"error": "File not found"}
        
        try:
            # Parse XML
            tree = ET.parse(domain_summary_fn)
            root = tree.getroot()
            
            # Log root tag and structure for debugging
            logger.debug(f"XML root tag: {root.tag}, attributes: {root.attrib}")
            
            # Initialize result structure
            result = {
                "sequence": "",
                "sequence_length": 0,
                "chain_id": "",
                "blast_hits": [],
                "domain_blast_hits": [],
                "hhsearch_hits": []
            }
            
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
            
            # Process chain-level BLAST hits using model
            chain_blast_hits = []
            for hit_elem in root.findall(".//chain_blast_run/hits/hit"):
                hit = BlastHit.from_xml(hit_elem)
                hit.hit_type = "chain_blast"
                chain_blast_hits.append(hit.to_dict())
            
            result["blast_hits"] = chain_blast_hits
            logger.debug(f"Found {len(chain_blast_hits)} chain BLAST hits")
            
            # Process domain-level BLAST hits using model
            domain_blast_hits = []
            for hit_elem in root.findall(".//blast_run/hits/hit"):
                hit = BlastHit.from_xml(hit_elem)
                hit.hit_type = "domain_blast"
                domain_blast_hits.append(hit.to_dict())
            
            result["domain_blast_hits"] = domain_blast_hits
            logger.debug(f"Found {len(domain_blast_hits)} domain BLAST hits")
            
            # Process HHSearch hits using model
            hhsearch_hits = []
            for hit_elem in root.findall(".//hh_run/hits/hit"):
                hit = HHSearchHit.from_xml(hit_elem)
                hhsearch_hits.append(hit.to_dict())
            
            result["hhsearch_hits"] = hhsearch_hits
            logger.debug(f"Found {len(hhsearch_hits)} HHSearch hits")
            
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
    
    def _analyze_domain_blast_hits(self, domain_blast_hits):
        """
        Analyze domain BLAST hits to identify multiple domains
        
        Args:
            domain_blast_hits: List of domain BLAST hits
            
        Returns:
            List of domain candidates with proper structure for _determine_domain_boundaries
        """
        self.logger.debug(f"Domain BLAST hits count: {len(domain_blast_hits) if domain_blast_hits else 0}")
        
        if not domain_blast_hits:
            return []
        
        # First, sort hits by e-value to prioritize the most significant hits
        sorted_hits = sorted(domain_blast_hits, key=lambda h: h.get('evalue', 999.0))
        
        # Log top hits for debugging
        for i, hit in enumerate(sorted_hits[:10]):  # Show first 10 hits
            if 'range' in hit:
                range_str = hit['range']
            elif 'query_regions' in hit:
                range_str = hit['query_regions']
            elif 'query_reg' in hit:
                range_str = hit['query_reg']
            else:
                range_str = "unknown"
                
            self.logger.debug(f"Top hit #{i+1}: domain_id={hit.get('domain_id', 'unknown')}, "
                           f"range={range_str}, evalue={hit.get('evalue', 'unknown')}")
        
        # Group hits by region (bin by roughly 50 residue windows)
        region_groups = {}
        
        for hit in sorted_hits:
            if 'range' not in hit:
                # Try alternate keys that might contain the range
                if 'query_regions' in hit:
                    range_str = hit['query_regions']
                elif 'query_reg' in hit:
                    range_str = hit['query_reg']
                else:
                    continue
            else:
                range_str = hit['range']
            
            if not range_str:
                continue
                
            # Process each range segment
            ranges = range_str.split(',')
            for r_str in ranges:
                if '-' not in r_str:
                    continue
                    
                try:
                    start, end = map(int, r_str.split('-'))
                    
                    # Create a bin key based on region center
                    bin_center = (start + end) // 2
                    bin_key = bin_center // 50
                    
                    if bin_key not in region_groups:
                        region_groups[bin_key] = []
                    
                    region_groups[bin_key].append({
                        'start': start,
                        'end': end,
                        'domain_id': hit.get('domain_id', ''),
                        'evalue': hit.get('evalue', 999.0),
                        'hit': hit
                    })
                except (ValueError, TypeError):
                    continue
        
        # Log each region group
        self.logger.debug(f"Found {len(region_groups)} region groups from domain BLAST hits")
        for bin_key in sorted(region_groups.keys()):
            hits = region_groups[bin_key]
            starts = [h['start'] for h in hits]
            ends = [h['end'] for h in hits]
            self.logger.debug(f"Region group {bin_key}: {len(hits)} hits, spanning approx. {min(starts)}-{max(ends)}")
        
        # Create domain candidates from region groups
        domain_candidates = []
        
        for bin_key, hits in sorted(region_groups.items()):
            if not hits:
                continue
                
            # Calculate consensus boundaries for this region
            starts = [h['start'] for h in hits]
            ends = [h['end'] for h in hits]
            
            # Use median to be more robust against outliers
            starts.sort()
            ends.sort()
            median_start = starts[len(starts) // 2]
            median_end = ends[len(ends) // 2]
            
            # Find the hit with the best e-value in this region
            best_hit = min(hits, key=lambda h: h.get('evalue', 999.0))
            
            # Create domain candidate
            domain_candidate = {
                "start": median_start,
                "end": median_end,
                "size": median_end - median_start + 1,
                "evidence": [{
                    "type": "domain_blast",
                    "domain_id": best_hit.get('domain_id', ''),
                    "query_range": f"{median_start}-{median_end}",
                    "evalue": best_hit.get('evalue', 999.0)
                }]
            }
            
            domain_candidates.append(domain_candidate)
            
            # Log detailed information about this domain candidate
            self.logger.debug(f"Domain candidate from region {bin_key}: {median_start}-{median_end} "
                          f"({median_end - median_start + 1} residues), "
                          f"best evalue: {best_hit.get('evalue', 'unknown')}, "
                          f"domain_id: {best_hit.get('domain_id', 'unknown')}")
        
        # Sort domain candidates by start position for clearer logs
        domain_candidates.sort(key=lambda d: d["start"])
        
        # Log consolidated list
        domains_str = ", ".join([f"{d['start']}-{d['end']}" for d in domain_candidates])
        self.logger.info(f"Found {len(domain_candidates)} domain candidates from domain BLAST hits: {domains_str}")
        
        return domain_candidates

    def _analyze_chainwise_hits_for_domains(self, chain_blast_hits):
        # Keep most of the improved version, but add field name mapping
        
        # Initialize counters for diagnostics (same as before)
        hits_total = len(chain_blast_hits) if chain_blast_hits else 0
        hits_processed = 0
        hits_with_valid_fields = 0
        hits_with_valid_regions = 0
        hits_with_reference_domains = 0
        hits_with_mapped_domains = 0
        domains_mapped_total = 0

        self.logger.info(f"Starting analysis of {hits_total} chain BLAST hits")

        if not chain_blast_hits:
            self.logger.warning("No chain BLAST hits provided to analyze")
            return []

        # Log first hit structure
        if hits_total > 0:
            self.logger.debug(f"First chain hit keys: {sorted(chain_blast_hits[0].keys())}")
        
        domain_candidates = []
        
        for hit_idx, hit in enumerate(chain_blast_hits):
            hits_processed += 1
            
            # Validate hit structure - MODIFIED to check for alternative field names
            missing_fields = []
            if "pdb_id" not in hit:
                missing_fields.append("pdb_id")
            if "chain_id" not in hit:
                missing_fields.append("chain_id")
                
            # Check for either query_regions OR range for query
            has_query_regions = "query_regions" in hit and hit["query_regions"]
            has_range = "range" in hit and hit["range"]
            
            # Check for either hit_regions OR hit_range for hit
            has_hit_regions = "hit_regions" in hit and hit["hit_regions"]
            has_hit_range = "hit_range" in hit and hit["hit_range"]
            
            if not (has_query_regions or has_range):
                missing_fields.append("query_regions/range")
            if not (has_hit_regions or has_hit_range):
                missing_fields.append("hit_regions/hit_range")
            
            if missing_fields:
                self.logger.warning(f"Hit #{hit_idx+1} missing required fields: {', '.join(missing_fields)}")
                continue
            
            hits_with_valid_fields += 1
            
            # Extract PDB and chain IDs
            hit_pdb_id = hit.get("pdb_id", "")
            hit_chain_id = hit.get("chain_id", "")
            source_id = f"{hit_pdb_id}_{hit_chain_id}"

            self.logger.debug(f"Processing hit #{hit_idx+1} for chain: {source_id}")
            
            # Get reference domains for this chain
            reference_domains = self._get_reference_chain_domains(source_id)
            
            if not reference_domains:
                self.logger.warning(f"No reference domains found for chain: {source_id}")
                continue
                
            hits_with_reference_domains += 1
            self.logger.debug(f"Found {len(reference_domains)} reference domains for {source_id}")
            
            # Log detailed information about reference domains
            for i, domain in enumerate(reference_domains[:3]):  # Log first 3
                self.logger.debug(f"Reference domain {i+1}: {domain.get('domain_id', 'unknown')}, "
                               f"range: {domain.get('range', 'unknown')}, "
                               f"t_group: {domain.get('t_group', 'unknown')}")
            
            # Log query and hit regions
            # Extract query and hit region strings - MODIFIED to use either field name
            query_regions_str = ""
            hit_regions_str = ""
            
            if has_query_regions:
                query_regions_str = hit.get("query_regions", "")
            elif has_range:
                query_regions_str = hit.get("range", "")
                
            if has_hit_regions:
                hit_regions_str = hit.get("hit_regions", "")
            elif has_hit_range:
                hit_regions_str = hit.get("hit_range", "")
                
            self.logger.debug(f"Query regions: {query_regions_str}")
            self.logger.debug(f"Hit regions: {hit_regions_str}")
            
            # Parse query and hit regions from alignment
            query_regions = []
            hit_regions = []

            # Check both naming conventions
            has_query_data = False
            has_hit_data = False

            if ("query_regions" in hit and hit["query_regions"]) or ("range" in hit and hit["range"]):
                has_query_data = True
                
            if ("hit_regions" in hit and hit["hit_regions"]) or ("hit_range" in hit and hit["hit_range"]):
                has_hit_data = True

            if has_query_data and has_hit_data:
                # Now use the strings we already extracted
                # This maintains the validation logic while using the extracted data
                query_region_strs = query_regions_str.split(",")
                hit_region_strs = hit_regions_str.split(",")
                
                if len(query_region_strs) != len(hit_region_strs):
                    self.logger.warning(f"Mismatched region counts in hit #{hit_idx+1}: "
                                      f"{len(query_region_strs)} query regions, "
                                      f"{len(hit_region_strs)} hit regions")
                
                for region_idx, (q_range, h_range) in enumerate(zip(query_region_strs, hit_region_strs)):
                    if not q_range or not h_range:
                        self.logger.warning(f"Empty region at index {region_idx} for hit #{hit_idx+1}")
                        continue
                    
                    try:
                        q_start, q_end = map(int, q_range.split("-"))
                        h_start, h_end = map(int, h_range.split("-"))
                        
                        query_regions.append((q_start, q_end))
                        hit_regions.append((h_start, h_end))
                    except (ValueError, AttributeError) as e:
                        self.logger.warning(f"Failed to parse region at index {region_idx} for hit #{hit_idx+1}: {e}")
                        self.logger.warning(f"  Query range: '{q_range}', Hit range: '{h_range}'")
                        continue
            else:
                self.logger.warning(f"Hit #{hit_idx+1} missing required region data for {source_id}")
                continue
            
            # Skip if no alignment regions were parsed
            if not query_regions or not hit_regions:
                self.logger.warning(f"Skipping hit #{hit_idx+1} for {source_id}: No valid alignment regions parsed")
                continue
            
            hits_with_valid_regions += 1
            self.logger.debug(f"Successfully parsed {len(query_regions)} alignment regions")
            
            # Track if any domains were mapped for this hit
            hit_mapped_domains = 0
            
            # Map reference domains to query coordinates
            for domain_idx, domain in enumerate(reference_domains):
                domain_id = domain.get("domain_id", "unknown")
                domain_range = domain.get("range", "")
                
                if not domain_range:
                    self.logger.warning(f"Domain #{domain_idx+1} ({domain_id}) has empty range")
                    continue
                    
                domain_ranges = self._parse_range(domain_range)
                
                # Skip domains without a valid range
                if not domain_ranges:
                    self.logger.warning(f"Failed to parse range '{domain_range}' for domain {domain_id}")
                    continue
                    
                self.logger.debug(f"Attempting to map domain {domain_id} with range {domain_range}")
                self.logger.debug(f"Parsed domain ranges: {domain_ranges}")
                
                # Find corresponding query positions
                mapped_ranges = []
                
                for d_range_idx, (d_start, d_end) in enumerate(domain_ranges):
                    segment_mapped = False
                    
                    for hit_range_idx, (h_start, h_end) in enumerate(hit_regions):
                        # Check for overlap
                        if h_end < d_start or h_start > d_end:
                            continue
                            
                        # Calculate overlap region
                        overlap_start = max(d_start, h_start)
                        overlap_end = min(d_end, h_end)
                        
                        if overlap_end >= overlap_start:
                            # Log overlap details
                            self.logger.debug(f"Found overlap between domain segment {d_start}-{d_end} "
                                            f"and hit region {h_start}-{h_end}: "
                                            f"overlap {overlap_start}-{overlap_end}")
                            
                            # Map to query coordinates
                            q_start, q_end = query_regions[hit_range_idx]
                            
                            # Calculate proportion
                            h_len = h_end - h_start
                            if h_len == 0:
                                self.logger.warning(f"Zero-length hit region at index {hit_range_idx}")
                                continue
                                
                            q_len = q_end - q_start
                            
                            # Map start position
                            start_ratio = (overlap_start - h_start) / h_len
                            q_mapped_start = int(q_start + start_ratio * q_len)
                            
                            # Map end position
                            end_ratio = (overlap_end - h_start) / h_len
                            q_mapped_end = int(q_start + end_ratio * q_len)
                            
                            # Validate mapped positions
                            if q_mapped_end < q_mapped_start:
                                self.logger.warning(f"Invalid mapped range: {q_mapped_start}-{q_mapped_end} "
                                                  f"(end < start)")
                                continue
                                
                            mapped_ranges.append((q_mapped_start, q_mapped_end))
                            
                            self.logger.debug(f"Mapped domain segment to query: {q_mapped_start}-{q_mapped_end}")
                            segment_mapped = True
                    
                    if not segment_mapped:
                        self.logger.debug(f"Could not map domain segment {d_start}-{d_end} to any hit region")
                
                # Only add domain if we could map it
                if not mapped_ranges:
                    self.logger.warning(f"Failed to map any segments of domain {domain_id}")
                    continue
                    
                # Calculate domain start and end
                # Use the earliest start and latest end from all mapped ranges
                start = min(r[0] for r in mapped_ranges)
                end = max(r[1] for r in mapped_ranges)
                size = end - start + 1
                
                # Validate size
                if size < 20:  # Minimum domain size threshold
                    self.logger.warning(f"Mapped domain too small ({size} residues), skipping: {start}-{end}")
                    continue
                    
                # Create domain candidate with structure expected by _determine_domain_boundaries
                domain_candidate = {
                    "start": start,
                    "end": end,
                    "size": size,
                    "t_group": domain.get("t_group", ""),
                    "h_group": domain.get("h_group", ""),
                    "x_group": domain.get("x_group", ""),
                    "a_group": domain.get("a_group", ""),
                    "evidence": [{
                        "type": "chain_blast",
                        "domain_id": domain_id,
                        "query_range": f"{start}-{end}",
                        "hit_range": domain_range
                    }]
                }
                
                domain_candidates.append(domain_candidate)
                hit_mapped_domains += 1
                domains_mapped_total += 1
                
                self.logger.info(f"Successfully mapped domain {domain_id} from {source_id}: "
                              f"{start}-{end} ({size} residues)")
            
            if hit_mapped_domains > 0:
                hits_with_mapped_domains += 1
                self.logger.info(f"Mapped {hit_mapped_domains} domains from hit #{hit_idx+1}")
            else:
                self.logger.warning(f"No domains could be mapped from hit #{hit_idx+1}")
        
        # Log summary statistics
        self.logger.info(f"Chain BLAST hit analysis summary:")
        self.logger.info(f"  Total hits: {hits_total}")
        self.logger.info(f"  Hits processed: {hits_processed}")
        self.logger.info(f"  Hits with valid fields: {hits_with_valid_fields}")
        self.logger.info(f"  Hits with valid regions: {hits_with_valid_regions}")
        self.logger.info(f"  Hits with reference domains: {hits_with_reference_domains}")
        self.logger.info(f"  Hits with mapped domains: {hits_with_mapped_domains}")
        self.logger.info(f"  Total domains mapped: {domains_mapped_total}")
        
        # Log domain candidate distribution
        if domain_candidates:
            self.logger.info(f"Domain size distribution:")
            
            # Group by size ranges
            size_ranges = [(0, 50), (51, 100), (101, 200), (201, 500), (501, float('inf'))]
            size_counts = {f"{r[0]}-{r[1] if r[1] != float('inf') else 'inf'}": 0 for r in size_ranges}
            
            for candidate in domain_candidates:
                size = candidate["size"]
                for r_start, r_end in size_ranges:
                    if r_start <= size <= r_end:
                        range_key = f"{r_start}-{r_end if r_end != float('inf') else 'inf'}"
                        size_counts[range_key] += 1
                        break
            
            for range_key, count in size_counts.items():
                if count > 0:
                    self.logger.info(f"  Size {range_key}: {count} domains")
        
        return domain_candidates

    def _determine_domain_boundaries(self, blast_data, sequence_length, pdb_chain):
        """
        Determine final domain boundaries using all available evidence
        """
        logger = logging.getLogger("ecod.domain_partition")
        logger.info(f"Determining domain boundaries for {pdb_chain} (length: {sequence_length})")

        # First, identify highly confident domain hits that should be preserved
        high_confidence_domains = self._respect_high_scoring_hits(blast_data, sequence_length, pdb_chain)
        
        # If we have high-confidence domains, they take precedence
        if high_confidence_domains:
            logger.info(f"Using {len(high_confidence_domains)} high-confidence domain hits as domain boundaries")
            
            # Check if high-confidence domains cover most of the sequence
            coverage = set()
            for domain in high_confidence_domains:
                for i in range(domain["start"], domain["end"] + 1):
                    coverage.add(i)
                    
            coverage_pct = len(coverage) / sequence_length if sequence_length > 0 else 0
            
            # If high-confidence domains cover at least 70% of sequence, use them directly
            if coverage_pct >= 0.7:
                logger.info(f"High-confidence domains cover {coverage_pct:.1%} of sequence, using them directly")
                
                # Check for overlaps and resolve if needed
                final_domains = self._resolve_domain_overlaps(high_confidence_domains)
                
                # Sort by position
                final_domains.sort(key=lambda d: d["start"])
                
                # Ensure domain IDs are set
                for i, domain in enumerate(final_domains):
                    if "domain_id" not in domain:
                        domain["domain_id"] = f"e{pdb_chain}{i+1}"
                        
                    # Convert range format
                    domain["range"] = f"{domain['start']}-{domain['end']}"
                
                return final_domains
        
        # Get domains from different sources
        #blast_domains = self._identify_domains_from_blast(blast_data.get("blast_hits", []), sequence_length)
        hhsearch_domains = self._identify_domains_from_hhsearch(blast_data.get("hhsearch_hits", []), sequence_length)
        chain_domain_candidates = self._analyze_chainwise_hits_for_domains(blast_data.get("blast_hits", []))
        domain_blast_candidates = self._analyze_domain_blast_hits(blast_data.get("domain_blast_hits", []))

        # Log domain counts
        logger.info(f"Domain candidates:  {len(hhsearch_domains)} from HHSearch, " 
            f"{len(chain_domain_candidates)} from chain BLAST, {len(domain_blast_candidates)} from domain BLAST")
        
        # Consolidate domains from different sources - excluding regular BLAST domains
        all_domains = []
        #all_domains.extend([{**d, "source": "blast"} for d in blast_domains])
        all_domains.extend([{**d, "source": "hhsearch"} for d in hhsearch_domains])
        all_domains.extend([{**d, "source": "chain_blast_domain"} for d in chain_domain_candidates])
        all_domains.extend([{**d, "source": "domain_blast"} for d in domain_blast_candidates])
        
        # If no domains found, use whole chain
        if not all_domains:
            logger.warning(f"No domain evidence found for {pdb_chain}, marking as unclassified")
        
        # Track which domains cover each position
        position_domain_coverage = [[] for _ in range(sequence_length + 1)]  # 1-indexed
        
        for domain_idx, domain in enumerate(all_domains):
            for i in range(max(1, domain["start"]), min(sequence_length + 1, domain["end"] + 1)):
                position_domain_coverage[i-1].append(domain_idx)
        
        # Find regions with consistent domain coverage
        regions = []
        current_region_start = None
        current_covering_domains = set()
        
        for i in range(sequence_length):
            covering_domains = set(position_domain_coverage[i])
            
            if covering_domains:
                if current_region_start is None:
                    # Start new region
                    current_region_start = i + 1
                    current_covering_domains = covering_domains
                elif covering_domains != current_covering_domains:
                    # Domain coverage changed, finalize current region and start new one
                    regions.append({
                        "start": current_region_start,
                        "end": i,
                        "size": i - current_region_start + 1,
                        "covering_domains": current_covering_domains,
                        "domain_count": len(current_covering_domains)
                    })
                    current_region_start = i + 1
                    current_covering_domains = covering_domains
            else:
                if current_region_start is not None:
                    # End of covered region
                    regions.append({
                        "start": current_region_start,
                        "end": i,
                        "size": i - current_region_start + 1,
                        "covering_domains": current_covering_domains,
                        "domain_count": len(current_covering_domains)
                    })
                    current_region_start = None
                    current_covering_domains = set()
        
        # Add final region if exists
        if current_region_start is not None:
            regions.append({
                "start": current_region_start,
                "end": sequence_length,
                "size": sequence_length - current_region_start + 1,
                "covering_domains": current_covering_domains,
                "domain_count": len(current_covering_domains)
            })
        
        # Log regions
        logger.debug(f"Found {len(regions)} regions with consistent domain coverage")
        for i, region in enumerate(regions):
            logger.debug(f"Region {i+1}: {region['start']}-{region['end']} ({region['size']} residues), "
                       f"covered by {region['domain_count']} domains")
        
        # Group regions by their covering domains
        domain_region_groups = {}
        for region in regions:
            domain_set = frozenset(region["covering_domains"])
            if domain_set not in domain_region_groups:
                domain_region_groups[domain_set] = []
            domain_region_groups[domain_set].append(region)
        
        # Combine adjacent regions in each group
        consolidated_regions = {}
        for domain_set, regions_list in domain_region_groups.items():
            # Sort regions by start position
            regions_list.sort(key=lambda r: r["start"])
            
            # Combine adjacent regions with small gaps
            gap_threshold = 20
            merged_regions = []
            current_region = None
            
            for region in regions_list:
                if current_region is None:
                    current_region = region.copy()
                elif region["start"] <= current_region["end"] + gap_threshold:
                    # Merge regions
                    current_region["end"] = region["end"]
                    current_region["size"] = current_region["end"] - current_region["start"] + 1
                else:
                    # Start new region
                    merged_regions.append(current_region)
                    current_region = region.copy()
            
            if current_region:
                merged_regions.append(current_region)
            
            consolidated_regions[domain_set] = merged_regions
        
        # Score each domain set based on:
        # 1. Coverage (total residues covered)
        # 2. Evidence quality (prioritize domain BLAST)
        # 3. Domain size (prefer larger domains)
        domain_set_scores = {}
        
        for domain_set, regions_list in consolidated_regions.items():
            total_coverage = sum(region["size"] for region in regions_list)
            
            # Calculate evidence quality score
            evidence_quality = 0
            domain_blast_count = 0
            chain_blast_count = 0
            hhsearch_count = 0
            
            for domain_idx in domain_set:
                source = all_domains[domain_idx]["source"]
                if source == "domain_blast":
                    domain_blast_count += 1
                    evidence_quality += 3  # Higher weight for domain BLAST
                elif source == "chain_blast_domain":
                    chain_blast_count += 1
                    evidence_quality += 2
                elif source == "hhsearch":
                    hhsearch_count += 1
                    evidence_quality += 1
            
            # Calculate size score (prefer regions > 50 residues)
            size_score = sum(1 for region in regions_list if region["size"] >= 50)
            
            # Combined score
            domain_set_scores[domain_set] = {
                "total_coverage": total_coverage,
                "evidence_quality": evidence_quality,
                "size_score": size_score,
                "combined_score": total_coverage * (1 + evidence_quality) * (1 + size_score),
                "domain_blast_count": domain_blast_count,
                "domain_indices": list(domain_set)
            }
        
        # Sort domain sets by score
        sorted_domain_sets = sorted(
            domain_set_scores.items(),
            key=lambda x: x[1]["combined_score"],
            reverse=True
        )
        
        # Log scoring results
        logger.debug(f"Scored {len(sorted_domain_sets)} domain sets")
        for i, (domain_set, score) in enumerate(sorted_domain_sets[:5]):  # Log top 5
            logger.debug(f"Domain set {i+1}: score={score['combined_score']:.1f}, "
                       f"coverage={score['total_coverage']}, "
                       f"quality={score['evidence_quality']}, "
                       f"size_score={score['size_score']}, "
                       f"domain_blast_count={score['domain_blast_count']}")
        
        # Create final domains
        final_domains = []
        assigned_positions = [False] * (sequence_length + 1)  # Track assigned positions
        min_domain_size = 30  # Minimum domain size
        
        # Process domain sets in score order
        for domain_set, score in sorted_domain_sets:
            # Skip if no domain BLAST evidence in this set
            if score["domain_blast_count"] == 0 and len(final_domains) > 0:
                continue
                
            # Get regions for this domain set
            regions_list = consolidated_regions[domain_set]
            
            for region in regions_list:
                # Skip if too small
                if region["size"] < min_domain_size:
                    continue
                    
                # Calculate overlap with existing domains
                overlap_count = sum(1 for i in range(region["start"], region["end"] + 1) 
                                  if i <= sequence_length and assigned_positions[i-1])
                overlap_percentage = overlap_count / region["size"]
                
                # Skip if too much overlap (>5%)
                if overlap_percentage > 0.05 and len(final_domains) > 0:
                    logger.debug(f"Skipping region {region['start']}-{region['end']} due to {overlap_percentage:.1%} overlap")
                    continue
                
                # Determine source
                sources = []
                evidence_items = []
                
                for domain_idx in domain_set:
                    domain = all_domains[domain_idx]
                    sources.append(domain["source"])
                    if "evidence" in domain:
                        evidence_items.extend(domain["evidence"])
                
                source_str = "+".join(set(sources))
                
                # Mark positions as assigned
                for i in range(region["start"], region["end"] + 1):
                    if i <= sequence_length:
                        assigned_positions[i-1] = True
                
                # Create domain
                final_domains.append({
                    "domain_num": len(final_domains) + 1,
                    "range": f"{region['start']}-{region['end']}",
                    "size": region["size"],
                    "confidence": "high" if "domain_blast" in sources else "medium",
                    "source": source_str,
                    "reason": f"Domain supported by {len(set(sources))} evidence types: {source_str}",
                    "evidence": evidence_items
                })
        
        # If no domains were created, fallback to whole chain
        if not final_domains:
            logger.info(f"No domains found for {pdb_chain}, returning empty domain list")
            return []

        # Sort domains by position
        final_domains.sort(key=lambda d: int(d["range"].split("-")[0]))
        
        # Log final domains
        logger.info(f"Final domain boundaries for {pdb_chain} after filtering: {len(final_domains)} domains")
        for domain in final_domains:
            logger.debug(f"Domain {domain['domain_num']}: {domain['range']} ({domain['size']} residues), "
                       f"confidence: {domain['confidence']}, source: {domain['source']}")
        
        return final_domains

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

        
        # Filter out any domains that have no evidence
        final_domains = [d for d in final_domains if d.get("evidence", [])]
        
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
        self.logger.debug(f"  Chain BLAST hits: {len(blast_data.get('chain_blast_hits', []))}")
        self.logger.debug(f"  Domain BLAST hits: {len(blast_data.get('domain_blast_hits', []))}")
        self.logger.debug(f"  HHSearch hits: {len(blast_data.get('hhsearch_hits', []))}")
        
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
        
        # Assign classifications to domains
        for i, domain in enumerate(domains):
            self.logger.debug(f"Assigning classification to domain {i+1}: {domain.get('range', 'unknown_range')}")
            
            # If domain has reference, use it directly
            if "reference" in domain and domain["reference"]:
                domain_id = domain.get("domain_id", "")
                if domain_id in reference_classifications:
                    domain.update(reference_classifications[domain_id])
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
                elif evidence["type"] == "domain_blast":
                    # NEW: Handle domain_blast evidence type properly
                    e_value = evidence.get("evalue", 999)
                    # Convert e-value to score (lower e-value = higher score)
                    score = 100.0 / (1.0 + e_value) if e_value < 10 else 0
                    # Apply a bonus for domain blast hits because they're more specific
                    score *= 1.5  # Give domain_blast evidence a 50% bonus
                    self.logger.debug(f"  Domain BLAST evidence - evalue: {e_value}, score: {score}")
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

                    for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
                        if cls_attr in classification and classification[cls_attr]:
                            best_evidence[cls_attr] = classification[cls_attr]

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
        
        # Define significance thresholds
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
            if "query_reg" not in hit or not hit["query_reg"].text:
                continue
                
            # Parse range and convert to coordinates
            query_range = hit["query_reg"].text
            range_parsed = self._parse_range(query_range)
            
            if not range_parsed:
                continue
                
            # Check significance (either condition)
            probability = float(hit.get("probability", 0))
            evalue = float(hit.get("evalue", 999))
            
            if (probability >= thresholds["probability"] or
                evalue <= thresholds["evalue"]):
                # Add parsed range to hit
                hit_copy = hit.attrib.copy()
                hit_copy["range_parsed"] = range_parsed
                significant_hits.append(hit_copy)
        
        logger.info(f"Found {len(significant_hits)}/{len(hhsearch_hits)} significant HHSearch hits")
        
        if not significant_hits:
            logger.warning("No significant HHSearch hits found")
            return []
        
        # Analyze coverage with significant hits
        position_coverage = [0] * (sequence_length + 1)  # 1-indexed
        hit_regions = []
        
        for hit in significant_hits:
            hit_region = {
                "domain_id": hit.get("domain_id", "unknown"),
                "probability": float(hit.get("probability", 0)),
                "evalue": float(hit.get("evalue", 999)),
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
        
        # Create domain candidates with structure expected by _determine_domain_boundaries
        domain_candidates = []
        
        for region in merged_regions:
            # Find best hit for this region
            best_hit = None
            best_evidence_quality = 0
            
            for hit in significant_hits:
                # Check overlap
                overlap = 0
                hit_coverage = 0
                
                for start, end in hit["range_parsed"]:
                    for i in range(start, end + 1):
                        if region["start"] <= i <= region["end"]:
                            overlap += 1
                    hit_coverage += end - start + 1
                
                # Calculate overlap percentage
                region_size = region["size"]
                hit_overlap_percentage = overlap / region_size if region_size > 0 else 0
                hit_coverage_percentage = overlap / hit_coverage if hit_coverage > 0 else 0
                
                # Calculate evidence quality (combination of probability and coverage)
                probability = float(hit.get("probability", 0))
                evidence_quality = probability * hit_overlap_percentage * hit_coverage_percentage
                
                if evidence_quality > best_evidence_quality:
                    best_hit = hit
                    best_evidence_quality = evidence_quality
            
            # Create domain candidate
            if best_hit:
                domain_candidate = {
                    "start": region["start"],
                    "end": region["end"],
                    "size": region["size"],
                    "source": "hhsearch",
                    "evidence": [{
                        "type": "hhsearch",
                        "domain_id": best_hit.get("domain_id", ""),
                        "query_range": f"{region['start']}-{region['end']}",
                        "probability": float(best_hit.get("probability", 0)),
                        "evalue": float(best_hit.get("evalue", 999))
                    }]
                }
                domain_candidates.append(domain_candidate)
        
        return domain_candidates

    def _extract_hhsearch_hits(self, blast_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Extract HHSearch hits from blast data
        
        Args:
            blast_data: Blast data dictionary
            
        Returns:
            List of HHSearch hits
        """
        # Check if HHSearch data is present
        hh_run = blast_data.get("blast_summ", {}).find(".//hh_run")
        if hh_run is None:
            return []
        
        # Extract hits
        hits = []
        for hit_elem in hh_run.findall(".//hit"):
            hits.append(hit_elem)
        
        return hits

    def _summarize_domain_evidence(self, evidence_list, domain_range=None):
        """
        Summarizes domain evidence to create a more compact representation,
        handling domains with non-overlapping ranges (discontinuous domains)
        
        Args:
            evidence_list: List of evidence dictionaries
            domain_range: String containing the domain range (e.g. "1-100,150-200")
            
        Returns:
            List of summarized evidence items with only essential fields
        """
        if not evidence_list:
            return []
        
        # Parse domain range if provided to identify domain segments
        domain_segments = []
        if domain_range:
            segments = domain_range.split(",")
            for segment in segments:
                if "-" in segment:
                    try:
                        start, end = map(int, segment.split("-"))
                        domain_segments.append((start, end))
                    except (ValueError, IndexError):
                        pass
        
        # Group evidence by type and by segment
        evidence_by_type_segment = {}
        for evidence in evidence_list:
            ev_type = evidence.get("type", "unknown")
            
            # Determine which segment this evidence corresponds to
            segment_idx = 0  # Default to first segment
            query_range = evidence.get("query_range", "")
            
            if query_range and domain_segments:
                try:
                    # Parse query range
                    if "-" in query_range:
                        q_start, q_end = map(int, query_range.split("-"))
                        
                        # Find matching segment
                        best_overlap = 0
                        for idx, (d_start, d_end) in enumerate(domain_segments):
                            # Calculate overlap
                            overlap_start = max(q_start, d_start)
                            overlap_end = min(q_end, d_end)
                            overlap = max(0, overlap_end - overlap_start + 1)
                            
                            if overlap > best_overlap:
                                best_overlap = overlap
                                segment_idx = idx
                except (ValueError, IndexError):
                    pass
            
            key = (ev_type, segment_idx)
            if key not in evidence_by_type_segment:
                evidence_by_type_segment[key] = []
            evidence_by_type_segment[key].append(evidence)
        
        # Summarize each type+segment combination
        summarized_evidence = []
        
        for (ev_type, segment_idx), items in evidence_by_type_segment.items():
            # Sort by quality (lower e-value or higher probability)
            if ev_type == "hhsearch":
                # For HHSearch, sort by probability (higher is better)
                items.sort(key=lambda x: float(x.get("probability", 0)), reverse=True)
                quality_field = "probability"
                quality_format = lambda x: f"{x:.1f}%"
            else:
                # For BLAST, sort by e-value (lower is better)
                items.sort(key=lambda x: float(x.get("evalue", 999)))
                quality_field = "evalue"
                quality_format = lambda x: f"{x:.2e}"
            
            # Take top 3 items
            top_items = items[:3]
            
            # Generate concise summary
            if len(items) == 1:
                # If only one item, include full details
                summarized_evidence.append(items[0])
            else:
                # Create a summary item
                best_item = top_items[0]
                summary = {
                    "type": ev_type,
                    "domain_id": best_item.get("domain_id", ""),
                    "query_range": best_item.get("query_range", ""),
                    "hit_range": best_item.get("hit_range", ""),
                    quality_field: best_item.get(quality_field, "")
                }
                
                # Add segment information for discontinuous domains
                if domain_segments and len(domain_segments) > 1:
                    summary["segment_idx"] = segment_idx
                    if segment_idx < len(domain_segments):
                        segment = domain_segments[segment_idx]
                        summary["segment_range"] = f"{segment[0]}-{segment[1]}"
                
                # Add summary of additional items
                additional_ids = [item.get("domain_id", "") for item in top_items[1:] 
                                 if item.get("domain_id") != best_item.get("domain_id")]
                
                if additional_ids:
                    summary["additional_matches"] = ", ".join(additional_ids)
                    summary["match_count"] = len(items)
                
                summarized_evidence.append(summary)
        
        # Sort by segment index for discontinuous domains
        if domain_segments and len(domain_segments) > 1:
            summarized_evidence.sort(key=lambda x: x.get("segment_idx", 0))
        
        return summarized_evidence

    def _respect_high_scoring_hits(self, blast_data, sequence_length, pdb_chain):
        """
        Pre-process domain boundaries by identifying high-confidence domain hits
        that should be preserved as single domains
        
        Args:
            blast_data: Dictionary containing BLAST results
            sequence_length: Length of the protein sequence
            pdb_chain: String in format "pdbid_chain" (e.g., "8c9i_A")
            
        Returns:
            List of domain dictionaries for high-confidence hits
        """
        logger = logging.getLogger("ecod.domain_partition")
        
        # Parse pdb_chain to get PDB ID and chain ID
        if "_" in pdb_chain:
            pdb_id, chain_id = pdb_chain.split("_")
        else:
            # Handle case where pdb_chain format is different
            logger.warning(f"Unexpected pdb_chain format: {pdb_chain}")
            pdb_id = pdb_chain[:4] if len(pdb_chain) >= 4 else pdb_chain
            chain_id = pdb_chain[4:5] if len(pdb_chain) >= 5 else "A"
        
        # Extract domain BLAST hits - these should be considered authoritative
        domain_blast_hits = blast_data.get("domain_blast_hits", [])
        
        # Filter for high-confidence hits
        high_confidence_domains = []
        domain_num = 1  # Initialize domain counter
        
        for hit in domain_blast_hits:
            # Check for required attributes
            if "evalue" not in hit or not hit.get("query_range"):
                continue
                
            try:
                evalue = float(hit.get("evalue", 999))
                # Use a stricter threshold for "definitive" hits
                if evalue < 1e-30:  # Very high confidence
                    # Parse query range exactly as provided in evidence
                    query_range = hit.get("query_range", "")
                    if "-" in query_range:
                        try:
                            start, end = map(int, query_range.split("-"))
                            # Check if it's a substantial domain (not a fragment)
                            if end - start + 1 >= 50:  # Minimum size threshold
                                # Generate proper domain ID
                                domain_id = f"e{pdb_id}{chain_id}{domain_num}"
                                domain_num += 1
                                
                                # Create domain with EXACT range from evidence
                                high_confidence_domains.append({
                                    "start": start,
                                    "end": end,
                                    "size": end - start + 1,
                                    "domain_id": domain_id,  # Use proper generated ID
                                    "hit_domain_id": hit.get("domain_id", ""),  # Store original hit ID separately
                                    "evalue": evalue,
                                    "protected": True,  # Mark as protected
                                    "range": f"{start}-{end}",  # Store range string format too
                                    "evidence": [{
                                        "type": "domain_blast",
                                        "domain_id": hit.get("domain_id", ""),  # Original hit domain ID
                                        "query_range": query_range,  # Exact query range
                                        "hit_range": hit.get("hit_range", ""),
                                        "evalue": evalue
                                    }]
                                })
                                #logger.info(f"Found high-confidence domain: {domain_id} (based on {hit.get('domain_id', '')}) "
                                #         f"at {query_range} with e-value {evalue}")
                        except (ValueError, TypeError):
                            continue
            except (ValueError, TypeError):
                continue
        
        # Handle overlapping high-confidence domains
        if len(high_confidence_domains) > 1:
            high_confidence_domains = self._resolve_domain_overlaps(high_confidence_domains)
            
            # Renumber domains after overlap resolution
            for i, domain in enumerate(sorted(high_confidence_domains, 
                                       key=lambda d: d["start"])):
                domain["domain_id"] = f"e{pdb_id}{chain_id}{i+1}"
        
        return high_confidence_domains

    def _resolve_domain_overlaps(self, domains):
        """
        Resolve overlaps between domains, preserving high-confidence domains
        """
        if not domains or len(domains) <= 1:
            return domains
            
        # Sort domains by confidence (lower e-value = higher confidence)
        sorted_domains = sorted(domains, key=lambda d: d.get("evalue", 999))
        
        # Track which positions are assigned
        max_pos = max(d["end"] for d in sorted_domains)
        assigned = [False] * (max_pos + 1)
        
        # Assign positions to domains in order of confidence
        final_domains = []
        
        for domain in sorted_domains:
            # Calculate overlap with already assigned positions
            overlap = 0
            for i in range(domain["start"], domain["end"] + 1):
                if i <= max_pos and assigned[i]:
                    overlap += 1
                    
            # Calculate overlap percentage
            overlap_pct = overlap / (domain["end"] - domain["start"] + 1)
            
            # If minimal overlap (<20%), include this domain
            if overlap_pct < 0.2:
                # Mark positions as assigned
                for i in range(domain["start"], domain["end"] + 1):
                    if i <= max_pos:
                        assigned[i] = True
                        
                final_domains.append(domain)
        
        return final_domains


    def _update_process_status(self, process_id: int, result: DomainPartitionResult) -> None:
        """Update process status in database"""
        # Implementation depends on your database setup
        pass

    def _find_domain_summary(self, batch_path: str, pdb_id: str, chain_id: str,
                           reference: str, blast_only: bool) -> Optional[str]:
        """Find domain summary file"""
        from ecod.utils.path_utils import get_all_evidence_paths

        # Get all evidence paths
        all_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

        # Choose correct summary type
        summary_type = 'blast_only_summary' if blast_only else 'domain_summary'

        if summary_type in all_paths and all_paths[summary_type]['exists_at']:
            path = all_paths[summary_type]['exists_at']
            logging.getLogger(__name__).info(f"Found domain summary at: {path}")
            return path

        return None
