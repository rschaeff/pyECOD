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
from typing import Dict, Any, List, Optional, Tuple, Set, Union

from ecod.exceptions import PipelineError, FileOperationError
from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.models.pipeline import BlastHit, HHSearchHit, DomainSummaryModel, PipelineResult
from ecod.models.domain_analysis import DomainPartitionResult
from ecod.models.domain_analysis.domain_candidate import DomainCandidate
from ecod.models.domain import Domain, DomainRange, DomainRangeSegment
from ecod.utils.xml_utils import ensure_dict, ensure_list_of_dicts
from ecod.utils.path_utils import get_standardized_paths, get_all_evidence_paths, resolve_file_path
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

    #########################################
    # Main processing methods
    #########################################

    def process_batch(self, batch_id: int, batch_path: str, reference: str,
                     blast_only: bool = False, limit: int = None,
                     reps_only: bool = False) -> List[DomainPartitionResult]:
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

        # Get proteins to process
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

    def _process_domains_internal(self, pdb_id: str, chain_id: str,
                                 domain_summary_path: str, domain_file: str,
                                 reference: str) -> Union[List[Dict[str, Any]], str]:
        """
        Process domains for a protein chain - internal implementation

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            domain_summary_path: Path to domain summary file
            domain_file: Path to output domain file
            reference: Reference version

        Returns:
            List of domain dictionaries if successful, error message string if failed
        """
        logger = logging.getLogger(__name__)
        logger.info(f"Processing domains for {pdb_id}_{chain_id}")

        try:
            # 1. Process domain summary
            blast_data = self._process_domain_summary(domain_summary_path)

            if not blast_data or "error" in blast_data:
                error_msg = f"Failed to parse domain summary for {pdb_id}_{chain_id}"
                if "error" in blast_data:
                    error_msg += f": {blast_data['error']}"
                logger.error(error_msg)
                return f"ERROR: {error_msg}"

            logger.info(f"Summary parsed successfully for {pdb_id}_{chain_id}")

            # Log summary stats
            sequence_length = blast_data.get("sequence_length", 0)
            is_peptide = blast_data.get("is_peptide", False)

            logger.info(f"Chain length: {sequence_length}")
            logger.info(f"Is peptide: {is_peptide}")
            logger.info(f"Chain BLAST hits: {len(blast_data.get('chain_blast_hits', []))}")
            logger.info(f"Domain BLAST hits: {len(blast_data.get('domain_blast_hits', []))}")
            logger.info(f"HHSearch hits: {len(blast_data.get('hhsearch_hits', []))}")

            # 2. Check if this is a peptide-length chain
            if is_peptide:
                logger.info(f"Chain {pdb_id}_{chain_id} marked as peptide")
                return []

            if sequence_length < 20:
                logger.info(f"Chain {pdb_id}_{chain_id} classified as peptide-length (length={sequence_length})")
                return []

            # 3. Determine domain boundaries
            logger.info(f"Detecting domain boundaries for {pdb_id}_{chain_id}")
            domains = self._determine_domain_boundaries(blast_data, sequence_length, f"{pdb_id}_{chain_id}")

            if not domains:
                logger.warning(f"No domain candidates found for {pdb_id}_{chain_id}")
                return []

            logger.info(f"Found {len(domains)} domain candidates")

            # 4. Assign domain classifications
            logger.info(f"Classifying domains for {pdb_id}_{chain_id}")
            self._assign_domain_classifications(domains, blast_data, f"{pdb_id}_{chain_id}")

            # 5. Create domain document
            domain_doc, domain_stats = self._create_domain_document(
                pdb_id, chain_id, reference, domains, sequence_length
            )

            # 6. Write output file
            os.makedirs(os.path.dirname(domain_file), exist_ok=True)
            tree = ET.ElementTree(domain_doc)
            tree.write(domain_file, encoding='utf-8', xml_declaration=True)

            logger.info(f"Created domain partition file with {len(domains)} domains")

            return domains

        except Exception as e:
            error_msg = f"Error processing domains for {pdb_id}_{chain_id}: {str(e)}"
            logger.error(error_msg)
            import traceback
            logger.error(traceback.format_exc())
            return f"ERROR: {str(e)}"

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

    #########################################
    # Reference data methods
    #########################################

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

    def _get_domain_classification(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """
        Get classification details for a domain

        Args:
            domain_id: Domain identifier

        Returns:
            Dictionary with classification details or None
        """
        try:
            query = """
            SELECT
                domain_id,
                t_group,
                h_group,
                x_group,
                a_group
            FROM
                pdb_analysis.domain
            WHERE
                domain_id = %s
            LIMIT 1
            """
            results = self.context.db.execute_dict_query(query, (domain_id,))

            if results:
                return {
                    'domain_id': results[0]['domain_id'],
                    't_group': results[0]['t_group'],
                    'h_group': results[0]['h_group'],
                    'x_group': results[0]['x_group'],
                    'a_group': results[0]['a_group']
                }
            return None
        except Exception as e:
            logging.getLogger(__name__).error(f"Error getting domain classification: {str(e)}")
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

    def _get_reference_domains(self, pdb_id: str, chain_id: str) -> List[Dict[str, Any]]:
        """
        Get reference domains for a PDB chain

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier

        Returns:
            List of domain dictionaries from reference database
        """
        # This would typically query your reference database
        # Simplified implementation for explanation purposes
        try:
            query = """
            SELECT
                d.domain_id,
                d.range,
                d.t_group,
                d.h_group,
                d.x_group,
                d.a_group
            FROM
                pdb_analysis.domain d
            JOIN
                pdb_analysis.protein p ON d.protein_id = p.id
            WHERE
                p.pdb_id = %s AND p.chain_id = %s
            """
            results = self.context.db.execute_dict_query(query, (pdb_id, chain_id))

            domains = []
            for row in results:
                domain = {
                    'domain_id': row['domain_id'],
                    'range': row['range'],
                    't_group': row['t_group'],
                    'h_group': row['h_group'],
                    'x_group': row['x_group'],
                    'a_group': row['a_group']
                }
                # Parse range
                ranges = self._parse_range(row['range'])
                if ranges:
                    domain['start'] = ranges[0][0]
                    domain['end'] = ranges[-1][1]
                    domains.append(domain)

            return domains
        except Exception as e:
            logging.getLogger(__name__).error(f"Error getting reference domains: {str(e)}")
            return []

    def _get_reference_domain_by_id(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """
        Get reference domain by ID

        Args:
            domain_id: Domain identifier

        Returns:
            Domain dictionary if found, None otherwise
        """
        try:
            query = """
            SELECT
                d.domain_id,
                d.t_group,
                d.h_group,
                d.x_group,
                d.a_group,
                d.is_manual_rep,
                d.is_f70,
                d.is_f40
            FROM
                pdb_analysis.domain d
            WHERE
                d.domain_id = %s
            LIMIT 1
            """
            results = self.context.db.execute_dict_query(query, (domain_id,))

            if results:
                return {
                    "domain_id": results[0]["domain_id"],
                    "t_group": results[0]["t_group"],
                    "h_group": results[0]["h_group"],
                    "x_group": results[0]["x_group"],
                    "a_group": results[0]["a_group"],
                    "is_manual_rep": results[0]["is_manual_rep"],
                    "is_f70": results[0]["is_f70"],
                    "is_f40": results[0]["is_f40"]
                }

            return None
        except Exception as e:
            logging.getLogger(__name__).error(f"Error getting reference domain: {str(e)}")
            return None

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

    #########################################
    # Domain determination and classification methods
    #########################################

    def detect_boundaries(self, summary: DomainSummaryModel) -> List[DomainCandidate]:
        """Detect domain boundaries from evidence

        Args:
            summary: Domain summary model with hit data

        Returns:
            List of domain candidates
        """
        logger = self.logger
        logger.info(f"Detecting domain boundaries for {summary.pdb_id}_{summary.chain_id}")

        # Step 1: Check for high-confidence domains from HHSearch
        high_confidence_candidates = []
        if summary.hhsearch_hits:
            for hit in summary.hhsearch_hits:
                if hit.probability >= 90.0:  # High confidence threshold
                    # Convert to evidence
                    evidence = Evidence(
                        type="hhsearch",
                        source_id=hit.domain_id or hit.hit_id,
                        query_range=hit.range,
                        hit_range=hit.hit_range,
                        confidence=hit.probability / 100.0  # Normalize to 0-1
                    )

                    # Extract positions using range utilities
                    from ecod.utils.range_utils import parse_range
                    for start, end in parse_range(hit.range):
                        candidate = DomainCandidate(
                            start=start,
                            end=end,
                            evidence=[evidence],
                            confidence=hit.probability / 100.0,
                            source="hhsearch",
                            protected=hit.probability >= 98.0  # Protected if very high confidence
                        )
                        high_confidence_candidates.append(candidate)

        # Step 2: Check for domains from chain BLAST
        chain_blast_candidates = []
        if summary.chain_blast_hits:
            # Group hits by source chain to find conserved domain architectures
            chain_groups = {}
            for hit in summary.chain_blast_hits:
                key = f"{hit.pdb_id}_{hit.chain_id}"
                if key not in chain_groups:
                    chain_groups[key] = []
                chain_groups[key].append(hit)

            # Process each chain group
            for source_chain, hits in chain_groups.items():
                # Create list to track source regions
                mapped_regions = []

                # Find reference domains for this chain
                ref_domains = self._get_reference_domains(hits[0].pdb_id, hits[0].chain_id)

                if not ref_domains:
                    continue

                # Map reference domains to query
                for hit in hits:
                    # Parse ranges using utility
                    from ecod.utils.range_utils import parse_range
                    query_ranges = parse_range(hit.range)
                    hit_ranges = parse_range(hit.hit_range)

                    # Map each reference domain to query
                    for ref_domain in ref_domains:
                        ref_start = ref_domain.get("start", 0)
                        ref_end = ref_domain.get("end", 0)

                        for (q_start, q_end), (h_start, h_end) in zip(query_ranges, hit_ranges):
                            # Check if hit region overlaps reference domain
                            if max(h_start, ref_start) <= min(h_end, ref_end):
                                # Calculate mapping ratio
                                h_length = h_end - h_start + 1
                                q_length = q_end - q_start + 1

                                # Map domain to query coordinates
                                mapped_start = q_start + round((ref_start - h_start) * q_length / h_length)
                                mapped_end = q_start + round((ref_end - h_start) * q_length / h_length)

                                # Create evidence
                                evidence = Evidence(
                                    type="chain_blast",
                                    source_id=ref_domain.get("domain_id", ""),
                                    query_range=f"{mapped_start}-{mapped_end}",
                                    hit_range=f"{ref_start}-{ref_end}",
                                    confidence=1.0 / (1.0 + hit.evalue)  # Convert evalue to confidence
                                )

                                # Add classification if available
                                if "t_group" in ref_domain and ref_domain["t_group"]:
                                    evidence.t_group = ref_domain["t_group"]
                                    evidence.h_group = ref_domain.get("h_group")
                                    evidence.x_group = ref_domain.get("x_group")
                                    evidence.a_group = ref_domain.get("a_group")

                                # Create domain candidate
                                candidate = DomainCandidate(
                                    start=mapped_start,
                                    end=mapped_end,
                                    evidence=[evidence],
                                    confidence=1.0 / (1.0 + hit.evalue),
                                    source="chain_blast"
                                )

                                # Store classification
                                candidate.t_group = evidence.t_group
                                candidate.h_group = evidence.h_group
                                candidate.x_group = evidence.x_group
                                candidate.a_group = evidence.a_group

                                chain_blast_candidates.append(candidate)

        # Step 3: Check for domains from domain BLAST
        domain_blast_candidates = []
        if summary.domain_blast_hits:
            # Create positional bins to group overlapping hits
            position_bins = {}

            # Bin hits by position
            for hit in summary.domain_blast_hits:
                from ecod.utils.range_utils import parse_range
                for start, end in parse_range(hit.range):
                    bin_key = start // 50  # Bin by ~50 residue windows

                    if bin_key not in position_bins:
                        position_bins[bin_key] = []

                    position_bins[bin_key].append({
                        "hit": hit,
                        "start": start,
                        "end": end,
                        "evalue": hit.evalue,
                        "domain_id": hit.domain_id
                    })

            # For each position bin, find consensus domain
            for bin_key, hits in position_bins.items():
                if len(hits) < 3:  # Require at least 3 hits for confidence
                    continue

                # Find median boundaries
                starts = sorted(h["start"] for h in hits)
                ends = sorted(h["end"] for h in hits)

                median_start = starts[len(starts) // 2]
                median_end = ends[len(ends) // 2]

                # Get best hit (lowest evalue)
                best_hit = min(hits, key=lambda h: h["evalue"])

                # Create evidence from best hit
                evidence = Evidence(
                    type="domain_blast",
                    source_id=best_hit["domain_id"],
                    query_range=f"{median_start}-{median_end}",
                    hit_range=best_hit["hit"]["hit_range"],
                    confidence=1.0 / (1.0 + best_hit["evalue"])
                )

                # Create domain candidate
                candidate = DomainCandidate(
                    start=median_start,
                    end=median_end,
                    evidence=[evidence],
                    confidence=1.0 / (1.0 + best_hit["evalue"]),
                    source="domain_blast"
                )

                domain_blast_candidates.append(candidate)

        # Step 4: Merge all candidates
        all_candidates = high_confidence_candidates + chain_blast_candidates + domain_blast_candidates

        # Step 5: Resolve overlaps
        final_candidates = []

        # Sort by confidence (descending)
        sorted_candidates = sorted(all_candidates, key=lambda c: c.confidence, reverse=True)

        # Track covered positions
        covered_positions = set()

        for candidate in sorted_candidates:
            from ecod.utils.range_utils import range_to_positions
            positions = range_to_positions(candidate.range)

            # Calculate overlap with existing domains
            overlap = len(positions.intersection(covered_positions))
            overlap_pct = overlap / len(positions) if positions else 0

            # If minimal overlap or this is a protected high-confidence domain, include it
            if overlap_pct < 0.3 or candidate.protected:
                final_candidates.append(candidate)
                covered_positions.update(positions)

        # Sort by position for clarity
        final_candidates.sort(key=lambda c: c.start)

        logger.info(f"Detected {len(final_candidates)} domain boundaries")
        return final_candidates

    def _determine_domain_boundaries(self, pdb_id: str, chain_id: str, summary: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Determine domain boundaries based on summary information

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            summary: Domain summary data (from _process_domain_summary)

        Returns:
            List of domain dictionaries with boundary information
        """
        logger = logging.getLogger(__name__)
        logger.info(f"Determining domain boundaries for {pdb_id}_{chain_id} (length: {summary.get('sequence_length', 0)})")

        # Initialize list for domain candidates
        domain_candidates = []

        # Get candidates from HHSearch hits (if available)
        hhsearch_candidates = self._identify_domains_from_hhsearch(
            pdb_id, chain_id, summary.get('hhsearch_hits', [])
        )
        if hhsearch_candidates:
            logger.info(f"Using {len(hhsearch_candidates)} high-confidence domain hits as domain boundaries")
            domain_candidates.extend(hhsearch_candidates)
        else:
            logger.warning("No HHSearch hits or invalid sequence length")

        # Get candidates from chain BLAST hits (if needed)
        if not domain_candidates or not self._is_fully_covered(domain_candidates, summary.get('sequence_length', 0)):
            chain_blast_candidates = self._identify_domains_from_chain_blast(
                pdb_id, chain_id, summary.get('chain_blast_hits', [])
            )
            if chain_blast_candidates:
                logger.info(f"Found {len(chain_blast_candidates)} domain candidates from chain BLAST hits")
                domain_candidates.extend(chain_blast_candidates)

        # Get candidates from domain BLAST hits (if needed)
        if not domain_candidates or not self._is_fully_covered(domain_candidates, summary.get('sequence_length', 0)):
            domain_blast_candidates = self._identify_domains_from_domain_blast(
                pdb_id, chain_id, summary.get('domain_blast_hits', [])
            )
            if domain_blast_candidates:
                logger.info(f"Found {len(domain_blast_candidates)} domain candidates from domain BLAST hits")
                domain_candidates.extend(domain_blast_candidates)

        # Log domain candidate sources
        logger.info(
            f"Domain candidates: "
            f" {len([d for d in domain_candidates if d.get('source') == 'hhsearch'])} from HHSearch, "
            f" {len([d for d in domain_candidates if d.get('source') == 'chain_blast'])} from chain BLAST, "
            f" {len([d for d in domain_candidates if d.get('source') == 'domain_blast'])} from domain BLAST"
        )

        # Resolve overlapping domain boundaries
        final_domains = self._resolve_domain_boundaries(
            domain_candidates,
            summary.get('sequence_length', 0)
        )

        logger.info(f"Final domain boundaries for {pdb_id}_{chain_id} after filtering: {len(final_domains)} domains")

        return final_domains

    def _identify_domains_from_hhsearch(self, pdb_id: str, chain_id: str, hhsearch_hits: List[Any]) -> List[Dict[str, Any]]:
        """
        Identify domains from HHSearch hits

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            hhsearch_hits: List of HHSearch hits from summary

        Returns:
            List of domain dictionaries from HHSearch
        """
        logger = logging.getLogger(__name__)

        # Filter hits by probability threshold
        high_prob_hits = [hit for hit in hhsearch_hits
                         if hasattr(hit, 'probability') and hit.probability >= 90.0]

        if not high_prob_hits:
            return []

        # Create domain candidates from high probability hits
        domains = []
        for hit in high_prob_hits:
            # Parse range
            ranges = []
            if hasattr(hit, 'range') and hit.range:
                ranges = self._parse_range(hit.range)
            else:
                continue

            if not ranges:
                continue

            # Create domain for each range segment
            for start, end in ranges:
                # Convert hit to dictionary for evidence
                hit_dict = {
                    'type': 'hhsearch',
                    'domain_id': hit.domain_id if hasattr(hit, 'domain_id') else '',
                    'query_range': hit.range if hasattr(hit, 'range') else '',
                    'hit_range': hit.hit_range if hasattr(hit, 'hit_range') else '',
                    'probability': hit.probability if hasattr(hit, 'probability') else 0.0,
                    'evalue': hit.evalue if hasattr(hit, 'evalue') else 999.0,
                    'score': hit.score if hasattr(hit, 'score') else 0.0
                }

                domain = {
                    'start': start,
                    'end': end,
                    'range': f"{start}-{end}",
                    'source': 'hhsearch',
                    'confidence': hit.probability / 100.0 if hasattr(hit, 'probability') else 0.0,
                    'source_id': hit.domain_id if hasattr(hit, 'domain_id') else
                                 (hit.hit_id if hasattr(hit, 'hit_id') else ''),
                    't_group': None,  # Will be set during classification
                    'h_group': None,
                    'x_group': None,
                    'a_group': None,
                    'evidence': [hit_dict]  # Store as dictionary for compatibility
                }
                domains.append(domain)

        return domains

    def _identify_domains_from_chain_blast(self, pdb_id: str, chain_id: str, chain_blast_hits: List[Any]) -> List[Dict[str, Any]]:
        """
        Identify domains from chain BLAST hits

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            chain_blast_hits: List of chain BLAST hits from summary

        Returns:
            List of domain dictionaries from chain BLAST
        """
        logger = logging.getLogger(__name__)

        if not chain_blast_hits:
            logger.warning("No chain BLAST hits provided to analyze")
            return []

        # Process each hit to find domains
        domains = []
        for i, hit in enumerate(chain_blast_hits[:30]):  # Limit to first 30 hits for performance
            # Skip hits without required fields
            if not (hasattr(hit, 'pdb_id') and hasattr(hit, 'chain_id')):
                continue

            # Parse range
            ranges = []
            if hasattr(hit, 'range') and hit.range:
                ranges = self._parse_range(hit.range)
            elif hasattr(hit, 'range_parsed') and hit.range_parsed:
                ranges = hit.range_parsed
            else:
                continue

            if not ranges:
                continue

            # Get reference domains for this hit
            hit_domains = self._get_reference_domains(
                hit.pdb_id if hasattr(hit, 'pdb_id') else '',
                hit.chain_id if hasattr(hit, 'chain_id') else ''
            )
            if not hit_domains:
                continue

            # Map reference domains to query sequence
            mapped_domains = self._map_domains_to_query(ranges, hit_domains)

            if mapped_domains:
                # Convert hit to dictionary for evidence
                hit_dict = {
                    'type': 'chain_blast',
                    'hit_id': hit.hit_id if hasattr(hit, 'hit_id') else '',
                    'domain_id': hit.domain_id if hasattr(hit, 'domain_id') else '',
                    'pdb_id': hit.pdb_id if hasattr(hit, 'pdb_id') else '',
                    'chain_id': hit.chain_id if hasattr(hit, 'chain_id') else '',
                    'evalue': hit.evalue if hasattr(hit, 'evalue') else 999.0,
                    'query_range': hit.range if hasattr(hit, 'range') else '',
                    'hit_range': hit.hit_range if hasattr(hit, 'hit_range') else ''
                }

                # Update evidence in mapped domains
                for domain in mapped_domains:
                    domain['evidence'] = [hit_dict]

                logger.info(f"Mapped {len(mapped_domains)} domains from hit #{i+1}")
                domains.extend(mapped_domains)

        # Log statistics
        logger.info(f"Chain BLAST hit analysis summary:")
        logger.info(f"  Total hits: {len(chain_blast_hits)}")
        logger.info(f"  Hits processed: {min(30, len(chain_blast_hits))}")
        logger.info(f"  Total domains mapped: {len(domains)}")

        return domains

    def _identify_domains_from_domain_blast(self, pdb_id: str, chain_id: str, domain_blast_hits: List[Any]) -> List[Dict[str, Any]]:
        """
        Identify domains from domain BLAST hits

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            domain_blast_hits: List of domain BLAST hits from summary

        Returns:
            List of domain dictionaries from domain BLAST
        """
        if not domain_blast_hits:
            return []

        # Group hits by query region
        region_hits = {}
        for hit in domain_blast_hits:
            # Parse range
            ranges = []
            if hasattr(hit, 'range') and hit.range:
                ranges = self._parse_range(hit.range)
            elif hasattr(hit, 'range_parsed') and hit.range_parsed:
                ranges = hit.range_parsed
            else:
                continue

            for start, end in ranges:
                region_key = f"{start}-{end}"
                if region_key not in region_hits:
                    region_hits[region_key] = []
                region_hits[region_key].append(hit)

        # Create domain candidates from regions with sufficient hits
        domains = []
        for region, hits in region_hits.items():
            if len(hits) < 3:  # Require at least 3 hits for confidence
                continue

            start, end = map(int, region.split('-'))

            # Find most common classification among hits
            t_groups = {}
            for hit in hits:
                if hasattr(hit, 'domain_id') and hit.domain_id:
                    domain_info = self._get_domain_classification_by_id(hit.domain_id)
                    if domain_info and 't_group' in domain_info:
                        t_group = domain_info['t_group']
                        t_groups[t_group] = t_groups.get(t_group, 0) + 1

            # Use most common t_group if available
            t_group = None
            if t_groups:
                t_group = max(t_groups.items(), key=lambda x: x[1])[0]

            # Convert top 5 hits to dictionaries for evidence
            evidence = []
            for hit in hits[:5]:  # Include top 5 hits as evidence
                hit_dict = {
                    'type': 'domain_blast',
                    'hit_id': hit.hit_id if hasattr(hit, 'hit_id') else '',
                    'domain_id': hit.domain_id if hasattr(hit, 'domain_id') else '',
                    'pdb_id': hit.pdb_id if hasattr(hit, 'pdb_id') else '',
                    'chain_id': hit.chain_id if hasattr(hit, 'chain_id') else '',
                    'evalue': hit.evalue if hasattr(hit, 'evalue') else 999.0,
                    'query_range': hit.range if hasattr(hit, 'range') else '',
                    'hit_range': hit.hit_range if hasattr(hit, 'hit_range') else ''
                }
                evidence.append(hit_dict)

            domain = {
                'start': start,
                'end': end,
                'range': region,
                'source': 'domain_blast',
                'confidence': 0.7,  # Medium confidence level
                'source_id': '',  # Multiple hits, no single source
                't_group': t_group,
                'h_group': None,  # Will be set during classification
                'x_group': None,
                'a_group': None,
                'evidence': evidence
            }
            domains.append(domain)

        return domains

    def _map_domains_to_query(self, query_ranges: List[Tuple[int, int]], reference_domains: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Map reference domains to query sequence

        Args:
            query_ranges: List of query ranges (start, end)
            reference_domains: List of reference domains

        Returns:
            List of mapped domain dictionaries
        """
        if not query_ranges or not reference_domains:
            return []

        logger = logging.getLogger(__name__)

        # Simplest case: single range mapping
        if len(query_ranges) == 1:
            q_start, q_end = query_ranges[0]
            q_length = q_end - q_start + 1

            mapped_domains = []
            for ref_domain in reference_domains:
                ref_start = ref_domain.get('start', 0)
                ref_end = ref_domain.get('end', 0)

                if ref_start <= 0 or ref_end <= 0:
                    continue

                ref_length = ref_end - ref_start + 1

                # Calculate mapping ratio
                ratio = q_length / ref_length if ref_length > 0 else 1.0

                # Map domain to query coordinates
                domain = ref_domain.copy()
                domain.update({
                    'start': max(1, round(q_start + (ref_start - 1) * ratio)),
                    'end': min(q_end, round(q_start + (ref_end - 1) * ratio)),
                    'source': 'chain_blast',
                    'confidence': 0.8,  # High confidence for chain blast
                    'evidence': [{'type': 'chain_blast', 'source_id': f"{ref_domain.get('domain_id', '')}"}]
                })

                # Validate mapped range
                if domain['start'] <= domain['end']:
                    domain['range'] = f"{domain['start']}-{domain['end']}"
                    mapped_domains.append(domain)
                    logger.info(f"Successfully mapped domain {ref_domain.get('domain_id', '')}: {domain['range']} ({domain['end']-domain['start']+1} residues)")

        # More complex case: multi-segment ranges would need more sophisticated mapping
        else:
            logger.warning("Multi-segment query ranges not fully supported for domain mapping")

        return mapped_domains

    def _resolve_domain_boundaries(self, candidates: List[Dict[str, Any]], sequence_length: int) -> List[Dict[str, Any]]:
        """
        Resolve overlapping domain boundaries

        Args:
            candidates: List of domain candidates
            sequence_length: Length of the protein sequence

        Returns:
            List of resolved domains with non-overlapping boundaries
        """
        if not candidates:
            return []

        logger = logging.getLogger(__name__)

        # Sort candidates by confidence and then by size (larger domains first)
        sorted_candidates = sorted(
            candidates,
            key=lambda d: (d.get('confidence', 0), d.get('end', 0) - d.get('start', 0) + 1),
            reverse=True
        )

        # Initialize coverage tracker
        covered = set()
        final_domains = []

        # Process candidates in order
        for candidate in sorted_candidates:
            start = candidate.get('start', 0)
            end = candidate.get('end', 0)

            if start <= 0 or end <= 0 or end < start:
                continue

            # Check overlap with existing domains
            domain_positions = set(range(start, end + 1))
            overlap = domain_positions.intersection(covered)

            # If no significant overlap, add domain
            if len(overlap) < 0.2 * len(domain_positions):
                final_domains.append(candidate)
                covered.update(domain_positions)

        # Sort final domains by position
        final_domains.sort(key=lambda d: d.get('start', 0))

        # Log coverage
        if sequence_length > 0:
            coverage_pct = (len(covered) / sequence_length) * 100
            logger.info(f"Domain coverage: {len(covered)}/{sequence_length} residues ({coverage_pct:.1f}%)")

        return final_domains

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

    def _is_fully_covered(self, domains: List[Dict[str, Any]], sequence_length: int, threshold: float = 0.9) -> bool:
        """
        Check if the domains cover the full sequence

        Args:
            domains: List of domain dictionaries
            sequence_length: Length of the protein sequence
            threshold: Coverage threshold (0.0-1.0)

        Returns:
            True if domains cover at least threshold% of the sequence
        """
        if not domains or sequence_length == 0:
            return False

        # Create a set of all positions covered
        covered = set()
        for domain in domains:
            start = domain.get('start', 0)
            end = domain.get('end', 0)

            if start > 0 and end > 0 and end >= start:
                covered.update(range(start, end + 1))

        # Calculate coverage percentage
        coverage = len(covered) / sequence_length

        return coverage >= threshold

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

    def _classify_domains(self, domains: List[Dict[str, Any]], reference: str) -> None:
        """
        Assign domain classifications

        Args:
            domains: List of domain dictionaries
            reference: Reference version
        """
        logger = logging.getLogger(__name__)

        # This would typically be a complex algorithm
        # Simplified implementation for explanation purposes
        unclassified_count = 0
        for domain in domains:
            if not domain.get('t_group'):
                unclassified_count += 1

        if unclassified_count:
            logger.warning(f"Could not classify {unclassified_count} domains")

    #########################################
    # File processing and I/O methods
    #########################################

    def _find_domain_summary(self, batch_path: str, pdb_id: str, chain_id: str,
                           reference: str, blast_only: bool = False) -> str:
        """Find domain summary file using path utilities to check standard and legacy paths

        Args:
            batch_path: Base directory for I/O operations
            pdb_id: PDB identifier
            chain_id: Chain identifier
            reference: Reference version
            blast_only: Whether to look for BLAST-only summary

        Returns:
            Path to domain summary file if found, empty string otherwise
        """
        from ecod.utils.path_utils import get_all_evidence_paths, resolve_file_path

        # Use path_utils to get standard and legacy paths
        file_type = 'blast_only_summary' if blast_only else 'domain_summary'

        try:
            # Use path_utils to check all possible paths
            evidence_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

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
                full_summ_path = resolve_file_path(batch_path, db_summ_path)
                if os.path.exists(full_summ_path):
                    self.logger.info(f"Found domain summary in database: {full_summ_path}")
                    return full_summ_path
        except Exception as e:
            self.logger.warning(f"Error querying database for domain summary: {e}")

        # If all else fails, return empty string
        self.logger.warning(f"Domain summary not found for {pdb_id}_{chain_id}")
        return ""


    def _find_fasta_file(self, pdb_id: str, chain_id: str, base_dir: str) -> Optional[str]:
        """
        Find FASTA file for a protein chain

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            base_dir: Base directory to search in

        Returns:
            Path to FASTA file if found, None otherwise
        """
        # Check canonical locations
        fasta_dir = os.path.join(os.path.dirname(base_dir), "fastas")

        # Try standard FASTA file pattern
        fasta_path = os.path.join(fasta_dir, f"{pdb_id}_{chain_id}.fa")
        if os.path.exists(fasta_path):
            return fasta_path

        # Try alternative patterns
        patterns = [
            f"{pdb_id}_{chain_id}.fasta",
            f"{pdb_id.lower()}_{chain_id}.fa",
            f"{pdb_id.lower()}_{chain_id}.fasta",
            f"{pdb_id}_{chain_id.lower()}.fa",
            f"{pdb_id}_{chain_id.lower()}.fasta"
        ]

        for pattern in patterns:
            path = os.path.join(fasta_dir, pattern)
            if os.path.exists(path):
                return path

        return None

    def _read_fasta_sequence(self, fasta_path: str) -> str:
        """
        Read sequence from FASTA file

        Args:
            fasta_path: Path to FASTA file

        Returns:
            Sequence string
        """
        sequence = ""

        try:
            with open(fasta_path, 'r') as f:
                lines = f.readlines()

            for line in lines:
                line = line.strip()
                if not line or line.startswith('>'):
                    continue
                sequence += line

            return sequence
        except Exception as e:
            logging.getLogger(__name__).error(f"Error reading FASTA file: {str(e)}")
            return ""

    def _get_sequence_length(self, pdb_id: str, chain_id: str, domain_summary_path: str) -> int:
        """
        Get sequence length for a protein chain

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            domain_summary_path: Path to domain summary file (for resolving fasta location)

        Returns:
            Sequence length if available, 0 otherwise
        """
        logger = logging.getLogger(__name__)
        self.logger.info(f"Getting sequence length for {pdb_id}_{chain_id}")

        # First try to get from database
        self.logger.info(f"Trying to get sequence length from database for {pdb_id}_{chain_id}")
        try:
            query = """
            SELECT length FROM ecod_schema.protein
            WHERE pdb_id = %s AND chain_id = %s
            LIMIT 1
            """
            results = self.context.db.execute_query(query, (pdb_id, chain_id))
            if results:
                if results[0][0]:
                    length = results[0][0]
                    self.logger.info(f"Got sequence length from database for {pdb_id}_{chain_id}: {length}")
                    return length
                else:
                    self.logger.warning(f"Database returned NULL length for {pdb_id}_{chain_id}")
            else:
                self.logger.warning(f"No records found in database for {pdb_id}_{chain_id}")
        except Exception as e:
            logger.warning(f"Error getting sequence length from database for {pdb_id}_{chain_id}: {str(e)}")

        # Try to get from FASTA file
        self.logger.info(f"Trying to get sequence length from FASTA file for {pdb_id}_{chain_id}")
        fasta_path = self._find_fasta_file(pdb_id, chain_id, os.path.dirname(domain_summary_path))
        if fasta_path:
            self.logger.info(f"Found FASTA file at {fasta_path}")
            sequence = self._read_fasta_sequence(fasta_path)
            if sequence:
                length = len(sequence)
                self.logger.info(f"Got sequence from FASTA for {pdb_id}_{chain_id}, length: {length}")
                return length
            else:
                self.logger.warning(f"Could not read sequence from FASTA file for {pdb_id}_{chain_id}")
        else:
            self.logger.warning(f"No FASTA file found for {pdb_id}_{chain_id}")

        # Check if there's a sequence in the domain summary XML itself
        self.logger.info(f"Checking domain summary XML for sequence length for {pdb_id}_{chain_id}")
        try:
            tree = ET.parse(domain_summary_path)
            root = tree.getroot()

            # Try to find sequence length in XML
            query_len_elem = root.find(".//query_len")
            if query_len_elem is not None and query_len_elem.text:
                try:
                    length = int(query_len_elem.text.strip())
                    self.logger.info(f"Found sequence length in domain summary XML: {length}")
                    return length
                except ValueError:
                    self.logger.warning(f"Invalid sequence length in XML: {query_len_elem.text}")
            else:
                self.logger.warning("No sequence length found in domain summary XML")

        except Exception as e:
            self.logger.warning(f"Error parsing domain summary XML: {str(e)}")

        # Default to 0
        self.logger.error(f"Failed to get sequence length for {pdb_id}_{chain_id}, defaulting to 0")
        return 0

    def _extract_domain_sequence(self, full_sequence: str, start: int, end: int) -> str:
        """
        Extract domain sequence from full protein sequence

        Args:
            full_sequence: Full protein sequence
            start: Domain start position (1-indexed)
            end: Domain end position (1-indexed)

        Returns:
            Domain sequence
        """
        if not full_sequence:
            return ""

        # Convert to 0-indexed
        start_idx = start - 1
        end_idx = end

        # Validate indices
        if start_idx < 0:
            start_idx = 0
        if end_idx > len(full_sequence):
            end_idx = len(full_sequence)

        # Extract sequence
        return full_sequence[start_idx:end_idx]

    def _process_domain_summary(self, domain_summary_path: str) -> Dict[str, Any]:
        """
        Process domain summary file

        Args:
            domain_summary_path: Path to domain summary file

        Returns:
            Dictionary with summary information
        """
        logger = logging.getLogger("ecod.domain_partition")
        logger.info(f"ENTRY: _process_domain_summary for {domain_summary_path}")

        # Check file existence
        if not os.path.exists(domain_summary_path):
            logger.error(f"Domain summary file not found: {domain_summary_path}")
            return {"error": "File not found"}

        try:
            logger.info("Starting to parse XML")
            # Parse XML
            tree = ET.parse(domain_summary_path)
            root = tree.getroot()

            # Log root tag and structure for debugging
            logger.info(f"XML parsed successfully, root tag: {root.tag}, attributes: {root.attrib}")

            # Get basic information
            summary_elem = root.find("blast_summ")
            if summary_elem is None:
                logger.error(f"Invalid domain summary format: missing blast_summ element")
                return {}

            pdb_id = summary_elem.get("pdb", "")
            chain_id = summary_elem.get("chain", "")
            logger.info(f"Found basic info: pdb_id={pdb_id}, chain_id={chain_id}")

            # Create summary dictionary
            logger.info("Creating summary dictionary")
            summary = {
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "is_peptide": summary_elem.get("is_peptide", "false").lower() == "true",
                "chain_blast_hits": [],
                "domain_blast_hits": [],
                "hhsearch_hits": []
            }

            # Log if is_peptide is true
            if summary["is_peptide"]:
                logger.info(f"Chain {pdb_id}_{chain_id} is explicitly marked as peptide in XML")

            logger.info("Processing chain BLAST hits")
            # Get chain BLAST hits
            chain_blast_run = root.find("chain_blast_run")
            if chain_blast_run is not None:
                logger.info("Found chain_blast_run element")
                hits_elem = chain_blast_run.find("hits")
                if hits_elem is not None:
                    hit_count = 0
                    for hit_elem in hits_elem.findall("hit"):
                        hit_count += 1
                        # Create BlastHit from element
                        from ecod.models.pipeline import BlastHit
                        try:
                            blast_hit = BlastHit.from_xml(hit_elem)
                            blast_hit.hit_type = "chain_blast"  # Set hit type

                            # Add hit to summary
                            summary["chain_blast_hits"].append(blast_hit)
                        except Exception as e:
                            logger.warning(f"Error parsing chain BLAST hit: {e}")

                    logger.info(f"Processed {hit_count} chain BLAST hits")
                else:
                    logger.info("No hits element found in chain_blast_run")
            else:
                logger.info("No chain_blast_run element found")

            logger.info("Processing domain BLAST hits")
            # Get domain BLAST hits
            blast_run = root.find("blast_run")
            if blast_run is not None:
                logger.info("Found blast_run element")
                hits_elem = blast_run.find("hits")
                if hits_elem is not None:
                    hit_count = 0
                    for hit_elem in hits_elem.findall("hit"):
                        hit_count += 1
                        # Create BlastHit from element
                        from ecod.models.pipeline import BlastHit
                        try:
                            blast_hit = BlastHit.from_xml(hit_elem)
                            blast_hit.hit_type = "domain_blast"  # Set hit type

                            # Add hit to summary
                            summary["domain_blast_hits"].append(blast_hit)
                        except Exception as e:
                            logger.warning(f"Error parsing domain BLAST hit: {e}")

                    logger.info(f"Processed {hit_count} domain BLAST hits")
                else:
                    logger.info("No hits element found in blast_run")
            else:
                logger.info("No blast_run element found")

            logger.info("Processing HHSearch hits")
            # Get HHSearch hits
            hh_run = root.find("hh_run")
            if hh_run is not None:
                logger.info("Found hh_run element")
                hits_elem = hh_run.find("hits")
                if hits_elem is not None:
                    hit_count = 0
                    for hit_elem in hits_elem.findall("hit"):
                        hit_count += 1
                        # Create HHSearchHit from element
                        from ecod.models.pipeline import HHSearchHit
                        try:
                            hhsearch_hit = HHSearchHit.from_xml(hit_elem)

                            # Add hit to summary
                            summary["hhsearch_hits"].append(hhsearch_hit)
                        except Exception as e:
                            logger.warning(f"Error parsing HHSearch hit: {e}")

                    logger.info(f"Processed {hit_count} HHSearch hits")
                else:
                    logger.info("No hits element found in hh_run")
            else:
                logger.info("No hh_run element found")

            # Get sequence length
            logger.info("CHECKPOINT: About to process sequence length")
            logger.info(f"ABOUT TO CALL _get_sequence_length for {pdb_id}_{chain_id}")
            sequence_length = self._get_sequence_length(pdb_id, chain_id, domain_summary_path)
            logger.info(f"AFTER CALLING _get_sequence_length, result: {sequence_length}")
            summary["sequence_length"] = sequence_length
            logger.info(f"Set sequence_length in summary to {sequence_length}")

            # Log hit counts after processing (verification)
            logger.info(f"Final hit counts in summary:")
            logger.info(f"  Chain BLAST hits: {len(summary['chain_blast_hits'])}")
            logger.info(f"  Domain BLAST hits: {len(summary['domain_blast_hits'])}")
            logger.info(f"  HHSearch hits: {len(summary['hhsearch_hits'])}")

            logger.info(f"COMPLETING _process_domain_summary with keys: {sorted(summary.keys())}")
            return summary

        except Exception as e:
            logger.error(f"EXCEPTION in _process_domain_summary: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return {"error": str(e)}

    def process_batch(self, batch_id: int, blast_only: bool = False,
                     limit: int = None, reps_only: bool = False) -> PipelineResult:
        """
        Process domains for batch with proper handling of DomainPartitionResult

        Args:
            batch_id: Batch ID
            blast_only: Whether to use only BLAST results (no HHSearch)
            limit: Maximum number of proteins to process
            reps_only: Whether to process only representative proteins

        Returns:
            PipelineResult with processing results
        """
        self.logger.info(f"Processing batch domains with ID {batch_id}")

        # Create result object
        result = PipelineResult(batch_id=batch_id)

        # Retrieve batch_info internally
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            result.success = False
            result.error = f"Batch ID {batch_id} not found"
            return result

        result.set_batch_info(batch_info)

        # Get partition component
        partition = DomainPartition(self.context)

        # Get batch path and reference from batch_info
        batch_path = batch_info.get('base_path')
        reference = batch_info.get('ref_version')

        # Run the domain partition process
        try:
            partition_results = partition.process_batch(
                batch_id,
                batch_path,
                reference,
                blast_only,
                limit,
                reps_only
            )

            # Process the results correctly based on the type
            if partition_results is None:
                self.logger.warning("process_batch returned None")
                result.success = False
                result.error = "No results returned from process_batch"
                result.partition_stats = {
                    "files_created": 0,
                    "total_proteins": 0,
                    "errors": ["No results returned from domain partition"]
                }
            elif isinstance(partition_results, list):
                # Handle list of DomainPartitionResult objects
                success_count = sum(1 for r in partition_results if r.success)
                failed_count = sum(1 for r in partition_results if not r.success)

                result.success = success_count > 0
                result.partition_stats = {
                    "files_created": success_count,
                    "total_proteins": len(partition_results),
                    "errors": [r.error for r in partition_results if not r.success and r.error]
                }
                self.logger.info(f"Created {success_count} domain partition files, {failed_count} failed")
            else:
                # Fallback for backward compatibility
                self.logger.warning(f"Unexpected result type from process_batch: {type(partition_results)}")
                result.success = bool(partition_results)
                result.partition_stats = {
                    "files_created": 0 if not result.success else 1,
                    "total_proteins": 0 if not result.success else 1
                }

        except Exception as e:
            self.logger.error(f"Exception in process_batch: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            result.success = False
            result.error = f"Exception in process_batch: {str(e)}"
            result.partition_stats = {
                "files_created": 0,
                "total_proteins": 0,
                "errors": [str(e)]
            }

        return result

    def _update_process_status(self, process_id: int, result: DomainPartitionResult) -> None:
        """Update process status in database"""
        if not process_id:
            return

        try:
            # Determine status based on result
            status = "success" if result.success else "error"
            stage = "domain_partition_complete" if result.success else "domain_partition_failed"

            # Update process status
            self.context.db.update(
                "ecod_schema.process_status",
                {
                    "current_stage": stage,
                    "status": status,
                    "error_message": result.error if not result.success else None
                },
                "id = %s",
                (process_id,)
            )

            # Add file record if successful
            if result.success and result.domain_file:
                # Check if file record already exists
                query = """
                SELECT id FROM ecod_schema.process_file
                WHERE process_id = %s AND file_type = 'domain_partition'
                """
                existing = self.context.db.execute_query(query, (process_id,))

                if existing:
                    # Update existing record
                    self.context.db.update(
                        "ecod_schema.process_file",
                        {
                            "file_path": result.domain_file,
                            "file_exists": True,
                            "file_size": os.path.getsize(result.domain_file) if os.path.exists(result.domain_file) else 0
                        },
                        "id = %s",
                        (existing[0][0],)
                    )
                else:
                    # Create new record
                    self.context.db.insert(
                        "ecod_schema.process_file",
                        {
                            "process_id": process_id,
                            "file_type": "domain_partition",
                            "file_path": result.domain_file,
                            "file_exists": True,
                            "file_size": os.path.getsize(result.domain_file) if os.path.exists(result.domain_file) else 0
                        }
                    )
        except Exception as e:
            self.logger.error(f"Error updating process status: {str(e)}")
