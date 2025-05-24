#!/usr/bin/env python3
"""
Domain partition module for the ECOD pipeline - Model-Based Implementation

This module has been completely rewritten to use the new consolidated models:
- Evidence: Unified evidence model
- DomainModel: Comprehensive domain model
- DomainPartitionResult: Enhanced result model

All dictionary-based processing has been replaced with model-based processing.
"""

import os
import re
import logging
import datetime
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Set, Union

from ecod.exceptions import PipelineError, FileOperationError
from ecod.core.context import ApplicationContext
from ecod.db import DBManager

# ONLY new gold standard models - NO legacy imports
from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult

#from ecod.utils.xml_utils import ensure_dict, ensure_list_of_dicts
from ecod.utils.path_utils import get_standardized_paths, get_all_evidence_paths, resolve_file_path
from ecod.utils.file import find_fasta_file, read_sequence_from_fasta
from ecod.utils.range_utils import parse_range


# MIGRATION NOTES:
# ================
# REMOVED IMPORTS (now using new models):
# - from ecod.models.pipeline import BlastHit, HHSearchHit, DomainSummaryModel, PipelineResult
#
# REPLACED WITH:
# - Evidence model handles all hit parsing (from_blast_xml, from_hhsearch_xml)
# - DomainModel handles all domain representation
# - DomainPartitionResult handles all result processing
#
# BENEFITS:
# - Single parsing path (no legacy model conversion)
# - Consistent data representation throughout pipeline
# - Reduced memory usage and processing overhead
# - Simplified debugging and maintenance


class DomainPartition:
    """Determine domain boundaries and classifications from search results using new models"""

    def __init__(self, context=None):
        """Initialize with configuration"""
        self.context = context or ApplicationContext()
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.partition")

        # Set default thresholds
        self.high_confidence_threshold = 0.95
        self.medium_confidence_threshold = 0.7
        self.overlap_threshold = 0.3
        self.gap_tolerance = 20

        # Initialize classification caches
        self.domain_classification_cache = {}
        self.domain_id_classification_cache = {}

    #########################################
    # Main processing methods
    #########################################

    def process_batch(self, batch_id: int, batch_path: str, reference: str,
                     blast_only: bool = False, limit: int = None,
                     reps_only: bool = False
    ) -> List[DomainPartitionResult]:
        """
        Process domains for a batch of proteins using model-based approach

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
        self.logger.info(f"Processing batch {batch_id} with model-based approach")
        results = []

        # Get proteins to process
        proteins = self._get_proteins_to_process(batch_id, limit, reps_only)

        if reps_only:
            self.logger.info("Processing representative proteins only")

        # Process each protein
        for protein in proteins:
            pdb_id = protein["pdb_id"]
            chain_id = protein["chain_id"]

            # Find domain summary file
            domain_summary_path = self._find_domain_summary(
                batch_path, pdb_id, chain_id, reference, blast_only
            )

            if not domain_summary_path:
                self.logger.warning(f"No domain summary found for {pdb_id}_{chain_id}")
                result = DomainPartitionResult(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    reference=reference,
                    success=False,
                    error="Domain summary not found",
                    domain_file=os.path.join(batch_path, "domains", f"{pdb_id}_{chain_id}.{reference}.domains.xml")
                )
                results.append(result)
                continue

            # Process domains using model-based approach
            result = self.process_protein_domains(
                pdb_id, chain_id, domain_summary_path, batch_path, reference
            )
            results.append(result)

            try:
                # Update database status
                self._update_process_status(protein["process_id"], result)
            except Exception as e:
                self.logger.error(f"Error updating process status: {str(e)}")

        self.logger.info(f"Processed {len(proteins)} proteins from batch {batch_id}")
        return results

    def process_specific_ids(self, batch_id: int, process_ids: List[int],
                            dump_dir: str, reference: str, blast_only: bool = False
        ) -> bool:
        """Process domain partition for specific process IDs using models"""

        # Get database connection
        db_config = self.context.config_manager.get_db_config()
        db = DBManager(db_config)

        # Get specific protein details by process IDs
        query = """
        SELECT ps.id, p.pdb_id, p.chain_id, ps.relative_path
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.id IN %s AND ps.batch_id = %s
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
                SELECT file_path FROM ecod_schema.process_file
                WHERE process_id = %s AND file_type = 'domain_summary'
                """
                summary_result = db.execute_query(summary_query, (process_id,))

                if not summary_result:
                    self.logger.warning(f"No domain summary found for {pdb_id}_{chain_id}")
                    continue

                summary_path = os.path.join(dump_dir, summary_result[0][0])
                if not os.path.exists(summary_path):
                    self.logger.warning(f"Domain summary file not found: {summary_path}")
                    continue

                # Process using model-based approach
                self.logger.info(f"Processing {pdb_id}_{chain_id} (process_id: {process_id})")

                result = self.process_protein_domains(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    domain_summary_path=summary_path,
                    output_dir=dump_dir,
                    reference=reference
                )

                if result.success:
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

                    if result.domain_file:
                        relative_path = os.path.relpath(result.domain_file, dump_dir)
                        self.register_domain_file(process_id, relative_path)
                        self.logger.info(f"Successfully processed {pdb_id}_{chain_id}")
                else:
                    self.logger.error(f"Failed to process {pdb_id}_{chain_id}: {result.error}")

                    db.update(
                        "ecod_schema.process_status",
                        {
                            "current_stage": "domain_partition_failed",
                            "status": "error",
                            "error_message": str(result.error)[:500]
                        },
                        "id = %s",
                        (process_id,)
                    )

            except Exception as e:
                self.logger.error(f"Error processing {pdb_id}_{chain_id}: {e}")
                import traceback
                self.logger.error(traceback.format_exc())

        self.logger.info(f"Processed {success_count}/{len(rows)} proteins successfully")
        return success_count > 0

    def process_protein_domains(self, pdb_id: str, chain_id: str,
                               domain_summary_path: str,
                               output_dir: str,
                               reference: str = "develop291"
        ) -> DomainPartitionResult:
        """Process domains for a protein using pure model-based approach WITH VALIDATION"""

        start_time = datetime.datetime.now()

        try:
            self.logger.info(f"Processing domains for {pdb_id}_{chain_id} with validation")

            # Create result with output file path
            domains_dir = os.path.join(output_dir, "domains")
            os.makedirs(domains_dir, exist_ok=True)

            domain_file = os.path.join(
                domains_dir,
                f"{pdb_id}_{chain_id}.{reference}.domains.xml"
            )

            result = DomainPartitionResult(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference=reference,
                domain_file=domain_file
            )

            # Parse domain summary into evidence models WITH VALIDATION
            summary_data = self._parse_domain_summary(domain_summary_path)

            if "error" in summary_data:
                result.success = False
                result.error = summary_data["error"]
                result.is_unclassified = True
                result.save()
                return result

            # Get sequence length and check for peptide
            sequence_length = summary_data.get("sequence_length", 0)
            result.sequence_length = sequence_length

            is_peptide = summary_data.get("is_peptide", False) or sequence_length < 20
            if is_peptide:
                self.logger.info(f"Chain {pdb_id}_{chain_id} classified as peptide")
                result.is_peptide = True
                result.is_unclassified = True
                result.success = True
                result.save()
                return result

            # Extract evidence from summary WITH VALIDATION
            all_evidence = self._extract_evidence_from_summary(summary_data)

            if not all_evidence:
                self.logger.warning(f"No valid evidence found for {pdb_id}_{chain_id}")
                result.is_unclassified = True
                result.success = True
                result.save()
                return result

            # Identify domain boundaries from evidence WITH VALIDATION
            self.logger.info(f"Identifying domain boundaries from {len(all_evidence)} validated evidence items")
            domain_models = self._identify_domain_boundaries(all_evidence, sequence_length, pdb_id, chain_id)

            if not domain_models:
                self.logger.warning(f"No valid domains identified for {pdb_id}_{chain_id}")
                result.is_unclassified = True
                result.success = True
                result.save()
                return result

            # Assign classifications to domains WITH VALIDATION
            self.logger.info(f"Assigning classifications to {len(domain_models)} validated domains")
            self._assign_domain_classifications(domain_models)

            # Set domains in result (domains are already validated)
            result.domains = domain_models
            result.is_classified = True
            result.success = True

            # Calculate processing time
            end_time = datetime.datetime.now()
            result.processing_time = (end_time - start_time).total_seconds()

            # Save result
            if result.save():
                self.logger.info(f"Successfully processed {pdb_id}_{chain_id} with {len(domain_models)} validated domains")
            else:
                result.error = "Failed to save domain partition file"
                result.success = False

            return result

        except Exception as e:
            self.logger.error(f"Error processing {pdb_id}_{chain_id}: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())

            result = DomainPartitionResult(
                pdb_id=pdb_id,
                chain_id=chain_id,
                reference=reference,
                success=False,
                error=str(e),
                domain_file=os.path.join(output_dir, "domains", f"{pdb_id}_{chain_id}.{reference}.domains.xml")
            )
            return result

    #########################################
    # Evidence processing methods
    #########################################

    def _parse_domain_summary(self, domain_summary_path: str) -> Dict[str, Any]:
        """Parse domain summary file directly into Evidence objects WITH VALIDATION"""

        if not os.path.exists(domain_summary_path):
            return {"error": "File not found"}

        try:
            tree = ET.parse(domain_summary_path)
            root = tree.getroot()

            # Get basic information
            summary_elem = root.find("blast_summ")
            if summary_elem is None:
                return {"error": "Invalid domain summary format"}

            pdb_id = summary_elem.get("pdb", "")
            chain_id = summary_elem.get("chain", "")
            is_peptide = summary_elem.get("is_peptide", "false").lower() == "true"

            # Create summary dictionary with Evidence objects
            summary = {
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "is_peptide": is_peptide,
                "chain_blast_evidence": [],
                "domain_blast_evidence": [],
                "hhsearch_evidence": []
            }

            # Parse chain BLAST hits directly to Evidence WITH VALIDATION
            chain_blast_run = root.find("chain_blast_run")
            if chain_blast_run is not None:
                hits_elem = chain_blast_run.find("hits")
                if hits_elem is not None:
                    for hit_elem in hits_elem.findall("hit"):
                        try:
                            evidence = Evidence.from_blast_xml(hit_elem, "chain_blast")

                            # VALIDATE EVIDENCE
                            if self._validate_evidence(evidence, "chain_blast_xml_parsing"):
                                summary["chain_blast_evidence"].append(evidence)
                                self.logger.debug(f"Added valid chain BLAST evidence: {evidence.source_id}")
                            else:
                                self.logger.warning(f"Skipped invalid chain BLAST evidence from {domain_summary_path}")

                        except Exception as e:
                            self.logger.warning(f"Error parsing chain BLAST hit: {e}")

            # Parse domain BLAST hits directly to Evidence WITH VALIDATION
            blast_run = root.find("blast_run")
            if blast_run is not None:
                hits_elem = blast_run.find("hits")
                if hits_elem is not None:
                    for hit_elem in hits_elem.findall("hit"):
                        try:
                            evidence = Evidence.from_blast_xml(hit_elem, "domain_blast")

                            # VALIDATE EVIDENCE
                            if self._validate_evidence(evidence, "domain_blast_xml_parsing"):
                                summary["domain_blast_evidence"].append(evidence)
                                self.logger.debug(f"Added valid domain BLAST evidence: {evidence.source_id}")
                            else:
                                self.logger.warning(f"Skipped invalid domain BLAST evidence from {domain_summary_path}")

                        except Exception as e:
                            self.logger.warning(f"Error parsing domain BLAST hit: {e}")

            # Parse HHSearch hits directly to Evidence WITH VALIDATION
            hh_run = root.find("hh_run")
            if hh_run is not None:
                hits_elem = hh_run.find("hits")
                if hits_elem is not None:
                    for hit_elem in hits_elem.findall("hit"):
                        try:
                            evidence = Evidence.from_hhsearch_xml(hit_elem)

                            # VALIDATE EVIDENCE
                            if self._validate_evidence(evidence, "hhsearch_xml_parsing"):
                                summary["hhsearch_evidence"].append(evidence)
                                self.logger.debug(f"Added valid HHSearch evidence: {evidence.source_id}")
                            else:
                                self.logger.warning(f"Skipped invalid HHSearch evidence from {domain_summary_path}")

                        except Exception as e:
                            self.logger.warning(f"Error parsing HHSearch hit: {e}")

            # Get sequence length
            sequence_length = self._get_sequence_length(pdb_id, chain_id, domain_summary_path)
            summary["sequence_length"] = sequence_length

            # Log validation summary
            total_evidence = (len(summary["chain_blast_evidence"]) +
                             len(summary["domain_blast_evidence"]) +
                             len(summary["hhsearch_evidence"]))

            self.logger.info(f"Parsed and validated {total_evidence} evidence items from {domain_summary_path}")

            return summary

        except Exception as e:
            self.logger.error(f"Error parsing domain summary: {str(e)}")
            return {"error": str(e)}


    def _extract_evidence_from_summary(self, summary_data: Dict[str, Any]) -> List[Evidence]:
        """Extract Evidence objects from parsed summary data WITH VALIDATION"""

        all_evidence = []
        validation_stats = {"valid": 0, "invalid": 0, "classification_enriched": 0}

        # Process HHSearch evidence (highest priority) - already Evidence objects
        for evidence in summary_data.get("hhsearch_evidence", []):
            try:
                # RE-VALIDATE (in case evidence was modified)
                if not self._validate_evidence(evidence, "hhsearch_processing"):
                    validation_stats["invalid"] += 1
                    continue

                # Get classification from database if available
                if evidence.domain_id:
                    classification = self._get_domain_classification_by_id(evidence.domain_id)
                    if classification:
                        evidence.t_group = classification.get("t_group")
                        evidence.h_group = classification.get("h_group")
                        evidence.x_group = classification.get("x_group")
                        evidence.a_group = classification.get("a_group")
                        validation_stats["classification_enriched"] += 1

                # VALIDATE AGAIN after classification enrichment
                if self._validate_evidence(evidence, "hhsearch_post_classification"):
                    all_evidence.append(evidence)
                    validation_stats["valid"] += 1
                else:
                    validation_stats["invalid"] += 1

            except Exception as e:
                self.logger.warning(f"Error processing HHSearch evidence: {e}")
                validation_stats["invalid"] += 1

        # Process domain BLAST evidence - already Evidence objects
        for evidence in summary_data.get("domain_blast_evidence", []):
            try:
                # RE-VALIDATE
                if not self._validate_evidence(evidence, "domain_blast_processing"):
                    validation_stats["invalid"] += 1
                    continue

                # Get classification from database if available
                if evidence.domain_id:
                    classification = self._get_domain_classification_by_id(evidence.domain_id)
                    if classification:
                        evidence.t_group = classification.get("t_group")
                        evidence.h_group = classification.get("h_group")
                        evidence.x_group = classification.get("x_group")
                        evidence.a_group = classification.get("a_group")
                        validation_stats["classification_enriched"] += 1

                # VALIDATE AGAIN after classification enrichment
                if self._validate_evidence(evidence, "domain_blast_post_classification"):
                    all_evidence.append(evidence)
                    validation_stats["valid"] += 1
                else:
                    validation_stats["invalid"] += 1

            except Exception as e:
                self.logger.warning(f"Error processing domain BLAST evidence: {e}")
                validation_stats["invalid"] += 1

        # Process chain BLAST evidence (map to domains) - already Evidence objects
        for evidence in summary_data.get("chain_blast_evidence", []):
            try:
                # RE-VALIDATE
                if not self._validate_evidence(evidence, "chain_blast_processing"):
                    validation_stats["invalid"] += 1
                    continue

                # Get reference domains for this chain
                pdb_id = evidence.extra_attributes.get("pdb_id", "")
                chain_id = evidence.extra_attributes.get("chain_id", "")

                if pdb_id and chain_id:
                    ref_domains = self._get_reference_chain_domains(f"{pdb_id}_{chain_id}")

                    if ref_domains:
                        # Map each reference domain to query
                        mapped_evidence = self._map_chain_blast_to_domains(evidence, ref_domains)

                        # VALIDATE each mapped evidence
                        for mapped_ev in mapped_evidence:
                            if self._validate_evidence(mapped_ev, "chain_blast_mapping"):
                                all_evidence.append(mapped_ev)
                                validation_stats["valid"] += 1
                            else:
                                validation_stats["invalid"] += 1
                    else:
                        # Use as-is if no reference domains found
                        if self._validate_evidence(evidence, "chain_blast_no_mapping"):
                            all_evidence.append(evidence)
                            validation_stats["valid"] += 1
                        else:
                            validation_stats["invalid"] += 1
                else:
                    if self._validate_evidence(evidence, "chain_blast_no_ids"):
                        all_evidence.append(evidence)
                        validation_stats["valid"] += 1
                    else:
                        validation_stats["invalid"] += 1

            except Exception as e:
                self.logger.warning(f"Error processing chain BLAST evidence: {e}")
                validation_stats["invalid"] += 1

        # Log validation statistics
        self.logger.info(f"Evidence validation: {validation_stats['valid']} valid, "
                        f"{validation_stats['invalid']} invalid, "
                        f"{validation_stats['classification_enriched']} classification-enriched")

        return all_evidence


    def _map_chain_blast_to_domains(self, chain_evidence: Evidence, ref_domains: List[Dict[str, Any]]) -> List[Evidence]:
        """Map chain BLAST evidence to domain evidence using reference domains WITH VALIDATION"""

        mapped_evidence = []

        try:
            # Parse query and hit ranges
            evidence_ranges = self._get_evidence_ranges(chain_evidence, "chain_blast_mapping")
            query_ranges = evidence_ranges["query"]
            hit_ranges = evidence_ranges["hit"]

            if not query_ranges or not hit_ranges:
                self.logger.warning("Invalid ranges in chain BLAST mapping")
                return []

            # For simplicity, use first range pair
            q_start, q_end = query_ranges[0]
            h_start, h_end = hit_ranges[0]

            q_length = q_end - q_start + 1
            h_length = h_end - h_start + 1

            for ref_domain in ref_domains:
                ref_start = ref_domain.get("start", 0)
                ref_end = ref_domain.get("end", 0)

                if ref_start <= 0 or ref_end <= 0:
                    self.logger.debug(f"Skipping reference domain with invalid coordinates: {ref_start}-{ref_end}")
                    continue

                # Check if reference domain overlaps with hit range
                if max(h_start, ref_start) <= min(h_end, ref_end):
                    # Map domain to query coordinates
                    ratio = q_length / h_length if h_length > 0 else 1.0

                    mapped_start = max(1, round(q_start + (ref_start - h_start) * ratio))
                    mapped_end = min(q_end, round(q_start + (ref_end - h_start) * ratio))

                    if mapped_start <= mapped_end:
                        try:
                            # Create evidence for mapped domain
                            evidence = Evidence(
                                type="chain_blast",
                                source_id=ref_domain.get("domain_id", ""),
                                domain_id=ref_domain.get("domain_id", ""),
                                query_range=f"{mapped_start}-{mapped_end}",
                                hit_range=f"{ref_start}-{ref_end}",
                                evalue=chain_evidence.evalue,
                                t_group=ref_domain.get("t_group"),
                                h_group=ref_domain.get("h_group"),
                                x_group=ref_domain.get("x_group"),
                                a_group=ref_domain.get("a_group"),
                                extra_attributes=chain_evidence.extra_attributes.copy()
                            )

                            # VALIDATE mapped evidence
                            if self._validate_evidence(evidence, "chain_blast_domain_mapping"):
                                mapped_evidence.append(evidence)
                                self.logger.debug(f"Mapped chain BLAST to domain {ref_domain.get('domain_id', 'unknown')}: {mapped_start}-{mapped_end}")
                            else:
                                self.logger.warning(f"Invalid mapped evidence for domain {ref_domain.get('domain_id', 'unknown')}")

                        except Exception as e:
                            self.logger.warning(f"Error creating mapped evidence: {e}")

        except Exception as e:
            self.logger.warning(f"Error mapping chain BLAST to domains: {e}")

        return mapped_evidence

    #########################################
    # Domain boundary identification
    #########################################

    def _identify_domain_boundaries(self, evidence_list: List[Evidence],
                                   sequence_length: int, pdb_id: str, chain_id: str) -> List[DomainModel]:
        """Identify domain boundaries from Evidence objects WITH VALIDATION"""

        # VALIDATE all input evidence first
        valid_evidence = []
        for evidence in evidence_list:
            if self._validate_evidence(evidence, "domain_boundary_identification"):
                valid_evidence.append(evidence)
            else:
                self.logger.warning(f"Skipping invalid evidence in boundary identification: {evidence.source_id}")

        if not valid_evidence:
            self.logger.warning("No valid evidence for domain boundary identification")
            return []

        # Group evidence by query ranges to identify domain candidates
        domain_candidates = {}

        for evidence in valid_evidence:
            if not evidence.query_range:
                continue

            try:
                evidence_ranges = self._get_evidence_ranges(evidence, "domain_boundary_identification")
                ranges = evidence_ranges["query"]
                for start, end in ranges:
                    # Validate coordinates
                    if start <= 0 or end <= 0 or start > end:
                        self.logger.warning(f"Invalid range {start}-{end} in evidence {evidence.source_id}")
                        continue

                    # Create a position-based key for grouping
                    position_key = (start // 50) * 50  # Group by 50-residue windows

                    if position_key not in domain_candidates:
                        domain_candidates[position_key] = []

                    domain_candidates[position_key].append({
                        "evidence": evidence,
                        "start": start,
                        "end": end
                    })
            except Exception as e:
                self.logger.warning(f"Error parsing range {evidence.query_range}: {e}")

        # Convert candidates to domain models WITH VALIDATION
        domain_models = []

        for position_key, candidates in domain_candidates.items():
            if len(candidates) < 1:  # Require at least one piece of evidence
                continue

            try:
                # Find consensus boundaries
                starts = [c["start"] for c in candidates]
                ends = [c["end"] for c in candidates]

                # Use median for consensus
                consensus_start = sorted(starts)[len(starts) // 2]
                consensus_end = sorted(ends)[len(ends) // 2]

                # Validate consensus boundaries
                if consensus_start <= 0 or consensus_end <= 0 or consensus_start > consensus_end:
                    self.logger.warning(f"Invalid consensus boundaries: {consensus_start}-{consensus_end}")
                    continue

                if sequence_length > 0 and (consensus_start > sequence_length or consensus_end > sequence_length):
                    self.logger.warning(f"Consensus boundaries exceed sequence length: {consensus_start}-{consensus_end} > {sequence_length}")
                    continue

                # Get best evidence (highest confidence)
                best_candidate = max(candidates, key=lambda c: c["evidence"].confidence or 0)
                best_evidence = best_candidate["evidence"]

                # Create domain model with robust ID generation
                domain_id = f"{pdb_id}_{chain_id}_d{consensus_start}_{consensus_end}"

                # Create domain model
                domain = DomainModel(
                    id=domain_id,
                    start=consensus_start,
                    end=consensus_end,
                    range=f"{consensus_start}-{consensus_end}",
                    source=best_evidence.type,
                    source_id=best_evidence.source_id or best_evidence.domain_id,
                    t_group=best_evidence.t_group,
                    h_group=best_evidence.h_group,
                    x_group=best_evidence.x_group,
                    a_group=best_evidence.a_group,
                    evidence=[c["evidence"] for c in candidates]
                )

                # VALIDATE domain model
                if self._validate_domain_model(domain, "domain_boundary_identification"):
                    domain_models.append(domain)
                    self.logger.debug(f"Created valid domain model: {domain_id} ({consensus_start}-{consensus_end})")
                else:
                    self.logger.warning(f"Invalid domain model created: {domain_id}")

            except Exception as e:
                self.logger.warning(f"Error creating domain model for position {position_key}: {e}")

        # Resolve overlaps (domains are already validated)
        final_domains = self._resolve_domain_overlaps(domain_models, sequence_length)

        # Sort by position
        final_domains.sort(key=lambda d: d.start)

        # Final validation of all domains
        validated_domains = []
        for domain in final_domains:
            if self._validate_domain_model(domain, "final_domain_validation"):
                validated_domains.append(domain)
            else:
                self.logger.warning(f"Domain failed final validation: {domain.id}")

        self.logger.info(f"Identified and validated {len(validated_domains)} domains from {len(evidence_list)} evidence items")
        return validated_domains

    def _resolve_domain_overlaps(self, domains: List[DomainModel], sequence_length: int) -> List[DomainModel]:
        """Resolve overlapping domains by confidence and protection status"""

        if len(domains) <= 1:
            return domains

        # Sort by confidence (descending) and protection status
        sorted_domains = sorted(
            domains,
            key=lambda d: (d.protected, d.confidence, d.end - d.start + 1),
            reverse=True
        )

        # Track covered positions
        covered_positions = set()
        final_domains = []

        for domain in sorted_domains:
            domain_positions = domain.get_positions()
            overlap = domain_positions.intersection(covered_positions)
            overlap_pct = len(overlap) / len(domain_positions) if domain_positions else 0

            # Include domain if minimal overlap or it's protected
            if overlap_pct < self.overlap_threshold or domain.protected:
                final_domains.append(domain)
                covered_positions.update(domain_positions)

                self.logger.debug(f"Added domain {domain.id} ({domain.range}) with {overlap_pct:.2f} overlap")

        return final_domains

    def _parse_range_with_fallback(self, range_string: str, context: str = "unknown") -> List[Tuple[int, int]]:
        """
        Parse range string using model methods when available, with fallback to utility

        Args:
            range_string: Range string like "1-100,150-200"
            context: Context for logging

        Returns:
            List of (start, end) tuples
        """
        if not range_string:
            return []

        try:
            # First check if we have a model that can parse ranges
            # This would be if Evidence or DomainModel objects have range parsing

            # For now, use the external utility but add validation and logging
            from ecod.utils.range_utils import parse_range

            ranges = parse_range(range_string)

            # Validate parsed ranges
            valid_ranges = []
            for start, end in ranges:
                if start > 0 and end > 0 and start <= end:
                    valid_ranges.append((start, end))
                else:
                    self.logger.warning(f"Invalid range segment {start}-{end} in {context}")

            if not valid_ranges and ranges:
                self.logger.warning(f"No valid ranges found in '{range_string}' for {context}")

            return valid_ranges

        except Exception as e:
            self.logger.warning(f"Error parsing range '{range_string}' in {context}: {e}")
            return []


    # =============================================================================
    # ENHANCED RANGE PARSING WITH MODEL INTEGRATION
    # =============================================================================

    def _get_evidence_ranges(self, evidence: Evidence, context: str = "unknown") -> Dict[str, List[Tuple[int, int]]]:
        """
        Get query and hit ranges from Evidence object using best available method

        Args:
            evidence: Evidence object
            context: Context for logging

        Returns:
            Dict with 'query' and 'hit' range lists
        """
        result = {"query": [], "hit": []}

        try:
            # Check if Evidence has built-in range parsing methods
            if hasattr(evidence, 'get_query_ranges'):
                result["query"] = evidence.get_query_ranges()
            elif hasattr(evidence, 'parse_query_range'):
                result["query"] = evidence.parse_query_range()
            else:
                # Fallback to utility
                result["query"] = self._parse_range_with_fallback(evidence.query_range, f"{context}_query")

            if hasattr(evidence, 'get_hit_ranges'):
                result["hit"] = evidence.get_hit_ranges()
            elif hasattr(evidence, 'parse_hit_range'):
                result["hit"] = evidence.parse_hit_range()
            else:
                # Fallback to utility
                result["hit"] = self._parse_range_with_fallback(evidence.hit_range, f"{context}_hit")

        except Exception as e:
            self.logger.warning(f"Error getting ranges from evidence in {context}: {e}")

        return result


    def _get_domain_ranges(self, domain: DomainModel, context: str = "unknown") -> List[Tuple[int, int]]:
        """
        Get ranges from DomainModel using best available method

        Args:
            domain: DomainModel object
            context: Context for logging

        Returns:
            List of (start, end) tuples
        """
        try:
            # Check if DomainModel has built-in range parsing
            if hasattr(domain, 'get_range_segments'):
                segments = domain.get_range_segments()
                return [(seg.start, seg.end) if hasattr(seg, 'start') else seg for seg in segments]
            elif hasattr(domain, 'parse_range'):
                return domain.parse_range()
            else:
                # Fallback to utility
                return self._parse_range_with_fallback(domain.range, context)

        except Exception as e:
            self.logger.warning(f"Error getting ranges from domain {getattr(domain, 'id', 'unknown')} in {context}: {e}")
            return []

    #########################################
    # Classification assignment
    #########################################

     def _assign_domain_classifications(self, domains: List[DomainModel]) -> None:
        """Assign classifications to domains using evidence WITH VALIDATION"""

        for domain in domains:
            try:
                # RE-VALIDATE domain before processing
                if not self._validate_domain_model(domain, "classification_assignment"):
                    self.logger.warning(f"Skipping classification for invalid domain: {domain.id}")
                    continue

                self.logger.debug(f"Assigning classification to domain {domain.id}")

                # If domain already has full classification, skip
                if domain.is_fully_classified():
                    continue

                # Find best evidence with classification
                best_evidence = None
                best_score = 0

                for evidence in domain.evidence:
                    # VALIDATE evidence before using
                    if not self._validate_evidence(evidence, f"classification_evidence_{domain.id}"):
                        continue

                    if not any([evidence.t_group, evidence.h_group, evidence.x_group, evidence.a_group]):
                        continue

                    # Score evidence by type and confidence
                    type_weights = {
                        "hhsearch": 3.0,
                        "domain_blast": 2.5,
                        "chain_blast": 2.0,
                        "blast": 1.5
                    }

                    weight = type_weights.get(evidence.type, 1.0)
                    confidence = evidence.confidence or 0.0
                    score = confidence * weight

                    if score > best_score:
                        best_score = score
                        best_evidence = evidence

                # Apply classification from best evidence
                if best_evidence:
                    # Only update if current value is None/empty
                    if not domain.t_group and best_evidence.t_group:
                        domain.t_group = best_evidence.t_group
                    if not domain.h_group and best_evidence.h_group:
                        domain.h_group = best_evidence.h_group
                    if not domain.x_group and best_evidence.x_group:
                        domain.x_group = best_evidence.x_group
                    if not domain.a_group and best_evidence.a_group:
                        domain.a_group = best_evidence.a_group

                    # VALIDATE domain after classification update
                    if self._validate_domain_model(domain, "post_classification"):
                        self.logger.debug(f"Applied classification from {best_evidence.type} evidence: "
                                        f"{domain.get_classification_level()}")
                    else:
                        self.logger.warning(f"Domain became invalid after classification update: {domain.id}")
                else:
                    self.logger.debug(f"No valid classification evidence found for domain {domain.id}")

            except Exception as e:
                self.logger.warning(f"Error assigning classification to domain {getattr(domain, 'id', 'unknown')}: {e}")


    #########################################
    # Database and reference data methods
    #########################################

    def _get_domain_classification_by_id(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """Get domain classification from database by domain ID with caching"""

        if not domain_id:
            return None

        # Check cache first
        if domain_id in self.domain_id_classification_cache:
            return self.domain_id_classification_cache[domain_id]

        # Query database
        try:

            query = """
            SELECT t_group, h_group, x_group, a_group,
                   is_manual_rep, is_f70, is_f40, is_f99
            FROM pdb_analysis.domain
            WHERE domain_id = %s
            """

            rows = self.context.db.execute_dict_query(query, (domain_id,))
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
                self.domain_id_classification_cache[domain_id] = classification
                return classification
        except Exception as e:
            self.logger.error(f"Error getting classification for {domain_id}: {e}")

        return None

    def _get_reference_chain_domains(self, source_id: str) -> List[Dict[str, Any]]:
        """Get reference domain information for a chain"""

        try:

            # Parse source_id to get pdb_id and chain_id
            if "_" in source_id:
                pdb_id, chain_id = source_id.split("_", 1)
            else:
                return []

            query = """
            SELECT d.domain_id, d.range, d.t_group, d.h_group, d.x_group, d.a_group
            FROM pdb_analysis.domain d
            JOIN pdb_analysis.protein p ON d.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            ORDER BY d.start_position
            """

            rows = self.context.db.execute_dict_query(query, (pdb_id, chain_id))

            domains = []
            for row in rows:
                # Parse range to get start/end positions
                range_str = row["range"]
                ranges = self._parse_range_with_fallback(range_str, f"reference_domain_{row['domain_id']}")

                if ranges:
                    start, end = ranges[0]  # Use first range segment
                    domain = {
                        "domain_id": row["domain_id"],
                        "range": range_str,
                        "start": start,
                        "end": end,
                        "t_group": row["t_group"],
                        "h_group": row["h_group"],
                        "x_group": row["x_group"],
                        "a_group": row["a_group"]
                    }
                    domains.append(domain)

            self.logger.debug(f"Found {len(domains)} reference domains for {source_id}")
            return domains

        except Exception as e:
            self.logger.error(f"Error getting reference domains for {source_id}: {e}")
            return []

    def _get_proteins_to_process(self, batch_id: int, limit: int = None, reps_only: bool = False) -> List[Dict[str, Any]]:
        """Get proteins to process for domain partitioning from a batch"""

        query = """
        SELECT p.id as protein_id, p.pdb_id, p.chain_id, ps.id as process_id, ps.is_representative
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN ecod_schema.process_file pf_summ ON (
            pf_summ.process_id = ps.id AND
            pf_summ.file_type = 'domain_summary' AND
            pf_summ.file_exists = TRUE
        )
        LEFT JOIN ecod_schema.process_file pf_part ON (
            pf_part.process_id = ps.id AND
            pf_part.file_type = 'domain_partition' AND
            pf_part.file_exists = TRUE
        )
        WHERE ps.batch_id = %s
          AND pf_summ.id IS NOT NULL
          AND pf_part.id IS NULL
        """

        params = [batch_id]

        if reps_only:
            query += " AND ps.is_representative = TRUE"

        query += " ORDER BY p.id"

        if limit is not None:
            query += f" LIMIT {limit}"

        try:
            results = self.context.db.execute_dict_query(query, tuple(params))
            self.logger.info(f"Found {len(results)} proteins to process for batch {batch_id}")
            return results
        except Exception as e:
            self.logger.error(f"Error fetching proteins to process: {str(e)}")
            return []

    #########################################
    # Utility methods
    #########################################

    def _find_domain_summary(self, batch_path: str, pdb_id: str, chain_id: str,
                           reference: str, blast_only: bool = False) -> str:
        """Find domain summary file using path utilities"""

        file_type = 'blast_only_summary' if blast_only else 'domain_summary'

        try:
            evidence_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, reference)

            if file_type in evidence_paths and evidence_paths[file_type]['exists_at']:
                summary_path = evidence_paths[file_type]['exists_at']
                self.logger.debug(f"Found domain summary: {summary_path}")
                return summary_path

        except Exception as e:
            self.logger.warning(f"Error using path_utils to find domain summary: {e}")

        return ""

    def _get_sequence_length(self, pdb_id: str, chain_id: str, domain_summary_path: str) -> int:
        """Get sequence length for a protein chain"""

        # Try database first
        try:
            query = "SELECT length FROM ecod_schema.protein WHERE pdb_id = %s AND chain_id = %s LIMIT 1"
            results = self.context.db.execute_query(query, (pdb_id, chain_id))
            if results and results[0][0]:
                return results[0][0]
        except Exception as e:
            self.logger.warning(f"Error getting sequence length from database: {e}")

        # Try XML file
        try:
            tree = ET.parse(domain_summary_path)
            root = tree.getroot()
            query_len_elem = root.find(".//query_len")
            if query_len_elem is not None and query_len_elem.text:
                return int(query_len_elem.text.strip())
        except Exception as e:
            self.logger.warning(f"Error getting sequence length from XML: {e}")

        return 0

    def _update_process_status(self, process_id: int, result: DomainPartitionResult) -> None:
        """Update process status in database"""

        if not process_id:
            return

        try:
            status = "success" if result.success else "error"
            stage = "domain_partition_complete" if result.success else "domain_partition_failed"

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
                self.register_domain_file(process_id, result.domain_file)

        except Exception as e:
            self.logger.error(f"Error updating process status: {str(e)}")

    def register_domain_file(self, process_id: int, file_path: str, db: DBManager = None) -> None:
        """Register domain partition file in database"""

        if db is None:
            db_config = self.context.config_manager.get_db_config()
            db = DBManager(db_config)

        try:
            # Check if record already exists
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = 'domain_partition'
            """

            existing = db.execute_query(query, (process_id,))

            file_size = os.path.getsize(file_path) if os.path.exists(file_path) else 0

            if existing:
                # Update existing record
                db.update(
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
                # Insert new record
                db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": "domain_partition",
                        "file_path": file_path,
                        "file_exists": True,
                        "file_size": file_size
                    }
                )
        except Exception as e:
            self.logger.warning(f"Error registering domain file: {e}")

    def _validate_evidence(self, evidence: Evidence, context: str = "unknown") -> bool:
        """
        Validate Evidence object with comprehensive error handling

        Args:
            evidence: Evidence object to validate
            context: Context string for logging (e.g., "blast_xml_parsing")

        Returns:
            bool: True if valid, False if invalid
        """
        try:
            # Check if Evidence has a validate method
            if hasattr(evidence, 'validate'):
                evidence.validate()
                return True

            # Manual validation if no validate method
            if not evidence.type:
                self.logger.warning(f"Evidence missing type in {context}")
                return False

            if not evidence.source_id and not evidence.domain_id:
                self.logger.warning(f"Evidence missing source_id and domain_id in {context}")
                return False

            # Validate confidence is in valid range
            if evidence.confidence is not None and not (0.0 <= evidence.confidence <= 1.0):
                self.logger.warning(f"Evidence confidence {evidence.confidence} out of range in {context}")
                return False

            # Validate evalue is positive if present
            if evidence.evalue is not None and evidence.evalue < 0:
                self.logger.warning(f"Evidence evalue {evidence.evalue} is negative in {context}")
                return False

            # Validate probability is in valid range if present
            if evidence.probability is not None and not (0.0 <= evidence.probability <= 100.0):
                self.logger.warning(f"Evidence probability {evidence.probability} out of range in {context}")
                return False

            return True

        except Exception as e:
            self.logger.warning(f"Evidence validation failed in {context}: {str(e)}")
            return False

    def _validate_domain_model(self, domain: DomainModel, context: str = "unknown") -> bool:
        """
        Validate DomainModel object with comprehensive error handling

        Args:
            domain: DomainModel object to validate
            context: Context string for logging

        Returns:
            bool: True if valid, False if invalid
        """
        try:
            # Check if DomainModel has a validate method
            if hasattr(domain, 'validate'):
                domain.validate()
                return True

            # Manual validation if no validate method
            if not domain.id:
                self.logger.warning(f"Domain missing ID in {context}")
                return False

            if domain.start <= 0 or domain.end <= 0:
                self.logger.warning(f"Domain {domain.id} has invalid coordinates: {domain.start}-{domain.end} in {context}")
                return False

            if domain.start > domain.end:
                self.logger.warning(f"Domain {domain.id} has start > end: {domain.start}-{domain.end} in {context}")
                return False

            if domain.confidence is not None and not (0.0 <= domain.confidence <= 1.0):
                self.logger.warning(f"Domain {domain.id} confidence {domain.confidence} out of range in {context}")
                return False

            return True

        except Exception as e:
            self.logger.warning(f"Domain validation failed for {domain.id if hasattr(domain, 'id') else 'unknown'} in {context}: {str(e)}")
            return False

    #########################################
    # Backward compatibility methods (deprecated)
    #########################################

    def partition_domains(self, pdb_id: str, chain_id: str, dump_dir: str, input_mode: str,
                         reference: str, blast_only: bool = False) -> Dict[str, Any]:
        """
        DEPRECATED: Backward compatibility wrapper

        This method is deprecated. Use process_protein_domains instead.
        Now uses the new Evidence-based parsing internally.
        """
        import warnings
        warnings.warn(
            "partition_domains is deprecated. Use process_protein_domains instead.",
            DeprecationWarning,
            stacklevel=2
        )

        # Find domain summary
        summary_path = self._find_domain_summary(dump_dir, pdb_id, chain_id, reference, blast_only)
        if not summary_path:
            return {
                "success": False,
                "file_path": "",
                "error": "Domain summary not found",
                "domains": [],
                "stats": {}
            }

        # Use new model-based method
        result = self.process_protein_domains(
            pdb_id=pdb_id,
            chain_id=chain_id,
            domain_summary_path=summary_path,
            output_dir=dump_dir,
            reference=reference
        )

        # Convert to legacy format
        legacy_domains = []
        for domain in result.domains:
            if hasattr(domain, 'to_dict'):
                legacy_domains.append(domain.to_dict())
            elif isinstance(domain, dict):
                legacy_domains.append(domain)
            else:
                # Fallback conversion
                legacy_domains.append({
                    "id": getattr(domain, 'id', 'unknown'),
                    "start": getattr(domain, 'start', 0),
                    "end": getattr(domain, 'end', 0),
                    "range": getattr(domain, 'range', ''),
                    "source": getattr(domain, 'source', ''),
                    "confidence": getattr(domain, 'confidence', 0.0)
                })

        return {
            "success": result.success,
            "file_path": result.domain_file or "",
            "error": result.error,
            "domains": legacy_domains,
            "stats": result.get_summary_stats() if hasattr(result, 'get_summary_stats') else {}
        }

    # Alias for the deprecated method name
    process_domains = process_protein_domains
