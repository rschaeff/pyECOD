#!/usr/bin/env python3
"""
Core domain partitioning processor.

This module contains the core algorithms for identifying domain boundaries,
resolving overlaps, and assigning classifications.
"""

import logging
from typing import Dict, Any, List, Optional, Tuple, Set
from collections import defaultdict

from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult
from ecod.db import DBManager

from .models import (
    PartitionOptions, EvidenceGroup, DomainCandidate,
    PartitionContext, ValidationLevel, ClassificationCache
)
from .analyzer import EvidenceAnalyzer


class PartitionProcessor:
    """
    Core processor for domain partitioning.

    This class implements the main algorithms for:
    - Identifying domain boundaries from evidence
    - Resolving overlapping domains
    - Assigning classifications to domains
    """

    def __init__(self, options: PartitionOptions,
                 analyzer: Optional[EvidenceAnalyzer] = None,
                 db_manager: Optional[DBManager] = None):
        """
        Initialize the partition processor.

        Args:
            options: Partitioning options
            analyzer: Evidence analyzer (created if not provided)
            db_manager: Database manager for classification lookups
        """
        self.options = options
        self.analyzer = analyzer or EvidenceAnalyzer(options)
        self.db = db_manager
        self.logger = logging.getLogger(__name__)

        # Classification cache
        self.classification_cache = ClassificationCache()

        # Processing statistics
        self.stats = {
            'domains_identified': 0,
            'overlaps_resolved': 0,
            'classifications_assigned': 0,
            'domains_validated': 0
        }

    def process_evidence(self, evidence_list: List[Evidence],
                        context: PartitionContext) -> DomainPartitionResult:
        """
        Process evidence to identify domains.

        Args:
            evidence_list: List of validated Evidence objects
            context: Partition context with protein information

        Returns:
            DomainPartitionResult with identified domains
        """
        self.logger.info(
            f"Processing {len(evidence_list)} evidence items for {context.protein_id}"
        )

        # Create result
        result = DomainPartitionResult(
            pdb_id=context.pdb_id,
            chain_id=context.chain_id,
            reference=context.reference,
            sequence_length=context.sequence_length,
            domain_file=str(context.output_file)
        )

        # Check if peptide
        if context.sequence_length < self.options.min_domain_size:
            self.logger.info(f"{context.protein_id} is a peptide (length: {context.sequence_length})")
            result.is_peptide = True
            result.is_unclassified = True
            result.success = True
            return result

        # No evidence means unclassified
        if not evidence_list:
            self.logger.info(f"No evidence for {context.protein_id} - marking as unclassified")
            result.is_unclassified = True
            result.success = True
            return result

        try:
            # Identify domain boundaries
            domain_candidates = self.identify_domain_boundaries(
                evidence_list, context.sequence_length
            )

            if not domain_candidates:
                self.logger.warning(f"No domain candidates identified for {context.protein_id}")
                result.is_unclassified = True
                result.success = True
                return result

            # Resolve overlaps
            resolved_domains = self.resolve_domain_overlaps(
                domain_candidates, context.sequence_length
            )

            # Convert candidates to domain models
            domain_models = []
            for i, candidate in enumerate(resolved_domains):
                domain_id = f"{context.protein_id}_d{i+1}"
                domain_model = candidate.to_domain_model(domain_id)

                # Validate domain
                validation = self.analyzer.validate_domain(
                    domain_model,
                    context=f"final_domain_{i+1}",
                    sequence_length=context.sequence_length
                )

                if validation.is_valid or self.options.validation_level != ValidationLevel.STRICT:
                    domain_models.append(domain_model)
                    self.stats['domains_validated'] += 1
                else:
                    self.logger.warning(
                        f"Domain {domain_id} failed validation: {validation.get_summary()}"
                    )

            # Assign classifications
            self.assign_domain_classifications(domain_models)

            # Set results
            result.domains = domain_models
            result.is_classified = len(domain_models) > 0
            result.success = True

            # Update statistics
            self.stats['domains_identified'] += len(domain_models)

            self.logger.info(
                f"Identified {len(domain_models)} domains for {context.protein_id}"
            )

            return result

        except Exception as e:
            self.logger.error(f"Error processing evidence: {str(e)}", exc_info=True)
            result.success = False
            result.error = str(e)
            return result

    def identify_domain_boundaries(self, evidence_list: List[Evidence],
                                 sequence_length: int) -> List[DomainCandidate]:
        """
        Identify domain boundaries from evidence.

        This implements the core algorithm for determining domain boundaries
        based on consensus from multiple evidence sources.

        Args:
            evidence_list: List of Evidence objects
            sequence_length: Length of the protein sequence

        Returns:
            List of DomainCandidate objects
        """
        # Group evidence by position windows
        evidence_groups = self.analyzer.group_evidence_by_position(
            evidence_list,
            window_size=self.options.merge_gap_tolerance * 2
        )

        if not evidence_groups:
            return []

        # Create domain candidates from evidence groups
        candidates = []

        for position_key, group in evidence_groups.items():
            # Skip groups with insufficient evidence
            if len(group.evidence_items) < 1:
                continue

            # Skip low-confidence groups
            if group.consensus_confidence < self.options.min_evidence_confidence:
                self.logger.debug(
                    f"Skipping low-confidence group at position {position_key}: "
                    f"{group.consensus_confidence:.3f}"
                )
                continue

            # Get consensus boundaries
            if group.consensus_start and group.consensus_end:
                # Validate boundaries
                if (group.consensus_start > 0 and
                    group.consensus_end > 0 and
                    group.consensus_start <= group.consensus_end):

                    # Check sequence bounds
                    if sequence_length > 0:
                        group.consensus_end = min(group.consensus_end, sequence_length)

                    # Create candidate
                    best_evidence = group.get_best_evidence()
                    candidate = DomainCandidate(
                        start=group.consensus_start,
                        end=group.consensus_end,
                        evidence_group=group,
                        source=best_evidence.type if best_evidence else "unknown",
                        confidence=group.consensus_confidence
                    )

                    # Check domain size
                    if candidate.size >= self.options.min_domain_size:
                        if not self.options.max_domain_size or candidate.size <= self.options.max_domain_size:
                            candidates.append(candidate)
                            self.logger.debug(
                                f"Created domain candidate: {candidate.range} "
                                f"(confidence: {candidate.confidence:.3f})"
                            )
                        else:
                            self.logger.debug(
                                f"Domain candidate too large: {candidate.size} > "
                                f"{self.options.max_domain_size}"
                            )
                    else:
                        self.logger.debug(
                            f"Domain candidate too small: {candidate.size} < "
                            f"{self.options.min_domain_size}"
                        )

        # Merge nearby candidates
        merged_candidates = self._merge_nearby_candidates(candidates)

        # Sort by position
        merged_candidates.sort(key=lambda c: c.start)

        return merged_candidates

    def _merge_nearby_candidates(self, candidates: List[DomainCandidate]) -> List[DomainCandidate]:
        """Merge domain candidates that are very close together"""
        if len(candidates) <= 1:
            return candidates

        # Sort by start position
        sorted_candidates = sorted(candidates, key=lambda c: c.start)
        merged = []
        current = sorted_candidates[0]

        for next_candidate in sorted_candidates[1:]:
            # Check if candidates should be merged
            gap = next_candidate.start - current.end - 1

            if gap <= self.options.merge_gap_tolerance:
                # Merge candidates
                self.logger.debug(
                    f"Merging candidates {current.range} and {next_candidate.range} "
                    f"(gap: {gap})"
                )

                # Combine evidence
                current.evidence_group.evidence_items.extend(
                    next_candidate.evidence_group.evidence_items
                )

                # Update boundaries
                current.end = max(current.end, next_candidate.end)

                # Update confidence (weighted average)
                total_weight = len(current.evidence_group.evidence_items) + len(next_candidate.evidence_group.evidence_items)
                current.confidence = (
                    (current.confidence * len(current.evidence_group.evidence_items) +
                     next_candidate.confidence * len(next_candidate.evidence_group.evidence_items)) /
                    total_weight
                )

                # Update source if next candidate has higher confidence
                if next_candidate.confidence > current.confidence:
                    current.source = next_candidate.source
            else:
                # No merge - add current and move to next
                merged.append(current)
                current = next_candidate

        # Add the last candidate
        merged.append(current)

        self.logger.info(
            f"Merged {len(candidates)} candidates into {len(merged)} domains"
        )

        return merged

    def resolve_domain_overlaps(self, candidates: List[DomainCandidate],
                              sequence_length: int) -> List[DomainCandidate]:
        """
        Resolve overlapping domains based on confidence and protection status.

        Args:
            candidates: List of domain candidates
            sequence_length: Length of the protein sequence

        Returns:
            List of non-overlapping domain candidates
        """
        if len(candidates) <= 1:
            return candidates

        # Mark protected domains (very high confidence)
        for candidate in candidates:
            if candidate.confidence >= 0.98:
                candidate.protected = True
                self.logger.debug(f"Protected domain: {candidate.range} (confidence: {candidate.confidence:.3f})")

        # Sort by confidence (descending), protection status, and size
        sorted_candidates = sorted(
            candidates,
            key=lambda c: (c.protected, c.confidence, c.size),
            reverse=True
        )

        # Track covered positions
        covered_positions = set()
        resolved = []

        for candidate in sorted_candidates:
            # Get positions covered by this candidate
            candidate_positions = set(range(candidate.start, candidate.end + 1))

            # Calculate overlap with already covered positions
            overlap = candidate_positions.intersection(covered_positions)
            overlap_pct = len(overlap) / len(candidate_positions) if candidate_positions else 0

            # Decision logic
            if candidate.protected:
                # Protected domains are always included
                resolved.append(candidate)
                covered_positions.update(candidate_positions)

                if overlap_pct > 0:
                    self.logger.debug(
                        f"Including protected domain {candidate.range} despite "
                        f"{overlap_pct:.1%} overlap"
                    )
                    self.stats['overlaps_resolved'] += 1

            elif overlap_pct < self.options.overlap_threshold:
                # Include if minimal overlap
                resolved.append(candidate)
                covered_positions.update(candidate_positions)

                self.logger.debug(
                    f"Including domain {candidate.range} with {overlap_pct:.1%} overlap"
                )
            else:
                # Skip due to excessive overlap
                self.logger.debug(
                    f"Skipping domain {candidate.range} due to {overlap_pct:.1%} overlap"
                )
                self.stats['overlaps_resolved'] += 1

        # Re-sort by position
        resolved.sort(key=lambda c: c.start)

        self.logger.info(
            f"Resolved {len(candidates)} candidates to {len(resolved)} non-overlapping domains"
        )

        return resolved

    def assign_domain_classifications(self, domains: List[DomainModel]) -> None:
        """
        Assign classifications to domains using their evidence.

        Args:
            domains: List of DomainModel objects to classify
        """
        for domain in domains:
            try:
                # Skip if already fully classified
                if domain.is_fully_classified():
                    self.logger.debug(f"Domain {domain.id} already fully classified")
                    continue

                # Find best evidence for classification
                best_evidence = self._find_best_classification_evidence(domain)

                if best_evidence:
                    # Apply classification from best evidence
                    updated = False

                    for attr in ['t_group', 'h_group', 'x_group', 'a_group']:
                        current_value = getattr(domain, attr)
                        evidence_value = getattr(best_evidence, attr, None)

                        if not current_value and evidence_value:
                            setattr(domain, attr, evidence_value)
                            updated = True

                    if updated:
                        self.stats['classifications_assigned'] += 1
                        self.logger.debug(
                            f"Applied classification to domain {domain.id} from "
                            f"{best_evidence.type} evidence: {domain.get_classification_level()}"
                        )

                    # If still not fully classified and DB available, try lookup
                    if not domain.is_fully_classified() and self.db and best_evidence.domain_id:
                        self._enrich_classification_from_db(domain, best_evidence.domain_id)
                else:
                    self.logger.debug(f"No classification evidence found for domain {domain.id}")

            except Exception as e:
                self.logger.warning(f"Error assigning classification to domain {domain.id}: {e}")

    def _find_best_classification_evidence(self, domain: DomainModel) -> Optional[Evidence]:
        """Find the best evidence for classification"""
        if not domain.evidence:
            return None

        # Score evidence based on type and classification completeness
        best_evidence = None
        best_score = 0

        for evidence in domain.evidence:
            # Skip evidence without classification
            if not any(getattr(evidence, attr, None) for attr in ['t_group', 'h_group', 'x_group', 'a_group']):
                continue

            # Calculate score
            type_weight = self.options.classification_confidence_weight.get(evidence.type, 1.0)
            confidence = evidence.confidence or 0.0

            # Count classification levels
            classification_levels = sum(1 for attr in ['t_group', 'h_group', 'x_group', 'a_group']
                                      if getattr(evidence, attr, None))

            # Score = type_weight * confidence * classification_completeness
            score = type_weight * confidence * (classification_levels / 4.0)

            # Prefer HHSearch if configured
            if self.options.prefer_hhsearch_classification and evidence.type == "hhsearch":
                score *= 1.5

            if score > best_score:
                best_score = score
                best_evidence = evidence

        return best_evidence

    def _enrich_classification_from_db(self, domain: DomainModel, domain_id: str) -> None:
        """Enrich domain classification from database"""
        try:
            # Check cache first
            classification = self.classification_cache.get_domain_classification(domain_id)

            if classification is None and self.db:
                # Query database
                query = """
                SELECT t_group, h_group, x_group, a_group,
                       is_manual_rep, is_f70, is_f40, is_f99
                FROM pdb_analysis.domain
                WHERE domain_id = %s
                """

                rows = self.db.execute_dict_query(query, (domain_id,))
                if rows:
                    classification = rows[0]
                    self.classification_cache.set_domain_classification(domain_id, classification)

            if classification:
                # Apply classification
                updated = False
                for attr in ['t_group', 'h_group', 'x_group', 'a_group']:
                    current_value = getattr(domain, attr)
                    db_value = classification.get(attr)

                    if not current_value and db_value:
                        setattr(domain, attr, db_value)
                        updated = True

                # Apply flags
                for flag in ['is_manual_rep', 'is_f70', 'is_f40', 'is_f99']:
                    if flag in classification:
                        setattr(domain, flag, classification[flag])

                if updated:
                    self.logger.debug(
                        f"Enriched classification for domain {domain.id} from database"
                    )

        except Exception as e:
            self.logger.warning(f"Error enriching classification from database: {e}")

    def get_statistics(self) -> Dict[str, Any]:
        """Get processing statistics"""
        stats = self.stats.copy()
        stats['cache_stats'] = self.classification_cache.get_stats()
        return stats

    def reset_statistics(self) -> None:
        """Reset processing statistics"""
        self.stats = {
            'domains_identified': 0,
            'overlaps_resolved': 0,
            'classifications_assigned': 0,
            'domains_validated': 0
        }

    def clear_cache(self) -> None:
        """Clear classification cache"""
        self.classification_cache.clear()
        self.logger.info("Cleared processor cache")
