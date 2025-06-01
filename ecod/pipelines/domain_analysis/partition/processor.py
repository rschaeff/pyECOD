#!/usr/bin/env python3
"""
Core domain partitioning processor with reference coverage validation.

This module contains the core algorithms for:
- Identifying domain boundaries from evidence
- Validating reference domain coverage
- Handling discontinuous domains properly
- Resolving overlapping domains
- Assigning classifications to domains
"""

import logging
from typing import Dict, Any, List, Optional, Tuple, Set
from collections import defaultdict
from dataclasses import dataclass, field

from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult
from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.core.sequence_range import SequenceRange

from .models import (
    PartitionOptions, EvidenceGroup, DomainCandidate,
    PartitionContext, ValidationLevel, ClassificationCache,
    PartitionStage
)
from .analyzer import EvidenceAnalyzer
from .reference_analyzer import ReferenceCoverageAnalyzer, EvidenceWithCoverage


@dataclass
class DiscontinuousDomainCandidate(DomainCandidate):
    """Domain candidate with multiple segments"""
    segments: List[Tuple[int, int]] = field(default_factory=list)

    def __init__(self, segments: List[Tuple[int, int]], evidence_group: EvidenceGroup,
                 source: str = "", confidence: float = 0.0, protected: bool = False):
        # Calculate overall start and end from segments
        if segments:
            start = min(s[0] for s in segments)
            end = max(s[1] for s in segments)
        else:
            start = end = 0

        # Set segments first
        self.segments = segments

        # Initialize parent with calculated start/end - but bypass property conflict
        # by setting the underlying attributes directly
        self.__dict__['start'] = start
        self.__dict__['end'] = end
        self.evidence_group = evidence_group
        self.source = source
        self.confidence = confidence
        self.protected = protected

    @property
    def start(self) -> int:
        """Get start of first segment"""
        return min(s[0] for s in self.segments) if self.segments else self.__dict__.get('start', 0)

    @property
    def end(self) -> int:
        """Get end of last segment"""
        return max(s[1] for s in self.segments) if self.segments else self.__dict__.get('end', 0)

    @property
    def size(self) -> int:
        """Get total size across all segments"""
        return sum(s[1] - s[0] + 1 for s in self.segments)

    @property
    def range(self) -> str:
        """Get discontinuous range string"""
        return ','.join(f"{s[0]}-{s[1]}" for s in self.segments)

    def to_domain_model(self, domain_id: str) -> DomainModel:
        """Convert to DomainModel with discontinuous range"""
        best_evidence = self.evidence_group.get_best_evidence()

        return DomainModel(
            id=domain_id,
            start=self.start,
            end=self.end,
            range=self.range,
            source=self.source or (best_evidence.type if best_evidence else "unknown"),
            confidence=self.confidence,
            source_id=best_evidence.source_id if best_evidence else "",
            evidence=self.evidence_group.evidence_items,
            protected=self.protected,
            is_discontinuous=True
        )


class PartitionProcessor:
    """
    Core processor for domain partitioning with reference coverage validation.

    This class implements the main algorithms for:
    - Identifying domain boundaries from evidence with coverage validation
    - Properly handling discontinuous domains
    - Resolving overlapping domains
    - Assigning classifications to domains
    """

    def __init__(self, options: PartitionOptions,
                 analyzer: Optional[EvidenceAnalyzer] = None,
                 context: Optional[ApplicationContext] = None):
        """
        Initialize the partition processor.

        Args:
            options: Partitioning options
            analyzer: Evidence analyzer (created if not provided)
            db_manager: Database manager for classification lookups
        """
        self.options = options
        self.analyzer = analyzer or EvidenceAnalyzer(options)
        self.logger = logging.getLogger(__name__)

        # Get database from context
        self.db = None
        if context and hasattr(context, 'db_manager'):
            self.db = context.db_manager
        elif context and hasattr(context, 'db'):
            self.db = context.db

        # Initialize reference coverage analyzer
        self.coverage_analyzer = ReferenceCoverageAnalyzer(
            self.db,
            min_coverage=options.min_reference_coverage
        ) if self.db else None

        # Classification cache
        self.classification_cache = ClassificationCache()

        # Processing statistics
        self.stats = {
            'domains_identified': 0,
            'overlaps_resolved': 0,
            'classifications_assigned': 0,
            'domains_validated': 0,
            'discontinuous_domains': 0,
            'low_coverage_rejected': 0,
            'domains_extended': 0
        }

    def process_evidence(self, evidence_list: List[Evidence],
                        context: PartitionContext) -> DomainPartitionResult:
        """
        Process evidence to identify domains with reference coverage validation.

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
            # Enhance evidence with coverage information
            context.record_stage_time(PartitionStage.VALIDATING_EVIDENCE)
            enhanced_evidence = self._enhance_evidence_with_coverage(evidence_list)

            if not enhanced_evidence:
                self.logger.warning(f"All evidence rejected due to low coverage for {context.protein_id}")
                result.is_unclassified = True
                result.success = True
                return result

            # Process evidence with coverage awareness
            context.record_stage_time(PartitionStage.IDENTIFYING_BOUNDARIES)
            domain_candidates = self._process_evidence_with_coverage(
                enhanced_evidence, context.sequence_length
            )

            if not domain_candidates:
                self.logger.warning(f"No domain candidates identified for {context.protein_id}")
                result.is_unclassified = True
                result.success = True
                return result

            # Resolve overlaps
            context.record_stage_time(PartitionStage.RESOLVING_OVERLAPS)
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

                    if isinstance(candidate, DiscontinuousDomainCandidate):
                        self.stats['discontinuous_domains'] += 1
                else:
                    self.logger.warning(
                        f"Domain {domain_id} failed validation: {validation.get_summary()}"
                    )

            # Assign classifications
            context.record_stage_time(PartitionStage.ASSIGNING_CLASSIFICATIONS)
            self.assign_domain_classifications(domain_models)

            # Set results
            result.domains = domain_models
            result.calculate_coverage()
            result.is_classified = len(domain_models) > 0
            result.success = True

            # Update statistics
            self.stats['domains_identified'] += len(domain_models)

            self.logger.info(
                f"Identified {len(domain_models)} domains for {context.protein_id} "
                f"({self.stats['discontinuous_domains']} discontinuous)"
            )

            return result

        except Exception as e:
            self.logger.error(f"Error processing evidence: {str(e)}", exc_info=True)
            result.success = False
            result.error = str(e)
            return result

    def _enhance_evidence_with_coverage(self, evidence_list: List[Evidence]) -> List[Evidence]:
        """Enhance evidence with reference coverage information"""
        if not self.coverage_analyzer:
            self.logger.warning("No coverage analyzer available - skipping coverage enhancement")
            return evidence_list

        # Evidence-type-specific coverage thresholds
        coverage_thresholds = {
            'hhsearch': 0.65,
            'chain_blast': 0.80,
            'domain_blast': 0.50,
        }

        enhanced_evidence = []
        coverage_stats = {'high': 0, 'medium': 0, 'low': 0, 'rejected': 0,
                          'rejected_by_type': {'hhsearch': 0, 'chain_blast': 0, 'domain_blast': 0}}

        for evidence in evidence_list:
            # DEBUG: Check what we're getting
            if evidence.type == 'hhsearch':
                self.logger.debug(f"HHSearch evidence: {evidence.source_id}, has coverage analyzer: {self.coverage_analyzer is not None}")

            enhanced = self.coverage_analyzer.analyze_evidence_coverage(evidence)

            # DEBUG: Log coverage info
            if evidence.type == 'hhsearch':
                self.logger.debug(f"HHSearch {evidence.source_id}: reference_coverage = {getattr(enhanced, 'reference_coverage', 'NONE')}")

            # Get evidence-type-specific threshold
            evidence_threshold = coverage_thresholds.get(
                evidence.type,
                self.options.min_reference_coverage
            )

            # Apply evidence-specific coverage threshold
            if hasattr(enhanced, 'reference_coverage') and enhanced.reference_coverage >= evidence_threshold:
                if enhanced.reference_coverage >= self.options.strict_reference_coverage:
                    coverage_stats['high'] += 1
                else:
                    coverage_stats['medium'] += 1
                enhanced_evidence.append(enhanced)
            elif hasattr(enhanced, 'reference_coverage') and enhanced.reference_coverage >= self.options.partial_coverage_threshold:
                coverage_stats['low'] += 1
                enhanced_evidence.append(enhanced)
            else:
                # Reject
                coverage_stats['rejected'] += 1
                if evidence.type in coverage_stats['rejected_by_type']:
                    coverage_stats['rejected_by_type'][evidence.type] += 1

                self.stats['low_coverage_rejected'] += 1

                if hasattr(enhanced, 'reference_coverage'):
                    self.logger.info(
                        f"Rejected {evidence.type} evidence with {enhanced.reference_coverage:.1%} "
                        f"coverage for {enhanced.domain_id} (threshold: {evidence_threshold:.1%})"
                    )
                else:
                    self.logger.warning(f"No reference coverage for {evidence.type} evidence {evidence.source_id}")

        # Log summary
        self.logger.info(
            f"Evidence coverage filtering complete: {len(enhanced_evidence)}/{len(evidence_list)} passed. "
            f"Rejected by type: {coverage_stats['rejected_by_type']}"
        )

        return enhanced_evidence

    def _process_evidence_with_coverage(self, evidence_list: List[Evidence],
                                      sequence_length: int) -> List[DomainCandidate]:
        """Process evidence with proper handling of discontinuous domains and coverage"""

        # Separate evidence by type
        discontinuous_evidence = []
        continuous_evidence = []

        for evidence in evidence_list:
            if hasattr(evidence, 'is_discontinuous') and evidence.is_discontinuous:
                discontinuous_evidence.append(evidence)
            else:
                continuous_evidence.append(evidence)

        domain_candidates = []

        # Process discontinuous evidence
        if discontinuous_evidence:
            discontinuous_domains = self._process_discontinuous_evidence(
                discontinuous_evidence, sequence_length
            )
            domain_candidates.extend(discontinuous_domains)
            self.logger.info(f"Created {len(discontinuous_domains)} discontinuous domain candidates")

        # Process continuous evidence
        if continuous_evidence:
            continuous_domains = self._process_continuous_evidence(
                continuous_evidence, sequence_length
            )
            domain_candidates.extend(continuous_domains)

        return domain_candidates

    def _process_discontinuous_evidence(self, evidence_list: List[Evidence],
                                      sequence_length: int) -> List[DomainCandidate]:
        """Process evidence with multi-segment domains"""
        candidates = []

        for evidence in evidence_list:
            # Parse multi-segment range
            if hasattr(evidence, 'query_range') and evidence.query_range:
                segments = self._parse_discontinuous_range(evidence.query_range)

                if segments and len(segments) > 1:
                    # Create a discontinuous candidate
                    # Apply high protection if high coverage
                    protected = False
                    if isinstance(evidence, EvidenceWithCoverage):
                        protected = evidence.reference_coverage >= self.options.strict_reference_coverage
                    else:
                        protected = evidence.confidence >= 0.95

                    candidate = DiscontinuousDomainCandidate(
                        segments=segments,
                        evidence_group=EvidenceGroup([evidence]),
                        source=f"{evidence.type}_discontinuous",
                        confidence=evidence.confidence or 0.0,
                        protected=protected
                    )
                    candidates.append(candidate)

                    self.logger.debug(
                        f"Created discontinuous domain candidate: {candidate.range} "
                        f"from {evidence.source_id}"
                    )
                elif segments and len(segments) == 1:
                    # Single segment, treat as continuous
                    start, end = segments[0]
                    candidate = DomainCandidate(
                        start=start,
                        end=end,
                        evidence_group=EvidenceGroup([evidence]),
                        source=evidence.type,
                        confidence=evidence.confidence or 0.0
                    )
                    candidates.append(candidate)

        return candidates

    def _process_continuous_evidence(self, evidence_list: List[Evidence],
                                   sequence_length: int) -> List[DomainCandidate]:
        """Process continuous evidence with coverage-based boundary decisions"""

        # Group evidence by position windows
        evidence_groups = self.analyzer.group_evidence_by_position(
            evidence_list,
            window_size=self.options.merge_gap_tolerance * 2
        )

        if not evidence_groups:
            return []

        candidates = []

        for position_key, evidence_items in evidence_groups.items():
            # Skip groups with insufficient evidence
            if len(evidence_items) < 1:
                continue

            # Create EvidenceGroup from list
            from ecod.pipelines.domain_analysis.partition.models import EvidenceGroup
            group = EvidenceGroup(evidence_items=evidence_items)

            # Get best evidence from group
            best_evidence = group.get_best_evidence()
            if not best_evidence:
                continue

            # For high coverage evidence, trust boundaries exactly
            if isinstance(best_evidence, EvidenceWithCoverage):
                if best_evidence.reference_coverage >= self.options.strict_reference_coverage:
                    # Trust boundaries exactly
                    if best_evidence.query_range:
                        start, end = self._parse_range(best_evidence.query_range)
                        candidate = DomainCandidate(
                            start=start,
                            end=end,
                            evidence_group=group,
                            source=f"{best_evidence.type}_high_coverage",
                            confidence=best_evidence.confidence or 0.0,
                            protected=True
                        )
                        candidates.append(candidate)
                        continue
                elif best_evidence.reference_coverage >= self.options.min_reference_coverage:
                    # Medium coverage - consider extension
                    if best_evidence.query_range:
                        start, end = self._parse_range(best_evidence.query_range)

                        # Optionally extend based on reference size
                        if (self.options.extend_to_reference_size and
                            self.coverage_analyzer and
                            best_evidence.reference_info):

                            extended_start, extended_end = self.coverage_analyzer.suggest_extended_boundaries(
                                best_evidence, start, end, sequence_length
                            )

                            if extended_start != start or extended_end != end:
                                self.stats['domains_extended'] += 1
                                self.logger.debug(
                                    f"Extended domain from {start}-{end} to "
                                    f"{extended_start}-{extended_end} based on reference"
                                )

                            start, end = extended_start, extended_end

                        candidate = DomainCandidate(
                            start=start,
                            end=end,
                            evidence_group=group,
                            source=best_evidence.type,
                            confidence=best_evidence.confidence or 0.0
                        )
                        candidates.append(candidate)
                        continue

            # Default processing for evidence without coverage info
            # FIXED: Calculate consensus from evidence directly since group may not have consensus attributes
            consensus_starts = []
            consensus_ends = []

            for evidence in group.evidence_items:
                if evidence.query_range:
                    try:
                        start, end = self._parse_range(evidence.query_range)
                        if start > 0 and end > start:
                            consensus_starts.append(start)
                            consensus_ends.append(end)
                    except:
                        continue

            if consensus_starts and consensus_ends:
                consensus_start = min(consensus_starts)  # Use most conservative start
                consensus_end = max(consensus_ends)      # Use most liberal end

                candidate = DomainCandidate(
                    start=consensus_start,
                    end=consensus_end,
                    evidence_group=group,
                    source=best_evidence.type if best_evidence else "unknown",
                    confidence=best_evidence.confidence or 0.0
                )

                # Check domain size
                if candidate.size >= self.options.min_domain_size:
                    if not self.options.max_domain_size or candidate.size <= self.options.max_domain_size:
                        candidates.append(candidate)

        # Merge nearby candidates
        merged_candidates = self._merge_nearby_candidates(candidates)

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
            # Don't merge protected domains
            if current.protected or next_candidate.protected:
                merged.append(current)
                current = next_candidate
                continue

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

                # Update confidence (weighted average) - FIX: Handle empty evidence lists
                current_evidence_count = len(current.evidence_group.evidence_items)
                next_evidence_count = len(next_candidate.evidence_group.evidence_items)
                total_evidence_count = current_evidence_count + next_evidence_count

                if total_evidence_count > 0:
                    current.confidence = (
                        (current.confidence * current_evidence_count +
                         next_candidate.confidence * next_evidence_count) /
                        total_evidence_count
                    )
                else:
                    # If no evidence, use simple average
                    current.confidence = (current.confidence + next_candidate.confidence) / 2

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
        Protected domains are always included and can coexist with overlapping domains.
        """
        if len(candidates) <= 1:
            return candidates

        # Sort by protection status first (protected first), then by confidence (descending)
        sorted_candidates = sorted(
            candidates,
            key=lambda c: (not c.protected, -c.confidence, -c.size)
        )

        # Track positions covered by NON-PROTECTED domains only
        # Protected domains can overlap with anything
        non_protected_positions = set()
        resolved = []

        for candidate in sorted_candidates:
            # Get positions covered by this candidate
            if isinstance(candidate, DiscontinuousDomainCandidate):
                candidate_positions = set()
                for start, end in candidate.segments:
                    candidate_positions.update(range(start, end + 1))
            else:
                candidate_positions = set(range(candidate.start, candidate.end + 1))

            if candidate.protected:
                # Protected domains are ALWAYS included, regardless of overlap
                resolved.append(candidate)
                self.logger.debug(f"Including protected domain {candidate.range}")

            else:
                # For non-protected domains, check overlap with other non-protected domains only
                overlap = candidate_positions.intersection(non_protected_positions)
                overlap_pct = len(overlap) / len(candidate_positions) if candidate_positions else 0

                if overlap_pct < self.options.overlap_threshold:
                    # Include if minimal overlap with existing non-protected domains
                    resolved.append(candidate)
                    non_protected_positions.update(candidate_positions)

                    self.logger.debug(
                        f"Including non-protected domain {candidate.range} with {overlap_pct:.1%} overlap"
                    )
                else:
                    # Skip due to excessive overlap with existing non-protected domains
                    self.logger.debug(
                        f"Skipping non-protected domain {candidate.range} due to {overlap_pct:.1%} overlap"
                    )
                    self.stats['overlaps_resolved'] += 1

        # Re-sort by position for final output
        resolved.sort(key=lambda c: c.start)

        self.logger.info(
            f"Resolved {len(candidates)} candidates to {len(resolved)} domains "
            f"({sum(1 for d in resolved if d.protected)} protected)"
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

            # Bonus for high coverage evidence
            if isinstance(evidence, EvidenceWithCoverage) and evidence.reference_coverage >= 0.9:
                score *= 1.2

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

    def _parse_discontinuous_range(self, range_str: str) -> List[Tuple[int, int]]:
        """Parse discontinuous range - NOW USING UNIFIED PARSER"""
        try:
            sequence_range = SequenceRange.parse(range_str)
            return [(seg.start, seg.end) for seg in sequence_range.segments]
        except ValueError as e:
            self.logger.warning(f"Invalid range string: {range_str} - {e}")
            return []

    def _parse_range(self, range_str: str) -> Tuple[int, int]:
        """Parse range string - handle discontinuous ranges gracefully"""
        if not range_str:
            return (0, 0)

        # Handle discontinuous ranges by taking overall span
        if ',' in range_str:
            try:
                seq_range = SequenceRange.parse(range_str)
                return seq_range.span  # Returns (overall_start, overall_end)
            except ValueError as e:
                self.logger.warning(f"Invalid range string: {range_str} - {e}")
                return (0, 0)

        # Handle simple ranges
        if '-' not in range_str:
            return (0, 0)

        try:
            parts = range_str.split('-', 1)  # Only split on first dash
            return (int(parts[0]), int(parts[1]))
        except (ValueError, IndexError):
            self.logger.warning(f"Invalid range string: {range_str}")
            return (0, 0)

    def _domains_overlap(self, domain1: Tuple[int, int], domain2: Tuple[int, int]) -> bool:
        """Check if two domains overlap"""
        start1, end1 = domain1
        start2, end2 = domain2

        return not (end1 < start2 or end2 < start1)

    def get_statistics(self) -> Dict[str, Any]:
        """Get processing statistics"""
        stats = self.stats.copy()
        stats['cache_stats'] = self.classification_cache.get_stats()
        if self.coverage_analyzer:
            stats['coverage_cache_stats'] = self.coverage_analyzer.get_cache_statistics()
        return stats

    def reset_statistics(self) -> None:
        """Reset processing statistics"""
        self.stats = {
            'domains_identified': 0,
            'overlaps_resolved': 0,
            'classifications_assigned': 0,
            'domains_validated': 0,
            'discontinuous_domains': 0,
            'low_coverage_rejected': 0,
            'domains_extended': 0
        }

    def clear_cache(self) -> None:
        """Clear all caches"""
        self.classification_cache.clear()
        if self.coverage_analyzer:
            self.coverage_analyzer._reference_cache.clear()
        self.logger.info("Cleared processor caches")
