#!/usr/bin/env python3
"""
Domain partition processor implementation.

Handles the core domain identification and overlap resolution logic.
"""

import logging
from typing import List, Dict, Any, Optional
from dataclasses import dataclass

from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.partition import DomainPartitionResult
from ecod.pipelines.domain_analysis.partition.models import (
    PartitionOptions, PartitionContext, DomainCandidate, EvidenceGroup,
    ValidationResult, ValidationLevel
)


class PartitionProcessor:
    """Processes evidence into domain boundaries"""

    def __init__(self, options: PartitionOptions, analyzer):
        """Initialize processor

        Args:
            options: Partition configuration options
            analyzer: Evidence analyzer instance
        """
        self.options = options
        self.analyzer = analyzer
        self.logger = logging.getLogger(__name__)

        # Initialize statistics
        self.stats = {
            'domains_identified': 0,
            'overlaps_resolved': 0,
            'classifications_assigned': 0,
            'domains_validated': 0,
            'discontinuous_domains': 0,
            'cache_stats': {}
        }

    def process_evidence(
        self,
        evidence_list: List[Evidence],
        context: PartitionContext
    ) -> DomainPartitionResult:
        """Process evidence into domain partition result

        Args:
            evidence_list: List of evidence to process
            context: Processing context

        Returns:
            DomainPartitionResult with identified domains
        """
        try:
            # Check for peptide
            if context.sequence_length > 0 and context.sequence_length < self.options.min_domain_size:
                return DomainPartitionResult(
                    pdb_id=context.pdb_id,
                    chain_id=context.chain_id,
                    reference=context.reference,
                    success=True,
                    is_peptide=True,
                    is_unclassified=True,
                    sequence_length=context.sequence_length
                )

            # If no evidence, mark as unclassified
            if not evidence_list:
                return DomainPartitionResult(
                    pdb_id=context.pdb_id,
                    chain_id=context.chain_id,
                    reference=context.reference,
                    success=True,
                    is_unclassified=True,
                    sequence_length=context.sequence_length
                )

            # Group evidence by position
            evidence_groups = self.analyzer.group_evidence_by_position(evidence_list)

            # Create domain candidates from evidence groups
            candidates = self._create_domain_candidates(evidence_groups, context)

            # Resolve overlaps
            resolved_candidates = self.resolve_domain_overlaps(candidates, context.sequence_length)

            # Convert to domain models
            domains = self._candidates_to_domains(resolved_candidates, context)

            # Assign classifications
            self.assign_domain_classifications(domains)

            # Create result
            result = DomainPartitionResult(
                pdb_id=context.pdb_id,
                chain_id=context.chain_id,
                reference=context.reference,
                domains=domains,
                success=True,
                is_classified=len(domains) > 0,
                is_unclassified=len(domains) == 0,
                sequence_length=context.sequence_length
            )

            self.stats['domains_identified'] += len(domains)
            return result

        except Exception as e:
            self.logger.error(f"Error processing evidence for {context.chain_identifier}: {e}")
            return DomainPartitionResult(
                pdb_id=context.pdb_id,
                chain_id=context.chain_id,
                reference=context.reference,
                success=False,
                error=str(e),
                sequence_length=context.sequence_length
            )

    def _create_domain_candidates(
        self,
        evidence_groups: Dict[int, EvidenceGroup],
        context: PartitionContext
    ) -> List[DomainCandidate]:
        """Create domain candidates from evidence groups"""
        candidates = []

        for position_key, group in evidence_groups.items():
            if group.consensus_start is not None and group.consensus_end is not None:
                # Check minimum size
                size = group.consensus_end - group.consensus_start + 1
                if size >= self.options.min_domain_size:
                    candidate = DomainCandidate(
                        start=group.consensus_start,
                        end=group.consensus_end,
                        evidence_group=group,
                        confidence=group.consensus_confidence,
                        source=group.get_dominant_source()
                    )
