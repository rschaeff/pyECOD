#!/usr/bin/env python3
"""
Evidence analysis and validation for domain partitioning.

Provides comprehensive evidence parsing, validation, and grouping
functionality for the domain partitioning pipeline.
"""

import logging
import xml.etree.ElementTree as ET
from typing import List, Dict, Any, Optional, Union, Callable
from pathlib import Path
import math

# Import required models
from ecod.models.pipeline.evidence import Evidence
from .models import (
    PartitionOptions, ValidationResult, EvidenceGroup, ValidationLevel,
    ClassificationCache
)


class EvidenceAnalyzer:
    """Analyzes and validates evidence for domain partitioning"""

    def __init__(self, options: PartitionOptions):
        """Initialize analyzer with options

        Args:
            options: Partition configuration options
        """
        self.options = options
        self.logger = logging.getLogger(__name__)

        # Initialize cache if enabled
        if options.use_cache:
            self.classification_cache = ClassificationCache()
        else:
            self.classification_cache = None

        # Initialize statistics
        self.stats = {
            'summaries_parsed': 0,
            'evidence_extracted': 0,
            'validation_errors': 0,
            'classification_lookups': 0,
            'domains_validated': 0
        }

    def parse_domain_summary(self, summary_file_path: str) -> Dict[str, Any]:
        """Parse domain summary XML file

        Args:
            summary_file_path: Path to domain summary XML file

        Returns:
            Dictionary containing parsed summary data or error information
        """
        try:
            if not Path(summary_file_path).exists():
                return {'error': f'Domain summary file not found: {summary_file_path}'}

            # Parse XML
            tree = ET.parse(summary_file_path)
            root = tree.getroot()

            # Extract basic metadata
            result = {
                'pdb_id': root.get('pdb_id', ''),
                'chain_id': root.get('chain_id', ''),
                'sequence_length': int(root.get('sequence_length', 0)),
                'is_peptide': False,
                'is_discontinuous': False,
                'blast_hits': [],
                'hhsearch_hits': []
            }

            # Check metadata section
            metadata = root.find('metadata')
            if metadata is not None:
                is_peptide_elem = metadata.find('is_peptide')
                if is_peptide_elem is not None:
                    result['is_peptide'] = is_peptide_elem.text.lower() == 'true'

                is_discontinuous_elem = metadata.find('is_discontinuous')
                if is_discontinuous_elem is not None:
                    result['is_discontinuous'] = is_discontinuous_elem.text.lower() == 'true'

            # Parse evidence
            evidence_list = root.find('evidence_list')
            if evidence_list is not None:
                # Parse BLAST hits
                blast_hits = evidence_list.find('blast_hits')
                if blast_hits is not None:
                    for hit in blast_hits.findall('hit'):
                        hit_data = dict(hit.attrib)

                        # Extract range elements
                        query_reg = hit.find('query_range')
                        if query_reg is not None:
                            hit_data['query_range'] = query_reg.text

                        hit_reg = hit.find('hit_range')
                        if hit_reg is not None:
                            hit_data['hit_range'] = hit_reg.text

                        result['blast_hits'].append(hit_data)

                # Parse HHSearch hits
                hhsearch_hits = evidence_list.find('hhsearch_hits')
                if hhsearch_hits is not None:
                    for hit in hhsearch_hits.findall('hit'):
                        hit_data = dict(hit.attrib)

                        # Extract range elements
                        query_reg = hit.find('query_range')
                        if query_reg is not None:
                            hit_data['query_range'] = query_reg.text

                        hit_reg = hit.find('hit_range')
                        if hit_reg is not None:
                            hit_data['hit_range'] = hit_reg.text

                        result['hhsearch_hits'].append(hit_data)

            self.stats['summaries_parsed'] += 1
            return result

        except ET.ParseError as e:
            return {'error': f'XML parsing error: {str(e)}'}
        except Exception as e:
            return {'error': f'Error parsing domain summary: {str(e)}'}

    def extract_evidence_with_classification(
        self,
        summary_data: Dict[str, Any],
        db_lookup_func: Optional[Callable[[str], Optional[Dict[str, Any]]]] = None
    ) -> List[Evidence]:
        """Extract evidence from summary data with classification lookup

        Args:
            summary_data: Parsed domain summary data
            db_lookup_func: Optional function to lookup domain classifications

        Returns:
            List of Evidence objects
        """
        evidence_list = []

        try:
            # Process BLAST hits if enabled
            if self.options.use_domain_blast and 'blast_hits' in summary_data:
                for hit_data in summary_data['blast_hits']:
                    evidence = self._create_evidence_from_blast_hit(hit_data)
                    if evidence and self._passes_confidence_filter(evidence):
                        # Add classification if available
                        if db_lookup_func and evidence.domain_id:
                            classification = db_lookup_func(evidence.domain_id)
                            if classification:
                                evidence.t_group = classification.get('t_group')
                                evidence.h_group = classification.get('h_group')
                                evidence.x_group = classification.get('x_group')
                                evidence.a_group = classification.get('a_group')
                                self.stats['classification_lookups'] += 1

                        evidence_list.append(evidence)

            # Process HHSearch hits if enabled
            if self.options.use_hhsearch and 'hhsearch_hits' in summary_data:
                for hit_data in summary_data['hhsearch_hits']:
                    evidence = self._create_evidence_from_hhsearch_hit(hit_data)
                    if evidence and self._passes_confidence_filter(evidence):
                        # Add classification if available
                        if db_lookup_func and evidence.domain_id:
                            classification = db_lookup_func(evidence.domain_id)
                            if classification:
                                evidence.t_group = classification.get('t_group')
                                evidence.h_group = classification.get('h_group')
                                evidence.x_group = classification.get('x_group')
                                evidence.a_group = classification.get('a_group')
                                self.stats['classification_lookups'] += 1

                        evidence_list.append(evidence)

            self.stats['evidence_extracted'] += len(evidence_list)
            return evidence_list

        except Exception as e:
            self.logger.error(f"Error extracting evidence: {e}")
            return []

    def _create_evidence_from_blast_hit(self, hit_data: Dict[str, Any]) -> Optional[Evidence]:
        """Create Evidence object from BLAST hit data"""
        try:
            # Extract evalue safely
            evalue_str = hit_data.get('evalues', '999.0')
            try:
                evalue = float(evalue_str)
                if not math.isfinite(evalue):
                    evalue = 999.0
            except (ValueError, TypeError):
                evalue = 999.0

            evidence = Evidence(
                type='domain_blast',
                source_id=hit_data.get('num', ''),
                domain_id=hit_data.get('domain_id', ''),
                query_range=hit_data.get('query_range', ''),
                hit_range=hit_data.get('hit_range', ''),
                evalue=evalue,
                identity=self._safe_float(hit_data.get('identity')),
                coverage=self._safe_float(hit_data.get('coverage')),
                hsp_count=self._safe_int(hit_data.get('hsp_count', 1))
            )

            return evidence

        except Exception as e:
            self.logger.warning(f"Error creating BLAST evidence: {e}")
            return None

    def _create_evidence_from_hhsearch_hit(self, hit_data: Dict[str, Any]) -> Optional[Evidence]:
        """Create Evidence object from HHSearch hit data"""
        try:
            # Extract probability and evalue safely
            probability = self._safe_float(hit_data.get('probability', 0.0))
            evalue = self._safe_float(hit_data.get('evalue', 999.0))
            score = self._safe_float(hit_data.get('score'))

            evidence = Evidence(
                type='hhsearch',
                source_id=hit_data.get('hit_id', ''),
                domain_id=hit_data.get('domain_id', ''),
                query_range=hit_data.get('query_range', ''),
                hit_range=hit_data.get('hit_range', ''),
                probability=probability,
                evalue=evalue,
                score=score
            )

            return evidence

        except Exception as e:
            self.logger.warning(f"Error creating HHSearch evidence: {e}")
            return None

    def _safe_float(self, value: Any) -> Optional[float]:
        """Safely convert value to float"""
        if value is None:
            return None
        try:
            result = float(value)
            return result if math.isfinite(result) else None
        except (ValueError, TypeError):
            return None

    def _safe_int(self, value: Any) -> Optional[int]:
        """Safely convert value to int"""
        if value is None:
            return None
        try:
            return int(value)
        except (ValueError, TypeError):
            return None

    def _passes_confidence_filter(self, evidence: Evidence) -> bool:
        """Check if evidence passes minimum confidence filter"""
        return evidence.confidence >= self.options.min_evidence_confidence

    def validate_evidence(self, evidence: Evidence, context: str) -> ValidationResult:
        """Validate evidence object

        Args:
            evidence: Evidence to validate
            context: Context string for error reporting

        Returns:
            ValidationResult with validation status and messages
        """
        result = ValidationResult(is_valid=True, context=context)

        try:
            # Check required fields
            if not evidence.type:
                result.add_error("Evidence type is required")

            if not evidence.domain_id:
                result.add_error("Evidence domain_id is required")

            # Check confidence
            if evidence.confidence < self.options.min_evidence_confidence:
                if self.options.validation_level == ValidationLevel.STRICT:
                    result.add_error(f"Evidence confidence {evidence.confidence} below minimum {self.options.min_evidence_confidence}")
                else:
                    result.add_warning(f"Low evidence confidence: {evidence.confidence}")

            # Validate coordinates if required
            if self.options.validate_coordinates and evidence.query_range:
                if not self._validate_range_format(evidence.query_range):
                    result.add_warning(f"Invalid query range format: {evidence.query_range}")

            if result.errors:
                self.stats['validation_errors'] += 1

        except Exception as e:
            result.add_error(f"Exception during validation: {str(e)}")
            self.stats['validation_errors'] += 1

        return result

    def _validate_range_format(self, range_str: str) -> bool:
        """Validate range string format (e.g., '10-50' or '10-20,30-40')"""
        if not range_str:
            return False

        try:
            # Split on commas for discontinuous ranges
            segments = range_str.split(',')
            for segment in segments:
                segment = segment.strip()
                if '-' not in segment:
                    return False

                parts = segment.split('-')
                if len(parts) != 2:
                    return False

                start, end = parts
                start_num = int(start.strip())
                end_num = int(end.strip())

                if start_num <= 0 or end_num <= 0 or start_num > end_num:
                    return False

            return True

        except (ValueError, AttributeError):
            return False

    def group_evidence_by_position(
        self,
        evidence_list: List[Evidence],
        window_size: int = 50
    ) -> Dict[int, EvidenceGroup]:
        """Group evidence by position on sequence

        Args:
            evidence_list: List of evidence to group
            window_size: Size of grouping window

        Returns:
            Dictionary mapping position keys to EvidenceGroup objects
        """
        groups = {}

        for evidence in evidence_list:
            # Parse position from evidence
            start_pos, end_pos = self._parse_evidence_position(evidence)

            if start_pos is None or end_pos is None:
                continue

            # Find or create group
            group_key = self._find_group_key(start_pos, end_pos, groups, window_size)

            if group_key not in groups:
                groups[group_key] = EvidenceGroup()

            groups[group_key].add_evidence(evidence, start_pos, end_pos)

        # Calculate consensus for each group
        for group in groups.values():
            group.consensus_confidence = group.calculate_consensus_confidence()

        return groups

    def _parse_evidence_position(self, evidence: Evidence) -> tuple[Optional[int], Optional[int]]:
        """Parse start and end positions from evidence"""
        if not evidence.query_range:
            return None, None

        try:
            # Handle simple range like "10-50"
            if '-' in evidence.query_range and ',' not in evidence.query_range:
                parts = evidence.query_range.split('-')
                if len(parts) == 2:
                    start = int(parts[0].strip())
                    end = int(parts[1].strip())
                    return start, end

            # Handle complex ranges like "10-20,30-40"
            segments = evidence.query_range.split(',')
            if segments:
                first_segment = segments[0].strip()
                if '-' in first_segment:
                    parts = first_segment.split('-')
                    if len(parts) == 2:
                        start = int(parts[0].strip())

                        # For end, use the last segment
                        last_segment = segments[-1].strip()
                        if '-' in last_segment:
                            last_parts = last_segment.split('-')
                            end = int(last_parts[1].strip())
                            return start, end

        except (ValueError, IndexError):
            pass

        return None, None

    def _find_group_key(
        self,
        start: int,
        end: int,
        existing_groups: Dict[int, EvidenceGroup],
        window_size: int
    ) -> int:
        """Find appropriate group key for position range"""
        center = (start + end) // 2

        # Check if this position fits into an existing group
        for key, group in existing_groups.items():
            if group.consensus_start is not None and group.consensus_end is not None:
                group_center = (group.consensus_start + group.consensus_end) // 2
                if abs(center - group_center) <= window_size:
                    return key

        # Create new group key
        return center

    def validate_domain(
        self,
        domain,
        context: str,
        sequence_length: int = 0
    ) -> ValidationResult:
        """Validate domain model

        Args:
            domain: Domain object to validate
            context: Context string for error reporting
            sequence_length: Sequence length for boundary checking

        Returns:
            ValidationResult with validation status and messages
        """
        result = ValidationResult(is_valid=True, context=context)

        try:
            # Extract domain properties (handle both dict and object)
            if hasattr(domain, 'start'):
                start = domain.start
                end = domain.end
                confidence = getattr(domain, 'confidence', 0.0)
            else:
                start = domain.get('start', 0)
                end = domain.get('end', 0)
                confidence = domain.get('confidence', 0.0)

            # Check domain size
            domain_size = end - start + 1
            if domain_size < self.options.min_domain_size:
                result.add_error(f"Domain size {domain_size} below minimum {self.options.min_domain_size}")

            if hasattr(self.options, 'max_domain_size') and self.options.max_domain_size and domain_size > self.options.max_domain_size:
                result.add_error(f"Domain size {domain_size} above maximum {self.options.max_domain_size}")

            # Check coordinates
            if start <= 0 or end <= 0:
                result.add_error("Domain coordinates must be positive")

            if start > end:
                result.add_error("Domain start must be <= end")

            # Check sequence boundaries
            if sequence_length > 0:
                if start > sequence_length or end > sequence_length:
                    result.add_error(f"Domain extends beyond sequence (length: {sequence_length})")

            # Check confidence
            if confidence < self.options.min_evidence_confidence:
                if self.options.validation_level == ValidationLevel.STRICT:
                    result.add_error(f"Domain confidence {confidence} below minimum")
                else:
                    result.add_warning(f"Low domain confidence: {confidence}")

            if result.errors:
                self.stats['validation_errors'] += 1

            self.stats['domains_validated'] += 1

        except Exception as e:
            result.add_error(f"Exception during domain validation: {str(e)}")
            self.stats['validation_errors'] += 1

        return result

    def get_statistics(self) -> Dict[str, Any]:
        """Get analyzer statistics"""
        return self.stats.copy()

    def get_cache_statistics(self) -> Dict[str, Any]:
        """Get cache statistics"""
        if self.classification_cache:
            return self.classification_cache.get_stats()
        else:
            return {
                'cache_hits': 0,
                'cache_misses': 0,
                'hit_rate': 0.0,
                'domain_classifications_cached': 0,
                'chain_domains_cached': 0
            }

    def clear_cache(self):
        """Clear all caches"""
        if self.classification_cache:
            self.classification_cache.clear()

    def reset_statistics(self):
        """Reset all statistics"""
        for key in self.stats:
            self.stats[key] = 0
