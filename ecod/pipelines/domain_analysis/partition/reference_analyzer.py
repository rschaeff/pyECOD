#!/usr/bin/env python3
"""
Reference domain coverage analysis for domain partitioning.
"""

import logging
from typing import Dict, Any, List, Optional, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict

from ecod.db import DBManager
from .models import ReferenceInfo, EvidenceWithCoverage, Evidence


class ReferenceCoverageAnalyzer:
    """Analyzes and enforces reference domain coverage requirements"""
    
    def __init__(self, db_manager: DBManager, min_coverage: float = 0.7):
        self.db = db_manager
        self.min_coverage = min_coverage
        self.logger = logging.getLogger(__name__)
        
        # Cache for reference domain information
        self._reference_cache: Dict[str, ReferenceInfo] = {}
        self._cache_stats = {'hits': 0, 'misses': 0}
    
    def analyze_evidence_coverage(self, evidence: Evidence) -> EvidenceWithCoverage:
        """
        Analyze reference coverage for a piece of evidence.

        Args:
            evidence: Evidence to analyze

        Returns:
            EvidenceWithCoverage with coverage information
        """
        # Create enhanced evidence object by copying all fields manually
        enhanced = EvidenceWithCoverage(
            type=getattr(evidence, 'type', ''),
            source_id=getattr(evidence, 'source_id', ''),
            domain_id=getattr(evidence, 'domain_id', ''),
            query_range=getattr(evidence, 'query_range', ''),
            hit_range=getattr(evidence, 'hit_range', ''),
            confidence=getattr(evidence, 'confidence', None),
            probability=getattr(evidence, 'probability', None),
            evalue=getattr(evidence, 'evalue', None),
            score=getattr(evidence, 'score', None),
            identity=getattr(evidence, 'identity', None),
            coverage=getattr(evidence, 'coverage', None),
            hsp_count=getattr(evidence, 'hsp_count', None),
            t_group=getattr(evidence, 't_group', None),
            h_group=getattr(evidence, 'h_group', None),
            x_group=getattr(evidence, 'x_group', None),
            a_group=getattr(evidence, 'a_group', None),
            extra_attributes=getattr(evidence, 'extra_attributes', {})
        )
        
        # Skip if no domain_id or hit_range
        if not evidence.domain_id or not evidence.hit_range:
            enhanced.reference_coverage = 0.0
            enhanced.coverage_warning = "No domain_id or hit_range"
            return enhanced
        
        # Get reference domain information
        ref_info = self.get_reference_info(evidence.domain_id)
        if not ref_info:
            enhanced.reference_coverage = 0.0
            enhanced.coverage_warning = f"Reference domain {evidence.domain_id} not found"
            return enhanced
        
        enhanced.reference_info = ref_info
        
        # Calculate coverage
        coverage_result = self._calculate_coverage(evidence.hit_range, ref_info)
        
        enhanced.reference_coverage = coverage_result['coverage']
        enhanced.hit_length = coverage_result['hit_length']
        enhanced.alignment_gaps = coverage_result['gaps']
        
        # Add warnings for low coverage
        if enhanced.reference_coverage < self.min_coverage:
            enhanced.coverage_warning = (
                f"Low reference coverage: {enhanced.reference_coverage:.1%} "
                f"(covers {enhanced.hit_length}/{ref_info.domain_length} residues)"
            )
        
        return enhanced
    
    def get_reference_info(self, domain_id: str) -> Optional[ReferenceInfo]:
        """Get reference domain information from database with caching"""

        # Check cache
        if domain_id in self._reference_cache:
            self._cache_stats['hits'] += 1
            return self._reference_cache[domain_id]

        self._cache_stats['misses'] += 1

        # Query database - join with protein table to get pdb_id
        query = """
        SELECT d.domain_id,
               d.range as domain_range,
               d.length as domain_length,
               d.t_group,
               d.chain_id,
               p.pdb_id,
               CASE
                   WHEN d.range LIKE '%,%' THEN true
                   ELSE false
               END as is_discontinuous
        FROM pdb_analysis.domain d
        JOIN pdb_analysis.protein p ON d.protein_id = p.id
        WHERE d.domain_id = %s
        """

        try:
            results = self.db.execute_dict_query(query, (domain_id,))
            if not results:
                self.logger.warning(f"Reference domain not found: {domain_id}")
                return None

            row = results[0]

            # Parse discontinuous ranges if needed
            discontinuous_ranges = []
            if row['is_discontinuous'] and ',' in row['domain_range']:
                for segment in row['domain_range'].split(','):
                    if '-' in segment:
                        start, end = segment.strip().split('-')
                        discontinuous_ranges.append((int(start), int(end)))

            ref_info = ReferenceInfo(
                domain_id=row['domain_id'],
                domain_range=row['domain_range'],
                domain_length=row['domain_length'],  # Already calculated!
                pdb_id=row['pdb_id'],
                chain_id=row['chain_id'],
                t_group=row.get('t_group'),
                is_discontinuous=row.get('is_discontinuous', False),
                discontinuous_ranges=discontinuous_ranges
            )
            
            # Cache it
            self._reference_cache[domain_id] = ref_info
            
            return ref_info
            
        except Exception as e:
            self.logger.error(f"Error fetching reference info for {domain_id}: {e}")
            return None
    
    def _calculate_coverage(self, hit_range: str, ref_info: ReferenceInfo) -> Dict[str, Any]:
        """Calculate detailed coverage statistics"""
        
        # Parse hit range (can be discontinuous)
        hit_segments = self._parse_range_segments(hit_range)
        if not hit_segments:
            return {'coverage': 0.0, 'hit_length': 0, 'gaps': 0}
        
        # Calculate total hit length
        hit_length = sum(end - start + 1 for start, end in hit_segments)
        
        # For discontinuous references, need more complex calculation
        if ref_info.is_discontinuous and ref_info.discontinuous_ranges:
            coverage = self._calculate_discontinuous_coverage(
                hit_segments, ref_info.discontinuous_ranges
            )
        else:
            # Simple coverage calculation
            coverage = hit_length / ref_info.domain_length if ref_info.domain_length > 0 else 0.0
        
        # Estimate alignment gaps (simplified - would need full alignment for accuracy)
        gaps = 0
        if len(hit_segments) > 1:
            for i in range(len(hit_segments) - 1):
                gap_size = hit_segments[i+1][0] - hit_segments[i][1] - 1
                if gap_size > 0:
                    gaps += gap_size
        
        return {
            'coverage': coverage,
            'hit_length': hit_length,
            'gaps': gaps
        }
    
    def _calculate_discontinuous_coverage(self, hit_segments: List[Tuple[int, int]], 
                                        ref_segments: List[Tuple[int, int]]) -> float:
        """Calculate coverage for discontinuous domains"""
        
        # Calculate overlap between hit and reference segments
        total_overlap = 0
        
        for hit_start, hit_end in hit_segments:
            for ref_start, ref_end in ref_segments:
                # Calculate overlap between this hit and ref segment
                overlap_start = max(hit_start, ref_start)
                overlap_end = min(hit_end, ref_end)
                
                if overlap_start <= overlap_end:
                    total_overlap += overlap_end - overlap_start + 1
        
        # Total reference length
        ref_length = sum(end - start + 1 for start, end in ref_segments)
        
        return total_overlap / ref_length if ref_length > 0 else 0.0
    
    def _parse_range_segments(self, range_str: str) -> List[Tuple[int, int]]:
        """Parse range string into segments"""
        if not range_str:
            return []
        
        segments = []
        
        # Handle discontinuous ranges
        for segment in range_str.split(','):
            segment = segment.strip()
            if '-' in segment:
                try:
                    start, end = segment.split('-')
                    segments.append((int(start), int(end)))
                except ValueError:
                    self.logger.warning(f"Invalid range segment: {segment}")
        
        return segments
    
    def group_evidence_by_reference(self, evidence_list: List[EvidenceWithCoverage]) -> Dict[str, List[EvidenceWithCoverage]]:
        """Group evidence by reference domain for potential combination"""
        
        grouped = defaultdict(list)
        
        for evidence in evidence_list:
            if evidence.domain_id:
                grouped[evidence.domain_id].append(evidence)
        
        return grouped
    
    def calculate_combined_coverage(self, evidence_group: List[EvidenceWithCoverage]) -> float:
        """Calculate combined coverage when multiple alignments hit the same reference"""
        
        if not evidence_group:
            return 0.0
        
        # Get reference info (should be same for all)
        ref_info = evidence_group[0].reference_info
        if not ref_info:
            return 0.0
        
        # Combine all hit ranges
        all_hit_positions = set()
        
        for evidence in evidence_group:
            if evidence.hit_range:
                segments = self._parse_range_segments(evidence.hit_range)
                for start, end in segments:
                    all_hit_positions.update(range(start, end + 1))
        
        # Calculate coverage
        combined_hit_length = len(all_hit_positions)
        combined_coverage = combined_hit_length / ref_info.domain_length
        
        return combined_coverage
    
    def suggest_extended_boundaries(self, evidence: EvidenceWithCoverage, 
                                  query_start: int, query_end: int,
                                  sequence_length: int) -> Tuple[int, int]:
        """
        Suggest extended domain boundaries based on reference domain size.
        
        Args:
            evidence: Evidence with coverage information
            query_start: Current start position on query
            query_end: Current end position on query  
            sequence_length: Total length of query sequence
            
        Returns:
            Tuple of (suggested_start, suggested_end)
        """
        if not evidence.reference_info:
            return (query_start, query_end)
        
        current_length = query_end - query_start + 1
        expected_length = evidence.reference_info.domain_length
        
        # If current assignment is already close to expected size, keep it
        if abs(current_length - expected_length) / expected_length < 0.15:
            return (query_start, query_end)
        
        # Calculate how much we need to extend
        extension_needed = expected_length - current_length
        
        if extension_needed <= 0:
            # Current domain is larger than expected - don't shrink
            return (query_start, query_end)
        
        # Try to extend symmetrically
        extend_left = extension_needed // 2
        extend_right = extension_needed - extend_left
        
        # Apply extensions with boundary checks
        new_start = max(1, query_start - extend_left)
        new_end = min(sequence_length, query_end + extend_right)
        
        # Adjust if we hit sequence boundaries
        actual_left_extension = query_start - new_start
        actual_right_extension = new_end - query_end
        
        # If we couldn't extend enough on one side, try to extend more on the other
        if actual_left_extension < extend_left:
            shortfall = extend_left - actual_left_extension
            new_end = min(sequence_length, new_end + shortfall)
        elif actual_right_extension < extend_right:
            shortfall = extend_right - actual_right_extension  
            new_start = max(1, new_start - shortfall)
        
        self.logger.debug(
            f"Extended domain from {query_start}-{query_end} ({current_length} res) "
            f"to {new_start}-{new_end} ({new_end - new_start + 1} res) "
            f"based on reference length {expected_length}"
        )
        
        return (new_start, new_end)
    
    def get_cache_statistics(self) -> Dict[str, Any]:
        """Get cache statistics"""
        total = self._cache_stats['hits'] + self._cache_stats['misses']
        hit_rate = (self._cache_stats['hits'] / total * 100) if total > 0 else 0
        
        return {
            'reference_cache_size': len(self._reference_cache),
            'cache_hits': self._cache_stats['hits'],
            'cache_misses': self._cache_stats['misses'],
            'cache_hit_rate': hit_rate
        }
