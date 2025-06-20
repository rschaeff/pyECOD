# mini/core/boundary_optimizer.py (simplified with SequenceRange integration)
"""
Boundary optimization with fragment merging rules

Simplified implementation using SequenceRange for all position operations.
Implements the boundary optimization logic from the Perl boundary_optimize routine.
Follows mini isolation rules - no external dependencies, pure sequence-based optimization.
"""

from typing import List, Set, Tuple, Optional, Protocol
from .models import Domain, UnassignedSegment, DomainLayout, SegmentType, FragmentSize
from .sequence_range import SequenceRange


class ContactAnalyzer(Protocol):
    """Protocol for contact analysis - allows future extensibility"""
    
    def contact_boundary_resolution(self, domain1: Domain, domain2: Domain, 
                                  segment: UnassignedSegment) -> Tuple[Set[int], Set[int]]:
        """
        Resolve boundary between two domains using contact analysis.
        
        Args:
            domain1: First domain
            domain2: Second domain  
            segment: Intervening segment to split
            
        Returns:
            Tuple of (positions_to_domain1, positions_to_domain2)
        """
        ...


class SequenceContactAnalyzer:
    """
    Stub implementation of contact analysis using sequence-based heuristics.
    
    This is a placeholder that follows mini isolation rules.
    Future implementations could add structural analysis.
    """
    
    def contact_boundary_resolution(self, domain1: Domain, domain2: Domain,
                                  segment: UnassignedSegment) -> Tuple[Set[int], Set[int]]:
        """
        Simple sequence-based boundary resolution - split segment in half.
        
        This is a fallback when structural contact analysis is not available.
        Real contact analysis would use CA-CA distances from 3D coordinates.
        """
        segment_positions = sorted(segment.positions)
        midpoint = len(segment_positions) // 2
        
        # Split roughly in half, with slight bias toward first domain
        to_domain1 = set(segment_positions[:midpoint + 1])
        to_domain2 = set(segment_positions[midpoint + 1:])
        
        # Handle single residue case
        if len(segment_positions) == 1:
            # Give to domain with higher confidence if available
            if (domain1.evidence_items and domain2.evidence_items):
                confidence1 = domain1.evidence_items[0].confidence
                confidence2 = domain2.evidence_items[0].confidence
                if confidence1 > confidence2:
                    to_domain1 = segment.positions
                    to_domain2 = set()
                else:
                    to_domain1 = set()
                    to_domain2 = segment.positions
            else:
                # Default to first domain
                to_domain1 = segment.positions
                to_domain2 = set()
        
        return to_domain1, to_domain2


class BoundaryOptimizer:
    """
    Boundary optimization using fragment merging rules.
    
    Simplified implementation using SequenceRange for all position operations.
    Implements the boundary optimization algorithm from the Perl code.
    """
    
    def __init__(self, contact_analyzer: Optional[ContactAnalyzer] = None):
        """
        Initialize boundary optimizer.
        
        Args:
            contact_analyzer: Optional contact analyzer for structural boundary resolution
        """
        self.contact_analyzer = contact_analyzer or SequenceContactAnalyzer()
        self.optimization_stats = {
            'nterm_merges': 0,
            'cterm_merges': 0,
            'inter_single_merges': 0,
            'inter_split_merges': 0,
            'singleton_merges': 0,
            'large_gaps_skipped': 0,
            'overlaps_resolved': 0
        }
    
    def optimize_boundaries(self, layout: DomainLayout, min_domain_size: int = 25,
                          neighbor_tolerance: int = 5, verbose: bool = False) -> DomainLayout:
        """
        Apply boundary optimization to a domain layout.
        
        Args:
            layout: Domain layout to optimize
            min_domain_size: Size threshold for fragment classification
            neighbor_tolerance: Distance tolerance for neighbor detection
            verbose: Whether to print optimization details
            
        Returns:
            Optimized domain layout
        """
        if verbose:
            print(f"\nSTEP: BOUNDARY OPTIMIZATION")
            print("=" * 40)
        
        # Reset stats
        self._reset_stats()
        
        # Step 1: Resolve small overlaps between domains
        overlaps_resolved = layout.resolve_small_overlaps(max_overlap_size=5)
        self.optimization_stats['overlaps_resolved'] = overlaps_resolved
        
        if overlaps_resolved > 0 and verbose:
            print(f"Resolved {overlaps_resolved} small domain overlaps")
        
        # Step 2: Analyze gaps after overlap resolution
        layout.analyze_gaps(min_domain_size)
        
        if not layout.unassigned_segments:
            if verbose:
                print("No unassigned segments found - optimization complete")
            return layout
        
        # Step 3: Apply fragment merging rules
        self._apply_fragment_merging_rules(layout, verbose)
        
        # Step 4: Final gap analysis to update statistics
        layout.analyze_gaps(min_domain_size)
        
        if verbose:
            self._print_optimization_summary(layout)
        
        return layout
    
    def _apply_fragment_merging_rules(self, layout: DomainLayout, verbose: bool = False) -> None:
        """Apply the fragment merging rules from the Perl boundary optimization"""
        
        # Process segments by type and size priority
        # Process small fragments first, then handle large gaps
        small_fragments = [s for s in layout.unassigned_segments if s.is_small_fragment]
        large_gaps = [s for s in layout.unassigned_segments if not s.is_small_fragment]
        
        if verbose:
            print(f"Processing {len(small_fragments)} small fragments and {len(large_gaps)} large gaps")
        
        # Apply rules to small fragments
        for segment in small_fragments:
            self._process_small_fragment(segment, layout, verbose)
        
        # Large gaps are left for external domain parser (not implemented in mini)
        for segment in large_gaps:
            self.optimization_stats['large_gaps_skipped'] += 1
            if verbose:
                print(f"Skipping large gap {segment.start}-{segment.end} (length {segment.length}) - "
                      f"would require external domain parser")
    
    def _process_small_fragment(self, segment: UnassignedSegment, layout: DomainLayout, 
                              verbose: bool = False) -> None:
        """Process a single small fragment according to optimization rules"""
        
        if segment.segment_type == SegmentType.NTERM:
            self._merge_nterm_fragment(segment, layout, verbose)
            
        elif segment.segment_type == SegmentType.CTERM:
            self._merge_cterm_fragment(segment, layout, verbose)
            
        elif segment.segment_type == SegmentType.INTERSTITIAL:
            self._merge_interstitial_fragment(segment, layout, verbose)
            
        elif segment.segment_type == SegmentType.SINGLETON:
            self._merge_singleton(segment, layout, verbose)
    
    def _merge_nterm_fragment(self, segment: UnassignedSegment, layout: DomainLayout,
                            verbose: bool = False) -> None:
        """Merge N-terminal fragment with first domain"""
        
        if not layout.domains:
            if verbose:
                print(f"Warning: No domains available for N-terminal fragment {segment.start}-{segment.end}")
            return
        
        # Find first domain (by start position)
        first_domain = min(layout.domains, key=lambda d: d.start_position)
        
        if verbose:
            print(f"Merging N-terminal fragment {segment.start}-{segment.end} with domain {first_domain.id}")

        first_domain.record_optimization_action('nterm_merge', f"{segment.start}-{segment.end}")
        
        layout.merge_segment_with_domain(segment, first_domain)
        self.optimization_stats['nterm_merges'] += 1
    
    def _merge_cterm_fragment(self, segment: UnassignedSegment, layout: DomainLayout,
                            verbose: bool = False) -> None:
        """Merge C-terminal fragment with last domain"""
        
        if not layout.domains:
            if verbose:
                print(f"Warning: No domains available for C-terminal fragment {segment.start}-{segment.end}")
            return
        
        # Find last domain (by end position)
        last_domain = max(layout.domains, key=lambda d: d.end_position)
        
        if verbose:
            print(f"Merging C-terminal fragment {segment.start}-{segment.end} with domain {last_domain.id}")

        last_domain.record_optimization_action("cterm_merge", f"{segment.start}-{segment.end}")
        
        layout.merge_segment_with_domain(segment, last_domain)
        self.optimization_stats['cterm_merges'] += 1
    
    def _merge_interstitial_fragment(self, segment: UnassignedSegment, layout: DomainLayout,
                                   verbose: bool = False) -> None:
        """Merge interstitial fragment according to neighbor rules"""
        
        preceding = segment.preceding_domain
        following = segment.following_domain
        
        if preceding and following:
            # Fragment between two domains - split using contact analysis
            if verbose:
                print(f"Splitting interstitial fragment {segment.start}-{segment.end} "
                      f"between domains {preceding.id} and {following.id}")
            
            split_positions = self.contact_analyzer.contact_boundary_resolution(
                preceding, following, segment)
            
            layout.split_segment_between_domains(segment, preceding, following, split_positions)
            self.optimization_stats['inter_split_merges'] += 1
            
        elif preceding:
            # Fragment follows a domain - merge with preceding domain
            if verbose:
                print(f"Merging interstitial fragment {segment.start}-{segment.end} "
                      f"with preceding domain {preceding.id}")
            
            layout.merge_segment_with_domain(segment, preceding)
            self.optimization_stats['inter_single_merges'] += 1
            
        elif following:
            # Fragment precedes a domain - merge with following domain
            if verbose:
                print(f"Merging interstitial fragment {segment.start}-{segment.end} "
                      f"with following domain {following.id}")
            
            layout.merge_segment_with_domain(segment, following)
            self.optimization_stats['inter_single_merges'] += 1
            
        else:
            # No adjacent domains within tolerance - treat as N-terminal
            if verbose:
                print(f"Interstitial fragment {segment.start}-{segment.end} has no close neighbors, "
                      f"treating as N-terminal")
            self._merge_nterm_fragment(segment, layout, verbose)
    
    def _merge_singleton(self, segment: UnassignedSegment, layout: DomainLayout,
                       verbose: bool = False) -> None:
        """Merge singleton residue with closest domain"""
        
        if not layout.domains:
            if verbose:
                print(f"Warning: No domains available for singleton {segment.start}")
            return
        
        singleton_pos = segment.start  # start == end for singletons
        
        # Find closest domain
        closest_domain = None
        min_distance = float('inf')
        
        for domain in layout.domains:
            distance = min(
                abs(singleton_pos - domain.start_position),
                abs(singleton_pos - domain.end_position)
            )
            
            # Bias toward N-terminal domain in case of ties
            if distance < min_distance or (distance == min_distance and 
                                         (closest_domain is None or 
                                          domain.start_position < closest_domain.start_position)):
                min_distance = distance
                closest_domain = domain
        
        if closest_domain:
            if verbose:
                print(f"Merging singleton {singleton_pos} with closest domain {closest_domain.id} "
                      f"(distance: {min_distance})")
            
            layout.merge_segment_with_domain(segment, closest_domain)
            self.optimization_stats['singleton_merges'] += 1
    
    def _reset_stats(self) -> None:
        """Reset optimization statistics"""
        for key in self.optimization_stats:
            self.optimization_stats[key] = 0
    
    def _print_optimization_summary(self, layout: DomainLayout) -> None:
        """Print summary of optimization results"""
        stats = layout.get_coverage_stats()
        
        print(f"\nBoundary optimization summary:")
        print(f"  Final coverage: {stats['assigned_residues']}/{stats['total_residues']} "
              f"residues ({stats['coverage_percent']:.1f}%)")
        print(f"  Final domains: {stats['num_domains']}")
        print(f"  Remaining gaps: {stats['num_gaps']}")
        
        print(f"\nOptimization actions:")
        for action, count in self.optimization_stats.items():
            if count > 0:
                action_name = action.replace('_', ' ').title()
                print(f"  {action_name}: {count}")
    
    def get_optimization_stats(self) -> dict:
        """Get optimization statistics"""
        return self.optimization_stats.copy()


# Convenience function for standalone boundary optimization
def optimize_domain_boundaries(domains: List[Domain], sequence_length: int,
                             min_domain_size: int = 25, neighbor_tolerance: int = 5,
                             contact_analyzer: Optional[ContactAnalyzer] = None,
                             verbose: bool = False) -> Tuple[List[Domain], dict]:
    """
    Optimize domain boundaries using fragment merging rules.
    
    Args:
        domains: List of domains to optimize
        sequence_length: Total protein sequence length
        min_domain_size: Size threshold for fragment classification
        neighbor_tolerance: Distance tolerance for neighbor detection
        contact_analyzer: Optional contact analyzer for boundary resolution
        verbose: Whether to print optimization details
        
    Returns:
        Tuple of (optimized_domains, optimization_stats)
    """
    # Create layout
    layout = DomainLayout.from_domains(domains, sequence_length)
    
    # Optimize
    optimizer = BoundaryOptimizer(contact_analyzer)
    optimized_layout = optimizer.optimize_boundaries(
        layout, min_domain_size, neighbor_tolerance, verbose)
    
    return optimized_layout.domains, optimizer.get_optimization_stats()


# Enhanced contact analyzer interface for future structural extensions
class StructuralContactAnalyzer:
    """
    Future interface for structural contact analysis.
    
    This would replace the sequence-based stub when structural analysis
    is added to mini (requires rule modification for external dependencies).
    """
    
    def __init__(self, structure_path: Optional[str] = None):
        """
        Initialize with optional structure file.
        
        Args:
            structure_path: Path to PDB/mmCIF structure file
        """
        self.structure_path = structure_path
        # Future: load structure and parse coordinates
    
    def contact_boundary_resolution(self, domain1: Domain, domain2: Domain,
                                  segment: UnassignedSegment) -> Tuple[Set[int], Set[int]]:
        """
        Structural contact-based boundary resolution.
        
        This would implement the CA-CA distance calculation from the Perl code.
        Currently returns sequence-based fallback.
        """
        # Future implementation would:
        # 1. Load coordinates for domain1, domain2, segment residues
        # 2. Calculate CA-CA distances < cutoff (typically 8-10 Å) 
        # 3. For each segment residue, count contacts to each domain
        # 4. Find optimal cut point that maximizes contact satisfaction
        # 5. Return split positions
        
        # For now, fall back to sequence-based split
        fallback = SequenceContactAnalyzer()
        return fallback.contact_boundary_resolution(domain1, domain2, segment)
    
    def load_structure(self, structure_path: str) -> bool:
        """
        Load structure file for contact analysis.
        
        Args:
            structure_path: Path to structure file
            
        Returns:
            True if successful, False otherwise
        """
        # Future implementation
        self.structure_path = structure_path
        return False  # Not implemented yet
