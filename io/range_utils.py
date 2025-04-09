#!/usr/bin/env python3
"""
Range utilities module for ECOD pipeline

Provides utility functions for working with sequence ranges,
calculating overlaps, and manipulating range strings.
"""

import re
import logging
from typing import List, Set, Dict, Tuple, Optional, Union

logger = logging.getLogger("ecod.io.range_utils")


class RangeUtils:
    """
    Utilities for working with sequence ranges
    
    This class provides static methods for manipulating protein sequence
    ranges in string format (e.g., "10-50,60-90") and converting between
    different representations.
    """
    
    @staticmethod
    def get_start_position(range_str: str) -> int:
        """
        Get the start position from a range string
        
        Args:
            range_str: Range string (e.g., "10-50" or "10-20,30-40")
            
        Returns:
            Start position as integer, or 0 if invalid
        """
        if not range_str:
            return 0
            
        if "-" in range_str:
            parts = range_str.split("-")
            try:
                return int(parts[0])
            except ValueError:
                logger.warning(f"Invalid range start in '{range_str}'")
                return 0
        elif "," in range_str:
            # Multi-segment range
            first_segment = range_str.split(",")[0]
            return RangeUtils.get_start_position(first_segment)
        else:
            try:
                return int(range_str)
            except ValueError:
                logger.warning(f"Invalid range format: '{range_str}'")
                return 0
    
    @staticmethod
    def get_end_position(range_str: str) -> int:
        """
        Get the end position from a range string
        
        Args:
            range_str: Range string (e.g., "10-50" or "10-20,30-40")
            
        Returns:
            End position as integer, or 0 if invalid
        """
        if not range_str:
            return 0
            
        if "-" in range_str:
            parts = range_str.split("-")
            try:
                return int(parts[1])
            except (ValueError, IndexError):
                logger.warning(f"Invalid range end in '{range_str}'")
                return 0
        elif "," in range_str:
            # Multi-segment range - get the last segment
            segments = range_str.split(",")
            return RangeUtils.get_end_position(segments[-1])
        else:
            try:
                return int(range_str)
            except ValueError:
                logger.warning(f"Invalid range format: '{range_str}'")
                return 0
    
    @staticmethod
    def range_to_positions(range_str: str) -> Set[int]:
        """
        Convert a range string to a set of positions
        
        Args:
            range_str: Range string (e.g., "10-50" or "10-20,30-40")
            
        Returns:
            Set of integer positions
        """
        positions = set()
        
        if not range_str:
            return positions
            
        # Handle multi-segment ranges
        segments = range_str.split(",")
        
        for segment in segments:
            segment = segment.strip()
            if "-" in segment:
                try:
                    start, end = map(int, segment.split("-"))
                    positions.update(range(start, end + 1))
                except ValueError:
                    logger.warning(f"Invalid range segment: '{segment}'")
                    continue
            else:
                try:
                    positions.add(int(segment))
                except ValueError:
                    logger.warning(f"Invalid position: '{segment}'")
                    continue
                    
        return positions
    
    @staticmethod
    def positions_to_ranges(positions: List[int]) -> str:
        """
        Convert a list of positions to a compact range string
        
        Args:
            positions: List of integer positions
            
        Returns:
            Range string (e.g., "10-20,30-40")
        """
        if not positions:
            return ""
        
        # Sort and remove duplicates
        positions = sorted(set(positions))
        ranges = []
        
        # Initialize with the first position
        start = positions[0]
        prev = start
        
        # Iterate through remaining positions
        for pos in positions[1:]:
            # If there's a gap, end the current range and start a new one
            if pos > prev + 1:
                ranges.append(f"{start}-{prev}")
                start = pos
            prev = pos
        
        # Add the final range
        ranges.append(f"{start}-{prev}")
        
        # Special case: single position ranges (convert "5-5" to "5")
        result_ranges = []
        for r in ranges:
            if r.startswith(r.split("-")[0] + "-" + r.split("-")[0]):
                result_ranges.append(r.split("-")[0])
            else:
                result_ranges.append(r)
        
        return ",".join(result_ranges)
    
    @staticmethod
    def calculate_overlap_percentage(range1: str, range2: str) -> float:
        """
        Calculate the percentage of overlap between two ranges
        
        Args:
            range1: First range string
            range2: Second range string
            
        Returns:
            Overlap percentage (0.0 to 1.0)
        """
        # Convert ranges to sets of positions
        positions1 = RangeUtils.range_to_positions(range1)
        positions2 = RangeUtils.range_to_positions(range2)
        
        # Calculate overlap
        if not positions1 or not positions2:
            return 0.0
            
        overlap = len(positions1.intersection(positions2))
        
        # Calculate coverage relative to the smaller range
        min_length = min(len(positions1), len(positions2))
        if min_length == 0:
            return 0.0
            
        return overlap / min_length
    
    @staticmethod
    def calculate_coverage_percentage(range_str: str, sequence_length: int) -> float:
        """
        Calculate the percentage of a sequence covered by a range
        
        Args:
            range_str: Range string
            sequence_length: Total length of the sequence
            
        Returns:
            Coverage percentage (0.0 to 1.0)
        """
        if not range_str or sequence_length <= 0:
            return 0.0
            
        positions = RangeUtils.range_to_positions(range_str)
        
        return len(positions) / sequence_length
    
    @staticmethod
    def merge_overlapping_ranges(ranges: List[str], gap_tolerance: int = 0) -> str:
        """
        Merge overlapping ranges into a single range string
        
        Args:
            ranges: List of range strings to merge
            gap_tolerance: Maximum gap size to consider as continuous (0 means exact adjacency)
            
        Returns:
            Merged range string
        """
        if not ranges:
            return ""
            
        # Convert all ranges to position sets and merge
        all_positions = set()
        for range_str in ranges:
            positions = RangeUtils.range_to_positions(range_str)
            all_positions.update(positions)
            
        # If gap tolerance > 0, fill in small gaps
        if gap_tolerance > 0:
            # Sort positions
            sorted_positions = sorted(all_positions)
            
            # Find gaps
            gaps = []
            for i in range(len(sorted_positions) - 1):
                current_pos = sorted_positions[i]
                next_pos = sorted_positions[i + 1]
                gap_size = next_pos - current_pos - 1
                
                if 0 < gap_size <= gap_tolerance:
                    # Fill gap
                    for pos in range(current_pos + 1, next_pos):
                        all_positions.add(pos)
        
        # Convert back to range string
        return RangeUtils.positions_to_ranges(sorted(all_positions))
    
    @staticmethod
    def split_ranges(range_str: str, max_length: int) -> List[str]:
        """
        Split a range into multiple ranges with maximum length
        
        Args:
            range_str: Range string to split
            max_length: Maximum length of each segment
            
        Returns:
            List of range strings
        """
        if not range_str or max_length <= 0:
            return []
            
        positions = sorted(RangeUtils.range_to_positions(range_str))
        
        # Split positions into chunks of max_length
        chunks = []
        current_chunk = []
        
        for pos in positions:
            current_chunk.append(pos)
            if len(current_chunk) >= max_length:
                chunks.append(current_chunk)
                current_chunk = []
                
        # Add remaining positions
        if current_chunk:
            chunks.append(current_chunk)
            
        # Convert each chunk back to range string
        return [RangeUtils.positions_to_ranges(chunk) for chunk in chunks]
    
    @staticmethod
    def intersect_ranges(range1: str, range2: str) -> str:
        """
        Find the intersection of two ranges
        
        Args:
            range1: First range string
            range2: Second range string
            
        Returns:
            Range string representing the intersection
        """
        positions1 = RangeUtils.range_to_positions(range1)
        positions2 = RangeUtils.range_to_positions(range2)
        
        intersection = positions1.intersection(positions2)
        
        return RangeUtils.positions_to_ranges(sorted(intersection))
    
    @staticmethod
    def subtract_ranges(range1: str, range2: str) -> str:
        """
        Subtract range2 from range1
        
        Args:
            range1: Range string to subtract from
            range2: Range string to subtract
            
        Returns:
            Range string representing range1 - range2
        """
        positions1 = RangeUtils.range_to_positions(range1)
        positions2 = RangeUtils.range_to_positions(range2)
        
        difference = positions1.difference(positions2)
        
        return RangeUtils.positions_to_ranges(sorted(difference))
    
    @staticmethod
    def validate_range(range_str: str, sequence_length: int) -> bool:
        """
        Validate that a range is well-formed and within sequence bounds
        
        Args:
            range_str: Range string to validate
            sequence_length: Length of the sequence
            
        Returns:
            True if valid, False otherwise
        """
        if not range_str:
            return False
            
        try:
            positions = RangeUtils.range_to_positions(range_str)
            
            # Check that positions are within bounds
            return all(1 <= pos <= sequence_length for pos in positions)
            
        except Exception as e:
            logger.warning(f"Range validation error: {e}")
            return False
    
    @staticmethod
    def range_length(range_str: str) -> int:
        """
        Calculate the number of positions in a range
        
        Args:
            range_str: Range string
            
        Returns:
            Number of positions
        """
        if not range_str:
            return 0
            
        positions = RangeUtils.range_to_positions(range_str)
        return len(positions)
    
    @staticmethod
    def normalize_range(range_str: str) -> str:
        """
        Normalize a range string (sort, merge overlaps, etc.)
        
        Args:
            range_str: Range string to normalize
            
        Returns:
            Normalized range string
        """
        if not range_str:
            return ""
            
        # Convert to positions and back to get a normalized representation
        positions = RangeUtils.range_to_positions(range_str)
        return RangeUtils.positions_to_ranges(sorted(positions))
    
    @staticmethod
    def format_chain_range(chain_id: str, range_str: str) -> str:
        """
        Format a range with chain identifier
        
        Args:
            chain_id: Chain identifier
            range_str: Range string
            
        Returns:
            Chain range string (e.g., "A:10-50")
        """
        if not chain_id or not range_str:
            return range_str
            
        return f"{chain_id}:{range_str}"
    
    @staticmethod
    def parse_chain_range(chain_range: str) -> Tuple[str, str]:
        """
        Parse a chain range string into chain_id and range
        
        Args:
            chain_range: Chain range string (e.g., "A:10-50")
            
        Returns:
            Tuple of (chain_id, range_str)
        """
        if not chain_range:
            return ("", "")
            
        if ":" not in chain_range:
            # No chain identifier
            return ("", chain_range)
            
        parts = chain_range.split(":", 1)
        return (parts[0], parts[1])