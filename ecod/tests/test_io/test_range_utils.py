#!/usr/bin/env python3
"""
Unit tests for Range Utilities module
"""
import unittest
from typing import Set, List

from ecod.io.range_utils import RangeUtils


class TestRangeUtils(unittest.TestCase):
    """Test cases for RangeUtils class"""
    
    def test_get_start_position(self):
        """Test get_start_position method"""
        # Test simple range
        self.assertEqual(RangeUtils.get_start_position("10-50"), 10)
        
        # Test multi-segment range
        self.assertEqual(RangeUtils.get_start_position("10-20,30-40"), 10)
        
        # Test single position
        self.assertEqual(RangeUtils.get_start_position("15"), 15)
        
        # Test empty string
        self.assertEqual(RangeUtils.get_start_position(""), 0)
        
        # Test invalid range
        self.assertEqual(RangeUtils.get_start_position("invalid"), 0)
    
    def test_get_end_position(self):
        """Test get_end_position method"""
        # Test simple range
        self.assertEqual(RangeUtils.get_end_position("10-50"), 50)
        
        # Test multi-segment range
        self.assertEqual(RangeUtils.get_end_position("10-20,30-40"), 40)
        
        # Test single position
        self.assertEqual(RangeUtils.get_end_position("15"), 15)
        
        # Test empty string
        self.assertEqual(RangeUtils.get_end_position(""), 0)
        
        # Test invalid range
        self.assertEqual(RangeUtils.get_end_position("invalid"), 0)
    
    def test_range_to_positions(self):
        """Test range_to_positions method"""
        # Test simple range
        self.assertEqual(RangeUtils.range_to_positions("10-15"), {10, 11, 12, 13, 14, 15})
        
        # Test multi-segment range
        self.assertEqual(RangeUtils.range_to_positions("10-12,14-16"), {10, 11, 12, 14, 15, 16})
        
        # Test single position
        self.assertEqual(RangeUtils.range_to_positions("10"), {10})
        
        # Test empty string
        self.assertEqual(RangeUtils.range_to_positions(""), set())
        
        # Test invalid range
        self.assertEqual(RangeUtils.range_to_positions("invalid"), set())
        
        # Test mixed valid and invalid segments
        self.assertEqual(RangeUtils.range_to_positions("10-12,invalid,14-16"), {10, 11, 12, 14, 15, 16})
    
    def test_positions_to_ranges(self):
        """Test positions_to_ranges method"""
        # Test simple consecutive positions
        self.assertEqual(RangeUtils.positions_to_ranges([10, 11, 12, 13, 14, 15]), "10-15")
        
        # Test non-consecutive positions
        self.assertEqual(RangeUtils.positions_to_ranges([10, 11, 12, 14, 15, 16]), "10-12,14-16")
        
        # Test single position
        self.assertEqual(RangeUtils.positions_to_ranges([10]), "10")
        
        # Test empty list
        self.assertEqual(RangeUtils.positions_to_ranges([]), "")
        
        # Test unsorted positions
        self.assertEqual(RangeUtils.positions_to_ranges([15, 10, 11, 14, 12, 16]), "10-12,14-16")
        
        # Test duplicate positions
        self.assertEqual(RangeUtils.positions_to_ranges([10, 11, 11, 12, 14, 15, 15, 16]), "10-12,14-16")
    
    def test_calculate_overlap_percentage(self):
        """Test calculate_overlap_percentage method"""
        # Test complete overlap
        self.assertEqual(RangeUtils.calculate_overlap_percentage("10-20", "10-20"), 1.0)
        
        # Test partial overlap
        self.assertEqual(RangeUtils.calculate_overlap_percentage("10-20", "15-25"), 0.5)
        
        # Test no overlap
        self.assertEqual(RangeUtils.calculate_overlap_percentage("10-20", "30-40"), 0.0)
        
        # Test one range as a subset of the other
        self.assertEqual(RangeUtils.calculate_overlap_percentage("10-30", "15-25"), 0.5)
        
        # Test with multi-segment ranges
        self.assertEqual(RangeUtils.calculate_overlap_percentage("10-20,30-40", "15-35"), 0.5)
        
        # Test with empty ranges
        self.assertEqual(RangeUtils.calculate_overlap_percentage("", "10-20"), 0.0)
        self.assertEqual(RangeUtils.calculate_overlap_percentage("10-20", ""), 0.0)
        self.assertEqual(RangeUtils.calculate_overlap_percentage("", ""), 0.0)
    
    def test_calculate_coverage_percentage(self):
        """Test calculate_coverage_percentage method"""
        # Test full coverage
        self.assertEqual(RangeUtils.calculate_coverage_percentage("1-100", 100), 1.0)
        
        # Test partial coverage
        self.assertEqual(RangeUtils.calculate_coverage_percentage("10-20", 100), 0.11)
        
        # Test multi-segment coverage
        self.assertEqual(RangeUtils.calculate_coverage_percentage("10-20,30-40", 100), 0.22)
        
        # Test empty range
        self.assertEqual(RangeUtils.calculate_coverage_percentage("", 100), 0.0)
        
        # Test invalid sequence length
        self.assertEqual(RangeUtils.calculate_coverage_percentage("10-20", 0), 0.0)
    
    def test_merge_overlapping_ranges(self):
        """Test merge_overlapping_ranges method"""
        # Test overlapping ranges
        self.assertEqual(RangeUtils.merge_overlapping_ranges(["10-20", "15-25"]), "10-25")
        
        # Test adjacent ranges
        self.assertEqual(RangeUtils.merge_overlapping_ranges(["10-20", "21-30"]), "10-30")
        
        # Test ranges with gaps
        self.assertEqual(RangeUtils.merge_overlapping_ranges(["10-20", "25-35"]), "10-20,25-35")
        
        # Test with gap tolerance
        self.assertEqual(RangeUtils.merge_overlapping_ranges(["10-20", "25-35"], gap_tolerance=5), "10-35")
        
        # Test empty list
        self.assertEqual(RangeUtils.merge_overlapping_ranges([]), "")
        
        # Test single range
        self.assertEqual(RangeUtils.merge_overlapping_ranges(["10-20"]), "10-20")
        
        # Test complex case
        self.assertEqual(
            RangeUtils.merge_overlapping_ranges(["10-20", "15-25", "30-40", "50-60"]),
            "10-25,30-40,50-60"
        )
    
    def test_split_ranges(self):
        """Test split_ranges method"""
        # Test splitting a range
        self.assertEqual(RangeUtils.split_ranges("1-10", 5), ["1-5", "6-10"])
        
        # Test range smaller than max_length
        self.assertEqual(RangeUtils.split_ranges("1-5", 10), ["1-5"])
        
        # Test multi-segment range
        self.assertEqual(RangeUtils.split_ranges("1-5,10-15", 5), ["1-5", "10-14", "15"])
        
        # Test empty range
        self.assertEqual(RangeUtils.split_ranges("", 5), [])
        
        # Test invalid max_length
        self.assertEqual(RangeUtils.split_ranges("1-10", 0), [])
    
    def test_intersect_ranges(self):
        """Test intersect_ranges method"""
        # Test overlapping ranges
        self.assertEqual(RangeUtils.intersect_ranges("10-20", "15-25"), "15-20")
        
        # Test no overlap
        self.assertEqual(RangeUtils.intersect_ranges("10-20", "30-40"), "")
        
        # Test exact match
        self.assertEqual(RangeUtils.intersect_ranges("10-20", "10-20"), "10-20")
        
        # Test partial overlap with multi-segment
        self.assertEqual(RangeUtils.intersect_ranges("10-20,30-40", "15-35"), "15-20,30-35")
        
        # Test empty ranges
        self.assertEqual(RangeUtils.intersect_ranges("", "10-20"), "")
        self.assertEqual(RangeUtils.intersect_ranges("10-20", ""), "")
    
    def test_subtract_ranges(self):
        """Test subtract_ranges method"""
        # Test partial overlap
        self.assertEqual(RangeUtils.subtract_ranges("10-20", "15-25"), "10-14")
        
        # Test no overlap
        self.assertEqual(RangeUtils.subtract_ranges("10-20", "30-40"), "10-20")
        
        # Test exact match
        self.assertEqual(RangeUtils.subtract_ranges("10-20", "10-20"), "")
        
        # Test middle section removal
        self.assertEqual(RangeUtils.subtract_ranges("10-30", "15-25"), "10-14,26-30")
        
        # Test with multi-segment ranges
        self.assertEqual(RangeUtils.subtract_ranges("10-20,30-40", "15-35"), "10-14,36-40")
        
        # Test empty ranges
        self.assertEqual(RangeUtils.subtract_ranges("", "10-20"), "")
        self.assertEqual(RangeUtils.subtract_ranges("10-20", ""), "10-20")
    
    def test_validate_range(self):
        """Test validate_range method"""
        # Test valid range
        self.assertTrue(RangeUtils.validate_range("10-20", 100))
        
        # Test multi-segment valid range
        self.assertTrue(RangeUtils.validate_range("10-20,30-40", 100))
        
        # Test range exceeding sequence length
        self.assertFalse(RangeUtils.validate_range("10-120", 100))
        
        # Test negative positions
        self.assertFalse(RangeUtils.validate_range("-10-20", 100))
        
        # Test invalid range format
        self.assertFalse(RangeUtils.validate_range("invalid", 100))
        
        # Test empty range
        self.assertFalse(RangeUtils.validate_range("", 100))
    
    def test_range_length(self):
        """Test range_length method"""
        # Test simple range
        self.assertEqual(RangeUtils.range_length("10-15"), 6)
        
        # Test multi-segment range
        self.assertEqual(RangeUtils.range_length("10-12,14-16"), 6)
        
        # Test single position
        self.assertEqual(RangeUtils.range_length("10"), 1)
        
        # Test empty range
        self.assertEqual(RangeUtils.range_length(""), 0)
        
        # Test invalid range
        self.assertEqual(RangeUtils.range_length("invalid"), 0)
    
    def test_normalize_range(self):
        """Test normalize_range method"""
        # Test already normalized range
        self.assertEqual(RangeUtils.normalize_range("10-20"), "10-20")
        
        # Test unsorted segments
        self.assertEqual(RangeUtils.normalize_range("30-40,10-20"), "10-20,30-40")
        
        # Test overlapping segments
        self.assertEqual(RangeUtils.normalize_range("10-25,20-30"), "10-30")
        
        # Test duplicate positions
        self.assertEqual(RangeUtils.normalize_range("10-20,15-25"), "10-25")
        
        # Test empty range
        self.assertEqual(RangeUtils.normalize_range(""), "")
    
    def test_format_chain_range(self):
        """Test format_chain_range method"""
        # Test simple case
        self.assertEqual(RangeUtils.format_chain_range("A", "10-20"), "A:10-20")
        
        # Test multi-segment range
        self.assertEqual(RangeUtils.format_chain_range("A", "10-20,30-40"), "A:10-20,30-40")
        
        # Test empty chain
        self.assertEqual(RangeUtils.format_chain_range("", "10-20"), "10-20")
        
        # Test empty range
        self.assertEqual(RangeUtils.format_chain_range("A", ""), "")
    
    def test_parse_chain_range(self):
        """Test parse_chain_range method"""
        # Test simple case
        self.assertEqual(RangeUtils.parse_chain_range("A:10-20"), ("A", "10-20"))
        
        # Test multi-segment range
        self.assertEqual(RangeUtils.parse_chain_range("A:10-20,30-40"), ("A", "10-20,30-40"))
        
        # Test no chain separator
        self.assertEqual(RangeUtils.parse_chain_range("10-20"), ("", "10-20"))
        
        # Test empty string
        self.assertEqual(RangeUtils.parse_chain_range(""), ("", ""))
        
        # Test multiple colons (should only split on first one)
        self.assertEqual(RangeUtils.parse_chain_range("A:10-20:30-40"), ("A", "10-20:30-40"))


if __name__ == '__main__':
    unittest.main()