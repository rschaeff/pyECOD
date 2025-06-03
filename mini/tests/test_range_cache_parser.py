#!/usr/bin/env python3
"""
Range cache parser tests for mini_pyecod

Tests parsing of ECOD range cache files.
"""

import pytest
from pathlib import Path
import csv

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.core.range_cache_parser import (
    parse_range_cache,
    create_domain_definitions_from_cache,
    extract_protein_lengths_from_cache,
    RangeCacheEntry
)
from ecod.core.sequence_range import SequenceRange


class TestRangeCacheParsing:
    """Test parsing ECOD range cache files"""
    
    @pytest.mark.unit
    def test_parse_basic_range_cache(self, tmp_path):
        """Test parsing basic range cache entries"""
        cache_file = tmp_path / "range_cache.txt"
        
        # Create test cache file (tab-separated)
        with open(cache_file, 'w') as f:
            f.write("001520984\te2ia4A1\tA:110-209\t2ia4\tA\n")
            f.write("001520985\te2ia4A2\tA:210-308\t2ia4\tA\n")
            f.write("001234567\te6dgvA1\tA:1-238\t6dgv\tA\n")
        
        entries = parse_range_cache(str(cache_file))
        
        assert len(entries) == 3
        
        # Check first entry
        assert 'e2ia4A1' in entries
        entry = entries['e2ia4A1']
        assert entry.cache_id == "001520984"
        assert entry.domain_id == "e2ia4A1"
        assert entry.range_spec == "A:110-209"
        assert entry.pdb_id == "2ia4"
        assert entry.chain_id == "A"
        assert str(entry.range) == "110-209"
        assert entry.length == 100
    
    @pytest.mark.unit
    def test_parse_discontinuous_ranges(self, tmp_path):
        """Test parsing discontinuous ranges"""
        cache_file = tmp_path / "disco_cache.txt"
        
        with open(cache_file, 'w') as f:
            # Discontinuous range
            f.write("123456\tediscoA1\tA:10-50,A:100-150\tdisco\tA\n")
        
        entries = parse_range_cache(str(cache_file))
        
        assert len(entries) == 1
        entry = entries['ediscoA1']
        assert entry.range.is_discontinuous
        assert str(entry.range) == "10-50,100-150"
        assert entry.length == 92  # 41 + 52
    
    @pytest.mark.unit
    def test_parse_without_chain_prefix(self, tmp_path):
        """Test parsing ranges without chain prefix"""
        cache_file = tmp_path / "no_prefix.txt"
        
        with open(cache_file, 'w') as f:
            # Range without chain prefix
            f.write("789\tesimpleA1\t1-100\tsimple\tA\n")
        
        entries = parse_range_cache(str(cache_file))
        
        assert len(entries) == 1
        entry = entries['esimpleA1']
        assert str(entry.range) == "1-100"
        assert entry.length == 100
    
    @pytest.mark.unit
    def test_parse_invalid_lines(self, tmp_path):
        """Test handling of invalid lines"""
        cache_file = tmp_path / "invalid.txt"
        
        with open(cache_file, 'w') as f:
            f.write("good1\tegoodA1\tA:1-50\tgood\tA\n")
            f.write("invalid\tline\twith\ttoo\tfew\tcolumns\textra\n")  # Too many columns
            f.write("short\tline\n")  # Too few columns
            f.write("bad2\tebadA1\tinvalid-range\tbad\tA\n")  # Invalid range
            f.write("good2\tegoodB1\tB:51-100\tgood\tB\n")
            f.write("\n")  # Empty line
        
        entries = parse_range_cache(str(cache_file), verbose=True)
        
        # Should only have the two good entries
        assert len(entries) == 2
        assert 'egoodA1' in entries
        assert 'egoodB1' in entries
        assert 'ebadA1' not in entries
    
    @pytest.mark.unit
    def test_empty_cache_file(self, tmp_path):
        """Test empty cache file"""
        cache_file = tmp_path / "empty.txt"
        cache_file.touch()  # Create empty file
        
        entries = parse_range_cache(str(cache_file))
        assert len(entries) == 0
    
    @pytest.mark.unit
    def test_missing_cache_file(self):
        """Test missing cache file"""
        entries = parse_range_cache("/nonexistent/cache.txt")
        assert len(entries) == 0


class TestDomainDefinitionGeneration:
    """Test generating domain definitions from cache"""
    
    @pytest.mark.unit
    def test_create_domain_definitions(self, tmp_path):
        """Test creating domain definitions CSV from cache"""
        # Create cache file
        cache_file = tmp_path / "cache.txt"
        with open(cache_file, 'w') as f:
            f.write("001\te1abcA1\tA:1-100\t1abc\tA\n")
            f.write("002\te1abcA2\tA:110-200\t1abc\tA\n")
            f.write("003\te2defB1\tB:5-95\t2def\tB\n")
        
        # Generate definitions
        output_file = tmp_path / "definitions.csv"
        create_domain_definitions_from_cache(str(cache_file), str(output_file))
        
        # Verify output
        assert output_file.exists()
        
        # Read and check content
        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        assert len(rows) == 3
        
        # Check first row
        assert rows[0]['domain_id'] == 'e1abcA1'
        assert rows[0]['pdb_id'] == '1abc'
        assert rows[0]['chain_id'] == 'A'
        assert rows[0]['range'] == '1-100'
        assert rows[0]['length'] == '100'
    
    @pytest.mark.unit
    def test_definitions_sorted_output(self, tmp_path):
        """Test that output is sorted by domain ID"""
        cache_file = tmp_path / "unsorted.txt"
        with open(cache_file, 'w') as f:
            # Write in unsorted order
            f.write("003\tezzzA1\tA:1-50\tzzz\tA\n")
            f.write("001\teaaaA1\tA:1-50\taaa\tA\n")
            f.write("002\temmmA1\tA:1-50\tmmm\tA\n")
        
        output_file = tmp_path / "sorted.csv"
        create_domain_definitions_from_cache(str(cache_file), str(output_file))
        
        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            domain_ids = [row['domain_id'] for row in reader]
        
        # Should be alphabetically sorted
        assert domain_ids == ['eaaaA1', 'emmmA1', 'ezzzA1']


class TestProteinLengthExtraction:
    """Test extracting protein lengths from cache"""
    
    @pytest.mark.unit
    def test_extract_protein_lengths(self, tmp_path):
        """Test extracting maximum protein extent"""
        cache_file = tmp_path / "cache.txt"
        with open(cache_file, 'w') as f:
            # Multiple domains for same chain
            f.write("001\te1abcA1\tA:1-100\t1abc\tA\n")
            f.write("002\te1abcA2\tA:150-250\t1abc\tA\n")  # Max position = 250
            f.write("003\te1abcA3\tA:50-125\t1abc\tA\n")
            f.write("004\te2defB1\tB:10-300\t2def\tB\n")   # Max position = 300
        
        output_file = tmp_path / "protein_lengths.csv"
        extract_protein_lengths_from_cache(str(cache_file), str(output_file))
        
        # Read output
        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = {(r['pdb_id'], r['chain_id']): int(r['length']) for r in reader}
        
        assert len(rows) == 2
        assert rows[('1abc', 'A')] == 250  # Maximum extent
        assert rows[('2def', 'B')] == 300
    
    @pytest.mark.unit
    def test_protein_lengths_discontinuous(self, tmp_path):
        """Test protein length with discontinuous domains"""
        cache_file = tmp_path / "disco.txt"
        with open(cache_file, 'w') as f:
            # Discontinuous domain - should still find max position
            f.write("001\tediscoA1\tA:10-50,A:200-300\tdisco\tA\n")
            f.write("002\tediscoA2\tA:60-150\tdisco\tA\n")
        
        output_file = tmp_path / "disco_lengths.csv"
        extract_protein_lengths_from_cache(str(cache_file), str(output_file))
        
        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        assert len(rows) == 1
        assert rows[0]['pdb_id'] == 'disco'
        assert rows[0]['chain_id'] == 'A'
        assert rows[0]['length'] == '300'  # Max position from discontinuous range


class TestRealCacheFile:
    """Test with real ECOD cache file format"""
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_real_cache_file_format(self):
        """Test parsing real ECOD cache file if available"""
        real_cache = "/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt"
        
        if not Path(real_cache).exists():
            pytest.skip("Real ECOD cache file not available")
        
        # Parse just first 100 entries to test
        entries = {}
        line_count = 0
        
        with open(real_cache, 'r') as f:
            for line in f:
                if line_count >= 100:
                    break
                
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) == 5:
                    cache_id, domain_id, range_spec, pdb_id, chain_id = parts
                    entries[domain_id] = {
                        'cache_id': cache_id,
                        'range_spec': range_spec,
                        'pdb_id': pdb_id,
                        'chain_id': chain_id
                    }
                    line_count += 1
        
        # Verify we parsed some entries
        assert len(entries) > 0
        
        # Check that entries have expected format
        for domain_id, entry in entries.items():
            assert entry['cache_id'].isdigit()  # Should be numeric ID
            assert len(entry['pdb_id']) == 4    # PDB IDs are 4 characters
            assert entry['chain_id']             # Should have chain


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
