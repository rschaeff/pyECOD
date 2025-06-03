#!/usr/bin/env python3
"""
Domain decomposer tests for mini_pyecod

Tests domain definition loading and blacklist functionality.
"""

import pytest
from pathlib import Path
import csv

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.decomposer import (
    load_domain_definitions,
    load_reference_blacklist,
    DomainReference,
    BlacklistEntry
)
from ecod.core.sequence_range import SequenceRange


class TestDomainDefinitionLoading:
    """Test loading domain definitions from CSV"""
    
    @pytest.mark.unit
    def test_load_basic_domain_definitions(self, tmp_path):
        """Test loading basic domain definitions"""
        csv_file = tmp_path / "domain_defs.csv"
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['domain_id', 'pdb_id', 'chain_id', 'range', 'length', 't_group', 'h_group'])
            writer.writerow(['e2ia4A1', '2ia4', 'A', '110-209', '100', '1234.1.1', '1234.1'])
            writer.writerow(['e2ia4A2', '2ia4', 'A', '210-308', '99', '1234.1.2', '1234.1'])
            writer.writerow(['e6dgvA1', '6dgv', 'A', '1-238', '238', '5678.2.1', '5678.2'])
        
        definitions = load_domain_definitions(str(csv_file))
        
        # Check we got the right structure
        assert len(definitions) == 2  # 2 unique chains
        assert ('2ia4', 'A') in definitions
        assert ('6dgv', 'A') in definitions
        
        # Check 2ia4 domains
        ia4_domains = definitions[('2ia4', 'A')]
        assert len(ia4_domains) == 2
        assert ia4_domains[0].domain_id == 'e2ia4A1'
        assert ia4_domains[0].pdb_id == '2ia4'
        assert ia4_domains[0].chain_id == 'A'
        assert str(ia4_domains[0].range) == '110-209'
        assert ia4_domains[0].length == 100
        assert ia4_domains[0].t_group == '1234.1.1'
        assert ia4_domains[0].h_group == '1234.1'
        
        # Check sorting by position
        assert ia4_domains[0].range.segments[0].start < ia4_domains[1].range.segments[0].start
    
    @pytest.mark.unit
    def test_load_with_invalid_ranges(self, tmp_path):
        """Test handling of invalid range formats"""
        csv_file = tmp_path / "bad_ranges.csv"
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['domain_id', 'pdb_id', 'chain_id', 'range', 'length'])
            writer.writerow(['good1', 'test', 'A', '1-100', '100'])
            writer.writerow(['bad1', 'test', 'A', 'invalid-range', '50'])  # Bad range
            writer.writerow(['good2', 'test', 'A', '150-200', '51'])
            writer.writerow(['bad2', 'test', 'B', '', '0'])  # Empty range
        
        definitions = load_domain_definitions(str(csv_file), verbose=True)
        
        # Should load only valid domains
        assert ('test', 'A') in definitions
        assert len(definitions[('test', 'A')]) == 2  # Only good1 and good2
        
        # Bad domain shouldn't be there
        domain_ids = [d.domain_id for d in definitions[('test', 'A')]]
        assert 'good1' in domain_ids
        assert 'good2' in domain_ids
        assert 'bad1' not in domain_ids
    
    @pytest.mark.unit
    def test_load_with_invalid_lengths(self, tmp_path):
        """Test handling of invalid lengths"""
        csv_file = tmp_path / "bad_lengths.csv"
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['domain_id', 'pdb_id', 'chain_id', 'range', 'length'])
            writer.writerow(['good', 'test', 'A', '1-100', '100'])
            writer.writerow(['zero', 'test', 'A', '101-150', '0'])     # Zero length
            writer.writerow(['negative', 'test', 'A', '151-200', '-50'])  # Negative
            writer.writerow(['missing', 'test', 'A', '201-250', ''])   # Missing
        
        definitions = load_domain_definitions(str(csv_file))
        
        # Should only have the good one
        assert ('test', 'A') in definitions
        assert len(definitions[('test', 'A')]) == 1
        assert definitions[('test', 'A')][0].domain_id == 'good'
    
    @pytest.mark.unit
    def test_case_normalization(self, tmp_path):
        """Test that PDB IDs are normalized to lowercase"""
        csv_file = tmp_path / "mixed_case.csv"
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['domain_id', 'pdb_id', 'chain_id', 'range', 'length'])
            writer.writerow(['e1ABCA1', '1ABC', 'A', '1-100', '100'])
            writer.writerow(['e2DefB1', '2DeF', 'B', '1-50', '50'])
        
        definitions = load_domain_definitions(str(csv_file))
        
        # Should normalize PDB IDs to lowercase
        assert ('1abc', 'A') in definitions
        assert ('2def', 'B') in definitions
        assert ('1ABC', 'A') not in definitions


class TestBlacklistLoading:
    """Test reference blacklist functionality"""
    
    @pytest.mark.unit
    def test_load_basic_blacklist(self, tmp_path):
        """Test loading a basic blacklist"""
        blacklist_file = tmp_path / "blacklist.csv"
        
        with open(blacklist_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['pdb_id', 'chain_id', 'reason', 'date_added', 'added_by'])
            writer.writerow(['1bad', 'A', 'incomplete structure', '2025-01-15', 'admin'])
            writer.writerow(['2bad', 'B', 'non-standard residues', '2025-01-16', 'curator'])
        
        blacklist = load_reference_blacklist(str(blacklist_file))
        
        assert len(blacklist) == 2
        assert ('1bad', 'A') in blacklist
        assert ('2bad', 'B') in blacklist
    
    @pytest.mark.unit
    def test_load_blacklist_with_verbose(self, tmp_path):
        """Test verbose blacklist loading"""
        blacklist_file = tmp_path / "blacklist_verbose.csv"
        
        with open(blacklist_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['pdb_id', 'chain_id', 'reason', 'date_added', 'added_by'])
            writer.writerow(['3xyz', 'C', 'test reason', '2025-01-17', 'tester'])
        
        # Verbose should print details but still work
        blacklist = load_reference_blacklist(str(blacklist_file), verbose=True)
        
        assert len(blacklist) == 1
        assert ('3xyz', 'C') in blacklist
    
    @pytest.mark.unit
    def test_empty_blacklist(self, tmp_path):
        """Test loading empty blacklist"""
        blacklist_file = tmp_path / "empty_blacklist.csv"
        
        with open(blacklist_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['pdb_id', 'chain_id', 'reason', 'date_added', 'added_by'])
            # No data rows
        
        blacklist = load_reference_blacklist(str(blacklist_file))
        assert len(blacklist) == 0
    
    @pytest.mark.unit
    def test_missing_blacklist_file(self):
        """Test handling of missing blacklist file"""
        blacklist = load_reference_blacklist("/nonexistent/blacklist.csv")
        assert len(blacklist) == 0
    
    @pytest.mark.unit
    def test_blacklist_with_empty_rows(self, tmp_path):
        """Test blacklist with empty rows"""
        blacklist_file = tmp_path / "blacklist_empty_rows.csv"
        
        with open(blacklist_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['pdb_id', 'chain_id', 'reason', 'date_added', 'added_by'])
            writer.writerow(['good', 'A', 'valid entry', '2025-01-18', 'admin'])
            writer.writerow(['', '', '', '', ''])  # Empty row
            writer.writerow(['also_good', 'B', 'another valid', '2025-01-19', 'admin'])
        
        blacklist = load_reference_blacklist(str(blacklist_file))
        
        # Should skip empty rows
        assert len(blacklist) == 2
        assert ('good', 'A') in blacklist
        assert ('also_good', 'B') in blacklist


class TestBlacklistIntegration:
    """Test blacklist integration with domain loading"""
    
    @pytest.mark.unit
    def test_domain_loading_with_blacklist(self, tmp_path):
        """Test that blacklisted chains are excluded from domain definitions"""
        # Create domain definitions
        defs_file = tmp_path / "domains.csv"
        with open(defs_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['domain_id', 'pdb_id', 'chain_id', 'range', 'length'])
            writer.writerow(['eGoodA1', 'good', 'A', '1-100', '100'])
            writer.writerow(['eBadB1', 'bad', 'B', '1-50', '50'])
            writer.writerow(['eGoodC1', 'good', 'C', '1-75', '75'])
        
        # Create blacklist
        blacklist_file = tmp_path / "blacklist.csv"
        with open(blacklist_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['pdb_id', 'chain_id', 'reason', 'date_added', 'added_by'])
            writer.writerow(['bad', 'B', 'problematic structure', '2025-01-20', 'admin'])
        
        # Load with blacklist
        definitions = load_domain_definitions(
            str(defs_file),
            blacklist_path=str(blacklist_file)
        )
        
        # Should have domains from 'good' but not 'bad'
        assert ('good', 'A') in definitions
        assert ('good', 'C') in definitions
        assert ('bad', 'B') not in definitions
        
        # Verify counts
        assert len(definitions) == 2  # Only good chains
        assert len(definitions[('good', 'A')]) == 1
        assert len(definitions[('good', 'C')]) == 1
    
    @pytest.mark.unit
    def test_case_insensitive_blacklist(self, tmp_path):
        """Test that blacklist matching is case-insensitive for PDB IDs"""
        defs_file = tmp_path / "domains_mixed.csv"
        with open(defs_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['domain_id', 'pdb_id', 'chain_id', 'range', 'length'])
            writer.writerow(['e1ABCA1', '1ABC', 'A', '1-100', '100'])  # Uppercase in file
        
        blacklist_file = tmp_path / "blacklist_lower.csv"
        with open(blacklist_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['pdb_id', 'chain_id', 'reason', 'date_added', 'added_by'])
            writer.writerow(['1abc', 'A', 'test', '2025-01-21', 'admin'])  # Lowercase in blacklist
        
        definitions = load_domain_definitions(
            str(defs_file),
            blacklist_path=str(blacklist_file)
        )
        
        # Should be blacklisted despite case difference
        assert ('1abc', 'A') not in definitions


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
