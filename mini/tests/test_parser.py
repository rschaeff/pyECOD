#!/usr/bin/env python3
"""
Domain summary XML parser tests for mini_pyecod

Tests the parsing of domain summary XML files and reference data loading.
"""

import pytest
from pathlib import Path
import csv

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.core.parser import (
    parse_domain_summary, 
    load_reference_lengths, 
    load_protein_lengths,
    get_evidence_summary
)
from mini.core.models import Evidence
from ecod.core.sequence_range import SequenceRange


class TestDomainSummaryParsing:
    """Test domain summary XML parsing"""
    
    @pytest.mark.unit
    def test_parse_valid_domain_summary(self, tmp_path):
        """Test parsing a valid domain summary XML"""
        xml_content = """<?xml version="1.0"?>
<blast_summ_doc>
  <blast_summ pdb="8ovp" chain="A"/>
  <chain_blast_run program="blastp">
    <hits>
      <hit num="1" pdb_id="6dgv" chain_id="A" hsp_count="1" evalues="1e-50">
        <query_reg>252-494</query_reg>
        <hit_reg>1-243</hit_reg>
      </hit>
    </hits>
  </chain_blast_run>
  <blast_run program="blastp">
    <hits>
      <hit domain_id="e6dgvA1" pdb_id="6dgv" chain_id="A" evalues="1e-30">
        <query_reg>260-480</query_reg>
        <hit_reg>8-228</hit_reg>
      </hit>
    </hits>
  </blast_run>
  <hh_run program="hhsearch">
    <hits>
      <hit hit_id="6dgv_A" domain_id="e6dgvA1" num="1" probability="99.5" evalue="1e-20" score="150.5">
        <query_reg>255-490</query_reg>
        <hit_reg>5-240</hit_reg>
      </hit>
    </hits>
  </hh_run>
</blast_summ_doc>"""
        
        xml_file = tmp_path / "test_summary.xml"
        xml_file.write_text(xml_content)
        
        # Create mock reference data
        reference_lengths = {"e6dgvA1": 238, "6dgv": 238}
        protein_lengths = {("6dgv", "A"): 238}
        
        # Parse
        evidence = parse_domain_summary(
            str(xml_file),
            reference_lengths=reference_lengths,
            protein_lengths=protein_lengths,
            require_reference_lengths=False
        )
        
        # Verify we got all evidence types
        assert len(evidence) == 3
        
        # Check chain BLAST
        chain_blast = [e for e in evidence if e.type == "chain_blast"]
        assert len(chain_blast) == 1
        assert chain_blast[0].source_pdb == "6dgv"
        assert str(chain_blast[0].query_range) == "252-494"
        assert chain_blast[0].evalue == 1e-50
        
        # Check domain BLAST
        domain_blast = [e for e in evidence if e.type == "domain_blast"]
        assert len(domain_blast) == 1
        assert domain_blast[0].domain_id == "e6dgvA1"
        assert str(domain_blast[0].query_range) == "260-480"
        
        # Check HHSearch
        hhsearch = [e for e in evidence if e.type == "hhsearch"]
        assert len(hhsearch) == 1
        assert hhsearch[0].confidence == 0.995  # 99.5% converted to 0-1
        assert str(hhsearch[0].query_range) == "255-490"
    
    @pytest.mark.unit
    def test_parse_empty_domain_summary(self, tmp_path):
        """Test parsing empty domain summary"""
        xml_content = """<?xml version="1.0"?>
<blast_summ_doc>
  <blast_summ pdb="test" chain="A"/>
  <chain_blast_run program="blastp">
    <hits></hits>
  </chain_blast_run>
  <blast_run program="blastp">
    <hits></hits>
  </blast_run>
  <hh_run program="hhsearch">
    <hits></hits>
  </hh_run>
</blast_summ_doc>"""
        
        xml_file = tmp_path / "empty_summary.xml"
        xml_file.write_text(xml_content)
        
        evidence = parse_domain_summary(str(xml_file))
        assert len(evidence) == 0
    
    @pytest.mark.unit
    def test_parse_missing_query_regions(self, tmp_path):
        """Test handling of hits without query regions"""
        xml_content = """<?xml version="1.0"?>
<blast_summ_doc>
  <blast_summ pdb="test" chain="A"/>
  <chain_blast_run program="blastp">
    <hits>
      <hit num="1" pdb_id="test" chain_id="B" evalues="0.01">
        <!-- Missing query_reg -->
      </hit>
    </hits>
  </chain_blast_run>
</blast_summ_doc>"""
        
        xml_file = tmp_path / "missing_regions.xml"
        xml_file.write_text(xml_content)
        
        evidence = parse_domain_summary(str(xml_file))
        # Should skip hits without query regions
        assert len(evidence) == 0
    
    @pytest.mark.unit
    def test_reference_length_filtering(self, tmp_path):
        """Test that evidence without reference lengths is filtered when required"""
        xml_content = """<?xml version="1.0"?>
<blast_summ_doc>
  <blast_summ pdb="test" chain="A"/>
  <domain_blast_run>
    <hits>
      <hit domain_id="eTestA1" evalues="1e-10">
        <query_reg>1-100</query_reg>
      </hit>
      <hit domain_id="eTestB1" evalues="1e-20">
        <query_reg>150-250</query_reg>
      </hit>
    </hits>
  </domain_blast_run>
</blast_summ_doc>"""
        
        xml_file = tmp_path / "test.xml"
        xml_file.write_text(xml_content)
        
        # Only provide reference for one domain
        reference_lengths = {"eTestA1": 100}
        
        # With require_reference_lengths=True (default)
        evidence = parse_domain_summary(
            str(xml_file),
            reference_lengths=reference_lengths,
            require_reference_lengths=True
        )
        
        # Should only have the one with reference length
        assert len(evidence) == 1
        assert evidence[0].domain_id == "eTestA1"
        assert evidence[0].reference_length == 100
        
        # With require_reference_lengths=False
        evidence = parse_domain_summary(
            str(xml_file),
            reference_lengths=reference_lengths,
            require_reference_lengths=False
        )
        
        # Should have both
        assert len(evidence) == 2
    
    @pytest.mark.unit
    def test_confidence_calculation(self, tmp_path):
        """Test confidence score calculation from e-values"""
        xml_content = """<?xml version="1.0"?>
<blast_summ_doc>
  <blast_summ pdb="test" chain="A"/>
  <blast_run program="blastp">
    <hits>
      <hit domain_id="d1" evalues="1e-50">
        <query_reg>1-100</query_reg>
      </hit>
      <hit domain_id="d2" evalues="1e-10">
        <query_reg>101-200</query_reg>
      </hit>
      <hit domain_id="d3" evalues="1e-5">
        <query_reg>201-300</query_reg>
      </hit>
      <hit domain_id="d4" evalues="0.01">
        <query_reg>301-400</query_reg>
      </hit>
    </hits>
  </blast_run>
</blast_summ_doc>"""
        
        xml_file = tmp_path / "confidence_test.xml"
        xml_file.write_text(xml_content)
        
        # Don't require reference lengths for this test
        evidence = parse_domain_summary(str(xml_file), require_reference_lengths=False)
        
        # Check confidence scores based on e-values
        assert len(evidence) == 4
        
        # Better e-value = higher confidence
        conf_by_domain = {e.domain_id: e.confidence for e in evidence}
        assert conf_by_domain["d1"] > conf_by_domain["d2"]
        assert conf_by_domain["d2"] > conf_by_domain["d3"]
        assert conf_by_domain["d3"] > conf_by_domain["d4"]
        
        # Check specific thresholds
        assert conf_by_domain["d1"] >= 0.9  # e-value < 1e-10
        assert conf_by_domain["d2"] >= 0.9  # e-value = 1e-10
        assert conf_by_domain["d3"] >= 0.7  # e-value < 1e-5
        assert conf_by_domain["d4"] >= 0.5  # e-value >= 0.001


class TestReferenceLengthLoading:
    """Test loading reference length files"""
    
    @pytest.mark.unit
    def test_load_domain_lengths(self, tmp_path):
        """Test loading domain reference lengths from CSV"""
        csv_file = tmp_path / "domain_lengths.csv"
        
        # Create test CSV
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['domain_id', 'length'])  # Header
            writer.writerow(['e6dgvA1', '238'])
            writer.writerow(['e2ia4A1', '98'])
            writer.writerow(['e2ia4A2', '156'])
        
        lengths = load_reference_lengths(str(csv_file))
        
        assert len(lengths) == 3
        assert lengths['e6dgvA1'] == 238
        assert lengths['e2ia4A1'] == 98
        assert lengths['e2ia4A2'] == 156
    
    @pytest.mark.unit
    def test_load_domain_lengths_no_header(self, tmp_path):
        """Test loading without header row"""
        csv_file = tmp_path / "lengths_no_header.csv"
        
        # No header, straight to data
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['test1', '100'])
            writer.writerow(['test2', '200'])
        
        lengths = load_reference_lengths(str(csv_file))
        
        assert len(lengths) == 2
        assert lengths['test1'] == 100
        assert lengths['test2'] == 200
    
    @pytest.mark.unit
    def test_load_protein_lengths(self, tmp_path):
        """Test loading protein lengths from CSV"""
        csv_file = tmp_path / "protein_lengths.csv"
        
        # Format: pdb_id,chain_id,length
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['pdb_id', 'chain_id', 'length'])
            writer.writerow(['6dgv', 'A', '238'])
            writer.writerow(['2ia4', 'A', '508'])
            writer.writerow(['8ovp', 'A', '569'])
        
        lengths = load_protein_lengths(str(csv_file))
        
        assert len(lengths) == 3
        assert lengths[('6dgv', 'A')] == 238
        assert lengths[('2ia4', 'A')] == 508
        assert lengths[('8ovp', 'A')] == 569
    
    @pytest.mark.unit
    def test_load_protein_lengths_combined_format(self, tmp_path):
        """Test loading protein lengths with pdb_chain format"""
        csv_file = tmp_path / "protein_lengths_alt.csv"
        
        # Format: pdb_chain,length
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['protein_id', 'length'])
            writer.writerow(['6dgv_A', '238'])
            writer.writerow(['2ia4_B', '508'])
        
        lengths = load_protein_lengths(str(csv_file))
        
        assert len(lengths) == 2
        assert lengths[('6dgv', 'A')] == 238
        assert lengths[('2ia4', 'B')] == 508
    
    @pytest.mark.unit
    def test_load_missing_file(self):
        """Test handling of missing files"""
        lengths = load_reference_lengths("/nonexistent/file.csv")
        assert lengths == {}
        
        protein_lengths = load_protein_lengths("/nonexistent/file.csv")
        assert protein_lengths == {}


class TestEvidenceSummary:
    """Test evidence summary statistics"""
    
    @pytest.mark.unit
    def test_evidence_summary_empty(self):
        """Test summary of empty evidence"""
        summary = get_evidence_summary([])
        
        assert summary['total'] == 0
        assert summary['by_type'] == {}
        assert summary['high_confidence'] == 0
        assert summary['unique_families'] == 0
    
    @pytest.mark.unit
    def test_evidence_summary_mixed(self):
        """Test summary of mixed evidence types"""
        evidence = [
            Evidence(
                type="chain_blast",
                source_pdb="pdb1",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.95,
                evalue=1e-50
            ),
            Evidence(
                type="chain_blast",
                source_pdb="pdb2",
                query_range=SequenceRange.parse("150-250"),
                confidence=0.6,
                evalue=1e-5
            ),
            Evidence(
                type="domain_blast",
                source_pdb="pdb1",
                query_range=SequenceRange.parse("10-90"),
                confidence=0.8,
                evalue=1e-20
            ),
            Evidence(
                type="hhsearch",
                source_pdb="pdb3",
                query_range=SequenceRange.parse("300-400"),
                confidence=0.85,
                evalue=1e-15
            )
        ]
        
        summary = get_evidence_summary(evidence)
        
        assert summary['total'] == 4
        assert summary['by_type']['chain_blast'] == 2
        assert summary['by_type']['domain_blast'] == 1
        assert summary['by_type']['hhsearch'] == 1
        assert summary['high_confidence'] == 3  # conf > 0.7 or evalue < 1e-10
        assert summary['unique_families'] == 3  # pdb1, pdb2, pdb3


class TestIntegrationWithRealData:
    """Integration tests with real domain summary files"""
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_real_domain_summary_8ovp(self, domain_summary_path, real_reference_data):
        """Test parsing real 8ovp_A domain summary"""
        evidence = parse_domain_summary(
            domain_summary_path,
            reference_lengths=real_reference_data['domain_lengths'],
            protein_lengths=real_reference_data['protein_lengths']
        )
        
        # Should have evidence
        assert len(evidence) > 0
        
        # Check evidence types
        summary = get_evidence_summary(evidence)
        assert summary['total'] > 0
        
        # Should have multiple evidence types
        assert len(summary['by_type']) >= 2
        
        # Check for expected families
        families = {e.source_pdb for e in evidence if e.source_pdb}
        expected_families = {"6dgv", "2ia4"}  # Known hits for 8ovp_A
        assert len(families.intersection(expected_families)) > 0


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
