#!/usr/bin/env python3
"""
BLAST XML parser tests for mini_pyecod

Tests the BLAST alignment parsing functionality.
"""

import pytest
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from mini.blast_parser import parse_blast_xml, load_chain_blast_alignments, BlastAlignment


class TestBlastXMLParsing:
    """Test BLAST XML parsing functionality"""
    
    @pytest.mark.unit
    def test_parse_valid_blast_xml(self, tmp_path):
        """Test parsing a valid BLAST XML file"""
        # Create a minimal valid BLAST XML
        blast_xml = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_def>2ia4 A</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_query-from>10</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>5</Hsp_hit-from>
              <Hsp_hit-to>95</Hsp_hit-to>
              <Hsp_qseq>ABCDEFGHIJ</Hsp_qseq>
              <Hsp_hseq>ABCDEFGHIJ</Hsp_hseq>
              <Hsp_evalue>1e-50</Hsp_evalue>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
        
        # Write to temp file
        xml_file = tmp_path / "test_blast.xml"
        xml_file.write_text(blast_xml)
        
        # Parse
        alignments = parse_blast_xml(str(xml_file))
        
        # Verify
        assert len(alignments) == 1
        assert ("2ia4", "A") in alignments
        
        alignment = alignments[("2ia4", "A")]
        assert alignment.query_start == 10
        assert alignment.query_end == 100
        assert alignment.hit_start == 5
        assert alignment.hit_end == 95
        assert alignment.query_seq == "ABCDEFGHIJ"
        assert alignment.hit_seq == "ABCDEFGHIJ"
        assert alignment.evalue == 1e-50
    
    @pytest.mark.unit
    def test_parse_empty_blast_xml(self, tmp_path):
        """Test parsing BLAST XML with no hits"""
        blast_xml = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
        
        xml_file = tmp_path / "empty_blast.xml"
        xml_file.write_text(blast_xml)
        
        alignments = parse_blast_xml(str(xml_file))
        assert len(alignments) == 0
    
    @pytest.mark.unit
    def test_parse_malformed_xml(self, tmp_path):
        """Test handling of malformed XML"""
        malformed_xml = """<?xml version="1.0"?>
<BlastOutput>
  <Unclosed tag
</BlastOutput>"""
        
        xml_file = tmp_path / "malformed.xml"
        xml_file.write_text(malformed_xml)
        
        # Should return empty dict on parse error
        alignments = parse_blast_xml(str(xml_file))
        assert len(alignments) == 0
    
    @pytest.mark.unit
    def test_parse_missing_file(self):
        """Test handling of missing file"""
        alignments = parse_blast_xml("/nonexistent/file.xml")
        assert len(alignments) == 0
    
    @pytest.mark.unit
    def test_multiple_hits(self, tmp_path):
        """Test parsing multiple hits"""
        blast_xml = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_def>2ia4 A</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>50</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>50</Hsp_hit-to>
              <Hsp_qseq>AAAAA</Hsp_qseq>
              <Hsp_hseq>AAAAA</Hsp_hseq>
              <Hsp_evalue>1e-20</Hsp_evalue>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_def>6dgv B</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_query-from>60</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>50</Hsp_hit-to>
              <Hsp_qseq>BBBBB</Hsp_qseq>
              <Hsp_hseq>BBBBB</Hsp_hseq>
              <Hsp_evalue>1e-30</Hsp_evalue>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
        
        xml_file = tmp_path / "multi_blast.xml"
        xml_file.write_text(blast_xml)
        
        alignments = parse_blast_xml(str(xml_file))
        
        assert len(alignments) == 2
        assert ("2ia4", "A") in alignments
        assert ("6dgv", "B") in alignments
        
        # Verify different e-values
        assert alignments[("2ia4", "A")].evalue == 1e-20
        assert alignments[("6dgv", "B")].evalue == 1e-30


class TestBlastAlignmentLoader:
    """Test the BLAST alignment loader"""
    
    @pytest.mark.unit
    def test_load_chain_blast_alignments(self, tmp_path):
        """Test loading alignments for a specific query"""
        # Create test BLAST directory and file
        blast_dir = tmp_path / "blast" / "chain"
        blast_dir.mkdir(parents=True)
        
        # Create BLAST XML for 8ovp_A
        blast_xml = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_def>2ia4 B</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_query-from>2</Hsp_query-from>
              <Hsp_query-to>248</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>247</Hsp_hit-to>
              <Hsp_qseq>QUERYSEQ</Hsp_qseq>
              <Hsp_hseq>HITSEQ</Hsp_hseq>
              <Hsp_evalue>1e-100</Hsp_evalue>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
        
        blast_file = blast_dir / "8ovp_A.develop291.xml"
        blast_file.write_text(blast_xml)
        
        # Load alignments
        alignments = load_chain_blast_alignments(str(blast_dir), "8ovp", "A")
        
        assert len(alignments) == 1
        assert ("2ia4", "B") in alignments
        
        alignment = alignments[("2ia4", "B")]
        assert alignment.query_start == 2
        assert alignment.query_end == 248
    
    @pytest.mark.unit
    def test_load_missing_blast_file(self, tmp_path):
        """Test loading when BLAST file doesn't exist"""
        blast_dir = tmp_path / "blast" / "chain"
        blast_dir.mkdir(parents=True)
        
        # Try to load non-existent file
        alignments = load_chain_blast_alignments(str(blast_dir), "9xyz", "Z")
        
        assert len(alignments) == 0
    
    @pytest.mark.unit
    def test_load_uppercase_filename(self, tmp_path):
        """Test loading with uppercase PDB ID in filename"""
        blast_dir = tmp_path / "blast" / "chain"
        blast_dir.mkdir(parents=True)
        
        # Create file with uppercase PDB ID
        blast_xml = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_def>test A</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>10</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>10</Hsp_hit-to>
              <Hsp_qseq>TEST</Hsp_qseq>
              <Hsp_hseq>TEST</Hsp_hseq>
              <Hsp_evalue>0.001</Hsp_evalue>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
        
        # Save with uppercase
        blast_file = blast_dir / "1ABC_A.develop291.xml"
        blast_file.write_text(blast_xml)
        
        # Load with lowercase query
        alignments = load_chain_blast_alignments(str(blast_dir), "1abc", "A")
        
        # Should find the uppercase file
        assert len(alignments) == 1


class TestBlastAlignmentIntegration:
    """Integration tests with real BLAST data"""
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_real_blast_file_8ovp(self, stable_batch_dir):
        """Test with real 8ovp_A BLAST file if available"""
        blast_dir = Path(stable_batch_dir) / "blast" / "chain"
        
        if not blast_dir.exists():
            pytest.skip("BLAST directory not available")
        
        alignments = load_chain_blast_alignments(str(blast_dir), "8ovp", "A")
        
        # Should have alignments
        assert len(alignments) > 0
        
        # Check for expected hit to 2ia4
        ia4_hits = [k for k in alignments.keys() if k[0] == "2ia4"]
        assert len(ia4_hits) > 0, "Expected to find 2ia4 hits for 8ovp_A"
        
        # Verify alignment has required fields
        for (pdb, chain), alignment in alignments.items():
            assert isinstance(alignment.query_start, int)
            assert isinstance(alignment.query_end, int)
            assert isinstance(alignment.hit_start, int)
            assert isinstance(alignment.hit_end, int)
            assert isinstance(alignment.query_seq, str)
            assert isinstance(alignment.hit_seq, str)
            assert isinstance(alignment.evalue, float)
            assert alignment.query_end > alignment.query_start
            assert alignment.hit_end > alignment.hit_start


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
