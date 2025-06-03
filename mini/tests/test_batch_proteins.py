#!/usr/bin/env python3
"""
Test cases for proteins discovered from batch analysis
"""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from mini.batch_test_proteins import BATCH_TEST_PROTEINS


class TestBatchProteins:
    """Test cases for proteins from batch"""

    @pytest.fixture
    def run_mini_algorithm(self, stable_batch_dir, real_reference_data, blast_alignments):
        """Run mini algorithm on a protein"""
        def _run(protein_id):
            from mini.parser import parse_domain_summary
            from mini.partitioner import partition_domains
            import os
            
            xml_path = os.path.join(stable_batch_dir, "domains", 
                                   f"{protein_id}.develop291.domain_summary.xml")
            
            if not os.path.exists(xml_path):
                return {"success": False, "error": "Domain summary not found"}
            
            try:
                evidence = parse_domain_summary(
                    xml_path,
                    reference_lengths=real_reference_data.get("domain_lengths", {}),
                    protein_lengths=real_reference_data.get("protein_lengths", {}),
                    blast_alignments=blast_alignments,
                    require_reference_lengths=False
                )
                
                if not evidence:
                    return {"success": False, "error": "No evidence found"}
                
                max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
                sequence_length = int(max_pos * 1.1)
                
                domains = partition_domains(
                    evidence,
                    sequence_length=sequence_length,
                    domain_definitions=real_reference_data.get("domain_definitions", {})
                )
                
                return {
                    "success": True,
                    "domains": domains,
                    "sequence_length": sequence_length,
                    "evidence_count": len(evidence)
                }
            except Exception as e:
                return {"success": False, "error": str(e)}
        
        return _run

    @pytest.mark.integration
    @pytest.mark.parametrize("protein_id", ['8oni_L', '8p6i_L', '8p2e_B', '8oz3_B', '8p12_L'])
    def test_chain_blast_multi_proteins(self, protein_id, run_mini_algorithm):
        """Test chain blast multi proteins"""
        info = BATCH_TEST_PROTEINS[protein_id]
        result = run_mini_algorithm(protein_id)
        
        # Basic success check
        assert result["success"], f"Failed: {result.get('error')}"
        assert len(result["domains"]) > 0, "Should find at least one domain"
        
        # Category-specific checks

    @pytest.mark.integration
    @pytest.mark.parametrize("protein_id", ['8oqj_A', '8oqh_A'])
    def test_single_domain_small_proteins(self, protein_id, run_mini_algorithm):
        """Test single domain small proteins"""
        info = BATCH_TEST_PROTEINS[protein_id]
        result = run_mini_algorithm(protein_id)
        
        # Basic success check
        assert result["success"], f"Failed: {result.get('error')}"
        assert len(result["domains"]) > 0, "Should find at least one domain"
        
        # Category-specific checks
        # Single domain proteins should have 1 domain
        assert len(result["domains"]) == 1, \
            f"Expected 1 domain, found {len(result['domains'])}"

    @pytest.mark.integration
    @pytest.mark.parametrize("protein_id", ['8oyx_A', '8oyw_A'])
    def test_single_domain_medium_proteins(self, protein_id, run_mini_algorithm):
        """Test single domain medium proteins"""
        info = BATCH_TEST_PROTEINS[protein_id]
        result = run_mini_algorithm(protein_id)
        
        # Basic success check
        assert result["success"], f"Failed: {result.get('error')}"
        assert len(result["domains"]) > 0, "Should find at least one domain"
        
        # Category-specific checks
        # Single domain proteins should have 1 domain
        assert len(result["domains"]) == 1, \
            f"Expected 1 domain, found {len(result['domains'])}"

    @pytest.mark.integration
    @pytest.mark.parametrize("protein_id", ['8p49_A', '8p8o_H'])
    def test_multi_domain_clear_proteins(self, protein_id, run_mini_algorithm):
        """Test multi domain clear proteins"""
        info = BATCH_TEST_PROTEINS[protein_id]
        result = run_mini_algorithm(protein_id)
        
        # Basic success check
        assert result["success"], f"Failed: {result.get('error')}"
        assert len(result["domains"]) > 0, "Should find at least one domain"
        
        # Category-specific checks
        # Multi-domain proteins should have multiple domains
        assert len(result["domains"]) >= 2, \
            f"Expected multiple domains, found {len(result['domains'])}"

    @pytest.mark.integration
    @pytest.mark.parametrize("protein_id", ['8oyu_A'])
    def test_large_complex_proteins(self, protein_id, run_mini_algorithm):
        """Test large complex proteins"""
        info = BATCH_TEST_PROTEINS[protein_id]
        result = run_mini_algorithm(protein_id)
        
        # Basic success check
        assert result["success"], f"Failed: {result.get('error')}"
        assert len(result["domains"]) > 0, "Should find at least one domain"
        
        # Category-specific checks

    @pytest.mark.integration
    @pytest.mark.parametrize("protein_id", ['8olg_A'])
    def test_minimal_evidence_proteins(self, protein_id, run_mini_algorithm):
        """Test minimal evidence proteins"""
        info = BATCH_TEST_PROTEINS[protein_id]
        result = run_mini_algorithm(protein_id)
        
        # Basic success check
        assert result["success"], f"Failed: {result.get('error')}"
        assert len(result["domains"]) > 0, "Should find at least one domain"
        
        # Category-specific checks
