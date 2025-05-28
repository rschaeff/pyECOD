#!/usr/bin/env python3
"""
Comprehensive tests for Evidence confidence calculation methods.

Tests critical confidence calculation logic including:
- HHSearch confidence calculation
- BLAST confidence calculation  
- Chain BLAST confidence calculation
- E-value to confidence conversion
- XML parsing edge cases
- Round-trip serialization fidelity
"""

import pytest
import math
import xml.etree.ElementTree as ET
from unittest.mock import Mock, patch

from ecod.models.pipeline.evidence import Evidence


class TestConfidenceCalculationCore:
    """Test core confidence calculation methods"""
    
    def test_explicit_confidence_preserved(self):
        """Test that explicitly set confidence is preserved"""
        evidence = Evidence(
            type='hhsearch',
            confidence=0.85  # Explicitly set
        )
        
        assert evidence.confidence == 0.85
        assert evidence._confidence_explicitly_set == True
    
    def test_auto_calculated_confidence(self):
        """Test that confidence is auto-calculated when None"""
        evidence = Evidence(
            type='hhsearch',
            probability=95.0,
            confidence=None  # Auto-calculate
        )
        
        assert evidence.confidence > 0.0
        assert evidence._confidence_explicitly_set == False
    
    def test_set_confidence_explicitly(self):
        """Test setting confidence explicitly after creation"""
        evidence = Evidence(type='blast', evalue=1e-10)
        
        original_confidence = evidence.confidence
        evidence.set_confidence(0.95)
        
        assert evidence.confidence == 0.95
        assert evidence._confidence_explicitly_set == True
        assert evidence.confidence != original_confidence
    
    def test_set_confidence_validation(self):
        """Test confidence validation during explicit setting"""
        evidence = Evidence(type='blast', evalue=1e-10)
        
        # Valid range
        evidence.set_confidence(0.5)
        assert evidence.confidence == 0.5
        
        # Invalid ranges should raise error
        with pytest.raises(ValueError, match="between 0.0 and 1.0"):
            evidence.set_confidence(-0.1)
        
        with pytest.raises(ValueError, match="between 0.0 and 1.0"):
            evidence.set_confidence(1.5)
    
    def test_recalculate_confidence(self):
        """Test forcing confidence recalculation"""
        evidence = Evidence(
            type='hhsearch',
            probability=95.0,
            confidence=0.5  # Explicitly set low
        )
        
        assert evidence.confidence == 0.5
        
        new_confidence = evidence.recalculate_confidence()
        
        assert new_confidence != 0.5
        assert evidence.confidence == new_confidence
        assert evidence._confidence_explicitly_set == False


class TestHHSearchConfidence:
    """Test HHSearch confidence calculation"""
    
    def test_hhsearch_high_probability(self):
        """Test high probability HHSearch hits"""
        evidence = Evidence(
            type='hhsearch',
            probability=99.5,  # Very high
            confidence=None
        )
        
        assert evidence.confidence >= 0.9
        assert evidence.confidence <= 1.0
    
    def test_hhsearch_medium_probability(self):
        """Test medium probability HHSearch hits"""
        evidence = Evidence(
            type='hhsearch',
            probability=85.0,  # Medium
            confidence=None
        )
        
        assert 0.7 <= evidence.confidence <= 0.95
    
    def test_hhsearch_low_probability(self):
        """Test low probability HHSearch hits"""
        evidence = Evidence(
            type='hhsearch',
            probability=50.0,  # Low
            confidence=None
        )
        
        assert 0.2 <= evidence.confidence <= 0.7
    
    def test_hhsearch_probability_scales(self):
        """Test both 0-1 and 0-100 probability scales"""
        # 0-100 scale (typical HHSearch)
        evidence_100 = Evidence(
            type='hhsearch',
            probability=95.0,
            confidence=None
        )
        
        # 0-1 scale
        evidence_1 = Evidence(
            type='hhsearch',
            probability=0.95,
            confidence=None
        )
        
        # Should give similar confidence
        assert abs(evidence_100.confidence - evidence_1.confidence) < 0.1
    
    def test_hhsearch_probability_over_100(self):
        """Test handling probability > 100"""
        evidence = Evidence(
            type='hhsearch',
            probability=150.0,  # Invalid, > 100
            confidence=None
        )
        
        # Should cap at 1.0 confidence
        assert evidence.confidence == 1.0
    
    def test_hhsearch_evalue_fallback(self):
        """Test fallback to evalue when probability missing"""
        evidence = Evidence(
            type='hhsearch',
            evalue=1e-50,  # Very good evalue
            confidence=None
        )
        
        assert evidence.confidence > 0.8
        # Should be slightly lower than probability-based due to penalty
        
        evidence_with_prob = Evidence(
            type='hhsearch',
            probability=99.0,
            confidence=None
        )
        
        assert evidence_with_prob.confidence >= evidence.confidence
    
    def test_hhsearch_score_boost(self):
        """Test score providing confidence boost"""
        evidence_no_score = Evidence(
            type='hhsearch',
            probability=85.0,
            confidence=None
        )
        
        evidence_with_score = Evidence(
            type='hhsearch',
            probability=85.0,
            score=75.0,  # Good score
            confidence=None
        )
        
        # Score should provide small boost
        assert evidence_with_score.confidence >= evidence_no_score.confidence
        assert evidence_with_score.confidence - evidence_no_score.confidence <= 0.1
    
    def test_hhsearch_no_metrics(self):
        """Test HHSearch with no probability or evalue"""
        evidence = Evidence(
            type='hhsearch',
            confidence=None
        )
        
        assert evidence.confidence == 0.0


class TestBlastConfidence:
    """Test BLAST confidence calculation"""
    
    def test_blast_excellent_evalue(self):
        """Test excellent BLAST evalue"""
        evidence = Evidence(
            type='domain_blast',
            evalue=1e-100,  # Excellent
            confidence=None
        )
        
        assert evidence.confidence >= 0.9
    
    def test_blast_good_evalue(self):
        """Test good BLAST evalue"""
        evidence = Evidence(
            type='domain_blast',
            evalue=1e-20,  # Good
            confidence=None
        )
        
        assert 0.7 <= evidence.confidence <= 0.9
    
    def test_blast_moderate_evalue(self):
        """Test moderate BLAST evalue"""
        evidence = Evidence(
            type='domain_blast',
            evalue=1e-5,  # Moderate
            confidence=None
        )
        
        assert 0.5 <= evidence.confidence <= 0.8
    
    def test_blast_poor_evalue(self):
        """Test poor BLAST evalue"""
        evidence = Evidence(
            type='domain_blast',
            evalue=1e-1,  # Poor
            confidence=None
        )
        
        assert evidence.confidence <= 0.5
    
    def test_blast_identity_boost(self):
        """Test identity providing confidence boost"""
        evidence_low_id = Evidence(
            type='domain_blast',
            evalue=1e-10,
            identity=30.0,  # Low identity
            confidence=None
        )
        
        evidence_high_id = Evidence(
            type='domain_blast',
            evalue=1e-10,
            identity=90.0,  # High identity
            confidence=None
        )
        
        # High identity should boost confidence
        assert evidence_high_id.confidence > evidence_low_id.confidence
    
    def test_blast_identity_penalty(self):
        """Test low identity penalizing confidence"""
        evidence = Evidence(
            type='domain_blast',
            evalue=1e-20,  # Good evalue
            identity=20.0,  # Very low identity
            confidence=None
        )
        
        # Should be significantly penalized
        assert evidence.confidence <= 0.6
    
    def test_blast_coverage_effects(self):
        """Test coverage affecting confidence"""
        evidence_poor_cov = Evidence(
            type='domain_blast',
            evalue=1e-10,
            coverage=25.0,  # Poor coverage
            confidence=None
        )
        
        evidence_good_cov = Evidence(
            type='domain_blast',
            evalue=1e-10,
            coverage=85.0,  # Good coverage
            confidence=None
        )
        
        # Good coverage should boost confidence
        assert evidence_good_cov.confidence > evidence_poor_cov.confidence
    
    def test_blast_probability_fallback(self):
        """Test fallback to probability for BLAST variants"""
        evidence = Evidence(
            type='blast',
            probability=0.95,  # Some BLAST variants provide this
            confidence=None
        )
        
        assert evidence.confidence == 0.95
    
    def test_blast_combined_metrics(self):
        """Test BLAST with multiple quality metrics"""
        evidence = Evidence(
            type='domain_blast',
            evalue=1e-30,    # Good evalue
            identity=80.0,   # Good identity
            coverage=90.0,   # Excellent coverage
            confidence=None
        )
        
        # Should get high confidence from combination
        assert evidence.confidence >= 0.8


class TestChainBlastConfidence:
    """Test chain BLAST confidence calculation"""
    
    def test_chain_blast_mapping_penalty(self):
        """Test chain BLAST gets mapping penalty"""
        domain_evidence = Evidence(
            type='domain_blast',
            evalue=1e-20,
            confidence=None
        )
        
        chain_evidence = Evidence(
            type='chain_blast',
            evalue=1e-20,  # Same evalue
            confidence=None
        )
        
        # Chain BLAST should have lower confidence due to mapping uncertainty
        assert chain_evidence.confidence < domain_evidence.confidence
        assert chain_evidence.confidence <= domain_evidence.confidence * 0.85
    
    def test_chain_blast_hsp_penalty(self):
        """Test additional penalty for complex chain alignments"""
        evidence_simple = Evidence(
            type='chain_blast',
            evalue=1e-20,
            hsp_count=1,  # Simple alignment
            confidence=None
        )
        
        evidence_complex = Evidence(
            type='chain_blast',
            evalue=1e-20,
            hsp_count=5,  # Complex multi-HSP
            confidence=None
        )
        
        # Complex alignment should have additional penalty
        assert evidence_complex.confidence < evidence_simple.confidence


class TestEvalueToConfidence:
    """Test E-value to confidence conversion"""
    
    def test_evalue_zero_handling(self):
        """Test handling E-value of 0"""
        evidence = Evidence(
            type='blast',
            evalue=0.0,  # Perfect score
            confidence=None
        )
        
        assert evidence.confidence == 1.0
    
    def test_evalue_negative_handling(self):
        """Test handling negative E-value (invalid)"""
        evidence = Evidence(
            type='blast',
            evalue=-1e-10,  # Invalid
            confidence=None
        )
        
        assert evidence.confidence == 1.0  # Treated as perfect
    
    def test_evalue_boundary_1_0(self):
        """Test E-value boundary at 1.0"""
        evidence = Evidence(
            type='blast',
            evalue=1.0,  # Boundary case
            confidence=None
        )
        
        # Should be in moderate range
        assert 0.2 <= evidence.confidence <= 0.5
    
    def test_evalue_very_large(self):
        """Test very large E-values"""
        evidence = Evidence(
            type='blast',
            evalue=1e10,  # Very poor
            confidence=None
        )
        
        assert evidence.confidence <= 0.2
    
    def test_evalue_progression(self):
        """Test E-value progression gives logical confidence progression"""
        evalues = [1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-2, 1.0, 10.0]
        confidences = []
        
        for evalue in evalues:
            evidence = Evidence(type='blast', evalue=evalue, confidence=None)
            confidences.append(evidence.confidence)
        
        # Confidence should decrease as E-value increases
        for i in range(len(confidences) - 1):
            assert confidences[i] >= confidences[i + 1], \
                f"Confidence not decreasing: {evalues[i]}→{confidences[i]}, {evalues[i+1]}→{confidences[i+1]}"
    
    def test_evalue_exceptional_boost(self):
        """Test boost for exceptionally good E-values"""
        evidence_excellent = Evidence(
            type='blast',
            evalue=1e-200,  # Exceptional
            confidence=None
        )
        
        evidence_very_good = Evidence(
            type='blast',
            evalue=1e-15,  # Very good but not exceptional
            confidence=None
        )
        
        # Exceptional should get boost
        assert evidence_excellent.confidence >= evidence_very_good.confidence
    
    def test_evalue_overflow_handling(self):
        """Test handling of overflow in E-value calculations"""
        evidence = Evidence(
            type='blast',
            evalue=float('inf'),  # Overflow
            confidence=None
        )
        
        # Should handle gracefully
        assert 0.0 <= evidence.confidence <= 1.0


class TestGenericConfidence:
    """Test generic confidence calculation for unknown types"""
    
    def test_generic_probability_priority(self):
        """Test probability has priority over other metrics"""
        evidence = Evidence(
            type='unknown_method',
            probability=0.9,
            evalue=1e-5,  # Also present
            score=50.0,   # Also present
            confidence=None
        )
        
        # Should primarily use probability
        assert evidence.confidence >= 0.8  # Close to probability but with penalty
    
    def test_generic_evalue_fallback(self):
        """Test E-value fallback for unknown types"""
        evidence = Evidence(
            type='unknown_method',
            evalue=1e-20,  # Good evalue
            score=50.0,    # Also present
            confidence=None
        )
        
        # Should use evalue with penalty
        assert 0.4 <= evidence.confidence <= 0.8
    
    def test_generic_score_fallback(self):
        """Test score fallback for unknown types"""
        evidence = Evidence(
            type='unknown_method',
            score=80.0,  # Only metric available
            confidence=None
        )
        
        # Should use score with significant penalty
        assert evidence.confidence <= 0.56  # 80/100 * 0.7 penalty
    
    def test_generic_conservative_penalty(self):
        """Test conservative penalty for unknown types"""
        evidence_known = Evidence(
            type='hhsearch',
            probability=90.0,
            confidence=None
        )
        
        evidence_unknown = Evidence(
            type='unknown_method',
            probability=90.0,
            confidence=None
        )
        
        # Unknown type should have lower confidence
        assert evidence_unknown.confidence < evidence_known.confidence


class TestXMLParsingEdgeCases:
    """Test XML parsing edge cases and error handling"""
    
    def test_blast_xml_malformed_evalue(self):
        """Test BLAST XML with malformed E-value"""
        xml_str = '''
        <hit num="1" domain_id="d1abcA1" evalues="not_a_number">
            <query_range>10-50</query_range>
        </hit>
        '''
        element = ET.fromstring(xml_str)
        
        evidence = Evidence.from_blast_xml(element)
        
        # Should use fallback value
        assert evidence.evalue == 999.0
        assert evidence.confidence <= 0.2  # Poor confidence for fallback
    
    def test_blast_xml_missing_evalue(self):
        """Test BLAST XML with missing E-value"""
        xml_str = '''
        <hit num="1" domain_id="d1abcA1">
            <query_range>10-50</query_range>
        </hit>
        '''
        element = ET.fromstring(xml_str)
        
        evidence = Evidence.from_blast_xml(element)
        
        assert evidence.evalue == 999.0
    
    def test_blast_xml_infinite_evalue(self):
        """Test BLAST XML with infinite E-value"""
        xml_str = '''
        <hit num="1" domain_id="d1abcA1" evalues="inf">
            <query_range>10-50</query_range>
        </hit>
        '''
        element = ET.fromstring(xml_str)
        
        evidence = Evidence.from_blast_xml(element)
        
        assert evidence.evalue == 999.0  # Fallback for non-finite
    
    def test_hhsearch_xml_malformed_probability(self):
        """Test HHSearch XML with malformed probability"""
        xml_str = '''
        <hit hit_id="h1" domain_id="d1abcA1" probability="invalid">
            <query_range>10-50</query_range>
        </hit>
        '''
        element = ET.fromstring(xml_str)
        
        evidence = Evidence.from_hhsearch_xml(element)
        
        assert evidence.probability == 0.0
        assert evidence.confidence == 0.0
    
    def test_hhsearch_xml_nan_values(self):
        """Test HHSearch XML with NaN values"""
        xml_str = '''
        <hit hit_id="h1" domain_id="d1abcA1" probability="nan" score="inf">
            <query_range>10-50</query_range>
        </hit>
        '''
        element = ET.fromstring(xml_str)
        
        evidence = Evidence.from_hhsearch_xml(element)
        
        assert evidence.probability == 0.0
        assert evidence.score == 0.0
    
    def test_xml_missing_ranges(self):
        """Test XML without range elements"""
        xml_str = '''
        <hit hit_id="h1" domain_id="d1abcA1" probability="95.0">
        </hit>
        '''
        element = ET.fromstring(xml_str)
        
        evidence = Evidence.from_hhsearch_xml(element)
        
        assert evidence.query_range == ""
        assert evidence.hit_range == ""
        assert evidence.confidence > 0.9  # Still good confidence from probability
    
    def test_xml_range_in_attributes(self):
        """Test XML with ranges as attributes instead of elements"""
        xml_str = '''
        <hit hit_id="h1" domain_id="d1abcA1" probability="95.0" 
             query_range="10-50" hit_range="5-45">
        </hit>
        '''
        element = ET.fromstring(xml_str)
        
        evidence = Evidence.from_hhsearch_xml(element)
        
        assert evidence.query_range == "10-50"
        assert evidence.hit_range == "5-45"


class TestXMLRoundTrip:
    """Test XML serialization round-trip fidelity"""
    
    def test_evidence_xml_round_trip_basic(self):
        """Test basic evidence XML round-trip"""
        original = Evidence(
            type='hhsearch',
            source_id='test_source',
            domain_id='d1abcA1',
            query_range='10-50',
            hit_range='5-45',
            probability=95.0,
            evalue=1e-30,
            score=75.0,
            confidence=0.95,
            t_group='1.1.1.1',
            h_group='1.1.1'
        )
        
        # Convert to XML and back
        xml_element = original.to_xml()
        reconstructed = Evidence.from_xml(xml_element)
        
        # Check core fields preserved
        assert reconstructed.type == original.type
        assert reconstructed.source_id == original.source_id
        assert reconstructed.domain_id == original.domain_id
        assert reconstructed.query_range == original.query_range
        assert reconstructed.hit_range == original.hit_range
        assert reconstructed.probability == original.probability
        assert reconstructed.evalue == original.evalue
        assert reconstructed.score == original.score
        assert abs(reconstructed.confidence - original.confidence) < 1e-9
        assert reconstructed.t_group == original.t_group
        assert reconstructed.h_group == original.h_group
    
    def test_evidence_xml_round_trip_none_values(self):
        """Test round-trip with None values"""
        original = Evidence(
            type='blast',
            domain_id='d1abcA1',
            evalue=1e-20,
            probability=None,  # None value
            score=None,       # None value
            confidence=0.85
        )
        
        xml_element = original.to_xml()
        reconstructed = Evidence.from_xml(xml_element)
        
        assert reconstructed.type == original.type
        assert reconstructed.domain_id == original.domain_id
        assert reconstructed.evalue == original.evalue
        assert reconstructed.probability is None
        assert reconstructed.score is None
        assert abs(reconstructed.confidence - original.confidence) < 1e-9
    
    def test_evidence_xml_confidence_state_preservation(self):
        """Test confidence state preservation in round-trip"""
        # Explicitly set confidence
        explicit = Evidence(
            type='hhsearch',
            probability=90.0,
            confidence=0.75  # Explicitly different from auto-calculated
        )
        
        xml_element = explicit.to_xml()
        reconstructed = Evidence.from_xml(xml_element)
        
        # Should preserve explicit confidence exactly
        assert abs(reconstructed.confidence - 0.75) < 1e-9
        assert reconstructed._confidence_explicitly_set == True
    
    def test_evidence_xml_auto_calculated_confidence(self):
        """Test auto-calculated confidence in round-trip"""
        # Auto-calculated confidence
        auto = Evidence(
            type='hhsearch',
            probability=90.0,
            confidence=None  # Auto-calculate
        )
        
        original_confidence = auto.confidence
        
        xml_element = auto.to_xml()
        reconstructed = Evidence.from_xml(xml_element)
        
        # Should preserve calculated confidence
        assert abs(reconstructed.confidence - original_confidence) < 1e-9
        assert reconstructed._confidence_explicitly_set == True  # From XML


class TestConfidenceExplanation:
    """Test confidence explanation functionality"""
    
    def test_confidence_explanation_explicit(self):
        """Test explanation for explicitly set confidence"""
        evidence = Evidence(
            type='hhsearch',
            probability=90.0,
            confidence=0.85  # Explicit
        )
        
        explanation = evidence.get_confidence_explanation()
        
        assert 'explicitly set' in explanation.lower()
        assert '0.85' in explanation
    
    def test_confidence_explanation_hhsearch(self):
        """Test explanation for HHSearch confidence"""
        evidence = Evidence(
            type='hhsearch',
            probability=95.0,
            confidence=None  # Auto-calculate
        )
        
        explanation = evidence.get_confidence_explanation()
        
        assert 'hhsearch probability' in explanation.lower()
        assert '95.0' in explanation
    
    def test_confidence_explanation_blast(self):
        """Test explanation for BLAST confidence"""
        evidence = Evidence(
            type='domain_blast',
            evalue=1e-30,
            identity=85.0,
            coverage=90.0,
            confidence=None
        )
        
        explanation = evidence.get_confidence_explanation()
        
        assert 'blast e-value' in explanation.lower()
        assert '1e-30' in explanation
        assert 'identity: 85.0%' in explanation
        assert 'coverage: 90.0%' in explanation
    
    def test_confidence_explanation_no_confidence(self):
        """Test explanation when no confidence calculated"""
        evidence = Evidence(type='unknown')
        evidence.confidence = None
        
        explanation = evidence.get_confidence_explanation()
        
        assert 'no confidence calculated' in explanation.lower()


class TestConfidenceEdgeCases:
    """Test edge cases in confidence calculation"""
    
    def test_evidence_no_metrics(self):
        """Test evidence with no quality metrics"""
        evidence = Evidence(
            type='blast',
            domain_id='d1abcA1',
            confidence=None
        )
        
        assert evidence.confidence == 0.0
    
    def test_evidence_all_invalid_metrics(self):
        """Test evidence with all invalid metrics"""
        evidence = Evidence(
            type='hhsearch',
            probability=float('nan'),
            evalue=float('inf'),
            score=float('-inf'),
            confidence=None
        )
        
        assert evidence.confidence == 0.0
    
    def test_evidence_mixed_valid_invalid(self):
        """Test evidence with mix of valid and invalid metrics"""
        evidence = Evidence(
            type='hhsearch',
            probability=float('nan'),  # Invalid
            evalue=1e-20,             # Valid
            confidence=None
        )
        
        # Should use valid evalue with fallback
        assert evidence.confidence > 0.0
    
    def test_confidence_calculation_exception_handling(self):
        """Test exception handling in confidence calculation"""
        # Mock the calculation to raise an exception
        with patch.object(Evidence, '_calculate_hhsearch_confidence', 
                         side_effect=Exception("Test error")):
            evidence = Evidence(
                type='hhsearch',
                probability=95.0,
                confidence=None
            )
            
            # Should default to 0.0 on exception
            assert evidence.confidence == 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
