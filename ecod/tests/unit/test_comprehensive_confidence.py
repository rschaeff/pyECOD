import pytest
import xml.etree.ElementTree as ET

from ecod.models.pipeline.evidence import Evidence


class TestEvidenceEvalueToConfidence:
    """Test the _evalue_to_confidence method"""

    def test_zero_evalue_returns_max_confidence(self):
        """E-value of 0 should return confidence of 1.0"""
        evidence = Evidence(type="test")
        result = evidence._evalue_to_confidence(0.0)
        assert result == 1.0

    def test_negative_evalue_returns_max_confidence(self):
        """Negative E-values should return confidence of 1.0"""
        evidence = Evidence(type="test")
        result = evidence._evalue_to_confidence(-1.0)
        assert result == 1.0

    def test_very_small_evalue_high_confidence(self):
        """Very small E-values (< 1e-10) should get high confidence with boost"""
        evidence = Evidence(type="test")
        result = evidence._evalue_to_confidence(1e-15)
        assert result > 0.8  # Should be high after fixes
        assert result <= 1.0  # But not exceed 1.0

    def test_good_evalue_reasonable_confidence(self):
        """E-values around 1e-5 should give reasonable confidence"""
        evidence = Evidence(type="test")
        result = evidence._evalue_to_confidence(1e-5)
        assert 0.5 < result < 1.0  # Should be reasonably high

    def test_poor_evalue_low_confidence(self):
        """Large E-values (> 10) should get significant penalty"""
        evidence = Evidence(type="test")
        result = evidence._evalue_to_confidence(100.0)
        assert result < 0.5  # Should be low due to penalty
        assert result >= 0.0

    def test_boundary_evalue_one(self):
        """E-value of 1.0 should give moderate confidence"""
        evidence = Evidence(type="test")
        result = evidence._evalue_to_confidence(1.0)
        assert 0.2 < result < 0.8  # Should be moderate after fixes

    def test_confidence_always_in_valid_range(self):
        """Test various E-values to ensure confidence is always 0.0-1.0"""
        evidence = Evidence(type="test")
        test_evalues = [0, 1e-100, 1e-50, 1e-10, 1e-5, 0.001, 0.1, 1.0, 10.0, 1000.0]

        for evalue in test_evalues:
            result = evidence._evalue_to_confidence(evalue)
            assert 0.0 <= result <= 1.0, f"E-value {evalue} produced invalid confidence {result}"

    def test_confidence_decreases_with_increasing_evalue(self):
        """Confidence should generally decrease as E-value increases"""
        evidence = Evidence(type="test")

        # Test with a sequence of increasing E-values
        evalues = [1e-10, 1e-8, 1e-6, 1e-4, 0.01, 1.0]
        confidences = [evidence._evalue_to_confidence(e) for e in evalues]

        # Check that confidence generally decreases (allowing for small variations due to boosts/penalties)
        for i in range(len(confidences) - 1):
            # Allow some tolerance since there might be boosts for very small evalues
            assert confidences[i] >= confidences[i + 1] * 0.8, \
                f"Confidence should generally decrease: E-value {evalues[i]} -> {confidences[i]}, " \
                f"E-value {evalues[i+1]} -> {confidences[i+1]}"

    def test_extreme_values_handled_gracefully(self):
        """Test that extreme values don't cause crashes"""
        evidence = Evidence(type="test")

        # Test very large numbers
        result = evidence._evalue_to_confidence(1e100)
        assert 0.0 <= result <= 1.0

        # Test very small numbers
        result = evidence._evalue_to_confidence(1e-100)
        assert 0.0 <= result <= 1.0

        # Test float('inf') - this might cause a ValueError or OverflowError
        # which should be handled by the fallback
        result = evidence._evalue_to_confidence(float('inf'))
        assert 0.0 <= result <= 1.0


class TestEvidenceConfidenceCalculation:
    """Test the overall confidence calculation logic"""

    def test_hhsearch_evidence_uses_probability(self):
        """HHSearch evidence should prioritize probability over evalue"""
        evidence = Evidence(
            type="hhsearch",
            probability=95.0,  # 95% on 0-100 scale
            evalue=1e-5
        )
        # Should convert 95% to 0.95 and potentially boost it
        assert evidence.confidence > 0.9

    def test_blast_evidence_uses_evalue(self):
        """BLAST evidence should use evalue for confidence"""
        evidence = Evidence(
            type="domain_blast",
            evalue=1e-10
        )
        # Should calculate confidence from evalue
        assert evidence.confidence >= 0.5

    def test_explicitly_set_confidence_preserved(self):
        """Explicitly set confidence should not be auto-calculated"""
        evidence = Evidence(
            type="blast",
            evalue=1e-5,
            confidence=0.8  # Explicitly set
        )
        assert evidence.confidence == 0.8
        assert evidence._confidence_explicitly_set == True

    def test_auto_calculated_confidence_can_be_recalculated(self):
        """Auto-calculated confidence should be recalculable"""
        evidence = Evidence(
            type="blast",
            evalue=1e-5,
            confidence=None  # Auto-calculate
        )
        original_confidence = evidence.confidence
        assert evidence._confidence_explicitly_set == False

        # Change the underlying data
        evidence.evalue = 1e-10

        # Recalculate
        new_confidence = evidence.recalculate_confidence()
        assert new_confidence != original_confidence
        assert new_confidence > original_confidence  # Better evalue should give higher confidence


class TestHHSearchConfidenceLogic:
    """Test HHSearch-specific confidence calculation logic"""

    def test_hhsearch_prefers_probability_over_evalue(self):
        """HHSearch should prioritize probability over e-value"""
        evidence = Evidence(
            type="hhsearch",
            probability=95.0,  # Excellent
            evalue=1.0,        # Poor
            confidence=None    # Auto-calculate
        )

        # Should be high confidence due to excellent probability
        assert evidence.confidence > 0.9, \
            f"HHSearch with 95% probability should give high confidence, got {evidence.confidence}"

    def test_hhsearch_probability_scale_detection(self):
        """Test both 0-1 and 0-100 probability scales"""

        # Test 0-100 scale (typical HHSearch format)
        evidence_100_scale = Evidence(type="hhsearch", probability=85.0, confidence=None)

        # Test 0-1 scale
        evidence_1_scale = Evidence(type="hhsearch", probability=0.85, confidence=None)

        # Both should give similar confidence (85% is the same as 0.85)
        conf_diff = abs(evidence_100_scale.confidence - evidence_1_scale.confidence)
        assert conf_diff < 0.1, \
            f"85.0 and 0.85 probability should give similar confidence: {evidence_100_scale.confidence} vs {evidence_1_scale.confidence}"

    def test_hhsearch_high_probability_boost(self):
        """Probabilities >= 90% should get confidence boost"""
        high_prob = Evidence(type="hhsearch", probability=95.0, confidence=None)
        medium_prob = Evidence(type="hhsearch", probability=85.0, confidence=None)

        # High probability should get boosted
        assert high_prob.confidence > medium_prob.confidence, \
            "95% probability should give higher confidence than 85%"

        # And should be very high
        assert high_prob.confidence > 0.9, \
            f"95% HHSearch probability should give >0.9 confidence, got {high_prob.confidence}"

    def test_hhsearch_evalue_fallback(self):
        """When probability missing, should fall back to e-value with penalty"""
        with_prob = Evidence(type="hhsearch", probability=80.0, evalue=1e-5, confidence=None)
        evalue_only = Evidence(type="hhsearch", probability=None, evalue=1e-5, confidence=None)

        # E-value only should get penalty
        assert with_prob.confidence > evalue_only.confidence, \
            "HHSearch with probability should beat e-value-only"


class TestBlastConfidenceLogic:
    """Test BLAST-specific confidence calculation logic"""

    def test_blast_uses_evalue_primarily(self):
        """BLAST should use e-value as primary confidence metric"""
        excellent_blast = Evidence(type="domain_blast", evalue=1e-10, confidence=None)
        poor_blast = Evidence(type="domain_blast", evalue=1.0, confidence=None)

        assert excellent_blast.confidence > poor_blast.confidence, \
            f"Better BLAST e-value should give higher confidence: {excellent_blast.confidence} vs {poor_blast.confidence}"

    def test_blast_identity_adjustment(self):
        """BLAST confidence should be adjusted by sequence identity"""
        high_identity = Evidence(
            type="domain_blast",
            evalue=1e-5,
            identity=90.0,  # 90% identity
            confidence=None
        )
        low_identity = Evidence(
            type="domain_blast",
            evalue=1e-5,
            identity=30.0,  # 30% identity
            confidence=None
        )

        # Same e-value, but high identity should win
        assert high_identity.confidence > low_identity.confidence, \
            f"Higher identity should boost confidence: {high_identity.confidence} vs {low_identity.confidence}"

        # Low identity should get significant penalty
        assert low_identity.confidence < 0.5, \
            f"30% identity should get major penalty, got {low_identity.confidence}"

    def test_blast_coverage_adjustment(self):
        """BLAST confidence should consider alignment coverage"""
        good_coverage = Evidence(
            type="domain_blast",
            evalue=1e-5,
            coverage=85.0,  # 85% coverage
            confidence=None
        )
        poor_coverage = Evidence(
            type="domain_blast",
            evalue=1e-5,
            coverage=20.0,  # 20% coverage
            confidence=None
        )

        # Good coverage should beat poor coverage
        assert good_coverage.confidence > poor_coverage.confidence, \
            f"Better coverage should increase confidence: {good_coverage.confidence} vs {poor_coverage.confidence}"

    def test_blast_combined_factors(self):
        """Test BLAST with multiple quality factors"""
        excellent_hit = Evidence(
            type="domain_blast",
            evalue=1e-10,
            identity=85.0,
            coverage=90.0,
            confidence=None
        )

        mediocre_hit = Evidence(
            type="domain_blast",
            evalue=1e-3,
            identity=45.0,
            coverage=30.0,
            confidence=None
        )

        # Excellent hit should significantly outperform mediocre
        assert excellent_hit.confidence > mediocre_hit.confidence + 0.2, \
            f"Excellent hit should significantly outperform mediocre: {excellent_hit.confidence} vs {mediocre_hit.confidence}"


class TestChainBlastConfidenceLogic:
    """Test chain BLAST confidence logic"""

    def test_chain_blast_penalty(self):
        """Chain BLAST should get penalty vs domain BLAST"""
        domain_blast = Evidence(type="domain_blast", evalue=1e-5, confidence=None)
        chain_blast = Evidence(type="chain_blast", evalue=1e-5, confidence=None)

        # Chain BLAST should get mapping uncertainty penalty
        assert domain_blast.confidence > chain_blast.confidence, \
            f"Domain BLAST should beat chain BLAST: {domain_blast.confidence} vs {chain_blast.confidence}"

        # Penalty should be meaningful but not devastating
        penalty = domain_blast.confidence - chain_blast.confidence
        assert 0.05 < penalty < 0.3, \
            f"Chain BLAST penalty should be 5-30%, got {penalty:.2f}"

    def test_chain_blast_hsp_count_penalty(self):
        """Multiple HSPs should increase uncertainty penalty"""
        simple_chain = Evidence(type="chain_blast", evalue=1e-5, hsp_count=1, confidence=None)
        complex_chain = Evidence(type="chain_blast", evalue=1e-5, hsp_count=5, confidence=None)

        # Complex alignments should get additional penalty
        assert simple_chain.confidence > complex_chain.confidence, \
            f"Simple chain alignment should beat complex: {simple_chain.confidence} vs {complex_chain.confidence}"


class TestGenericConfidenceLogic:
    """Test generic confidence calculation"""

    def test_generic_conservative_penalty(self):
        """Generic/unknown types should get conservative penalty"""
        blast_evidence = Evidence(type="domain_blast", evalue=1e-5, confidence=None)
        generic_evidence = Evidence(type="unknown", evalue=1e-5, confidence=None)

        # Generic should be penalized for unknown calibration
        assert blast_evidence.confidence > generic_evidence.confidence, \
            f"Known method should beat generic: {blast_evidence.confidence} vs {generic_evidence.confidence}"

    def test_generic_metric_preference_order(self):
        """Generic should prefer probability > evalue > score"""

        # Test with probability
        with_prob = Evidence(type="unknown", probability=80.0, evalue=1.0, score=50.0, confidence=None)

        # Test with evalue only
        with_evalue = Evidence(type="unknown", probability=None, evalue=1e-5, score=50.0, confidence=None)

        # Test with score only
        with_score = Evidence(type="unknown", probability=None, evalue=None, score=80.0, confidence=None)

        # Should prefer probability over evalue over score
        assert with_prob.confidence > with_evalue.confidence, \
            "Should prefer probability over evalue"
        assert with_evalue.confidence > with_score.confidence, \
            "Should prefer evalue over score"


class TestConfidenceCalculationIntegration:
    """Test overall confidence calculation system"""

    def test_explicit_confidence_preserved(self):
        """Explicitly set confidence should not be auto-calculated"""
        evidence = Evidence(
            type="blast",
            evalue=1e-10,  # Would normally give high confidence
            confidence=0.3  # Explicitly set to low
        )

        assert evidence.confidence == 0.3, "Explicit confidence should be preserved"
        assert evidence._confidence_explicitly_set == True, "Should track explicit setting"

    def test_confidence_recalculation(self):
        """Test confidence recalculation after data changes"""
        evidence = Evidence(type="blast", evalue=1e-3, confidence=None)
        original = evidence.confidence

        # Change underlying data
        evidence.evalue = 1e-10
        new_confidence = evidence.recalculate_confidence()

        # Should recalculate and improve
        assert new_confidence > original, \
            f"Better evalue should improve confidence: {original} -> {new_confidence}"
        assert evidence._confidence_explicitly_set == False, \
            "Recalculation should reset explicit flag"

    def test_confidence_explanation_accuracy(self):
        """Test that confidence explanations match actual calculation"""

        # HHSearch with probability
        hhsearch = Evidence(type="hhsearch", probability=90.0, confidence=None)
        explanation = hhsearch.get_confidence_explanation()
        assert "probability" in explanation.lower(), "Should mention probability"
        assert "90" in explanation, "Should include actual probability value"

        # BLAST with evalue
        blast = Evidence(type="domain_blast", evalue=1e-8, confidence=None)
        explanation = blast.get_confidence_explanation()
        assert "evalue" in explanation.lower() or "e-value" in explanation.lower(), "Should mention evalue"
        assert "1e-08" in explanation or "1.00e-08" in explanation, "Should include actual evalue"


class TestDomainLogicConsistency:
    """Test that confidence calculations match bioinformatics domain knowledge"""

    def test_biologically_meaningful_thresholds(self):
        """Test that confidence thresholds align with field standards"""

        # These are standard significance thresholds in bioinformatics
        test_cases = [
            # (evidence_params, expected_min_confidence, description)
            ({"type": "blast", "evalue": 1e-10}, 0.8, "Highly significant BLAST"),
            ({"type": "blast", "evalue": 1e-5}, 0.6, "Significant BLAST"),
            ({"type": "blast", "evalue": 0.01}, 0.4, "Marginally significant BLAST"),
            ({"type": "blast", "evalue": 1.0}, 0.1, "Non-significant BLAST"),

            ({"type": "hhsearch", "probability": 95.0}, 0.9, "Highly confident HHSearch"),
            ({"type": "hhsearch", "probability": 80.0}, 0.7, "Confident HHSearch"),
            ({"type": "hhsearch", "probability": 60.0}, 0.5, "Moderate HHSearch"),
        ]

        for params, min_conf, description in test_cases:
            evidence = Evidence(confidence=None, **params)
            assert evidence.confidence >= min_conf, \
                f"{description} should have confidence >= {min_conf}, got {evidence.confidence}"

    def test_relative_method_performance(self):
        """Test that method-specific calibrations make sense"""

        # Similar "quality" hits across methods
        excellent_blast = Evidence(type="blast", evalue=1e-10, confidence=None)
        excellent_hhsearch = Evidence(type="hhsearch", probability=95.0, confidence=None)

        # Should both be high confidence
        assert excellent_blast.confidence > 0.8, "Excellent BLAST should be high confidence"
        assert excellent_hhsearch.confidence > 0.8, "Excellent HHSearch should be high confidence"

        # Should be reasonably comparable
        conf_diff = abs(excellent_blast.confidence - excellent_hhsearch.confidence)
        assert conf_diff < 0.2, \
            f"Excellent hits should have similar confidence: {excellent_blast.confidence} vs {excellent_hhsearch.confidence}"
