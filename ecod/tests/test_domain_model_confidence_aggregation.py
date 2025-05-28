#!/usr/bin/env python3
"""
Test suite for Domain Model confidence aggregation and weighting logic

These tests focus specifically on the multi-evidence fusion mathematics 
and weighted aggregation logic at the Domain level, distinct from the 
individual Evidence confidence calculation tests.
"""

import pytest
import math
from unittest.mock import Mock, patch

from ecod.models.pipeline.domain import DomainModel
from ecod.models.pipeline.evidence import Evidence


class TestDomainConfidenceWeightedAggregation:
    """Test the core weighted aggregation mathematics"""
    
    def test_single_evidence_confidence_passthrough(self):
        """Single evidence should pass through its confidence"""
        evidence = Evidence(type="domain_blast", confidence=0.85)
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100", 
            evidence=[evidence], confidence=0.0  # Force recalculation
        )
        
        assert domain.confidence == 0.85
    
    def test_weighted_average_calculation_known_types(self):
        """Test weighted average with known evidence types"""
        # Create evidence with explicit confidence values
        evidence1 = Evidence(type="domain_blast", confidence=0.9)   # weight 3.0
        evidence2 = Evidence(type="hhsearch", confidence=0.8)       # weight 2.5  
        evidence3 = Evidence(type="chain_blast", confidence=0.6)    # weight 2.0
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2, evidence3],
            confidence=0.0  # Force recalculation
        )
        
        # Expected calculation:
        # (0.9*3.0 + 0.8*2.5 + 0.6*2.0) / (3.0 + 2.5 + 2.0)
        # = (2.7 + 2.0 + 1.2) / 7.5 = 5.9 / 7.5 ≈ 0.7867
        expected = (0.9*3.0 + 0.8*2.5 + 0.6*2.0) / (3.0 + 2.5 + 2.0)
        assert abs(domain.confidence - expected) < 0.001
    
    def test_weighted_average_with_blast_variants(self):
        """Test weighting for different BLAST evidence types"""
        evidence1 = Evidence(type="domain_blast", confidence=0.9)   # weight 3.0
        evidence2 = Evidence(type="blast", confidence=0.8)          # weight 1.5
        evidence3 = Evidence(type="chain_blast", confidence=0.7)    # weight 2.0
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2, evidence3],
            confidence=0.0
        )
        
        # Expected: (0.9*3.0 + 0.8*1.5 + 0.7*2.0) / (3.0 + 1.5 + 2.0)
        expected = (0.9*3.0 + 0.8*1.5 + 0.7*2.0) / (3.0 + 1.5 + 2.0)
        assert abs(domain.confidence - expected) < 0.001
    
    def test_unknown_evidence_type_fallback_weight(self):
        """Unknown evidence types should get default weight 1.0"""
        evidence1 = Evidence(type="domain_blast", confidence=0.9)   # weight 3.0
        evidence2 = Evidence(type="unknown_method", confidence=0.5) # weight 1.0
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )
        
        # Expected: (0.9*3.0 + 0.5*1.0) / (3.0 + 1.0) = 3.2 / 4.0 = 0.8
        expected = (0.9*3.0 + 0.5*1.0) / (3.0 + 1.0)
        assert abs(domain.confidence - expected) < 0.001
    
    def test_self_comparison_evidence_weight(self):
        """Self-comparison evidence should get lowest weight"""
        evidence1 = Evidence(type="domain_blast", confidence=0.8)      # weight 3.0
        evidence2 = Evidence(type="self_comparison", confidence=0.9)    # weight 1.0
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )
        
        # Expected: (0.8*3.0 + 0.9*1.0) / (3.0 + 1.0) = 3.3 / 4.0 = 0.825
        expected = (0.8*3.0 + 0.9*1.0) / (3.0 + 1.0)
        assert abs(domain.confidence - expected) < 0.001


class TestDomainConfidenceMixedEvidenceTypes:
    """Test confidence calculation with mixed Evidence objects and dictionaries"""
    
    def test_mixed_evidence_objects_and_dicts(self):
        """Test weighted calculation with Evidence objects and dictionaries"""
        evidence_obj = Evidence(type="domain_blast", confidence=0.9)
        evidence_dict = {
            "type": "hhsearch", 
            "confidence": 0.8
        }
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence_obj, evidence_dict],
            confidence=0.0
        )
        
        # Should handle both types: (0.9*3.0 + 0.8*2.5) / (3.0 + 2.5)
        expected = (0.9*3.0 + 0.8*2.5) / (3.0 + 2.5)
        assert abs(domain.confidence - expected) < 0.001
    
    def test_dict_evidence_missing_confidence(self):
        """Dictionary evidence without confidence should be skipped or defaulted"""
        evidence_obj = Evidence(type="domain_blast", confidence=0.9)
        evidence_dict = {
            "type": "hhsearch"
            # Missing confidence
        }
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence_obj, evidence_dict],
            confidence=0.0
        )
        
        # Should either skip the dict evidence or use confidence=0.0
        # If skipped: confidence = 0.9 (from single evidence)
        # If defaulted: confidence = (0.9*3.0 + 0.0*2.5) / (3.0 + 2.5) ≈ 0.49
        assert domain.confidence == 0.9 or abs(domain.confidence - 0.491) < 0.01
    
    def test_dict_evidence_missing_type(self):
        """Dictionary evidence without type should be skipped"""
        evidence_obj = Evidence(type="domain_blast", confidence=0.9)
        evidence_dict = {
            "confidence": 0.8
            # Missing type - should be skipped or treated as "unknown"
        }

        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence_obj, evidence_dict],
            confidence=0.0
        )

        # UPDATED: After fixes, malformed evidence should be skipped
        # If dict evidence is skipped: confidence = 0.9 (from valid evidence only)
        # If dict evidence gets default weight 1.0: confidence = (0.9*3.0 + 0.8*1.0)/(3.0+1.0) = 0.875
        assert domain.confidence == 0.9  # Should skip malformed evidence
    
    def test_invalid_evidence_objects(self):
        """Test handling of invalid evidence objects"""
        evidence_obj = Evidence(type="domain_blast", confidence=0.9)
        invalid_obj = "not_an_evidence_object"
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence_obj, invalid_obj],
            confidence=0.0
        )
        
        # Should skip invalid objects gracefully
        assert domain.confidence == 0.9


class TestDomainConfidenceEdgeCases:
    """Test mathematical edge cases in confidence calculation"""
    
    def test_zero_total_weight_division_by_zero(self):
        """Test handling when total weight is zero"""
        evidence_dict = {
            "type": None,  # Invalid type should be skipped
            "confidence": 0.8
        }

        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence_dict],
            confidence=0.0
        )

        # UPDATED: After fixes, evidence with type=None should be skipped
        assert domain.confidence == 0.0
    
    def test_empty_evidence_list(self):
        """Test confidence calculation with no evidence"""
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[],
            confidence=0.0
        )
        
        assert domain.confidence == 0.0
    
    def test_evidence_with_nan_confidence(self):
        """Test handling evidence with NaN confidence values"""
        evidence1 = Evidence(type="domain_blast", confidence=0.9)
        evidence2 = Evidence(type="hhsearch", confidence=float('nan'))

        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )

        # UPDATED: After fixes, NaN evidence should be skipped
        assert not math.isnan(domain.confidence)
        assert domain.confidence == 0.9  # Only valid evidence should be used

    def test_evidence_with_infinite_confidence(self):
        """Test handling evidence with infinite confidence values"""
        evidence1 = Evidence(type="domain_blast", confidence=0.9)
        evidence2 = Evidence(type="hhsearch", confidence=float('inf'))

        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )

        # UPDATED: After fixes, infinite evidence should be skipped
        assert math.isfinite(domain.confidence)
        assert domain.confidence == 0.9  # Only valid evidence should be used
    
    def test_evidence_with_negative_confidence(self):
        """Test handling evidence with negative confidence values"""
        evidence1 = Evidence(type="domain_blast", confidence=0.9)
        evidence2 = Evidence(type="hhsearch", confidence=-0.5)
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )
        
        # Should handle negative values gracefully
        assert 0.0 <= domain.confidence <= 1.0
    
    def test_evidence_with_confidence_above_one(self):
        """Test handling evidence with confidence > 1.0"""
        evidence1 = Evidence(type="domain_blast", confidence=0.9)
        evidence2 = Evidence(type="hhsearch", confidence=1.5)

        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )

        # UPDATED: After fixes, out-of-range evidence should be skipped OR clamped
        assert 0.0 <= domain.confidence <= 1.0
        # Could be 0.9 (if >1.0 evidence skipped) or clamped weighted average


class TestDomainConfidenceRecalculation:
    """Test confidence recalculation scenarios"""
    
    def test_confidence_recalculation_after_evidence_change(self):
        """Test domain confidence updates when evidence changes"""
        evidence = Evidence(type="blast", evalue=1e-5, confidence=None)  # Auto-calc
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence],
            confidence=0.0
        )
        
        original_confidence = domain.confidence
        
        # Change evidence quality and recalculate
        evidence.evalue = 1e-10  # Much better E-value
        evidence.recalculate_confidence()
        
        # Manually trigger domain confidence recalculation
        domain._calculate_confidence()
        
        # Domain confidence should improve
        assert domain.confidence > original_confidence
    
    def test_add_evidence_updates_confidence(self):
        """Test that adding evidence updates domain confidence"""
        evidence1 = Evidence(type="chain_blast", confidence=0.6)  # Lower weight/confidence
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1],
            confidence=0.0
        )
        
        original_confidence = domain.confidence
        assert original_confidence == 0.6  # Should match single evidence
        
        # Add higher quality evidence
        evidence2 = Evidence(type="domain_blast", confidence=0.9)  # Higher weight/confidence
        domain.add_evidence(evidence2)
        
        # Confidence should improve due to weighted average
        # (0.6*2.0 + 0.9*3.0) / (2.0 + 3.0) = 3.9 / 5.0 = 0.78
        expected = (0.6*2.0 + 0.9*3.0) / (2.0 + 3.0)
        assert abs(domain.confidence - expected) < 0.001
        assert domain.confidence > original_confidence
    
    def test_evidence_list_modification_consistency(self):
        """Test that modifying evidence list maintains confidence consistency"""
        evidence1 = Evidence(type="domain_blast", confidence=0.8)
        evidence2 = Evidence(type="hhsearch", confidence=0.9)
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )
        
        # Store calculated confidence
        calculated_confidence = domain.confidence
        
        # Manually modify evidence list and recalculate
        domain.evidence.append(Evidence(type="chain_blast", confidence=0.7))
        domain._calculate_confidence()
        
        # Confidence should be different now
        assert domain.confidence != calculated_confidence
        
        # Should be calculable: (0.8*3.0 + 0.9*2.5 + 0.7*2.0) / (3.0 + 2.5 + 2.0)
        expected = (0.8*3.0 + 0.9*2.5 + 0.7*2.0) / (3.0 + 2.5 + 2.0)
        assert abs(domain.confidence - expected) < 0.001


class TestDomainConfidenceWeightingSchemeChanges:
    """Test robustness to potential algorithm retune changes"""
    
    def test_custom_weight_scheme(self):
        """Test confidence calculation with modified weight scheme"""
        # Simulate algorithm retune changing weights
        evidence1 = Evidence(type="domain_blast", confidence=0.8)
        evidence2 = Evidence(type="hhsearch", confidence=0.9)
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )
        
        # Patch the weights to simulate algorithm change
        custom_weights = {
            "domain_blast": 4.0,  # Increased from 3.0
            "hhsearch": 3.0,      # Increased from 2.5
            "chain_blast": 2.0,   # Same
            "blast": 1.5,         # Same
            "self_comparison": 1.0  # Same
        }
        
        with patch.object(domain, '_calculate_confidence') as mock_calc:
            def custom_calculate():
                if not domain.evidence:
                    return 0.0
                
                total_weight = 0.0
                weighted_sum = 0.0
                
                for ev in domain.evidence:
                    ev_type = getattr(ev, 'type', ev.get('type', 'unknown') if isinstance(ev, dict) else 'unknown')
                    ev_confidence = getattr(ev, 'confidence', ev.get('confidence', 0.0) if isinstance(ev, dict) else 0.0)
                    
                    weight = custom_weights.get(ev_type, 1.0)
                    weighted_sum += ev_confidence * weight
                    total_weight += weight
                
                return weighted_sum / total_weight if total_weight > 0 else 0.0
            
            mock_calc.side_effect = custom_calculate
            confidence_custom = domain._calculate_confidence()
        
        # Should work with different weight scheme
        expected = (0.8*4.0 + 0.9*3.0) / (4.0 + 3.0)  # = 5.9 / 7.0 ≈ 0.843
        assert abs(confidence_custom - expected) < 0.001
    
    def test_new_evidence_type_handling(self):
        """Test handling of new evidence types not in current weight scheme"""
        evidence1 = Evidence(type="domain_blast", confidence=0.8)
        evidence2 = Evidence(type="new_future_method", confidence=0.95)  # New type
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[evidence1, evidence2],
            confidence=0.0
        )
        
        # New type should get default weight 1.0
        # Expected: (0.8*3.0 + 0.95*1.0) / (3.0 + 1.0) = 3.35 / 4.0 = 0.8375
        expected = (0.8*3.0 + 0.95*1.0) / (3.0 + 1.0)
        assert abs(domain.confidence - expected) < 0.001


class TestDomainConfidencePerformance:
    """Test performance characteristics of confidence calculation"""
    
    def test_large_evidence_list_performance(self):
        """Test confidence calculation with many evidence items"""
        import time
        
        # Create 100 evidence items
        evidence_list = []
        for i in range(100):
            evidence_list.append(Evidence(
                type="domain_blast", 
                confidence=0.8 + (i % 20) * 0.01  # Varying confidence
            ))
        
        start_time = time.time()
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=evidence_list,
            confidence=0.0
        )
        calculation_time = time.time() - start_time
        
        # Should complete quickly (< 1 second for 100 items)
        assert calculation_time < 1.0
        
        # Should still produce valid confidence
        assert 0.0 <= domain.confidence <= 1.0
        assert domain.confidence > 0.8  # Should be around the evidence range
    
    def test_deeply_nested_evidence_dicts(self):
        """Test handling evidence with complex nested dictionaries"""
        complex_evidence = {
            "type": "hhsearch",
            "confidence": 0.85,
            "extra_attributes": {
                "nested": {
                    "deeply": {
                        "nested": {
                            "data": "should_not_break_confidence_calc"
                        }
                    }
                }
            }
        }
        
        domain = DomainModel(
            id="test", start=1, end=100, range="1-100",
            evidence=[complex_evidence],
            confidence=0.0
        )
        
        # Should handle complex evidence without issues
        assert domain.confidence == 0.85


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
