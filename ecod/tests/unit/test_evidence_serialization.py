import pytest
import xml.etree.ElementTree as ET

from ecod.models.pipeline.evidence import Evidence


class TestXMLSerialization:
    """Test XML serialization/deserialization"""

    def test_to_xml_basic_structure(self):
        """Test that to_xml creates proper XML structure"""
        evidence = Evidence(
            type="blast",
            source_id="test_hit",
            domain_id="PF00001",
            evalue=1e-5,
            confidence=0.8
        )

        xml_element = evidence.to_xml()

        # Check basic structure
        assert xml_element.tag == "evidence"
        assert xml_element.get("type") == "blast"
        assert xml_element.get("source_id") == "test_hit"
        assert xml_element.get("domain_id") == "PF00001"
        assert xml_element.get("evalue") == "1e-05"
        assert xml_element.get("confidence") == "0.8000000000"

    def test_xml_round_trip_preserves_data(self):
        """Test that XML round-trip preserves all data"""
        original = Evidence(
            type="hhsearch",
            source_id="2ABC_A",
            domain_id="PF00123",
            query_range="10-50",
            hit_range="5-45",
            probability=85.5,
            evalue=1e-8,
            score=42.3,
            confidence=0.9,
            identity=78.5,
            coverage=65.2,
            hsp_count=3
        )

        # Convert to XML and back
        xml_element = original.to_xml()
        reconstructed = Evidence.from_xml(xml_element)

        # Check all core fields are preserved
        assert reconstructed.type == original.type
        assert reconstructed.source_id == original.source_id
        assert reconstructed.domain_id == original.domain_id
        assert reconstructed.query_range == original.query_range
        assert reconstructed.hit_range == original.hit_range
        assert reconstructed.probability == original.probability
        assert reconstructed.evalue == original.evalue
        assert reconstructed.score == original.score
        assert abs(reconstructed.confidence - original.confidence) < 0.001
        assert reconstructed.identity == original.identity
        assert reconstructed.coverage == original.coverage
        assert reconstructed.hsp_count == original.hsp_count

    def test_xml_handles_none_values(self):
        """Test XML serialization with None values"""
        evidence = Evidence(
            type="blast",
            source_id="test",
            evalue=None,  # None values should be handled gracefully
            probability=None,
            score=None
        )

        xml_element = evidence.to_xml()
        reconstructed = Evidence.from_xml(xml_element)

        # None values should remain None
        assert reconstructed.evalue is None
        assert reconstructed.probability is None
        assert reconstructed.score is None

    def test_xml_numeric_precision(self):
        """Test that numeric precision is preserved in XML"""
        evidence = Evidence(
            type="blast",
            evalue=1.23456789e-10,  # High precision E-value
            probability=95.123456,   # High precision probability
            confidence=0.987654321   # High precision confidence
        )

        xml_element = evidence.to_xml()
        reconstructed = Evidence.from_xml(xml_element)

        # Check precision is reasonably preserved
        assert abs(reconstructed.evalue - evidence.evalue) < 1e-15
        assert abs(reconstructed.probability - evidence.probability) < 1e-6
        assert abs(reconstructed.confidence - evidence.confidence) < 1e-6


class TestBLASTXMLParsing:
    """Test BLAST XML parsing specifically"""

    def create_blast_xml_element(self, **kwargs):
        """Helper to create BLAST XML test elements"""
        defaults = {
            "num": "1",
            "domain_id": "PF00001",
            "pdb_id": "1ABC",
            "chain_id": "A",
            "evalues": "1e-5",
            "hsp_count": "2"
        }
        defaults.update(kwargs)

        element = ET.Element("hit")
        for key, value in defaults.items():
            if value is not None:
                element.set(key, str(value))

        # Add query and hit regions
        if "query_range" in kwargs:
            query_reg = ET.SubElement(element, "query_reg")
            query_reg.text = kwargs["query_range"]

        if "hit_range" in kwargs:
            hit_reg = ET.SubElement(element, "hit_reg")
            hit_reg.text = kwargs["hit_range"]

        return element

    def test_blast_xml_basic_parsing(self):
        """Test basic BLAST XML parsing"""
        xml_element = self.create_blast_xml_element(
            domain_id="PF00123",
            evalues="1e-8",
            hsp_count="3",
            query_range="10-50",
            hit_range="15-55"
        )

        evidence = Evidence.from_blast_xml(xml_element)

        assert evidence.type == "domain_blast"
        assert evidence.domain_id == "PF00123"
        assert evidence.evalue == 1e-8
        assert evidence.hsp_count == 3
        assert evidence.query_range == "10-50"
        assert evidence.hit_range == "15-55"

    def test_blast_xml_comma_separated_evalues(self):
        """Test BLAST XML with comma-separated E-values"""
        xml_element = self.create_blast_xml_element(
            evalues="1e-8,2e-7,3e-6"  # Multiple E-values
        )

        evidence = Evidence.from_blast_xml(xml_element)

        # Should take the first (best) E-value
        assert evidence.evalue == 1e-8

    def test_blast_xml_malformed_evalues(self):
        """Test BLAST XML with malformed E-values"""
        xml_element = self.create_blast_xml_element(
            evalues="not_a_number"
        )

        evidence = Evidence.from_blast_xml(xml_element)

        # Should fall back to default E-value
        assert evidence.evalue == 999.0

    def test_blast_xml_missing_attributes(self):
        """Test BLAST XML with missing attributes"""
        # Create minimal XML element
        element = ET.Element("hit")
        element.set("num", "1")
        # Missing most attributes

        evidence = Evidence.from_blast_xml(element)

        # Should handle missing attributes gracefully
        assert evidence.type == "domain_blast"
        assert evidence.source_id == "1"  # Falls back to num
        assert evidence.domain_id == ""
        assert evidence.evalue == 999.0
        assert evidence.hsp_count == 0

    def test_blast_xml_extra_attributes(self):
        """Test BLAST XML preserves extra attributes"""
        xml_element = self.create_blast_xml_element(
            pdb_id="2DEF",
            chain_id="B",
            discontinuous="true"
        )

        evidence = Evidence.from_blast_xml(xml_element)

        # Extra attributes should be preserved
        assert evidence.extra_attributes["pdb_id"] == "2DEF"
        assert evidence.extra_attributes["chain_id"] == "B"
        assert evidence.extra_attributes["discontinuous"] == True


class TestHHSearchXMLParsing:
    """Test HHSearch XML parsing specifically"""

    def create_hhsearch_xml_element(self, **kwargs):
        """Helper to create HHSearch XML test elements"""
        defaults = {
            "hit_id": "2ABC_A",
            "domain_id": "PF00001",
            "probability": "85.5",
            "evalue": "1e-6",
            "score": "42.3",
            "num": "1"
        }
        defaults.update(kwargs)

        element = ET.Element("hit")
        for key, value in defaults.items():
            if value is not None:
                element.set(key, str(value))

        # Add query and hit regions
        if "query_range" in kwargs:
            query_reg = ET.SubElement(element, "query_reg")
            query_reg.text = kwargs["query_range"]

        if "hit_range" in kwargs:
            hit_reg = ET.SubElement(element, "hit_reg")
            hit_reg.text = kwargs["hit_range"]

        return element

    def test_hhsearch_xml_basic_parsing(self):
        """Test basic HHSearch XML parsing"""
        xml_element = self.create_hhsearch_xml_element(
            domain_id="PF00456",
            probability="92.1",
            evalue="1e-12",
            score="55.7",
            query_range="5-45",
            hit_range="10-50"
        )

        evidence = Evidence.from_hhsearch_xml(xml_element)

        assert evidence.type == "hhsearch"
        assert evidence.domain_id == "PF00456"
        assert evidence.probability == 92.1
        assert evidence.evalue == 1e-12
        assert evidence.score == 55.7
        assert evidence.query_range == "5-45"
        assert evidence.hit_range == "10-50"

    def test_hhsearch_xml_malformed_numbers(self):
        """Test HHSearch XML with malformed numeric values"""
        xml_element = self.create_hhsearch_xml_element(
            probability="not_a_number",
            evalue="also_not_a_number",
            score="still_not_a_number"
        )

        evidence = Evidence.from_hhsearch_xml(xml_element)

        # Should fall back to default values
        assert evidence.probability == 0.0
        assert evidence.evalue == 999.0
        assert evidence.score == 0.0

    def test_hhsearch_xml_missing_regions(self):
        """Test HHSearch XML with missing query/hit regions"""
        xml_element = self.create_hhsearch_xml_element()
        # No query_reg or hit_reg child elements

        evidence = Evidence.from_hhsearch_xml(xml_element)

        assert evidence.query_range == ""
        assert evidence.hit_range == ""


class TestDictionarySerialization:
    """Test dictionary serialization/deserialization"""

    def test_to_dict_completeness(self):
        """Test that to_dict includes all relevant fields"""
        evidence = Evidence(
            type="blast",
            source_id="test_hit",
            domain_id="PF00001",
            query_range="10-50",
            hit_range="5-45",
            evalue=1e-5,
            probability=85.0,
            score=42.0,
            confidence=0.8,
            identity=75.5,
            coverage=80.2,
            hsp_count=2,
            t_group="Hydrolase",
            h_group="Enzyme",
            extra_attributes={"custom_field": "custom_value"}
        )

        data_dict = evidence.to_dict()

        # Check all fields are present
        assert data_dict["type"] == "blast"
        assert data_dict["source_id"] == "test_hit"
        assert data_dict["domain_id"] == "PF00001"
        assert data_dict["query_range"] == "10-50"
        assert data_dict["hit_range"] == "5-45"
        assert data_dict["evalue"] == 1e-5
        assert data_dict["probability"] == 85.0
        assert data_dict["score"] == 42.0
        assert data_dict["confidence"] == 0.8
        assert data_dict["identity"] == 75.5
        assert data_dict["coverage"] == 80.2
        assert data_dict["hsp_count"] == 2
        assert data_dict["t_group"] == "Hydrolase"
        assert data_dict["h_group"] == "Enzyme"
        assert data_dict["custom_field"] == "custom_value"

    def test_dict_round_trip_preserves_data(self):
        """Test dictionary round-trip preserves all data"""
        original = Evidence(
            type="hhsearch",
            source_id="2ABC_A",
            evalue=1e-10,
            probability=95.7,
            confidence=0.95,
            extra_attributes={"pdb_id": "2ABC", "resolution": 1.5}
        )

        # Convert to dict and back
        data_dict = original.to_dict()
        reconstructed = Evidence.from_dict(data_dict)

        # Check key fields preserved
        assert reconstructed.type == original.type
        assert reconstructed.source_id == original.source_id
        assert reconstructed.evalue == original.evalue
        assert reconstructed.probability == original.probability
        assert reconstructed.confidence == original.confidence
        assert reconstructed.extra_attributes == original.extra_attributes

    def test_from_dict_handles_missing_fields(self):
        """Test from_dict handles missing fields gracefully"""
        minimal_dict = {
            "type": "blast",
            "source_id": "test"
            # Missing most fields
        }

        evidence = Evidence.from_dict(minimal_dict)

        # Should use default values for missing fields
        assert evidence.type == "blast"
        assert evidence.source_id == "test"
        assert evidence.domain_id == ""
        assert evidence.evalue is None
        assert evidence.confidence is not None  # Will auto-calculate

    def test_from_dict_preserves_types(self):
        """Test that from_dict preserves correct data types"""
        data_dict = {
            "type": "blast",
            "evalue": 1e-8,        # float
            "hsp_count": 5,        # int
            "confidence": 0.85,    # float
            "source_id": "hit_1"   # string
        }

        evidence = Evidence.from_dict(data_dict)

        # Types should be preserved
        assert isinstance(evidence.evalue, float)
        assert isinstance(evidence.hsp_count, int)
        assert isinstance(evidence.confidence, float)
        assert isinstance(evidence.source_id, str)


class TestAutoDetectionFromXML:
    """Test automatic evidence type detection from XML"""

    def test_auto_detect_hhsearch(self):
        """Test auto-detection of HHSearch evidence"""
        element = ET.Element("hit")
        element.set("probability", "85.0")
        element.set("score", "42.0")

        evidence = Evidence.from_xml(element)

        assert evidence.type == "hhsearch"
        assert evidence.probability == 85.0
        assert evidence.score == 42.0

    def test_auto_detect_domain_blast(self):
        """Test auto-detection of domain BLAST evidence"""
        element = ET.Element("hit")
        element.set("evalues", "1e-5")
        element.set("domain_id", "PF00001")

        evidence = Evidence.from_xml(element)

        assert evidence.type == "domain_blast"
        assert evidence.evalue == 1e-5
        assert evidence.domain_id == "PF00001"

    def test_auto_detect_chain_blast(self):
        """Test auto-detection of chain BLAST evidence"""
        element = ET.Element("hit")
        element.set("evalues", "1e-5")
        element.set("hsp_count", "3")
        # No domain_id - should be chain BLAST

        evidence = Evidence.from_xml(element)

        assert evidence.type == "chain_blast"
        assert evidence.evalue == 1e-5
        assert evidence.hsp_count == 3

    def test_auto_detect_unknown(self):
        """Test detection of unknown evidence type"""
        element = ET.Element("hit")
        element.set("some_custom_field", "value")
        # No recognizable attributes

        evidence = Evidence.from_xml(element)

        assert evidence.type == "unknown"


class TestSerializationErrorHandling:
    """Test error handling in serialization methods"""

    def test_malformed_xml_elements(self):
        """Test handling of malformed XML elements"""
        # Empty element
        empty_element = ET.Element("hit")
        evidence = Evidence.from_xml(empty_element)
        assert evidence.type == "unknown"

        # Element with only text, no attributes
        text_element = ET.Element("hit")
        text_element.text = "some random text"
        evidence = Evidence.from_xml(text_element)
        assert evidence.type == "unknown"

    def test_xml_with_invalid_numeric_values(self):
        """Test XML parsing with invalid numeric strings"""
        element = ET.Element("hit")
        element.set("evalues", "infinity")
        element.set("probability", "not_a_percentage")
        element.set("hsp_count", "lots")

        # Should not crash, should use fallback values
        evidence = Evidence.from_blast_xml(element)
        assert evidence.evalue == 999.0  # fallback
        assert evidence.hsp_count == 0   # fallback

    def test_confidence_recalculation_after_serialization(self):
        """Test that confidence is recalculated correctly after serialization"""
        # Create evidence without explicit confidence
        original = Evidence(type="blast", evalue=1e-8, confidence=None)
        original_confidence = original.confidence

        # Serialize and deserialize
        xml_element = original.to_xml()
        reconstructed = Evidence.from_xml(xml_element)

        # Confidence should be recalculated and match
        assert abs(reconstructed.confidence - original_confidence) < 0.001
        assert reconstructed._confidence_explicitly_set == True  # Because it was stored in XML
