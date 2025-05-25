#!/usr/bin/env python3
"""
HHSearch evidence processor for domain analysis.

This module provides processors for HHSearch results in various formats:
- HHSearch XML format (.xml)
- HHR format (.hhr)
"""

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Set

from .base import EvidenceProcessor, XMLProcessor, ValidationResult, ProcessingResult
from ecod.models.pipeline import Evidence
from ecod.exceptions import ValidationError


@dataclass
class HHSearchHitData:
    """Data for a single HHSearch hit"""
    hit_num: int
    hit_id: str
    domain_id: str = ""
    probability: float = 0.0
    evalue: float = 999.0
    score: float = 0.0
    query_start: int = 0
    query_end: int = 0
    hit_start: int = 0
    hit_end: int = 0
    query_ali: str = ""
    template_ali: str = ""
    query_range: str = ""
    hit_range: str = ""
    cols: int = 0
    identities: float = 0.0
    similarity: float = 0.0
    sum_probs: float = 0.0

    def calculate_confidence(self) -> float:
        """Calculate confidence score from probability"""
        # HHSearch probability is 0-100
        return self.probability / 100.0 if self.probability <= 100 else 1.0


class HHSearchEvidenceProcessor(XMLProcessor):
    """Process HHSearch results into Evidence objects"""

    def __init__(self, probability_threshold: float = 0.0,
                 max_hits: int = 100,
                 logger=None):
        """
        Initialize HHSearch processor.

        Args:
            probability_threshold: Minimum probability to accept hits
            max_hits: Maximum number of hits to process
            logger: Logger instance
        """
        super().__init__(logger)
        self.probability_threshold = probability_threshold
        self.max_hits = max_hits
        self.expected_root_tag = "hh_summ_doc"

    @property
    def supported_formats(self) -> Set[str]:
        """HHSearch processor supports .xml and .hhr files"""
        return {'.xml', '.hhr'}

    @property
    def evidence_type(self) -> str:
        """Return evidence type"""
        return "hhsearch"

    def validate_content(self, content: Any) -> ValidationResult:
        """Validate HHSearch content based on format"""
        if isinstance(content, ET.Element):
            # XML format validation
            return self._validate_xml_content(content)
        elif isinstance(content, dict):
            # HHR parsed format validation
            return self._validate_hhr_content(content)
        else:
            result = ValidationResult(valid=False)
            result.error = "Unknown content format"
            return result

    def _validate_xml_content(self, content: ET.Element) -> ValidationResult:
        """Validate HHSearch XML content"""
        result = super().validate_content(content)

        if not result.valid:
            return result

        # Check for hit list
        hit_list = content.find(".//hh_hit_list")
        if hit_list is None:
            result.error = "No hit list found in HHSearch XML"
            result.valid = False
            return result

        # Count hits
        hits = hit_list.findall("hh_hit")
        result.has_content = len(hits) > 0

        if len(hits) == 0:
            result.warnings.append("No hits found in HHSearch output")

        # Extract metadata
        result.metadata['hit_count'] = len(hits)
        result.metadata['format'] = 'xml'

        # Extract query info if available
        query_info = content.find(".//query")
        if query_info is not None:
            result.metadata['query_id'] = query_info.get("id", "")
            result.metadata['query_length'] = int(query_info.get("length", "0"))

        result.valid = True
        return result

    def _validate_hhr_content(self, content: Dict[str, Any]) -> ValidationResult:
        """Validate parsed HHR content"""
        result = ValidationResult(valid=False)

        # Check required keys
        if 'header' not in content or 'hits' not in content:
            result.error = "Invalid HHR content structure"
            return result

        result.format_valid = True

        # Check hits
        hits = content.get('hits', [])
        result.has_content = len(hits) > 0

        if len(hits) == 0:
            result.warnings.append("No hits found in HHR output")

        # Extract metadata
        header = content.get('header', {})
        result.metadata['hit_count'] = len(hits)
        result.metadata['format'] = 'hhr'
        result.metadata['query_id'] = header.get('query_id', '')
        result.metadata['match_columns'] = header.get('match_columns', 0)
        result.metadata['no_of_seqs'] = header.get('no_of_seqs', 0)

        result.valid = True
        return result

    def extract_evidence(self, content: Any, file_path: Path) -> List[Evidence]:
        """Extract evidence from HHSearch results"""
        if isinstance(content, ET.Element):
            return self._extract_from_xml(content, file_path)
        elif isinstance(content, dict):
            return self._extract_from_hhr(content, file_path)
        else:
            self.logger.error(f"Unknown content format for {file_path}")
            return []

    def _extract_from_xml(self, content: ET.Element, file_path: Path) -> List[Evidence]:
        """Extract evidence from HHSearch XML"""
        evidence_list = []

        # Find hit list
        hit_list = content.find(".//hh_hit_list")
        if hit_list is None:
            return []

        # Process hits
        for i, hit_elem in enumerate(hit_list.findall("hh_hit")):
            if i >= self.max_hits:
                break

            hit_data = self._parse_xml_hit(hit_elem)

            # Apply probability threshold
            if hit_data.probability < self.probability_threshold:
                continue

            # Create Evidence object
            evidence = self._create_evidence_from_hit(hit_data)
            evidence_list.append(evidence)

        self.logger.info(
            f"Extracted {len(evidence_list)} evidence items from {file_path.name}"
        )

        return evidence_list

    def _parse_xml_hit(self, hit_elem: ET.Element) -> HHSearchHitData:
        """Parse HHSearch XML hit element"""
        hit_data = HHSearchHitData(
            hit_num=int(hit_elem.get("num", "0")),
            hit_id=hit_elem.get("hit_id", ""),
            domain_id=hit_elem.get("domain_id", ""),
            probability=float(hit_elem.get("probability", "0")),
            evalue=float(hit_elem.get("evalue", "999")),
            score=float(hit_elem.get("score", "0"))
        )

        # Extract query regions
        query_reg = hit_elem.find("query_reg")
        if query_reg is not None and query_reg.text:
            hit_data.query_range = query_reg.text.strip()

        # Extract hit regions
        hit_reg = hit_elem.find("hit_reg")
        if hit_reg is not None and hit_reg.text:
            hit_data.hit_range = hit_reg.text.strip()

        # Extract additional metrics if available
        hit_data.cols = int(hit_elem.get("cols", "0"))
        hit_data.identities = float(hit_elem.get("identities", "0"))
        hit_data.similarity = float(hit_elem.get("similarity", "0"))
        hit_data.sum_probs = float(hit_elem.get("sum_probs", "0"))

        return hit_data

    def _extract_from_hhr(self, content: Dict[str, Any], file_path: Path) -> List[Evidence]:
        """Extract evidence from parsed HHR content"""
        evidence_list = []
        hits = content.get('hits', [])

        for i, hit_dict in enumerate(hits):
            if i >= self.max_hits:
                break

            hit_data = self._parse_hhr_hit(hit_dict)

            # Apply probability threshold
            if hit_data.probability < self.probability_threshold:
                continue

            # Create Evidence object
            evidence = self._create_evidence_from_hit(hit_data)
            evidence_list.append(evidence)

        self.logger.info(
            f"Extracted {len(evidence_list)} evidence items from {file_path.name}"
        )

        return evidence_list

    def _parse_hhr_hit(self, hit_dict: Dict[str, Any]) -> HHSearchHitData:
        """Parse HHR hit dictionary"""
        hit_data = HHSearchHitData(
            hit_num=hit_dict.get('hit_num', 0),
            hit_id=hit_dict.get('hit_id', ''),
            probability=hit_dict.get('probability', 0.0),
            evalue=hit_dict.get('e_value', 999.0),
            score=hit_dict.get('score', 0.0),
            query_start=hit_dict.get('query_start', 0),
            hit_start=hit_dict.get('template_start', 0),
            query_ali=hit_dict.get('query_ali', ''),
            template_ali=hit_dict.get('template_ali', '')
        )

        # Extract domain ID from hit_id (e.g., "e4tm9c1")
        domain_pattern = r'[ed]\d\w{3}\w\d+'
        match = re.search(domain_pattern, hit_data.hit_id)
        if match:
            hit_data.domain_id = match.group(0)

        # Calculate ranges from alignments
        if hit_data.query_ali and hit_data.query_start:
            hit_data.query_range = self._calculate_range(
                hit_data.query_ali,
                hit_data.query_start
            )

        if hit_data.template_ali and hit_data.hit_start:
            hit_data.hit_range = self._calculate_range(
                hit_data.template_ali,
                hit_data.hit_start
            )

        return hit_data

    def _create_evidence_from_hit(self, hit_data: HHSearchHitData) -> Evidence:
        """Create Evidence object from HHSearch hit data"""
        evidence = Evidence(
            type=self.evidence_type,
            source_id=hit_data.domain_id or hit_data.hit_id,
            domain_id=hit_data.domain_id,
            query_range=hit_data.query_range,
            hit_range=hit_data.hit_range,
            probability=hit_data.probability,
            evalue=hit_data.evalue,
            score=hit_data.score,
            confidence=None,  # Auto-calculated by Evidence model
            extra_attributes={
                "hit_id": hit_data.hit_id,
                "hit_num": hit_data.hit_num,
                "cols": hit_data.cols,
                "identities": hit_data.identities,
                "similarity": hit_data.similarity,
                "sum_probs": hit_data.sum_probs
            }
        )

        # Add alignment data if available
        if hit_data.query_ali:
            evidence.extra_attributes["query_ali"] = hit_data.query_ali
        if hit_data.template_ali:
            evidence.extra_attributes["template_ali"] = hit_data.template_ali

        return evidence

    def _parse_file(self, file_path: Path) -> Any:
        """Parse file based on extension"""
        if file_path.suffix.lower() == '.xml':
            return self._parse_xml(file_path)
        elif file_path.suffix.lower() == '.hhr':
            return self._parse_hhr_file(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_path.suffix}")

    def _parse_hhr_file(self, hhr_file: Path) -> Dict[str, Any]:
        """Parse HHR file to extract structured data"""
        cache_key = str(hhr_file)
        if cache_key in self._file_cache:
            return self._file_cache[cache_key]

        try:
            with open(hhr_file, 'r') as f:
                content = f.read()

            # Parse header information
            header = self._parse_hhr_header(content)

            # Find the beginning of the hit table
            lines = content.split('\n')
            hit_table_start = None

            for i, line in enumerate(lines):
                if line.startswith(' No Hit'):
                    hit_table_start = i + 1
                    break

            if not hit_table_start:
                return {'header': header, 'hits': []}

            # Process hits
            hits = self._parse_hhr_hits(lines[hit_table_start:])

            result = {
                'header': header,
                'hits': hits
            }

            # Cache result
            self._file_cache[cache_key] = result
            return result

        except Exception as e:
            self.logger.error(f"Error parsing HHR file {hhr_file}: {e}")
            return {'header': {}, 'hits': []}

    def _parse_hhr_header(self, content: str) -> Dict[str, Any]:
        """Parse HHR header information"""
        header = {}
        lines = content.split('\n')

        for line in lines:
            if line.startswith('Query '):
                parts = line.strip().split()
                if len(parts) > 1:
                    header['query_id'] = parts[1]
            elif line.startswith('Match_columns '):
                parts = line.strip().split()
                if len(parts) > 1:
                    header['match_columns'] = int(parts[1])
            elif line.startswith('No_of_seqs '):
                parts = line.strip().split()
                if len(parts) > 1:
                    header['no_of_seqs'] = int(parts[1])

            # Stop at the beginning of the hit table
            if line.startswith(' No Hit'):
                break

        return header

    def _parse_hhr_hits(self, lines: List[str]) -> List[Dict[str, Any]]:
        """Parse HHR hit section"""
        hits = []
        current_hit = None
        i = 0

        while i < len(lines):
            line = lines[i].strip()

            # New hit begins with a line like " 1 e4tm9c1 etc"
            match = re.match(r'^\s*(\d+)\s+(\S+)', line)
            if match:
                # Store previous hit if exists
                if current_hit and 'query_ali' in current_hit and 'template_ali' in current_hit:
                    hits.append(current_hit)

                # Parse hit line
                hit_num = int(match.group(1))
                hit_id = match.group(2)

                # Create new hit
                current_hit = {
                    'hit_num': hit_num,
                    'hit_id': hit_id,
                    'probability': None,
                    'e_value': None,
                    'score': None,
                    'query_ali': "",
                    'template_ali': ""
                }

                # Find probability, e-value, score in following lines
                j = i + 1
                while j < len(lines) and not lines[j].startswith('>'):
                    if 'Probab=' in lines[j]:
                        prob_match = re.search(r'Probab=(\d+\.\d+)', lines[j])
                        if prob_match:
                            current_hit['probability'] = float(prob_match.group(1))

                    if 'E-value=' in lines[j]:
                        eval_match = re.search(r'E-value=(\S+)', lines[j])
                        if eval_match:
                            try:
                                current_hit['e_value'] = float(eval_match.group(1))
                            except ValueError:
                                pass

                    if 'Score=' in lines[j]:
                        score_match = re.search(r'Score=(\d+\.\d+)', lines[j])
                        if score_match:
                            current_hit['score'] = float(score_match.group(1))

                    j += 1

            # Process alignment
            elif line.startswith('Q ') and current_hit:
                parts = line.split()
                if len(parts) >= 4:
                    if 'query_start' not in current_hit:
                        try:
                            current_hit['query_start'] = int(parts[2])
                        except ValueError:
                            pass
                    current_hit['query_ali'] += parts[3]

            elif line.startswith('T ') and current_hit:
                parts = line.split()
                if len(parts) >= 4:
                    if 'template_start' not in current_hit:
                        try:
                            current_hit['template_start'] = int(parts[2])
                        except ValueError:
                            pass
                    current_hit['template_ali'] += parts[3]

            i += 1

        # Add the last hit
        if current_hit and 'query_ali' in current_hit and 'template_ali' in current_hit:
            hits.append(current_hit)

        return hits

    def _calculate_range(self, alignment: str, start_pos: int) -> str:
        """
        Calculate range from alignment string.

        Args:
            alignment: Alignment string (with gaps)
            start_pos: Starting position

        Returns:
            Range string in format "start-end,start-end,..."
        """
        ranges = []
        current_range_start = None
        current_pos = start_pos

        for char in alignment:
            if char != '-':  # Not a gap
                if current_range_start is None:
                    current_range_start = current_pos
                current_pos += 1
            else:  # Gap
                if current_range_start is not None:
                    ranges.append(f"{current_range_start}-{current_pos-1}")
                    current_range_start = None

        # Add the last range if exists
        if current_range_start is not None:
            ranges.append(f"{current_range_start}-{current_pos-1}")

        return ",".join(ranges)

    def _calculate_coverage(self, query_ali: str, template_ali: str) -> float:
        """
        Calculate coverage between query and template alignments.

        Args:
            query_ali: Query alignment string
            template_ali: Template alignment string

        Returns:
            Coverage as a float between 0 and 1
        """
        if not query_ali or not template_ali or len(query_ali) != len(template_ali):
            return 0.0

        # Count aligned (non-gap) positions
        aligned_positions = sum(1 for q, t in zip(query_ali, template_ali) 
                              if q != '-' and t != '-')
        total_template_positions = sum(1 for t in template_ali if t != '-')

        if total_template_positions == 0:
            return 0.0

        return aligned_positions / total_template_positions


class HHSearchXMLProcessor(HHSearchEvidenceProcessor):
    """Specialized processor for HHSearch XML format only"""

    @property
    def supported_formats(self) -> Set[str]:
        """Only supports XML format"""
        return {'.xml'}

    def _parse_file(self, file_path: Path) -> Any:
        """Parse XML file only"""
        return self._parse_xml(file_path)


class HHRProcessor(HHSearchEvidenceProcessor):
    """Specialized processor for HHR format only"""

    @property
    def supported_formats(self) -> Set[str]:
        """Only supports HHR format"""
        return {'.hhr'}

    def _parse_file(self, file_path: Path) -> Any:
        """Parse HHR file only"""
        return self._parse_hhr_file(file_path)
