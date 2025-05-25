#!/usr/bin/env python3
"""
BLAST evidence processors for domain analysis.

This module provides processors for different types of BLAST results:
- Standard BLAST (domain vs domain database)
- Chain BLAST (full chain vs chain database)
- PSI-BLAST iterations
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
class HSPData:
    """Data for a single High-scoring Segment Pair"""
    query_from: int
    query_to: int
    hit_from: int
    hit_to: int
    evalue: float
    score: float = 0.0
    bit_score: float = 0.0
    identity: float = 0.0
    positive: float = 0.0
    gaps: int = 0
    align_len: int = 0
    query_seq: str = ""
    hit_seq: str = ""
    midline: str = ""

    @property
    def query_coverage(self) -> float:
        """Calculate query coverage percentage"""
        if self.align_len == 0:
            return 0.0
        return ((self.query_to - self.query_from + 1) / self.align_len) * 100

    @property
    def query_range(self) -> str:
        """Get query range as string"""
        return f"{self.query_from}-{self.query_to}"

    @property
    def hit_range(self) -> str:
        """Get hit range as string"""
        return f"{self.hit_from}-{self.hit_to}"


@dataclass
class BlastHitData:
    """Data for a complete BLAST hit (may contain multiple HSPs)"""
    hit_num: int
    hit_id: str
    hit_def: str
    hit_accession: str = ""
    hit_len: int = 0
    hsps: List[HSPData] = field(default_factory=list)

    # Parsed identifiers
    pdb_id: str = ""
    chain_id: str = ""
    domain_id: str = ""

    def parse_identifiers(self):
        """Parse PDB/chain/domain IDs from hit definition"""
        # Try to parse PDB and chain from hit_def
        if " " in self.hit_def:
            parts = self.hit_def.split()
            if len(parts) >= 2:
                self.pdb_id = parts[0][:4]  # First 4 chars
                self.chain_id = parts[1] if len(parts[1]) == 1 else parts[0][5:]

        # Try to parse domain ID (e.g., e4hluA1)
        domain_pattern = r'[ed]\d\w{3}[A-Z]\d+'
        match = re.search(domain_pattern, self.hit_def)
        if match:
            self.domain_id = match.group(0)

    @property
    def best_evalue(self) -> float:
        """Get best (lowest) e-value from all HSPs"""
        if not self.hsps:
            return 999.0
        return min(hsp.evalue for hsp in self.hsps)

    @property
    def total_score(self) -> float:
        """Get sum of all HSP scores"""
        return sum(hsp.score for hsp in self.hsps)


class BlastEvidenceProcessor(XMLProcessor):
    """Process BLAST XML results into Evidence objects"""

    def __init__(self, blast_type: str = "domain_blast",
                 hsp_evalue_threshold: float = 0.005,
                 hit_coverage_threshold: float = 0.7,
                 logger=None):
        """
        Initialize BLAST processor.

        Args:
            blast_type: Type of BLAST (domain_blast, chain_blast, etc.)
            hsp_evalue_threshold: Maximum e-value for HSPs
            hit_coverage_threshold: Minimum coverage for hits
            logger: Logger instance
        """
        super().__init__(logger)
        self.blast_type = blast_type
        self.hsp_evalue_threshold = hsp_evalue_threshold
        self.hit_coverage_threshold = hit_coverage_threshold
        self.expected_root_tag = "BlastOutput"

    @property
    def evidence_type(self) -> str:
        """Return evidence type based on BLAST type"""
        return self.blast_type

    def validate_content(self, content: ET.Element) -> ValidationResult:
        """Validate BLAST XML content"""
        result = super().validate_content(content)

        if not result.valid:
            return result

        # Check for required BLAST elements
        required = [
            ".//BlastOutput_program",
            ".//BlastOutput_version",
            ".//BlastOutput_query-def",
            ".//BlastOutput_query-len"
        ]

        for element_path in required:
            if content.find(element_path) is None:
                result.warnings.append(f"Missing BLAST element: {element_path}")

        # Check for iterations
        iterations = content.findall(".//Iteration")
        if not iterations:
            result.error = "No iterations found in BLAST output"
            result.valid = False
            return result

        # Check for hits
        total_hits = sum(len(iter.findall(".//Hit")) for iter in iterations)
        result.has_content = total_hits > 0

        if total_hits == 0:
            result.warnings.append("No hits found in BLAST output")

        # Extract metadata
        result.metadata['program'] = content.findtext(".//BlastOutput_program", "")
        result.metadata['version'] = content.findtext(".//BlastOutput_version", "")
        result.metadata['query_def'] = content.findtext(".//BlastOutput_query-def", "")
        result.metadata['query_len'] = int(content.findtext(".//BlastOutput_query-len", "0"))
        result.metadata['database'] = content.findtext(".//BlastOutput_db", "")
        result.metadata['iteration_count'] = len(iterations)
        result.metadata['total_hits'] = total_hits

        result.valid = True
        return result

    def extract_evidence(self, content: ET.Element, file_path: Path) -> List[Evidence]:
        """Extract evidence from BLAST XML"""
        evidence_list = []

        # Get query information
        query_def = content.findtext(".//BlastOutput_query-def", "")
        query_len = int(content.findtext(".//BlastOutput_query-len", "0"))

        # Extract query identifiers if available
        query_pdb, query_chain = self._parse_query_identifiers(query_def)

        # Process each iteration
        for iteration in content.findall(".//Iteration"):
            iter_num = int(iteration.findtext("Iteration_iter-num", "1"))

            # Process hits in this iteration
            for hit_elem in iteration.findall(".//Hit"):
                hit_data = self._parse_hit(hit_elem)

                # Process HSPs and create evidence
                evidence = self._process_hit_to_evidence(
                    hit_data, query_len, query_pdb, query_chain
                )

                evidence_list.extend(evidence)

        # Log summary
        self.logger.info(
            f"Extracted {len(evidence_list)} evidence items from {file_path.name}"
        )

        return evidence_list

    def _parse_hit(self, hit_elem: ET.Element) -> BlastHitData:
        """Parse a BLAST hit element"""
        hit_data = BlastHitData(
            hit_num=int(hit_elem.findtext("Hit_num", "0")),
            hit_id=hit_elem.findtext("Hit_id", ""),
            hit_def=hit_elem.findtext("Hit_def", ""),
            hit_accession=hit_elem.findtext("Hit_accession", ""),
            hit_len=int(hit_elem.findtext("Hit_len", "0"))
        )

        # Parse identifiers from hit definition
        hit_data.parse_identifiers()

        # Parse HSPs
        for hsp_elem in hit_elem.findall(".//Hsp"):
            hsp = self._parse_hsp(hsp_elem)
            if hsp:
                hit_data.hsps.append(hsp)

        return hit_data

    def _parse_hsp(self, hsp_elem: ET.Element) -> Optional[HSPData]:
        """Parse a single HSP element"""
        try:
            hsp = HSPData(
                query_from=int(hsp_elem.findtext("Hsp_query-from", "0")),
                query_to=int(hsp_elem.findtext("Hsp_query-to", "0")),
                hit_from=int(hsp_elem.findtext("Hsp_hit-from", "0")),
                hit_to=int(hsp_elem.findtext("Hsp_hit-to", "0")),
                evalue=float(hsp_elem.findtext("Hsp_evalue", "999")),
                score=float(hsp_elem.findtext("Hsp_score", "0")),
                bit_score=float(hsp_elem.findtext("Hsp_bit-score", "0")),
                identity=float(hsp_elem.findtext("Hsp_identity", "0")),
                positive=float(hsp_elem.findtext("Hsp_positive", "0")),
                gaps=int(hsp_elem.findtext("Hsp_gaps", "0")),
                align_len=int(hsp_elem.findtext("Hsp_align-len", "0")),
                query_seq=hsp_elem.findtext("Hsp_qseq", ""),
                hit_seq=hsp_elem.findtext("Hsp_hseq", ""),
                midline=hsp_elem.findtext("Hsp_midline", "")
            )

            # Calculate identity percentage
            if hsp.align_len > 0:
                hsp.identity = (hsp.identity / hsp.align_len) * 100

            return hsp

        except (ValueError, TypeError) as e:
            self.logger.warning(f"Error parsing HSP: {e}")
            return None

    def _process_hit_to_evidence(self, hit_data: BlastHitData, query_len: int,
                                query_pdb: str = "", query_chain: str = "") -> List[Evidence]:
        """Convert BLAST hit to Evidence objects"""
        if not hit_data.hsps:
            return []

        # Filter HSPs by e-value threshold
        valid_hsps = [
            hsp for hsp in hit_data.hsps
            if hsp.evalue < self.hsp_evalue_threshold
        ]

        if not valid_hsps:
            return []

        # For domain BLAST, stitch HSPs if needed
        if self.blast_type == "domain_blast" and len(valid_hsps) > 1:
            return self._create_stitched_evidence(hit_data, valid_hsps, query_len)
        else:
            # Create individual evidence for each HSP
            return self._create_individual_evidence(hit_data, valid_hsps)

    def _create_individual_evidence(self, hit_data: BlastHitData,
                                  hsps: List[HSPData]) -> List[Evidence]:
        """Create individual Evidence objects for each HSP"""
        evidence_list = []

        for hsp in hsps:
            evidence = Evidence(
                type=self.blast_type,
                source_id=hit_data.domain_id or hit_data.hit_id,
                domain_id=hit_data.domain_id,
                query_range=hsp.query_range,
                hit_range=hsp.hit_range,
                evalue=hsp.evalue,
                score=hsp.score,
                identity=hsp.identity,
                coverage=hsp.query_coverage,
                confidence=None,  # Auto-calculated
                extra_attributes={
                    "pdb_id": hit_data.pdb_id,
                    "chain_id": hit_data.chain_id,
                    "hit_def": hit_data.hit_def,
                    "query_seq": hsp.query_seq,
                    "hit_seq": hsp.hit_seq,
                    "align_len": hsp.align_len,
                    "bit_score": hsp.bit_score,
                    "gaps": hsp.gaps
                }
            )

            evidence_list.append(evidence)

        return evidence_list

    def _create_stitched_evidence(self, hit_data: BlastHitData,
                                hsps: List[HSPData], query_len: int) -> List[Evidence]:
        """Create stitched Evidence for discontinuous domains"""
        # Group HSPs that should be stitched together
        hsp_groups = self._group_hsps_for_stitching(hsps, hit_data.domain_id)

        evidence_list = []

        for group in hsp_groups:
            if len(group) == 1:
                # Single HSP - create normal evidence
                evidence_list.extend(
                    self._create_individual_evidence(hit_data, group)
                )
            else:
                # Multiple HSPs - create stitched evidence
                evidence = self._create_merged_evidence(hit_data, group)
                if evidence:
                    evidence_list.append(evidence)

        return evidence_list

    def _group_hsps_for_stitching(self, hsps: List[HSPData],
                                 domain_id: str = "") -> List[List[HSPData]]:
        """Group HSPs that represent parts of the same discontinuous domain"""
        # Sort HSPs by query start position
        sorted_hsps = sorted(hsps, key=lambda x: x.query_from)

        # Parameters for grouping
        max_gap = 30  # Maximum gap between HSPs to consider stitching

        groups = []
        current_group = [sorted_hsps[0]]

        for i in range(1, len(sorted_hsps)):
            hsp = sorted_hsps[i]
            prev_hsp = current_group[-1]

            # Calculate gap between HSPs
            query_gap = hsp.query_from - prev_hsp.query_to

            # Check if HSPs should be grouped
            if query_gap <= max_gap and self._hsps_compatible(prev_hsp, hsp):
                current_group.append(hsp)
            else:
                # Start new group
                groups.append(current_group)
                current_group = [hsp]

        # Add final group
        if current_group:
            groups.append(current_group)

        self.logger.debug(
            f"Grouped {len(hsps)} HSPs into {len(groups)} groups for {domain_id}"
        )

        return groups

    def _hsps_compatible(self, hsp1: HSPData, hsp2: HSPData) -> bool:
        """Check if two HSPs are compatible for stitching"""
        # Check that hit coordinates also make sense
        hit_gap = abs(hsp2.hit_from - hsp1.hit_to)

        # HSPs should maintain relative order in both query and hit
        query_order_ok = hsp2.query_from > hsp1.query_to
        hit_order_ok = hsp2.hit_from > hsp1.hit_to

        # E-values should be reasonably similar (within an order of magnitude)
        evalue_ratio = max(hsp1.evalue, hsp2.evalue) / min(hsp1.evalue, hsp2.evalue)
        evalue_similar = evalue_ratio < 10

        return query_order_ok and hit_order_ok and evalue_similar and hit_gap < 100

    def _create_merged_evidence(self, hit_data: BlastHitData,
                               hsp_group: List[HSPData]) -> Optional[Evidence]:
        """Create a single Evidence object from multiple HSPs"""
        if not hsp_group:
            return None

        # Combine ranges
        query_ranges = [hsp.query_range for hsp in hsp_group]
        hit_ranges = [hsp.hit_range for hsp in hsp_group]

        # Use best (lowest) e-value
        best_evalue = min(hsp.evalue for hsp in hsp_group)

        # Calculate combined statistics
        total_identity = sum(hsp.identity * hsp.align_len for hsp in hsp_group)
        total_align_len = sum(hsp.align_len for hsp in hsp_group)
        avg_identity = total_identity / total_align_len if total_align_len > 0 else 0

        # Calculate coverage
        query_positions = set()
        for hsp in hsp_group:
            query_positions.update(range(hsp.query_from, hsp.query_to + 1))
        coverage = len(query_positions)  # Will be converted to percentage by Evidence

        # Create merged evidence
        evidence = Evidence(
            type=self.blast_type,
            source_id=hit_data.domain_id or hit_data.hit_id,
            domain_id=hit_data.domain_id,
            query_range=",".join(query_ranges),
            hit_range=",".join(hit_ranges),
            evalue=best_evalue,
            score=sum(hsp.score for hsp in hsp_group),
            identity=avg_identity,
            coverage=coverage,
            hsp_count=len(hsp_group),
            confidence=None,  # Auto-calculated
            extra_attributes={
                "pdb_id": hit_data.pdb_id,
                "chain_id": hit_data.chain_id,
                "hit_def": hit_data.hit_def,
                "discontinuous": True,
                "segment_count": len(hsp_group),
                "gap_info": self._calculate_gap_info(hsp_group)
            }
        )

        return evidence

    def _calculate_gap_info(self, hsp_group: List[HSPData]) -> str:
        """Calculate gap information for stitched HSPs"""
        gaps = []

        for i in range(len(hsp_group) - 1):
            current = hsp_group[i]
            next_hsp = hsp_group[i + 1]
            gap_size = next_hsp.query_from - current.query_to - 1
            gaps.append(str(gap_size))

        return ",".join(gaps)

    def _parse_query_identifiers(self, query_def: str) -> Tuple[str, str]:
        """Parse PDB and chain IDs from query definition"""
        pdb_id = ""
        chain_id = ""

        # Try common patterns
        # Pattern 1: "1ABC_A" or "1ABC:A"
        match = re.match(r'(\w{4})[_:](\w)', query_def)
        if match:
            pdb_id = match.group(1)
            chain_id = match.group(2)
        # Pattern 2: "pdb|1ABC|A"
        else:
            match = re.search(r'pdb\|(\w{4})\|(\w)', query_def)
            if match:
                pdb_id = match.group(1)
                chain_id = match.group(2)

        return pdb_id, chain_id


class ChainBlastProcessor(BlastEvidenceProcessor):
    """
    Specialized processor for chain BLAST results.

    Chain BLAST requires additional processing:
    - HSP coverage validation
    - Residue usage tracking
    - Special handling for multi-domain chains
    """

    def __init__(self, hsp_evalue_threshold: float = 0.005,
                 hit_coverage_threshold: float = 0.7,
                 hit_diff_tolerance: int = 50,
                 query_diff_tolerance: int = 50,
                 logger=None):
        """
        Initialize chain BLAST processor.

        Args:
            hsp_evalue_threshold: Maximum e-value for HSPs
            hit_coverage_threshold: Minimum coverage for hits
            hit_diff_tolerance: Maximum difference between hit length and alignment
            query_diff_tolerance: Maximum difference between query length and hit length
            logger: Logger instance
        """
        super().__init__(
            blast_type="chain_blast",
            hsp_evalue_threshold=hsp_evalue_threshold,
            hit_coverage_threshold=hit_coverage_threshold,
            logger=logger
        )
        self.hit_diff_tolerance = hit_diff_tolerance
        self.query_diff_tolerance = query_diff_tolerance

    def _process_hit_to_evidence(self, hit_data: BlastHitData, query_len: int,
                                query_pdb: str = "", query_chain: str = "") -> List[Evidence]:
        """Process chain BLAST hit with coverage validation"""
        if not hit_data.hsps:
            return []

        # Track residue usage
        used_query_positions = set()
        used_hit_positions = set()

        # Filter and validate HSPs
        valid_hsps = []

        for hsp in hit_data.hsps:
            # Check e-value
            if hsp.evalue >= self.hsp_evalue_threshold:
                continue

            # Check for overlap with previously used regions
            hsp_query_positions = set(range(hsp.query_from, hsp.query_to + 1))
            hsp_hit_positions = set(range(hsp.hit_from, hsp.hit_to + 1))

            query_overlap = len(hsp_query_positions & used_query_positions)
            hit_overlap = len(hsp_hit_positions & used_hit_positions)

            # Allow small overlaps but not large ones
            if query_overlap > 10 or hit_overlap > 10:
                self.logger.debug(
                    f"Skipping HSP due to overlap: query={query_overlap}, hit={hit_overlap}"
                )
                continue

            # Add to valid HSPs
            valid_hsps.append(hsp)

            # Update used positions
            used_query_positions.update(hsp_query_positions)
            used_hit_positions.update(hsp_hit_positions)

        if not valid_hsps:
            return []

        # Check overall coverage requirements
        total_align_len = sum(hsp.align_len for hsp in valid_hsps)
        hit_coverage = len(used_hit_positions) / hit_data.hit_len if hit_data.hit_len > 0 else 0

        # Apply coverage validation for chain BLAST
        if abs(total_align_len - hit_data.hit_len) > self.hit_diff_tolerance:
            self.logger.debug(
                f"Hit coverage mismatch: align_len={total_align_len}, "
                f"hit_len={hit_data.hit_len}"
            )
            return []

        if abs(query_len - hit_data.hit_len) > self.query_diff_tolerance:
            self.logger.debug(
                f"Query/hit length mismatch: query_len={query_len}, "
                f"hit_len={hit_data.hit_len}"
            )
            return []

        # Create evidence with aggregated HSP data
        evidence = Evidence(
            type=self.blast_type,
            source_id=f"{hit_data.pdb_id}_{hit_data.chain_id}",
            domain_id="",  # Chain BLAST doesn't have domain IDs
            query_range=",".join([hsp.query_range for hsp in valid_hsps]),
            hit_range=",".join([hsp.hit_range for hsp in valid_hsps]),
            evalue=min(hsp.evalue for hsp in valid_hsps),
            score=sum(hsp.score for hsp in valid_hsps),
            identity=self._calculate_avg_identity(valid_hsps),
            coverage=hit_coverage * 100,  # Convert to percentage
            hsp_count=len(valid_hsps),
            confidence=None,  # Auto-calculated
            extra_attributes={
                "pdb_id": hit_data.pdb_id,
                "chain_id": hit_data.chain_id,
                "hit_def": hit_data.hit_def,
                "query_coverage": (len(used_query_positions) / query_len * 100) if query_len > 0 else 0,
                "hit_coverage": hit_coverage * 100,
                "query_seqs": [hsp.query_seq for hsp in valid_hsps],
                "hit_seqs": [hsp.hit_seq for hsp in valid_hsps]
            }
        )

        return [evidence]

    def _calculate_avg_identity(self, hsps: List[HSPData]) -> float:
        """Calculate weighted average identity across HSPs"""
        if not hsps:
            return 0.0

        total_identity = sum(hsp.identity * hsp.align_len for hsp in hsps)
        total_align_len = sum(hsp.align_len for hsp in hsps)

        return total_identity / total_align_len if total_align_len > 0 else 0.0


class PSIBlastProcessor(BlastEvidenceProcessor):
    """Process PSI-BLAST results with iteration handling"""

    def __init__(self, min_iteration: int = 1, max_iteration: int = 5, **kwargs):
        """
        Initialize PSI-BLAST processor.

        Args:
            min_iteration: Minimum iteration to consider
            max_iteration: Maximum iteration to consider
            **kwargs: Additional arguments for parent class
        """
        super().__init__(blast_type="psi_blast", **kwargs)
        self.min_iteration = min_iteration
        self.max_iteration = max_iteration

    def extract_evidence(self, content: ET.Element, file_path: Path) -> List[Evidence]:
        """Extract evidence from PSI-BLAST iterations"""
        evidence_list = []

        # Get query information
        query_len = int(content.findtext(".//BlastOutput_query-len", "0"))

        # Process specified iterations
        for iteration in content.findall(".//Iteration"):
            iter_num = int(iteration.findtext("Iteration_iter-num", "1"))

            # Skip iterations outside range
            if iter_num < self.min_iteration or iter_num > self.max_iteration:
                continue

            # Process hits in this iteration
            for hit_elem in iteration.findall(".//Hit"):
                hit_data = self._parse_hit(hit_elem)

                # Add iteration info
                for evidence in self._process_hit_to_evidence(hit_data, query_len):
                    evidence.extra_attributes["iteration"] = iter_num
                    evidence_list.append(evidence)

        return evidence_list


# Factory function
def create_blast_processor(blast_type: str = "domain_blast", **kwargs) -> BlastEvidenceProcessor:
    """
    Create appropriate BLAST processor based on type.

    Args:
        blast_type: Type of BLAST processor to create
        **kwargs: Additional arguments for processor

    Returns:
        BlastEvidenceProcessor instance
    """
    if blast_type == "chain_blast":
        return ChainBlastProcessor(**kwargs)
    elif blast_type == "psi_blast":
        return PSIBlastProcessor(**kwargs)
    else:
        return BlastEvidenceProcessor(blast_type=blast_type, **kwargs)
