#!/usr/bin/env python3
"""
Enhanced Mini PyECOD Results Importer v2.0

Updated to handle the new comprehensive provenance tracking and evidence metrics
from the enhanced writer module. Imports results as mini_pyecod_v2.0.

Key enhancements:
- Handles new XML structure with comprehensive metadata
- Extracts and stores git version information and file hashes
- Processes enhanced evidence information with quality metrics
- Stores boundary optimization data
- Collision-safe import with version tracking
- Separated propagation logic (recommend SQL-only approach)

Usage:
    python enhanced_mini_importer.py --assess-and-import --tier-filter excellent,good
    python enhanced_mini_importer.py --comparative-analysis
"""

import os
import sys
import argparse
import yaml
import psycopg2
import psycopg2.extras
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass, asdict
from datetime import datetime
import logging
from collections import defaultdict
import hashlib

# Import your existing quality assessment
sys.path.append(str(Path(__file__).parent))
from assess_quality_v2 import EnhancedQualityAssessment, ProteinQuality

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class EnhancedMiniResult:
    """Enhanced representation of a mini PyECOD result with full provenance"""
    protein_id: str
    pdb_id: str
    chain_id: str
    is_classified: bool
    
    # Enhanced metadata from new XML format
    algorithm_version: Optional[str] = None
    git_commit_hash: Optional[str] = None
    processing_timestamp: Optional[datetime] = None
    source_xml_path: Optional[str] = None
    source_xml_hash: Optional[str] = None
    batch_id: Optional[str] = None
    
    # Sequence and coverage statistics
    sequence_length: Optional[int] = None
    total_coverage: Optional[float] = None
    residues_assigned: Optional[int] = None
    domains_optimized: Optional[int] = None
    
    # Processing parameters
    process_parameters: Optional[Dict] = None
    
    # Domain data
    domains: List['EnhancedDomain'] = None
    
    def __post_init__(self):
        if self.domains is None:
            self.domains = []
        if self.process_parameters is None:
            self.process_parameters = {}

@dataclass
class EnhancedEvidence:
    """Enhanced evidence information with quality metrics"""
    source_type: str
    source_id: str
    domain_id: Optional[str] = None
    confidence: Optional[float] = None
    
    # Type-specific metrics
    hh_probability: Optional[float] = None  # Original HHsearch probability
    evalue: Optional[float] = None
    reference_coverage: Optional[float] = None
    reference_length: Optional[int] = None
    hsp_count: Optional[int] = None
    
    # Quality flags
    quality_flags: List[str] = None
    
    # Ranges
    evidence_range: Optional[str] = None
    hit_range: Optional[str] = None
    discontinuous: bool = False
    
    def __post_init__(self):
        if self.quality_flags is None:
            self.quality_flags = []

@dataclass
class EnhancedDomain:
    """Enhanced domain with boundary optimization and evidence details"""
    domain_id: str
    range_str: str
    family: str
    source: str
    evidence_count: int
    is_discontinuous: bool
    
    # Classification
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    x_group: Optional[str] = None
    confidence: Optional[float] = None
    reference_ecod_domain_id: Optional[str] = None
    
    # Enhanced evidence
    primary_evidence: Optional[EnhancedEvidence] = None
    supporting_evidence_count: Optional[int] = None
    average_confidence: Optional[float] = None
    
    # Boundary optimization
    was_optimized: bool = False
    original_range: Optional[str] = None
    optimization_actions: List[str] = None
    position_change: Optional[int] = None
    
    def __post_init__(self):
        if self.optimization_actions is None:
            self.optimization_actions = []

class EnhancedMiniImporter:
    """Enhanced importer for mini PyECOD v2.0 results with comprehensive provenance"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.db_conn = psycopg2.connect(**self.config["database"])
        self.quality_assessor = EnhancedQualityAssessment(config_path)
        
        # Version for this import
        self.process_version = "mini_pyecod_v2.0"
        
    def parse_enhanced_mini_xml(self, xml_file: Path) -> EnhancedMiniResult:
        """Parse enhanced mini XML with comprehensive metadata"""
        
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # Basic identification
        pdb_id = root.get("pdb_id")
        chain_id = root.get("chain_id")
        protein_id = f"{pdb_id}_{chain_id}"
        is_classified = root.get("is_classified", "false").lower() == "true"
        
        result = EnhancedMiniResult(
            protein_id=protein_id,
            pdb_id=pdb_id,
            chain_id=chain_id,
            is_classified=is_classified,
            source_xml_path=str(xml_file)
        )
        
        # Calculate XML file hash
        result.source_xml_hash = self._calculate_file_hash(xml_file)
        
        # Parse enhanced metadata
        metadata_elem = root.find("metadata")
        if metadata_elem is not None:
            self._parse_metadata(metadata_elem, result)
        
        # Parse domains with enhanced information
        domains_elem = root.find("domains")
        if domains_elem is not None:
            result.domains = self._parse_enhanced_domains(domains_elem)
        
        return

    def _calculate_file_hash(self, file_path: Path) -> str:
        """Calculate SHA256 hash of file"""
        hash_sha256 = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_sha256.update(chunk)
        return hash_sha256.hexdigest()

    def _parse_metadata(self, metadata_elem: ET.Element, result: EnhancedMiniResult):
        """Parse comprehensive metadata from XML"""

        # Version information
        version_elem = metadata_elem.find("version")
        if version_elem is not None:
            result.algorithm_version = version_elem.get("algorithm")
            result.git_commit_hash = version_elem.get("git_commit")
            timestamp_str = version_elem.get("timestamp")
            if timestamp_str:
                result.processing_timestamp = datetime.fromisoformat(timestamp_str)

        # Source provenance
        source_elem = metadata_elem.find("source")
        if source_elem is not None:
            result.batch_id = source_elem.get("batch_id")

        # Processing parameters
        params_elem = metadata_elem.find("parameters")
        if params_elem is not None:
            result.process_parameters = {}
            for param in params_elem.findall("parameter"):
                name = param.get("name")
                value = param.get("value")
                # Try to convert to appropriate type
                try:
                    if value.lower() in ('true', 'false'):
                        value = value.lower() == 'true'
                    elif '.' in value:
                        value = float(value)
                    else:
                        value = int(value)
                except (ValueError, AttributeError):
                    pass  # Keep as string
                result.process_parameters[name] = value

        # Statistics
        stats_elem = metadata_elem.find("statistics")
        if stats_elem is not None:
            result.sequence_length = self._get_int_attr(stats_elem, "sequence_length")
            result.total_coverage = self._get_float_attr(stats_elem, "total_coverage")
            result.residues_assigned = self._get_int_attr(stats_elem, "residues_assigned")
            result.domains_optimized = self._get_int_attr(stats_elem, "domains_optimized")

    def _parse_enhanced_domains(self, domains_elem: ET.Element) -> List[EnhancedDomain]:
        """Parse domains with enhanced evidence and optimization data"""

        domains = []

        for domain_elem in domains_elem.findall("domain"):
            domain = EnhancedDomain(
                domain_id=domain_elem.get("id"),
                range_str=domain_elem.get("range"),
                family=domain_elem.get("family"),
                source=domain_elem.get("source"),
                evidence_count=self._get_int_attr(domain_elem, "evidence_count", 0),
                is_discontinuous=domain_elem.get("is_discontinuous", "false").lower() == "true"
            )

            # Classification hierarchy
            domain.t_group = domain_elem.get("t_group")
            domain.h_group = domain_elem.get("h_group")
            domain.x_group = domain_elem.get("x_group")
            domain.confidence = self._get_float_attr(domain_elem, "confidence")
            domain.reference_ecod_domain_id = domain_elem.get("reference_ecod_domain_id")

            # Parse primary evidence
            evidence_elem = domain_elem.find("primary_evidence")
            if evidence_elem is not None:
                domain.primary_evidence = self._parse_enhanced_evidence(evidence_elem)

            # Parse supporting evidence summary
            support_elem = domain_elem.find("supporting_evidence")
            if support_elem is not None:
                domain.supporting_evidence_count = self._get_int_attr(support_elem, "count")
                domain.average_confidence = self._get_float_attr(support_elem, "average_confidence")

            # Parse boundary optimization
            opt_elem = domain_elem.find("boundary_optimization")
            if opt_elem is not None:
                domain.was_optimized = True
                domain.original_range = opt_elem.get("original_range")
                domain.position_change = self._get_int_attr(opt_elem, "position_change")
                actions = opt_elem.get("actions")
                if actions:
                    domain.optimization_actions = actions.split(",")

            domains.append(domain)

        return domains

    def _parse_enhanced_evidence(self, evidence_elem: ET.Element) -> EnhancedEvidence:
        """Parse enhanced evidence with quality metrics"""

        evidence = EnhancedEvidence(
            source_type=evidence_elem.get("source_type"),
            source_id=evidence_elem.get("source_id"),
            domain_id=evidence_elem.get("domain_id"),
            confidence=self._get_float_attr(evidence_elem, "confidence")
        )

        # Type-specific metrics
        evidence.hh_probability = self._get_float_attr(evidence_elem, "hh_probability")
        evidence.evalue = self._get_float_attr(evidence_elem, "evalue")
        evidence.reference_coverage = self._get_float_attr(evidence_elem, "reference_coverage")
        evidence.reference_length = self._get_int_attr(evidence_elem, "reference_length")
        evidence.hsp_count = self._get_int_attr(evidence_elem, "hsp_count")

        # Ranges
        evidence.evidence_range = evidence_elem.get("evidence_range")
        evidence.hit_range = evidence_elem.get("hit_range")
        evidence.discontinuous = evidence_elem.get("discontinuous", "false").lower() == "true"

        # Quality flags
        quality_flags = evidence_elem.get("quality_flags")
        if quality_flags:
            evidence.quality_flags = quality_flags.split(",")

        return evidence

    def _get_int_attr(self, elem: ET.Element, attr: str, default: Optional[int] = None) -> Optional[int]:
        """Safely get integer attribute"""
        value = elem.get(attr)
        if value is None:
            return default
        try:
            return int(value)
        except ValueError:
            return default

    def _get_float_attr(self, elem: ET.Element, attr: str, default: Optional[float] = None) -> Optional[float]:
        """Safely get float attribute"""
        value = elem.get(attr)
        if value is None:
            return default
        try:
            return float(value)
        except ValueError:
            return default

    def import_enhanced_result(self, result: EnhancedMiniResult,
                             collision_strategy: str = "separate") -> bool:
        """Import enhanced mini result with comprehensive provenance tracking"""

        try:
            with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                # Check for existing imports
                if collision_strategy == "skip":
                    cursor.execute("""
                        SELECT id FROM pdb_analysis.partition_proteins
                        WHERE pdb_id = %s AND chain_id = %s AND process_version = %s
                    """, (result.pdb_id, result.chain_id, self.process_version))

                    if cursor.fetchone():
                        logger.info(f"Skipping {result.protein_id} - already exists")
                        return True

                # Look up batch_id from batch_name
                batch_id_int = None
                batch_name = result.batch_id  # This is actually batch_name from XML

                if batch_name:
                    # Look up the batch in ecod_schema.batch table
                    cursor.execute("""
                        SELECT id FROM ecod_schema.batch
                        WHERE batch_name = %s
                    """, (batch_name,))

                    batch_row = cursor.fetchone()
                    if batch_row:
                        batch_id_int = batch_row['id']
                        logger.debug(f"Found batch_id {batch_id_int} for batch_name '{batch_name}'")
                    else:
                        # Batch doesn't exist in the table
                        logger.warning(f"Batch '{batch_name}' not found in ecod_schema.batch table")
                        # For now, we'll leave batch_id as NULL and store batch_name in parameters
                        batch_id_int = None

                # Prepare enhanced process parameters with batch name
                enhanced_params = result.process_parameters.copy() if result.process_parameters else {}
                if batch_name:
                    enhanced_params['source_batch_name'] = batch_name
                enhanced_params['import_timestamp'] = datetime.now().isoformat()
                enhanced_params['import_version'] = self.process_version

                # Insert partition_proteins record with enhanced metadata
                cursor.execute("""
                    INSERT INTO pdb_analysis.partition_proteins (
                        pdb_id, chain_id, batch_id, reference_version, is_classified,
                        sequence_length, coverage, residues_assigned, domains_with_evidence,
                        fully_classified_domains, process_version, algorithm_version,
                        git_commit_hash, source_xml_path, source_xml_hash,
                        domain_summary_path, domain_summary_hash, process_parameters
                    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    RETURNING id
                """, (
                    result.pdb_id,
                    result.chain_id,
                    batch_id_int,  # Use looked-up batch_id (may be NULL)
                    "mini_pyecod_v2.0",  # reference_version
                    result.is_classified,
                    result.sequence_length,
                    result.total_coverage,
                    result.residues_assigned,
                    len([d for d in result.domains if d.primary_evidence]),
                    len([d for d in result.domains if d.t_group]),
                    self.process_version,
                    result.algorithm_version,
                    result.git_commit_hash,
                    result.source_xml_path,
                    result.source_xml_hash,
                    None,  # domain_summary_path (not applicable for v2.0)
                    None,  # domain_summary_hash
                    psycopg2.extras.Json(enhanced_params)
                ))

                partition_protein_id = cursor.fetchone()[0]

                # Insert domains with enhanced data
                for i, domain in enumerate(result.domains, 1):
                    self._import_enhanced_domain(cursor, partition_protein_id, i, domain, result)

                self.db_conn.commit()
                return True

        except Exception as e:
            logger.error(f"Error importing {result.protein_id}: {e}")
            self.db_conn.rollback()
            return False

    def create_missing_batch_if_needed(self, batch_name: str) -> Optional[int]:
        """Create a batch record if it doesn't exist and return its ID"""

        if not batch_name:
            return None

        try:
            with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                # Check if batch already exists
                cursor.execute("""
                    SELECT id FROM ecod_schema.batch WHERE batch_name = %s
                """, (batch_name,))

                existing = cursor.fetchone()
                if existing:
                    return existing['id']

                # Create new batch record
                cursor.execute("""
                    INSERT INTO ecod_schema.batch (
                        batch_name, base_path, type, ref_version,
                        total_items, completed_items, status, created_at
                    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                    RETURNING id
                """, (
                    batch_name,
                    f"/imported/mini_v2/{batch_name}",  # Placeholder path
                    "mini_import_v2",  # Type to distinguish imported batches
                    "mini_pyecod_v2.0",  # Reference version
                    0,  # total_items (unknown for imported batches)
                    0,  # completed_items
                    "imported",  # Status
                    datetime.now()
                ))

                new_batch_id = cursor.fetchone()['id']
                self.db_conn.commit()

                logger.info(f"Created new batch record: {batch_name} -> ID {new_batch_id}")
                return new_batch_id

        except Exception as e:
            logger.error(f"Error creating batch record for '{batch_name}': {e}")
            self.db_conn.rollback()
            return None

    def import_enhanced_result_with_batch_creation(self, result: EnhancedMiniResult,
                                                 collision_strategy: str = "separate",
                                                 create_missing_batches: bool = True) -> bool:
        """Import with option to create missing batch records"""

        # First try normal import
        success = self.import_enhanced_result(result, collision_strategy)

        # If it failed and we allow batch creation, try creating the batch
        if not success and create_missing_batches and result.batch_id:
            logger.info(f"Attempting to create missing batch: {result.batch_id}")
            batch_id = self.create_missing_batch_if_needed(result.batch_id)
            if batch_id:
                # Try import again
                success = self.import_enhanced_result(result, collision_strategy)

        return success

    def _import_enhanced_domain(self, cursor, partition_protein_id: int, domain_number: int,
                              domain: EnhancedDomain, result: EnhancedMiniResult):
        """Import domain with enhanced evidence and optimization data"""

        # Parse range for start/end positions
        start_pos, end_pos = self._parse_range_positions(domain.range_str)

        # Insert partition_domains record
        cursor.execute("""
            INSERT INTO pdb_analysis.partition_domains (
                protein_id, domain_number, domain_id, start_pos, end_pos, range,
                source, source_id, confidence, t_group, h_group, x_group, a_group,
                is_manual_rep, is_f70, is_f40, is_f99, created_at,
                primary_evidence_type, primary_evidence_id, evidence_evalue,
                evidence_query_range, evidence_hit_range, original_range,
                optimization_actions
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id
        """, (
            partition_protein_id,
            domain_number,
            domain.domain_id,
            start_pos,
            end_pos,
            domain.range_str,
            domain.source,
            domain.reference_ecod_domain_id,
            domain.confidence,
            domain.t_group,
            domain.h_group,
            domain.x_group,
            None,  # a_group
            False,  # is_manual_rep
            False,  # is_f70
            False,  # is_f40
            False,  # is_f99
            result.processing_timestamp,
            domain.primary_evidence.source_type if domain.primary_evidence else None,
            domain.primary_evidence.source_id if domain.primary_evidence else None,
            domain.primary_evidence.evalue if domain.primary_evidence else None,
            domain.primary_evidence.evidence_range if domain.primary_evidence else None,
            domain.primary_evidence.hit_range if domain.primary_evidence else None,
            domain.original_range if domain.was_optimized else None,
            domain.optimization_actions if domain.optimization_actions else None
        ))

        domain_id = cursor.fetchone()[0]

        # Insert enhanced evidence record
        if domain.primary_evidence:
            self._import_enhanced_evidence(cursor, domain_id, domain.primary_evidence)

    def _import_enhanced_evidence(self, cursor, domain_id: int, evidence: EnhancedEvidence):
        """Import enhanced evidence with quality metrics"""

        cursor.execute("""
            INSERT INTO pdb_analysis.domain_evidence (
                domain_id, evidence_type, source_id, domain_ref_id, hit_id,
                confidence, probability, evalue, score, hsp_count,
                is_discontinuous, query_range, hit_range, created_at
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """, (
            domain_id,
            evidence.source_type,
            evidence.source_id,
            evidence.domain_id,
            evidence.source_id,  # hit_id
            evidence.confidence,
            evidence.hh_probability,  # Use as probability
            evidence.evalue,
            None,  # score (not available in v2.0 format)
            evidence.hsp_count,
            evidence.discontinuous,
            evidence.evidence_range,
            evidence.hit_range,
            datetime.now()
        ))

    def _parse_range_positions(self, range_str: str) -> Tuple[int, int]:
        """Parse range string to get start and end positions"""
        # Handle discontinuous ranges by taking first and last positions
        parts = range_str.split(',')

        # Get first segment start
        first_segment = parts[0].strip()
        start_pos = int(first_segment.split('-')[0])

        # Get last segment end
        last_segment = parts[-1].strip()
        end_pos = int(last_segment.split('-')[-1])

        return start_pos, end_pos

    def assess_and_filter_results(self, tier_filter: List[str] = None) -> Dict[str, List[ProteinQuality]]:
        """Assess all results and filter by quality tier (unchanged from original)"""

        if tier_filter is None:
            tier_filter = ['excellent', 'good']

        logger.info(f"🔍 Assessing quality and filtering for tiers: {tier_filter}")

        # Get all quality assessments
        all_results = self.quality_assessor.assess_all_batches()

        # Filter by tier
        filtered_by_batch = defaultdict(list)
        total_assessed = len(all_results)
        total_filtered = 0

        for result in all_results:
            if result.tier in tier_filter:
                filtered_by_batch[result.batch_name].append(result)
                total_filtered += 1

        logger.info(f"✓ Quality filtering: {total_filtered}/{total_assessed} results pass tier filter")

        return dict(filtered_by_batch)

    def import_quality_filtered_results(self, tier_filter: List[str] = None,
                                      limit: Optional[int] = None,
                                      collision_strategy: str = "separate",
                                      create_missing_batches: bool = False) -> Dict[str, any]:
        """Import only high-quality enhanced results"""

        if tier_filter is None:
            tier_filter = ['excellent', 'good']

        logger.info(f"🚀 Enhanced quality-filtered import v2.0: tiers {tier_filter}")

        # Get filtered results
        filtered_results = self.assess_and_filter_results(tier_filter)

        # Import each filtered result
        stats = {
            'total_assessed': 0,
            'quality_filtered': 0,
            'imported': 0,
            'failed': 0,
            'by_batch': {},
            'by_tier': defaultdict(int)
        }

        imported_count = 0

        for batch_name, quality_results in filtered_results.items():
            logger.info(f"Processing batch {batch_name}: {len(quality_results)} quality results")

            batch_stats = {'imported': 0, 'failed': 0, 'total': len(quality_results)}

            for quality_result in quality_results:
                if limit and imported_count >= limit:
                    break

                try:
                    # Look for enhanced XML files
                    xml_file = Path(self.config["paths"]["batch_base_dir"]) / batch_name / "mini_domains" / f"{quality_result.protein_id}.mini.domains.xml"

                    if not xml_file.exists():
                        logger.warning(f"Enhanced XML file not found: {xml_file}")
                        batch_stats['failed'] += 1
                        continue

                    # Parse enhanced XML
                    enhanced_result = self.parse_enhanced_mini_xml(xml_file)

                    # Import with collision handling and batch creation
                    if create_missing_batches:
                        success = self.import_enhanced_result_with_batch_creation(
                            enhanced_result, collision_strategy, create_missing_batches=True)
                    else:
                        success = self.import_enhanced_result(enhanced_result, collision_strategy)

                    if success:
                        batch_stats['imported'] += 1
                        stats['imported'] += 1
                        stats['by_tier'][quality_result.tier] += 1
                        imported_count += 1

                        logger.info(f"✓ Imported {quality_result.protein_id} v2.0 (tier: {quality_result.tier}, "
                                  f"domains: {len(enhanced_result.domains)})")
                    else:
                        batch_stats['failed'] += 1
                        stats['failed'] += 1

                except Exception as e:
                    logger.error(f"Error importing {quality_result.protein_id}: {e}")
                    batch_stats['failed'] += 1
                    stats['failed'] += 1

            stats['by_batch'][batch_name] = batch_stats
            logger.info(f"✓ Batch {batch_name}: {batch_stats['imported']}/{batch_stats['total']} imported")

        stats['total_assessed'] = sum(len(results) for results in filtered_results.values())
        stats['quality_filtered'] = stats['total_assessed']

        logger.info(f"🎉 Enhanced v2.0 import complete:")
        logger.info(f"  Quality filtered: {stats['quality_filtered']}")
        logger.info(f"  Successfully imported: {stats['imported']}")
        logger.info(f"  Failed: {stats['failed']}")

        return stats

    def generate_v2_analysis(self) -> Dict[str, any]:
        """Generate analysis of v2.0 imports with enhanced metrics"""

        logger.info("📊 Generating v2.0 analysis...")

        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # v2.0 specific statistics
            cursor.execute("""
                SELECT
                    COUNT(*) as protein_count,
                    COUNT(*) FILTER (WHERE is_classified = true) as classified_count,
                    AVG(domains_with_evidence) as avg_domains,
                    SUM(domains_with_evidence) as total_domains,
                    COUNT(*) FILTER (WHERE algorithm_version IS NOT NULL) as with_version_info,
                    COUNT(*) FILTER (WHERE git_commit_hash IS NOT NULL) as with_git_hash,
                    COUNT(*) FILTER (WHERE process_parameters IS NOT NULL) as with_parameters
                FROM pdb_analysis.partition_proteins
                WHERE process_version = 'mini_pyecod_v2.0'
            """)

            v2_stats = cursor.fetchone()

            # Enhanced evidence quality metrics
            cursor.execute("""
                SELECT
                    de.evidence_type,
                    COUNT(*) as evidence_count,
                    AVG(de.confidence) as avg_confidence,
                    COUNT(*) FILTER (WHERE de.evalue IS NOT NULL) as with_evalue,
                    COUNT(*) FILTER (WHERE de.probability IS NOT NULL) as with_probability
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                JOIN pdb_analysis.domain_evidence de ON pd.id = de.domain_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                GROUP BY de.evidence_type
                ORDER BY COUNT(*) DESC
            """)

            evidence_quality = cursor.fetchall()

            # Boundary optimization statistics
            cursor.execute("""
                SELECT
                    COUNT(*) as total_domains,
                    COUNT(*) FILTER (WHERE original_range IS NOT NULL) as optimized_domains,
                    COUNT(*) FILTER (WHERE optimization_actions IS NOT NULL) as with_optimization_actions
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
            """)

            optimization_stats = cursor.fetchone()

        analysis = {
            'v2_statistics': dict(v2_stats) if v2_stats else {},
            'evidence_quality_breakdown': [dict(row) for row in evidence_quality],
            'boundary_optimization': dict(optimization_stats) if optimization_stats else {},
            'analysis_timestamp': datetime.now().isoformat()
        }

        return analysis


def main():
    """Command line interface for enhanced v2.0 importer"""
    parser = argparse.ArgumentParser(
        description='Enhanced Mini PyECOD v2.0 Importer with Comprehensive Provenance'
    )

    parser.add_argument('--assess-and-import', action='store_true',
                       help='Assess quality and import enhanced v2.0 results')
    parser.add_argument('--tier-filter', type=str, default='excellent,good',
                       help='Comma-separated quality tiers to import (default: excellent,good)')
    parser.add_argument('--limit', type=int,
                       help='Maximum results to import')
    parser.add_argument('--collision-strategy', type=str, default='separate',
                       choices=['separate', 'skip'],
                       help='Collision handling strategy (default: separate)')
    parser.add_argument('--create-missing-batches', action='store_true',
                       help='Create batch records for missing batch names')

    parser.add_argument('--link-existing-batches', action='store_true',
                       help='Link existing partition records to batch table')
    parser.add_argument('--diagnose-batches', action='store_true',
                       help='Diagnose batch relationship issues')

    parser.add_argument('--comparative-analysis', action='store_true',
                       help='Generate v2.0 analysis report')

    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')

    args = parser.parse_args()

    # Initialize enhanced importer
    importer = EnhancedMiniImporter(args.config)

    if args.link_existing_batches:
        importer = EnhancedMiniImporter(args.config)

        # Use raw SQL to call the linking function
        with importer.db_conn.cursor() as cursor:
            cursor.execute("SELECT pdb_analysis.link_partitions_to_batches()")
            result = cursor.fetchone()[0]
            print(f"✅ Batch Linking Results: {result}")
        return

    if args.diagnose_batches:
        importer = EnhancedMiniImporter(args.config)

        print("\n🔍 Batch Relationship Diagnosis")
        print("=" * 50)

        with importer.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            cursor.execute("SELECT * FROM pdb_analysis.diagnose_batch_relationships()")
            results = cursor.fetchall()

            for result in results:
                print(f"\n📊 {result['issue_type'].replace('_', ' ').title()}:")
                print(f"  Count: {result['count']:,}")
                if result['example_details']:
                    print(f"  Details: {result['example_details']}")
                print(f"  Action: {result['suggested_action']}")
        return result

    if args.assess_and_import:
        tier_filter = [tier.strip() for tier in args.tier_filter.split(',')]
        stats = importer.import_quality_filtered_results(
            tier_filter=tier_filter,
            limit=args.limit,
            collision_strategy=args.collision_strategy,
            create_missing_batches=args.create_missing_batches
        )

        print(f"\n✅ Enhanced v2.0 Import Results:")
        print(f"  Process version: mini_pyecod_v2.0")
        print(f"  Tier filter: {tier_filter}")
        print(f"  Create missing batches: {args.create_missing_batches}")
        print(f"  Quality filtered: {stats['quality_filtered']}")
        print(f"  Successfully imported: {stats['imported']}")
        print(f"  Failed: {stats['failed']}")

        if stats['failed'] > 0:
            print(f"\n💡 If imports failed due to missing batches, try:")
            print(f"  --create-missing-batches    (create batch records)")
            print(f"  --diagnose-batches         (check batch issues)")
            print(f"  --link-existing-batches    (link to existing batches)")
        return
    
    if args.comparative_analysis:
        analysis = importer.generate_v2_analysis()
        
        print("\n🔍 Enhanced Mini PyECOD v2.0 Analysis Report")
        print("=" * 60)
        print(f"Generated: {analysis['analysis_timestamp']}")
        print()
        
        # v2.0 Statistics
        v2_stats = analysis['v2_statistics']
        if v2_stats:
            print("📊 v2.0 Import Statistics:")
            print(f"  Total proteins imported: {v2_stats['protein_count']:,}")
            print(f"  Classified proteins: {v2_stats['classified_count']:,}")
            print(f"  Total domains: {v2_stats['total_domains']:,}")
            print(f"  With version info: {v2_stats['with_version_info']:,}")
            print(f"  With git hash: {v2_stats['with_git_hash']:,}")
            print(f"  With parameters: {v2_stats['with_parameters']:,}")
            print()
        
        # Evidence Quality
        if analysis['evidence_quality_breakdown']:
            print("🎯 Evidence Quality Breakdown:")
            for ev in analysis['evidence_quality_breakdown']:
                print(f"  {ev['evidence_type']:<20} {ev['evidence_count']:>6} domains "
                      f"(conf: {ev['avg_confidence']:.3f})")
            print()
        
        # Optimization
        opt_stats = analysis['boundary_optimization']
        if opt_stats:
            opt_rate = opt_stats['optimized_domains'] / max(1, opt_stats['total_domains']) * 100
            print("🔧 Boundary Optimization:")
            print(f"  Total domains: {opt_stats['total_domains']:,}")
            print(f"  Optimized domains: {opt_stats['optimized_domains']:,} ({opt_rate:.1f}%)")
            print(f"  With optimization actions: {opt_stats['with_optimization_actions']:,}")
        
        return
    
    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
