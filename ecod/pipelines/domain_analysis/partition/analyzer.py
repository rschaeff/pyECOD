# ecod/pipelines/domain_analysis/partition/analyzer.py
"""
Comprehensive Evidence Analyzer for Domain Analysis Pipeline

This analyzer provides advanced evidence processing, classification integration,
quality control, and performance optimization for domain analysis workflows.
"""

import logging
import os
import json
import time
import threading
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Set, Tuple, Union, Iterator
from pathlib import Path
import xml.etree.ElementTree as ET
import re
from datetime import datetime, timedelta
import math
import statistics

from ecod.models.pipeline.evidence import Evidence
from ecod.models.pipeline.domain import DomainModel
from ecod.pipelines.domain_analysis.partition.models import (
    PartitionOptions, ValidationResult, ValidationLevel, EvidenceGroup, ProcessingStats
)


@dataclass
class EvidenceQualityMetrics:
    """Quality metrics for evidence analysis"""
    total_evidence_count: int = 0
    valid_evidence_count: int = 0
    blast_hit_count: int = 0
    hhsearch_hit_count: int = 0
    self_comparison_count: int = 0

    # Coverage metrics
    sequence_coverage: float = 0.0
    gap_count: int = 0
    largest_gap_size: int = 0
    overlap_count: int = 0
    largest_overlap_size: int = 0

    # Quality scores
    average_confidence: float = 0.0
    confidence_std: float = 0.0
    high_confidence_count: int = 0  # > 0.8
    low_confidence_count: int = 0   # < 0.3

    # Classification metrics
    classified_evidence_count: int = 0
    classification_consistency: float = 0.0  # How consistent are classifications

    # Performance metrics
    processing_time_seconds: float = 0.0
    memory_usage_mb: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'total_evidence': self.total_evidence_count,
            'valid_evidence': self.valid_evidence_count,
            'blast_hits': self.blast_hit_count,
            'hhsearch_hits': self.hhsearch_hit_count,
            'self_comparison': self.self_comparison_count,
            'sequence_coverage': self.sequence_coverage,
            'gap_count': self.gap_count,
            'largest_gap': self.largest_gap_size,
            'overlap_count': self.overlap_count,
            'largest_overlap': self.largest_overlap_size,
            'average_confidence': self.average_confidence,
            'confidence_std': self.confidence_std,
            'high_confidence_count': self.high_confidence_count,
            'low_confidence_count': self.low_confidence_count,
            'classified_evidence': self.classified_evidence_count,
            'classification_consistency': self.classification_consistency,
            'processing_time': self.processing_time_seconds,
            'memory_usage_mb': self.memory_usage_mb
        }


@dataclass
class ClassificationCache:
    """Advanced caching system for domain classifications"""

    def __init__(self, max_size: int = 10000, ttl_hours: int = 24):
        self.max_size = max_size
        self.ttl = timedelta(hours=ttl_hours)
        self._cache: Dict[str, Dict[str, Any]] = {}
        self._timestamps: Dict[str, datetime] = {}
        self._access_count: Dict[str, int] = defaultdict(int)
        self._lock = threading.RLock()
        self.logger = logging.getLogger(__name__)

        # Statistics
        self.hits = 0
        self.misses = 0
        self.evictions = 0

    def get_domain_classification(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """Get cached classification with TTL and LRU eviction"""
        with self._lock:
            if domain_id not in self._cache:
                self.misses += 1
                return None

            # Check TTL
            if datetime.now() - self._timestamps[domain_id] > self.ttl:
                self._evict(domain_id)
                self.misses += 1
                return None

            # Update access tracking
            self._access_count[domain_id] += 1
            self.hits += 1
            return self._cache[domain_id].copy()

    def set_domain_classification(self, domain_id: str, classification: Dict[str, Any]) -> bool:
        """Set classification with automatic eviction"""
        with self._lock:
            try:
                # Evict if at capacity
                if len(self._cache) >= self.max_size and domain_id not in self._cache:
                    self._evict_lru()

                self._cache[domain_id] = classification.copy()
                self._timestamps[domain_id] = datetime.now()
                self._access_count[domain_id] = 1
                return True

            except Exception as e:
                self.logger.warning(f"Error setting cache for {domain_id}: {str(e)}")
                return False

    def _evict(self, domain_id: str) -> None:
        """Evict specific entry"""
        self._cache.pop(domain_id, None)
        self._timestamps.pop(domain_id, None)
        self._access_count.pop(domain_id, None)
        self.evictions += 1

    def _evict_lru(self) -> None:
        """Evict least recently used entry"""
        if not self._access_count:
            return

        lru_domain = min(self._access_count.keys(), key=self._access_count.get)
        self._evict(lru_domain)

    def cleanup_expired(self) -> int:
        """Clean up expired entries"""
        with self._lock:
            now = datetime.now()
            expired = [
                domain_id for domain_id, timestamp in self._timestamps.items()
                if now - timestamp > self.ttl
            ]

            for domain_id in expired:
                self._evict(domain_id)

            return len(expired)

    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        total_requests = self.hits + self.misses
        hit_rate = (self.hits / total_requests * 100) if total_requests > 0 else 0

        return {
            'size': len(self._cache),
            'max_size': self.max_size,
            'hits': self.hits,
            'misses': self.misses,
            'hit_rate_percent': hit_rate,
            'evictions': self.evictions,
            'ttl_hours': self.ttl.total_seconds() / 3600
        }


class EvidenceAnalyzer:
    """
    Comprehensive evidence analyzer for domain analysis pipeline.

    Provides advanced evidence processing, classification integration,
    quality control, and performance optimization.
    """

    def __init__(self, options: PartitionOptions, db_manager=None):
        self.options = options
        self.db_manager = db_manager
        self.logger = logging.getLogger(__name__)

        # Initialize classification cache
        if options.use_cache:
            cache_size = getattr(options, 'cache_size', 10000)
            cache_ttl = getattr(options, 'cache_ttl_hours', 24)
            self.classification_cache = ClassificationCache(cache_size, cache_ttl)
        else:
            self.classification_cache = None

        # Performance tracking
        self.stats = ProcessingStats()
        self._processing_times = []
        self._memory_usage = []

        # Evidence processing patterns
        self.range_patterns = [
            re.compile(r'(\d+)-(\d+)'),  # Standard range: 123-456
            re.compile(r'(\d+):(\d+)'),  # Colon range: 123:456
            re.compile(r'(\d+)\.\.(\d+)'),  # Dot range: 123..456
        ]

        # Classification integration
        self.ecod_hierarchy = self._load_ecod_hierarchy()
        self.representative_domains = self._load_representative_domains()

        # Thread pool for parallel processing
        if options.parallel_processing:
            self.executor = ThreadPoolExecutor(max_workers=options.max_workers)
        else:
            self.executor = None

        # Quality thresholds
        self.quality_thresholds = {
            'min_sequence_coverage': 0.7,
            'max_gap_size': 50,
            'min_evidence_confidence': options.min_evidence_confidence,
            'max_classification_conflicts': 3,
            'min_evidence_density': 0.1  # evidence per 100 residues
        }

        self.logger.info(f"EvidenceAnalyzer initialized with {options.validation_level.value} validation")

    def analyze_domain_summary(self, file_path: str, protein_id: str = None,
                             sequence_length: int = 0) -> Dict[str, Any]:
        """
        Comprehensive analysis of domain summary file.

        Args:
            file_path: Path to domain summary XML file
            protein_id: Protein identifier for context
            sequence_length: Protein sequence length for coverage analysis

        Returns:
            dict: Comprehensive analysis results
        """
        start_time = time.time()

        try:
            self.logger.debug(f"Starting comprehensive analysis: {file_path}")

            # Parse domain summary file
            summary_data = self.parse_domain_summary(file_path)

            if 'error' in summary_data:
                return {
                    'success': False,
                    'error': summary_data['error'],
                    'file_path': file_path,
                    'protein_id': protein_id
                }

            # Extract and process evidence
            evidence_list = self.extract_evidence_with_classification(summary_data)

            # Validate evidence
            validation_results = self.validate_evidence_list(evidence_list, protein_id)

            # Filter evidence based on validation and options
            filtered_evidence = self.filter_evidence(evidence_list, validation_results)

            # Group evidence by position and type
            evidence_groups = self.group_evidence_comprehensive(filtered_evidence, sequence_length)

            # Resolve conflicts between evidence sources
            resolved_evidence = self.resolve_evidence_conflicts(evidence_groups)

            # Calculate quality metrics
            quality_metrics = self.calculate_quality_metrics(
                filtered_evidence, sequence_length, time.time() - start_time
            )

            # Generate domain suggestions
            domain_suggestions = self.generate_domain_suggestions(resolved_evidence, sequence_length)

            # Classification analysis
            classification_analysis = self.analyze_classifications(filtered_evidence)

            processing_time = time.time() - start_time
            self._processing_times.append(processing_time)

            return {
                'success': True,
                'file_path': file_path,
                'protein_id': protein_id,
                'sequence_length': sequence_length,
                'evidence_count': len(evidence_list),
                'filtered_evidence_count': len(filtered_evidence),
                'evidence_groups': evidence_groups,
                'resolved_evidence': resolved_evidence,
                'domain_suggestions': domain_suggestions,
                'quality_metrics': quality_metrics.to_dict(),
                'classification_analysis': classification_analysis,
                'validation_summary': self._summarize_validation_results(validation_results),
                'processing_time_seconds': processing_time,
                'cache_stats': self.classification_cache.get_stats() if self.classification_cache else None
            }

        except Exception as e:
            error_msg = f"Error in comprehensive analysis: {str(e)}"
            self.logger.error(error_msg)
            return {
                'success': False,
                'error': error_msg,
                'file_path': file_path,
                'protein_id': protein_id,
                'processing_time_seconds': time.time() - start_time
            }

    def parse_domain_summary(self, file_path: str) -> Dict[str, Any]:
        """
        Enhanced domain summary parsing with comprehensive error handling.

        Returns dict with either parsed data or error information.
        """
        try:
            # Validate file existence and permissions
            if not os.path.exists(file_path):
                return {"error": f"Domain summary file not found: {file_path}"}

            if not os.access(file_path, os.R_OK):
                return {"error": f"Permission denied reading file: {file_path}"}

            # Check file size (warn for very large files)
            file_size = os.path.getsize(file_path)
            if file_size > 50 * 1024 * 1024:  # 50MB
                self.logger.warning(f"Large domain summary file: {file_size / 1024 / 1024:.1f}MB")

            # Parse XML with multiple fallback strategies
            try:
                tree = ET.parse(file_path)
                root = tree.getroot()
            except ET.ParseError as e:
                # Try to repair common XML issues
                return self._parse_with_repair(file_path, str(e))
            except UnicodeDecodeError:
                # Try different encodings
                return self._parse_with_encoding_fallback(file_path)

            # Validate XML structure
            if root.tag != "blast_summ_doc":
                return {"error": f"Invalid XML structure: expected 'blast_summ_doc' root, got '{root.tag}'"}

            # Extract metadata
            metadata = self._extract_metadata(root)

            # Parse evidence sections with parallel processing if enabled
            if self.executor:
                evidence_results = self._parse_evidence_parallel(root)
            else:
                evidence_results = self._parse_evidence_sequential(root)

            # Combine results
            result = {
                **metadata,
                **evidence_results,
                'file_size_bytes': file_size,
                'parsing_method': 'parallel' if self.executor else 'sequential'
            }

            self.logger.debug(f"Successfully parsed domain summary: {file_path}")
            return result

        except Exception as e:
            error_msg = f"Unexpected error parsing {file_path}: {str(e)}"
            self.logger.error(error_msg)
            return {"error": error_msg}

    # Replace the _parse_with_repair method in ecod/pipelines/domain_analysis/partition/analyzer.py

    def _parse_with_repair(self, file_path: str, parse_error: str) -> Dict[str, Any]:
        """Attempt to repair and parse XML with common issues - FIXED to be less aggressive"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # Check for severe structural problems that can't be repaired
            severe_issues = [
                r'<[^/>]*$',  # Unclosed tag at end
                r'<[^>]*></[^>]*>',  # Mismatched tags
                r'<[^>]*>[^<]*</[^>]+[^>]*$',  # Malformed closing tag
            ]

            for pattern in severe_issues:
                if re.search(pattern, content):
                    return {"error": f"Severely malformed XML cannot be repaired: {parse_error}"}

            # Only attempt simple repairs
            simple_repairs = [
                (r'&(?!amp;|lt;|gt;|quot;|apos;)', '&amp;'),  # Fix unescaped ampersands
            ]

            repaired_content = content
            repairs_made = 0

            for pattern, replacement in simple_repairs:
                new_content = re.sub(pattern, replacement, repaired_content)
                if new_content != repaired_content:
                    repairs_made += 1
                    repaired_content = new_content

            # If no repairs were possible, return error
            if repairs_made == 0:
                return {"error": f"XML parsing failed and no repairs possible: {parse_error}"}

            # Try parsing repaired content
            root = ET.fromstring(repaired_content)

            self.logger.warning(f"Successfully repaired XML with {repairs_made} fixes: {file_path}")
            return self._extract_all_data(root)

        except Exception as e:
            return {"error": f"XML parsing failed even after repair attempts: {str(e)}"}

    def _parse_with_encoding_fallback(self, file_path: str) -> Dict[str, Any]:
        """Try parsing with different encodings"""
        encodings = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']

        for encoding in encodings:
            try:
                with open(file_path, 'r', encoding=encoding) as f:
                    tree = ET.parse(f)
                    root = tree.getroot()

                self.logger.warning(f"Successfully parsed with {encoding} encoding: {file_path}")
                return self._extract_all_data(root)

            except (UnicodeDecodeError, ET.ParseError):
                continue

        return {"error": f"Could not parse XML with any encoding: {file_path}"}

    def _extract_metadata(self, root: ET.Element) -> Dict[str, Any]:
        """Extract metadata from XML root"""
        metadata = {}

        # Basic identification
        blast_summ = root.find("blast_summ")
        if blast_summ is not None:
            metadata['pdb_id'] = blast_summ.get("pdb", "")
            metadata['chain_id'] = blast_summ.get("chain", "")

        # Additional metadata
        metadata['xml_structure_version'] = root.get("version", "unknown")
        metadata['creation_date'] = root.get("created", "")

        return metadata

    def _parse_evidence_parallel(self, root: ET.Element) -> Dict[str, Any]:
        """Parse evidence sections using parallel processing"""
        futures = []

        # Submit parsing tasks
        futures.append(self.executor.submit(self._parse_blast_hits, root))
        futures.append(self.executor.submit(self._parse_hhsearch_hits, root))
        futures.append(self.executor.submit(self._parse_domain_suggestions, root))
        futures.append(self.executor.submit(self._parse_self_comparison, root))

        # Collect results
        results = {}
        for future in as_completed(futures):
            try:
                section_results = future.result(timeout=30)
                results.update(section_results)
            except Exception as e:
                self.logger.warning(f"Error in parallel evidence parsing: {str(e)}")

        return results

    def _parse_evidence_sequential(self, root: ET.Element) -> Dict[str, Any]:
        """Parse evidence sections sequentially"""
        results = {}

        try:
            results.update(self._parse_blast_hits(root))
        except Exception as e:
            self.logger.warning(f"Error parsing BLAST hits: {str(e)}")
            results['blast_hits'] = []

        try:
            results.update(self._parse_hhsearch_hits(root))
        except Exception as e:
            self.logger.warning(f"Error parsing HHSearch hits: {str(e)}")
            results['hhsearch_hits'] = []

        try:
            results.update(self._parse_domain_suggestions(root))
        except Exception as e:
            self.logger.warning(f"Error parsing domain suggestions: {str(e)}")
            results['domain_suggestions'] = []

        try:
            results.update(self._parse_self_comparison(root))
        except Exception as e:
            self.logger.warning(f"Error parsing self comparison: {str(e)}")
            results['self_comparison'] = []

        return results

    def _extract_all_data(self, root: ET.Element) -> Dict[str, Any]:
        """Extract all data from parsed XML root"""
        metadata = self._extract_metadata(root)
        evidence = self._parse_evidence_sequential(root)
        return {**metadata, **evidence}

    def _parse_blast_hits(self, root: ET.Element) -> Dict[str, Any]:
        """Parse BLAST hits with enhanced processing"""
        hits = []

        # Chain BLAST hits
        chain_blast = root.find("chain_blast_run")
        if chain_blast is not None:
            hits_elem = chain_blast.find("hits")
            if hits_elem is not None:
                for hit_elem in hits_elem.findall("hit"):
                    hit_data = self._parse_blast_hit_element(hit_elem, "chain_blast")
                    if hit_data:
                        hits.append(hit_data)

        # Domain BLAST hits
        domain_blast = root.find("blast_run")
        if domain_blast is not None:
            hits_elem = domain_blast.find("hits")
            if hits_elem is not None:
                for hit_elem in hits_elem.findall("hit"):
                    hit_data = self._parse_blast_hit_element(hit_elem, "domain_blast")
                    if hit_data:
                        hits.append(hit_data)

        return {'blast_hits': hits}

    def _parse_blast_hit_element(self, hit_elem: ET.Element, hit_type: str) -> Optional[Dict[str, Any]]:
        """Parse individual BLAST hit with comprehensive data extraction"""
        try:
            hit_data = {
                "type": hit_type,
                "hit_id": hit_elem.get("num", ""),
                "domain_id": hit_elem.get("domain_id", ""),
                "pdb_id": hit_elem.get("pdb_id", ""),
                "chain_id": hit_elem.get("chain_id", ""),
            }

            # Parse numeric values with robust error handling
            evalues = hit_elem.get("evalues", "")
            if evalues:
                try:
                    # Handle multiple e-values (comma-separated)
                    evalue_list = [float(e.strip()) for e in evalues.split(",") if e.strip()]
                    hit_data["evalue"] = min(evalue_list) if evalue_list else 999.0
                    hit_data["evalue_list"] = evalue_list
                except (ValueError, IndexError):
                    hit_data["evalue"] = 999.0

            # Additional BLAST-specific fields
            for field in ["hsp_count", "identity", "coverage", "score"]:
                value = hit_elem.get(field, "")
                if value:
                    try:
                        if field == "hsp_count":
                            hit_data[field] = int(value)
                        else:
                            hit_data[field] = float(value)
                    except ValueError:
                        pass

            # Extract ranges with multiple format support
            hit_data["query_range"] = self._extract_range_text(hit_elem, "query_reg", "query_range")
            hit_data["hit_range"] = self._extract_range_text(hit_elem, "hit_reg", "hit_range")

            # Parse ranges into structured format
            hit_data["query_range_parsed"] = self._parse_range_comprehensive(hit_data["query_range"])
            hit_data["hit_range_parsed"] = self._parse_range_comprehensive(hit_data["hit_range"])

            # Additional metadata
            hit_data["discontinuous"] = hit_elem.get("discontinuous", "false").lower() == "true"
            hit_data["aligned_length"] = self._calculate_aligned_length(hit_data["query_range_parsed"])

            return hit_data

        except Exception as e:
            self.logger.warning(f"Error parsing BLAST hit element: {str(e)}")
            return None

    def _parse_hhsearch_hits(self, root: ET.Element) -> Dict[str, Any]:
        """Parse HHSearch hits with comprehensive processing"""
        hits = []

        hh_run = root.find("hh_run")
        if hh_run is not None:
            hits_elem = hh_run.find("hits")
            if hits_elem is not None:
                for hit_elem in hits_elem.findall("hit"):
                    hit_data = self._parse_hhsearch_hit_element(hit_elem)
                    if hit_data:
                        hits.append(hit_data)

        return {'hhsearch_hits': hits}

    def _parse_hhsearch_hit_element(self, hit_elem: ET.Element) -> Optional[Dict[str, Any]]:
        """Parse individual HHSearch hit with comprehensive data extraction"""
        try:
            hit_data = {
                "type": "hhsearch",
                "hit_id": hit_elem.get("hit_id", ""),
                "domain_id": hit_elem.get("domain_id", ""),
                "source_id": hit_elem.get("source_id", ""),
            }

            # Parse numeric values with validation
            for attr in ["probability", "evalue", "score", "ss_score", "cols", "query_length", "template_length"]:
                value_str = hit_elem.get(attr, "")
                if value_str:
                    try:
                        value = float(value_str)
                        if math.isfinite(value):  # Check for NaN/inf
                            hit_data[attr] = value
                    except ValueError:
                        pass

            # Extract alignment information
            hit_data["query_range"] = self._extract_range_text(hit_elem, "query_reg", "query_range")
            hit_data["hit_range"] = self._extract_range_text(hit_elem, "hit_reg", "hit_range")

            # Parse ranges
            hit_data["query_range_parsed"] = self._parse_range_comprehensive(hit_data["query_range"])
            hit_data["hit_range_parsed"] = self._parse_range_comprehensive(hit_data["hit_range"])

            # Extract additional HHSearch-specific data
            ss_elem = hit_elem.find("secondary_structure")
            if ss_elem is not None:
                hit_data["secondary_structure"] = ss_elem.text

            # Calculate coverage metrics
            if "query_length" in hit_data and hit_data["query_range_parsed"]:
                hit_data["query_coverage"] = self._calculate_coverage(
                    hit_data["query_range_parsed"], hit_data["query_length"]
                )

            return hit_data

        except Exception as e:
            self.logger.warning(f"Error parsing HHSearch hit element: {str(e)}")
            return None

    def _parse_domain_suggestions(self, root: ET.Element) -> Dict[str, Any]:
        """Parse domain suggestions if present"""
        suggestions = []

        domain_suggestions = root.find("domain_suggestions")
        if domain_suggestions is not None:
            for domain_elem in domain_suggestions.findall("domain"):
                try:
                    suggestion = {
                        "id": domain_elem.get("id", ""),
                        "range": domain_elem.get("range", ""),
                        "confidence": float(domain_elem.get("confidence", "0.0")),
                        "source": domain_elem.get("source", "predicted")
                    }

                    # Parse range
                    suggestion["range_parsed"] = self._parse_range_comprehensive(suggestion["range"])

                    suggestions.append(suggestion)

                except Exception as e:
                    self.logger.warning(f"Error parsing domain suggestion: {str(e)}")

        return {'domain_suggestions': suggestions}

    def _parse_self_comparison(self, root: ET.Element) -> Dict[str, Any]:
        """Parse self-comparison results if present"""
        self_comp = []

        self_comp_elem = root.find("self_comparison")
        if self_comp_elem is not None:
            for hit_elem in self_comp_elem.findall("hit"):
                try:
                    hit_data = {
                        "type": "self_comparison",
                        "range1": self._extract_range_text(hit_elem, "range1"),
                        "range2": self._extract_range_text(hit_elem, "range2"),
                        "score": float(hit_elem.get("score", "0.0")),
                        "identity": float(hit_elem.get("identity", "0.0"))
                    }

                    self_comp.append(hit_data)

                except Exception as e:
                    self.logger.warning(f"Error parsing self-comparison hit: {str(e)}")

        return {'self_comparison': self_comp}

    def _extract_range_text(self, element: ET.Element, *possible_tags: str) -> str:
        """Extract range text from multiple possible element names"""
        for tag in possible_tags:
            # Try as child element
            child = element.find(tag)
            if child is not None and child.text:
                return child.text.strip()

            # Try as attribute
            attr_value = element.get(tag, "")
            if attr_value:
                return attr_value

        return ""

    def _parse_range_comprehensive(self, range_str: str) -> List[Tuple[int, int]]:
        """Parse range string into list of (start, end) tuples with multiple format support"""
        if not range_str:
            return []

        ranges = []

        # Split by comma for multiple ranges
        for segment in range_str.split(","):
            segment = segment.strip()
            if not segment:
                continue

            # Try different range patterns
            for pattern in self.range_patterns:
                match = pattern.search(segment)
                if match:
                    try:
                        start, end = int(match.group(1)), int(match.group(2))
                        if start <= end and start > 0:  # Validate range
                            ranges.append((start, end))
                        break
                    except ValueError:
                        continue
            else:
                # Try single number
                try:
                    pos = int(segment)
                    if pos > 0:
                        ranges.append((pos, pos))
                except ValueError:
                    pass

        return ranges

    def _calculate_aligned_length(self, ranges: List[Tuple[int, int]]) -> int:
        """Calculate total aligned length from ranges"""
        return sum(end - start + 1 for start, end in ranges)

    def _calculate_coverage(self, ranges: List[Tuple[int, int]], total_length: int) -> float:
        """Calculate coverage percentage"""
        if total_length <= 0:
            return 0.0

        aligned_length = self._calculate_aligned_length(ranges)
        return (aligned_length / total_length) * 100.0

    def extract_evidence_with_classification(self, summary_data: Dict[str, Any]) -> List[Evidence]:
        """Extract Evidence objects with integrated classification lookup"""
        evidence_list = []

        try:
            # Process BLAST hits
            for hit_data in summary_data.get("blast_hits", []):
                evidence = self._create_evidence_from_blast(hit_data)
                if evidence:
                    # Enhance with classification
                    self._enhance_evidence_with_classification(evidence)
                    evidence_list.append(evidence)

            # Process HHSearch hits
            for hit_data in summary_data.get("hhsearch_hits", []):
                evidence = self._create_evidence_from_hhsearch(hit_data)
                if evidence:
                    # Enhance with classification
                    self._enhance_evidence_with_classification(evidence)
                    evidence_list.append(evidence)

            # Process self-comparison results
            for hit_data in summary_data.get("self_comparison", []):
                evidence = self._create_evidence_from_self_comparison(hit_data)
                if evidence:
                    evidence_list.append(evidence)

            self.logger.debug(f"Extracted {len(evidence_list)} evidence items")

        except Exception as e:
            self.logger.error(f"Error extracting evidence: {str(e)}")

        return evidence_list

    def _create_evidence_from_blast(self, hit_data: Dict[str, Any]) -> Optional[Evidence]:
        """Create Evidence object from BLAST hit data"""
        try:
            return Evidence(
                type=hit_data.get("type", "blast"),
                source_id=hit_data.get("domain_id", "") or hit_data.get("hit_id", ""),
                domain_id=hit_data.get("domain_id", ""),
                query_range=hit_data.get("query_range", ""),
                hit_range=hit_data.get("hit_range", ""),
                evalue=hit_data.get("evalue"),
                hsp_count=hit_data.get("hsp_count"),
                identity=hit_data.get("identity"),
                coverage=hit_data.get("coverage"),
                confidence=None,  # Will be auto-calculated
                extra_attributes={
                    "pdb_id": hit_data.get("pdb_id", ""),
                    "chain_id": hit_data.get("chain_id", ""),
                    "discontinuous": hit_data.get("discontinuous", False),
                    "aligned_length": hit_data.get("aligned_length", 0)
                }
            )
        except Exception as e:
            self.logger.warning(f"Error creating BLAST evidence: {str(e)}")
            return None

    def _create_evidence_from_hhsearch(self, hit_data: Dict[str, Any]) -> Optional[Evidence]:
        """Create Evidence object from HHSearch hit data"""
        try:
            return Evidence(
                type="hhsearch",
                source_id=hit_data.get("source_id", "") or hit_data.get("domain_id", ""),
                domain_id=hit_data.get("domain_id", ""),
                query_range=hit_data.get("query_range", ""),
                hit_range=hit_data.get("hit_range", ""),
                probability=hit_data.get("probability"),
                evalue=hit_data.get("evalue"),
                score=hit_data.get("score"),
                confidence=None,  # Will be auto-calculated
                extra_attributes={
                    "hit_id": hit_data.get("hit_id", ""),
                    "ss_score": hit_data.get("ss_score"),
                    "cols": hit_data.get("cols"),
                    "query_coverage": hit_data.get("query_coverage"),
                    "secondary_structure": hit_data.get("secondary_structure", "")
                }
            )
        except Exception as e:
            self.logger.warning(f"Error creating HHSearch evidence: {str(e)}")
            return None

    def _create_evidence_from_self_comparison(self, hit_data: Dict[str, Any]) -> Optional[Evidence]:
        """Create Evidence object from self-comparison data"""
        try:
            return Evidence(
                type="self_comparison",
                source_id="self",
                query_range=hit_data.get("range1", ""),
                hit_range=hit_data.get("range2", ""),
                score=hit_data.get("score"),
                identity=hit_data.get("identity"),
                confidence=None,  # Will be auto-calculated
                extra_attributes={
                    "comparison_type": "internal_repeat"
                }
            )
        except Exception as e:
            self.logger.warning(f"Error creating self-comparison evidence: {str(e)}")
            return None

    def _enhance_evidence_with_classification(self, evidence: Evidence) -> None:
        """Enhance evidence with classification information"""
        if not evidence.domain_id:
            return

        try:
            # Try cache first
            classification = None
            if self.classification_cache:
                classification = self.classification_cache.get_domain_classification(evidence.domain_id)

            # Look up in database if not cached
            if not classification and self.db_manager:
                classification = self._lookup_domain_classification(evidence.domain_id)

                # Cache the result
                if classification and self.classification_cache:
                    self.classification_cache.set_domain_classification(evidence.domain_id, classification)

            # Apply classification to evidence
            if classification:
                evidence.t_group = classification.get("t_group")
                evidence.h_group = classification.get("h_group")
                evidence.x_group = classification.get("x_group")
                evidence.a_group = classification.get("a_group")

                # Add representative flags
                evidence.extra_attributes.update({
                    "is_manual_rep": classification.get("is_manual_rep", False),
                    "is_f70": classification.get("is_f70", False),
                    "is_f40": classification.get("is_f40", False),
                    "is_f99": classification.get("is_f99", False)
                })

        except Exception as e:
            self.logger.warning(f"Error enhancing evidence with classification: {str(e)}")

    def _lookup_domain_classification(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """Look up domain classification in database"""
        if not self.db_manager:
            return None

        try:
            query = """
                SELECT t_group, h_group, x_group, a_group,
                       is_manual_rep, is_f70, is_f40, is_f99
                FROM domain_classification
                WHERE domain_id = %s
            """

            if hasattr(self.db_manager, 'execute_dict_query'):
                results = self.db_manager.execute_dict_query(query, (domain_id,))
                return results[0] if results else None
            else:
                cursor = self.db_manager.execute(query, (domain_id,))
                row = cursor.fetchone()
                if row:
                    columns = [desc[0] for desc in cursor.description]
                    return dict(zip(columns, row))
                return None

        except Exception as e:
            self.logger.warning(f"Error looking up classification for {domain_id}: {str(e)}")
            return None

    def validate_evidence_list(self, evidence_list: List[Evidence], context: str = "") -> List[ValidationResult]:
        """Validate a list of evidence with comprehensive checking"""
        validation_results = []

        for i, evidence in enumerate(evidence_list):
            context_str = f"{context}_evidence_{i}" if context else f"evidence_{i}"
            result = self.validate_evidence(evidence, context_str)
            validation_results.append(result)

        return validation_results

# Replace the validate_evidence method in ecod/pipelines/domain_analysis/partition/analyzer.py

    def validate_evidence(self, evidence: Evidence, context: str) -> ValidationResult:
        """Comprehensive evidence validation - FIXED to check domain_id"""
        try:
            errors = []
            warnings = []

            # Basic validation
            if not evidence.type or evidence.type.strip() == "":
                errors.append("Evidence type is empty")

            # FIXED: Add domain_id validation
            if evidence.type in ("domain_blast", "hhsearch") and not evidence.domain_id:
                errors.append("Domain ID is empty for domain-level evidence")

            # Type-specific validation
            if evidence.type == "hhsearch":
                if evidence.probability is None and evidence.evalue is None:
                    errors.append("HHSearch evidence missing both probability and e-value")

                if evidence.probability is not None:
                    if not (0 <= evidence.probability <= 100):
                        warnings.append(f"HHSearch probability out of expected range: {evidence.probability}")

            elif evidence.type in ("blast", "domain_blast", "chain_blast"):
                if evidence.evalue is None:
                    warnings.append("BLAST evidence missing e-value")
                elif evidence.evalue < 0:
                    errors.append(f"Invalid negative e-value: {evidence.evalue}")

                if evidence.identity is not None and not (0 <= evidence.identity <= 100):
                    warnings.append(f"Identity percentage out of range: {evidence.identity}")

            # Range validation
            if evidence.query_range:
                if not self._is_valid_range_format(evidence.query_range):
                    warnings.append(f"Invalid query range format: {evidence.query_range}")

            # Confidence validation
            if evidence.confidence is not None:
                if not (0.0 <= evidence.confidence <= 1.0):
                    errors.append(f"Confidence out of range [0,1]: {evidence.confidence}")

            # Classification validation
            if evidence.t_group and not self._is_valid_t_group(evidence.t_group):
                warnings.append(f"Unrecognized T-group format: {evidence.t_group}")

            # Determine validation level
            level = self.options.validation_level
            is_valid = len(errors) == 0

            # In lenient mode, treat some errors as warnings
            if level == ValidationLevel.LENIENT and is_valid == False:
                lenient_errors = []
                for error in errors:
                    if "empty" in error.lower() or "missing" in error.lower():
                        warnings.append(f"Warning: {error}")
                    else:
                        lenient_errors.append(error)

                errors = lenient_errors
                is_valid = len(errors) == 0

            return ValidationResult(
                is_valid=is_valid,
                errors=errors,
                warnings=warnings,
                context=context,
                validation_level=level
            )

        except Exception as e:
            return ValidationResult(
                is_valid=False,
                errors=[f"Validation error: {str(e)}"],
                warnings=[],
                context=context,
                validation_level=self.options.validation_level
            )

    def _is_valid_range_format(self, range_str: str) -> bool:
        """Validate range string format"""
        try:
            ranges = self._parse_range_comprehensive(range_str)
            return len(ranges) > 0
        except:
            return False

    def _is_valid_t_group(self, t_group: str) -> bool:
        """Validate T-group format (simplified check)"""
        # Basic pattern: number.number.number
        return bool(re.match(r'^\d+\.\d+\.\d+$', t_group))

    def filter_evidence(self, evidence_list: List[Evidence],
                       validation_results: List[ValidationResult]) -> List[Evidence]:
        """Filter evidence based on validation results and quality criteria"""
        filtered_evidence = []

        for evidence, validation in zip(evidence_list, validation_results):
            # Skip invalid evidence in strict mode
            if self.options.validation_level == ValidationLevel.STRICT and not validation.is_valid:
                continue

            # Apply confidence filter
            if (evidence.confidence is not None and
                evidence.confidence < self.options.min_evidence_confidence):
                continue

            # Apply type-specific filters
            if evidence.type == "hhsearch":
                # Filter very low probability hits
                if evidence.probability is not None and evidence.probability < 10.0:
                    continue

            elif evidence.type in ("blast", "domain_blast", "chain_blast"):
                # Filter very high e-values
                if evidence.evalue is not None and evidence.evalue > 10.0:
                    continue

                # Filter very low identity
                if evidence.identity is not None and evidence.identity < 20.0:
                    continue

            filtered_evidence.append(evidence)

        self.logger.debug(f"Filtered {len(evidence_list)} -> {len(filtered_evidence)} evidence items")
        return filtered_evidence

    def group_evidence_comprehensive(self, evidence_list: List[Evidence],
                                   sequence_length: int = 0) -> List[EvidenceGroup]:
        """Group evidence by position and type with advanced clustering"""
        position_groups = defaultdict(list)

        # Group by overlapping positions
        for evidence in evidence_list:
            ranges = self._parse_range_comprehensive(evidence.query_range)
            for start, end in ranges:
                # Use range midpoint as group key
                midpoint = (start + end) // 2
                window = max(50, (end - start) * 2)  # Adaptive window size

                group_key = midpoint // window
                position_groups[group_key].append(evidence)

        # Create EvidenceGroup objects
        evidence_groups = []
        for group_id, group_evidence in position_groups.items():
            group = EvidenceGroup(
                group_id=f"group_{group_id}",
                evidence_items=group_evidence
            )

            # Update group statistics
            group._update_group_stats()
            evidence_groups.append(group)

        # Merge nearby groups if they have similar evidence
        merged_groups = self._merge_similar_groups(evidence_groups)

        return merged_groups

    def _merge_similar_groups(self, groups: List[EvidenceGroup]) -> List[EvidenceGroup]:
        """Merge evidence groups that are similar and nearby"""
        # Simple implementation - could be enhanced with sophisticated clustering
        merged = []
        used = set()

        for i, group1 in enumerate(groups):
            if i in used:
                continue

            current_group = group1
            used.add(i)

            for j, group2 in enumerate(groups[i+1:], i+1):
                if j in used:
                    continue

                # Check if groups should be merged
                if self._should_merge_groups(current_group, group2):
                    # Merge group2 into current_group
                    current_group.evidence_items.extend(group2.evidence_items)
                    current_group._update_group_stats()
                    used.add(j)

            merged.append(current_group)

        return merged

    def _should_merge_groups(self, group1: EvidenceGroup, group2: EvidenceGroup) -> bool:
        """Determine if two evidence groups should be merged"""
        # Simple criteria - could be enhanced

        # Check if groups have similar confidence
        confidence_diff = abs(group1.consensus_confidence - group2.consensus_confidence)
        if confidence_diff > 0.3:
            return False

        # Check if groups have similar types
        type1 = group1.group_type
        type2 = group2.group_type

        if type1 != type2 and type1 != "mixed" and type2 != "mixed":
            return False

        return True

    def resolve_evidence_conflicts(self, evidence_groups: List[EvidenceGroup]) -> List[Evidence]:
        """Resolve conflicts between different evidence sources"""
        resolved_evidence = []

        for group in evidence_groups:
            if len(group.evidence_items) == 1:
                # No conflicts
                resolved_evidence.extend(group.evidence_items)
                continue

            # Group by evidence type
            by_type = defaultdict(list)
            for evidence in group.evidence_items:
                by_type[evidence.type].append(evidence)

            # Resolve conflicts within each type
            for evidence_type, type_evidence in by_type.items():
                if len(type_evidence) == 1:
                    resolved_evidence.extend(type_evidence)
                else:
                    # Multiple evidence of same type - pick best
                    best_evidence = self._select_best_evidence(type_evidence)
                    resolved_evidence.append(best_evidence)

        return resolved_evidence

    def _select_best_evidence(self, evidence_list: List[Evidence]) -> Evidence:
        """Select the best evidence from a list of conflicting evidence"""
        # Score each evidence item
        scored_evidence = []

        for evidence in evidence_list:
            score = evidence.confidence or 0.0

            # Boost score for classified evidence
            if evidence.t_group:
                score += 0.1

            # Boost score for representative domains
            if evidence.extra_attributes.get("is_manual_rep"):
                score += 0.2

            scored_evidence.append((score, evidence))

        # Return highest scoring evidence
        scored_evidence.sort(key=lambda x: x[0], reverse=True)
        return scored_evidence[0][1]

    def calculate_quality_metrics(self, evidence_list: List[Evidence],
                                sequence_length: int, processing_time: float) -> EvidenceQualityMetrics:
        """Calculate comprehensive quality metrics"""
        metrics = EvidenceQualityMetrics()
        metrics.total_evidence_count = len(evidence_list)
        metrics.processing_time_seconds = processing_time

        if not evidence_list:
            return metrics

        # Count by type
        type_counts = Counter(e.type for e in evidence_list)
        metrics.blast_hit_count = type_counts.get("blast", 0) + type_counts.get("domain_blast", 0) + type_counts.get("chain_blast", 0)
        metrics.hhsearch_hit_count = type_counts.get("hhsearch", 0)
        metrics.self_comparison_count = type_counts.get("self_comparison", 0)

        # Valid evidence (has confidence score)
        valid_evidence = [e for e in evidence_list if e.confidence is not None]
        metrics.valid_evidence_count = len(valid_evidence)

        if valid_evidence:
            # Confidence statistics
            confidences = [e.confidence for e in valid_evidence]
            metrics.average_confidence = statistics.mean(confidences)
            metrics.confidence_std = statistics.stdev(confidences) if len(confidences) > 1 else 0.0
            metrics.high_confidence_count = sum(1 for c in confidences if c > 0.8)
            metrics.low_confidence_count = sum(1 for c in confidences if c < 0.3)

        # Coverage analysis
        if sequence_length > 0:
            metrics.sequence_coverage = self._calculate_sequence_coverage(evidence_list, sequence_length)
            gaps = self._find_coverage_gaps(evidence_list, sequence_length)
            metrics.gap_count = len(gaps)
            metrics.largest_gap_size = max(gaps) if gaps else 0

            overlaps = self._find_overlaps(evidence_list)
            metrics.overlap_count = len(overlaps)
            metrics.largest_overlap_size = max(overlaps) if overlaps else 0

        # Classification metrics
        classified_evidence = [e for e in evidence_list if e.t_group]
        metrics.classified_evidence_count = len(classified_evidence)

        if classified_evidence:
            # Check classification consistency
            t_groups = [e.t_group for e in classified_evidence if e.t_group]
            unique_t_groups = set(t_groups)
            metrics.classification_consistency = 1.0 - (len(unique_t_groups) - 1) / max(len(t_groups), 1)

        return metrics

    def _calculate_sequence_coverage(self, evidence_list: List[Evidence], sequence_length: int) -> float:
        """Calculate what percentage of sequence is covered by evidence"""
        covered_positions = set()

        for evidence in evidence_list:
            ranges = self._parse_range_comprehensive(evidence.query_range)
            for start, end in ranges:
                covered_positions.update(range(start, end + 1))

        return len(covered_positions) / sequence_length * 100.0 if sequence_length > 0 else 0.0

    def _find_coverage_gaps(self, evidence_list: List[Evidence], sequence_length: int) -> List[int]:
        """Find gaps in evidence coverage"""
        covered_positions = set()

        for evidence in evidence_list:
            ranges = self._parse_range_comprehensive(evidence.query_range)
            for start, end in ranges:
                covered_positions.update(range(start, end + 1))

        # Find gaps
        gaps = []
        current_gap = 0

        for pos in range(1, sequence_length + 1):
            if pos not in covered_positions:
                current_gap += 1
            else:
                if current_gap > 0:
                    gaps.append(current_gap)
                    current_gap = 0

        if current_gap > 0:
            gaps.append(current_gap)

        return gaps

    def _find_overlaps(self, evidence_list: List[Evidence]) -> List[int]:
        """Find overlaps between evidence items"""
        overlaps = []

        # Get all ranges
        evidence_ranges = []
        for evidence in evidence_list:
            ranges = self._parse_range_comprehensive(evidence.query_range)
            evidence_ranges.extend(ranges)

        # Check for overlaps between all pairs
        for i, (start1, end1) in enumerate(evidence_ranges):
            for j, (start2, end2) in enumerate(evidence_ranges[i+1:], i+1):
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)

                if overlap_start <= overlap_end:
                    overlap_size = overlap_end - overlap_start + 1
                    overlaps.append(overlap_size)

        return overlaps

    def generate_domain_suggestions(self, evidence_list: List[Evidence],
                                  sequence_length: int) -> List[Dict[str, Any]]:
        """Generate domain boundary suggestions from evidence"""
        suggestions = []

        # Group evidence by position
        position_map = defaultdict(list)

        for evidence in evidence_list:
            ranges = self._parse_range_comprehensive(evidence.query_range)
            for start, end in ranges:
                for pos in range(start, end + 1):
                    position_map[pos].append(evidence)

        # Find consensus regions
        consensus_regions = self._find_consensus_regions(position_map, sequence_length)

        # Convert regions to domain suggestions
        for i, (start, end) in enumerate(consensus_regions):
            # Get evidence supporting this region
            supporting_evidence = []
            for pos in range(start, end + 1):
                supporting_evidence.extend(position_map[pos])

            # Remove duplicates
            unique_evidence = list({e.source_id: e for e in supporting_evidence}.values())

            # Calculate confidence
            confidences = [e.confidence for e in unique_evidence if e.confidence is not None]
            avg_confidence = statistics.mean(confidences) if confidences else 0.0

            # Get most common classification
            classifications = [e.t_group for e in unique_evidence if e.t_group]
            most_common_classification = Counter(classifications).most_common(1)
            t_group = most_common_classification[0][0] if most_common_classification else None

            suggestion = {
                "id": f"domain_{i+1}",
                "start": start,
                "end": end,
                "range": f"{start}-{end}",
                "confidence": avg_confidence,
                "evidence_count": len(unique_evidence),
                "t_group": t_group,
                "source": "evidence_consensus"
            }

            suggestions.append(suggestion)

        return suggestions

    def _find_consensus_regions(self, position_map: Dict[int, List[Evidence]],
                              sequence_length: int) -> List[Tuple[int, int]]:
        """Find regions with strong evidence consensus"""
        regions = []
        current_start = None
        min_evidence_count = 2  # Minimum evidence items to consider a consensus

        for pos in range(1, sequence_length + 1):
            evidence_count = len(position_map.get(pos, []))

            if evidence_count >= min_evidence_count:
                if current_start is None:
                    current_start = pos
            else:
                if current_start is not None:
                    # End current region
                    regions.append((current_start, pos - 1))
                    current_start = None

        # Close final region if needed
        if current_start is not None:
            regions.append((current_start, sequence_length))

        # Filter out very small regions
        min_domain_size = self.options.min_domain_size
        regions = [(start, end) for start, end in regions if end - start + 1 >= min_domain_size]

        return regions

    def analyze_classifications(self, evidence_list: List[Evidence]) -> Dict[str, Any]:
        """Analyze classification patterns in evidence"""
        analysis = {
            "total_evidence": len(evidence_list),
            "classified_evidence": 0,
            "t_group_distribution": Counter(),
            "h_group_distribution": Counter(),
            "x_group_distribution": Counter(),
            "a_group_distribution": Counter(),
            "representative_counts": {
                "manual_rep": 0,
                "f70": 0,
                "f40": 0,
                "f99": 0
            },
            "classification_conflicts": []
        }

        classified_evidence = []

        for evidence in evidence_list:
            if evidence.t_group:
                classified_evidence.append(evidence)
                analysis["t_group_distribution"][evidence.t_group] += 1

                if evidence.h_group:
                    analysis["h_group_distribution"][evidence.h_group] += 1

                if evidence.x_group:
                    analysis["x_group_distribution"][evidence.x_group] += 1

                if evidence.a_group:
                    analysis["a_group_distribution"][evidence.a_group] += 1

                # Count representative flags
                if evidence.extra_attributes.get("is_manual_rep"):
                    analysis["representative_counts"]["manual_rep"] += 1
                if evidence.extra_attributes.get("is_f70"):
                    analysis["representative_counts"]["f70"] += 1
                if evidence.extra_attributes.get("is_f40"):
                    analysis["representative_counts"]["f40"] += 1
                if evidence.extra_attributes.get("is_f99"):
                    analysis["representative_counts"]["f99"] += 1

        analysis["classified_evidence"] = len(classified_evidence)

        # Detect classification conflicts
        t_groups = list(analysis["t_group_distribution"].keys())
        if len(t_groups) > 1:
            analysis["classification_conflicts"] = t_groups

        return analysis

    def _summarize_validation_results(self, validation_results: List[ValidationResult]) -> Dict[str, Any]:
        """Summarize validation results"""
        total = len(validation_results)
        valid = sum(1 for r in validation_results if r.is_valid)

        all_errors = []
        all_warnings = []

        for result in validation_results:
            all_errors.extend(result.errors)
            all_warnings.extend(result.warnings)

        return {
            "total_evidence": total,
            "valid_evidence": valid,
            "invalid_evidence": total - valid,
            "validation_rate": (valid / total * 100) if total > 0 else 0,
            "total_errors": len(all_errors),
            "total_warnings": len(all_warnings),
            "common_errors": Counter(all_errors).most_common(5),
            "common_warnings": Counter(all_warnings).most_common(5)
        }

    def _load_ecod_hierarchy(self) -> Dict[str, Any]:
        """Load ECOD hierarchy from database"""
        # Placeholder - implement based on your database schema
        return {}

    def _load_representative_domains(self) -> Set[str]:
        """Load set of representative domain IDs"""
        # Placeholder - implement based on your database schema
        return set()

    def get_processing_stats(self) -> Dict[str, Any]:
        """Get processing performance statistics"""
        return {
            "total_files_processed": len(self._processing_times),
            "average_processing_time": statistics.mean(self._processing_times) if self._processing_times else 0,
            "median_processing_time": statistics.median(self._processing_times) if self._processing_times else 0,
            "cache_stats": self.classification_cache.get_stats() if self.classification_cache else None,
            "parallel_processing": self.executor is not None,
            "max_workers": self.options.max_workers if self.executor else 1
        }

    def get_cache_statistics(self) -> Dict[str, Any]:
        """Get cache statistics"""
        if self.classification_cache:
            return self.classification_cache.get_stats()
        return {}

    def clear_cache(self) -> None:
        """Clear classification cache"""
        if self.classification_cache:
            self.classification_cache.cleanup_expired()

    def validate_domain(self, domain: DomainModel, context: str = "",
                   sequence_length: int = 0) -> ValidationResult:
        """
        Validate a domain model.

        Args:
            domain: Domain to validate
            context: Context for validation
            sequence_length: Sequence length for range validation

        Returns:
            ValidationResult
        """
        try:
            errors = []
            warnings = []

            # Basic validation
            if not domain.id:
                errors.append("Domain ID is empty")

            if domain.start <= 0:
                errors.append(f"Invalid start position: {domain.start}")

            if domain.end <= 0:
                errors.append(f"Invalid end position: {domain.end}")

            if domain.start > domain.end:
                errors.append(f"Start position ({domain.start}) > end position ({domain.end})")

            # Size validation
            domain_size = domain.end - domain.start + 1
            if domain_size < self.options.min_domain_size:
                if self.options.validation_level == ValidationLevel.STRICT:
                    errors.append(f"Domain too small: {domain_size} < {self.options.min_domain_size}")
                else:
                    warnings.append(f"Domain size ({domain_size}) below minimum ({self.options.min_domain_size})")

            if (self.options.max_domain_size and
                domain_size > self.options.max_domain_size):
                warnings.append(f"Domain size ({domain_size}) above maximum ({self.options.max_domain_size})")

            # Sequence length validation
            if sequence_length > 0:
                if domain.end > sequence_length:
                    errors.append(f"Domain end ({domain.end}) beyond sequence length ({sequence_length})")

            # Confidence validation
            if domain.confidence is not None and not (0.0 <= domain.confidence <= 1.0):
                warnings.append(f"Invalid confidence value: {domain.confidence}")

            is_valid = len(errors) == 0

            return ValidationResult(
                is_valid=is_valid,
                errors=errors,
                warnings=warnings,
                context=context,
                validation_level=self.options.validation_level
            )

        except Exception as e:
            return ValidationResult(
                is_valid=False,
                errors=[f"Domain validation error: {str(e)}"],
                warnings=[],
                context=context,
                validation_level=self.options.validation_level
            )

    def group_evidence_by_position(self, evidence_list: List[Evidence],
                              window_size: int = 50) -> Dict[int, List[Evidence]]:
        """
        Group evidence by position windows.

        Args:
            evidence_list: List of evidence items to group
            window_size: Size of position windows for grouping

        Returns:
            dict: Position groups with evidence items
        """
        try:
            position_groups = defaultdict(list)

            for evidence in evidence_list:
                # Parse ranges from evidence
                ranges = self._parse_range_comprehensive(evidence.query_range)

                if not ranges:
                    # If no range, use a default group
                    position_groups[0].append(evidence)
                    continue

                # Group by window based on range midpoint
                for start, end in ranges:
                    midpoint = (start + end) // 2
                    window = midpoint // window_size
                    position_groups[window].append(evidence)

            # Convert defaultdict to regular dict
            return dict(position_groups)

        except Exception as e:
            self.logger.error(f"Error grouping evidence by position: {str(e)}")
            return {}

    def cleanup_resources(self) -> None:
        """Clean up resources"""
        if self.executor:
            self.executor.shutdown(wait=True)

        if self.classification_cache:
            expired_count = self.classification_cache.cleanup_expired()
            if expired_count > 0:
                self.logger.info(f"Cleaned up {expired_count} expired cache entries")

    def __enter__(self):
        """Context manager entry"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.cleanup_resources()
        return False
