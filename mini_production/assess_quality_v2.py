#!/usr/bin/env python3
"""
Enhanced Quality Assessment for Mini PyECOD - Tier-Specific Analysis

Provides detailed breakdown of excellent and good tier results to validate
production readiness and catch logical inconsistencies.

Key enhancements:
- Tier-specific quality metrics
- Issue analysis by tier (catch bugs like no_classification in excellent)
- Coverage analysis for top tiers
- Production readiness validation

Usage:
    python enhanced_quality_assessment.py --scan-all
    python enhanced_quality_assessment.py --batch-name batch_031
"""

import os
import sys
import argparse
import yaml
import xml.etree.ElementTree as ET
import psycopg2
import psycopg2.extras
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
from datetime import datetime
import json
import logging
from collections import defaultdict

# Import SequenceRange for proper domain handling
from ecod.core.sequence_range import SequenceRange

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class DomainQuality:
    """Quality metrics for a single domain"""
    domain_id: str
    range: str
    length: int
    source: str
    has_classification: bool
    t_group: Optional[str] = None
    h_group: Optional[str] = None

@dataclass
class ProteinQuality:
    """Quality assessment for a complete protein"""
    protein_id: str
    pdb_id: str
    chain_id: str
    batch_name: str

    # Basic metrics
    total_domains: int
    sequence_length: int
    coverage_percentage: float
    residues_covered: int

    # Quality metrics
    classified_domains: int
    evidence_diversity: int

    # Quality scores (0-1)
    coverage_score: float = 0.0
    classification_score: float = 0.0
    parsimony_score: float = 0.0
    overall_score: float = 0.0

    # Quality tier
    tier: str = "unknown"

    # Issues (separated into quality vs system errors)
    quality_issues: List[str] = None
    system_errors: List[str] = None

    # Domain details
    domains: List[DomainQuality] = None

    def __post_init__(self):
        if self.quality_issues is None:
            self.quality_issues = []
        if self.system_errors is None:
            self.system_errors = []
        if self.domains is None:
            self.domains = []

    @property
    def issues(self):
        """Backward compatibility - all issues combined"""
        return self.quality_issues + self.system_errors

    @property
    def has_system_errors(self):
        """Check if this result has system errors"""
        return len(self.system_errors) > 0

class EnhancedQualityAssessment:
    """Enhanced quality assessment with tier-specific analysis"""

    def __init__(self, config_path: str = "config.local.yml"):
        self.config = self._load_config(config_path)
        self.db_conn = self._init_db_connection() if self.config.get('database') else None

    def _load_config(self, config_path: str) -> Dict:
        """Load configuration"""
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            return {"paths": {"batch_base_dir": "/data/ecod/pdb_updates/batches"}}

    def _init_db_connection(self):
        """Initialize database connection for sequence length lookup"""
        try:
            return psycopg2.connect(**self.config['database'])
        except Exception as e:
            logger.warning(f"Database connection failed: {e}")
            return None

    def get_sequence_length(self, pdb_id: str, chain_id: str) -> int:
        """Get sequence length from database or estimate from domains"""
        if not self.db_conn:
            return 0

        try:
            with self.db_conn.cursor() as cursor:
                cursor.execute("""
                    SELECT length FROM ecod_schema.protein
                    WHERE pdb_id = %s AND chain_id = %s
                    LIMIT 1
                """, (pdb_id, chain_id))

                result = cursor.fetchone()
                return result[0] if result else 0
        except Exception:
            return 0

    def parse_mini_result(self, xml_file: Path) -> ProteinQuality:
        """Parse and assess a single mini result"""

        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()

            # Extract protein info
            protein_id = xml_file.stem.replace('.mini.domains', '')
            pdb_id = protein_id.split('_')[0]
            chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
            batch_name = xml_file.parent.parent.name

            # Parse domains
            domains = []
            domain_elements = root.findall(".//domain")

            evidence_types = set()
            classified_count = 0
            total_length = 0

            for i, domain_elem in enumerate(domain_elements):
                # Extract domain data
                range_str = domain_elem.get('range', '')

                # Use SequenceRange to parse discontinuous domains
                try:
                    seq_range = SequenceRange.parse(range_str)
                    start_pos, end_pos = seq_range.span
                    length = seq_range.total_length
                except (ValueError, Exception) as e:
                    logger.warning(f"Could not parse range '{range_str}' for {protein_id}: {e}")
                    start_pos, end_pos, length = 0, 0, 0

                source = domain_elem.get('source', 'unknown')
                t_group = domain_elem.get('t_group')
                has_classification = t_group is not None and t_group.strip() != ''

                domain_quality = DomainQuality(
                    domain_id=f"{protein_id}_d{i+1}",
                    range=range_str,
                    length=length,
                    source=source,
                    has_classification=has_classification,
                    t_group=t_group,
                    h_group=domain_elem.get('h_group')
                )

                domains.append(domain_quality)
                evidence_types.add(source)
                total_length += length

                if has_classification:
                    classified_count += 1

            # Get sequence length
            sequence_length = self.get_sequence_length(pdb_id, chain_id)
            if sequence_length == 0:
                # Estimate from domain coverage
                try:
                    max_end = max((SequenceRange.parse(d.range).span[1] for d in domains), default=0)
                    sequence_length = int(max_end * 1.1)  # Add 10% buffer
                except Exception:
                    sequence_length = total_length + 50  # Fallback estimation

            # Calculate metrics
            total_domains = len(domains)
            coverage_percentage = (total_length / max(1, sequence_length)) * 100
            evidence_diversity = len(evidence_types)

            # Calculate quality scores
            coverage_score = min(1.0, coverage_percentage / 80.0)  # 80% coverage = perfect
            classification_score = classified_count / max(1, total_domains) if total_domains > 0 else 0.0
            
            # Parsimony score (replaces confidence)
            parsimony_score = self._calculate_parsimony_score(domains, sequence_length)

            # Overall score
            overall_score = (
                classification_score * 0.5 +
                coverage_score * 0.3 +
                parsimony_score * 0.2
            )

            # Determine tier
            tier = self._determine_tier_no_confidence(overall_score, coverage_percentage, 
                                                    classified_count, total_domains)

            # Identify issues (separated into quality vs system errors)
            quality_issues, system_errors = self._identify_issues(domains, coverage_percentage, sequence_length)

            return ProteinQuality(
                protein_id=protein_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                batch_name=batch_name,
                total_domains=total_domains,
                sequence_length=sequence_length,
                coverage_percentage=coverage_percentage,
                residues_covered=total_length,
                classified_domains=classified_count,
                evidence_diversity=evidence_diversity,
                coverage_score=coverage_score,
                classification_score=classification_score,
                parsimony_score=parsimony_score,
                overall_score=overall_score,
                tier=tier,
                quality_issues=quality_issues,
                system_errors=system_errors,
                domains=domains
            )

        except Exception as e:
            logger.error(f"Error parsing {xml_file}: {e}")
            # Return minimal quality object for failed parse
            protein_id = xml_file.stem.replace('.mini.domains', '')
            return ProteinQuality(
                protein_id=protein_id,
                pdb_id=protein_id.split('_')[0],
                chain_id=protein_id.split('_')[1] if '_' in protein_id else 'A',
                batch_name=xml_file.parent.parent.name,
                total_domains=0,
                sequence_length=0,
                coverage_percentage=0,
                residues_covered=0,
                classified_domains=0,
                evidence_diversity=0,
                tier="failed",
                quality_issues=[],
                system_errors=["parse_error"]
            )

    def _calculate_parsimony_score(self, domains: List[DomainQuality], sequence_length: int) -> float:
        """Calculate parsimony score - penalizes over-segmentation and tiny domains"""
        
        if not domains:
            return 0.0
        
        score = 1.0  # Start perfect
        
        # Penalty for tiny domains (likely spurious)
        tiny_domains = sum(1 for d in domains if d.length < 30)
        if tiny_domains > 0:
            score -= tiny_domains * 0.2
        
        # Penalty for excessive segmentation
        domain_density = len(domains) / max(1, sequence_length / 100)
        if domain_density > 0.5:
            score -= (domain_density - 0.5) * 0.3
        
        # Penalty for overlapping domains
        overlaps = self._count_overlapping_domains(domains)
        if overlaps > 0:
            score -= overlaps * 0.15
        
        return max(0.0, min(1.0, score))

    def _count_overlapping_domains(self, domains: List[DomainQuality]) -> int:
        """Count overlapping domains"""
        overlaps = 0
        
        try:
            # Parse all domain ranges
            domain_ranges = []
            for domain in domains:
                try:
                    seq_range = SequenceRange.parse(domain.range)
                    domain_ranges.append(seq_range)
                except Exception:
                    continue

            # Check for overlaps
            for i in range(len(domain_ranges)):
                for j in range(i + 1, len(domain_ranges)):
                    if domain_ranges[i].overlaps(domain_ranges[j]):
                        overlaps += 1
                        
        except Exception as e:
            logger.warning(f"Could not check domain overlaps: {e}")
        
        return overlaps

    def _determine_tier_no_confidence(self, overall_score: float, coverage: float,
                                    classified: int, total: int) -> str:
        """Determine quality tier WITHOUT confidence requirement"""

        if overall_score >= 0.75 and coverage >= 50 and classified >= total * 0.8:
            return "excellent"
        elif overall_score >= 0.60 and coverage >= 40 and classified >= total * 0.6:
            return "good"
        elif overall_score >= 0.40 and coverage >= 20:
            return "acceptable"
        else:
            return "poor"

    def _identify_issues(self, domains: List[DomainQuality], coverage: float,
                        sequence_length: int) -> Tuple[List[str], List[str]]:
        """Identify quality issues vs system errors"""
        quality_issues = []
        system_errors = []

        if len(domains) == 0:
            quality_issues.append("no_domains")
        elif len(domains) > 10:
            quality_issues.append("too_many_domains")

        if coverage < 20:
            quality_issues.append("low_coverage")
        elif coverage > 150:
            quality_issues.append("overcoverage")

        # CRITICAL: no_classification should be impossible in this toolchain
        # If we see it, it's a system error (parse/data sync), not quality issue
        if not any(d.has_classification for d in domains) and len(domains) > 0:
            system_errors.append("no_classification_SYSTEM_ERROR")

        if sequence_length < 50:
            quality_issues.append("short_protein")

        # Check for overlapping domains
        if self._count_overlapping_domains(domains) > 0:
            quality_issues.append("overlapping_domains")
        
        # Check for tiny domains
        tiny_domains = sum(1 for d in domains if d.length < 30)
        if tiny_domains > 0:
            quality_issues.append("tiny_domains")

        return quality_issues, system_errors
    
    def assess_batch(self, batch_name: str) -> List[ProteinQuality]:
        """Assess all results in a specific batch"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        batch_dir = batch_base / batch_name
        mini_domains_dir = batch_dir / "mini_domains"
        
        if not mini_domains_dir.exists():
            logger.error(f"No mini results found in {batch_name}")
            return []
        
        xml_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
        logger.info(f"Assessing {len(xml_files)} results in {batch_name}")
        
        results = []
        for xml_file in xml_files:
            quality = self.parse_mini_result(xml_file)
            results.append(quality)
        
        return results
    
    def assess_all_batches(self) -> List[ProteinQuality]:
        """Assess all available mini results"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        all_results = []
        
        for batch_dir in batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
                
            mini_domains_dir = batch_dir / "mini_domains"
            if not mini_domains_dir.exists():
                continue
            
            batch_results = self.assess_batch(batch_dir.name)
            all_results.extend(batch_results)
            
            logger.info(f"Batch {batch_dir.name}: {len(batch_results)} results assessed")
        
        logger.info(f"Total: {len(all_results)} results assessed")
        return all_results
    
    def generate_tier_specific_report(self, results: List[ProteinQuality]) -> Dict:
        """Generate detailed tier-specific analysis"""
        
        # Group results by tier
        by_tier = defaultdict(list)
        for result in results:
            by_tier[result.tier].append(result)
        
        tier_analysis = {}
        
        for tier, tier_results in by_tier.items():
            if not tier_results:
                continue
                
            # Calculate tier-specific metrics
            avg_coverage = sum(r.coverage_percentage for r in tier_results) / len(tier_results)
            avg_domains = sum(r.total_domains for r in tier_results) / len(tier_results)
            avg_classified = sum(r.classified_domains for r in tier_results) / len(tier_results)
            avg_classification_rate = sum(r.classification_score for r in tier_results) / len(tier_results)
            avg_parsimony = sum(r.parsimony_score for r in tier_results) / len(tier_results)
            avg_overall_score = sum(r.overall_score for r in tier_results) / len(tier_results)
            
            # Coverage distribution
            coverage_ranges = {
                "0-25%": sum(1 for r in tier_results if r.coverage_percentage < 25),
                "25-50%": sum(1 for r in tier_results if 25 <= r.coverage_percentage < 50),
                "50-75%": sum(1 for r in tier_results if 50 <= r.coverage_percentage < 75),
                "75-100%": sum(1 for r in tier_results if 75 <= r.coverage_percentage < 100),
                "100%+": sum(1 for r in tier_results if r.coverage_percentage >= 100)
            }
            
            # Domain count distribution
            domain_ranges = {
                "1": sum(1 for r in tier_results if r.total_domains == 1),
                "2-3": sum(1 for r in tier_results if 2 <= r.total_domains <= 3),
                "4-5": sum(1 for r in tier_results if 4 <= r.total_domains <= 5),
                "6+": sum(1 for r in tier_results if r.total_domains >= 6)
            }
            
            tier_analysis[tier] = {
                "count": len(tier_results),
                "metrics": {
                    "avg_coverage": avg_coverage,
                    "avg_domains": avg_domains,
                    "avg_classified": avg_classified,
                    "avg_classification_rate": avg_classification_rate,
                    "avg_parsimony": avg_parsimony,
                    "avg_overall_score": avg_overall_score
                },
                "coverage_distribution": coverage_ranges,
                "domain_count_distribution": domain_ranges
            }
        
        return tier_analysis
    
    def diagnose_no_classification_errors(self, results: List[ProteinQuality]) -> Dict:
        """Diagnose the cause of no_classification system errors"""
        
        no_class_errors = [r for r in results if "no_classification_SYSTEM_ERROR" in r.system_errors]
        
        if not no_class_errors:
            return {"no_errors": True}
        
        diagnostics = []
        
        for result in no_class_errors:
            # Examine the raw XML to diagnose the issue
            batch_base = Path(self.config["paths"]["batch_base_dir"])
            xml_file = batch_base / result.batch_name / "mini_domains" / f"{result.protein_id}.mini.domains.xml"
            
            diagnosis = {
                "protein_id": result.protein_id,
                "domains_found": result.total_domains,
                "xml_exists": xml_file.exists(),
                "diagnosis": "unknown"
            }
            
            if xml_file.exists():
                try:
                    tree = ET.parse(xml_file)
                    root = tree.getroot()
                    domain_elements = root.findall(".//domain")
                    
                    # Check what classification fields are actually present
                    classification_fields = []
                    for domain_elem in domain_elements:
                        fields = {}
                        for attr in ['t_group', 'h_group', 'x_group', 'family']:
                            value = domain_elem.get(attr)
                            fields[attr] = value if value and value.strip() else None
                        classification_fields.append(fields)
                    
                    diagnosis["classification_fields"] = classification_fields
                    
                    # Determine likely cause
                    if not domain_elements:
                        diagnosis["diagnosis"] = "no_domain_elements_in_xml"
                    elif all(not any(f.values()) for f in classification_fields):
                        diagnosis["diagnosis"] = "domains_present_but_no_classification_data"
                    elif all(f.get('t_group') is None for f in classification_fields):
                        diagnosis["diagnosis"] = "missing_t_group_field_specifically"
                    else:
                        diagnosis["diagnosis"] = "partial_classification_data"
                        
                except Exception as e:
                    diagnosis["diagnosis"] = f"xml_parse_error: {e}"
            else:
                diagnosis["diagnosis"] = "xml_file_missing"
            
            diagnostics.append(diagnosis)
        
        # Summarize common patterns
        diagnosis_counts = {}
        for d in diagnostics:
            diag = d["diagnosis"]
            diagnosis_counts[diag] = diagnosis_counts.get(diag, 0) + 1
        
        return {
            "total_errors": len(no_class_errors),
            "error_rate": len(no_class_errors) / len(results) * 100,
            "diagnosis_summary": diagnosis_counts,
            "detailed_diagnostics": diagnostics
        }

    def analyze_issues_by_tier(self, results: List[ProteinQuality]) -> Dict:
        """Analyze which issues occur in which tiers"""

        quality_issue_by_tier = defaultdict(lambda: defaultdict(int))
        system_error_by_tier = defaultdict(lambda: defaultdict(int))
        logical_inconsistencies = []

        for result in results:
            tier = result.tier

            # Track quality issues
            for issue in result.quality_issues:
                quality_issue_by_tier[tier][issue] += 1

            # Track system errors
            for error in result.system_errors:
                system_error_by_tier[tier][error] += 1

                # Any system error in excellent/good tier is a logical inconsistency
                if tier in ["excellent", "good"]:
                    logical_inconsistencies.append({
                        "protein_id": result.protein_id,
                        "tier": tier,
                        "error_type": "system_error",
                        "error": error,
                        "classification_score": result.classification_score,
                        "classified_domains": result.classified_domains,
                        "total_domains": result.total_domains
                    })

        return {
            "quality_issue_by_tier": dict(quality_issue_by_tier),
            "system_error_by_tier": dict(system_error_by_tier),
            "logical_inconsistencies": logical_inconsistencies
        }

    def generate_tier_specific_report(self, results: List[ProteinQuality]) -> Dict:
        """Generate detailed tier-specific analysis"""

        # Group results by tier
        by_tier = defaultdict(list)
        for result in results:
            by_tier[result.tier].append(result)

        tier_analysis = {}

        for tier, tier_results in by_tier.items():
            if not tier_results:
                continue

            # Calculate tier-specific metrics
            avg_coverage = sum(r.coverage_percentage for r in tier_results) / len(tier_results)
            avg_domains = sum(r.total_domains for r in tier_results) / len(tier_results)
            avg_classified = sum(r.classified_domains for r in tier_results) / len(tier_results)
            avg_classification_rate = sum(r.classification_score for r in tier_results) / len(tier_results)
            avg_parsimony = sum(r.parsimony_score for r in tier_results) / len(tier_results)
            avg_overall_score = sum(r.overall_score for r in tier_results) / len(tier_results)

            # Coverage distribution
            coverage_ranges = {
                "0-25%": sum(1 for r in tier_results if r.coverage_percentage < 25),
                "25-50%": sum(1 for r in tier_results if 25 <= r.coverage_percentage < 50),
                "50-75%": sum(1 for r in tier_results if 50 <= r.coverage_percentage < 75),
                "75-100%": sum(1 for r in tier_results if 75 <= r.coverage_percentage < 100),
                "100%+": sum(1 for r in tier_results if r.coverage_percentage >= 100)
            }

            # Domain count distribution
            domain_ranges = {
                "1": sum(1 for r in tier_results if r.total_domains == 1),
                "2-3": sum(1 for r in tier_results if 2 <= r.total_domains <= 3),
                "4-5": sum(1 for r in tier_results if 4 <= r.total_domains <= 5),
                "6+": sum(1 for r in tier_results if r.total_domains >= 6)
            }

            tier_analysis[tier] = {
                "count": len(tier_results),
                "metrics": {
                    "avg_coverage": avg_coverage,
                    "avg_domains": avg_domains,
                    "avg_classified": avg_classified,
                    "avg_classification_rate": avg_classification_rate,
                    "avg_parsimony": avg_parsimony,
                    "avg_overall_score": avg_overall_score
                },
                "coverage_distribution": coverage_ranges,
                "domain_count_distribution": domain_ranges
            }

        return tier_analysis
        """Analyze which issues occur in which tiers"""
        
        issue_by_tier = defaultdict(lambda: defaultdict(int))
        logical_inconsistencies = []
        
        for result in results:
            tier = result.tier
            for issue in result.issues:
                issue_by_tier[tier][issue] += 1
                
                # Check for logical inconsistencies
                if tier in ["excellent", "good"] and issue == "no_classification":
                    logical_inconsistencies.append({
                        "protein_id": result.protein_id,
                        "tier": tier,
                        "issue": issue,
                        "classification_score": result.classification_score,
                        "classified_domains": result.classified_domains,
                        "total_domains": result.total_domains
                    })
        
        return {
            "issue_by_tier": dict(issue_by_tier),
            "logical_inconsistencies": logical_inconsistencies
        }
    
    def generate_quality_report(self, results: List[ProteinQuality]) -> Dict:
        """Generate comprehensive quality report with tier analysis"""
        
        if not results:
            return {"error": "No results to analyze"}
        
        # Basic report (existing)
        tier_counts = {}
        for result in results:
            tier_counts[result.tier] = tier_counts.get(result.tier, 0) + 1
        
        total_results = len(results)
        successful_results = [r for r in results if r.tier != "failed"]
        
        if successful_results:
            avg_coverage = sum(r.coverage_percentage for r in successful_results) / len(successful_results)
            avg_domains = sum(r.total_domains for r in successful_results) / len(successful_results)
            avg_classified = sum(r.classified_domains for r in successful_results) / len(successful_results)
            avg_parsimony = sum(r.parsimony_score for r in successful_results) / len(successful_results)
        else:
            avg_coverage = avg_domains = avg_classified = avg_parsimony = 0
        
        issue_counts = {}
        system_error_counts = {}
        for result in results:
            for issue in result.quality_issues:
                issue_counts[issue] = issue_counts.get(issue, 0) + 1
            for error in result.system_errors:
                system_error_counts[error] = system_error_counts.get(error, 0) + 1
        
        excellent = tier_counts.get("excellent", 0)
        good = tier_counts.get("good", 0)
        production_ready = excellent + good
        
        # Enhanced analysis
        tier_analysis = self.generate_tier_specific_report(results)
        issue_analysis = self.analyze_issues_by_tier(results)
        
        # CRITICAL: Diagnose no_classification system errors
        no_class_diagnosis = self.diagnose_no_classification_errors(results)
        
        return {
            "summary": {
                "total_results": total_results,
                "successful_parses": len(successful_results),
                "production_ready": production_ready,
                "production_rate": production_ready / max(1, total_results) * 100,
                "system_errors_detected": len([r for r in results if r.has_system_errors])
            },
            "tier_distribution": tier_counts,
            "quality_metrics": {
                "average_coverage": avg_coverage,
                "average_domains_per_protein": avg_domains,
                "average_classified_domains": avg_classified,
                "average_parsimony": avg_parsimony
            },
            "quality_issue_analysis": issue_counts,
            "system_error_analysis": system_error_counts,
            "tier_specific_analysis": tier_analysis,
            "issue_by_tier_analysis": issue_analysis,
            "no_classification_diagnosis": no_class_diagnosis,
            "scoring_method": "no_confidence_parsimony_based"
        }
    
    def print_quality_summary(self, results: List[ProteinQuality]):
        """Print comprehensive quality assessment summary"""
        
        report = self.generate_quality_report(results)
        
        print("🔍 Enhanced Mini PyECOD Quality Assessment")
        print("=" * 70)
        print(f"Assessment Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Scoring Method: Classification + Coverage + Parsimony (NO confidence)\n")
        
        # Summary
        summary = report["summary"]
        print("📊 Overall Summary:")
        print(f"  Total results:        {summary['total_results']:>8,}")
        print(f"  Successful parses:    {summary['successful_parses']:>8,}")
        print(f"  Production ready:     {summary['production_ready']:>8,}")
        print(f"  Production rate:      {summary['production_rate']:>7.1f}%")
        print(f"  System errors:        {summary['system_errors_detected']:>8,}")
        print()
        
        # CRITICAL: System error analysis first (these are bugs, not quality issues)
        system_errors = report["system_error_analysis"]
        if system_errors:
            print("🚨 SYSTEM ERRORS DETECTED (These are bugs, not quality issues):")
            for error, count in sorted(system_errors.items(), key=lambda x: -x[1]):
                pct = count / summary['total_results'] * 100
                print(f"  {error:<30} {count:>6,} ({pct:>5.1f}%)")
            
            # Special diagnosis for no_classification errors
            no_class_diag = report["no_classification_diagnosis"]
            if not no_class_diag.get("no_errors", False):
                print(f"\n🔍 NO_CLASSIFICATION ERROR DIAGNOSIS:")
                print(f"  Total affected proteins: {no_class_diag['total_errors']}")
                print(f"  Error rate: {no_class_diag['error_rate']:.1f}%")
                print(f"  Root causes identified:")
                for cause, count in no_class_diag["diagnosis_summary"].items():
                    print(f"    {cause:<35} {count:>4}")
                
                # Show specific examples for debugging
                if no_class_diag["detailed_diagnostics"]:
                    print(f"\n  Examples for debugging:")
                    for diag in no_class_diag["detailed_diagnostics"][:3]:  # Show first 3
                        print(f"    {diag['protein_id']}: {diag['diagnosis']}")
            print()
        else:
            print("✅ No system errors detected")
            print()
        
        # Tier distribution
        print("🏆 Quality Tier Distribution:")
        tier_counts = report["tier_distribution"]
        for tier in ["excellent", "good", "acceptable", "poor", "failed"]:
            count = tier_counts.get(tier, 0)
            pct = count / max(1, summary['total_results']) * 100
            print(f"  {tier.capitalize():<12} {count:>8,} ({pct:>5.1f}%)")
        print()
        
        # ENHANCED: Tier-specific analysis
        print("🎯 Tier-Specific Analysis:")
        tier_analysis = report["tier_specific_analysis"]
        
        for tier in ["excellent", "good"]:
            if tier not in tier_analysis:
                continue
                
            data = tier_analysis[tier]
            metrics = data["metrics"]
            
            print(f"\n  {tier.upper()} Tier ({data['count']} proteins):")
            print(f"    Average coverage:       {metrics['avg_coverage']:>7.1f}%")
            print(f"    Average domains:        {metrics['avg_domains']:>7.1f}")
            print(f"    Average classified:     {metrics['avg_classified']:>7.1f}")
            print(f"    Classification rate:    {metrics['avg_classification_rate']:>7.1%}")
            print(f"    Parsimony score:        {metrics['avg_parsimony']:>7.3f}")
            print(f"    Overall score:          {metrics['avg_overall_score']:>7.3f}")
            
            # Coverage distribution for this tier
            print(f"    Coverage distribution:")
            for range_name, count in data["coverage_distribution"].items():
                if count > 0:
                    pct = count / data['count'] * 100
                    print(f"      {range_name:<8} {count:>4} ({pct:>5.1f}%)")
        
        print()
        
        # ENHANCED: Issue analysis by tier (QUALITY ISSUES ONLY)
        print("⚠️  Quality Issue Analysis by Tier:")
        issue_analysis = report["issue_by_tier_analysis"]
        quality_issue_by_tier = issue_analysis["quality_issue_by_tier"]
        
        for tier in ["excellent", "good", "acceptable", "poor"]:
            if tier not in quality_issue_by_tier:
                continue
                
            tier_issues = quality_issue_by_tier[tier]
            if not tier_issues:
                print(f"  {tier.capitalize():<12} No quality issues detected ✅")
                continue
                
            print(f"  {tier.capitalize():<12} Quality Issues:")
            for issue, count in sorted(tier_issues.items(), key=lambda x: -x[1]):
                tier_count = tier_analysis.get(tier, {}).get('count', 1)
                pct = count / tier_count * 100
                print(f"    {issue:<20} {count:>4} ({pct:>5.1f}%)")
        
        # Show system errors by tier separately
        system_error_by_tier = issue_analysis["system_error_by_tier"]
        if any(system_error_by_tier.values()):
            print(f"\n🚨 System Errors by Tier:")
            for tier in ["excellent", "good", "acceptable", "poor"]:
                if tier not in system_error_by_tier:
                    continue
                    
                tier_errors = system_error_by_tier[tier]
                if not tier_errors:
                    continue
                    
                print(f"  {tier.capitalize():<12} System Errors:")
                for error, count in sorted(tier_errors.items(), key=lambda x: -x[1]):
                    tier_count = tier_analysis.get(tier, {}).get('count', 1)
                    pct = count / tier_count * 100
                    print(f"    {error:<30} {count:>4} ({pct:>5.1f}%)")
        
        print()
        
        # CRITICAL: Check for logical inconsistencies (system errors in high tiers)
        inconsistencies = issue_analysis["logical_inconsistencies"]
        if inconsistencies:
            print("🚨 CRITICAL BUGS DETECTED:")
            print("    (System errors should never appear in excellent/good tiers)")
            for inc in inconsistencies:
                print(f"    {inc['protein_id']}: {inc['tier']} tier but has system error '{inc['error']}'")
                print(f"      - Classification score: {inc['classification_score']:.3f}")
                print(f"      - Classified domains: {inc['classified_domains']}/{inc['total_domains']}")
            print(f"    → These {len(inconsistencies)} cases indicate bugs in partition-classify integration")
            print()
        else:
            print("✅ No logical inconsistencies detected (system errors properly isolated)")
            print()
        
        # Overall quality metrics  
        metrics = report["quality_metrics"]
        print("📈 Overall Quality Metrics:")
        print(f"  Avg coverage:         {metrics['average_coverage']:>7.1f}%")
        print(f"  Avg domains/protein:  {metrics['average_domains_per_protein']:>7.1f}")
        print(f"  Avg classified/protein: {metrics['average_classified_domains']:>5.1f}")
        print(f"  Avg parsimony score:  {metrics['average_parsimony']:>7.3f}")
        print()
        
        # Quality issues summary (separate from system errors)
        quality_issues = report["quality_issue_analysis"]
        if quality_issues:
            print("📋 Quality Issues Summary (algorithmic, not bugs):")
            sorted_issues = sorted(quality_issues.items(), key=lambda x: -x[1])[:5]
            for issue, count in sorted_issues:
                pct = count / max(1, summary['total_results']) * 100
                print(f"  {issue:<20} {count:>6,} ({pct:>5.1f}%)")
            print()
        
        # Production readiness assessment
        excellent_count = tier_counts.get("excellent", 0)
        good_count = tier_counts.get("good", 0)
        system_error_count = summary['system_errors_detected']
        
        print("🚀 Production Readiness Assessment:")
        
        # System errors are production blockers
        if system_error_count > 0:
            system_error_rate = system_error_count / summary['total_results'] * 100
            print("  ❌ NOT READY FOR PRODUCTION - SYSTEM ERRORS DETECTED")
            print(f"     {system_error_count} proteins have system errors ({system_error_rate:.1f}%)")
            print("     → Fix partition-classify integration bugs before export")
            
            # Show what needs to be fixed
            system_errors = report["system_error_analysis"]
            print("     Priority fixes needed:")
            for error, count in sorted(system_errors.items(), key=lambda x: -x[1]):
                print(f"       - {error}: {count} cases")
                
        elif excellent_count + good_count > summary['total_results'] * 0.7:
            print("  ✅ READY FOR PRODUCTION")
            print(f"     {excellent_count + good_count} high-quality results ({(excellent_count + good_count)/summary['total_results']*100:.1f}%)")
            print("     No system errors detected")
        elif excellent_count + good_count > summary['total_results'] * 0.5:
            print("  ⚠️  CONDITIONALLY READY")
            print(f"     {excellent_count + good_count} high-quality results ({(excellent_count + good_count)/summary['total_results']*100:.1f}%)")
            print("     No system errors, but consider additional quality validation")
        else:
            print("  ❌ NOT READY FOR PRODUCTION - LOW QUALITY")
            print(f"     Only {excellent_count + good_count} high-quality results ({(excellent_count + good_count)/summary['total_results']*100:.1f}%)")
            print("     Algorithm needs quality improvement")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Enhanced Quality Assessment for Mini PyECOD Results'
    )
    
    parser.add_argument('--scan-all', action='store_true',
                       help='Assess all available results')
    parser.add_argument('--batch-name', type=str,
                       help='Assess specific batch')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize assessment
    assessor = EnhancedQualityAssessment(args.config)
    
    # Assess results
    if args.batch_name:
        results = assessor.assess_batch(args.batch_name)
    elif args.scan_all:
        results = assessor.assess_all_batches()
    else:
        parser.print_help()
        return
    
    if not results:
        print("No results found to assess.")
        return
    
    # Print enhanced summary
    assessor.print_quality_summary(results)


if __name__ == "__main__":
    main()
