#!/usr/bin/env python3
"""
Quality Assessment Without Confidence - FIXED VERSION

Removes confidence scoring to see quality distribution without the artificial cap.
Uses binary classification success and consistency metrics instead.

Usage:
    python assess_quality_no_confidence.py --scan-all
    python assess_quality_no_confidence.py --batch-name batch_031
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

# Import SequenceRange for proper domain handling
from ecod.core.sequence_range import SequenceRange

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class DomainQuality:
    """Quality metrics for a single domain (NO CONFIDENCE)"""
    domain_id: str
    range: str
    length: int
    source: str
    has_classification: bool
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    # REMOVED: confidence field entirely

@dataclass
class ProteinQuality:
    """Quality assessment for a complete protein (NO CONFIDENCE)"""
    protein_id: str
    pdb_id: str
    chain_id: str
    batch_name: str

    # Basic metrics
    total_domains: int
    sequence_length: int
    coverage_percentage: float
    residues_covered: int

    # Quality metrics (NO CONFIDENCE)
    classified_domains: int
    evidence_diversity: int  # Number of different evidence types

    # Quality scores (0-1) - REWEIGHTED WITHOUT CONFIDENCE
    coverage_score: float = 0.0
    classification_score: float = 0.0
    parsimony_score: float = 0.0  # NEW: replaces confidence
    overall_score: float = 0.0

    # Quality tier
    tier: str = "unknown"

    # Issues
    issues: List[str] = None

    # Domain details
    domains: List[DomainQuality] = None

    def __post_init__(self):
        if self.issues is None:
            self.issues = []
        if self.domains is None:
            self.domains = []

class QualityAssessmentNoConfidence:
    """Quality assessment without confidence scoring"""

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
        """Parse and assess a single mini result - NO CONFIDENCE VERSION"""

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
                    # NO CONFIDENCE!
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

            # Calculate quality scores WITHOUT CONFIDENCE
            coverage_score = min(1.0, coverage_percentage / 80.0)  # 80% coverage = perfect
            classification_score = classified_count / max(1, total_domains) if total_domains > 0 else 0.0
            
            # NEW: Parsimony score (replaces confidence)
            parsimony_score = self._calculate_parsimony_score(domains, sequence_length)

            # REWEIGHTED: Overall score without confidence
            overall_score = (
                classification_score * 0.5 +    # Classification most important (was 0.4)
                coverage_score * 0.3 +          # Coverage important (was 0.4) 
                parsimony_score * 0.2           # Parsimony replaces confidence (was 0.3)
            )

            # Determine tier with ADJUSTED THRESHOLDS
            tier = self._determine_tier_no_confidence(overall_score, coverage_percentage, 
                                                    classified_count, total_domains)

            # Identify issues
            issues = self._identify_issues(domains, coverage_percentage, sequence_length)

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
                issues=issues,
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
                issues=["parse_error"]
            )

    def _calculate_parsimony_score(self, domains: List[DomainQuality], sequence_length: int) -> float:
        """Calculate parsimony score - penalizes over-segmentation and tiny domains"""
        
        if not domains:
            return 0.0
        
        score = 1.0  # Start perfect
        
        # Penalty for tiny domains (likely spurious)
        tiny_domains = sum(1 for d in domains if d.length < 30)
        if tiny_domains > 0:
            score -= tiny_domains * 0.2  # Heavy penalty for tiny domains
        
        # Penalty for excessive segmentation
        domain_density = len(domains) / max(1, sequence_length / 100)  # domains per 100 residues
        if domain_density > 0.5:  # More than 1 domain per 200 residues is suspicious
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

        # ADJUSTED THRESHOLDS - more achievable without confidence cap
        if overall_score >= 0.75 and coverage >= 50 and classified >= total * 0.8:
            return "excellent"
        elif overall_score >= 0.60 and coverage >= 40 and classified >= total * 0.6:
            return "good"
        elif overall_score >= 0.40 and coverage >= 20:
            return "acceptable"
        else:
            return "poor"

    def _identify_issues(self, domains: List[DomainQuality], coverage: float,
                        sequence_length: int) -> List[str]:
        """Identify potential issues with the result"""
        issues = []

        if len(domains) == 0:
            issues.append("no_domains")
        elif len(domains) > 10:
            issues.append("too_many_domains")

        if coverage < 20:
            issues.append("low_coverage")
        elif coverage > 150:
            issues.append("overcoverage")

        if not any(d.has_classification for d in domains):
            issues.append("no_classification")

        if sequence_length < 50:
            issues.append("short_protein")

        # Check for overlapping domains
        if self._count_overlapping_domains(domains) > 0:
            issues.append("overlapping_domains")
        
        # Check for tiny domains
        tiny_domains = sum(1 for d in domains if d.length < 30)
        if tiny_domains > 0:
            issues.append("tiny_domains")

        return issues
    
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
    
    def generate_quality_report(self, results: List[ProteinQuality]) -> Dict:
        """Generate comprehensive quality report"""
        
        if not results:
            return {"error": "No results to analyze"}
        
        # Tier distribution
        tier_counts = {}
        for result in results:
            tier_counts[result.tier] = tier_counts.get(result.tier, 0) + 1
        
        # Quality metrics
        total_results = len(results)
        successful_results = [r for r in results if r.tier != "failed"]
        
        if successful_results:
            avg_coverage = sum(r.coverage_percentage for r in successful_results) / len(successful_results)
            avg_domains = sum(r.total_domains for r in successful_results) / len(successful_results)
            avg_classified = sum(r.classified_domains for r in successful_results) / len(successful_results)
            avg_parsimony = sum(r.parsimony_score for r in successful_results) / len(successful_results)
        else:
            avg_coverage = avg_domains = avg_classified = avg_parsimony = 0
        
        # Issue analysis
        issue_counts = {}
        for result in results:
            for issue in result.issues:
                issue_counts[issue] = issue_counts.get(issue, 0) + 1
        
        # Production readiness
        excellent = tier_counts.get("excellent", 0)
        good = tier_counts.get("good", 0)
        production_ready = excellent + good
        
        return {
            "summary": {
                "total_results": total_results,
                "successful_parses": len(successful_results),
                "production_ready": production_ready,
                "production_rate": production_ready / max(1, total_results) * 100
            },
            "tier_distribution": tier_counts,
            "quality_metrics": {
                "average_coverage": avg_coverage,
                "average_domains_per_protein": avg_domains,
                "average_classified_domains": avg_classified,
                "average_parsimony": avg_parsimony  # REPLACES confidence
            },
            "issue_analysis": issue_counts,
            "scoring_method": "no_confidence_parsimony_based"
        }
    
    def print_quality_summary(self, results: List[ProteinQuality]):
        """Print quality assessment summary"""
        
        report = self.generate_quality_report(results)
        
        print("🔍 Mini PyECOD Quality Assessment (NO CONFIDENCE)")
        print("=" * 60)
        print(f"Assessment Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Scoring Method: Classification + Coverage + Parsimony (NO confidence)\n")
        
        # Summary
        summary = report["summary"]
        print("📊 Summary:")
        print(f"  Total results:        {summary['total_results']:>8,}")
        print(f"  Successful parses:    {summary['successful_parses']:>8,}")
        print(f"  Production ready:     {summary['production_ready']:>8,}")
        print(f"  Production rate:      {summary['production_rate']:>7.1f}%")
        print()
        
        # Tier distribution
        print("🏆 Quality Tier Distribution:")
        tier_counts = report["tier_distribution"]
        for tier in ["excellent", "good", "acceptable", "poor", "failed"]:
            count = tier_counts.get(tier, 0)
            pct = count / max(1, summary['total_results']) * 100
            print(f"  {tier.capitalize():<12} {count:>8,} ({pct:>5.1f}%)")
        print()
        
        # Quality metrics  
        metrics = report["quality_metrics"]
        print("📈 Quality Metrics:")
        print(f"  Avg coverage:         {metrics['average_coverage']:>7.1f}%")
        print(f"  Avg domains/protein:  {metrics['average_domains_per_protein']:>7.1f}")
        print(f"  Avg classified/protein: {metrics['average_classified_domains']:>5.1f}")
        print(f"  Avg parsimony score:  {metrics['average_parsimony']:>7.3f}")
        print()
        
        # Top issues
        print("⚠️  Top Issues:")
        issue_counts = report["issue_analysis"]
        sorted_issues = sorted(issue_counts.items(), key=lambda x: -x[1])[:5]
        for issue, count in sorted_issues:
            pct = count / max(1, summary['total_results']) * 100
            print(f"  {issue:<20} {count:>6,} ({pct:>5.1f}%)")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Assess Quality of Mini PyECOD Results (NO CONFIDENCE)'
    )
    
    parser.add_argument('--scan-all', action='store_true',
                       help='Assess all available results')
    parser.add_argument('--batch-name', type=str,
                       help='Assess specific batch')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize assessment
    assessor = QualityAssessmentNoConfidence(args.config)
    
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
    
    # Print summary
    assessor.print_quality_summary(results)


if __name__ == "__main__":
    main()
