#!/usr/bin/env python3
"""
Quality Assessment for Mini PyECOD Results - COMPLETE FIXED VERSION

Multi-tier assessment strategy:
1. Automated quality filtering (all results)
2. Coverage-based ranking
3. Evidence strength analysis
4. Production readiness scoring

FIXES:
- Properly handles discontinuous domains using SequenceRange
- Fixes all range parsing issues
- Correct overlap detection for discontinuous domains

Usage:
    python assess_quality.py --scan-all                    # Assess all results
    python assess_quality.py --batch-name batch_031       # Assess specific batch
    python assess_quality.py --export-production          # Export high-quality results
    python assess_quality.py --dashboard-prep             # Prepare data for dashboard
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

# FIXED: Import SequenceRange for proper discontinuous domain handling
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
    confidence: float
    has_classification: bool
    t_group: Optional[str] = None
    h_group: Optional[str] = None
    evidence_score: float = 0.0  # Computed score based on source and confidence

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
    average_confidence: float
    evidence_diversity: int  # Number of different evidence types

    # Quality scores (0-1)
    coverage_score: float = 0.0
    classification_score: float = 0.0
    confidence_score: float = 0.0
    overall_score: float = 0.0

    # Quality tier
    tier: str = "unknown"  # "excellent", "good", "acceptable", "poor"

    # Issues
    issues: List[str] = None

    # Domain details
    domains: List[DomainQuality] = None

    def __post_init__(self):
        if self.issues is None:
            self.issues = []
        if self.domains is None:
            self.domains = []

class QualityAssessment:
    """Assess quality of mini PyECOD results"""

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
        """Parse and assess a single mini result - FIXED for discontinuous domains"""

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
            total_confidence = 0.0
            classified_count = 0
            total_length = 0

            for i, domain_elem in enumerate(domain_elements):
                # Extract domain data
                range_str = domain_elem.get('range', '')

                # FIXED: Use SequenceRange to parse discontinuous domains
                try:
                    seq_range = SequenceRange.parse(range_str)
                    start_pos, end_pos = seq_range.span  # Overall start and end
                    length = seq_range.total_length      # Correct total length
                except (ValueError, Exception) as e:
                    logger.warning(f"Could not parse range '{range_str}' for {protein_id}: {e}")
                    start_pos, end_pos, length = 0, 0, 0

                source = domain_elem.get('source', 'unknown')
                confidence = float(domain_elem.get('confidence', '0.0'))
                t_group = domain_elem.get('t_group')
                has_classification = t_group is not None and t_group.strip() != ''

                # Calculate evidence score
                evidence_score = self._calculate_evidence_score(source, confidence)

                domain_quality = DomainQuality(
                    domain_id=f"{protein_id}_d{i+1}",
                    range=range_str,
                    length=length,
                    source=source,
                    confidence=confidence,
                    has_classification=has_classification,
                    t_group=t_group,
                    h_group=domain_elem.get('h_group'),
                    evidence_score=evidence_score
                )

                domains.append(domain_quality)
                evidence_types.add(source)
                total_confidence += confidence
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
            average_confidence = total_confidence / max(1, total_domains)
            evidence_diversity = len(evidence_types)

            # Calculate quality scores
            coverage_score = min(1.0, coverage_percentage / 80.0)  # 80% coverage = perfect
            classification_score = classified_count / max(1, total_domains)
            confidence_score = min(1.0, average_confidence)

            # Overall score (weighted combination)
            overall_score = (
                coverage_score * 0.4 +
                classification_score * 0.3 +
                confidence_score * 0.3
            )

            # Determine tier
            tier = self._determine_tier(overall_score, coverage_percentage, classified_count, total_domains)

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
                average_confidence=average_confidence,
                evidence_diversity=evidence_diversity,
                coverage_score=coverage_score,
                classification_score=classification_score,
                confidence_score=confidence_score,
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
                average_confidence=0,
                evidence_diversity=0,
                tier="failed",
                issues=["parse_error"]
            )

    def _calculate_evidence_score(self, source: str, confidence: float) -> float:
        """Calculate evidence strength score"""
        source_weights = {
            'hhsearch': 1.0,
            'domain_blast': 0.8,
            'chain_blast': 0.6,
            'chain_blast_decomposed': 0.7,  # ADDED: Handle mini's decomposed chains
            'mini_pyecod': 0.9,
            'unknown': 0.3
        }

        weight = source_weights.get(source, 0.5)
        return weight * confidence

    def _determine_tier(self, overall_score: float, coverage: float,
                       classified: int, total: int) -> str:
        """Determine quality tier"""

        if overall_score >= 0.8 and coverage >= 60 and classified >= total * 0.8:
            return "excellent"
        elif overall_score >= 0.6 and coverage >= 40 and classified >= total * 0.5:
            return "good"
        elif overall_score >= 0.4 and coverage >= 20:
            return "acceptable"
        else:
            return "poor"

    def _identify_issues(self, domains: List[DomainQuality], coverage: float,
                        sequence_length: int) -> List[str]:
        """Identify potential issues with the result - FIXED overlap detection"""
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

        # FIXED: Check for overlapping domains using SequenceRange
        if len(domains) > 1:
            try:
                # Parse all domain ranges
                domain_ranges = []
                for domain in domains:
                    try:
                        seq_range = SequenceRange.parse(domain.range)
                        domain_ranges.append(seq_range)
                    except Exception:
                        continue  # Skip unparseable ranges

                # Check for overlaps
                for i in range(len(domain_ranges)):
                    for j in range(i + 1, len(domain_ranges)):
                        if domain_ranges[i].overlaps(domain_ranges[j]):
                            issues.append("overlapping_domains")
                            break
                    if "overlapping_domains" in issues:
                        break

            except Exception as e:
                logger.warning(f"Could not check domain overlaps: {e}")
        
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
            avg_confidence = sum(r.average_confidence for r in successful_results) / len(successful_results)
        else:
            avg_coverage = avg_domains = avg_classified = avg_confidence = 0
        
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
                "average_confidence": avg_confidence
            },
            "issue_analysis": issue_counts,
            "recommendations": self._generate_recommendations(tier_counts, issue_counts, total_results)
        }
    
    def _generate_recommendations(self, tier_counts: Dict, issue_counts: Dict, total: int) -> List[str]:
        """Generate recommendations based on quality analysis"""
        recommendations = []
        
        production_ready = tier_counts.get("excellent", 0) + tier_counts.get("good", 0)
        production_rate = production_ready / max(1, total) * 100
        
        if production_rate >= 70:
            recommendations.append("EXCELLENT: >70% production ready - proceed with full import")
        elif production_rate >= 50:
            recommendations.append("GOOD: 50-70% production ready - import excellent+good tiers")
        elif production_rate >= 30:
            recommendations.append("MODERATE: 30-50% production ready - import excellent tier only")
        else:
            recommendations.append("POOR: <30% production ready - investigate algorithm issues")
        
        # Issue-specific recommendations
        if issue_counts.get("no_domains", 0) > total * 0.2:
            recommendations.append("HIGH: >20% proteins have no domains - check input data quality")
        
        if issue_counts.get("low_coverage", 0) > total * 0.3:
            recommendations.append("COVERAGE: >30% have low coverage - consider coverage thresholds")
        
        if issue_counts.get("no_classification", 0) > total * 0.4:
            recommendations.append("CLASSIFICATION: >40% lack classification - check reference databases")
        
        return recommendations
    
    def export_production_ready(self, results: List[ProteinQuality], 
                               output_file: str = "production_ready.json") -> Dict:
        """Export high-quality results for production import"""
        
        production_results = [
            r for r in results 
            if r.tier in ["excellent", "good"] and "parse_error" not in r.issues
        ]
        
        # Convert to serializable format
        export_data = {
            "export_timestamp": datetime.now().isoformat(),
            "total_assessed": len(results),
            "production_ready": len(production_results),
            "quality_threshold": "excellent + good tiers",
            "proteins": []
        }
        
        for result in production_results:
            export_data["proteins"].append({
                "protein_id": result.protein_id,
                "pdb_id": result.pdb_id,
                "chain_id": result.chain_id,
                "batch_name": result.batch_name,
                "tier": result.tier,
                "overall_score": result.overall_score,
                "coverage_percentage": result.coverage_percentage,
                "total_domains": result.total_domains,
                "classified_domains": result.classified_domains,
                "average_confidence": result.average_confidence,
                "issues": result.issues
            })
        
        # Write to file
        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)
        
        logger.info(f"Exported {len(production_results)} production-ready results to {output_file}")
        return export_data
    
    def print_quality_summary(self, results: List[ProteinQuality]):
        """Print quality assessment summary"""
        
        report = self.generate_quality_report(results)
        
        print("🔍 Mini PyECOD Quality Assessment")
        print("=" * 60)
        print(f"Assessment Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
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
        print(f"  Avg confidence:       {metrics['average_confidence']:>7.1f}")
        print()
        
        # Top issues
        print("⚠️  Top Issues:")
        issue_counts = report["issue_analysis"]
        sorted_issues = sorted(issue_counts.items(), key=lambda x: -x[1])[:5]
        for issue, count in sorted_issues:
            pct = count / max(1, summary['total_results']) * 100
            print(f"  {issue:<20} {count:>6,} ({pct:>5.1f}%)")
        print()
        
        # Recommendations
        print("💡 Recommendations:")
        for rec in report["recommendations"]:
            print(f"  • {rec}")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Assess Quality of Mini PyECOD Results'
    )
    
    parser.add_argument('--scan-all', action='store_true',
                       help='Assess all available results')
    parser.add_argument('--batch-name', type=str,
                       help='Assess specific batch')
    parser.add_argument('--export-production', action='store_true',
                       help='Export production-ready results')
    parser.add_argument('--output-file', type=str, default='production_ready.json',
                       help='Output file for export')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize assessment
    assessor = QualityAssessment(args.config)
    
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
    
    # Export production ready if requested
    if args.export_production:
        export_data = assessor.export_production_ready(results, args.output_file)
        print(f"\n✅ Exported {len(export_data['proteins'])} production-ready results")


if __name__ == "__main__":
    main()
