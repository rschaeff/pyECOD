#!/usr/bin/env python3
"""
Compare Mini vs Main PyECOD Results

Comparative analysis between mini results and existing main results
to demonstrate improvement and validate mini algorithm performance.

Usage:
    python compare_results.py --sample-comparison 100    # Compare 100 proteins
    python compare_results.py --batch-comparison batch_031   # Compare full batch  
    python compare_results.py --win-loss-analysis        # Show wins/losses/ties
"""

import argparse
import yaml
import psycopg2
import psycopg2.extras
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class ComparisonResult:
    """Comparison between mini and main results for one protein"""
    protein_id: str
    pdb_id: str
    chain_id: str
    
    # Mini results
    mini_domains: int
    mini_classified: int
    mini_coverage: float
    mini_confidence: float
    mini_tier: str
    
    # Main results  
    main_domains: int
    main_classified: int
    main_coverage: float
    main_confidence: float
    main_tier: str
    
    # Comparison
    winner: str  # "mini", "main", "tie", "both_poor"
    improvement_type: str  # "coverage", "classification", "confidence", "domains", "overall"
    score_difference: float
    
    # Details
    mini_issues: List[str]
    main_issues: List[str]
    notes: str = ""

class ResultComparator:
    """Compare mini vs main PyECOD results"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        self.config = self._load_config(config_path)
        self.db_conn = self._init_db_connection()
        
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration"""
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _init_db_connection(self):
        """Initialize database connection"""
        try:
            return psycopg2.connect(**self.config['database'])
        except Exception as e:
            logger.error(f"Database connection failed: {e}")
            raise
    
    def get_main_results(self, pdb_id: str, chain_id: str) -> Optional[Dict]:
        """Get main PyECOD results from database"""
        
        try:
            with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                # Get protein record (excluding mini results)
                cursor.execute("""
                    SELECT pp.*, 
                           COUNT(pd.id) as actual_domain_count,
                           COALESCE(AVG(pd.confidence), 0) as avg_confidence,
                           COUNT(pd.id) FILTER (WHERE pd.t_group IS NOT NULL) as classified_count
                    FROM pdb_analysis.partition_proteins pp
                    LEFT JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                    WHERE pp.pdb_id = %s AND pp.chain_id = %s 
                    AND (pp.process_version IS NULL OR pp.process_version != 'mini_pyecod_1.0')
                    GROUP BY pp.id
                    ORDER BY pp.timestamp DESC
                    LIMIT 1
                """, (pdb_id, chain_id))
                
                result = cursor.fetchone()
                
                if result:
                    return {
                        'domains': result['actual_domain_count'],
                        'classified': result['classified_count'],
                        'coverage': result['coverage'] or 0,
                        'confidence': result['avg_confidence'] or 0,
                        'is_classified': result['is_classified'],
                        'process_version': result['process_version'] or 'unknown'
                    }
                return None
                
        except Exception as e:
            logger.error(f"Error getting main results for {pdb_id}_{chain_id}: {e}")
            return None
    
    def parse_mini_result(self, xml_file: Path) -> Dict:
        """Parse mini result from XML file"""
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            domains = root.findall(".//domain")
            classified_count = 0
            total_confidence = 0.0
            total_coverage = 0
            
            for domain in domains:
                if domain.get('t_group'):
                    classified_count += 1
                
                confidence = float(domain.get('confidence', '0'))
                total_confidence += confidence
                
                # Calculate coverage (simple approach)
                start = int(domain.get('start', '0'))
                end = int(domain.get('end', '0'))
                total_coverage += (end - start + 1)
            
            avg_confidence = total_confidence / max(1, len(domains))
            
            return {
                'domains': len(domains),
                'classified': classified_count,
                'coverage': total_coverage,  # Raw residue count
                'confidence': avg_confidence,
                'is_classified': len(domains) > 0
            }
            
        except Exception as e:
            logger.error(f"Error parsing mini result {xml_file}: {e}")
            return {
                'domains': 0,
                'classified': 0,
                'coverage': 0,
                'confidence': 0,
                'is_classified': False
            }
    
    def compare_protein(self, protein_id: str, mini_xml_path: Path) -> Optional[ComparisonResult]:
        """Compare mini vs main results for one protein"""
        
        pdb_id = protein_id.split('_')[0]
        chain_id = protein_id.split('_')[1] if '_' in protein_id else 'A'
        
        # Get results
        mini_result = self.parse_mini_result(mini_xml_path)
        main_result = self.get_main_results(pdb_id, chain_id)
        
        if not main_result:
            # No main result to compare against
            return None
        
        # Determine tiers
        mini_tier = self._determine_tier(mini_result)
        main_tier = self._determine_tier(main_result)
        
        # Calculate scores (simple scoring system)
        mini_score = self._calculate_score(mini_result)
        main_score = self._calculate_score(main_result)
        
        # Determine winner
        winner, improvement_type = self._determine_winner(mini_result, main_result, mini_score, main_score)
        
        return ComparisonResult(
            protein_id=protein_id,
            pdb_id=pdb_id,
            chain_id=chain_id,
            mini_domains=mini_result['domains'],
            mini_classified=mini_result['classified'],
            mini_coverage=mini_result['coverage'],
            mini_confidence=mini_result['confidence'],
            mini_tier=mini_tier,
            main_domains=main_result['domains'],
            main_classified=main_result['classified'],
            main_coverage=main_result['coverage'],
            main_confidence=main_result['confidence'],
            main_tier=main_tier,
            winner=winner,
            improvement_type=improvement_type,
            score_difference=mini_score - main_score,
            mini_issues=[],  # Could add issue detection
            main_issues=[]
        )
    
    def _determine_tier(self, result: Dict) -> str:
        """Determine quality tier for a result"""
        if result['domains'] == 0:
            return "no_domains"
        elif result['classified'] == 0:
            return "unclassified"
        elif result['classified'] >= result['domains'] * 0.8 and result['confidence'] > 0.7:
            return "excellent"
        elif result['classified'] >= result['domains'] * 0.5 and result['confidence'] > 0.5:
            return "good"
        else:
            return "poor"
    
    def _calculate_score(self, result: Dict) -> float:
        """Calculate overall quality score"""
        if result['domains'] == 0:
            return 0.0
        
        classification_rate = result['classified'] / result['domains']
        confidence_score = min(1.0, result['confidence'])
        domain_bonus = min(1.0, result['domains'] / 3.0)  # Bonus for having domains
        
        return (classification_rate * 0.5 + confidence_score * 0.3 + domain_bonus * 0.2)
    
    def _determine_winner(self, mini: Dict, main: Dict, mini_score: float, main_score: float) -> Tuple[str, str]:
        """Determine which result is better and why"""
        
        score_diff = mini_score - main_score
        
        if abs(score_diff) < 0.1:
            return "tie", "similar_quality"
        elif score_diff > 0:
            # Mini wins - determine why
            if mini['classified'] > main['classified']:
                return "mini", "classification"
            elif mini['confidence'] > main['confidence']:
                return "mini", "confidence"
            elif mini['domains'] > main['domains']:
                return "mini", "domains"
            else:
                return "mini", "overall"
        else:
            # Main wins
            if main['classified'] > mini['classified']:
                return "main", "classification"
            elif main['confidence'] > mini['confidence']:
                return "main", "confidence"
            elif main['domains'] > mini['domains']:
                return "main", "domains"
            else:
                return "main", "overall"
    
    def compare_batch(self, batch_name: str) -> List[ComparisonResult]:
        """Compare all proteins in a batch"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        mini_domains_dir = batch_base / batch_name / "mini_domains"
        
        if not mini_domains_dir.exists():
            logger.error(f"No mini results found for batch {batch_name}")
            return []
        
        xml_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
        logger.info(f"Comparing {len(xml_files)} proteins in batch {batch_name}")
        
        comparisons = []
        for xml_file in xml_files:
            protein_id = xml_file.stem.replace('.mini.domains', '')
            comparison = self.compare_protein(protein_id, xml_file)
            if comparison:
                comparisons.append(comparison)
        
        logger.info(f"Successfully compared {len(comparisons)} proteins")
        return comparisons
    
    def sample_comparison(self, sample_size: int = 100) -> List[ComparisonResult]:
        """Compare a random sample of proteins across all batches"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        all_xml_files = []
        
        # Collect all mini result files
        for batch_dir in batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
            mini_domains_dir = batch_dir / "mini_domains"
            if mini_domains_dir.exists():
                xml_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
                all_xml_files.extend(xml_files)
        
        # Sample random files
        import random
        sample_files = random.sample(all_xml_files, min(sample_size, len(all_xml_files)))
        
        logger.info(f"Comparing random sample of {len(sample_files)} proteins")
        
        comparisons = []
        for xml_file in sample_files:
            protein_id = xml_file.stem.replace('.mini.domains', '')
            comparison = self.compare_protein(protein_id, xml_file)
            if comparison:
                comparisons.append(comparison)
        
        return comparisons
    
    def generate_comparison_report(self, comparisons: List[ComparisonResult]) -> Dict:
        """Generate comprehensive comparison report"""
        
        if not comparisons:
            return {"error": "No comparisons to analyze"}
        
        total = len(comparisons)
        
        # Win/loss statistics
        winners = {}
        improvements = {}
        
        for comp in comparisons:
            winners[comp.winner] = winners.get(comp.winner, 0) + 1
            improvements[comp.improvement_type] = improvements.get(comp.improvement_type, 0) + 1
        
        # Quality tier changes
        tier_changes = {}
        for comp in comparisons:
            change = f"{comp.main_tier} -> {comp.mini_tier}"
            tier_changes[change] = tier_changes.get(change, 0) + 1
        
        # Score differences
        score_diffs = [comp.score_difference for comp in comparisons]
        avg_improvement = sum(score_diffs) / len(score_diffs)
        
        # Significant improvements (score diff > 0.2)
        significant_improvements = [c for c in comparisons if c.score_difference > 0.2]
        significant_degradations = [c for c in comparisons if c.score_difference < -0.2]
        
        return {
            "summary": {
                "total_comparisons": total,
                "mini_wins": winners.get("mini", 0),
                "main_wins": winners.get("main", 0), 
                "ties": winners.get("tie", 0),
                "mini_win_rate": winners.get("mini", 0) / total * 100,
                "average_improvement": avg_improvement,
                "significant_improvements": len(significant_improvements),
                "significant_degradations": len(significant_degradations)
            },
            "win_loss_breakdown": winners,
            "improvement_types": improvements,
            "tier_changes": tier_changes,
            "top_improvements": sorted(comparisons, key=lambda x: x.score_difference, reverse=True)[:10],
            "top_degradations": sorted(comparisons, key=lambda x: x.score_difference)[:10]
        }
    
    def print_comparison_summary(self, comparisons: List[ComparisonResult]):
        """Print comparison analysis summary"""
        
        report = self.generate_comparison_report(comparisons)
        
        print("⚔️  Mini vs Main PyECOD Comparison")
        print("=" * 60)
        
        summary = report["summary"]
        print("📊 Overall Results:")
        print(f"  Total comparisons:      {summary['total_comparisons']:>8,}")
        print(f"  Mini wins:              {summary['mini_wins']:>8,} ({summary['mini_win_rate']:>5.1f}%)")
        print(f"  Main wins:              {summary['main_wins']:>8,}")
        print(f"  Ties:                   {summary['ties']:>8,}")
        print(f"  Average improvement:    {summary['average_improvement']:>+8.3f}")
        print(f"  Significant improvements: {summary['significant_improvements']:>6,}")
        print(f"  Significant degradations: {summary['significant_degradations']:>6,}")
        print()
        
        # Improvement types
        print("🎯 Improvement Types:")
        for imp_type, count in sorted(report["improvement_types"].items(), key=lambda x: -x[1]):
            pct = count / summary['total_comparisons'] * 100
            print(f"  {imp_type:<15} {count:>6,} ({pct:>5.1f}%)")
        print()
        
        # Top improvements
        print("🏆 Top 5 Improvements (Mini over Main):")
        for i, comp in enumerate(report["top_improvements"][:5]):
            print(f"  {i+1}. {comp.protein_id:<12} +{comp.score_difference:.3f} ({comp.improvement_type})")
            print(f"     Mini: {comp.mini_domains}d, {comp.mini_classified}c, {comp.mini_confidence:.2f}conf")
            print(f"     Main: {comp.main_domains}d, {comp.main_classified}c, {comp.main_confidence:.2f}conf")
        print()
        
        # Assessment
        win_rate = summary['mini_win_rate']
        if win_rate >= 70:
            print("✅ EXCELLENT: Mini algorithm significantly outperforms main")
        elif win_rate >= 55:
            print("✅ GOOD: Mini algorithm outperforms main")
        elif win_rate >= 45:
            print("⚖️  MIXED: Mini and main perform similarly")
        else:
            print("❌ POOR: Main algorithm outperforms mini")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Compare Mini vs Main PyECOD Results'
    )
    
    parser.add_argument('--sample-comparison', type=int, metavar='N',
                       help='Compare random sample of N proteins')
    parser.add_argument('--batch-comparison', type=str, metavar='BATCH',
                       help='Compare all proteins in specified batch')
    parser.add_argument('--win-loss-analysis', action='store_true',
                       help='Show detailed win/loss analysis')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize comparator
    comparator = ResultComparator(args.config)
    
    # Perform comparison
    if args.sample_comparison:
        comparisons = comparator.sample_comparison(args.sample_comparison)
    elif args.batch_comparison:
        comparisons = comparator.compare_batch(args.batch_comparison)
    else:
        parser.print_help()
        return
    
    if not comparisons:
        print("No comparisons found.")
        return
    
    # Print results
    comparator.print_comparison_summary(comparisons)


if __name__ == "__main__":
    main()
