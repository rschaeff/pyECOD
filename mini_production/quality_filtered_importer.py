#!/usr/bin/env python3
"""
Quality-Filtered Mini PyECOD Results Importer with Sequence Propagation

Imports only high-quality mini results (Excellent/Good tier) and propagates
classifications to redundant protein chains using sequence MD5 matching.

Key features:
- Quality tier filtering (only import excellent/good results)
- Collision-safe import with comparative analysis support
- Sequence-based propagation to redundant chains
- Comprehensive reporting and validation

Usage:
    python quality_filtered_importer.py --assess-and-import --tier-filter excellent,good
    python quality_filtered_importer.py --propagate-sequences --dry-run
    python quality_filtered_importer.py --comparative-analysis
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

# Import your existing quality assessment
sys.path.append(str(Path(__file__).parent))
from assess_quality_v2 import EnhancedQualityAssessment, ProteinQuality
from import_results import MiniResultsImporter, MiniProteinResult

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class PropagationCandidate:
    """A protein chain that could receive propagated classification"""
    protein_id: int
    source_id: str
    pdb_id: str
    chain_id: str
    sequence_md5: str
    sequence_length: int
    has_existing_classification: bool
    donor_protein_id: str  # The mini result that would donate classification

@dataclass
class PropagationResult:
    """Result of sequence-based propagation"""
    candidate: PropagationCandidate
    domains_propagated: int
    success: bool
    error_message: Optional[str] = None

class QualityFilteredImporter(MiniResultsImporter):
    """Enhanced importer with quality filtering and sequence propagation"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        super().__init__(config_path)
        self.quality_assessor = EnhancedQualityAssessment(config_path)
        
    def assess_and_filter_results(self, tier_filter: List[str] = None) -> Dict[str, List[ProteinQuality]]:
        """Assess all results and filter by quality tier
        
        Args:
            tier_filter: List of tiers to include (e.g., ['excellent', 'good'])
                        If None, includes all tiers
        
        Returns:
            Dict mapping batch_name to list of filtered quality results
        """
        
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
        
        # Print tier breakdown
        tier_counts = defaultdict(int)
        for result in all_results:
            tier_counts[result.tier] += 1
        
        logger.info("📊 Tier distribution:")
        for tier, count in sorted(tier_counts.items()):
            status = "✓" if tier in tier_filter else "✗"
            logger.info(f"  {status} {tier:<12} {count:>6} ({count/total_assessed*100:5.1f}%)")
        
        return dict(filtered_by_batch)
    
    def import_quality_filtered_results(self, tier_filter: List[str] = None, 
                                      limit: Optional[int] = None,
                                      collision_strategy: str = "separate") -> Dict[str, any]:
        """Import only high-quality results based on tier filtering
        
        Args:
            tier_filter: Quality tiers to include (default: ['excellent', 'good'])
            limit: Maximum results to import
            collision_strategy: How to handle collisions (default: 'separate')
        """
        
        if tier_filter is None:
            tier_filter = ['excellent', 'good']
        
        logger.info(f"🚀 Quality-filtered import: tiers {tier_filter}, strategy '{collision_strategy}'")
        
        # Get filtered results
        filtered_results = self.assess_and_filter_results(tier_filter)
        
        # Import each filtered result
        stats = {
            'total_assessed': 0,
            'quality_filtered': 0,
            'imported': 0,
            'failed': 0,
            'skipped': 0,
            'by_batch': {},
            'by_tier': defaultdict(int)
        }
        
        imported_count = 0
        
        for batch_name, quality_results in filtered_results.items():
            logger.info(f"Processing batch {batch_name}: {len(quality_results)} quality results")
            
            batch_stats = {'imported': 0, 'failed': 0, 'total': len(quality_results)}
            
            for quality_result in quality_results:
                # Check limit
                if limit and imported_count >= limit:
                    break
                
                try:
                    # Parse the XML file
                    xml_file = Path(self.config["paths"]["batch_base_dir"]) / batch_name / "mini_domains" / f"{quality_result.protein_id}.mini.domains.xml"
                    
                    if not xml_file.exists():
                        logger.warning(f"XML file not found: {xml_file}")
                        batch_stats['failed'] += 1
                        continue
                    
                    mini_result = self.parse_mini_xml(xml_file)
                    
                    # Import with collision handling
                    success = self.import_protein_result(mini_result, collision_strategy)
                    
                    if success:
                        batch_stats['imported'] += 1
                        stats['imported'] += 1
                        stats['by_tier'][quality_result.tier] += 1
                        imported_count += 1
                        
                        logger.info(f"✓ Imported {quality_result.protein_id} (tier: {quality_result.tier}, "
                                  f"score: {quality_result.overall_score:.3f}, domains: {quality_result.total_domains})")
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
        
        logger.info(f"🎉 Quality-filtered import complete:")
        logger.info(f"  Quality filtered: {stats['quality_filtered']}")
        logger.info(f"  Successfully imported: {stats['imported']}")
        logger.info(f"  Failed: {stats['failed']}")
        
        return stats
    
    def find_propagation_candidates(self, limit_per_sequence: int = 10) -> List[PropagationCandidate]:
        """Find protein chains that could receive propagated classifications
        
        Args:
            limit_per_sequence: Maximum chains per unique sequence (to control redundancy)
        
        Returns:
            List of propagation candidates
        """
        
        logger.info(f"🔍 Finding propagation candidates (max {limit_per_sequence} per sequence)")
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Find all imported mini results with their sequence MD5s
            cursor.execute("""
                SELECT DISTINCT
                    pp.pdb_id as donor_pdb_id,
                    pp.chain_id as donor_chain_id,
                    ps.sequence_md5,
                    ps.protein_id as donor_protein_db_id,
                    pp.domains_with_evidence,
                    pp.is_classified
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
                JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                WHERE pp.process_version = 'mini_pyecod_1.0'
                  AND pp.is_classified = true
                  AND pp.domains_with_evidence > 0
                  AND ps.sequence_md5 IS NOT NULL
            """)
            
            mini_donors = cursor.fetchall()
            
            if not mini_donors:
                logger.warning("No mini donor sequences found")
                return []
            
            logger.info(f"Found {len(mini_donors)} mini donor sequences")
            
            # For each donor sequence, find other proteins with same MD5 that lack classifications
            candidates = []
            
            for donor in mini_donors:
                donor_id = f"{donor['donor_pdb_id']}_{donor['donor_chain_id']}"
                
                # Find proteins with same sequence MD5 but no classification
                cursor.execute("""
                    SELECT 
                        p.id as protein_id,
                        p.source_id,
                        p.pdb_id,
                        p.chain_id,
                        ps.sequence_md5,
                        p.length as sequence_length,
                        CASE WHEN pp.id IS NOT NULL THEN true ELSE false END as has_existing_classification
                    FROM pdb_analysis.protein p
                    JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                    LEFT JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id AND p.chain_id = pp.chain_id
                    WHERE ps.sequence_md5 = %s
                      AND p.source_id != %s  -- Exclude the donor itself
                      AND (pp.id IS NULL OR pp.is_classified = false)  -- No existing classification or poor classification
                    ORDER BY p.id
                    LIMIT %s
                """, (donor['sequence_md5'], donor_id, limit_per_sequence))
                
                sequence_candidates = cursor.fetchall()
                
                for candidate_row in sequence_candidates:
                    candidate = PropagationCandidate(
                        protein_id=candidate_row['protein_id'],
                        source_id=candidate_row['source_id'],
                        pdb_id=candidate_row['pdb_id'],
                        chain_id=candidate_row['chain_id'],
                        sequence_md5=candidate_row['sequence_md5'],
                        sequence_length=candidate_row['sequence_length'],
                        has_existing_classification=candidate_row['has_existing_classification'],
                        donor_protein_id=donor_id
                    )
                    candidates.append(candidate)
        
        logger.info(f"✓ Found {len(candidates)} propagation candidates")
        
        # Group by sequence for reporting
        by_sequence = defaultdict(list)
        for candidate in candidates:
            by_sequence[candidate.sequence_md5].append(candidate)
        
        logger.info(f"📊 Propagation potential: {len(by_sequence)} unique sequences, "
                   f"avg {len(candidates)/len(by_sequence):.1f} targets per sequence")
        
        return candidates
    
    def propagate_sequence_classifications(self, candidates: List[PropagationCandidate],
                                         dry_run: bool = True) -> List[PropagationResult]:
        """Propagate classifications to sequence-identical proteins
        
        Args:
            candidates: List of propagation candidates
            dry_run: If True, don't actually write to database
        
        Returns:
            List of propagation results
        """
        
        action = "DRY RUN" if dry_run else "EXECUTING"
        logger.info(f"🔄 Propagating classifications to {len(candidates)} candidates ({action})")
        
        results = []
        
        for candidate in candidates:
            try:
                result = self._propagate_single_candidate(candidate, dry_run)
                results.append(result)
                
                if result.success:
                    logger.info(f"✓ {'[DRY RUN] ' if dry_run else ''}Propagated {result.domains_propagated} domains to {candidate.source_id}")
                else:
                    logger.warning(f"✗ Failed to propagate to {candidate.source_id}: {result.error_message}")
                    
            except Exception as e:
                error_result = PropagationResult(
                    candidate=candidate,
                    domains_propagated=0,
                    success=False,
                    error_message=str(e)
                )
                results.append(error_result)
                logger.error(f"Error propagating to {candidate.source_id}: {e}")
        
        # Summary
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        total_domains = sum(r.domains_propagated for r in successful)
        
        logger.info(f"🎉 Propagation complete ({action}):")
        logger.info(f"  Successful: {len(successful)}")
        logger.info(f"  Failed: {len(failed)}")
        logger.info(f"  Total domains propagated: {total_domains}")
        
        return results
    
    def _propagate_single_candidate(self, candidate: PropagationCandidate, 
                                   dry_run: bool) -> PropagationResult:
        """Propagate classification to a single candidate"""
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Get the donor's classification data
            donor_pdb_id, donor_chain_id = candidate.donor_protein_id.split('_', 1)
            
            cursor.execute("""
                SELECT pp.*, 
                       array_agg(
                           json_build_object(
                               'domain_number', pd.domain_number,
                               'start_pos', pd.start_pos,
                               'end_pos', pd.end_pos,
                               'range', pd.range,
                               'source', pd.source,
                               'source_id', pd.source_id,
                               'confidence', pd.confidence,
                               't_group', pd.t_group,
                               'h_group', pd.h_group,
                               'x_group', pd.x_group,
                               'a_group', pd.a_group,
                               'is_manual_rep', pd.is_manual_rep,
                               'is_f70', pd.is_f70,
                               'is_f40', pd.is_f40,
                               'is_f99', pd.is_f99
                           ) ORDER BY pd.domain_number
                       ) as domains
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.pdb_id = %s AND pp.chain_id = %s 
                  AND pp.process_version = 'mini_pyecod_1.0'
                GROUP BY pp.id, pp.pdb_id, pp.chain_id, pp.batch_id, pp.reference_version,
                         pp.is_classified, pp.sequence_length, pp.coverage, pp.residues_assigned,
                         pp.domains_with_evidence, pp.fully_classified_domains, pp.process_version,
                         pp.timestamp
            """, (donor_pdb_id, donor_chain_id))
            
            donor_data = cursor.fetchone()
            
            if not donor_data:
                return PropagationResult(
                    candidate=candidate,
                    domains_propagated=0,
                    success=False,
                    error_message="Donor classification not found"
                )
            
            if dry_run:
                # Just count domains for dry run
                domains_count = len(donor_data['domains'])
                return PropagationResult(
                    candidate=candidate,
                    domains_propagated=domains_count,
                    success=True
                )
            
            # Create partition_proteins record for candidate
            cursor.execute("""
                INSERT INTO pdb_analysis.partition_proteins 
                (pdb_id, chain_id, batch_id, reference_version, is_classified, 
                 sequence_length, coverage, residues_assigned, domains_with_evidence, 
                 fully_classified_domains, process_version)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                RETURNING id
            """, (
                candidate.pdb_id,
                candidate.chain_id,
                donor_data['batch_id'],
                donor_data['reference_version'],
                donor_data['is_classified'],
                candidate.sequence_length,
                donor_data['coverage'],
                donor_data['residues_assigned'],
                donor_data['domains_with_evidence'],
                donor_data['fully_classified_domains'],
                'mini_pyecod_propagated_1.0'  # Distinguish propagated results
            ))
            
            new_protein_id = cursor.fetchone()[0]
            
            # Create domain records
            domains_created = 0
            for domain_data in donor_data['domains']:
                cursor.execute("""
                    INSERT INTO pdb_analysis.partition_domains
                    (protein_id, domain_number, domain_id, start_pos, end_pos, range,
                     source, source_id, confidence, t_group, h_group, x_group, a_group,
                     is_manual_rep, is_f70, is_f40, is_f99)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                """, (
                    new_protein_id,
                    domain_data['domain_number'],
                    f"{candidate.source_id}_d{domain_data['domain_number']}_propagated",
                    domain_data['start_pos'],
                    domain_data['end_pos'],
                    domain_data['range'],
                    'mini_pyecod_propagated',  # Mark as propagated
                    domain_data['source_id'],
                    domain_data['confidence'],
                    domain_data['t_group'],
                    domain_data['h_group'],
                    domain_data['x_group'],
                    domain_data['a_group'],
                    domain_data['is_manual_rep'],
                    domain_data['is_f70'],
                    domain_data['is_f40'],
                    domain_data['is_f99']
                ))
                domains_created += 1
            
            self.db_conn.commit()
            
            return PropagationResult(
                candidate=candidate,
                domains_propagated=domains_created,
                success=True
            )
    
    def generate_comparative_analysis(self) -> Dict[str, any]:
        """Generate comparative analysis between original, mini, and propagated results"""
        
        logger.info("📊 Generating comparative analysis...")
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Overall statistics by process version
            cursor.execute("""
                SELECT 
                    COALESCE(process_version, 'original') as source_type,
                    COUNT(*) as protein_count,
                    COUNT(*) FILTER (WHERE is_classified = true) as classified_count,
                    AVG(domains_with_evidence) as avg_domains,
                    SUM(domains_with_evidence) as total_domains
                FROM pdb_analysis.partition_proteins
                GROUP BY process_version
                ORDER BY 
                    CASE process_version 
                        WHEN 'mini_pyecod_1.0' THEN 1
                        WHEN 'mini_pyecod_propagated_1.0' THEN 2
                        ELSE 0
                    END, process_version
            """)
            
            overall_stats = cursor.fetchall()
            
            # Classification quality comparison
            cursor.execute("""
                SELECT 
                    COALESCE(pp.process_version, 'original') as source_type,
                    COUNT(DISTINCT pd.t_group) as unique_t_groups,
                    COUNT(DISTINCT pd.h_group) as unique_h_groups,
                    AVG(pd.confidence) as avg_confidence,
                    COUNT(*) FILTER (WHERE pd.t_group IS NOT NULL) as classified_domains,
                    COUNT(*) as total_domains
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                GROUP BY pp.process_version
                ORDER BY 
                    CASE pp.process_version 
                        WHEN 'mini_pyecod_1.0' THEN 1
                        WHEN 'mini_pyecod_propagated_1.0' THEN 2
                        ELSE 0
                    END, pp.process_version
            """)
            
            classification_stats = cursor.fetchall()
            
            # Sequence propagation effectiveness
            cursor.execute("""
                SELECT 
                    COUNT(DISTINCT ps.sequence_md5) as total_unique_sequences,
                    COUNT(DISTINCT 
                        CASE WHEN pp.process_version = 'mini_pyecod_1.0' 
                             THEN ps.sequence_md5 END
                    ) as mini_sequences,
                    COUNT(DISTINCT 
                        CASE WHEN pp.process_version = 'mini_pyecod_propagated_1.0' 
                             THEN ps.sequence_md5 END
                    ) as propagated_sequences
                FROM pdb_analysis.protein p
                JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                LEFT JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id AND p.chain_id = pp.chain_id
                WHERE ps.sequence_md5 IS NOT NULL
            """)
            
            propagation_stats = cursor.fetchone()
            
            # Top families discovered
            cursor.execute("""
                SELECT 
                    pd.t_group,
                    COUNT(*) as domain_count,
                    COUNT(DISTINCT pp.id) as protein_count,
                    STRING_AGG(DISTINCT pp.process_version, ', ') as found_in_sources
                FROM pdb_analysis.partition_domains pd
                JOIN pdb_analysis.partition_proteins pp ON pd.protein_id = pp.id
                WHERE pd.t_group IS NOT NULL
                  AND pp.process_version IN ('mini_pyecod_1.0', 'mini_pyecod_propagated_1.0')
                GROUP BY pd.t_group
                ORDER BY COUNT(*) DESC
                LIMIT 20
            """)
            
            top_families = cursor.fetchall()
        
        analysis = {
            'overall_statistics': [dict(row) for row in overall_stats],
            'classification_quality': [dict(row) for row in classification_stats],
            'propagation_effectiveness': dict(propagation_stats) if propagation_stats else {},
            'top_families_discovered': [dict(row) for row in top_families],
            'analysis_timestamp': datetime.now().isoformat()
        }
        
        return analysis
    
    def print_comparative_analysis(self, analysis: Dict[str, any]):
        """Print formatted comparative analysis report"""
        
        print("\n🔍 Mini PyECOD Comparative Analysis Report")
        print("=" * 80)
        print(f"Generated: {analysis['analysis_timestamp']}")
        print()
        
        # Overall Statistics
        print("📊 Overall Statistics by Source:")
        print(f"{'Source Type':<30} {'Proteins':<10} {'Classified':<12} {'Avg Domains':<12} {'Total Domains':<15}")
        print("-" * 85)
        
        for stat in analysis['overall_statistics']:
            source = stat['source_type'] or 'original'
            classified_pct = stat['classified_count'] / max(1, stat['protein_count']) * 100
            avg_domains = stat['avg_domains'] or 0
            
            print(f"{source:<30} {stat['protein_count']:<10,} "
                  f"{stat['classified_count']:<7,} ({classified_pct:4.1f}%) "
                  f"{avg_domains:<12.1f} {stat['total_domains']:<15,}")
        print()
        
        # Classification Quality
        print("🎯 Classification Quality Comparison:")
        print(f"{'Source Type':<30} {'T-Groups':<10} {'H-Groups':<10} {'Confidence':<12} {'Class Rate':<12}")
        print("-" * 80)
        
        for stat in analysis['classification_quality']:
            source = stat['source_type'] or 'original'
            class_rate = stat['classified_domains'] / max(1, stat['total_domains']) * 100
            confidence = stat['avg_confidence'] or 0
            
            print(f"{source:<30} {stat['unique_t_groups']:<10} {stat['unique_h_groups']:<10} "
                  f"{confidence:<12.3f} {class_rate:<11.1f}%")
        print()
        
        # Propagation Effectiveness
        prop = analysis['propagation_effectiveness']
        if prop:
            print("🔄 Sequence Propagation Effectiveness:")
            print(f"  Total unique sequences:       {prop['total_unique_sequences']:>8,}")
            print(f"  Mini classified sequences:    {prop['mini_sequences']:>8,}")
            print(f"  Propagated sequences:         {prop['propagated_sequences']:>8,}")
            
            if prop['mini_sequences'] > 0:
                propagation_factor = prop['propagated_sequences'] / prop['mini_sequences']
                coverage_improvement = (prop['mini_sequences'] + prop['propagated_sequences']) / prop['total_unique_sequences'] * 100
                print(f"  Propagation factor:           {propagation_factor:>8.1f}x")
                print(f"  Total sequence coverage:      {coverage_improvement:>7.1f}%")
            print()
        
        # Top Families
        if analysis['top_families_discovered']:
            print("🏆 Top Families Discovered (by Mini + Propagation):")
            print(f"{'T-Group':<15} {'Domains':<8} {'Proteins':<10} {'Sources'}")
            print("-" * 60)
            
            for family in analysis['top_families_discovered'][:15]:
                print(f"{family['t_group']:<15} {family['domain_count']:<8} "
                      f"{family['protein_count']:<10} {family['found_in_sources']}")
            print()
        
        # Summary Assessment
        mini_stats = next((s for s in analysis['overall_statistics'] if s['source_type'] == 'mini_pyecod_1.0'), {})
        prop_stats = next((s for s in analysis['overall_statistics'] if s['source_type'] == 'mini_pyecod_propagated_1.0'), {})
        
        mini_proteins = mini_stats.get('protein_count', 0)
        prop_proteins = prop_stats.get('protein_count', 0)
        total_new_classifications = mini_proteins + prop_proteins
        
        print("🎉 Impact Summary:")
        print(f"  Direct mini classifications:   {mini_proteins:>8,}")
        print(f"  Propagated classifications:    {prop_proteins:>8,}")
        print(f"  Total new classifications:     {total_new_classifications:>8,}")
        
        if mini_proteins > 0 and prop_proteins > 0:
            amplification = (mini_proteins + prop_proteins) / mini_proteins
            print(f"  Classification amplification:  {amplification:>7.1f}x")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Quality-Filtered Mini PyECOD Importer with Sequence Propagation'
    )
    
    parser.add_argument('--assess-and-import', action='store_true',
                       help='Assess quality and import filtered results')
    parser.add_argument('--tier-filter', type=str, default='excellent,good',
                       help='Comma-separated quality tiers to import (default: excellent,good)')
    parser.add_argument('--limit', type=int,
                       help='Maximum results to import')
    parser.add_argument('--collision-strategy', type=str, default='separate',
                       choices=['separate', 'skip', 'update'],
                       help='Collision handling strategy (default: separate)')
    
    parser.add_argument('--propagate-sequences', action='store_true',
                       help='Propagate classifications to sequence-identical chains')
    parser.add_argument('--limit-per-sequence', type=int, default=10,
                       help='Maximum propagation targets per sequence (default: 10)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Dry run - analyze propagation without writing to database')
    
    parser.add_argument('--comparative-analysis', action='store_true',
                       help='Generate comparative analysis report')
    
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize enhanced importer
    importer = QualityFilteredImporter(args.config)
    
    if args.assess_and_import:
        tier_filter = [tier.strip() for tier in args.tier_filter.split(',')]
        stats = importer.import_quality_filtered_results(
            tier_filter=tier_filter,
            limit=args.limit,
            collision_strategy=args.collision_strategy
        )
        
        print(f"\n✅ Quality-Filtered Import Results:")
        print(f"  Tier filter: {tier_filter}")
        print(f"  Quality filtered: {stats['quality_filtered']}")
        print(f"  Successfully imported: {stats['imported']}")
        print(f"  Failed: {stats['failed']}")
        print(f"\n  By tier:")
        for tier, count in stats['by_tier'].items():
            print(f"    {tier}: {count}")
        return
    
    if args.propagate_sequences:
        candidates = importer.find_propagation_candidates(args.limit_per_sequence)
        
        if not candidates:
            print("No propagation candidates found.")
            return
        
        results = importer.propagate_sequence_classifications(candidates, args.dry_run)
        
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        
        action = "DRY RUN" if args.dry_run else "EXECUTED"
        print(f"\n🔄 Sequence Propagation Results ({action}):")
        print(f"  Candidates processed: {len(candidates)}")
        print(f"  Successful: {len(successful)}")
        print(f"  Failed: {len(failed)}")
        print(f"  Total domains propagated: {sum(r.domains_propagated for r in successful)}")
        
        if not args.dry_run:
            print(f"  ✅ Propagation complete - {len(successful)} new classified protein chains")
        return
    
    if args.comparative_analysis:
        analysis = importer.generate_comparative_analysis()
        importer.print_comparative_analysis(analysis)
        return
    
    # Default: show help
    parser.print_help()


if __name__ == "__main__":
    main()
