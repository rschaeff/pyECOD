#!/usr/bin/env python3
"""
Mini PyECOD v2.0 Data Validation Script

Comprehensive validation of imported mini PyECOD v2.0 data to ensure:
- Data integrity and completeness
- Proper batch relationships
- Evidence quality metrics
- Provenance tracking
- Classification coverage
- Boundary optimization data

Usage:
    python validate_mini_v2_data.py --comprehensive
    python validate_mini_v2_data.py --quick-check
    python validate_mini_v2_data.py --export-report
"""

import psycopg2
import psycopg2.extras
import yaml
import argparse
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Any
import json

class MiniV2DataValidator:
    """Comprehensive validator for mini PyECOD v2.0 imported data"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.db_conn = psycopg2.connect(**self.config["database"])
        self.validation_results = {}
        
    def run_comprehensive_validation(self) -> Dict[str, Any]:
        """Run all validation checks"""
        
        print("🔍 Starting Comprehensive Mini PyECOD v2.0 Data Validation")
        print("=" * 70)
        
        self.validation_results = {
            'timestamp': datetime.now().isoformat(),
            'import_statistics': self.check_import_statistics(),
            'batch_relationships': self.check_batch_relationships(),
            'data_integrity': self.check_data_integrity(),
            'evidence_quality': self.check_evidence_quality(),
            'provenance_tracking': self.check_provenance_tracking(),
            'classification_coverage': self.check_classification_coverage(),
            'boundary_optimization': self.check_boundary_optimization(),
            'version_comparison': self.check_version_comparison(),
            'propagation_effectiveness': self.check_propagation_effectiveness(),
            'potential_issues': self.identify_potential_issues()
        }
        
        return self.validation_results
    
    def check_import_statistics(self) -> Dict[str, Any]:
        """Check basic import statistics"""
        
        print("\n📊 Import Statistics")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Overall v2.0 statistics
            cursor.execute("""
                SELECT 
                    COUNT(*) as total_proteins,
                    COUNT(*) FILTER (WHERE is_classified = true) as classified_proteins,
                    SUM(domains_with_evidence) as total_domains,
                    AVG(domains_with_evidence) as avg_domains_per_protein,
                    AVG(coverage) as avg_coverage,
                    MIN(timestamp) as earliest_import,
                    MAX(timestamp) as latest_import
                FROM pdb_analysis.partition_proteins 
                WHERE process_version = 'mini_pyecod_v2.0'
            """)
            
            overall_stats = dict(cursor.fetchone() or {})
            
            # By batch statistics
            cursor.execute("""
                SELECT 
                    b.batch_name,
                    b.type as batch_type,
                    COUNT(pp.*) as protein_count,
                    SUM(pp.domains_with_evidence) as domain_count,
                    AVG(pp.coverage) as avg_coverage,
                    MIN(pp.timestamp) as import_time
                FROM pdb_analysis.partition_proteins pp
                JOIN ecod_schema.batch b ON pp.batch_id = b.id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                GROUP BY b.id, b.batch_name, b.type
                ORDER BY COUNT(pp.*) DESC
            """)
            
            batch_stats = [dict(row) for row in cursor.fetchall()]
            
            # Process version distribution
            cursor.execute("""
                SELECT 
                    process_version,
                    COUNT(*) as protein_count,
                    SUM(domains_with_evidence) as domain_count
                FROM pdb_analysis.partition_proteins 
                WHERE process_version LIKE '%mini%'
                GROUP BY process_version
                ORDER BY COUNT(*) DESC
            """)
            
            version_distribution = [dict(row) for row in cursor.fetchall()]
        
        stats = {
            'overall': overall_stats,
            'by_batch': batch_stats,
            'version_distribution': version_distribution
        }
        
        # Print summary
        if overall_stats:
            print(f"  Total v2.0 proteins: {overall_stats['total_proteins']:,}")
            print(f"  Classified proteins: {overall_stats['classified_proteins']:,} "
                  f"({overall_stats['classified_proteins']/max(1,overall_stats['total_proteins'])*100:.1f}%)")
            print(f"  Total domains: {overall_stats['total_domains']:,}")
            print(f"  Avg domains/protein: {overall_stats['avg_domains_per_protein']:.2f}")
            print(f"  Avg coverage: {overall_stats['avg_coverage']:.3f}")
            print(f"  Import timespan: {overall_stats['earliest_import']} to {overall_stats['latest_import']}")
        
        print(f"  Batches represented: {len(batch_stats)}")
        
        return stats
    
    def check_batch_relationships(self) -> Dict[str, Any]:
        """Validate batch relationships and referential integrity"""
        
        print("\n🔗 Batch Relationships")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Check for orphaned batch references
            cursor.execute("""
                SELECT COUNT(*) as orphaned_count
                FROM pdb_analysis.partition_proteins pp
                LEFT JOIN ecod_schema.batch b ON pp.batch_id = b.id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pp.batch_id IS NOT NULL
                  AND b.id IS NULL
            """)
            
            orphaned_batches = cursor.fetchone()['orphaned_count']
            
            # Check for NULL batch_ids
            cursor.execute("""
                SELECT COUNT(*) as null_batch_count
                FROM pdb_analysis.partition_proteins pp
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pp.batch_id IS NULL
            """)
            
            null_batches = cursor.fetchone()['null_batch_count']
            
            # Check batch names in process_parameters
            cursor.execute("""
                SELECT 
                    COUNT(*) as with_batch_names,
                    COUNT(DISTINCT process_parameters->>'source_batch_name') as unique_batch_names
                FROM pdb_analysis.partition_proteins pp
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pp.process_parameters ? 'source_batch_name'
            """)
            
            batch_names_info = dict(cursor.fetchone())
            
            # Batch type distribution
            cursor.execute("""
                SELECT 
                    b.type,
                    COUNT(pp.*) as protein_count
                FROM pdb_analysis.partition_proteins pp
                JOIN ecod_schema.batch b ON pp.batch_id = b.id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                GROUP BY b.type
                ORDER BY COUNT(pp.*) DESC
            """)
            
            batch_types = [dict(row) for row in cursor.fetchall()]
        
        batch_health = {
            'orphaned_batch_refs': orphaned_batches,
            'null_batch_ids': null_batches,
            'batch_names_preserved': batch_names_info,
            'batch_type_distribution': batch_types,
            'health_status': 'healthy' if orphaned_batches == 0 and null_batches == 0 else 'issues_found'
        }
        
        print(f"  Orphaned batch references: {orphaned_batches}")
        print(f"  NULL batch_ids: {null_batches}")
        print(f"  Batch names preserved: {batch_names_info.get('with_batch_names', 0)}")
        print(f"  Batch health: {batch_health['health_status']}")
        
        return batch_health
    
    def check_data_integrity(self) -> Dict[str, Any]:
        """Check data integrity and consistency"""
        
        print("\n🔍 Data Integrity")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Check for missing protein records
            cursor.execute("""
                SELECT COUNT(*) as missing_protein_records
                FROM pdb_analysis.partition_proteins pp
                LEFT JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND p.id IS NULL
            """)
            
            missing_proteins = cursor.fetchone()['missing_protein_records']
            
            # Check domain count consistency
            cursor.execute("""
                SELECT 
                    COUNT(*) as inconsistent_domain_counts
                FROM pdb_analysis.partition_proteins pp
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pp.domains_with_evidence != (
                      SELECT COUNT(*) 
                      FROM pdb_analysis.partition_domains pd 
                      WHERE pd.protein_id = pp.id
                  )
            """)
            
            inconsistent_domains = cursor.fetchone()['inconsistent_domain_counts']
            
            # Check for domains without evidence
            cursor.execute("""
                SELECT COUNT(*) as domains_without_evidence
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                LEFT JOIN pdb_analysis.domain_evidence de ON pd.id = de.domain_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND de.id IS NULL
            """)
            
            domains_no_evidence = cursor.fetchone()['domains_without_evidence']
            
            # Check range format consistency
            cursor.execute("""
                SELECT COUNT(*) as invalid_ranges
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND (pd.range IS NULL OR pd.range = '' OR NOT pd.range ~ '^\d+(-\d+)?(,\d+(-\d+)?)*$')
            """)
            
            invalid_ranges = cursor.fetchone()['invalid_ranges']
            
            # Check classification consistency
            cursor.execute("""
                SELECT COUNT(*) as classification_inconsistencies
                FROM pdb_analysis.partition_proteins pp
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pp.is_classified = true
                  AND pp.fully_classified_domains = 0
            """)
            
            classification_issues = cursor.fetchone()['classification_inconsistencies']
        
        integrity = {
            'missing_protein_records': missing_proteins,
            'inconsistent_domain_counts': inconsistent_domains,
            'domains_without_evidence': domains_no_evidence,
            'invalid_range_formats': invalid_ranges,
            'classification_inconsistencies': classification_issues,
            'overall_integrity': 'good' if all([
                missing_proteins == 0,
                inconsistent_domains == 0,
                domains_no_evidence == 0,
                invalid_ranges == 0,
                classification_issues == 0
            ]) else 'issues_found'
        }
        
        print(f"  Missing protein records: {missing_proteins}")
        print(f"  Inconsistent domain counts: {inconsistent_domains}")
        print(f"  Domains without evidence: {domains_no_evidence}")
        print(f"  Invalid range formats: {invalid_ranges}")
        print(f"  Classification inconsistencies: {classification_issues}")
        print(f"  Overall integrity: {integrity['overall_integrity']}")
        
        return integrity
    
    def check_evidence_quality(self) -> Dict[str, Any]:
        """Check enhanced evidence quality metrics"""
        
        print("\n🎯 Evidence Quality")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Evidence type distribution
            cursor.execute("""
                SELECT 
                    de.evidence_type,
                    COUNT(*) as evidence_count,
                    AVG(de.confidence) as avg_confidence,
                    COUNT(*) FILTER (WHERE de.evalue IS NOT NULL) as with_evalue,
                    COUNT(*) FILTER (WHERE de.probability IS NOT NULL) as with_probability,
                    AVG(de.evalue) FILTER (WHERE de.evalue IS NOT NULL) as avg_evalue
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                JOIN pdb_analysis.domain_evidence de ON pd.id = de.domain_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                GROUP BY de.evidence_type
                ORDER BY COUNT(*) DESC
            """)
            
            evidence_types = [dict(row) for row in cursor.fetchall()]
            
            # Quality metrics distribution
            cursor.execute("""
                SELECT 
                    COUNT(*) as total_evidence,
                    COUNT(*) FILTER (WHERE de.confidence >= 0.8) as high_confidence,
                    COUNT(*) FILTER (WHERE de.confidence < 0.5) as low_confidence,
                    COUNT(*) FILTER (WHERE de.evalue < 1e-10) as excellent_evalue,
                    COUNT(*) FILTER (WHERE de.evalue > 1.0) as poor_evalue,
                    AVG(de.confidence) as overall_avg_confidence
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                JOIN pdb_analysis.domain_evidence de ON pd.id = de.domain_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
            """)
            
            quality_metrics = dict(cursor.fetchone())
            
            # Primary evidence tracking
            cursor.execute("""
                SELECT 
                    COUNT(*) as domains_with_primary_evidence,
                    COUNT(*) FILTER (WHERE pd.primary_evidence_type IS NOT NULL) as with_evidence_type,
                    COUNT(*) FILTER (WHERE pd.evidence_evalue IS NOT NULL) as with_evidence_evalue
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
            """)
            
            primary_evidence = dict(cursor.fetchone())
        
        evidence_quality = {
            'evidence_type_distribution': evidence_types,
            'quality_metrics': quality_metrics,
            'primary_evidence_tracking': primary_evidence
        }
        
        print(f"  Evidence types: {len(evidence_types)}")
        for ev_type in evidence_types[:3]:  # Show top 3
            print(f"    {ev_type['evidence_type']}: {ev_type['evidence_count']:,} "
                  f"(conf: {ev_type['avg_confidence']:.3f})")
        
        if quality_metrics:
            total = quality_metrics['total_evidence']
            print(f"  High confidence (≥0.8): {quality_metrics['high_confidence']:,}/{total:,} "
                  f"({quality_metrics['high_confidence']/max(1,total)*100:.1f}%)")
            print(f"  Low confidence (<0.5): {quality_metrics['low_confidence']:,}/{total:,} "
                  f"({quality_metrics['low_confidence']/max(1,total)*100:.1f}%)")
        
        return evidence_quality
    
    def check_provenance_tracking(self) -> Dict[str, Any]:
        """Check provenance and metadata tracking"""
        
        print("\n📋 Provenance Tracking")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Git and version tracking
            cursor.execute("""
                SELECT 
                    COUNT(*) as total_proteins,
                    COUNT(*) FILTER (WHERE algorithm_version IS NOT NULL) as with_algorithm_version,
                    COUNT(*) FILTER (WHERE git_commit_hash IS NOT NULL) as with_git_hash,
                    COUNT(*) FILTER (WHERE source_xml_path IS NOT NULL) as with_source_path,
                    COUNT(*) FILTER (WHERE source_xml_hash IS NOT NULL) as with_source_hash,
                    COUNT(DISTINCT algorithm_version) as unique_algorithm_versions,
                    COUNT(DISTINCT git_commit_hash) as unique_git_commits
                FROM pdb_analysis.partition_proteins pp
                WHERE pp.process_version = 'mini_pyecod_v2.0'
            """)
            
            provenance_stats = dict(cursor.fetchone())
            
            # Process parameters tracking
            cursor.execute("""
                SELECT 
                    COUNT(*) as with_process_parameters,
                    COUNT(*) FILTER (WHERE process_parameters ? 'source_batch_name') as with_batch_names,
                    COUNT(*) FILTER (WHERE process_parameters ? 'import_timestamp') as with_import_timestamps,
                    array_agg(DISTINCT algorithm_version) as algorithm_versions,
                    array_agg(DISTINCT git_commit_hash) as git_commits
                FROM pdb_analysis.partition_proteins pp
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND process_parameters IS NOT NULL
            """)
            
            params_stats = dict(cursor.fetchone())
            
            # File integrity (source XML hashes)
            cursor.execute("""
                SELECT 
                    COUNT(DISTINCT source_xml_hash) as unique_source_files,
                    COUNT(*) FILTER (WHERE source_xml_hash IS NOT NULL) as with_file_hashes
                FROM pdb_analysis.partition_proteins pp
                WHERE pp.process_version = 'mini_pyecod_v2.0'
            """)
            
            file_integrity = dict(cursor.fetchone())
        
        provenance = {
            'version_tracking': provenance_stats,
            'process_parameters': params_stats,
            'file_integrity': file_integrity,
            'provenance_completeness': provenance_stats['with_git_hash'] / max(1, provenance_stats['total_proteins'])
        }
        
        total = provenance_stats['total_proteins']
        print(f"  Algorithm versions: {provenance_stats['with_algorithm_version']:,}/{total:,}")
        print(f"  Git commit hashes: {provenance_stats['with_git_hash']:,}/{total:,}")
        print(f"  Source file hashes: {provenance_stats['with_source_hash']:,}/{total:,}")
        print(f"  Unique algorithms: {provenance_stats['unique_algorithm_versions']}")
        print(f"  Unique git commits: {provenance_stats['unique_git_commits']}")
        print(f"  Provenance completeness: {provenance['provenance_completeness']:.1%}")
        
        return provenance
    
    def check_classification_coverage(self) -> Dict[str, Any]:
        """Check classification coverage and family distribution"""
        
        print("\n🏷️ Classification Coverage")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # T-group distribution
            cursor.execute("""
                SELECT 
                    pd.t_group,
                    COUNT(*) as domain_count,
                    COUNT(DISTINCT pp.id) as protein_count
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pd.t_group IS NOT NULL
                GROUP BY pd.t_group
                ORDER BY COUNT(*) DESC
                LIMIT 20
            """)
            
            top_tgroups = [dict(row) for row in cursor.fetchall()]
            
            # Classification completeness
            cursor.execute("""
                SELECT 
                    COUNT(*) as total_domains,
                    COUNT(*) FILTER (WHERE pd.t_group IS NOT NULL) as with_tgroup,
                    COUNT(*) FILTER (WHERE pd.h_group IS NOT NULL) as with_hgroup,
                    COUNT(*) FILTER (WHERE pd.x_group IS NOT NULL) as with_xgroup,
                    COUNT(DISTINCT pd.t_group) as unique_tgroups,
                    COUNT(DISTINCT pd.h_group) as unique_hgroups
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
            """)
            
            classification_stats = dict(cursor.fetchone())
            
            # Confidence distribution
            cursor.execute("""
                SELECT 
                    CASE 
                        WHEN pd.confidence >= 0.8 THEN 'high'
                        WHEN pd.confidence >= 0.5 THEN 'medium'
                        WHEN pd.confidence > 0 THEN 'low'
                        ELSE 'none'
                    END as confidence_tier,
                    COUNT(*) as domain_count
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                GROUP BY 1
                ORDER BY 2 DESC
            """)
            
            confidence_distribution = [dict(row) for row in cursor.fetchall()]
        
        classification = {
            'top_families': top_tgroups,
            'classification_completeness': classification_stats,
            'confidence_distribution': confidence_distribution
        }
        
        if classification_stats:
            total = classification_stats['total_domains']
            print(f"  Total domains: {total:,}")
            print(f"  With T-group: {classification_stats['with_tgroup']:,} "
                  f"({classification_stats['with_tgroup']/max(1,total)*100:.1f}%)")
            print(f"  Unique T-groups: {classification_stats['unique_tgroups']:,}")
            print(f"  Unique H-groups: {classification_stats['unique_hgroups']:,}")
        
        print(f"  Top families: {len(top_tgroups)} shown")
        
        return classification
    
    def check_boundary_optimization(self) -> Dict[str, Any]:
        """Check boundary optimization data"""
        
        print("\n🔧 Boundary Optimization")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Optimization statistics
            cursor.execute("""
                SELECT 
                    COUNT(*) as total_domains,
                    COUNT(*) FILTER (WHERE pd.original_range IS NOT NULL) as with_original_range,
                    COUNT(*) FILTER (WHERE pd.optimization_actions IS NOT NULL) as with_optimization_actions,
                    COUNT(*) FILTER (WHERE pd.original_range IS NOT NULL AND pd.original_range != pd.range) as actually_optimized
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
            """)
            
            optimization_stats = dict(cursor.fetchone())
            
            # Optimization actions distribution
            cursor.execute("""
                SELECT 
                    unnest(pd.optimization_actions) as action,
                    COUNT(*) as frequency
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pd.optimization_actions IS NOT NULL
                GROUP BY 1
                ORDER BY 2 DESC
            """)
            
            optimization_actions = [dict(row) for row in cursor.fetchall()]
        
        optimization = {
            'optimization_statistics': optimization_stats,
            'action_distribution': optimization_actions
        }
        
        if optimization_stats:
            total = optimization_stats['total_domains']
            print(f"  Total domains: {total:,}")
            print(f"  With original range: {optimization_stats['with_original_range']:,}")
            print(f"  With optimization actions: {optimization_stats['with_optimization_actions']:,}")
            print(f"  Actually optimized: {optimization_stats['actually_optimized']:,}")
            
            if optimization_stats['with_original_range'] > 0:
                opt_rate = optimization_stats['actually_optimized'] / optimization_stats['with_original_range'] * 100
                print(f"  Optimization rate: {opt_rate:.1f}%")
        
        return optimization
    
    def check_version_comparison(self) -> Dict[str, Any]:
        """Compare v2.0 with other mini versions"""
        
        print("\n📈 Version Comparison")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Compare mini versions
            cursor.execute("""
                SELECT 
                    process_version,
                    COUNT(*) as protein_count,
                    SUM(domains_with_evidence) as domain_count,
                    AVG(coverage) as avg_coverage,
                    COUNT(*) FILTER (WHERE is_classified = true) as classified_count
                FROM pdb_analysis.partition_proteins 
                WHERE process_version LIKE '%mini%'
                GROUP BY process_version
                ORDER BY 
                    CASE 
                        WHEN process_version = 'mini_pyecod_v2.0' THEN 1
                        WHEN process_version = 'mini_pyecod_v2.0_propagated' THEN 2
                        ELSE 3
                    END,
                    process_version
            """)
            
            version_comparison = [dict(row) for row in cursor.fetchall()]
        
        comparison = {
            'version_statistics': version_comparison
        }
        
        for version in version_comparison:
            print(f"  {version['process_version']}: {version['protein_count']:,} proteins, "
                  f"{version['domain_count']:,} domains")
        
        return comparison
    
    def check_propagation_effectiveness(self) -> Dict[str, Any]:
        """Check propagation effectiveness (only meaningful post-propagation)"""
        
        print("\n🔄 Propagation Effectiveness")
        print("-" * 30)
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Check if propagated data exists
            cursor.execute("""
                SELECT COUNT(*) as propagated_count
                FROM pdb_analysis.partition_proteins 
                WHERE process_version = 'mini_pyecod_v2.0_propagated'
            """)
            
            propagated_count = cursor.fetchone()['propagated_count']
            
            if propagated_count == 0:
                print("  No propagated data found - run pre-propagation validation")
                return {'status': 'no_propagated_data'}
            
            # Propagation amplification factor
            cursor.execute("""
                SELECT 
                    COUNT(*) FILTER (WHERE process_version = 'mini_pyecod_v2.0') as original_count,
                    COUNT(*) FILTER (WHERE process_version = 'mini_pyecod_v2.0_propagated') as propagated_count,
                    COUNT(DISTINCT ps.sequence_md5) as unique_sequences_covered
                FROM pdb_analysis.partition_proteins pp
                LEFT JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
                LEFT JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                WHERE pp.process_version IN ('mini_pyecod_v2.0', 'mini_pyecod_v2.0_propagated')
            """)
            
            amplification = dict(cursor.fetchone())
            
            # Quality comparison between original and propagated
            cursor.execute("""
                SELECT 
                    pp.process_version,
                    AVG(pd.confidence) as avg_confidence,
                    COUNT(*) FILTER (WHERE pd.confidence >= 0.8) as high_confidence_domains,
                    COUNT(*) as total_domains
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version IN ('mini_pyecod_v2.0', 'mini_pyecod_v2.0_propagated')
                GROUP BY pp.process_version
            """)
            
            quality_comparison = [dict(row) for row in cursor.fetchall()]
            
            # Sequence propagation coverage
            cursor.execute("""
                SELECT 
                    COUNT(DISTINCT ps.sequence_md5) as total_unique_sequences,
                    COUNT(DISTINCT CASE WHEN pp.process_version = 'mini_pyecod_v2.0' THEN ps.sequence_md5 END) as original_sequences,
                    COUNT(DISTINCT CASE WHEN pp.process_version = 'mini_pyecod_v2.0_propagated' THEN ps.sequence_md5 END) as propagated_sequences
                FROM pdb_analysis.protein p
                JOIN pdb_analysis.protein_sequence ps ON p.id = ps.protein_id
                LEFT JOIN pdb_analysis.partition_proteins pp ON p.pdb_id = pp.pdb_id AND p.chain_id = pp.chain_id
                WHERE pp.process_version IN ('mini_pyecod_v2.0', 'mini_pyecod_v2.0_propagated')
                  AND ps.sequence_md5 IS NOT NULL
            """)
            
            sequence_coverage = dict(cursor.fetchone())
        
        propagation = {
            'status': 'propagated_data_found',
            'amplification_metrics': amplification,
            'quality_comparison': quality_comparison,
            'sequence_coverage': sequence_coverage
        }
        
        # Print results
        if amplification:
            original = amplification['original_count']
            propagated = amplification['propagated_count']
            total = original + propagated
            amplification_factor = propagated / max(1, original)
            
            print(f"  Original proteins: {original:,}")
            print(f"  Propagated proteins: {propagated:,}")
            print(f"  Total proteins: {total:,}")
            print(f"  Amplification factor: {amplification_factor:.1f}x")
        
        if sequence_coverage:
            seq_total = sequence_coverage['total_unique_sequences']
            seq_original = sequence_coverage['original_sequences']
            seq_propagated = sequence_coverage['propagated_sequences']
            coverage_pct = (seq_original + seq_propagated) / max(1, seq_total) * 100
            
            print(f"  Unique sequences covered: {seq_original + seq_propagated:,}/{seq_total:,} ({coverage_pct:.1f}%)")
        
        return propagation
    
    def identify_potential_issues(self) -> List[Dict[str, Any]]:
        """Identify potential data issues that need attention"""
        
        print("\n⚠️ Potential Issues")
        print("-" * 30)
        
        issues = []
        
        with self.db_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            # Check for domains with very low confidence
            cursor.execute("""
                SELECT COUNT(*) as count
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pd.confidence < 0.3
            """)
            
            low_confidence = cursor.fetchone()['count']
            if low_confidence > 0:
                issues.append({
                    'type': 'low_confidence_domains',
                    'count': low_confidence,
                    'description': f'{low_confidence} domains with confidence < 0.3'
                })
            
            # Check for very large domains (potential issues)
            cursor.execute("""
                SELECT COUNT(*) as count
                FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pd.length > 1000
            """)
            
            large_domains = cursor.fetchone()['count']
            if large_domains > 0:
                issues.append({
                    'type': 'unusually_large_domains',
                    'count': large_domains,
                    'description': f'{large_domains} domains with length > 1000 residues'
                })
            
            # Check for proteins with very high coverage (>1.0)
            cursor.execute("""
                SELECT COUNT(*) as count
                FROM pdb_analysis.partition_proteins pp
                WHERE pp.process_version = 'mini_pyecod_v2.0'
                  AND pp.coverage > 1.0
            """)
            
            high_coverage = cursor.fetchone()['count']
            if high_coverage > 0:
                issues.append({
                    'type': 'coverage_over_100_percent',
                    'count': high_coverage,
                    'description': f'{high_coverage} proteins with coverage > 100%'
                })
        
        if not issues:
            print("  No significant issues detected ✓")
        else:
            for issue in issues:
                print(f"  {issue['type']}: {issue['description']}")
        
        return issues
    
    def print_validation_summary(self):
        """Print a comprehensive validation summary"""
        
        if not self.validation_results:
            print("No validation results available. Run validation first.")
            return
        
        print(f"\n🎯 Validation Summary")
        print("=" * 50)
        
        # Overall health
        overall_health = "HEALTHY"
        health_issues = []
        
        # Check each component
        if self.validation_results['batch_relationships']['health_status'] != 'healthy':
            health_issues.append("Batch relationship issues")
            overall_health = "ISSUES FOUND"
        
        if self.validation_results['data_integrity']['overall_integrity'] != 'good':
            health_issues.append("Data integrity issues")
            overall_health = "ISSUES FOUND"
        
        if self.validation_results['potential_issues']:
            health_issues.append(f"{len(self.validation_results['potential_issues'])} data quality concerns")
            overall_health = "REVIEW NEEDED"
        
        # Print status
        print(f"Overall Status: {overall_health}")
        if health_issues:
            print("Issues Found:")
            for issue in health_issues:
                print(f"  - {issue}")
        
        # Key metrics
        stats = self.validation_results['import_statistics']['overall']
        if stats:
            print(f"\nKey Metrics:")
            print(f"  Imported proteins: {stats['total_proteins']:,}")
            print(f"  Classified proteins: {stats['classified_proteins']:,}")
            print(f"  Total domains: {stats['total_domains']:,}")
            print(f"  Average coverage: {stats['avg_coverage']:.3f}")
        
        provenance = self.validation_results['provenance_tracking']
        print(f"  Provenance completeness: {provenance['provenance_completeness']:.1%}")
        
        print(f"\nValidation completed at: {self.validation_results['timestamp']}")
    
    def export_validation_report(self, output_path: str = "mini_v2_validation_report.json"):
        """Export validation results to JSON file"""
        
        if not self.validation_results:
            print("No validation results to export. Run validation first.")
            return
        
        with open(output_path, 'w') as f:
            json.dump(self.validation_results, f, indent=2, default=str)
        
        print(f"✅ Validation report exported to: {output_path}")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Mini PyECOD v2.0 Data Validation',
        epilog='''
Recommended workflow:
  1. python validate_mini_v2_data.py --pre-propagation
  2. Run propagation: SELECT * FROM pdb_analysis.propagate_mini_v2_classifications(false);
  3. python validate_mini_v2_data.py --post-propagation
        '''
    )
    
    parser.add_argument('--comprehensive', action='store_true',
                       help='Run comprehensive validation (default)')
    parser.add_argument('--quick-check', action='store_true',
                       help='Run quick validation check')
    parser.add_argument('--pre-propagation', action='store_true',
                       help='Pre-propagation validation (imported data only)')
    parser.add_argument('--post-propagation', action='store_true',
                       help='Post-propagation validation (includes propagation analysis)')
    parser.add_argument('--export-report', action='store_true',
                       help='Export validation report to JSON')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    parser.add_argument('--output', type=str, default='mini_v2_validation_report.json',
                       help='Output file for validation report')
    
    args = parser.parse_args()
    
    # Initialize validator
    validator = MiniV2DataValidator(args.config)
    
    if args.pre_propagation:
        print("🔍 PRE-PROPAGATION VALIDATION")
        print("Validating imported mini PyECOD v2.0 data before propagation")
        print("=" * 70)
        
        # Run validation but skip propagation effectiveness
        validator.validation_results = {
            'timestamp': datetime.now().isoformat(),
            'validation_type': 'pre_propagation',
            'import_statistics': validator.check_import_statistics(),
            'batch_relationships': validator.check_batch_relationships(),
            'data_integrity': validator.check_data_integrity(),
            'evidence_quality': validator.check_evidence_quality(),
            'provenance_tracking': validator.check_provenance_tracking(),
            'classification_coverage': validator.check_classification_coverage(),
            'boundary_optimization': validator.check_boundary_optimization(),
            'version_comparison': validator.check_version_comparison(),
            'potential_issues': validator.identify_potential_issues()
        }
        validator.print_validation_summary()
        
        print(f"\n🚀 Next Steps:")
        print(f"  1. If validation passed, run propagation:")
        print(f"     SELECT * FROM pdb_analysis.propagate_mini_v2_classifications(false);")
        print(f"  2. Then run post-propagation validation:")
        print(f"     python validate_mini_v2_data.py --post-propagation")
        
    elif args.post_propagation:
        print("🔍 POST-PROPAGATION VALIDATION")
        print("Validating complete dataset after propagation")
        print("=" * 70)
        
        # Run full validation including propagation analysis
        validator.run_comprehensive_validation()
        validator.validation_results['validation_type'] = 'post_propagation'
        validator.print_validation_summary()
        
        # Show propagation effectiveness if available
        if validator.validation_results['propagation_effectiveness']['status'] == 'propagated_data_found':
            print(f"\n🎉 Propagation Analysis Complete!")
            print(f"  Check 'Propagation Effectiveness' section above for amplification metrics")
        else:
            print(f"\n⚠️ No propagated data found. Run propagation first.")
        
    elif args.quick_check:
        # Just run import statistics and summary
        validator.validation_results = {
            'timestamp': datetime.now().isoformat(),
            'import_statistics': validator.check_import_statistics(),
            'batch_relationships': validator.check_batch_relationships(),
            'data_integrity': validator.check_data_integrity()
        }
        validator.print_validation_summary()
        
    else:
        # Run comprehensive validation (default)
        validator.run_comprehensive_validation()
        validator.print_validation_summary()
    
    if args.export_report:
        validator.export_validation_report(args.output)


if __name__ == "__main__":
    main()
