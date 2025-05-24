#!/usr/bin/env python3
"""
Identify Problematic Batches Script

This script analyzes batch completion status and identifies batches that need repair
based on various criteria like missing partitions, low success rates, and gaps.

Usage:
    python identify_problematic_batches.py --config config.yml [options]

Options:
    --config CONFIG              Path to configuration file
    --output-format FORMAT       Output format: table, json, csv (default: table)
    --severity-threshold LEVEL   Minimum severity to report: low, medium, high, critical (default: medium)
    --include-non-rep            Include non-representative proteins in analysis
    --export-repair-list FILE    Export list of batch IDs needing repair
    --auto-repair                Automatically trigger repair for critical batches
    --verbose                    Enable verbose output
"""

import os
import sys
import logging
import argparse
import json
import csv
import psycopg2
from psycopg2.extras import RealDictCursor
import yaml
from datetime import datetime
from typing import List, Dict, Any, Optional, Set
from enum import Enum
from dataclasses import dataclass


class Severity(Enum):
    """Severity levels for batch issues."""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


@dataclass
class BatchIssue:
    """Represents an issue found with a batch."""
    batch_id: int
    batch_name: str
    issue_type: str
    severity: Severity
    description: str
    metric_value: float
    recommendation: str
    auto_fixable: bool = False


def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format_str)
    return logging.getLogger(__name__)


def parse_config(config_path):
    """Parse configuration file with optional local overrides."""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        # Check for local config
        local_config_path = os.path.join(
            os.path.dirname(config_path),
            'config.local.yml'
        )

        if os.path.exists(local_config_path):
            with open(local_config_path, 'r') as f:
                local_config = yaml.safe_load(f)

            if local_config:
                config = deep_merge(config, local_config)

        return config
    except Exception as e:
        logging.error(f"Error parsing config file: {str(e)}")
        raise


def deep_merge(dict1, dict2):
    """Deep merge two dictionaries."""
    result = dict1.copy()
    for key, value in dict2.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def get_db_connection(config):
    """Create database connection from config."""
    db_config = config.get('database', {})

    try:
        conn = psycopg2.connect(
            host=db_config.get('host', 'dione'),
            port=db_config.get('port', 45000),
            dbname=db_config.get('name', 'ecod_protein'),
            user=db_config.get('user', 'ecod'),
            password=db_config.get('password', '')
        )
        return conn
    except psycopg2.Error as e:
        logging.error(f"Database connection error: {e}")
        raise


def get_batch_audit_data(conn, include_non_rep=False):
    """Get comprehensive batch audit data similar to the API endpoint."""
    
    rep_filter = "" if include_non_rep else "AND ps.is_representative = true"
    
    query = f"""
    WITH batch_overview AS (
        SELECT
            b.id as batch_id,
            b.batch_name,
            b.type as batch_type,
            b.total_items as reported_total,
            b.completed_items as reported_completed,
            b.status as batch_status,
            b.ref_version,
            b.created_at as batch_created
        FROM ecod_schema.batch b
        WHERE b.type IN ('pdb_hhsearch', 'domain_analysis')
    ),
    batch_proteins AS (
        SELECT
            ps.batch_id,
            COUNT(*) as proteins_in_batch,
            COUNT(CASE WHEN ps.current_stage IN ('completed', 'classified', 'domain_partition_complete') THEN 1 END) as proteins_reported_done,
            COUNT(CASE WHEN ps.is_representative = true THEN 1 END) as representative_count,
            COUNT(CASE WHEN ps.is_representative = false THEN 1 END) as non_representative_count,
            COUNT(CASE WHEN ps.status = 'error' THEN 1 END) as error_count,
            array_agg(DISTINCT ps.current_stage) as stages_present,
            array_agg(DISTINCT ps.status) as statuses_present
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
        WHERE 1=1 {rep_filter}
        GROUP BY ps.batch_id
    ),
    batch_protein_list AS (
        SELECT
            ps.batch_id,
            ep.source_id,
            ep.pdb_id,
            ep.chain_id,
            ps.is_representative
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein ep ON ps.protein_id = ep.id
        WHERE 1=1 {rep_filter}
    ),
    partition_results AS (
        SELECT
            bpl.batch_id,
            COUNT(pp.id) as partitions_attempted,
            COUNT(CASE WHEN pp.is_classified = true THEN 1 END) as partitions_classified,
            COUNT(CASE WHEN pp.is_classified = false THEN 1 END) as partitions_unclassified,
            COUNT(CASE WHEN pp.sequence_length < 30 THEN 1 END) as partitions_peptide,
            COUNT(pd.id) as total_domains_found,
            COUNT(de.id) as total_evidence_items,
            AVG(pp.coverage) as avg_coverage,
            AVG(CASE WHEN pp.is_classified THEN pp.coverage END) as avg_classified_coverage
        FROM batch_protein_list bpl
        LEFT JOIN pdb_analysis.partition_proteins pp ON (
            bpl.source_id = (pp.pdb_id || '_' || pp.chain_id)
            AND pp.batch_id = bpl.batch_id
        )
        LEFT JOIN pdb_analysis.partition_domains pd ON pp.id = pd.protein_id
        LEFT JOIN pdb_analysis.domain_evidence de ON pd.id = de.domain_id
        GROUP BY bpl.batch_id
    ),
    missing_analysis AS (
        SELECT
            bpl.batch_id,
            COUNT(*) as proteins_missing_partitions
        FROM batch_protein_list bpl
        LEFT JOIN pdb_analysis.partition_proteins pp ON (
            bpl.source_id = (pp.pdb_id || '_' || pp.chain_id)
            AND pp.batch_id = bpl.batch_id
        )
        WHERE pp.id IS NULL
        GROUP BY bpl.batch_id
    ),
    file_analysis AS (
        SELECT
            ps.batch_id,
            COUNT(CASE WHEN pf.file_type = 'fasta' AND pf.file_exists = true THEN 1 END) as fasta_files_exist,
            COUNT(CASE WHEN pf.file_type LIKE '%blast%' AND pf.file_exists = true THEN 1 END) as blast_files_exist,
            COUNT(CASE WHEN pf.file_type LIKE '%hhsearch%' AND pf.file_exists = true THEN 1 END) as hhsearch_files_exist,
            COUNT(CASE WHEN pf.file_type LIKE '%partition%' AND pf.file_exists = true THEN 1 END) as partition_files_exist,
            COUNT(CASE WHEN pf.file_type = 'domain_summary' AND pf.file_exists = true THEN 1 END) as summary_files_exist
        FROM ecod_schema.process_status ps
        LEFT JOIN ecod_schema.process_file pf ON ps.id = pf.process_id
        WHERE 1=1 {rep_filter}
        GROUP BY ps.batch_id
    )
    SELECT
        bo.batch_id,
        bo.batch_name,
        bo.batch_type,
        bo.ref_version,
        bo.batch_status,
        bo.batch_created,
        
        -- Expectations vs Reality
        bo.reported_total as batch_reported_total,
        bo.reported_completed as batch_reported_completed,
        COALESCE(bp.proteins_in_batch, 0) as actual_proteins_in_batch,
        COALESCE(bp.proteins_reported_done, 0) as proteins_reported_done,
        COALESCE(bp.representative_count, 0) as representative_count,
        COALESCE(bp.non_representative_count, 0) as non_representative_count,
        COALESCE(bp.error_count, 0) as error_count,
        
        -- Partition Results
        COALESCE(pr.partitions_attempted, 0) as partitions_attempted,
        COALESCE(pr.partitions_classified, 0) as partitions_classified,
        COALESCE(pr.partitions_unclassified, 0) as partitions_unclassified,
        COALESCE(pr.partitions_peptide, 0) as partitions_peptide,
        COALESCE(pr.total_domains_found, 0) as total_domains_found,
        COALESCE(pr.total_evidence_items, 0) as total_evidence_items,
        COALESCE(pr.avg_coverage, 0) as avg_coverage,
        COALESCE(pr.avg_classified_coverage, 0) as avg_classified_coverage,
        
        -- Gaps and Issues
        COALESCE(ma.proteins_missing_partitions, 0) as proteins_missing_partitions,
        (COALESCE(bp.proteins_in_batch, 0) - COALESCE(pr.partitions_attempted, 0)) as partition_gap,
        (bo.reported_total - COALESCE(bp.proteins_in_batch, 0)) as batch_definition_gap,
        
        -- File Status
        COALESCE(fa.fasta_files_exist, 0) as fasta_files_exist,
        COALESCE(fa.blast_files_exist, 0) as blast_files_exist,
        COALESCE(fa.hhsearch_files_exist, 0) as hhsearch_files_exist,
        COALESCE(fa.partition_files_exist, 0) as partition_files_exist,
        COALESCE(fa.summary_files_exist, 0) as summary_files_exist,
        
        -- Status Arrays
        COALESCE(bp.stages_present, ARRAY[]::text[]) as stages_present,
        COALESCE(bp.statuses_present, ARRAY[]::text[]) as statuses_present,
        
        -- Calculated Quality Metrics
        CASE
            WHEN COALESCE(bp.proteins_in_batch, 0) > 0
            THEN ROUND((COALESCE(pr.partitions_attempted, 0)::numeric / bp.proteins_in_batch::numeric) * 100, 1)
            ELSE 0.0
        END as partition_attempt_rate,
        
        CASE
            WHEN COALESCE(pr.partitions_attempted, 0) > 0
            THEN ROUND((COALESCE(pr.partitions_classified, 0)::numeric / pr.partitions_attempted::numeric) * 100, 1)
            ELSE 0.0
        END as classification_success_rate,
        
        CASE
            WHEN COALESCE(bp.proteins_in_batch, 0) > 0
            THEN ROUND((COALESCE(pr.partitions_classified, 0)::numeric / bp.proteins_in_batch::numeric) * 100, 1)
            ELSE 0.0
        END as overall_success_rate
        
    FROM batch_overview bo
    LEFT JOIN batch_proteins bp ON bo.batch_id = bp.batch_id
    LEFT JOIN partition_results pr ON bo.batch_id = pr.batch_id
    LEFT JOIN missing_analysis ma ON bo.batch_id = ma.batch_id
    LEFT JOIN file_analysis fa ON bo.batch_id = fa.batch_id
    ORDER BY
        ma.proteins_missing_partitions DESC NULLS LAST,
        (COALESCE(bp.proteins_in_batch, 0) - COALESCE(pr.partitions_attempted, 0)) DESC,
        bo.batch_id DESC
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query)
        return cur.fetchall()


def analyze_batch_issues(batch_data: Dict[str, Any]) -> List[BatchIssue]:
    """Analyze a single batch and identify issues."""
    issues = []
    batch_id = batch_data['batch_id']
    batch_name = batch_data['batch_name'] or f"Batch {batch_id}"
    
    # Issue 1: Missing partitions
    missing_partitions = batch_data['proteins_missing_partitions']
    total_proteins = batch_data['actual_proteins_in_batch']
    
    if missing_partitions > 0 and total_proteins > 0:
        missing_rate = (missing_partitions / total_proteins) * 100
        
        if missing_rate >= 90:
            severity = Severity.CRITICAL
            description = f"{missing_partitions}/{total_proteins} proteins missing partitions ({missing_rate:.1f}%)"
            recommendation = "Batch may need complete reprocessing from domain summary stage"
            auto_fixable = True
        elif missing_rate >= 50:
            severity = Severity.HIGH
            description = f"{missing_partitions}/{total_proteins} proteins missing partitions ({missing_rate:.1f}%)"
            recommendation = "Run import script in repair mode, then reprocess missing proteins"
            auto_fixable = True
        elif missing_rate >= 20:
            severity = Severity.MEDIUM
            description = f"{missing_partitions}/{total_proteins} proteins missing partitions ({missing_rate:.1f}%)"
            recommendation = "Run targeted repair for missing proteins"
            auto_fixable = True
        else:
            severity = Severity.LOW
            description = f"{missing_partitions}/{total_proteins} proteins missing partitions ({missing_rate:.1f}%)"
            recommendation = "Consider running import script to capture missed partitions"
            auto_fixable = True
        
        issues.append(BatchIssue(
            batch_id=batch_id,
            batch_name=batch_name,
            issue_type="missing_partitions",
            severity=severity,
            description=description,
            metric_value=missing_rate,
            recommendation=recommendation,
            auto_fixable=auto_fixable
        ))
    
    # Issue 2: Low classification success rate
    classification_rate = batch_data['classification_success_rate']
    partitions_attempted = batch_data['partitions_attempted']
    
    if partitions_attempted > 0:
        if classification_rate < 30:
            severity = Severity.HIGH
            description = f"Very low classification rate: {classification_rate:.1f}% of {partitions_attempted} attempts"
            recommendation = "Check reference database version and reprocess with full pipeline"
        elif classification_rate < 50:
            severity = Severity.MEDIUM
            description = f"Low classification rate: {classification_rate:.1f}% of {partitions_attempted} attempts"
            recommendation = "Review unclassified cases and consider parameter tuning"
        elif classification_rate < 70:
            severity = Severity.LOW
            description = f"Below-average classification rate: {classification_rate:.1f}% of {partitions_attempted} attempts"
            recommendation = "Monitor and compare with similar batches"
        
        if classification_rate < 70:
            issues.append(BatchIssue(
                batch_id=batch_id,
                batch_name=batch_name,
                issue_type="low_classification_rate",
                severity=severity,
                description=description,
                metric_value=classification_rate,
                recommendation=recommendation,
                auto_fixable=False
            ))
    
    # Issue 3: Large partition gap
    partition_gap = batch_data['partition_gap']
    
    if partition_gap > 0:
        gap_rate = (partition_gap / total_proteins) * 100 if total_proteins > 0 else 0
        
        if gap_rate >= 50:
            severity = Severity.HIGH
            description = f"Large gap between expected and actual partitions: {partition_gap} proteins ({gap_rate:.1f}%)"
            recommendation = "Run import script in repair mode to capture existing partition files"
            auto_fixable = True
        elif gap_rate >= 20:
            severity = Severity.MEDIUM
            description = f"Moderate partition gap: {partition_gap} proteins ({gap_rate:.1f}%)"
            recommendation = "Check for unimported partition files"
            auto_fixable=True
        elif gap_rate >= 10:
            severity = Severity.LOW
            description = f"Small partition gap: {partition_gap} proteins ({gap_rate:.1f}%)"
            recommendation = "Review import logs for any missed files"
            auto_fixable=True
            
        if gap_rate >= 10:
            issues.append(BatchIssue(
                batch_id=batch_id,
                batch_name=batch_name,
                issue_type="partition_gap",
                severity=severity,
                description=description,
                metric_value=gap_rate,
                recommendation=recommendation,
                auto_fixable=auto_fixable
            ))
    
    # Issue 4: High error count
    error_count = batch_data['error_count']
    
    if error_count > 0:
        error_rate = (error_count / total_proteins) * 100 if total_proteins > 0 else 0
        
        if error_rate >= 20:
            severity = Severity.HIGH
            description = f"High error rate: {error_count} proteins failed ({error_rate:.1f}%)"
            recommendation = "Review error logs and reset failed processes"
            auto_fixable = False
        elif error_rate >= 10:
            severity = Severity.MEDIUM
            description = f"Moderate error rate: {error_count} proteins failed ({error_rate:.1f}%)"
            recommendation = "Reset failed processes and retry"
            auto_fixable = False
        elif error_rate >= 5:
            severity = Severity.LOW
            description = f"Some processing errors: {error_count} proteins failed ({error_rate:.1f}%)"
            recommendation = "Review individual failure cases"
            auto_fixable = False
            
        if error_rate >= 5:
            issues.append(BatchIssue(
                batch_id=batch_id,
                batch_name=batch_name,
                issue_type="high_error_rate",
                severity=severity,
                description=description,
                metric_value=error_rate,
                recommendation=recommendation,
                auto_fixable=auto_fixable
            ))
    
    # Issue 5: Batch definition mismatch
    definition_gap = batch_data['batch_definition_gap']
    
    if abs(definition_gap) > 10:
        if abs(definition_gap) >= 100:
            severity = Severity.MEDIUM
            description = f"Large mismatch between reported total ({batch_data['batch_reported_total']}) and actual proteins ({total_proteins})"
            recommendation = "Update batch metadata or verify protein assignments"
        else:
            severity = Severity.LOW
            description = f"Small mismatch between reported total ({batch_data['batch_reported_total']}) and actual proteins ({total_proteins})"
            recommendation = "Update batch metadata"
            
        issues.append(BatchIssue(
            batch_id=batch_id,
            batch_name=batch_name,
            issue_type="batch_definition_mismatch",
            severity=severity,
            description=description,
            metric_value=abs(definition_gap),
            recommendation=recommendation,
            auto_fixable=False
        ))
    
    # Issue 6: Missing files
    files_missing = []
    if batch_data['summary_files_exist'] == 0 and total_proteins > 0:
        files_missing.append("domain summaries")
    if batch_data['partition_files_exist'] == 0 and total_proteins > 0:
        files_missing.append("partition files")
    
    if files_missing:
        severity = Severity.HIGH if "partition files" in files_missing else Severity.MEDIUM
        description = f"Missing critical files: {', '.join(files_missing)}"
        recommendation = "Rerun domain analysis pipeline from appropriate stage"
        
        issues.append(BatchIssue(
            batch_id=batch_id,
            batch_name=batch_name,
            issue_type="missing_files",
            severity=severity,
            description=description,
            metric_value=len(files_missing),
            recommendation=recommendation,
            auto_fixable=False
        ))
    
    return issues


def filter_issues_by_severity(issues: List[BatchIssue], min_severity: Severity) -> List[BatchIssue]:
    """Filter issues by minimum severity level."""
    severity_order = {
        Severity.LOW: 1,
        Severity.MEDIUM: 2,
        Severity.HIGH: 3,
        Severity.CRITICAL: 4
    }
    
    min_level = severity_order[min_severity]
    return [issue for issue in issues if severity_order[issue.severity] >= min_level]


def output_table_format(issues: List[BatchIssue], logger):
    """Output issues in table format."""
    if not issues:
        logger.info("No issues found matching the specified criteria.")
        return
    
    # Group by batch
    batch_issues = {}
    for issue in issues:
        if issue.batch_id not in batch_issues:
            batch_issues[issue.batch_id] = []
        batch_issues[issue.batch_id].append(issue)
    
    print(f"\n{'='*120}")
    print(f"PROBLEMATIC BATCHES ANALYSIS - {len(batch_issues)} batches with {len(issues)} total issues")
    print(f"{'='*120}")
    
    for batch_id in sorted(batch_issues.keys()):
        batch_issue_list = batch_issues[batch_id]
        batch_name = batch_issue_list[0].batch_name
        
        print(f"\nBatch {batch_id}: {batch_name}")
        print(f"{'-'*80}")
        
        for issue in sorted(batch_issue_list, key=lambda x: ['low', 'medium', 'high', 'critical'].index(x.severity.value), reverse=True):
            severity_color = {
                'low': '',
                'medium': '⚠️ ',
                'high': '🔥 ',
                'critical': '🚨 '
            }
            
            print(f"  {severity_color.get(issue.severity.value, '')}{issue.severity.value.upper()}: {issue.issue_type}")
            print(f"    {issue.description}")
            print(f"    → {issue.recommendation}")
            if issue.auto_fixable:
                print(f"    ✅ Auto-fixable")
            print()


def output_json_format(issues: List[BatchIssue], logger):
    """Output issues in JSON format."""
    issues_data = []
    for issue in issues:
        issues_data.append({
            'batch_id': issue.batch_id,
            'batch_name': issue.batch_name,
            'issue_type': issue.issue_type,
            'severity': issue.severity.value,
            'description': issue.description,
            'metric_value': issue.metric_value,
            'recommendation': issue.recommendation,
            'auto_fixable': issue.auto_fixable
        })
    
    print(json.dumps({
        'timestamp': datetime.now().isoformat(),
        'total_issues': len(issues),
        'issues': issues_data
    }, indent=2))


def output_csv_format(issues: List[BatchIssue], logger):
    """Output issues in CSV format."""
    import sys
    
    writer = csv.writer(sys.stdout)
    writer.writerow([
        'batch_id', 'batch_name', 'issue_type', 'severity', 
        'description', 'metric_value', 'recommendation', 'auto_fixable'
    ])
    
    for issue in issues:
        writer.writerow([
            issue.batch_id, issue.batch_name, issue.issue_type, issue.severity.value,
            issue.description, issue.metric_value, issue.recommendation, issue.auto_fixable
        ])


def export_repair_list(issues: List[BatchIssue], filename: str, logger):
    """Export list of batch IDs that need repair."""
    # Get unique batch IDs with auto-fixable issues
    repair_batches = set()
    for issue in issues:
        if issue.auto_fixable and issue.severity in [Severity.HIGH, Severity.CRITICAL]:
            repair_batches.add(issue.batch_id)
    
    with open(filename, 'w') as f:
        for batch_id in sorted(repair_batches):
            f.write(f"{batch_id}\n")
    
    logger.info(f"Exported {len(repair_batches)} batch IDs needing repair to {filename}")


def trigger_auto_repair(issues: List[BatchIssue], config, logger):
    """Automatically trigger repair for critical auto-fixable issues."""
    critical_batches = set()
    
    for issue in issues:
        if issue.auto_fixable and issue.severity == Severity.CRITICAL:
            critical_batches.add(issue.batch_id)
    
    if not critical_batches:
        logger.info("No critical auto-fixable issues found")
        return
    
    logger.info(f"Found {len(critical_batches)} batches with critical auto-fixable issues")
    
    # Here you would integrate with your domain partition pipeline
    # For now, just log the commands that should be run
    
    for batch_id in sorted(critical_batches):
        logger.info(f"Batch {batch_id} needs repair:")
        logger.info(f"  1. Run import script: python import_domain_partitions_with_repair.py --config {config} --batch-id {batch_id} --repair-mode --conflict-strategy update")
        logger.info(f"  2. If still missing partitions: python domain_partition_run.py process batch --batch-id {batch_id}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Identify problematic batches needing repair')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--output-format', choices=['table', 'json', 'csv'], default='table',
                       help='Output format')
    parser.add_argument('--severity-threshold', choices=['low', 'medium', 'high', 'critical'], 
                       default='medium', help='Minimum severity to report')
    parser.add_argument('--include-non-rep', action='store_true',
                       help='Include non-representative proteins in analysis')
    parser.add_argument('--export-repair-list', help='Export batch IDs needing repair to file')
    parser.add_argument('--auto-repair', action='store_true',
                       help='Automatically trigger repair for critical batches')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.verbose)

    # Parse config
    try:
        config = parse_config(args.config)
    except Exception as e:
        logger.error(f"Error parsing config file: {str(e)}")
        sys.exit(1)

    # Convert severity threshold
    try:
        min_severity = Severity(args.severity_threshold)
    except ValueError:
        logger.error(f"Invalid severity threshold: {args.severity_threshold}")
        sys.exit(1)

    # Get database connection
    try:
        conn = get_db_connection(config)
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        sys.exit(1)

    try:
        # Get batch audit data
        logger.info("Analyzing batch completion status...")
        batch_data = get_batch_audit_data(conn, args.include_non_rep)
        logger.info(f"Retrieved data for {len(batch_data)} batches")

        # Analyze each batch for issues
        all_issues = []
        for batch in batch_data:
            batch_issues = analyze_batch_issues(batch)
            all_issues.extend(batch_issues)

        # Filter by severity
        filtered_issues = filter_issues_by_severity(all_issues, min_severity)
        
        logger.info(f"Found {len(filtered_issues)} issues (>= {min_severity.value}) across {len(set(issue.batch_id for issue in filtered_issues))} batches")

        # Output results
        if args.output_format == 'table':
            output_table_format(filtered_issues, logger)
        elif args.output_format == 'json':
            output_json_format(filtered_issues, logger)
        elif args.output_format == 'csv':
            output_csv_format(filtered_issues, logger)

        # Export repair list if requested
        if args.export_repair_list:
            export_repair_list(filtered_issues, args.export_repair_list, logger)

        # Auto-repair if requested
        if args.auto_repair:
            trigger_auto_repair(filtered_issues, args.config, logger)

    except Exception as e:
        logger.error(f"Error during analysis: {str(e)}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
