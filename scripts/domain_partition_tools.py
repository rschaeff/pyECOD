#!/usr/bin/env python3
"""
Domain Partition Tools - Diagnostic utilities using service-based architecture

This script provides diagnostic tools for domain partitioning issues,
leveraging the new service components for better integration.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.core.context import ApplicationContext
from ecod.db import DBManager

# Service-based imports
from ecod.pipelines.domain_analysis.partition import (
    DomainPartitionService,
    StatusTracker,
    EvidenceAnalyzer,
    PartitionOptions
)
from ecod.models.pipeline.partition import DomainPartitionResult

# Utility imports
from ecod.utils.path_utils import (
    get_standardized_paths,
    get_all_evidence_paths,
    find_files_with_legacy_paths
)


class DomainPartitionDiagnostics:
    """Diagnostic tools for domain partition issues using service architecture"""

    def __init__(self, context: ApplicationContext):
        """Initialize diagnostics with application context"""
        self.context = context
        self.logger = logging.getLogger("domain_partition.diagnostics")

        # Initialize service components
        self.service = DomainPartitionService(context)
        self.tracker = self.service.tracker
        self.analyzer = self.service.analyzer

    def diagnose_process(self, process_id: int, fix_in_db: bool = False) -> bool:
        """
        Diagnose domain partition issues for a specific process

        Args:
            process_id: Process ID to diagnose
            fix_in_db: Whether to fix file paths in database

        Returns:
            True if diagnosis successful
        """
        self.logger.info(f"Diagnosing process {process_id}")

        # Get process details
        process_info = self._get_process_info(process_id)
        if not process_info:
            self.logger.error(f"Process {process_id} not found")
            return False

        pdb_id = process_info['pdb_id']
        chain_id = process_info['chain_id']
        batch_path = process_info['base_path']
        ref_version = process_info['ref_version']

        self.logger.info(f"Protein: {pdb_id}_{chain_id}")
        self.logger.info(f"Batch: {process_info['batch_id']} ({batch_path})")
        self.logger.info(f"Reference: {ref_version}")
        self.logger.info(f"Status: {process_info['current_stage']} / {process_info['status']}")

        # Run diagnostic checks
        diagnostics = {
            'process_info': process_info,
            'database_files': self._check_database_files(process_id, batch_path),
            'filesystem_files': self._check_filesystem_files(batch_path, pdb_id, chain_id, ref_version),
            'evidence_analysis': self._analyze_evidence_availability(batch_path, pdb_id, chain_id, ref_version),
            'partition_readiness': self._check_partition_readiness(process_id, batch_path, pdb_id, chain_id, ref_version),
            'recommendations': []
        }

        # Generate recommendations
        diagnostics['recommendations'] = self._generate_recommendations(diagnostics)

        # Display diagnostics
        self._display_diagnostics(diagnostics)

        # Fix database if requested
        if fix_in_db:
            self._fix_database_paths(process_id, diagnostics)

        return True

    def diagnose_batch(self, batch_id: int, sample_size: int = 10) -> bool:
        """
        Diagnose domain partition issues for an entire batch

        Args:
            batch_id: Batch ID to diagnose
            sample_size: Number of proteins to sample

        Returns:
            True if diagnosis successful
        """
        self.logger.info(f"Diagnosing batch {batch_id}")

        # Get batch statistics using service's tracker
        progress = self.tracker.get_batch_progress(batch_id)

        if not progress:
            self.logger.error(f"Batch {batch_id} not found")
            return False

        # Display batch overview
        self.logger.info("\nBatch Overview:")
        self.logger.info(f"  Total proteins: {progress['total']}")
        self.logger.info(f"  Completed: {progress['complete']} ({progress['complete_pct']:.1f}%)")
        self.logger.info(f"  Errors: {progress['errors']} ({progress['error_pct']:.1f}%)")
        self.logger.info(f"  Processing: {progress['processing']}")
        self.logger.info(f"  Skipped: {progress['skipped']}")
        self.logger.info(f"  Files created: {progress['files_created']}")

        # Sample problematic proteins
        problematic_proteins = self._get_problematic_proteins(batch_id, sample_size)

        if problematic_proteins:
            self.logger.info(f"\nSampling {len(problematic_proteins)} problematic proteins:")

            for i, protein in enumerate(problematic_proteins[:5]):  # Show first 5
                self.logger.info(f"\n  {i+1}. {protein['pdb_id']}_{protein['chain_id']} "
                               f"(Process {protein['process_id']})")
                self.logger.info(f"     Stage: {protein['current_stage']}")
                self.logger.info(f"     Status: {protein['status']}")
                if protein['error_message']:
                    self.logger.info(f"     Error: {protein['error_message'][:100]}...")

        # Check file consistency
        self.logger.info("\nChecking file consistency...")
        consistency = self.tracker.verify_file_consistency(batch_id, self._get_batch_path(batch_id))

        if consistency:
            self.logger.info(f"  Total files: {consistency['total']}")
            self.logger.info(f"  Consistent: {consistency['consistent']}")
            self.logger.info(f"  Missing: {consistency['missing']}")
            self.logger.info(f"  Size mismatches: {consistency['size_mismatch']}")

        return True

    def analyze_failures(self, batch_id: Optional[int] = None, limit: int = 20) -> bool:
        """
        Analyze common failure patterns

        Args:
            batch_id: Optional batch ID to limit analysis
            limit: Maximum number of failures to analyze

        Returns:
            True if analysis successful
        """
        self.logger.info("Analyzing domain partition failures")

        # Get failed processes
        query = """
        SELECT ps.id as process_id, p.pdb_id, p.chain_id, ps.batch_id,
               ps.current_stage, ps.error_message, b.base_path
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        JOIN ecod_schema.batch b ON ps.batch_id = b.id
        WHERE ps.status = 'error'
          AND ps.current_stage LIKE '%domain_partition%'
        """

        params = []
        if batch_id:
            query += " AND ps.batch_id = %s"
            params.append(batch_id)

        query += f" ORDER BY ps.updated_at DESC LIMIT {limit}"

        failures = self.context.db.execute_dict_query(query, tuple(params))

        if not failures:
            self.logger.info("No domain partition failures found")
            return True

        # Analyze failure patterns
        failure_patterns = {}

        for failure in failures:
            error_msg = failure['error_message'] or 'Unknown error'

            # Categorize errors
            if 'summary not found' in error_msg.lower():
                category = 'Missing domain summary'
            elif 'no evidence' in error_msg.lower():
                category = 'No evidence found'
            elif 'validation' in error_msg.lower():
                category = 'Validation error'
            elif 'xml' in error_msg.lower() or 'parse' in error_msg.lower():
                category = 'XML parsing error'
            elif 'permission' in error_msg.lower() or 'access' in error_msg.lower():
                category = 'File access error'
            else:
                category = 'Other'

            if category not in failure_patterns:
                failure_patterns[category] = []

            failure_patterns[category].append(failure)

        # Display analysis
        self.logger.info(f"\nAnalyzed {len(failures)} failures:")

        for category, cases in sorted(failure_patterns.items(), key=lambda x: len(x[1]), reverse=True):
            self.logger.info(f"\n{category}: {len(cases)} cases")

            # Show examples
            for i, case in enumerate(cases[:3]):  # Show first 3 examples
                self.logger.info(f"  {i+1}. {case['pdb_id']}_{case['chain_id']} "
                               f"(Batch {case['batch_id']})")
                if case['error_message']:
                    error_preview = case['error_message'][:150]
                    if len(case['error_message']) > 150:
                        error_preview += "..."
                    self.logger.info(f"     {error_preview}")

        # Provide recommendations
        self.logger.info("\nRecommendations:")

        if 'Missing domain summary' in failure_patterns:
            count = len(failure_patterns['Missing domain summary'])
            self.logger.info(f"  - {count} proteins need domain summary generation")
            self.logger.info("    Run: python scripts/domain_summary_run.py process missing ...")

        if 'No evidence found' in failure_patterns:
            count = len(failure_patterns['No evidence found'])
            self.logger.info(f"  - {count} proteins have no valid evidence")
            self.logger.info("    These may be truly unclassifiable or need BLAST re-run")

        if 'XML parsing error' in failure_patterns:
            count = len(failure_patterns['XML parsing error'])
            self.logger.info(f"  - {count} proteins have corrupted XML files")
            self.logger.info("    Run: python scripts/domain_partition_run.py repair missing ...")

        return True

    def check_readiness(self, batch_id: int, mode: str = 'summary') -> bool:
        """
        Check batch readiness for domain partition

        Args:
            batch_id: Batch ID to check
            mode: 'summary' to show summary, 'detailed' for protein list

        Returns:
            True if batch is ready
        """
        self.logger.info(f"Checking domain partition readiness for batch {batch_id}")

        # Check various readiness criteria
        query = """
        SELECT
            COUNT(*) as total,
            SUM(CASE WHEN ps.is_representative = true THEN 1 ELSE 0 END) as representatives,
            SUM(CASE WHEN EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                  AND pf.file_type = 'domain_summary'
                  AND pf.file_exists = true
            ) THEN 1 ELSE 0 END) as has_summary,
            SUM(CASE WHEN EXISTS (
                SELECT 1 FROM ecod_schema.process_file pf
                WHERE pf.process_id = ps.id
                  AND pf.file_type = 'domain_partition'
                  AND pf.file_exists = true
            ) THEN 1 ELSE 0 END) as has_partition,
            SUM(CASE WHEN ps.status = 'error' THEN 1 ELSE 0 END) as errors
        FROM ecod_schema.process_status ps
        WHERE ps.batch_id = %s
        """

        result = self.context.db.execute_dict_query(query, (batch_id,))[0]

        # Calculate readiness metrics
        total = result['total']
        ready_for_partition = result['has_summary'] - result['has_partition']
        readiness_pct = (result['has_summary'] / total * 100) if total > 0 else 0
        completion_pct = (result['has_partition'] / total * 100) if total > 0 else 0

        # Display summary
        self.logger.info("\nReadiness Summary:")
        self.logger.info(f"  Total proteins: {total}")
        self.logger.info(f"  Representatives: {result['representatives']}")
        self.logger.info(f"  Have domain summary: {result['has_summary']} ({readiness_pct:.1f}%)")
        self.logger.info(f"  Have domain partition: {result['has_partition']} ({completion_pct:.1f}%)")
        self.logger.info(f"  Ready for partition: {ready_for_partition}")
        self.logger.info(f"  Errors: {result['errors']}")

        # Check if batch is ready
        is_ready = readiness_pct >= 90.0  # 90% threshold

        if is_ready:
            self.logger.info("\n✓ Batch is ready for domain partition processing")
        else:
            self.logger.warning(f"\n✗ Batch is NOT ready (only {readiness_pct:.1f}% have summaries)")

        # Show detailed list if requested
        if mode == 'detailed' and ready_for_partition > 0:
            query = """
            SELECT p.pdb_id, p.chain_id, ps.id as process_id
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE ps.batch_id = %s
              AND EXISTS (
                  SELECT 1 FROM ecod_schema.process_file pf
                  WHERE pf.process_id = ps.id
                    AND pf.file_type = 'domain_summary'
                    AND pf.file_exists = true
              )
              AND NOT EXISTS (
                  SELECT 1 FROM ecod_schema.process_file pf
                  WHERE pf.process_id = ps.id
                    AND pf.file_type = 'domain_partition'
                    AND pf.file_exists = true
              )
            LIMIT 20
            """

            proteins = self.context.db.execute_dict_query(query, (batch_id,))

            if proteins:
                self.logger.info(f"\nSample proteins ready for partition (showing {len(proteins)}):")
                for i, protein in enumerate(proteins):
                    self.logger.info(f"  {i+1}. {protein['pdb_id']}_{protein['chain_id']} "
                                   f"(Process {protein['process_id']})")

        return is_ready

    def _get_process_info(self, process_id: int) -> Optional[Dict[str, Any]]:
        """Get detailed process information"""
        query = """
        SELECT
            ps.id as process_id,
            p.pdb_id,
            p.chain_id,
            ps.batch_id,
            ps.relative_path,
            ps.current_stage,
            ps.status,
            ps.error_message,
            ps.is_representative,
            b.base_path,
            b.ref_version,
            b.batch_name
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE
            ps.id = %s
        """

        results = self.context.db.execute_dict_query(query, (process_id,))
        return results[0] if results else None

    def _check_database_files(self, process_id: int, batch_path: str) -> Dict[str, Any]:
        """Check file records in database"""
        query = """
        SELECT pf.id, pf.file_type, pf.file_path, pf.file_exists, pf.file_size
        FROM ecod_schema.process_file pf
        WHERE pf.process_id = %s
        ORDER BY pf.file_type
        """

        file_records = self.context.db.execute_dict_query(query, (process_id,))

        results = {}
        for record in file_records:
            file_path = record['file_path']
            full_path = os.path.join(batch_path, file_path) if not os.path.isabs(file_path) else file_path
            fs_exists = os.path.exists(full_path)
            fs_size = os.path.getsize(full_path) if fs_exists else 0

            results[record['file_type']] = {
                'db_path': file_path,
                'full_path': full_path,
                'db_exists': record['file_exists'],
                'db_size': record['file_size'],
                'fs_exists': fs_exists,
                'fs_size': fs_size,
                'consistent': record['file_exists'] == fs_exists and record['file_size'] == fs_size
            }

        return results

    def _check_filesystem_files(self, batch_path: str, pdb_id: str, chain_id: str,
                               ref_version: str) -> Dict[str, Any]:
        """Check files on filesystem using path utilities"""
        results = {}

        # Check standard paths
        standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)
        for file_type, path in standard_paths.items():
            results[f"{file_type}_standard"] = {
                'path': path,
                'exists': os.path.exists(path),
                'size': os.path.getsize(path) if os.path.exists(path) else 0
            }

        # Check legacy paths
        legacy_files = find_files_with_legacy_paths(batch_path, pdb_id, chain_id, ref_version)
        for file_type, info in legacy_files.items():
            if info['exists_at']:
                results[f"{file_type}_legacy"] = {
                    'path': info['exists_at'],
                    'exists': True,
                    'size': os.path.getsize(info['exists_at'])
                }

        # Check comprehensive paths
        evidence_paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, ref_version)
        for file_type, info in evidence_paths.items():
            if info['exists_at'] and f"{file_type}_standard" not in results:
                results[f"{file_type}_found"] = {
                    'path': info['exists_at'],
                    'exists': True,
                    'size': os.path.getsize(info['exists_at'])
                }

        return results

    def _analyze_evidence_availability(self, batch_path: str, pdb_id: str,
                                      chain_id: str, ref_version: str) -> Dict[str, Any]:
        """Analyze what evidence is available for partition"""

        # Find domain summary
        summary_path = None
        for summary_type in ['domain_summary', 'blast_only_summary']:
            paths = get_all_evidence_paths(batch_path, pdb_id, chain_id, ref_version)
            if summary_type in paths and paths[summary_type]['exists_at']:
                summary_path = paths[summary_type]['exists_at']
                break

        if not summary_path:
            return {
                'summary_found': False,
                'summary_path': None,
                'evidence_count': 0,
                'evidence_types': []
            }

        # Parse summary to check evidence
        try:
            summary_data = self.analyzer.parse_domain_summary(summary_path)

            if 'error' in summary_data:
                return {
                    'summary_found': True,
                    'summary_path': summary_path,
                    'parse_error': summary_data['error'],
                    'evidence_count': 0,
                    'evidence_types': []
                }

            # Count evidence
            evidence_types = []
            evidence_count = 0

            if summary_data.get('chain_blast_evidence'):
                count = len(summary_data['chain_blast_evidence'])
                evidence_count += count
                evidence_types.append(f"chain_blast ({count})")

            if summary_data.get('domain_blast_evidence'):
                count = len(summary_data['domain_blast_evidence'])
                evidence_count += count
                evidence_types.append(f"domain_blast ({count})")

            if summary_data.get('hhsearch_evidence'):
                count = len(summary_data['hhsearch_evidence'])
                evidence_count += count
                evidence_types.append(f"hhsearch ({count})")

            return {
                'summary_found': True,
                'summary_path': summary_path,
                'evidence_count': evidence_count,
                'evidence_types': evidence_types,
                'is_peptide': summary_data.get('is_peptide', False),
                'sequence_length': summary_data.get('sequence_length', 0),
                'validation_stats': summary_data.get('validation_stats', {})
            }

        except Exception as e:
            return {
                'summary_found': True,
                'summary_path': summary_path,
                'parse_error': str(e),
                'evidence_count': 0,
                'evidence_types': []
            }

    def _check_partition_readiness(self, process_id: int, batch_path: str,
                                  pdb_id: str, chain_id: str, ref_version: str) -> Dict[str, Any]:
        """Check if process is ready for partition"""
        evidence = self._analyze_evidence_availability(batch_path, pdb_id, chain_id, ref_version)

        ready = (
            evidence['summary_found'] and
            'parse_error' not in evidence and
            (evidence['evidence_count'] > 0 or evidence.get('is_peptide', False))
        )

        reasons = []
        if not evidence['summary_found']:
            reasons.append("No domain summary file found")
        elif 'parse_error' in evidence:
            reasons.append(f"Summary parsing error: {evidence['parse_error']}")
        elif evidence['evidence_count'] == 0 and not evidence.get('is_peptide', False):
            reasons.append("No evidence found in summary (and not a peptide)")

        return {
            'ready': ready,
            'reasons': reasons,
            'evidence': evidence
        }

    def _generate_recommendations(self, diagnostics: Dict[str, Any]) -> List[str]:
        """Generate recommendations based on diagnostics"""
        recommendations = []

        readiness = diagnostics['partition_readiness']

        if readiness['ready']:
            recommendations.append("✓ Ready for domain partition - run:")
            recommendations.append(f"  python scripts/domain_partition_run.py process single "
                                 f"--pdb-id {diagnostics['process_info']['pdb_id']} "
                                 f"--chain-id {diagnostics['process_info']['chain_id']} "
                                 f"--batch-id {diagnostics['process_info']['batch_id']}")
        else:
            recommendations.append("✗ Not ready for domain partition:")
            for reason in readiness['reasons']:
                recommendations.append(f"  - {reason}")

            # Specific recommendations
            if not readiness['evidence']['summary_found']:
                # Check if BLAST files exist
                fs_files = diagnostics['filesystem_files']
                has_blast = any(
                    key.startswith(('chain_blast', 'domain_blast')) and info['exists']
                    for key, info in fs_files.items()
                )

                if has_blast:
                    recommendations.append("\n  Action: Generate domain summary")
                    recommendations.append("  python scripts/domain_summary_run.py process single "
                                         f"--pdb-id {diagnostics['process_info']['pdb_id']} "
                                         f"--chain-id {diagnostics['process_info']['chain_id']}")
                else:
                    recommendations.append("\n  Action: Run BLAST analysis first")
                    recommendations.append("  python scripts/run_blast.py ...")

        return recommendations

    def _display_diagnostics(self, diagnostics: Dict[str, Any]) -> None:
        """Display diagnostic results"""
        info = diagnostics['process_info']

        # Database files
        self.logger.info("\nDatabase File Records:")
        db_files = diagnostics['database_files']
        if not db_files:
            self.logger.warning("  No file records found in database")
        else:
            for file_type, file_info in sorted(db_files.items()):
                status = "✓" if file_info['consistent'] else "✗"
                self.logger.info(f"  [{status}] {file_type}:")
                self.logger.info(f"      Path: {file_info['db_path']}")
                self.logger.info(f"      DB: exists={file_info['db_exists']}, "
                               f"size={file_info['db_size']}")
                self.logger.info(f"      FS: exists={file_info['fs_exists']}, "
                               f"size={file_info['fs_size']}")

        # Filesystem files
        self.logger.info("\nFilesystem Files:")
        fs_files = diagnostics['filesystem_files']
        for file_type, file_info in sorted(fs_files.items()):
            if file_info['exists']:
                self.logger.info(f"  [✓] {file_type}: {file_info['path']}")
            else:
                self.logger.debug(f"  [✗] {file_type}: {file_info['path']}")

        # Evidence analysis
        self.logger.info("\nEvidence Analysis:")
        evidence = diagnostics['partition_readiness']['evidence']
        if evidence['summary_found']:
            self.logger.info(f"  Summary found: {evidence['summary_path']}")
            if 'parse_error' in evidence:
                self.logger.error(f"  Parse error: {evidence['parse_error']}")
            else:
                self.logger.info(f"  Evidence count: {evidence['evidence_count']}")
                self.logger.info(f"  Evidence types: {', '.join(evidence['evidence_types'])}")
                self.logger.info(f"  Is peptide: {evidence.get('is_peptide', False)}")
                self.logger.info(f"  Sequence length: {evidence.get('sequence_length', 0)}")
        else:
            self.logger.warning("  No domain summary found")

        # Recommendations
        self.logger.info("\nRecommendations:")
        for rec in diagnostics['recommendations']:
            self.logger.info(rec)

    def _fix_database_paths(self, process_id: int, diagnostics: Dict[str, Any]) -> None:
        """Fix database paths if legacy files found"""
        self.logger.info("\nFixing database paths...")

        db_files = diagnostics['database_files']
        fs_files = diagnostics['filesystem_files']

        fixed_count = 0

        for file_type in ['domain_summary', 'blast_only_summary']:
            if file_type in db_files and not db_files[file_type]['fs_exists']:
                # Check if we found it elsewhere
                legacy_key = f"{file_type}_legacy"
                found_key = f"{file_type}_found"

                new_path = None
                if legacy_key in fs_files and fs_files[legacy_key]['exists']:
                    new_path = fs_files[legacy_key]['path']
                elif found_key in fs_files and fs_files[found_key]['exists']:
                    new_path = fs_files[found_key]['path']

                if new_path:
                    # Update database
                    batch_path = diagnostics['process_info']['base_path']
                    rel_path = os.path.relpath(new_path, batch_path)

                    self.tracker.register_domain_file(process_id, new_path, batch_path)
                    self.logger.info(f"  Fixed {file_type}: {rel_path}")
                    fixed_count += 1

        if fixed_count > 0:
            self.logger.info(f"Fixed {fixed_count} file paths in database")
        else:
            self.logger.info("No paths needed fixing")

    def _get_problematic_proteins(self, batch_id: int, limit: int) -> List[Dict[str, Any]]:
        """Get proteins with problems in the batch"""
        query = """
        SELECT ps.id as process_id, p.pdb_id, p.chain_id,
               ps.current_stage, ps.status, ps.error_message
        FROM ecod_schema.process_status ps
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE ps.batch_id = %s
          AND (ps.status = 'error' OR
               (ps.current_stage LIKE '%summary%' AND ps.status = 'success' AND NOT EXISTS (
                   SELECT 1 FROM ecod_schema.process_file pf
                   WHERE pf.process_id = ps.id
                     AND pf.file_type = 'domain_partition'
                     AND pf.file_exists = true
               )))
        ORDER BY ps.status DESC, ps.updated_at DESC
        LIMIT %s
        """

        return self.context.db.execute_dict_query(query, (batch_id, limit))

    def _get_batch_path(self, batch_id: int) -> str:
        """Get batch base path"""
        query = "SELECT base_path FROM ecod_schema.batch WHERE id = %s"
        result = self.context.db.execute_query(query, (batch_id,))
        return result[0][0] if result else ""


def setup_logging(verbose: bool = False, log_file: Optional[str] = None) -> None:
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO

    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Domain Partition Diagnostic Tools - Service-Based Architecture'
    )
    parser.add_argument('--config', type=str, default='config/config.yml',
                      help='Path to configuration file')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    # Create subparsers for different diagnostic modes
    subparsers = parser.add_subparsers(dest="mode", help="Diagnostic mode")

    # Diagnose process
    process_parser = subparsers.add_parser("process", help="Diagnose a specific process")
    process_parser.add_argument('--process-id', type=int, required=True,
                              help='Process ID to diagnose')
    process_parser.add_argument('--fix-in-db', action='store_true',
                              help='Fix file paths in database if legacy paths are found')

    # Diagnose batch
    batch_parser = subparsers.add_parser("batch", help="Diagnose an entire batch")
    batch_parser.add_argument('--batch-id', type=int, required=True,
                            help='Batch ID to diagnose')
    batch_parser.add_argument('--sample-size', type=int, default=10,
                            help='Number of problematic proteins to sample')

    # Analyze failures
    failures_parser = subparsers.add_parser("failures", help="Analyze failure patterns")
    failures_parser.add_argument('--batch-id', type=int,
                               help='Limit to specific batch')
    failures_parser.add_argument('--limit', type=int, default=20,
                               help='Maximum failures to analyze')

    # Check readiness
    readiness_parser = subparsers.add_parser("readiness", help="Check batch readiness")
    readiness_parser.add_argument('--batch-id', type=int, required=True,
                                help='Batch ID to check')
    readiness_parser.add_argument('--detailed', action='store_true',
                                help='Show detailed protein list')

    args = parser.parse_args()

    if not args.mode:
        parser.print_help()
        return 1

    # Setup logging
    setup_logging(args.verbose, args.log_file)

    # Use absolute path for config if provided
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)

    # Initialize context and diagnostics
    context = ApplicationContext(args.config)
    diagnostics = DomainPartitionDiagnostics(context)

    # Run appropriate diagnostic
    try:
        if args.mode == "process":
            success = diagnostics.diagnose_process(args.process_id, args.fix_in_db)
        elif args.mode == "batch":
            success = diagnostics.diagnose_batch(args.batch_id, args.sample_size)
        elif args.mode == "failures":
            success = diagnostics.analyze_failures(args.batch_id, args.limit)
        elif args.mode == "readiness":
            mode = 'detailed' if args.detailed else 'summary'
            success = diagnostics.check_readiness(args.batch_id, mode)
        else:
            logging.getLogger().error(f"Unknown mode: {args.mode}")
            success = False

        return 0 if success else 1

    except Exception as e:
        logging.getLogger().error(f"Error running diagnostics: {str(e)}")
        import traceback
        logging.getLogger().error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
