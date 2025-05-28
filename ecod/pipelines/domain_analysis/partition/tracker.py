#!/usr/bin/env python3
"""
Status tracker for domain partitioning.

This module handles all database status tracking and file registration
for the domain partition service.
"""

import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, List

from ecod.db import DBManager


class StatusTracker:
    """
    Tracks processing status in the database.

    This class encapsulates all database operations related to
    status tracking, file registration, and batch management.
    """

    def __init__(self, db_manager: DBManager):
        """
        Initialize the status tracker.

        Args:
            db_manager: Database manager instance
        """
        self.db = db_manager
        self.logger = logging.getLogger(__name__)

    def update_process_status(self, process_id: int, stage: str, status: str,
                            error_message: Optional[str] = None) -> bool:
        """
        Update process status in database.

        Args:
            process_id: Process ID to update
            stage: Current stage name
            status: Status ("processing", "success", "error", etc.)
            error_message: Optional error message

        Returns:
            bool: True if update succeeded, False if failed
        """
        try:
            update_query = """
            UPDATE ecod_schema.process_status
            SET current_stage = %s,
                status = %s,
                last_updated = CURRENT_TIMESTAMP
            """

            params = [stage, status]

            if error_message:
                update_query += ", error_message = %s"
                params.append(error_message)

            update_query += " WHERE id = %s"
            params.append(process_id)

            # Execute the update with error handling
            self.db.execute_query(update_query, tuple(params))

            # Log successful update
            self.logger.debug(f"Updated process {process_id} to {stage}/{status}")
            return True

        except Exception as e:
            # Log the error but don't raise - return False to indicate failure
            self.logger.error(f"Failed to update process status for {process_id}: {e}")
            return False

    def register_domain_file(self, process_id: int, file_path: str,
                            output_dir: str) -> bool:
        """
        Register domain partition file in database.

        Args:
            process_id: Process ID
            file_path: Path to domain file
            output_dir: Output directory

        Returns:
            bool: True if registration succeeded, False if failed
        """
        try:
            # Register the file
            query = """
            INSERT INTO ecod_schema.process_file
            (process_id, file_type, file_path, file_exists, created_at)
            VALUES (%s, %s, %s, %s, CURRENT_TIMESTAMP)
            ON CONFLICT (process_id, file_type)
            DO UPDATE SET
                file_path = EXCLUDED.file_path,
                file_exists = EXCLUDED.file_exists,
                created_at = EXCLUDED.created_at
            """

            import os
            file_exists = os.path.exists(file_path)

            self.db.execute_query(query, (process_id, 'domain_partition', file_path, file_exists))

            self.logger.debug(f"Registered domain file for process {process_id}: {file_path}")
            return True

        except Exception as e:
            self.logger.error(f"Failed to register domain file for process {process_id}: {e}")
            return False

    def update_batch_completion_status(self, batch_id: int,
                                     representatives_only: bool = False) -> bool:
        """
        Update batch completion status.

        Args:
            batch_id: Batch ID
            representatives_only: Whether processing was limited to representatives

        Returns:
            bool: True if update succeeded, False if failed
        """
        try:
            # Your existing batch completion logic here
            # Just wrap it in try/catch and return bool

            # Example logic (replace with actual implementation):
            query = """
            UPDATE ecod_schema.batch_status
            SET status = 'completed',
                completed_at = CURRENT_TIMESTAMP,
                representatives_only = %s
            WHERE batch_id = %s
            """

            self.db.execute_query(query, (representatives_only, batch_id))

            self.logger.info(f"Updated batch {batch_id} completion status")
            return True

        except Exception as e:
            self.logger.error(f"Failed to update batch completion status for {batch_id}: {e}")
            return False

    def update_non_representative_status(self, batch_id: int) -> bool:
        """
        Update non-representative protein status.

        Args:
            batch_id: Batch ID

        Returns:
            bool: True if update succeeded, False if failed
        """
        try:
            # Your existing non-representative update logic here
            # Wrap in try/catch and return bool

            query = """
            UPDATE ecod_schema.process_status
            SET status = 'skipped',
                current_stage = 'non_representative_skipped',
                last_updated = CURRENT_TIMESTAMP
            WHERE batch_id = %s
              AND is_representative = FALSE
              AND status = 'pending'
            """

            self.db.execute_query(query, (batch_id,))

            self.logger.info(f"Updated non-representative status for batch {batch_id}")
            return True

        except Exception as e:
            self.logger.error(f"Failed to update non-representative status for batch {batch_id}: {e}")
            return False

    def get_batch_progress(self, batch_id: int) -> Dict[str, Any]:
        """
        Get detailed batch progress information.

        Args:
            batch_id: Batch ID

        Returns:
            Dictionary with progress information
        """
        try:
            query = """
            SELECT
                COUNT(*) as total,
                COUNT(CASE WHEN ps.current_stage = 'domain_partition_complete'
                           AND ps.status = 'success' THEN 1 END) as complete,
                COUNT(CASE WHEN ps.status = 'error' THEN 1 END) as errors,
                COUNT(CASE WHEN ps.status = 'processing' THEN 1 END) as processing,
                COUNT(CASE WHEN ps.status = 'skipped' THEN 1 END) as skipped,
                COUNT(CASE WHEN ps.is_representative = true THEN 1 END) as representatives,
                COUNT(CASE WHEN pf.id IS NOT NULL THEN 1 END) as files_created
            FROM ecod_schema.process_status ps
            LEFT JOIN ecod_schema.process_file pf ON (
                pf.process_id = ps.id AND
                pf.file_type = 'domain_partition' AND
                pf.file_exists = TRUE
            )
            WHERE ps.batch_id = %s
            GROUP BY ps.batch_id
            """

            result = self.db.execute_dict_query(query, (batch_id,))

            if result:
                progress = result[0]

                # Calculate percentages
                total = progress['total']
                if total > 0:
                    progress['complete_pct'] = (progress['complete'] / total) * 100
                    progress['error_pct'] = (progress['errors'] / total) * 100
                    progress['progress_pct'] = ((progress['complete'] + progress['errors']) / total) * 100
                else:
                    progress['complete_pct'] = 0
                    progress['error_pct'] = 0
                    progress['progress_pct'] = 0

                return progress
            else:
                return {
                    'total': 0,
                    'complete': 0,
                    'errors': 0,
                    'processing': 0,
                    'skipped': 0,
                    'representatives': 0,
                    'files_created': 0,
                    'complete_pct': 0,
                    'error_pct': 0,
                    'progress_pct': 0
                }

        except Exception as e:
            self.logger.error(f"Error getting batch progress: {e}")
            return {}

    def get_process_history(self, process_id: int) -> List[Dict[str, Any]]:
        """
        Get processing history for a specific process.

        Args:
            process_id: Process ID

        Returns:
            List of status history records
        """
        try:
            # This would require a process_history table to track all changes
            # For now, return current status
            query = """
            SELECT ps.*, p.pdb_id, p.chain_id
            FROM ecod_schema.process_status ps
            JOIN ecod_schema.protein p ON ps.protein_id = p.id
            WHERE ps.id = %s
            """

            result = self.db.execute_dict_query(query, (process_id,))

            if result:
                return result
            else:
                return []

        except Exception as e:
            self.logger.error(f"Error getting process history: {e}")
            return []

    def cleanup_orphaned_files(self, batch_id: int, base_path: str) -> int:
        """
        Clean up orphaned domain files that don't have corresponding process records.

        Args:
            batch_id: Batch ID
            base_path: Base path for files

        Returns:
            Number of files cleaned up
        """
        try:
            # Get all registered files
            query = """
            SELECT pf.file_path
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            WHERE ps.batch_id = %s AND pf.file_type = 'domain_partition'
            """

            registered_files = set()
            for row in self.db.execute_query(query, (batch_id,)):
                registered_files.add(row[0])

            # Check actual files in domains directory
            domains_dir = Path(base_path) / "domains"
            if not domains_dir.exists():
                return 0

            cleaned = 0
            for file_path in domains_dir.glob("*.domains.xml"):
                relative_path = os.path.relpath(file_path, base_path)

                if relative_path not in registered_files:
                    # Orphaned file
                    try:
                        file_path.unlink()
                        cleaned += 1
                        self.logger.info(f"Cleaned up orphaned file: {relative_path}")
                    except Exception as e:
                        self.logger.error(f"Error deleting orphaned file {file_path}: {e}")

            return cleaned

        except Exception as e:
            self.logger.error(f"Error cleaning up orphaned files: {e}")
            return 0

    def verify_file_consistency(self, batch_id: int, base_path: str) -> Dict[str, int]:
        """
        Verify consistency between database records and actual files.

        Args:
            batch_id: Batch ID
            base_path: Base path for files

        Returns:
            Dictionary with verification results
        """
        try:
            query = """
            SELECT pf.id, pf.file_path, pf.file_exists, pf.file_size
            FROM ecod_schema.process_file pf
            JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
            WHERE ps.batch_id = %s AND pf.file_type = 'domain_partition'
            """

            results = {
                'total': 0,
                'consistent': 0,
                'missing': 0,
                'size_mismatch': 0,
                'status_mismatch': 0
            }

            for row in self.db.execute_dict_query(query, (batch_id,)):
                results['total'] += 1

                file_id = row['id']
                relative_path = row['file_path']
                db_exists = row['file_exists']
                db_size = row['file_size']

                # Check actual file
                full_path = os.path.join(base_path, relative_path)
                actual_exists = os.path.exists(full_path)
                actual_size = os.path.getsize(full_path) if actual_exists else 0

                # Compare
                if actual_exists == db_exists and actual_size == db_size:
                    results['consistent'] += 1
                else:
                    if not actual_exists:
                        results['missing'] += 1
                    if actual_exists != db_exists:
                        results['status_mismatch'] += 1
                    if actual_exists and actual_size != db_size:
                        results['size_mismatch'] += 1

                    # Update database with correct info
                    self.db.update(
                        "ecod_schema.process_file",
                        {
                            "file_exists": actual_exists,
                            "file_size": actual_size,
                            "last_checked": datetime.now()
                        },
                        "id = %s",
                        (file_id,)
                    )

            self.logger.info(
                f"File consistency check for batch {batch_id}: "
                f"{results['consistent']}/{results['total']} consistent, "
                f"{results['missing']} missing, "
                f"{results['size_mismatch']} size mismatches"
            )

            return results

        except Exception as e:
            self.logger.error(f"Error verifying file consistency: {e}")
            return {}
