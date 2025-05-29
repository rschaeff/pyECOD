# ecod/pipelines/domain_analysis/partition/tracker.py
"""
Comprehensive Status Tracker for Domain Analysis Pipeline

This tracker provides detailed monitoring of all pipeline stages with
comprehensive error handling, performance metrics, and recovery capabilities.
"""

import logging
import json
import time
from typing import Optional, Dict, Any, List, Union
from datetime import datetime, timedelta
from dataclasses import dataclass, field
from enum import Enum
from collections import defaultdict, deque


class ProcessStage(Enum):
    """Pipeline processing stages"""
    QUEUED = "queued"
    PARSING = "parsing"
    EVIDENCE_EXTRACTION = "evidence_extraction"
    DOMAIN_PARTITIONING = "domain_partitioning"
    CLASSIFICATION = "classification"
    VALIDATION = "validation"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


class ProcessStatus(Enum):
    """Processing status within a stage"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    RETRYING = "retrying"
    CANCELLED = "cancelled"


@dataclass
class ProcessInfo:
    """Information about a single process"""
    process_id: int
    pdb_id: str
    chain_id: str
    batch_id: Optional[int] = None
    stage: ProcessStage = ProcessStage.QUEUED
    status: ProcessStatus = ProcessStatus.PENDING
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    progress: float = 0.0  # 0-100
    error_message: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    files_created: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    retry_count: int = 0
    last_updated: datetime = field(default_factory=datetime.now)

    @property
    def processing_time(self) -> Optional[timedelta]:
        """Calculate processing time"""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        elif self.start_time:
            return datetime.now() - self.start_time
        return None

    @property
    def is_active(self) -> bool:
        """Check if process is currently active"""
        return self.status in (ProcessStatus.PENDING, ProcessStatus.RUNNING, ProcessStatus.RETRYING)

    @property
    def is_finished(self) -> bool:
        """Check if process is finished (completed, failed, or cancelled)"""
        return self.status in (ProcessStatus.COMPLETED, ProcessStatus.FAILED, ProcessStatus.CANCELLED)


@dataclass
class BatchInfo:
    """Information about a batch of processes"""
    batch_id: int
    name: str
    total_processes: int = 0
    completed_processes: int = 0
    failed_processes: int = 0
    skipped_processes: int = 0
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    estimated_completion: Optional[datetime] = None

    @property
    def success_rate(self) -> float:
        """Calculate success rate"""
        if self.total_processes == 0:
            return 0.0
        return (self.completed_processes / self.total_processes) * 100

    @property
    def is_completed(self) -> bool:
        """Check if batch is completed"""
        return (self.completed_processes + self.failed_processes + self.skipped_processes) >= self.total_processes


@dataclass
class PerformanceMetrics:
    """Performance tracking metrics"""
    processes_per_hour: float = 0.0
    average_processing_time: timedelta = timedelta()
    peak_memory_usage: float = 0.0  # MB
    database_operations_per_second: float = 0.0
    error_rate: float = 0.0  # percentage
    active_connections: int = 0
    queue_depth: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            'processes_per_hour': self.processes_per_hour,
            'average_processing_time_seconds': self.average_processing_time.total_seconds(),
            'peak_memory_usage_mb': self.peak_memory_usage,
            'database_operations_per_second': self.database_operations_per_second,
            'error_rate_percent': self.error_rate,
            'active_connections': self.active_connections,
            'queue_depth': self.queue_depth
        }


class StatusTracker:
    """
    Comprehensive status tracker for domain analysis pipeline.

    Provides detailed monitoring, error handling, performance metrics,
    and recovery capabilities for all pipeline operations.
    """

    def __init__(self, db_manager=None, config: Optional[Dict[str, Any]] = None):
        self.db_manager = db_manager
        self.config = config or {}
        self.logger = logging.getLogger(__name__)

        # Connection status
        self.db_available = self._test_database_connection()
        self.last_db_check = datetime.now()
        self.db_check_interval = timedelta(minutes=5)

        # In-memory tracking (for when database is unavailable)
        self.processes: Dict[int, ProcessInfo] = {}
        self.batches: Dict[int, BatchInfo] = {}

        # Performance tracking
        self.metrics = PerformanceMetrics()
        self.recent_operations = deque(maxlen=1000)  # Track recent operations for metrics
        self.error_history = deque(maxlen=100)  # Track recent errors

        # Configuration
        self.max_retry_attempts = self.config.get('max_retry_attempts', 3)
        self.retry_delay = self.config.get('retry_delay_seconds', 60)
        self.checkpoint_interval = self.config.get('checkpoint_interval_seconds', 300)
        self.cleanup_completed_after = timedelta(hours=self.config.get('cleanup_completed_hours', 24))

        # Performance monitoring
        self._last_metrics_update = datetime.now()
        self._operation_times = deque(maxlen=100)

        self.logger.info(f"StatusTracker initialized. Database available: {self.db_available}")

    def _test_database_connection(self) -> bool:
        """Test if database is available with comprehensive checks"""
        if not self.db_manager:
            return False

        try:
            # Try connection test if available
            if hasattr(self.db_manager, 'test_connection'):
                result = self.db_manager.test_connection()
                if result:
                    self.logger.debug("Database connection test passed")
                return result

            # Try a simple query
            if hasattr(self.db_manager, 'execute_query'):
                self.db_manager.execute_query("SELECT 1")
                self.logger.debug("Database query test passed")
                return True

            # Try basic database operation
            if hasattr(self.db_manager, 'execute'):
                cursor = self.db_manager.execute("SELECT 1")
                cursor.fetchone()
                self.logger.debug("Database execute test passed")
                return True

            return False

        except Exception as e:
            self.logger.warning(f"Database connection test failed: {str(e)}")
            return False

    def _maybe_reconnect_database(self) -> None:
        """Periodically try to reconnect to database if it's unavailable"""
        now = datetime.now()
        if (not self.db_available and
            now - self.last_db_check > self.db_check_interval):

            self.logger.info("Attempting database reconnection...")
            self.db_available = self._test_database_connection()
            self.last_db_check = now

            if self.db_available:
                self.logger.info("Database reconnection successful")
                self._sync_in_memory_to_database()
            else:
                self.logger.warning("Database reconnection failed")

    def _sync_in_memory_to_database(self) -> None:
        """Sync in-memory process data to database when connection is restored"""
        if not self.db_available:
            return

        try:
            synced_count = 0
            for process_info in self.processes.values():
                if self._write_process_to_database(process_info):
                    synced_count += 1

            self.logger.info(f"Synced {synced_count} processes to database")

        except Exception as e:
            self.logger.error(f"Error syncing to database: {str(e)}")
            self.db_available = False

    def start_process(self, process_id: int, pdb_id: str, chain_id: str,
                     batch_id: Optional[int] = None, **metadata) -> bool:
        """
        Start tracking a new process.

        Args:
            process_id: Unique process identifier
            pdb_id: PDB structure ID
            chain_id: Chain identifier
            batch_id: Optional batch identifier
            **metadata: Additional metadata to track

        Returns:
            bool: True if process started successfully
        """
        try:
            process_info = ProcessInfo(
                process_id=process_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                batch_id=batch_id,
                start_time=datetime.now(),
                stage=ProcessStage.QUEUED,
                status=ProcessStatus.PENDING,
                metadata=metadata
            )

            # Store in memory
            self.processes[process_id] = process_info

            # Try to store in database
            if self.db_available:
                self._write_process_to_database(process_info)

            self.logger.debug(f"Started tracking process {process_id} ({pdb_id}_{chain_id})")
            return True

        except Exception as e:
            self.logger.error(f"Error starting process {process_id}: {str(e)}")
            return False

    def update_process_stage(self, process_id: int, stage: ProcessStage,
                           status: ProcessStatus = ProcessStatus.RUNNING,
                           progress: Optional[float] = None,
                           metadata: Optional[Dict[str, Any]] = None) -> bool:
        """
        Update process stage and status.

        Args:
            process_id: Process identifier
            stage: New processing stage
            status: Processing status
            progress: Progress percentage (0-100)
            metadata: Additional metadata to update

        Returns:
            bool: True if update successful
        """
        try:
            self._maybe_reconnect_database()

            # Update in-memory record
            if process_id in self.processes:
                process_info = self.processes[process_id]
                old_stage = process_info.stage

                process_info.stage = stage
                process_info.status = status
                process_info.last_updated = datetime.now()

                if progress is not None:
                    process_info.progress = progress

                if metadata:
                    process_info.metadata.update(metadata)

                # Set end time if completed or failed
                if status in (ProcessStatus.COMPLETED, ProcessStatus.FAILED, ProcessStatus.CANCELLED):
                    process_info.end_time = datetime.now()

                self.logger.debug(f"Process {process_id} stage: {old_stage} -> {stage}, status: {status}")
            else:
                self.logger.warning(f"Process {process_id} not found for stage update")
                return False

            # Try to update database
            if self.db_available:
                return self._update_process_in_database(process_id, stage, status, progress, metadata)
            else:
                self.logger.debug(f"Database unavailable, cached stage update for process {process_id}")
                return True

        except Exception as e:
            self.logger.error(f"Error updating process stage {process_id}: {str(e)}")
            return False

    def update_process_status(self, process_id: int, stage: str, status: str) -> bool:
        """
        Legacy method for compatibility with existing code.

        Args:
            process_id: Process identifier
            stage: Stage name as string
            status: Status name as string

        Returns:
            bool: True if update successful
        """
        try:
            # Convert string values to enums
            stage_enum = ProcessStage(stage.lower())
            status_enum = ProcessStatus(status.lower())

            return self.update_process_stage(process_id, stage_enum, status_enum)

        except ValueError as e:
            self.logger.warning(f"Invalid stage/status values: {stage}/{status}: {str(e)}")
            return False
        except Exception as e:
            self.logger.error(f"Error in legacy status update: {str(e)}")
            return False

    def add_process_error(self, process_id: int, error_message: str,
                         retry: bool = True) -> bool:
        """
        Add an error to a process and optionally retry.

        Args:
            process_id: Process identifier
            error_message: Error description
            retry: Whether to attempt retry

        Returns:
            bool: True if error recorded successfully
        """
        try:
            if process_id not in self.processes:
                self.logger.warning(f"Process {process_id} not found for error recording")
                return False

            process_info = self.processes[process_id]
            process_info.error_message = error_message
            process_info.last_updated = datetime.now()

            # Track error for metrics
            self.error_history.append({
                'process_id': process_id,
                'error': error_message,
                'timestamp': datetime.now(),
                'stage': process_info.stage.value
            })

            # Determine if we should retry
            if retry and process_info.retry_count < self.max_retry_attempts:
                process_info.retry_count += 1
                process_info.status = ProcessStatus.RETRYING

                self.logger.info(f"Process {process_id} error (retry {process_info.retry_count}): {error_message}")

                # Schedule retry (in a real implementation, this would go to a queue)
                process_info.metadata['retry_scheduled'] = datetime.now() + timedelta(seconds=self.retry_delay)

            else:
                process_info.status = ProcessStatus.FAILED
                process_info.end_time = datetime.now()

                self.logger.error(f"Process {process_id} failed permanently: {error_message}")

            # Update metrics
            self._update_performance_metrics()

            return True

        except Exception as e:
            self.logger.error(f"Error recording process error {process_id}: {str(e)}")
            return False

    def add_process_warning(self, process_id: int, warning_message: str) -> bool:
        """
        Add a warning to a process.

        Args:
            process_id: Process identifier
            warning_message: Warning description

        Returns:
            bool: True if warning recorded successfully
        """
        try:
            if process_id not in self.processes:
                self.logger.warning(f"Process {process_id} not found for warning recording")
                return False

            process_info = self.processes[process_id]
            process_info.warnings.append(warning_message)
            process_info.last_updated = datetime.now()

            self.logger.warning(f"Process {process_id} warning: {warning_message}")
            return True

        except Exception as e:
            self.logger.error(f"Error recording process warning {process_id}: {str(e)}")
            return False

    def register_domain_file(self, process_id: int, file_path: str,
                           base_path: Optional[str] = None) -> bool:
        """
        Register a domain file created by a process.

        Args:
            process_id: Process identifier
            file_path: Full path to created file
            base_path: Base path for calculating relative path

        Returns:
            bool: True if file registered successfully
        """
        try:
            self._maybe_reconnect_database()

            # Calculate relative path
            relative_path = file_path
            if base_path and file_path.startswith(base_path):
                relative_path = file_path[len(base_path):].lstrip('/')

            # Update in-memory record
            if process_id in self.processes:
                process_info = self.processes[process_id]
                process_info.files_created.append(file_path)
                process_info.metadata['domain_file'] = relative_path
                process_info.last_updated = datetime.now()

                self.logger.debug(f"Registered domain file for process {process_id}: {relative_path}")

            # Try to register in database
            if self.db_available:
                return self._register_file_in_database(process_id, file_path, relative_path)
            else:
                self.logger.debug(f"Database unavailable, cached file registration for process {process_id}")
                return True

        except Exception as e:
            self.logger.error(f"Error registering domain file for {process_id}: {str(e)}")
            return False

    def get_process_status(self, process_id: int) -> Optional[Dict[str, Any]]:
        """
        Get detailed status information for a process.

        Args:
            process_id: Process identifier

        Returns:
            dict: Process status information, or None if not found
        """
        try:
            # Check in-memory first
            if process_id in self.processes:
                process_info = self.processes[process_id]
                return {
                    'process_id': process_info.process_id,
                    'pdb_id': process_info.pdb_id,
                    'chain_id': process_info.chain_id,
                    'batch_id': process_info.batch_id,
                    'stage': process_info.stage.value,
                    'status': process_info.status.value,
                    'progress': process_info.progress,
                    'start_time': process_info.start_time.isoformat() if process_info.start_time else None,
                    'end_time': process_info.end_time.isoformat() if process_info.end_time else None,
                    'processing_time_seconds': process_info.processing_time.total_seconds() if process_info.processing_time else None,
                    'error_message': process_info.error_message,
                    'warnings': process_info.warnings,
                    'files_created': process_info.files_created,
                    'retry_count': process_info.retry_count,
                    'metadata': process_info.metadata,
                    'last_updated': process_info.last_updated.isoformat()
                }

            # Fall back to database if available
            if self.db_available:
                return self._get_process_from_database(process_id)

            return None

        except Exception as e:
            self.logger.error(f"Error getting process status {process_id}: {str(e)}")
            return None

    def get_batch_status(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """
        Get status information for a batch.

        Args:
            batch_id: Batch identifier

        Returns:
            dict: Batch status information, or None if not found
        """
        try:
            batch_processes = [p for p in self.processes.values() if p.batch_id == batch_id]

            if not batch_processes and self.db_available:
                # Try to get from database
                return self._get_batch_from_database(batch_id)

            if not batch_processes:
                return None

            # Calculate batch statistics
            total = len(batch_processes)
            completed = sum(1 for p in batch_processes if p.status == ProcessStatus.COMPLETED)
            failed = sum(1 for p in batch_processes if p.status == ProcessStatus.FAILED)
            running = sum(1 for p in batch_processes if p.is_active)

            start_times = [p.start_time for p in batch_processes if p.start_time]
            end_times = [p.end_time for p in batch_processes if p.end_time and p.is_finished]

            return {
                'batch_id': batch_id,
                'total_processes': total,
                'completed': completed,
                'failed': failed,
                'running': running,
                'success_rate': (completed / total * 100) if total > 0 else 0,
                'start_time': min(start_times).isoformat() if start_times else None,
                'end_time': max(end_times).isoformat() if end_times and len(end_times) == (completed + failed) else None,
                'processes': [p.process_id for p in batch_processes]
            }

        except Exception as e:
            self.logger.error(f"Error getting batch status {batch_id}: {str(e)}")
            return None

    def get_performance_metrics(self) -> Dict[str, Any]:
        """Get current performance metrics"""
        try:
            self._update_performance_metrics()
            return self.metrics.to_dict()
        except Exception as e:
            self.logger.error(f"Error getting performance metrics: {str(e)}")
            return {}

    def get_error_summary(self, hours: int = 24) -> Dict[str, Any]:
        """
        Get error summary for the specified time period.

        Args:
            hours: Number of hours to look back

        Returns:
            dict: Error summary statistics
        """
        try:
            cutoff = datetime.now() - timedelta(hours=hours)
            recent_errors = [e for e in self.error_history if e['timestamp'] > cutoff]

            # Group errors by stage
            errors_by_stage = defaultdict(list)
            for error in recent_errors:
                errors_by_stage[error['stage']].append(error)

            # Count error types (simplified)
            error_types = defaultdict(int)
            for error in recent_errors:
                # Simple error classification
                error_msg = error['error'].lower()
                if 'file not found' in error_msg or 'no such file' in error_msg:
                    error_types['file_not_found'] += 1
                elif 'xml' in error_msg or 'parsing' in error_msg:
                    error_types['xml_parsing'] += 1
                elif 'database' in error_msg or 'connection' in error_msg:
                    error_types['database'] += 1
                else:
                    error_types['other'] += 1

            return {
                'time_period_hours': hours,
                'total_errors': len(recent_errors),
                'errors_by_stage': {stage: len(errors) for stage, errors in errors_by_stage.items()},
                'error_types': dict(error_types),
                'error_rate_per_hour': len(recent_errors) / hours if hours > 0 else 0
            }

        except Exception as e:
            self.logger.error(f"Error getting error summary: {str(e)}")
            return {}

    def cleanup_old_processes(self, older_than: Optional[timedelta] = None) -> int:
        """
        Clean up old completed processes from memory.

        Args:
            older_than: Remove processes older than this, defaults to config value

        Returns:
            int: Number of processes cleaned up
        """
        try:
            if older_than is None:
                older_than = self.cleanup_completed_after

            cutoff = datetime.now() - older_than

            to_remove = []
            for process_id, process_info in self.processes.items():
                if (process_info.is_finished and
                    process_info.last_updated < cutoff):
                    to_remove.append(process_id)

            for process_id in to_remove:
                del self.processes[process_id]

            if to_remove:
                self.logger.info(f"Cleaned up {len(to_remove)} old processes")

            return len(to_remove)

        except Exception as e:
            self.logger.error(f"Error cleaning up old processes: {str(e)}")
            return 0

    def _update_performance_metrics(self) -> None:
        """Update performance metrics based on recent activity"""
        try:
            now = datetime.now()

            # Only update metrics periodically
            if now - self._last_metrics_update < timedelta(seconds=30):
                return

            # Calculate processes per hour
            recent_completed = [p for p in self.processes.values()
                             if p.status == ProcessStatus.COMPLETED and
                                p.end_time and p.end_time > now - timedelta(hours=1)]

            self.metrics.processes_per_hour = len(recent_completed)

            # Calculate average processing time
            if self._operation_times:
                avg_seconds = sum(self._operation_times) / len(self._operation_times)
                self.metrics.average_processing_time = timedelta(seconds=avg_seconds)

            # Calculate error rate
            recent_processes = [p for p in self.processes.values()
                             if p.is_finished and p.last_updated > now - timedelta(hours=1)]

            if recent_processes:
                failed_count = sum(1 for p in recent_processes if p.status == ProcessStatus.FAILED)
                self.metrics.error_rate = (failed_count / len(recent_processes)) * 100

            # Update queue depth
            self.metrics.queue_depth = sum(1 for p in self.processes.values() if p.is_active)

            self._last_metrics_update = now

        except Exception as e:
            self.logger.error(f"Error updating performance metrics: {str(e)}")

    # Database interaction methods (implement based on your schema)
    def _write_process_to_database(self, process_info: ProcessInfo) -> bool:
        """Write process info to database"""
        if not self.db_available:
            return False

        try:
            query = """
                INSERT INTO process_status
                (process_id, pdb_id, chain_id, batch_id, stage, status, start_time,
                 progress, error_message, retry_count, metadata, last_updated)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (process_id) DO UPDATE SET
                    stage = EXCLUDED.stage,
                    status = EXCLUDED.status,
                    progress = EXCLUDED.progress,
                    error_message = EXCLUDED.error_message,
                    retry_count = EXCLUDED.retry_count,
                    metadata = EXCLUDED.metadata,
                    last_updated = EXCLUDED.last_updated
            """

            params = (
                process_info.process_id, process_info.pdb_id, process_info.chain_id,
                process_info.batch_id, process_info.stage.value, process_info.status.value,
                process_info.start_time, process_info.progress, process_info.error_message,
                process_info.retry_count, json.dumps(process_info.metadata),
                process_info.last_updated
            )

            self.db_manager.execute_query(query, params)
            return True

        except Exception as e:
            self.logger.warning(f"Error writing process to database: {str(e)}")
            self.db_available = False
            return False

    def _update_process_in_database(self, process_id: int, stage: ProcessStage,
                                  status: ProcessStatus, progress: Optional[float],
                                  metadata: Optional[Dict[str, Any]]) -> bool:
        """Update process in database"""
        if not self.db_available:
            return False

        try:
            query = """
                UPDATE process_status
                SET stage = %s, status = %s, last_updated = %s
            """
            params = [stage.value, status.value, datetime.now()]

            if progress is not None:
                query += ", progress = %s"
                params.append(progress)

            if metadata:
                query += ", metadata = metadata || %s"
                params.append(json.dumps(metadata))

            query += " WHERE process_id = %s"
            params.append(process_id)

            self.db_manager.execute_query(query, params)
            return True

        except Exception as e:
            self.logger.warning(f"Error updating process in database: {str(e)}")
            self.db_available = False
            return False

    def _register_file_in_database(self, process_id: int, file_path: str, relative_path: str) -> bool:
        """Register file in database"""
        if not self.db_available:
            return False

        try:
            query = """
                INSERT INTO domain_files (process_id, file_path, relative_path, created_at)
                VALUES (%s, %s, %s, %s)
                ON CONFLICT (process_id) DO UPDATE SET
                    file_path = EXCLUDED.file_path,
                    relative_path = EXCLUDED.relative_path,
                    updated_at = %s
            """

            now = datetime.now()
            params = (process_id, file_path, relative_path, now, now)

            self.db_manager.execute_query(query, params)
            return True

        except Exception as e:
            self.logger.warning(f"Error registering file in database: {str(e)}")
            self.db_available = False
            return False

    def _get_process_from_database(self, process_id: int) -> Optional[Dict[str, Any]]:
        """Get process from database"""
        if not self.db_available:
            return None

        try:
            query = """
                SELECT process_id, pdb_id, chain_id, batch_id, stage, status,
                       start_time, end_time, progress, error_message, retry_count,
                       metadata, last_updated
                FROM process_status
                WHERE process_id = %s
            """

            if hasattr(self.db_manager, 'execute_dict_query'):
                results = self.db_manager.execute_dict_query(query, (process_id,))
                return results[0] if results else None
            else:
                cursor = self.db_manager.execute(query, (process_id,))
                row = cursor.fetchone()
                if row:
                    columns = [desc[0] for desc in cursor.description]
                    return dict(zip(columns, row))
                return None

        except Exception as e:
            self.logger.warning(f"Error getting process from database: {str(e)}")
            return None

    def _get_batch_from_database(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch info from database"""
        if not self.db_available:
            return None

        try:
            query = """
                SELECT batch_id,
                       COUNT(*) as total_processes,
                       COUNT(*) FILTER (WHERE status = 'completed') as completed,
                       COUNT(*) FILTER (WHERE status = 'failed') as failed,
                       COUNT(*) FILTER (WHERE status IN ('pending', 'running', 'retrying')) as running,
                       MIN(start_time) as start_time,
                       MAX(CASE WHEN status IN ('completed', 'failed') THEN last_updated END) as end_time
                FROM process_status
                WHERE batch_id = %s
                GROUP BY batch_id
            """

            if hasattr(self.db_manager, 'execute_dict_query'):
                results = self.db_manager.execute_dict_query(query, (batch_id,))
                return results[0] if results else None
            else:
                cursor = self.db_manager.execute(query, (batch_id,))
                row = cursor.fetchone()
                if row:
                    columns = [desc[0] for desc in cursor.description]
                    return dict(zip(columns, row))
                return None

        except Exception as e:
            self.logger.warning(f"Error getting batch from database: {str(e)}")
            return None

    def is_database_available(self) -> bool:
        """Check if database is currently available"""
        return self.db_available

    def reconnect_database(self) -> bool:
        """Attempt to reconnect to database"""
        self.logger.info("Manual database reconnection attempt...")
        self.db_available = self._test_database_connection()
        self.last_db_check = datetime.now()

        if self.db_available:
            self.logger.info("Database reconnection successful")
            self._sync_in_memory_to_database()
        else:
            self.logger.warning("Database reconnection failed")

        return self.db_available

    def get_summary_stats(self) -> Dict[str, Any]:
        """Get overall summary statistics"""
        try:
            active_processes = sum(1 for p in self.processes.values() if p.is_active)
            completed_today = sum(1 for p in self.processes.values()
                                if p.status == ProcessStatus.COMPLETED and
                                   p.end_time and p.end_time > datetime.now() - timedelta(days=1))

            return {
                'total_processes_tracked': len(self.processes),
                'active_processes': active_processes,
                'completed_today': completed_today,
                'database_available': self.db_available,
                'recent_errors': len([e for e in self.error_history
                                    if e['timestamp'] > datetime.now() - timedelta(hours=1)]),
                'performance_metrics': self.get_performance_metrics(),
                'last_updated': datetime.now().isoformat()
            }

        except Exception as e:
            self.logger.error(f"Error getting summary stats: {str(e)}")
            return {}
