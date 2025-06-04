#!/usr/bin/env python3
"""
PostgreSQL Tracking Modifications for filesystem_batch_processor.py

Replace these methods in your existing filesystem_batch_processor.py:
"""

import psycopg2
import psycopg2.extras

class FilesystemBatchProcessor:
    """Modified to use PostgreSQL instead of SQLite for tracking"""

    def __init__(self, config_path: str = "config/config.local.yml"):
        self.config = self._load_config(config_path)
        self.mini_executable = self._find_mini_executable()

        # Initialize PostgreSQL connections
        self.db_conn = self._init_db_connection()
        self.tracking_conn = self._init_tracking_connection()

        # Initialize tracking table
        self._init_tracking_table()

    def _init_tracking_connection(self):
        """Initialize separate PostgreSQL connection for tracking"""
        try:
            db_config = self.config['database']
            conn = psycopg2.connect(**db_config)
            conn.autocommit = True  # Auto-commit for tracking operations
            print("✓ Connected to PostgreSQL tracking")
            return conn
        except Exception as e:
            print(f"❌ Tracking connection failed: {e}")
            raise

    def _init_tracking_table(self):
        """Ensure tracking table exists and is ready"""
        try:
            with self.tracking_conn.cursor() as cursor:
                # Verify table exists
                cursor.execute("""
                    SELECT COUNT(*) FROM mini_production.tracking_status
                """)
                print(f"✓ Tracking table ready ({cursor.fetchone()[0]} existing records)")
        except Exception as e:
            print(f"⚠️  Tracking table issue: {e}")

    def _save_protein_job(self, protein_job: ProteinJob):
        """Save protein job to PostgreSQL tracking table"""
        try:
            with self.tracking_conn.cursor() as cursor:
                cursor.execute("""
                    INSERT INTO mini_production.tracking_status
                    (protein_id, slurm_job_id, status, submitted_at, error_message)
                    VALUES (%s, %s, %s, NOW(), %s)
                    ON CONFLICT (protein_id) DO UPDATE SET
                        slurm_job_id = EXCLUDED.slurm_job_id,
                        status = EXCLUDED.status,
                        submitted_at = EXCLUDED.submitted_at,
                        error_message = EXCLUDED.error_message
                """, (
                    protein_job.protein_id,
                    protein_job.slurm_job_id,
                    protein_job.status,
                    protein_job.error_message
                ))
                # No need to commit - using autocommit=True
        except Exception as e:
            print(f"⚠️  Failed to save tracking for {protein_job.protein_id}: {e}")

    def mark_job_completed(self, protein_id: str, domains_found: int = 0):
        """Mark job as completed with domain count"""
        try:
            with self.tracking_conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE mini_production.tracking_status
                    SET status = 'completed',
                        completed_at = NOW(),
                        domains_found = %s
                    WHERE protein_id = %s
                """, (domains_found, protein_id))
        except Exception as e:
            print(f"⚠️  Failed to mark completed {protein_id}: {e}")

    def mark_job_failed(self, protein_id: str, error_message: str):
        """Mark job as failed with error message"""
        try:
            with self.tracking_conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE mini_production.tracking_status
                    SET status = 'failed',
                        completed_at = NOW(),
                        error_message = %s
                    WHERE protein_id = %s
                """, (error_message, protein_id))
        except Exception as e:
            print(f"⚠️  Failed to mark failed {protein_id}: {e}")

    def get_statistics(self) -> Dict:
        """Get current processing statistics from PostgreSQL"""
        try:
            with self.tracking_conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                cursor.execute("""
                    SELECT
                        status,
                        COUNT(*) as count,
                        AVG(EXTRACT(EPOCH FROM (completed_at - submitted_at))) as avg_runtime_seconds
                    FROM mini_production.tracking_status
                    WHERE submitted_at >= NOW() - INTERVAL '24 hours'  -- Last 24 hours only
                    GROUP BY status
                """)

                results = cursor.fetchall()
                status_counts = {row['status']: row['count'] for row in results}

                # Get average runtime for completed jobs
                avg_runtime = next((row['avg_runtime_seconds'] for row in results
                                  if row['status'] == 'completed'), 0) or 0

                total = sum(status_counts.values())

                return {
                    'total': total,
                    'completed': status_counts.get('completed', 0),
                    'submitted': status_counts.get('submitted', 0),
                    'failed': status_counts.get('failed', 0),
                    'pending': 0,  # No pending in PostgreSQL tracking
                    'avg_runtime': avg_runtime
                }
        except Exception as e:
            print(f"⚠️  Failed to get statistics: {e}")
            return {'total': 0, 'completed': 0, 'submitted': 0, 'failed': 0, 'pending': 0}

    def check_already_processed(self, protein_id: str) -> bool:
        """Check if protein was already successfully processed"""
        try:
            with self.tracking_conn.cursor() as cursor:
                cursor.execute("""
                    SELECT status FROM mini_production.tracking_status
                    WHERE protein_id = %s AND status = 'completed'
                """, (protein_id,))
                return cursor.fetchone() is not None
        except Exception as e:
            print(f"⚠️  Failed to check processed status for {protein_id}: {e}")
            return False

    def update_job_statuses_from_slurm(self):
        """Update job statuses by checking SLURM job states"""
        try:
            # Get all submitted jobs from tracking
            with self.tracking_conn.cursor() as cursor:
                cursor.execute("""
                    SELECT protein_id, slurm_job_id
                    FROM mini_production.tracking_status
                    WHERE status = 'submitted' AND slurm_job_id IS NOT NULL
                """)
                submitted_jobs = cursor.fetchall()

            if not submitted_jobs:
                return

            # Check SLURM status for these jobs
            job_ids = [job[1] for job in submitted_jobs]
            job_ids_str = ','.join(job_ids)

            # Query SLURM accounting database
            result = subprocess.run([
                "sacct", "-j", job_ids_str, "--format=JobID,State",
                "--noheader", "--parsable2"
            ], capture_output=True, text=True)

            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if '|' in line:
                        job_id, state = line.strip().split('|')

                        # Find corresponding protein
                        protein_id = next((job[0] for job in submitted_jobs if job[1] == job_id), None)

                        if protein_id:
                            if state in ['COMPLETED']:
                                self.mark_job_completed(protein_id)
                            elif state in ['FAILED', 'CANCELLED', 'TIMEOUT']:
                                self.mark_job_failed(protein_id, f"SLURM job {state}")

        except Exception as e:
            print(f"⚠️  Failed to update job statuses: {e}")

    def cleanup_tracking(self, older_than_hours: int = 48):
        """Clean up old tracking records"""
        try:
            with self.tracking_conn.cursor() as cursor:
                cursor.execute("""
                    DELETE FROM mini_production.tracking_status
                    WHERE submitted_at < NOW() - INTERVAL '%s hours'
                    AND status IN ('completed', 'failed')
                """, (older_than_hours,))
                print(f"✓ Cleaned up tracking records older than {older_than_hours} hours")
        except Exception as e:
            print(f"⚠️  Failed to cleanup tracking: {e}")


# Additional monitoring functions
def show_tracking_summary():
    """Standalone function to show tracking summary"""
    import yaml

    with open("config/config.local.yml", 'r') as f:
        config = yaml.safe_load(f)

    conn = psycopg2.connect(**config['database'])

    with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
        # Overall summary
        cursor.execute("""
            SELECT
                status,
                COUNT(*) as count,
                MIN(submitted_at) as first_submitted,
                MAX(submitted_at) as last_submitted
            FROM mini_production.tracking_status
            GROUP BY status
            ORDER BY count DESC
        """)

        print("📊 Tracking Summary:")
        print("-" * 60)
        for row in cursor.fetchall():
            print(f"  {row['status']:<12} {row['count']:>8,} jobs")

        # Recent activity
        cursor.execute("""
            SELECT
                COUNT(*) as recent_submissions,
                COUNT(*) FILTER (WHERE status = 'completed') as recent_completions
            FROM mini_production.tracking_status
            WHERE submitted_at >= NOW() - INTERVAL '1 hour'
        """)

        recent = cursor.fetchone()
        print(f"\n🕐 Last Hour Activity:")
        print(f"  Submitted: {recent['recent_submissions']}")
        print(f"  Completed: {recent['recent_completions']}")


if __name__ == "__main__":
    # Test the tracking
    show_tracking_summary()
