s# ecod/db/migration_manager.py
import os
import logging
import psycopg2
from typing import Dict, Any, List, Optional

class MigrationManager:
    def __init__(self, db_config: Dict[str, Any], migrations_dir: str):
        self.db_config = db_config
        self.migrations_dir = migrations_dir
        self.logger = logging.getLogger("ecod.migration")
        
    def _create_migration_table(self, conn) -> None:
        """Create migration tracking table if it doesn't exist"""
        with conn.cursor() as cursor:
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS ecod_schema.schema_migrations (
                id SERIAL PRIMARY KEY,
                migration_name VARCHAR(255) UNIQUE NOT NULL,
                applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """)
            conn.commit()
    
    def _get_applied_migrations(self, conn) -> List[str]:
        """Get list of already applied migrations"""
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT migration_name FROM ecod_schema.schema_migrations
            ORDER BY id
            """)
            return [row[0] for row in cursor.fetchall()]
    
    def _get_migration_files(self) -> List[str]:
        """Get sorted list of migration files"""
        files = [f for f in os.listdir(self.migrations_dir) 
               if f.endswith('.sql') and f[0].isdigit()]
        return sorted(files)
    
    def _apply_migration(self, conn, migration_file: str) -> None:
        """Apply a single migration file"""
        file_path = os.path.join(self.migrations_dir, migration_file)
        self.logger.info(f"Applying migration: {migration_file}")
        
        try:
            with open(file_path, 'r') as f:
                sql = f.read()
                
            with conn.cursor() as cursor:
                cursor.execute(sql)
                
            # Record the migration
            with conn.cursor() as cursor:
                cursor.execute(
                    "INSERT INTO ecod_schema.schema_migrations (migration_name) VALUES (%s)",
                    (migration_file,)
                )
                
            conn.commit()
            self.logger.info(f"Migration applied successfully: {migration_file}")
            
        except Exception as e:
            conn.rollback()
            self.logger.error(f"Error applying migration {migration_file}: {str(e)}")
            raise
    
    def apply_migrations(self) -> None:
        """Apply all pending migrations"""
        conn = None
        try:
            conn = psycopg2.connect(**self.db_config)
            
            # Create schema and migration table
            with conn.cursor() as cursor:
                cursor.execute("CREATE SCHEMA IF NOT EXISTS ecod_schema")
                conn.commit()
                
            self._create_migration_table(conn)
            
            # Get applied and available migrations
            applied = self._get_applied_migrations(conn)
            available = self._get_migration_files()
            
            # Apply pending migrations
            for migration in available:
                if migration not in applied:
                    self._apply_migration(conn, migration)
            
        except Exception as e:
            self.logger.error(f"Migration failed: {str(e)}")
            raise
        finally:
            if conn:
                conn.close()