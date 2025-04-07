#!/usr/bin/env python3
"""
db_interaction.py - Script to interact with the PyECOD database
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

# Add parent directory to path if needed
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from ecod.core.config import ConfigManager
from ecod.core.db_manager import DBManager
from ecod.core.models import Protein, Batch, ProcessStatus

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

class DatabaseInteractor:
    """Class to interact with the PyECOD database"""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize with configuration"""
        self.config_manager = ConfigManager(config_path)
        self.db_config = self.config_manager.get_db_config()
        self.db = DBManager(self.db_config)
        self.logger = logging.getLogger("ecod.db_interaction")
        
    def list_schemas(self) -> List[str]:
        """List all schemas in the database"""
        query = """
        SELECT schema_name 
        FROM information_schema.schemata
        WHERE schema_name NOT LIKE 'pg_%' 
          AND schema_name != 'information_schema'
        ORDER BY schema_name
        """
        rows = self.db.execute_query(query)
        return [row[0] for row in rows]
    
    def list_tables(self, schema: str) -> List[str]:
        """List all tables in a schema"""
        query = """
        SELECT table_name
        FROM information_schema.tables
        WHERE table_schema = %s
        ORDER BY table_name
        """
        rows = self.db.execute_query(query, (schema,))
        return [row[0] for row in rows]
    
    def list_proteins(self, limit: int = 10) -> List[Protein]:
        """List proteins in the database"""
        query = """
        SELECT 
            p.id, p.pdb_id, p.chain_id, p.source_id, p.length, ps.sequence
        FROM 
            ecod_schema.protein p
        LEFT JOIN
            ecod_schema.protein_sequence ps ON p.id = ps.protein_id
        ORDER BY 
            p.id
        LIMIT %s
        """
        rows = self.db.execute_dict_query(query, (limit,))
        return [
            Protein(
                id=row['id'],
                pdb_id=row['pdb_id'],
                chain_id=row['chain_id'],
                source_id=row['source_id'],
                length=row['length'],
                sequence=row.get('sequence')
            )
            for row in rows
        ]
    
    def list_batches(self) -> List[Batch]:
        """List processing batches"""
        query = """
        SELECT 
            id, batch_name, base_path, type, ref_version, 
            total_items, completed_items, status,
            created_at, completed_at
        FROM 
            ecod_schema.batch
        ORDER BY 
            id DESC
        """
        rows = self.db.execute_dict_query(query)
        return [Batch.from_db_row(row) for row in rows]
    
    def count_items_by_status(self) -> Dict[str, int]:
        """Count process items by status"""
        query = """
        SELECT 
            status, COUNT(*) as count
        FROM 
            ecod_schema.process_status
        GROUP BY 
            status
        ORDER BY 
            count DESC
        """
        rows = self.db.execute_dict_query(query)
        return {row['status']: row['count'] for row in rows}
    
    def test_connection(self) -> bool:
        """Test database connection"""
        try:
            with self.db.get_connection() as conn:
                self.logger.info("Database connection successful!")
                return True
        except Exception as e:
            self.logger.error(f"Database connection failed: {e}")
            return False

def main():
    parser = argparse.ArgumentParser(description='PyECOD Database Interaction Tool')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Connection test command
    subparsers.add_parser('test', help='Test database connection')
    
    # List schemas command
    subparsers.add_parser('schemas', help='List all schemas')
    
    # List tables command
    tables_parser = subparsers.add_parser('tables', help='List tables in schema')
    tables_parser.add_argument('schema', type=str, help='Schema name')
    
    # List proteins command
    proteins_parser = subparsers.add_parser('proteins', help='List proteins')
    proteins_parser.add_argument('--limit', type=int, default=10, help='Maximum proteins to list')
    
    # List batches command
    subparsers.add_parser('batches', help='List processing batches')
    
    # Status count command
    subparsers.add_parser('status', help='Count items by status')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    db_interactor = DatabaseInteractor(args.config)
    
    if args.command == 'test':
        db_interactor.test_connection()
    
    elif args.command == 'schemas':
        schemas = db_interactor.list_schemas()
        print(f"Found {len(schemas)} schemas:")
        for schema in schemas:
            print(f"  - {schema}")
    
    elif args.command == 'tables':
        tables = db_interactor.list_tables(args.schema)
        print(f"Found {len(tables)} tables in schema '{args.schema}':")
        for table in tables:
            print(f"  - {table}")
    
    elif args.command == 'proteins':
        proteins = db_interactor.list_proteins(args.limit)
        print(f"Found {len(proteins)} proteins:")
        for protein in proteins:
            print(f"  - PDB: {protein.pdb_id}, Chain: {protein.chain_id}, Length: {protein.length}")
    
    elif args.command == 'batches':
        batches = db_interactor.list_batches()
        print(f"Found {len(batches)} processing batches:")
        for batch in batches:
            progress = f"{batch.completed_items}/{batch.total_items}"
            print(f"  - {batch.id}: {batch.batch_name} ({batch.type}) - {progress} completed, Status: {batch.status}")
    
    elif args.command == 'status':
        status_counts = db_interactor.count_items_by_status()
        print("Process status counts:")
        for status, count in status_counts.items():
            print(f"  - {status}: {count}")
    
    else:
        parser.print_help()

if __name__ == "__main__":
    main()