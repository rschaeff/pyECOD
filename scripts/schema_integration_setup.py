#!/usr/bin/env python3
"""
Schema Manager Integration Setup for pyECOD

This script shows how to integrate schema awareness into your existing pyECOD architecture.
"""

import os
import sys
from pathlib import Path

# Fix the import issue and create the schema manager module
def create_schema_manager_module():
    """Create the schema manager module with proper imports"""
    
    module_content = '''#!/usr/bin/env python3
"""
Schema Manager for pyECOD - Database Schema Awareness

This module provides database schema introspection and management
without relying on pg_dump outputs.
"""

import os
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Union
from dataclasses import dataclass, asdict

try:
    from sqlalchemy import create_engine, MetaData, Table, text
    from sqlalchemy.engine import Engine
    from sqlalchemy.exc import SQLAlchemyError
    SQLALCHEMY_AVAILABLE = True
except ImportError:
    SQLALCHEMY_AVAILABLE = False
    print("SQLAlchemy not available - install with: pip install sqlalchemy")

@dataclass
class TableInfo:
    """Information about a database table"""
    schema: str
    name: str
    columns: Dict[str, Dict[str, Any]]
    primary_keys: List[str]
    foreign_keys: List[str]
    indexes: List[str]
    constraints: List[str] = None

    def __post_init__(self):
        if self.constraints is None:
            self.constraints = []

@dataclass
class SchemaInfo:
    """Information about a database schema"""
    name: str
    tables: Dict[str, TableInfo]
    views: Dict[str, str]
    functions: Dict[str, Dict[str, Any]]
    sequences: List[str]

class SchemaManager:
    """
    Schema-aware manager for pyECOD database schemas
    """
    
    def __init__(self, db_config: Dict[str, Any], cache_dir: str = "schema_cache"):
        if not SQLALCHEMY_AVAILABLE:
            raise ImportError("SQLAlchemy is required for schema management")
            
        self.db_config = db_config
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
        # Create SQLAlchemy engine for introspection
        self.engine = self._create_engine()
        self.metadata = MetaData()
        
        # Cache for schema information
        self._schema_cache = {}
        
        self.logger = logging.getLogger(__name__)
        
        # Load cached schema info
        self._load_schema_cache()
    
    def _create_engine(self) -> Engine:
        """Create SQLAlchemy engine from config"""
        connection_string = (
            f"postgresql://{self.db_config['user']}:{self.db_config['password']}"
            f"@{self.db_config['host']}:{self.db_config['port']}/{self.db_config['database']}"
        )
        return create_engine(connection_string, echo=False)
    
    def _load_schema_cache(self):
        """Load schema information from cache files"""
        for cache_file in self.cache_dir.glob("*.json"):
            schema_name = cache_file.stem
            try:
                with open(cache_file, 'r') as f:
                    cached_data = json.load(f)
                    
                # Convert back to SchemaInfo object
                schema_info = SchemaInfo(
                    name=schema_name,
                    tables={},
                    views=cached_data.get('views', {}),
                    functions=cached_data.get('functions', {}),
                    sequences=cached_data.get('sequences', [])
                )
                
                # Rebuild table info
                for table_name, table_data in cached_data.get('tables', {}).items():
                    schema_info.tables[table_name] = TableInfo(
                        schema=schema_name,
                        name=table_name,
                        columns=table_data['columns'],
                        primary_keys=table_data['primary_keys'],
                        foreign_keys=table_data['foreign_keys'],
                        indexes=table_data['indexes'],
                        constraints=table_data.get('constraints', [])
                    )
                
                self._schema_cache[schema_name] = schema_info
                self.logger.info(f"Loaded cached schema info for {schema_name}")
                
            except Exception as e:
                self.logger.warning(f"Failed to load cache for {schema_name}: {e}")
    
    def refresh_schema(self, schema_name: str, force: bool = False) -> SchemaInfo:
        """Refresh schema information from database"""
        if not force and schema_name in self._schema_cache:
            return self._schema_cache[schema_name]
        
        self.logger.info(f"Refreshing schema info for {schema_name}")
        
        try:
            # Introspect schema
            metadata = MetaData()
            metadata.reflect(bind=self.engine, schema=schema_name)
            
            schema_info = SchemaInfo(
                name=schema_name,
                tables={},
                views=self._get_views(schema_name),
                functions=self._get_functions(schema_name),
                sequences=self._get_sequences(schema_name)
            )
            
            # Process tables
            for table_key, table in metadata.tables.items():
                if '.' in table_key:
                    _, table_name = table_key.split('.', 1)
                else:
                    table_name = table_key
                
                # Get column information
                columns = {}
                for column in table.columns:
                    columns[column.name] = {
                        'type': str(column.type),
                        'nullable': column.nullable,
                        'default': str(column.default) if column.default else None,
                        'primary_key': column.primary_key,
                        'foreign_keys': [str(fk.target_fullname) for fk in column.foreign_keys],
                        'autoincrement': getattr(column, 'autoincrement', False)
                    }
                
                # Create table info
                table_info = TableInfo(
                    schema=schema_name,
                    name=table_name,
                    columns=columns,
                    primary_keys=[col.name for col in table.primary_key.columns],
                    foreign_keys=[str(fk) for fk in table.foreign_keys],
                    indexes=[idx.name for idx in table.indexes],
                    constraints=[]
                )
                
                schema_info.tables[table_name] = table_info
            
            # Cache the result
            self._schema_cache[schema_name] = schema_info
            self._save_schema_cache(schema_name)
            
            return schema_info
            
        except Exception as e:
            self.logger.error(f"Failed to refresh schema {schema_name}: {e}")
            # Return cached version if available
            return self._schema_cache.get(schema_name)
    
    def _save_schema_cache(self, schema_name: str):
        """Save schema information to cache"""
        if schema_name not in self._schema_cache:
            return
        
        schema_info = self._schema_cache[schema_name]
        cache_file = self.cache_dir / f"{schema_name}.json"
        
        # Convert to JSON-serializable format
        cache_data = {
            'tables': {},
            'views': schema_info.views,
            'functions': schema_info.functions,
            'sequences': schema_info.sequences
        }
        
        for table_name, table_info in schema_info.tables.items():
            cache_data['tables'][table_name] = asdict(table_info)
        
        try:
            with open(cache_file, 'w') as f:
                json.dump(cache_data, f, indent=2, default=str)
            self.logger.debug(f"Cached schema info for {schema_name}")
        except Exception as e:
            self.logger.error(f"Failed to cache schema info for {schema_name}: {e}")
    
    def _get_views(self, schema_name: str) -> Dict[str, str]:
        """Get view definitions"""
        query = text("""
            SELECT table_name, view_definition 
            FROM information_schema.views 
            WHERE table_schema = :schema_name
        """)
        
        try:
            with self.engine.connect() as conn:
                result = conn.execute(query, {'schema_name': schema_name})
                return {row[0]: row[1] for row in result}
        except Exception as e:
            self.logger.warning(f"Failed to get views for {schema_name}: {e}")
            return {}
    
    def _get_functions(self, schema_name: str) -> Dict[str, Dict[str, Any]]:
        """Get function information"""
        query = text("""
            SELECT 
                routine_name,
                routine_type,
                data_type,
                routine_definition
            FROM information_schema.routines 
            WHERE routine_schema = :schema_name
        """)
        
        try:
            with self.engine.connect() as conn:
                result = conn.execute(query, {'schema_name': schema_name})
                return {
                    row[0]: {
                        'type': row[1],
                        'return_type': row[2],
                        'definition': row[3]
                    } for row in result
                }
        except Exception as e:
            self.logger.warning(f"Failed to get functions for {schema_name}: {e}")
            return {}
    
    def _get_sequences(self, schema_name: str) -> List[str]:
        """Get sequence names"""
        query = text("""
            SELECT sequence_name 
            FROM information_schema.sequences 
            WHERE sequence_schema = :schema_name
        """)
        
        try:
            with self.engine.connect() as conn:
                result = conn.execute(query, {'schema_name': schema_name})
                return [row[0] for row in result]
        except Exception as e:
            self.logger.warning(f"Failed to get sequences for {schema_name}: {e}")
            return []
    
    def get_table_info(self, schema_name: str, table_name: str) -> Optional[TableInfo]:
        """Get information about a specific table"""
        schema_info = self.get_schema_info(schema_name)
        if not schema_info:
            return None
        return schema_info.tables.get(table_name)
    
    def get_schema_info(self, schema_name: str) -> Optional[SchemaInfo]:
        """Get schema information (cached or fresh)"""
        if schema_name not in self._schema_cache:
            return self.refresh_schema(schema_name)
        return self._schema_cache[schema_name]
    
    def table_exists(self, schema_name: str, table_name: str) -> bool:
        """Check if a table exists"""
        table_info = self.get_table_info(schema_name, table_name)
        return table_info is not None
    
    def column_exists(self, schema_name: str, table_name: str, column_name: str) -> bool:
        """Check if a column exists"""
        table_info = self.get_table_info(schema_name, table_name)
        if not table_info:
            return False
        return column_name in table_info.columns
    
    def get_table_columns(self, schema_name: str, table_name: str) -> List[str]:
        """Get list of column names for a table"""
        table_info = self.get_table_info(schema_name, table_name)
        return list(table_info.columns.keys()) if table_info else []
    
    def export_schema_summary(self, schema_name: str) -> Dict[str, Any]:
        """Export schema summary as dictionary"""
        schema_info = self.get_schema_info(schema_name)
        if not schema_info:
            return {}
        
        return {
            'schema_name': schema_name,
            'table_count': len(schema_info.tables),
            'view_count': len(schema_info.views),
            'function_count': len(schema_info.functions),
            'tables': {
                name: {
                    'column_count': len(table.columns),
                    'primary_keys': table.primary_keys,
                    'has_foreign_keys': len(table.foreign_keys) > 0
                } for name, table in schema_info.tables.items()
            }
        }

# Integration helpers
def integrate_with_existing_db_manager():
    """Integration example with existing DBManager"""
    
    integration_code = '''
# Add this to your existing DBManager class

from ecod.db.schema_manager import SchemaManager

class EnhancedDBManager(DBManager):  # Your existing DBManager
    """Enhanced DBManager with schema awareness"""
    
    def __init__(self, config):
        super().__init__(config)
        
        # Add schema manager
        self.schema_manager = SchemaManager(config['database'])
    
    def execute_query(self, query: str, params: tuple = None, validate_schema: bool = False):
        """Execute query with optional schema validation"""
        
        if validate_schema:
            # Check if query references valid tables
            issues = self._validate_query_schema(query)
            if issues:
                self.logger.warning(f"Schema validation issues: {issues}")
        
        # Your existing query execution
        return super().execute_query(query, params)
    
    def _validate_query_schema(self, query: str) -> List[str]:
        """Basic schema validation for queries"""
        issues = []
        
        # Simple table validation (could be more sophisticated)
        import re
        table_refs = re.findall(r'\\b(?:FROM|JOIN|UPDATE|INTO)\\s+(\\w+)\\.(\\w+)', query, re.IGNORECASE)
        
        for schema, table in table_refs:
            if not self.schema_manager.table_exists(schema, table):
                issues.append(f"Table {schema}.{table} does not exist")
        
        return issues
    
    def get_table_info(self, schema: str, table: str):
        """Get table metadata"""
        return self.schema_manager.get_table_info(schema, table)
    
    def refresh_schema_cache(self):
        """Refresh schema metadata cache"""
        for schema in ['ecod_schema', 'pdb_analysis']:
            self.schema_manager.refresh_schema(schema, force=True)
    '''
    
    return integration_code

if __name__ == "__main__":
    print("Schema Manager module ready for integration")
'''
    
    # Write the module
    module_path = Path("ecod/db/schema_manager.py")
    module_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(module_path, 'w') as f:
        f.write(module_content)
    
    print(f"✅ Created schema manager module: {module_path}")
    return module_path

def create_cli_tool():
    """Create CLI tool for schema management"""
    
    cli_content = '''#!/usr/bin/env python3
"""
CLI tool for pyECOD schema management
"""

import sys
import argparse
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from ecod.core.context import ApplicationContext
from ecod.db.schema_manager import SchemaManager

def main():
    parser = argparse.ArgumentParser(description="pyECOD Schema Manager")
    parser.add_argument('--config', default='config/config.yml', help='Config file')
    
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Refresh command
    refresh_parser = subparsers.add_parser('refresh', help='Refresh schema cache')
    refresh_parser.add_argument('schema', help='Schema name (ecod_schema, pdb_analysis)')
    
    # Info command
    info_parser = subparsers.add_parser('info', help='Show schema info')
    info_parser.add_argument('schema', help='Schema name')
    info_parser.add_argument('--table', help='Specific table name')
    
    # List command
    list_parser = subparsers.add_parser('list', help='List tables in schema')
    list_parser.add_argument('schema', help='Schema name')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate table exists')
    validate_parser.add_argument('schema', help='Schema name')
    validate_parser.add_argument('table', help='Table name')
    
    # Summary command
    summary_parser = subparsers.add_parser('summary', help='Schema summary')
    summary_parser.add_argument('schema', help='Schema name')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    try:
        # Initialize
        context = ApplicationContext(args.config)
        schema_manager = SchemaManager(context.config['database'])
        
        if args.command == 'refresh':
            schema_info = schema_manager.refresh_schema(args.schema, force=True)
            print(f"✅ Refreshed schema cache for {args.schema}")
            print(f"   Tables: {len(schema_info.tables)}")
            print(f"   Views: {len(schema_info.views)}")
            print(f"   Functions: {len(schema_info.functions)}")
        
        elif args.command == 'info':
            schema_info = schema_manager.get_schema_info(args.schema)
            if not schema_info:
                print(f"❌ Schema {args.schema} not found")
                return 1
            
            if args.table:
                table_info = schema_info.tables.get(args.table)
                if not table_info:
                    print(f"❌ Table {args.table} not found in {args.schema}")
                    return 1
                
                print(f"📋 Table: {args.schema}.{args.table}")
                print(f"   Columns: {len(table_info.columns)}")
                print(f"   Primary Keys: {table_info.primary_keys}")
                print(f"   Foreign Keys: {len(table_info.foreign_keys)}")
                print(f"   Indexes: {table_info.indexes}")
                
                print(f"\\n📊 Columns:")
                for col_name, col_info in table_info.columns.items():
                    pk = " (PK)" if col_info['primary_key'] else ""
                    nullable = " NULL" if col_info['nullable'] else " NOT NULL"
                    print(f"   - {col_name}: {col_info['type']}{nullable}{pk}")
            else:
                print(f"📋 Schema: {args.schema}")
                print(f"   Tables: {len(schema_info.tables)}")
                print(f"   Views: {len(schema_info.views)}")
                print(f"   Functions: {len(schema_info.functions)}")
        
        elif args.command == 'list':
            schema_info = schema_manager.get_schema_info(args.schema)
            if not schema_info:
                print(f"❌ Schema {args.schema} not found")
                return 1
            
            print(f"📋 Tables in {args.schema}:")
            for table_name in sorted(schema_info.tables.keys()):
                table_info = schema_info.tables[table_name]
                print(f"   - {table_name} ({len(table_info.columns)} columns)")
        
        elif args.command == 'validate':
            exists = schema_manager.table_exists(args.schema, args.table)
            if exists:
                print(f"✅ Table {args.schema}.{args.table} exists")
                table_info = schema_manager.get_table_info(args.schema, args.table)
                print(f"   Columns: {len(table_info.columns)}")
            else:
                print(f"❌ Table {args.schema}.{args.table} does not exist")
                return 1
        
        elif args.command == 'summary':
            summary = schema_manager.export_schema_summary(args.schema)
            if not summary:
                print(f"❌ Schema {args.schema} not found")
                return 1
            
            print(f"📊 Schema Summary: {args.schema}")
            print(f"   Tables: {summary['table_count']}")
            print(f"   Views: {summary['view_count']}")
            print(f"   Functions: {summary['function_count']}")
            
            print(f"\\n📋 Table Details:")
            for table_name, table_data in sorted(summary['tables'].items()):
                fk_status = " (has FKs)" if table_data['has_foreign_keys'] else ""
                print(f"   - {table_name}: {table_data['column_count']} columns{fk_status}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
'''
    
    cli_path = Path("scripts/schema_manager.py")
    with open(cli_path, 'w') as f:
        f.write(cli_content)
    
    # Make executable
    os.chmod(cli_path, 0o755)
    
    print(f"✅ Created CLI tool: {cli_path}")
    return cli_path

def create_integration_example():
    """Create integration example showing how to use schema awareness"""
    
    example_content = '''#!/usr/bin/env python3
"""
Example: Using Schema Awareness in pyECOD

This shows how to integrate schema awareness into your existing code.
"""

from ecod.core.context import ApplicationContext
from ecod.db.schema_manager import SchemaManager

def example_schema_aware_query():
    """Example of schema-aware database operations"""
    
    # Initialize
    context = ApplicationContext("config/config.yml")
    schema_manager = SchemaManager(context.config['database'])
    
    # Check if table exists before querying
    if schema_manager.table_exists('ecod_schema', 'algorithm_version'):
        print("✅ algorithm_version table exists")
        
        # Get table info
        table_info = schema_manager.get_table_info('ecod_schema', 'algorithm_version')
        print(f"   Columns: {list(table_info.columns.keys())}")
        
        # Safe query construction
        columns = schema_manager.get_table_columns('ecod_schema', 'algorithm_version')
        if 'version_id' in columns and 'name' in columns:
            query = "SELECT version_id, name FROM ecod_schema.algorithm_version"
            # Execute with your existing DBManager
            # results = db_manager.execute_dict_query(query)
            print(f"   Safe to execute: {query}")
        else:
            print("   ❌ Required columns not found")
    else:
        print("❌ algorithm_version table not found")

def example_validate_migration():
    """Example of validating schema before migration"""
    
    context = ApplicationContext("config/config.yml")
    schema_manager = SchemaManager(context.config['database'])
    
    # Validate expected tables exist
    required_tables = [
        ('ecod_schema', 'algorithm_version'),
        ('ecod_schema', 'algorithm_run'),
        ('pdb_analysis', 'protein'),
        ('pdb_analysis', 'domain')
    ]
    
    missing_tables = []
    for schema, table in required_tables:
        if not schema_manager.table_exists(schema, table):
            missing_tables.append(f"{schema}.{table}")
    
    if missing_tables:
        print(f"❌ Missing tables: {missing_tables}")
        return False
    else:
        print("✅ All required tables exist")
        return True

def example_dynamic_query_building():
    """Example of building queries dynamically based on schema"""
    
    context = ApplicationContext("config/config.yml")
    schema_manager = SchemaManager(context.config['database'])
    
    # Get available columns for a table
    columns = schema_manager.get_table_columns('pdb_analysis', 'protein')
    
    # Build query based on available columns
    select_columns = []
    if 'pdb_id' in columns:
        select_columns.append('pdb_id')
    if 'chain_id' in columns:
        select_columns.append('chain_id')
    if 'source_id' in columns:
        select_columns.append('source_id')
    
    if select_columns:
        query = f"SELECT {', '.join(select_columns)} FROM pdb_analysis.protein LIMIT 10"
        print(f"Dynamic query: {query}")
    else:
        print("❌ No expected columns found")

if __name__ == "__main__":
    print("🔍 Schema Awareness Examples")
    print("=" * 40)
    
    try:
        print("\\n1️⃣ Schema-Aware Query Example")
        example_schema_aware_query()
        
        print("\\n2️⃣ Migration Validation Example")
        example_validate_migration()
        
        print("\\n3️⃣ Dynamic Query Building Example")
        example_dynamic_query_building()
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
'''
    
    example_path = Path("examples/schema_awareness_example.py")
    example_path.parent.mkdir(exist_ok=True)
    
    with open(example_path, 'w') as f:
        f.write(example_content)
    
    print(f"✅ Created integration example: {example_path}")
    return example_path

def main():
    """Setup schema awareness for pyECOD"""
    print("🔧 Setting up Schema Awareness for pyECOD")
    print("=" * 50)
    
    # Check if SQLAlchemy is available
    try:
        import sqlalchemy
        print("✅ SQLAlchemy is available")
    except ImportError:
        print("❌ SQLAlchemy not found")
        print("   Install with: pip install sqlalchemy")
        print("   Optional: pip install psycopg2-binary")
        return 1
    
    # Create schema manager module
    print("\\n1️⃣ Creating Schema Manager Module")
    module_path = create_schema_manager_module()
    
    # Create CLI tool
    print("\\n2️⃣ Creating CLI Tool")
    cli_path = create_cli_tool()
    
    # Create integration example
    print("\\n3️⃣ Creating Integration Example")
    example_path = create_integration_example()
    
    # Create __init__.py files
    print("\\n4️⃣ Creating Package Structure")
    init_files = [
        "ecod/db/__init__.py"
    ]
    
    for init_file in init_files:
        Path(init_file).parent.mkdir(parents=True, exist_ok=True)
        if not Path(init_file).exists():
            with open(init_file, 'w') as f:
                f.write('"""Database utilities for pyECOD"""\n')
            print(f"✅ Created {init_file}")
    
    print("\\n🎯 Setup Complete!")
    print("=" * 30)
    
    print("\\n📋 Next Steps:")
    print("1. Test the schema manager:")
    print(f"   python {cli_path} refresh ecod_schema")
    print(f"   python {cli_path} info ecod_schema")
    print(f"   python {cli_path} list pdb_analysis")
    
    print("\\n2. Integrate with your existing DBManager:")
    print("   - Add schema_manager as a property")
    print("   - Use table_exists() before queries")
    print("   - Use get_table_columns() for dynamic queries")
    
    print("\\n3. Run the integration example:")
    print(f"   python {example_path}")
    
    print("\\n🔧 Benefits:")
    print("✅ No more pg_dump schema outputs needed")
    print("✅ Real-time schema validation")
    print("✅ Cached schema information")
    print("✅ Integration with existing code")
    print("✅ CLI tools for schema management")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
