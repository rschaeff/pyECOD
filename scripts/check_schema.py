#!/usr/bin/env python3
"""
Schema validation CLI tool for pyECOD
Provides command-line interface for schema checking and validation
"""

import sys
import argparse
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def main():
    parser = argparse.ArgumentParser(description="pyECOD Schema Validation Tool")
    parser.add_argument('--config', default='config/config.yml', help='Config file path')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Columns command
    cols_parser = subparsers.add_parser('columns', help='Check table columns')
    cols_parser.add_argument('schema', help='Schema name (e.g., ecod_schema)')
    cols_parser.add_argument('table', help='Table name')
    cols_parser.add_argument('--check', nargs='+', help='Columns to validate')
    
    # Tables command
    tables_parser = subparsers.add_parser('tables', help='List tables in schema')
    tables_parser.add_argument('schema', help='Schema name')
    
    # Refresh command
    refresh_parser = subparsers.add_parser('refresh', help='Refresh schema cache')
    refresh_parser.add_argument('--schema', help='Specific schema to refresh')
    refresh_parser.add_argument('--table', help='Specific table to refresh (requires --schema)')
    
    # Query command
    query_parser = subparsers.add_parser('query', help='Validate query safety')
    query_parser.add_argument('query_file', help='File containing SQL query')
    
    # Schemas command
    schemas_parser = subparsers.add_parser('schemas', help='List available schemas')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    try:
        from ecod.core.context import ApplicationContext
        from ecod.db import DBManager
        
        # Initialize context and database manager
        context = ApplicationContext(args.config)
        db = DBManager(context.config['database'])
        
        if args.command == 'columns':
            return handle_columns_command(db, args)
        elif args.command == 'tables':
            return handle_tables_command(db, args)
        elif args.command == 'refresh':
            return handle_refresh_command(db, args)
        elif args.command == 'query':
            return handle_query_command(db, args)
        elif args.command == 'schemas':
            return handle_schemas_command(db, args)
        
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


def handle_columns_command(db, args):
    """Handle the columns command"""
    if not db.table_exists(args.schema, args.table):
        print(f"❌ Table {args.schema}.{args.table} does not exist")
        return 1
    
    # Get all columns
    columns = db.get_table_columns_cached(args.schema, args.table)
    print(f"📋 Columns in {args.schema}.{args.table}:")
    for i, col in enumerate(columns, 1):
        print(f"   {i:2d}. {col}")
    
    print(f"\nTotal: {len(columns)} columns")
    
    # Validate specific columns if requested
    if args.check:
        print(f"\n✅ Column Validation:")
        validation = db.validate_columns_exist(args.schema, args.table, args.check)
        
        for col, exists in validation.items():
            status = "✅" if exists else "❌"
            print(f"   {status} {col}")
        
        # Show safe columns
        safe_columns = db.get_safe_columns(args.schema, args.table, args.check)
        missing_columns = [col for col in args.check if col not in safe_columns]
        
        if safe_columns:
            print(f"\n✅ Safe columns for queries: {', '.join(safe_columns)}")
        if missing_columns:
            print(f"❌ Missing columns: {', '.join(missing_columns)}")
    
    return 0


def handle_tables_command(db, args):
    """Handle the tables command"""
    if not db.schema_exists(args.schema):
        print(f"❌ Schema {args.schema} does not exist")
        return 1
    
    # Get tables using information_schema
    query = """
    SELECT table_name, table_type
    FROM information_schema.tables 
    WHERE table_schema = %s 
    ORDER BY table_type, table_name
    """
    
    try:
        table_rows = db.execute_dict_query(query, (args.schema,))
        
        if not table_rows:
            print(f"📋 No tables found in schema {args.schema}")
            return 0
        
        print(f"📋 Tables in {args.schema}:")
        
        # Group by table type
        base_tables = [row for row in table_rows if row['table_type'] == 'BASE TABLE']
        views = [row for row in table_rows if row['table_type'] == 'VIEW']
        
        if base_tables:
            print(f"\n  📊 Base Tables ({len(base_tables)}):")
            for i, table in enumerate(base_tables, 1):
                # Get column count
                col_count = len(db.get_table_columns_cached(args.schema, table['table_name']))
                print(f"     {i:2d}. {table['table_name']} ({col_count} columns)")
        
        if views:
            print(f"\n  👁️  Views ({len(views)}):")
            for i, view in enumerate(views, 1):
                print(f"     {i:2d}. {view['table_name']}")
        
        print(f"\nTotal: {len(table_rows)} objects")
        
    except Exception as e:
        print(f"❌ Error querying tables: {e}")
        return 1
    
    return 0


def handle_refresh_command(db, args):
    """Handle the refresh command"""
    if args.schema and args.table:
        # Refresh specific table
        if not db.table_exists(args.schema, args.table):
            print(f"❌ Table {args.schema}.{args.table} does not exist")
            return 1
        
        db.refresh_column_cache(args.schema, args.table)
        print(f"✅ Refreshed cache for {args.schema}.{args.table}")
        
    elif args.schema:
        # Refresh specific schema
        if not db.schema_exists(args.schema):
            print(f"❌ Schema {args.schema} does not exist")
            return 1
        
        db.refresh_column_cache(args.schema)
        print(f"✅ Refreshed cache for schema {args.schema}")
        
    else:
        # Refresh all common schemas
        db.refresh_column_cache()
        print("✅ Refreshed cache for all common schemas")
    
    return 0


def handle_query_command(db, args):
    """Handle the query command"""
    if not Path(args.query_file).exists():
        print(f"❌ Query file not found: {args.query_file}")
        return 1
    
    try:
        with open(args.query_file, 'r') as f:
            query = f.read().strip()
        
        if not query:
            print("❌ Query file is empty")
            return 1
        
        print(f"🔍 Validating query from {args.query_file}")
        print(f"Query preview: {query[:100]}{'...' if len(query) > 100 else ''}")
        
        # Validate query safety
        safety = db.validate_query_safety(query)
        
        if safety['safe']:
            print("✅ Query validation passed")
            print("   All referenced tables exist")
        else:
            print("❌ Query validation failed")
            for issue in safety['issues']:
                print(f"   - {issue}")
        
        if safety['warnings']:
            print("⚠️  Warnings:")
            for warning in safety['warnings']:
                print(f"   - {warning}")
        
        return 0 if safety['safe'] else 1
        
    except Exception as e:
        print(f"❌ Error validating query: {e}")
        return 1


def handle_schemas_command(db, args):
    """Handle the schemas command"""
    try:
        # Get all schemas
        query = """
        SELECT schema_name 
        FROM information_schema.schemata 
        WHERE schema_name NOT IN ('information_schema', 'pg_catalog', 'pg_toast')
        ORDER BY schema_name
        """
        
        schema_rows = db.execute_query(query)
        schemas = [row[0] for row in schema_rows]
        
        if not schemas:
            print("📋 No user schemas found")
            return 0
        
        print("📋 Available Schemas:")
        
        for i, schema in enumerate(schemas, 1):
            # Get table count for each schema
            table_query = """
            SELECT COUNT(*) 
            FROM information_schema.tables 
            WHERE table_schema = %s AND table_type = 'BASE TABLE'
            """
            
            table_count_result = db.execute_query(table_query, (schema,))
            table_count = table_count_result[0][0] if table_count_result else 0
            
            # Check if schema exists and is accessible
            exists_status = "✅" if db.schema_exists(schema) else "❌"
            
            print(f"   {i:2d}. {exists_status} {schema} ({table_count} tables)")
        
        print(f"\nTotal: {len(schemas)} schemas")
        
    except Exception as e:
        print(f"❌ Error querying schemas: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
