#!/usr/bin/env python3
"""
Database manager for the ECOD pipeline
Handles database connections and queries with schema awareness
"""
import psycopg2
import psycopg2.extras
import logging
import json
from pathlib import Path
from contextlib import contextmanager
from typing import Dict, Any, List, Tuple, Optional, Generator, TypeVar, Union, Callable

from ecod.exceptions import ConnectionError, QueryError, DatabaseError

T = TypeVar('T')

class DBManager:
    """Database manager for the ECOD pipeline with schema awareness"""

    def __init__(self, config: Dict[str, Any]):
        """Initialize database manager

        Args:
            config: Database configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger("ecod.db")

        # Validate required configuration
        required_fields = ['host', 'port', 'database', 'user']
        for field in required_fields:
            if field not in config:
                raise ConnectionError(f"Missing required database configuration field: {field}")

        # Add schema awareness
        self.cache_dir = Path("schema_cache")
        self.cache_dir.mkdir(exist_ok=True)
        self._column_cache: Dict[str, List[str]] = {}
        self._load_column_cache()

    def _load_column_cache(self):
        """Load cached column information"""
        cache_file = self.cache_dir / "columns.json"
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    self._column_cache = json.load(f)
                self.logger.debug("Loaded column cache")
            except Exception as e:
                self.logger.warning(f"Failed to load column cache: {e}")

    def _save_column_cache(self):
        """Save column cache to disk"""
        cache_file = self.cache_dir / "columns.json"
        try:
            with open(cache_file, 'w') as f:
                json.dump(self._column_cache, f, indent=2)
        except Exception as e:
            self.logger.error(f"Failed to save column cache: {e}")

    @contextmanager
    def get_connection(self) -> Generator[psycopg2.extensions.connection, None, None]:
        """Context manager for database connections

        Yields:
            Database connection

        Raises:
            ConnectionError: If connection fails
        """
        conn = None
        try:
            conn = psycopg2.connect(**self.config)
            yield conn
            conn.commit()
        except psycopg2.Error as e:
            if conn:
                conn.rollback()
            error_msg = f"Database connection error: {str(e)}"
            self.logger.error(error_msg)
            raise ConnectionError(error_msg, {"code": e.pgcode} if hasattr(e, 'pgcode') else None) from e
        finally:
            if conn:
                conn.close()

    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict[str, Any]]] = None) -> List[Tuple]:
        """Execute a query and return results

        Args:
            query: SQL query
            params: Query parameters

        Returns:
            List of result tuples

        Raises:
            QueryError: If query execution fails
        """
        try:
            with self.get_connection() as conn:
                with conn.cursor() as cursor:
                    cursor.execute(query, params or ())
                    if cursor.description:  # If the query returns rows
                        return cursor.fetchall()
                    return []
        except ConnectionError:
            raise  # Re-raise connection errors
        except psycopg2.Error as e:
            error_msg = f"Query execution error: {str(e)}"
            self.logger.error(f"{error_msg}\nQuery: {query}\nParams: {params}")
            raise QueryError(error_msg, {"query": query, "params": params, "code": e.pgcode} if hasattr(e, 'pgcode') else None) from e

    def execute_dict_query(self, query: str, params: Optional[Union[Tuple, Dict[str, Any]]] = None) -> List[Dict[str, Any]]:
        """Execute a query and return results as dictionaries

        Args:
            query: SQL query
            params: Query parameters

        Returns:
            List of result dictionaries

        Raises:
            QueryError: If query execution fails
        """
        try:
            with self.get_connection() as conn:
                with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                    cursor.execute(query, params or ())
                    if cursor.description:  # If the query returns rows
                        return cursor.fetchall()
                    return []
        except ConnectionError:
            raise  # Re-raise connection errors
        except psycopg2.Error as e:
            error_msg = f"Query execution error: {str(e)}"
            self.logger.error(f"{error_msg}\nQuery: {query}\nParams: {params}")
            raise QueryError(error_msg, {"query": query, "params": params, "code": e.pgcode} if hasattr(e, 'pgcode') else None) from e

    def insert(self, table: str, data: Dict[str, Any], returning: Optional[str] = None) -> Optional[Any]:
        """Insert a record and optionally return a value

        Args:
            table: Table name
            data: Column data
            returning: Column to return (optional)

        Returns:
            Value of the returning column if specified

        Raises:
            QueryError: If insert fails
        """
        columns = list(data.keys())
        values = [data[col] for col in columns]
        placeholders = [f"%s" for _ in columns]

        query = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({', '.join(placeholders)})"

        if returning:
            query += f" RETURNING {returning}"

        try:
            with self.get_connection() as conn:
                with conn.cursor() as cursor:
                    cursor.execute(query, values)
                    if returning:
                        result = cursor.fetchone()
                        return result[0] if result else None
                    return None
        except (ConnectionError, psycopg2.Error) as e:
            if isinstance(e, psycopg2.Error):
                error_msg = f"Insert error: {str(e)}"
                self.logger.error(f"{error_msg}\nTable: {table}\nData: {data}")
                raise QueryError(error_msg, {"table": table, "data": data, "code": e.pgcode} if hasattr(e, 'pgcode') else None) from e
            raise

    def update(self, table: str, data: Dict[str, Any], condition: str,
              condition_params: Union[Tuple, List, Dict[str, Any]]) -> int:
        """Update records in a table

        Args:
            table: Table name
            data: Column data to update
            condition: WHERE clause condition
            condition_params: Parameters for the condition

        Returns:
            Number of rows updated

        Raises:
            QueryError: If update fails
        """
        set_items = [f"{k} = %s" for k in data.keys()]
        values = list(data.values())

        if isinstance(condition_params, (list, tuple)):
            values.extend(condition_params)
        elif isinstance(condition_params, dict):
            # For named parameters
            values.extend(condition_params.values())

        query = f"UPDATE {table} SET {', '.join(set_items)} WHERE {condition}"

        try:
            with self.get_connection() as conn:
                with conn.cursor() as cursor:
                    cursor.execute(query, values)
                    return cursor.rowcount
        except (ConnectionError, psycopg2.Error) as e:
            if isinstance(e, psycopg2.Error):
                error_msg = f"Update error: {str(e)}"
                self.logger.error(f"{error_msg}\nTable: {table}\nData: {data}\nCondition: {condition}")
                raise QueryError(error_msg, {"table": table, "data": data, "condition": condition, "code": e.pgcode} if hasattr(e, 'pgcode') else None) from e
            raise

    def delete(self, table: str, condition: str,
              condition_params: Union[Tuple, List, Dict[str, Any]]) -> int:
        """Delete records from a table

        Args:
            table: Table name
            condition: WHERE clause condition
            condition_params: Parameters for the condition

        Returns:
            Number of rows deleted

        Raises:
            QueryError: If delete fails
        """
        query = f"DELETE FROM {table} WHERE {condition}"

        try:
            with self.get_connection() as conn:
                with conn.cursor() as cursor:
                    cursor.execute(query, condition_params)
                    return cursor.rowcount
        except (ConnectionError, psycopg2.Error) as e:
            if isinstance(e, psycopg2.Error):
                error_msg = f"Delete error: {str(e)}"
                self.logger.error(f"{error_msg}\nTable: {table}\nCondition: {condition}")
                raise QueryError(error_msg, {"table": table, "condition": condition, "code": e.pgcode} if hasattr(e, 'pgcode') else None) from e
            raise

    def execute_transaction(self, callback: Callable[[psycopg2.extensions.cursor], T]) -> T:
        """Execute operations in a transaction

        Args:
            callback: Function that takes a cursor and performs operations

        Returns:
            Result of the callback function

        Raises:
            DatabaseError: If transaction fails
        """
        try:
            with self.get_connection() as conn:
                with conn.cursor() as cursor:
                    result = callback(cursor)
                    return result
        except (ConnectionError, psycopg2.Error) as e:
            if isinstance(e, psycopg2.Error):
                error_msg = f"Transaction error: {str(e)}"
                self.logger.error(error_msg)
                raise DatabaseError(error_msg, {"code": e.pgcode} if hasattr(e, 'pgcode') else None) from e
            raise

    def table_exists(self, schema: str, table: str) -> bool:
        """Check if a table exists

        Args:
            schema: Schema name
            table: Table name

        Returns:
            True if table exists
        """
        query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables
            WHERE table_schema = %s
            AND table_name = %s
        )
        """

        try:
            result = self.execute_query(query, (schema, table))
            return result[0][0] if result else False
        except QueryError:
            return False

    def schema_exists(self, schema: str) -> bool:
        """Check if a schema exists

        Args:
            schema: Schema name

        Returns:
            True if schema exists
        """
        query = """
        SELECT EXISTS (
            SELECT FROM information_schema.schemata
            WHERE schema_name = %s
        )
        """

        try:
            result = self.execute_query(query, (schema,))
            return result[0][0] if result else False
        except QueryError:
            return False

    def get_column_names(self, schema: str, table: str) -> List[str]:
        """Get column names for a table

        Args:
            schema: Schema name
            table: Table name

        Returns:
            List of column names
        """
        query = """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = %s
        AND table_name = %s
        ORDER BY ordinal_position
        """

        try:
            result = self.execute_query(query, (schema, table))
            return [row[0] for row in result]
        except QueryError:
            return []

    def get_table_columns_cached(self, schema: str, table: str) -> List[str]:
        """Get column names with caching

        Args:
            schema: Schema name
            table: Table name

        Returns:
            List of column names
        """
        cache_key = f"{schema}.{table}"

        if cache_key in self._column_cache:
            return self._column_cache[cache_key]

        # Use existing method
        columns = self.get_column_names(schema, table)

        # Cache the result
        self._column_cache[cache_key] = columns
        self._save_column_cache()

        return columns

    def validate_columns_exist(self, schema: str, table: str, columns: List[str]) -> Dict[str, bool]:
        """Check which columns actually exist in the table

        Args:
            schema: Schema name
            table: Table name
            columns: List of column names to check

        Returns:
            Dict mapping column names to existence status
        """
        available = set(self.get_table_columns_cached(schema, table))
        return {col: col in available for col in columns}

    def get_safe_columns(self, schema: str, table: str, requested_columns: List[str]) -> List[str]:
        """Get only the columns that actually exist

        Args:
            schema: Schema name
            table: Table name
            requested_columns: List of requested column names

        Returns:
            List of columns that actually exist
        """
        validation = self.validate_columns_exist(schema, table, requested_columns)
        safe_columns = [col for col, exists in validation.items() if exists]

        if len(safe_columns) != len(requested_columns):
            missing = [col for col, exists in validation.items() if not exists]
            self.logger.warning(f"Missing columns in {schema}.{table}: {missing}")

        return safe_columns

    def build_safe_select(self, schema: str, table: str, columns: List[str],
                         where_clause: str = None, limit: int = None) -> Optional[str]:
        """Build a SELECT query with validated columns

        Args:
            schema: Schema name
            table: Table name
            columns: List of column names
            where_clause: Optional WHERE clause (without WHERE keyword)
            limit: Optional LIMIT value

        Returns:
            Query string or None if no valid columns
        """
        safe_columns = self.get_safe_columns(schema, table, columns)

        if not safe_columns:
            self.logger.error(f"No valid columns found for {schema}.{table}")
            return None

        query = f"SELECT {', '.join(safe_columns)} FROM {schema}.{table}"

        if where_clause:
            query += f" WHERE {where_clause}"

        if limit:
            query += f" LIMIT {limit}"

        return query

    def refresh_column_cache(self, schema: str = None, table: str = None):
        """Refresh cached column information

        Args:
            schema: Optional schema to refresh (all if None)
            table: Optional specific table to refresh
        """
        if schema and table:
            # Refresh specific table
            cache_key = f"{schema}.{table}"
            columns = self.get_column_names(schema, table)
            self._column_cache[cache_key] = columns
            self.logger.info(f"Refreshed cache for {cache_key}")
        elif schema:
            # Refresh all tables in schema
            try:
                tables_query = """
                SELECT table_name
                FROM information_schema.tables
                WHERE table_schema = %s AND table_type = 'BASE TABLE'
                """
                table_rows = self.execute_query(tables_query, (schema,))

                for (table_name,) in table_rows:
                    cache_key = f"{schema}.{table_name}"
                    columns = self.get_column_names(schema, table_name)
                    self._column_cache[cache_key] = columns

                self.logger.info(f"Refreshed cache for schema {schema}")
            except Exception as e:
                self.logger.error(f"Failed to refresh schema cache for {schema}: {e}")
        else:
            # Refresh common schemas
            for schema_name in ['ecod_schema', 'pdb_analysis']:
                if self.schema_exists(schema_name):
                    self.refresh_column_cache(schema_name)

        self._save_column_cache()

    def validate_query_safety(self, query: str) -> Dict[str, Any]:
        """Basic validation of query safety

        Args:
            query: SQL query to validate

        Returns:
            Dict with validation results
        """
        import re

        results = {
            'safe': True,
            'issues': [],
            'warnings': []
        }

        # Find table references (simple regex)
        table_pattern = r'\b(\w+)\.(\w+)\b'
        table_matches = re.findall(table_pattern, query)

        for schema, table in table_matches:
            if not self.table_exists(schema, table):
                results['safe'] = False
                results['issues'].append(f"Table {schema}.{table} does not exist")

        return results

    def get_primary_key(self, schema: str, table: str) -> Optional[str]:
        """Get primary key column for a table

        Args:
            schema: Schema name
            table: Table name

        Returns:
            Primary key column name or None
        """
        query = """
        SELECT a.attname
        FROM pg_index i
        JOIN pg_attribute a ON a.attrelid = i.indrelid
        AND a.attnum = ANY(i.indkey)
        JOIN pg_namespace n ON n.oid = i.indrelid::regclass::oid
        JOIN pg_class c ON c.oid = i.indrelid
        WHERE i.indisprimary
        AND n.nspname = %s
        AND c.relname = %s
        LIMIT 1
        """

        try:
            result = self.execute_query(query, (schema, table))
            return result[0][0] if result else None
        except QueryError:
            return None

    def export_query_to_csv(self, query: str, output_path: str,
                          params: Optional[Tuple] = None, header: bool = True) -> bool:
        """Export query results to CSV file
        
        Args:
            query: SQL query
            params: Query parameters
            output_path: Output file path
            header: Include header row
            
        Returns:
            True if successful
        """
        import csv
        
        try:
            results = self.execute_query(query, params)
            
            with open(output_path, 'w', newline='') as f:
                writer = csv.writer(f)
                
                # Write header if requested
                if header and results:
                    with self.get_connection() as conn:
                        with conn.cursor() as cursor:
                            cursor.execute(query, params or ())
                            writer.writerow([desc[0] for desc in cursor.description])
                
                # Write data rows
                writer.writerows(results)
                
            self.logger.debug(f"Exported query results to {output_path}")
            return True
        except Exception as e:
            error_msg = f"Error exporting query to CSV: {str(e)}"
            self.logger.error(error_msg)
            return False
