# Enhanced ecod/core/db_manager.py
import psycopg2
import logging
from contextlib import contextmanager
from typing import Dict, Any, List, Tuple, Optional, Generator
import psycopg2.extras
import time
import random
from .exceptions import ConnectionError, QueryError

class DBManager:
    """Enhanced Database Manager with improved error handling"""
    
    MAX_RETRIES = 3
    RETRY_BACKOFF = [1, 3, 7]  # Seconds to wait between retries
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.logger = logging.getLogger("ecod.db")
        
        # Validate configuration
        self._validate_config()
        
    def _validate_config(self) -> None:
        """Validate database configuration"""
        required_keys = ['host', 'database', 'user']
        missing_keys = [key for key in required_keys if key not in self.config]
        
        if missing_keys:
            err_msg = f"Missing required database configuration keys: {', '.join(missing_keys)}"
            self.logger.error(err_msg)
            raise ConnectionError(err_msg)
        
    @contextmanager
    def get_connection(self, retry: bool = True) -> Generator[psycopg2.extensions.connection, None, None]:
        """Context manager for database connections with retry logic
        
        Args:
            retry: Whether to retry on connection failures
            
        Yields:
            Database connection
            
        Raises:
            ConnectionError: If connection fails after retries
        """
        conn = None
        retries = self.MAX_RETRIES if retry else 1
        last_exception = None
        
        for attempt in range(retries):
            try:
                self.logger.debug(f"Establishing database connection (attempt {attempt+1}/{retries})")
                conn = psycopg2.connect(**self.config)
                self.logger.debug("Database connection established successfully")
                yield conn
                conn.commit()
                break
            except psycopg2.OperationalError as e:
                last_exception = e
                if attempt < retries - 1:
                    wait_time = self.RETRY_BACKOFF[min(attempt, len(self.RETRY_BACKOFF)-1)]
                    # Add some jitter to avoid thundering herd
                    wait_time = wait_time + random.uniform(0, 0.5)
                    self.logger.warning(f"Database connection failed (attempt {attempt+1}/{retries}): {str(e)}. Retrying in {wait_time:.1f} seconds...")
                    time.sleep(wait_time)
                else:
                    self.logger.error(f"Database connection failed after {retries} attempts: {str(e)}")
                    raise ConnectionError(f"Failed to connect to database after {retries} attempts: {str(e)}") from e
            except Exception as e:
                last_exception = e
                self.logger.error(f"Unexpected database error: {str(e)}")
                if conn:
                    conn.rollback()
                raise QueryError(f"Unexpected database error: {str(e)}") from e
            finally:
                if conn:
                    conn.close()
                    self.logger.debug("Database connection closed")
    
    def execute_query(self, query: str, params: Optional[Tuple] = None, retry: bool = True) -> List[Tuple]:
        """Execute a query and return results with enhanced error handling
        
        Args:
            query: SQL query to execute
            params: Query parameters
            retry: Whether to retry on connection failures
            
        Returns:
            Query results as a list of tuples
            
        Raises:
            QueryError: If query execution fails
        """
        try:
            with self.get_connection(retry) as conn:
                with conn.cursor() as cursor:
                    self.logger.debug(f"Executing query: {query[:100]}{'...' if len(query) > 100 else ''}")
                    cursor.execute(query, params or ())
                    if cursor.description:  # If the query returns rows
                        results = cursor.fetchall()
                        self.logger.debug(f"Query returned {len(results)} rows")
                        return results
                    self.logger.debug("Query executed successfully (no results)")
                    return []
        except ConnectionError as e:
            # Re-raise connection errors as is
            raise
        except Exception as e:
            error_msg = f"Error executing query: {str(e)}"
            self.logger.error(error_msg)
            raise QueryError(error_msg) from e
    
    def execute_dict_query(self, query: str, params: Optional[Tuple] = None, retry: bool = True) -> List[Dict[str, Any]]:
        """Execute a query and return results as dictionaries with enhanced error handling"""
        try:
            with self.get_connection(retry) as conn:
                with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                    self.logger.debug(f"Executing dict query: {query[:100]}{'...' if len(query) > 100 else ''}")
                    cursor.execute(query, params or ())
                    if cursor.description:  # If the query returns rows
                        results = cursor.fetchall()
                        self.logger.debug(f"Dict query returned {len(results)} rows")
                        return results
                    self.logger.debug("Dict query executed successfully (no results)")
                    return []
        except ConnectionError as e:
            # Re-raise connection errors as is
            raise
        except Exception as e:
            error_msg = f"Error executing dict query: {str(e)}"
            self.logger.error(error_msg)
            raise QueryError(error_msg) from e
    
    def insert(self, table: str, data: Dict[str, Any], returning: Optional[str] = None, retry: bool = True) -> Optional[Any]:
        """Insert a record with enhanced error handling and logging"""
        if not data:
            self.logger.warning(f"Attempted to insert empty data into {table}")
            return None
            
        columns = list(data.keys())
        values = [data[col] for col in columns]
        placeholders = [f"%s" for _ in columns]
        
        query = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({', '.join(placeholders)})"
        
        if returning:
            query += f" RETURNING {returning}"
        
        try:
            with self.get_connection(retry) as conn:
                with conn.cursor() as cursor:
                    self.logger.debug(f"Executing insert into {table} ({len(columns)} columns)")
                    cursor.execute(query, values)
                    result = None
                    if returning:
                        result = cursor.fetchone()[0]
                        self.logger.debug(f"Insert returned: {result}")
                    else:
                        self.logger.debug(f"Insert successful (no return value)")
                    return result
        except Exception as e:
            error_msg = f"Error inserting into {table}: {str(e)}"
            self.logger.error(error_msg)
            raise QueryError(error_msg) from e
    
    def update(self, table: str, data: Dict[str, Any], condition: str, 
              condition_params: Tuple, retry: bool = True) -> int:
        """Update records in a table with enhanced error handling"""
        if not data:
            self.logger.warning(f"Attempted to update with empty data in {table}")
            return 0
            
        set_items = [f"{k} = %s" for k in data.keys()]
        values = list(data.values()) + list(condition_params)
        
        query = f"UPDATE {table} SET {', '.join(set_items)} WHERE {condition}"
        
        try:
            with self.get_connection(retry) as conn:
                with conn.cursor() as cursor:
                    self.logger.debug(f"Executing update on {table} ({len(set_items)} columns)")
                    cursor.execute(query, values)
                    affected = cursor.rowcount
                    self.logger.debug(f"Update affected {affected} rows")
                    return affected
        except Exception as e:
            error_msg = f"Error updating {table}: {str(e)}"
            self.logger.error(error_msg)
            raise QueryError(error_msg) from e
    
    def delete(self, table: str, condition: str, condition_params: Tuple, retry: bool = True) -> int:
        """Delete records from a table with enhanced error handling"""
        query = f"DELETE FROM {table} WHERE {condition}"
        
        try:
            with self.get_connection(retry) as conn:
                with conn.cursor() as cursor:
                    self.logger.debug(f"Executing delete from {table}")
                    cursor.execute(query, condition_params)
                    affected = cursor.rowcount
                    self.logger.debug(f"Delete affected {affected} rows")
                    return affected
        except Exception as e:
            error_msg = f"Error deleting from {table}: {str(e)}"
            self.logger.error(error_msg)
            raise QueryError(error_msg) from e
    
    def transaction(self, callback, retry: bool = True):
        """Execute a callback within a transaction with error handling
        
        Args:
            callback: Function that takes a connection and returns a result
            retry: Whether to retry on connection failures
            
        Returns:
            Result of the callback function
        """
        try:
            with self.get_connection(retry) as conn:
                self.logger.debug("Starting database transaction")
                result = callback(conn)
                self.logger.debug("Transaction completed successfully")
                return result
        except Exception as e:
            error_msg = f"Transaction failed: {str(e)}"
            self.logger.error(error_msg)
            raise QueryError(error_msg) from e