# ecod/core/db_manager.py
import psycopg2
import logging
from contextlib import contextmanager
from typing import Dict, Any, List, Tuple, Optional, Generator
import psycopg2.extras

class DBManager:
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.logger = logging.getLogger("ecod.db")
        
    @contextmanager
    def get_connection(self) -> Generator[psycopg2.extensions.connection, None, None]:
        """Context manager for database connections"""
        conn = None
        try:
            conn = psycopg2.connect(**self.config)
            yield conn
            conn.commit()
        except Exception as e:
            if conn:
                conn.rollback()
            self.logger.error(f"Database error: {str(e)}")
            raise
        finally:
            if conn:
                conn.close()
                
    def execute_query(self, query: str, params: Optional[Tuple] = None) -> List[Tuple]:
        """Execute a query and return results"""
        with self.get_connection() as conn:
            with conn.cursor() as cursor:
                cursor.execute(query, params or ())
                if cursor.description:  # If the query returns rows
                    return cursor.fetchall()
                return []
                
    def execute_dict_query(self, query: str, params: Optional[Tuple] = None) -> List[Dict[str, Any]]:
        """Execute a query and return results as dictionaries"""
        with self.get_connection() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                cursor.execute(query, params or ())
                if cursor.description:  # If the query returns rows
                    return cursor.fetchall()
                return []
                
    def insert(self, table: str, data: Dict[str, Any], returning: Optional[str] = None) -> Optional[Any]:
        """Insert a record and optionally return a value"""
        columns = list(data.keys())
        values = [data[col] for col in columns]
        placeholders = [f"%s" for _ in columns]
        
        query = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({', '.join(placeholders)})"
        
        if returning:
            query += f" RETURNING {returning}"
            
        with self.get_connection() as conn:
            with conn.cursor() as cursor:
                cursor.execute(query, values)
                if returning:
                    return cursor.fetchone()[0]
                return None
                
    def update(self, table: str, data: Dict[str, Any], condition: str, 
              condition_params: Tuple) -> int:
        """Update records in a table"""
        set_items = [f"{k} = %s" for k in data.keys()]
        values = list(data.values()) + list(condition_params)
        
        query = f"UPDATE {table} SET {', '.join(set_items)} WHERE {condition}"
        
        with self.get_connection() as conn:
            with conn.cursor() as cursor:
                cursor.execute(query, values)
                return cursor.rowcount