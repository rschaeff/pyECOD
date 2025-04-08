# ecod/tests/test_db_manager.py
import unittest
import psycopg2
from unittest.mock import patch, MagicMock

from ecod.db import DBManager

class TestDBManager(unittest.TestCase):
    def setUp(self):
        # Create test configuration
        self.test_config = {
            'host': 'localhost',
            'port': 5432,
            'database': 'test_ecod',
            'user': 'test_user',
            'password': 'test_password'
        }
        
        # Create DB manager with mock connection
        self.db_manager = DBManager(self.test_config)
        
    @patch('psycopg2.connect')
    def test_get_connection(self, mock_connect):
        # Setup mock
        mock_conn = MagicMock()
        mock_connect.return_value = mock_conn
        
        # Test connection context manager
        with self.db_manager.get_connection() as conn:
            self.assertEqual(conn, mock_conn)
            
        # Verify connection was committed and closed
        mock_conn.commit.assert_called_once()
        mock_conn.close.assert_called_once()
        
    @patch('psycopg2.connect')
    def test_execute_query(self, mock_connect):
        # Setup mocks
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_conn.cursor.return_value.__enter__.return_value = mock_cursor
        mock_connect.return_value = mock_conn
        
        # Setup return value
        mock_cursor.fetchall.return_value = [(1, 'test'), (2, 'test2')]
        
        # Execute test query
        result = self.db_manager.execute_query("SELECT * FROM test", ('param',))
        
        # Verify query execution
        mock_cursor.execute.assert_called_once_with("SELECT * FROM test", ('param',))
        self.assertEqual(result, [(1, 'test'), (2, 'test2')])
        
    @patch('psycopg2.connect')
    def test_insert(self, mock_connect):
        # Setup mocks
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_conn.cursor.return_value.__enter__.return_value = mock_cursor
        mock_connect.return_value = mock_conn
        
        # Setup return value for RETURNING clause
        mock_cursor.fetchone.return_value = (123,)
        
        # Test insert with returning
        result = self.db_manager.insert(
            "ecod_schema.test_table",
            {"column1": "value1", "column2": 42},
            "id"
        )
        
        # Verify correct SQL generation and execution
        self.assertEqual(result, 123)
        mock_cursor.execute.assert_called_once()
        call_args = mock_cursor.execute.call_args[0]
        self.assertTrue("INSERT INTO ecod_schema.test_table" in call_args[0])
        self.assertTrue("RETURNING id" in call_args[0])