#!/usr/bin/env python3
"""
Standalone Database Connection Test

Run this script to verify your test database connection works
before running the integration tests.

Usage: python test_db_connection.py
"""

import os
import sys
from pathlib import Path

def test_database_connection():
    """Test the database connection with your credentials"""
    
    # Your test database configuration
    config = {
        'host': 'localhost',
        'port': 5432,
        'database': 'ecod_test',
        'user': 'test_user',
        'password': 'test_pass'
    }
    
    print("=" * 60)
    print("ECOD Integration Test Database Connection Test")
    print("=" * 60)
    
    print("\nTesting database connection with:")
    print(f"  Host: {config['host']}")
    print(f"  Port: {config['port']}")
    print(f"  Database: {config['database']}")
    print(f"  User: {config['user']}")
    print(f"  Password: {'*' * len(config['password'])}")
    
    # Check environment variables
    print("\nEnvironment variables:")
    env_vars = ['PGHOST', 'PGPORT', 'PGDATABASE', 'PGUSER', 'PGPASSWORD', 
                'TEST_DB_HOST', 'TEST_DB_NAME', 'TEST_DB_USER', 'TEST_DB_PASSWORD']
    
    for var in env_vars:
        value = os.getenv(var)
        if value:
            print(f"  {var}: {value if 'PASSWORD' not in var else '*' * len(value)}")
        else:
            print(f"  {var}: NOT SET")
    
    try:
        import psycopg2
        print(f"\n✓ psycopg2 is installed (version: {psycopg2.__version__})")
        
    except ImportError:
        print("\n✗ psycopg2 not installed")
        print("  Install with: pip install psycopg2-binary")
        return False
    
    print("\nTesting connection...")
    
    try:
        # Test connection
        conn = psycopg2.connect(**config)
        print("✓ Connection successful!")
        
        # Test basic query
        with conn.cursor() as cur:
            cur.execute("SELECT version()")
            version = cur.fetchone()[0]
            print(f"✓ PostgreSQL version: {version[:70]}...")
            
            # Test current database and user
            cur.execute("SELECT current_database(), current_user")
            db_name, user = cur.fetchone()
            print(f"✓ Connected to database '{db_name}' as user '{user}'")
            
            # Test table creation permissions
            try:
                cur.execute("CREATE TEMP TABLE test_permissions (id INTEGER)")
                cur.execute("INSERT INTO test_permissions VALUES (1)")
                cur.execute("SELECT COUNT(*) FROM test_permissions")
                count = cur.fetchone()[0]
                cur.execute("DROP TABLE test_permissions")
                print(f"✓ Database permissions: OK (created/inserted/queried/dropped temp table)")
            except Exception as e:
                print(f"⚠ Limited database permissions: {e}")
            
            # Test schema access
            try:
                cur.execute("SELECT schema_name FROM information_schema.schemata WHERE schema_name IN ('ecod_schema', 'pdb_analysis')")
                schemas = [row[0] for row in cur.fetchall()]
                if schemas:
                    print(f"✓ Found schemas: {', '.join(schemas)}")
                else:
                    print("ℹ No ecod_schema or pdb_analysis schemas found (this is OK for a fresh database)")
                    
            except Exception as e:
                print(f"⚠ Could not check schemas: {e}")
        
        conn.close()
        
        print("\n" + "=" * 60)
        print("DATABASE CONNECTION TEST: SUCCESS")
        print("=" * 60)
        print("\nYou can now run integration tests with:")
        print("  export PGHOST=localhost")
        print("  export PGPORT=5432") 
        print("  export PGDATABASE=ecod_test")
        print("  export PGUSER=test_user")
        print("  export PGPASSWORD=test_pass")
        print("  pytest integration/test_domain_partition_e2e.py -v")
        
        return True
        
    except psycopg2.OperationalError as e:
        print(f"\n✗ Connection failed: {e}")
        
        if "does not exist" in str(e):
            print("\nTo fix: Create the database")
            print(f"  createdb -U {config['user']} -h {config['host']} {config['database']}")
            print("  OR")
            print(f"  psql -U {config['user']} -h {config['host']} postgres -c \"CREATE DATABASE {config['database']};\"")
            
        elif "authentication failed" in str(e):
            print("\nTo fix: Check authentication")
            print("  - Verify username and password")
            print("  - Check pg_hba.conf for authentication settings")
            print(f"  - Test with: psql -U {config['user']} -h {config['host']} -d postgres")
            
        elif "could not connect" in str(e):
            print("\nTo fix: Check PostgreSQL server")
            print(f"  - Is PostgreSQL running on {config['host']}:{config['port']}?")
            print("  - Check with: pg_isready -h localhost -p 5432")
            print("  - Or: sudo systemctl status postgresql")
            
        else:
            print(f"\nUnexpected connection error: {e}")
            
        return False
        
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        return False


def test_manual_connection():
    """Test with manual psql command"""
    print("\n" + "-" * 40)
    print("Manual Connection Test")
    print("-" * 40)
    
    cmd = "psql -U test_user -d ecod_test -h localhost -c 'SELECT current_database();'"
    print(f"Testing command: {cmd}")
    print("(You may be prompted for password: test_pass)")
    
    import subprocess
    try:
        result = subprocess.run(
            ['psql', '-U', 'test_user', '-d', 'ecod_test', '-h', 'localhost', '-c', 'SELECT current_database();'],
            capture_output=True,
            text=True,
            input='test_pass\n'  # This might not work depending on setup
        )
        
        if result.returncode == 0:
            print("✓ Manual psql connection successful")
            print(f"Output: {result.stdout.strip()}")
        else:
            print(f"✗ Manual psql connection failed: {result.stderr}")
            
    except FileNotFoundError:
        print("✗ psql command not found")
        print("  Make sure PostgreSQL client tools are installed")
    except Exception as e:
        print(f"✗ Error running psql: {e}")


if __name__ == "__main__":
    print("Starting database connection test...")
    
    success = test_database_connection()
    
    if success:
        test_manual_connection()
    
    print(f"\nTest result: {'PASSED' if success else 'FAILED'}")
    sys.exit(0 if success else 1)
