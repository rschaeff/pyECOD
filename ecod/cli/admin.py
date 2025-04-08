"""
    What is this? Target for removal.
"""

def _update_status(args: argparse.Namespace, db: DBManager) -> int:
    """Update status of batches"""
    try:
        # Build query for batches to update
        query = """
        SELECT id, batch_name, total_items
        FROM ecod_schema.batch
        WHERE status != 'completed'
        """
        
        params = []
        if args.batch_id:
            query += " AND id = %s"
            params.append(args.batch_id)
        
        rows = db.execute_dict_query(query, tuple(params) if params else None)
        
        if not rows:
            logger.info("No batches to update")
            return 0
        
        # Update each batch
        for batch in rows:
            batch_id = batch['id']
            
            # Count completed items
            count_query = """
            SELECT COUNT(*)
            FROM ecod_schema.process_status
            WHERE batch_id = %s
            AND status IN ('success', 'completed')
            """
            
            count_result = db.execute_query(count_query, (batch_id,))
            completed_items = count_result[0][0] if count_result else 0
            
            # Determine if batch is completed
            is_completed = completed_items >= batch['total_items']
            status = 'completed' if is_completed else 'processing'
            
            # Update batch status
            update_query = """
            UPDATE ecod_schema.batch
            SET completed_items = %s,
                status = %s,
                completed_at = CASE WHEN %s THEN CURRENT_TIMESTAMP ELSE completed_at END
            WHERE id = %s
            """
            
            db.execute_query(update_query, (completed_items, status, is_completed, batch_id))
            
            logger.info(f"Updated batch {batch_id}: {completed_items}/{batch['total_items']} completed, status: {status}")
        
        return 0
    except Exception as e:
        logger.error(f"Error updating batch status: {str(e)}")
        return 1

def _test_connection(args: argparse.Namespace, db: DBManager) -> int:
    """Test database connection"""
    try:
        # Try a simple query
        result = db.execute_query("SELECT 1")
        if result and result[0][0] == 1:
            logger.info("Database connection successful!")
            
            # Get server version
            version_result = db.execute_query("SELECT version()")
            if version_result:
                logger.info(f"Database server: {version_result[0][0]}")
            
            # Check schemas
            schema_query = """
            SELECT schema_name 
            FROM information_schema.schemata
            WHERE schema_name NOT LIKE 'pg_%' 
              AND schema_name != 'information_schema'
            ORDER BY schema_name
            """
            schema_rows = db.execute_query(schema_query)
            schemas = [row[0] for row in schema_rows if row and len(row) > 0]
            
            print("Database connection successful!")
            print(f"Found {len(schemas)} schemas: {', '.join(schemas)}")
            return 0
        else:
            logger.error("Database connection test failed")
            return 1
    except Exception as e:
        logger.error(f"Database connection error: {str(e)}")
        return 1