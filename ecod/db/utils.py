def update_protein_processing_state(context, process_id, new_state, error_message=None):
    """Update processing state with proper tracking and history"""
    try:
        # Get current state
        current_state = context.db.execute_query(
            "SELECT current_stage, status FROM ecod_schema.process_status WHERE id = %s", 
            (process_id,)
        )[0]
        
        # Record state transition in history table
        context.db.insert(
            "ecod_schema.process_history", 
            {
                "process_id": process_id,
                "previous_state": current_state[0],
                "new_state": new_state,
                "transition_time": "NOW()",
                "notes": error_message
            }
        )
        
        # Update current state
        context.db.update(
            "ecod_schema.process_status",
            {
                "current_stage": new_state,
                "status": "error" if error_message else "processing",
                "error_message": error_message
            },
            "id = %s",
            (process_id,)
        )
    except Exception as e:
        logger.error(f"Failed to update processing state: {str(e)}")