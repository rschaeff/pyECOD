# ecod/repositories/process_repository.py
class ProcessRepository:
    """Repository for process-related database operations"""
    
    def __init__(self, db_manager):
        self.db = db_manager
    
    def get_processes_for_batch(self, batch_id, filters=None):
        """Get processes for a batch with optional filters"""
        # Implementation...
    
    def update_process_status(self, process_id, stage, status="success", error_message=None):
        """Update process status"""
        # Implementation...
