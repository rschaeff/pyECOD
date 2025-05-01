# ecod/repositories/file_repository.py
class FileRepository:
    """Repository for file-related database operations"""
    
    def __init__(self, db_manager):
        self.db = db_manager
    
    def register_file(self, process_id, file_type, file_path, file_exists=True):
        """Register or update file in database"""
        # Implementation...
    
    def get_files_for_process(self, process_id, file_types=None):
        """Get files for a process with optional file type filter"""
        # Implementation...
