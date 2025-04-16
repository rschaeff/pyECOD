# ecod/core/base_pipeline.py
import logging
from typing import Optional, Dict, Any

from ecod.core.context import ApplicationContext
from ecod.exceptions import ConfigurationError

class BasePipeline:
    """Base class for all pipeline components"""
    
    def __init__(self, context=None, logger_name=None):
        """
        Initialize with application context
        
        Args:
            context: Application context
            logger_name: Name for the logger
        """
        self.context = context or ApplicationContext()
        self.db = self.context.db
        self.job_manager = self.context.job_manager
        self.config = self.context.config_manager.config
        self.logger = logging.getLogger(logger_name or "ecod.pipeline")
        
        # Initialize configuration
        self._load_configuration()
    
    def _load_configuration(self) -> None:
        """Load component-specific configuration"""
        pass  # To be implemented by subclasses
    
    def _validate_config(self) -> None:
        """Validate component configuration"""
        pass  # To be implemented by subclasses
    
    def register_file(self, process_id: int, file_type: str, file_path: str, file_exists: bool = True) -> None:
        """Register a file with proper duplicate handling"""
        try:
            # Calculate file size if file exists
            file_size = 0
            if file_exists and os.path.exists(file_path):
                file_size = os.path.getsize(file_path)
            
            # Check if record already exists
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """
            
            existing = self.db.execute_query(query, (process_id, file_type))
            
            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": file_path,
                        "file_exists": file_exists,
                        "file_size": file_size
                    },
                    "id = %s",
                    (existing[0][0],)
                )
                self.logger.debug(f"Updated existing {file_type} file record for process {process_id}")
            else:
                # Insert new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": file_path,
                        "file_exists": file_exists,
                        "file_size": file_size
                    }
                )
                self.logger.debug(f"Created new {file_type} file record for process {process_id}")
        except Exception as e:
            self.logger.warning(f"Error registering {file_type} file for process {process_id}: {e}")