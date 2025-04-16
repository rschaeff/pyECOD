"""
context.py -- Provide application context for pyECOD
"""
import threading
import os
import logging

from typing import Optional
from ecod.db import DBManager
from ecod.config import ConfigManager
from ecod.jobs.factory import create_job_manager

class ApplicationContext:
    """
    Core application context that initializes and provides access to shared resources.
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize application context with shared resources
        
        Args:
            config_path: Optional path to configuration file
        """
        # Initialize configuration
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        
        # Create logger
        self.logger = logging.getLogger("ecod.core")
        
        # Initialize database connection
        self.db = DBManager(self.config_manager.get_db_config())
        
        # Initialize job manager
        self.job_manager = create_job_manager(self.config)
        
        # Add a flag for forcing overwrite of files
        self._force_overwrite = self.config.get('force_overwrite', False)
    
    def is_force_overwrite(self) -> bool:
        """Check if force overwrite is enabled"""
        return self._force_overwrite
    
    def set_force_overwrite(self, value: bool) -> None:
        """Set force overwrite flag"""
        self._force_overwrite = value
    
    @property
    def db(self) -> DBManager:
        """Get database manager
        
        Returns:
            Database manager
        """
        return self.db_manager

    def update_config(self, section: str, key: str, value: Any) -> None:
        """Update configuration value
        
        Args:
            section: Configuration section
            key: Configuration key
            value: New value
        """
        if section not in self.config_manager.config:
            self.config_manager.config[section] = {}
        
        self.config_manager.config[section][key] = value
        self.logger.debug(f"Updated config {section}.{key} = {value}")

    def set_force_overwrite(self, value: bool) -> None:
        """Set the force_overwrite flag in pipeline config
        
        Args:
            value: Flag value to set
        """
        if 'pipeline' not in self.config_manager.config:
            self.config_manager.config['pipeline'] = {}
            
        self.config_manager.config['pipeline']['force_overwrite'] = value
        self.logger.info(f"Force overwrite {'enabled' if value else 'disabled'}")
        
    def is_force_overwrite(self) -> bool:
        """Check if force_overwrite is enabled
        
        Returns:
            True if force_overwrite is enabled
        """
        return self.config_manager.config.get('pipeline', {}).get('force_overwrite', False)