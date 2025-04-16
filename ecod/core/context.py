"""
context.py -- Provide application context for pyECOD
"""
import threading
import os
import logging
from typing import Dict, Any, Optional

from ecod.db.manager import DBManager
from ecod.config import ConfigManager

class ApplicationContext:
    """Application context for the ECOD pipeline"""
    
    _instance = None
    _lock = threading.Lock()
    
    def __new__(cls, config_path: Optional[str] = None):
        """Singleton instance
        
        Args:
            config_path: Path to configuration file
            
        Returns:
            ApplicationContext instance
        """
        with cls._lock:
            if cls._instance is None:
                cls._instance = super(ApplicationContext, cls).__new__(cls)
                cls._instance._initialized = False
            return cls._instance
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize application context
        
        Args:
            config_path: Path to configuration file
        """
        if not hasattr(self, '_initialized') or not self._initialized:
            self.logger = logging.getLogger("ecod.context")
            
            # Initialize configuration
            self.config_manager = ConfigManager(config_path)
            self.logger.info("Configuration initialized")
            
            # Initialize database manager
            self.db_manager = DBManager(self.config_manager.get_db_config())
            self.logger.info("Database manager initialized")
            
            # More initialization...
            
            self._initialized = True
        elif config_path is not None and config_path != self.config_manager.config_path:
            # Re-initialize with new config path
            self.logger.info(f"Re-initializing context with new config: {config_path}")
            self.config_manager = ConfigManager(config_path)
            self.db_manager = DBManager(self.config_manager.get_db_config())
    
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