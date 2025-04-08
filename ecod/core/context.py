#!/usr/bin/env python3
"""
Application context for the ECOD pipeline
Manages shared resources like configuration and database connections
"""
import logging
from typing import Dict, Any, Optional, Set
import threading

from ecod.config import ConfigManager
from ecod.db import DBManager
from ecod.exceptions import ConfigurationError
from ecod.jobs import DatabaseJobManager

class ApplicationContext:
    """Application context for the ECOD pipeline"""
    
    _instance = None
    _lock = threading.Lock()
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize application context
        
        Args:
            config_path: Path to configuration file
        """
        if not self._initialized:
            self.logger = logging.getLogger("ecod.context")
            
            # Initialize configuration
            self.config_manager = ConfigManager(config_path)
            self.logger.info("Configuration initialized")
            
            # Initialize database manager
            self.db_manager = DBManager(self.config_manager.get_db_config())
            self.logger.info("Database manager initialized")

            # Initialize database job manager
            self.job_manager = DatabaseJobManager(self.config_manager.config_path)
            
            # Resource cache
            self._resource_cache: Dict[str, Any] = {}
            
            # Initialize tracking sets
            self._active_batches: Set[int] = set()
            
            self._initialized = True
            self.logger.info("Application context initialized")
        elif config_path is not None and config_path != self.config_manager.config_path:
            # Re-initialize with new config path
            self.logger.info(f"Re-initializing context with new config: {config_path}")
            self.config_manager = ConfigManager(config_path)
            self.db_manager = DBManager(self.config_manager.get_db_config())
    
    @property
    def config(self) -> Dict[str, Any]:
        """Get configuration dictionary
        
        Returns:
            Configuration dictionary
        """
        return self.config_manager.config
    
    @property
    def db(self) -> DBManager:
        """Get database manager
        
        Returns:
            Database manager
        """
        return self.db_manager
    
    def get_resource(self, name: str) -> Any:
        """Get a cached resource
        
        Args:
            name: Resource name
            
        Returns:
            Cached resource or None
        """
        return self._resource_cache.get(name)
    
    def set_resource(self, name: str, resource: Any) -> None:
        """Set a cached resource
        
        Args:
            name: Resource name
            resource: Resource object
        """
        self._resource_cache[name] = resource
    
    def clear_resource(self, name: str) -> None:
        """Clear a cached resource
        
        Args:
            name: Resource name
        """
        if name in self._resource_cache:
            del self._resource_cache[name]
    
    def register_active_batch(self, batch_id: int) -> None:
        """Register an active batch
        
        Args:
            batch_id: Batch ID
        """
        self._active_batches.add(batch_id)
    
    def unregister_active_batch(self, batch_id: int) -> None:
        """Unregister an active batch
        
        Args:
            batch_id: Batch ID
        """
        if batch_id in self._active_batches:
            self._active_batches.remove(batch_id)
    
    def is_batch_active(self, batch_id: int) -> bool:
        """Check if a batch is active
        
        Args:
            batch_id: Batch ID
            
        Returns:
            True if batch is active
        """
        return batch_id in self._active_batches
    
    def get_active_batches(self) -> Set[int]:
        """Get set of active batch IDs
        
        Returns:
            Set of active batch IDs
        """
        return self._active_batches.copy()
    
    def check_database_connection(self) -> bool:
        """Check database connection
        
        Returns:
            True if connection successful
        """
        try:
            # Try a simple query to check connection
            self.db.execute_query("SELECT 1")
            return True
        except Exception as e:
            self.logger.error(f"Database connection check failed: {str(e)}")
            return False
    
    def get_config_value(self, key: str, default: Any = None) -> Any:
        """Get a configuration value
        
        Args:
            key: Configuration key (using dot notation for nested values)
            default: Default value if key not found
            
        Returns:
            Configuration value or default
        """
        return self.config_manager.get(key, default)
    
    def get_path(self, path_name: str, default: str = "") -> str:
        """Get a path from configuration
        
        Args:
            path_name: Path name
            default: Default value if path not found
            
        Returns:
            Path string
        """
        return self.config_manager.get_path(path_name, default)
    
    def get_tool_path(self, tool_name: str, default: str = "") -> str:
        """Get a tool path from configuration
        
        Args:
            tool_name: Tool name
            default: Default value if tool path not found
            
        Returns:
            Tool path string
        """
        return self.config_manager.get_tool_path(tool_name, default)
    
    def get_reference_path(self, ref_name: str, default: str = "") -> str:
        """Get a reference path from configuration
        
        Args:
            ref_name: Reference name
            default: Default value if reference path not found
            
        Returns:
            Reference path string
        """
        return self.config_manager.get_reference_path(ref_name, default)
    
    def is_feature_enabled(self, feature: str) -> bool:
        """Check if a feature is enabled
        
        Args:
            feature: Feature name
            
        Returns:
            True if feature is enabled

            What is this?
        """
        return self.config_manager.is_enabled(feature) __new__(cls, config_path: Optional[str] = None):
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
    