# ecod/core/config.py
import os
import yaml
import logging
from typing import Dict, Any, Optional

class ConfigManager:
    def __init__(self, config_path: Optional[str] = None):
        self.logger = logging.getLogger("ecod.config")
        self.config: Dict[str, Any] = {}
        
        # Load from config file if provided
        if config_path and os.path.exists(config_path):
            self._load_from_file(config_path)
            
        # Override with environment variables
        self._load_from_env()
        
    def _load_from_file(self, config_path: str) -> None:
        """Load configuration from YAML file"""
        try:
            with open(config_path, 'r') as f:
                self.config = yaml.safe_load(f) or {}
            self.logger.info(f"Loaded configuration from {config_path}")
        except Exception as e:
            self.logger.error(f"Error loading config file: {str(e)}")
            
    def _load_from_env(self) -> None:
        """Override config with environment variables"""
        # Database configuration
        db_config = self.config.setdefault('database', {})
        if os.environ.get('ECOD_DB_HOST'):
            db_config['host'] = os.environ.get('ECOD_DB_HOST')
        if os.environ.get('ECOD_DB_PORT'):
            db_config['port'] = int(os.environ.get('ECOD_DB_PORT', '5432'))
        if os.environ.get('ECOD_DB_NAME'):
            db_config['database'] = os.environ.get('ECOD_DB_NAME')
        if os.environ.get('ECOD_DB_USER'):
            db_config['user'] = os.environ.get('ECOD_DB_USER')
        if os.environ.get('ECOD_DB_PASSWORD'):
            db_config['password'] = os.environ.get('ECOD_DB_PASSWORD')
            
        # Path configuration
        paths = self.config.setdefault('paths', {})
        if os.environ.get('ECOD_OUTPUT_DIR'):
            paths['output_dir'] = os.environ.get('ECOD_OUTPUT_DIR')
            
        # Tool configuration
        tools = self.config.setdefault('tools', {})
        if os.environ.get('ECOD_BLAST_PATH'):
            tools['blast_path'] = os.environ.get('ECOD_BLAST_PATH')
        if os.environ.get('ECOD_HHSUITE_PATH'):
            tools['hhsuite_path'] = os.environ.get('ECOD_HHSUITE_PATH')
            
    def get_db_config(self) -> Dict[str, Any]:
        """Get database configuration"""
        return self.config.get('database', {})
        
    def get_path(self, path_name: str, default: str = "") -> str:
        """Get a path from configuration"""
        return self.config.get('paths', {}).get(path_name, default)
        
    def get_tool_path(self, tool_name: str, default: str = "") -> str:
        """Get a tool path from configuration"""
        return self.config.get('tools', {}).get(tool_name, default)