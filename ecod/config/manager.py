#!/usr/bin/env python3
"""
Configuration manager for the ECOD pipeline
Handles loading and accessing configuration from various sources.
"""
import os
import yaml
import json
import logging
from typing import Dict, Any, Optional, List, Union

from .schema import ConfigSchema
from .defaults import DEFAULT_CONFIG

class ConfigManager:
    """Configuration manager for the ECOD pipeline"""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize configuration manager
        
        Args:
            config_path: Path to configuration file (optional)
        """
        self.logger = logging.getLogger("ecod.config")
        self.config_path = config_path
        self.config: Dict[str, Any] = {}
        
        # Load configuration
        self._load_defaults()
        
        # Load from config file if provided
        if config_path and os.path.exists(config_path):
            self._load_from_file(config_path)
            
            # Try to load local config if it exists
            local_config_path = self._get_local_config_path(config_path)
            if os.path.exists(local_config_path):
                self._load_from_file(local_config_path)
                self.logger.info(f"Merged local configuration from {local_config_path}")
            
        # Override with environment variables
        self._load_from_env()
        
        # Validate configuration
        self._validate_config()

    def _get_local_config_path(self, config_path: str) -> str:
        """Get path to local configuration file based on main config path"""
        config_dir = os.path.dirname(config_path)
        config_name = os.path.basename(config_path)
        
        # Split filename and extension
        name_parts = os.path.splitext(config_name)
        
        # Format: <filename>.local.<extension>
        local_name = f"{name_parts[0]}.local{name_parts[1]}"
        local_path = os.path.join(config_dir, local_name)
        
        self.logger.debug(f"Looking for local config at: {local_path}")
        if os.path.exists(local_path):
            self.logger.debug(f"Found local config file")
        else:
            self.logger.debug(f"No local config file found at {local_path}")
        
        return local_path

    def _load_defaults(self) -> None:
        """Load default configuration values"""
        self.config = DEFAULT_CONFIG.copy()
        self.logger.debug("Loaded default configuration")
        
    def _load_from_file(self, config_path: str) -> None:
        """Load configuration from YAML or JSON file
        
        Args:
            config_path: Path to configuration file
        """
        try:
            with open(config_path, 'r') as f:
                if config_path.endswith('.json'):
                    file_config = json.load(f)
                else:  # Assume YAML otherwise
                    file_config = yaml.safe_load(f) or {}
                
            # Merge with current config (nested update)
            self._deep_update(self.config, file_config)
            self.logger.info(f"Loaded configuration from {config_path}")
            
        except Exception as e:
            error_msg = f"Error loading config file: {str(e)}"
            self.logger.error(error_msg)
            
    def _load_from_env(self) -> None:
        """Override config with environment variables
        
        Environment variables should be prefixed with ECOD_
        and use double underscore __ for nesting.
        Example: ECOD_DATABASE__HOST for database.host
        """
        prefix = "ECOD_"
        for key, value in os.environ.items():
            if key.startswith(prefix):
                config_key = key[len(prefix):].lower()
                
                # Handle nested keys
                if "__" in config_key:
                    parts = config_key.split("__")
                    self._set_nested_value(self.config, parts, value)
                else:
                    self.config[config_key] = self._convert_value(value)
        
        self.logger.debug("Applied environment variable overrides")
                    
    def _set_nested_value(self, config: Dict[str, Any], 
                         key_parts: List[str], value: str) -> None:
        """Set a nested value in the configuration dictionary
        
        Args:
            config: Configuration dictionary
            key_parts: List of nested key parts
            value: Value to set
        """
        current = config
        for part in key_parts[:-1]:
            if part not in current:
                current[part] = {}
            current = current[part]
        
        current[key_parts[-1]] = self._convert_value(value)
        
    def _convert_value(self, value: str) -> Any:
        """Convert string environment variable to appropriate type
        
        Args:
            value: String value to convert
            
        Returns:
            Converted value with appropriate type
        """
        if value.lower() in ('true', 'yes', '1'):
            return True
        elif value.lower() in ('false', 'no', '0'):
            return False
        try:
            # Try converting to int
            return int(value)
        except ValueError:
            try:
                # Try converting to float
                return float(value)
            except ValueError:
                # Keep as string
                return value
    
    def _deep_update(self, target: Dict[str, Any], source: Dict[str, Any]) -> None:
        """Recursively update target dictionary with values from source
        
        Args:
            target: Target dictionary to update
            source: Source dictionary with new values
        """
        for key, value in source.items():
            if key in target and isinstance(target[key], dict) and isinstance(value, dict):
                self._deep_update(target[key], value)
            else:
                target[key] = value
    
    def _validate_config(self) -> None:
        """Validate the configuration against the schema"""
        errors = ConfigSchema.validate(self.config)
        
        if errors:
            for error in errors:
                self.logger.error(f"Configuration error: {error}")
            self.logger.warning("Using configuration with validation errors")
        else:
            self.logger.debug("Configuration validated successfully")
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get a configuration value using dot notation
        
        Args:
            key: Configuration key (can use dot notation for nested values)
            default: Default value if key not found
            
        Returns:
            Configuration value or default
        """
        if '.' in key:
            parts = key.split('.')
            current = self.config
            for part in parts:
                if part not in current:
                    return default
                current = current[part]
            return current
        return self.config.get(key, default)
    
    def get_db_config(self) -> Dict[str, Any]:
        """Get database configuration as a dictionary
        
        Returns:
            Dictionary with database configuration
        """
        return self.config.get('database', {})
        
    def get_path(self, path_name: str, default: str = "") -> str:
        """Get a path from configuration
        
        Args:
            path_name: Name of the path in configuration
            default: Default value if path not found
            
        Returns:
            Path string
        """
        return self.config.get('paths', {}).get(path_name, default)
        
    def get_tool_path(self, tool_name: str, default: str = "") -> str:
        """Get a tool path from configuration
        
        Args:
            tool_name: Name of the tool in configuration
            default: Default value if tool path not found
            
        Returns:
            Tool path string
        """
        return self.config.get('tools', {}).get(f"{tool_name}_path", default)
    
    def get_reference_path(self, ref_name: str, default: str = "") -> str:
        """Get a reference resource path
        
        Args:
            ref_name: Name of the reference in configuration
            default: Default value if reference not found
            
        Returns:
            Reference path string
        """
        return self.config.get('reference', {}).get(ref_name, default)
    
    def is_enabled(self, feature: str) -> bool:
        """Check if a feature is enabled
        
        Args:
            feature: Feature name to check
            
        Returns:
            True if feature is enabled
        """
        return bool(self.config.get('features', {}).get(feature, False))