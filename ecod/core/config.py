# ecod/core/config.py
import os
import yaml
import json
import logging
from typing import Dict, Any, Optional, List, Union, Set
import os.path
from ecod.exceptions import ConfigurationError

class ConfigManager:
    """Enhanced configuration management with validation and multiple sources"""
    
    # Configuration schema with defaults and types
    CONFIG_SCHEMA = {
        'database': {
            'host': {'type': str, 'required': True, 'default': 'localhost'},
            'port': {'type': int, 'required': False, 'default': 5432},
            'database': {'type': str, 'required': True},
            'user': {'type': str, 'required': True},
            'password': {'type': str, 'required': False}
        },
        'paths': {
            'output_dir': {'type': str, 'required': True},
            'blast_db': {'type': str, 'required': False}
        },
        'tools': {
            'blast_path': {'type': str, 'required': False, 'default': 'blastp'},
            'hhblits_path': {'type': str, 'required': False, 'default': 'hhblits'},
            'hhsearch_path': {'type': str, 'required': False, 'default': 'hhsearch'},
            'hhmake_path': {'type': str, 'required': False, 'default': 'hhmake'}
        },
        'reference': {
            'uniclust_db': {'type': str, 'required': False},
            'current_version': {'type': str, 'required': True},
            'chain_db': {'type': str, 'required': False},
            'domain_db': {'type': str, 'required': False}
        },
        'logging': {
            'level': {'type': str, 'required': False, 'default': 'INFO'},
            'format': {'type': str, 'required': False, 
                     'default': '%(asctime)s - %(name)s - %(levelname)s - %(message)s'},
            'log_dir': {'type': str, 'required': False}
        },
        'pipeline': {
            'batch_size': {'type': int, 'required': False, 'default': 100},
            'threads': {'type': int, 'required': False, 'default': 8},
            'memory': {'type': str, 'required': False, 'default': '8G'},
            'timeout': {'type': int, 'required': False, 'default': 3600},
            'check_interval': {'type': int, 'required': False, 'default': 300}
        }
    }
    
    def __init__(self, config_path: Optional[str] = None):
        self.logger = logging.getLogger("ecod.config")
        self.config_path = config_path
        self.config: Dict[str, Any] = {}
        
        # Load configuration from multiple sources
        self._load_config()
        
        # Validate the configuration
        self._validate_config()
        
    def _load_config(self) -> None:
        """Load configuration from file and environment variables"""
        # Start with empty config
        self.config = {}
        
        # Load defaults from schema
        self._load_defaults()
        
        # Load from config file if provided
        if self.config_path and os.path.exists(self.config_path):
            self._load_from_file(self.config_path)
            
        # Override with environment variables
        self._load_from_env()
        
        # Create derived configurations
        self._create_derived_configs()
        
        self.logger.debug(f"Loaded configuration with {len(self.config)} top-level keys")
    
    def _load_defaults(self) -> None:
        """Load default values from schema"""
        for section, fields in self.CONFIG_SCHEMA.items():
            self.config.setdefault(section, {})
            for field, props in fields.items():
                if 'default' in props:
                    self.config[section][field] = props['default']
    
    def _load_from_file(self, config_path: str) -> None:
        """Load configuration from YAML or JSON file"""
        try:
            with open(config_path, 'r') as f:
                if config_path.endswith('.json'):
                    file_config = json.load(f)
                else:  # Assume YAML otherwise
                    file_config = yaml.safe_load(f) or {}
                
                # Merge with current config
                self._merge_config(file_config)
                
            self.logger.info(f"Loaded configuration from {config_path}")
        except Exception as e:
            error_msg = f"Error loading config file: {str(e)}"
            self.logger.error(error_msg)
            raise ConfigurationError(error_msg) from e
            
    def _load_from_env(self) -> None:
        """Override config with environment variables"""
        # Define prefix for ECOD environment variables
        prefix = "ECOD_"
        env_vars_used = []
        
        # Process environment variables
        for key, value in os.environ.items():
            if key.startswith(prefix):
                config_key = key[len(prefix):]  # Remove prefix
                
                # Split by double underscore to determine sections
                parts = config_key.split('__')
                
                if len(parts) == 1:
                    # Top-level config
                    self.config[parts[0].lower()] = self._parse_env_value(value)
                elif len(parts) == 2:
                    # Section config
                    section, field = parts
                    section = section.lower()
                    field = field.lower()
                    
                    # Create section if it doesn't exist
                    if section not in self.config:
                        self.config[section] = {}
                    
                    # Set value with appropriate type conversion
                    if section in self.CONFIG_SCHEMA and field in self.CONFIG_SCHEMA[section]:
                        value_type = self.CONFIG_SCHEMA[section][field]['type']
                        self.config[section][field] = self._convert_value(value, value_type)
                    else:
                        self.config[section][field] = self._parse_env_value(value)
                
                env_vars_used.append(key)
        
        if env_vars_used:
            self.logger.info(f"Applied configuration from {len(env_vars_used)} environment variables")
            self.logger.debug(f"Environment variables used: {', '.join(env_vars_used)}")
    
    def _parse_env_value(self, value: str) -> Any:
        """Parse environment variable value with type detection"""
        # Try to parse as JSON first
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            pass
        
        # Try to convert to appropriate Python type
        if value.lower() == 'true':
            return True
        elif value.lower() == 'false':
            return False
        elif value.isdigit():
            return int(value)
        elif self._is_float(value):
            return float(value)
        else:
            return value
    
    def _is_float(self, value: str) -> bool:
        """Check if string can be converted to float"""
        try:
            float(value)
            return True
        except ValueError:
            return False
    
    def _convert_value(self, value: str, value_type: type) -> Any:
        """Convert string value to specified type"""
        try:
            if value_type == bool:
                return value.lower() in ('true', 'yes', '1', 'y')
            return value_type(value)
        except ValueError:
            self.logger.warning(f"Could not convert value '{value}' to {value_type.__name__}")
            return value
    
    def _merge_config(self, new_config: Dict[str, Any]) -> None:
        """Merge new configuration into existing configuration"""
        for section, section_config in new_config.items():
            if isinstance(section_config, dict):
                # Merge section dictionary
                if section not in self.config:
                    self.config[section] = {}
                for key, value in section_config.items():
                    self.config[section][key] = value
            else:
                # Set top-level value
                self.config[section] = section_config
    
    def _create_derived_configs(self) -> None:
        """Create any derived configuration values"""
        # Example: Create full paths from relative paths
        if 'paths' in self.config:
            paths = self.config['paths']
            
            # Make output directory absolute if relative
            if 'output_dir' in paths and not os.path.isabs(paths['output_dir']):
                paths['output_dir'] = os.path.abspath(paths['output_dir'])
    
    def _validate_config(self) -> None:
        """Validate configuration against schema"""
        errors = []
        
        # Check required fields
        for section, fields in self.CONFIG_SCHEMA.items():
            # Check if required section exists
            if any(props.get('required', False) for props in fields.values()):
                if section not in self.config:
                    errors.append(f"Missing required configuration section: {section}")
                    continue
            
            # Skip if section doesn't exist
            if section not in self.config:
                continue
            
            # Check required fields in section
            section_config = self.config[section]
            for field, props in fields.items():
                if props.get('required', False) and field not in section_config:
                    errors.append(f"Missing required configuration field: {section}.{field}")
        
        # Validate field types
        for section, fields in self.CONFIG_SCHEMA.items():
            if section not in self.config:
                continue
                
            section_config = self.config[section]
            for field, props in fields.items():
                if field in section_config and props.get('type') is not None:
                    # Check type compatibility
                    expected_type = props['type']
                    if not isinstance(section_config[field], expected_type):
                        errors.append(
                            f"Invalid type for {section}.{field}: expected {expected_type.__name__}, "
                            f"got {type(section_config[field]).__name__}"
                        )
        
        # Report all validation errors
        if errors:
            error_msg = f"Configuration validation failed:\n- " + "\n- ".join(errors)
            self.logger.error(error_msg)
            raise ConfigurationError(error_msg)
        
        self.logger.debug("Configuration validated successfully")
    
    def get_db_config(self) -> Dict[str, Any]:
        """Get database configuration"""
        return self.config.get('database', {})
        
    def get_path(self, path_name: str, default: str = "") -> str:
        """Get a path from configuration"""
        return self.config.get('paths', {}).get(path_name, default)
        
    def get_tool_path(self, tool_name: str, default: str = "") -> str:
        """Get a tool path from configuration"""
        return self.config.get('tools', {}).get(tool_name, default)
    
    def get_reference_path(self, ref_name: str, default: str = "") -> str:
        """Get a reference path from configuration"""
        return self.config.get('reference', {}).get(ref_name, default)
    
    def get_logging_config(self) -> Dict[str, Any]:
        """Get logging configuration"""
        return self.config.get('logging', {})
    
    def get_pipeline_config(self) -> Dict[str, Any]:
        """Get pipeline configuration"""
        return self.config.get('pipeline', {})
    
    def get(self, path: str, default: Any = None) -> Any:
        """Get a configuration value using dot notation path
        
        Args:
            path: Configuration path (e.g., 'database.host')
            default: Default value if path doesn't exist
            
        Returns:
            Configuration value or default
        """
        parts = path.split('.')
        config = self.config
        
        for part in parts:
            if isinstance(config, dict) and part in config:
                config = config[part]
            else:
                return default
                
        return config
    
    def set(self, path: str, value: Any) -> None:
        """Set a configuration value using dot notation path
        
        Args:
            path: Configuration path (e.g., 'database.host')
            value: Value to set
        """
        parts = path.split('.')
        config = self.config
        
        # Navigate to parent section
        for part in parts[:-1]:
            if part not in config:
                config[part] = {}
            config = config[part]
        
        # Set value
        config[parts[-1]] = value