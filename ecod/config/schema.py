#!/usr/bin/env python3
"""
Configuration schema definition for validation
"""
from typing import Dict, Any, List, Optional, Type

class ConfigSchema:
    """Configuration schema for validation"""
    
    SCHEMA = {
        'database': {
            'host': {'type': str, 'required': True},
            'port': {'type': int, 'required': True},
            'database': {'type': str, 'required': True},
            'user': {'type': str, 'required': True},
            'password': {'type': str, 'required': False},
        },
        'paths': {
            'output_dir': {'type': str, 'required': True},
            'blast_db': {'type': str, 'required': False},
        },
        'tools': {
            'blast_path': {'type': str, 'required': False},
            'hhblits_path': {'type': str, 'required': False},
            'hhsearch_path': {'type': str, 'required': False},
            'hhmake_path': {'type': str, 'required': False},
        },
        'reference': {
            'uniclust_db': {'type': str, 'required': False},
            'current_version': {'type': str, 'required': True},
            'chain_db': {'type': str, 'required': False},
            'domain_db': {'type': str, 'required': False},
        },
        'logging': {
            'level': {'type': str, 'required': False},
            'format': {'type': str, 'required': False},
            'log_dir': {'type': str, 'required': False},
        },
        'pipeline': {
            'batch_size': {'type': int, 'required': False},
            'threads': {'type': int, 'required': False},
            'memory': {'type': str, 'required': False},
            'timeout': {'type': int, 'required': False},
            'check_interval': {'type': int, 'required': False},
        }
    }
    
    @classmethod
    def validate(cls, config: Dict[str, Any]) -> List[str]:
        """Validate configuration against schema
        
        Args:
            config: Configuration to validate
            
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Check required fields
        for section, fields in cls.SCHEMA.items():
            # Check if required section exists
            if any(props.get('required', False) for _, props in fields.items()):
                if section not in config:
                    errors.append(f"Missing required configuration section: {section}")
                    continue
            
            # Skip if section doesn't exist
            if section not in config:
                continue
                
            # Check required fields in section
            section_config = config[section]
            for field, props in fields.items():
                if props.get('required', False) and field not in section_config:
                    errors.append(f"Missing required configuration field: {section}.{field}")
        
        # Validate field types
        for section, fields in cls.SCHEMA.items():
            if section not in config:
                continue
                
            section_config = config[section]
            for field, props in fields.items():
                if field in section_config and 'type' in props:
                    expected_type = props['type']
                    if not isinstance(section_config[field], expected_type):
                        errors.append(
                            f"Invalid type for {section}.{field}: expected {expected_type.__name__}, "
                            f"got {type(section_config[field]).__name__}"
                        )
        
        return errors