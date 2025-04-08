#!/usr/bin/env python3
"""
Validation utilities for the ECOD pipeline
Functions for validating data and formats
"""
import re
import logging
from typing import Dict, Any, List, Optional, Tuple, Set, Union

from ..core.errors import ValidationError

logger = logging.getLogger("ecod.utils.validation")

def validate_pdb_id(pdb_id: str) -> bool:
    """Validate PDB ID format
    
    Args:
        pdb_id: PDB identifier
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not pdb_id:
        raise ValidationError("PDB ID cannot be empty")
    
    # PDB IDs are typically 4 characters
    if not re.match(r'^[A-Za-z0-9]{4}', pdb_id):
        raise ValidationError(f"Invalid PDB ID format: {pdb_id}")
    
    return True

def validate_chain_id(chain_id: str) -> bool:
    """Validate chain ID format
    
    Args:
        chain_id: Chain identifier
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not chain_id:
        raise ValidationError("Chain ID cannot be empty")
    
    # Chain IDs are typically a single character
    if not re.match(r'^[A-Za-z0-9], chain_id):
        raise ValidationError(f"Invalid chain ID format: {chain_id}")
    
    return True

def validate_domain_id(domain_id: str) -> bool:
    """Validate ECOD domain ID format
    
    Args:
        domain_id: Domain identifier (e.g., "e1abcA1")
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not domain_id:
        raise ValidationError("Domain ID cannot be empty")
    
    # ECOD domain IDs follow a specific pattern
    if not re.match(r'^[dge]\d\w{3}\w\d+, domain_id):
        raise ValidationError(f"Invalid domain ID format: {domain_id}")
    
    return True

def validate_range_format(range_str: str) -> bool:
    """Validate range string format
    
    Args:
        range_str: Range string (e.g., "1-100,120-150")
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not range_str:
        raise ValidationError("Range string cannot be empty")
    
    # Split by comma for multi-segment ranges
    for segment in range_str.split(','):
        segment = segment.strip()
        
        # Skip empty segments
        if not segment:
            continue
        
        # Handle single residue case
        if segment.isdigit():
            continue
        
        # Handle range
        if '-' in segment:
            start_str, end_str = segment.split('-')
            
            try:
                start = int(start_str.strip())
                end = int(end_str.strip())
                
                if start > end:
                    raise ValidationError(f"Invalid range: {segment} (start > end)")
            except ValueError:
                raise ValidationError(f"Invalid range format: {segment}")
        else:
            raise ValidationError(f"Invalid range format: {segment}")
    
    return True

def validate_fasta_format(fasta_content: str) -> bool:
    """Validate FASTA format
    
    Args:
        fasta_content: FASTA content
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not fasta_content:
        raise ValidationError("FASTA content cannot be empty")
    
    # Check for header
    lines = fasta_content.strip().split('\n')
    if not lines[0].startswith('>'):
        raise ValidationError("Invalid FASTA format: Missing header")
    
    # Check for sequence
    if len(lines) < 2:
        raise ValidationError("Invalid FASTA format: Missing sequence")
    
    # Check for valid sequence characters
    sequence = ''.join(lines[1:])
    invalid_chars = set(sequence.upper()) - set('ACDEFGHIKLMNPQRSTVWYX-')
    if invalid_chars:
        raise ValidationError(f"Invalid FASTA format: Invalid characters in sequence: {', '.join(invalid_chars)}")
    
    return True

def validate_path_safety(path: str) -> bool:
    """Validate path safety (prevent path traversal)
    
    Args:
        path: File or directory path
        
    Returns:
        True if safe
        
    Raises:
        ValidationError: If path is unsafe
    """
    if not path:
        raise ValidationError("Path cannot be empty")
    
    # Check for path traversal attempts
    if '..' in path or '~' in path:
        raise ValidationError(f"Unsafe path: {path}")
    
    # Check for absolute paths
    if path.startswith('/') or (len(path) > 1 and path[1] == ':'):
        raise ValidationError(f"Absolute paths not allowed: {path}")
    
    return True

def validate_email(email: str) -> bool:
    """Validate email format
    
    Args:
        email: Email address
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not email:
        raise ValidationError("Email cannot be empty")
    
    # Simple email validation
    if not re.match(r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}, email):
        raise ValidationError(f"Invalid email format: {email}")
    
    return True

def validate_url(url: str) -> bool:
    """Validate URL format
    
    Args:
        url: URL
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not url:
        raise ValidationError("URL cannot be empty")
    
    # Simple URL validation
    if not re.match(r'^https?://', url):
        raise ValidationError(f"Invalid URL format: {url}")
    
    return True

def validate_dict_keys(data: Dict[str, Any], required_keys: List[str],
                     optional_keys: Optional[List[str]] = None) -> bool:
    """Validate dictionary keys
    
    Args:
        data: Dictionary to validate
        required_keys: List of required keys
        optional_keys: List of optional keys
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not data:
        raise ValidationError("Data dictionary cannot be empty")
    
    # Check for required keys
    missing_keys = [key for key in required_keys if key not in data]
    if missing_keys:
        raise ValidationError(f"Missing required keys: {', '.join(missing_keys)}")
    
    # Check for unknown keys
    if optional_keys is not None:
        allowed_keys = set(required_keys + optional_keys)
        unknown_keys = [key for key in data if key not in allowed_keys]
        if unknown_keys:
            raise ValidationError(f"Unknown keys: {', '.join(unknown_keys)}")
    
    return True

def validate_date_format(date_str: str, format_str: str = '%Y-%m-%d') -> bool:
    """Validate date format
    
    Args:
        date_str: Date string
        format_str: Expected date format
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if not date_str:
        raise ValidationError("Date string cannot be empty")
    
    try:
        import datetime
        datetime.datetime.strptime(date_str, format_str)
        return True
    except ValueError:
        raise ValidationError(f"Invalid date format: {date_str} (expected {format_str})")

def validate_float_range(value: float, min_value: Optional[float] = None,
                       max_value: Optional[float] = None) -> bool:
    """Validate float value within range
    
    Args:
        value: Float value
        min_value: Minimum allowed value
        max_value: Maximum allowed value
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if min_value is not None and value < min_value:
        raise ValidationError(f"Value {value} is below minimum {min_value}")
    
    if max_value is not None and value > max_value:
        raise ValidationError(f"Value {value} is above maximum {max_value}")
    
    return True

def validate_int_range(value: int, min_value: Optional[int] = None,
                     max_value: Optional[int] = None) -> bool:
    """Validate integer value within range
    
    Args:
        value: Integer value
        min_value: Minimum allowed value
        max_value: Maximum allowed value
        
    Returns:
        True if valid
        
    Raises:
        ValidationError: If validation fails
    """
    if min_value is not None and value < min_value:
        raise ValidationError(f"Value {value} is below minimum {min_value}")
    
    if max_value is not None and value > max_value:
        raise ValidationError(f"Value {value} is above maximum {max_value}")
    
    return True