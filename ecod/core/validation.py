# ecod/core/validation.py
import re
import logging
from typing import Dict, Any, List, Optional, Union, Tuple
from ecod.exceptions import DataValidationError

logger = logging.getLogger("ecod.validation")

def validate_protein_data(data: Dict[str, Any]) -> Tuple[bool, Optional[str]]:
    """Validate protein data
    
    Args:
        data: Protein data dictionary
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    required_fields = ['pdb_id', 'chain_id', 'sequence']
    missing_fields = [field for field in required_fields if field not in data]
    
    if missing_fields:
        error_msg = f"Missing required fields: {', '.join(missing_fields)}"
        logger.error(error_msg)
        return False, error_msg
    
    # Validate PDB ID format (typically 4 characters)
    if not re.match(r'^[A-Za-z0-9]{4}$', data['pdb_id']):
        error_msg = f"Invalid PDB ID format: {data['pdb_id']}"
        logger.error(error_msg)
        return False, error_msg
    
    # Validate chain ID (typically a single character)
    if not re.match(r'^[A-Za-z0-9]$', data['chain_id']):
        error_msg = f"Invalid chain ID format: {data['chain_id']}"
        logger.error(error_msg)
        return False, error_msg
    
    # Validate sequence (must be valid amino acids)
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    invalid_chars = set(data['sequence'].upper()) - valid_aa
    if invalid_chars:
        error_msg = f"Invalid amino acids in sequence: {', '.join(invalid_chars)}"
        logger.error(error_msg)
        return False, error_msg
    
    # Validate sequence length
    if 'length' in data and data['length'] != len(data['sequence']):
        error_msg = f"Sequence length mismatch: provided {data['length']}, actual {len(data['sequence'])}"
        logger.error(error_msg)
        return False, error_msg
    
    return True, None

def validate_batch_config(config: Dict[str, Any]) -> Tuple[bool, Optional[str]]:
    """Validate batch configuration
    
    Args:
        config: Batch configuration dictionary
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    required_fields = ['type', 'output_dir', 'reference_version']
    missing_fields = [field for field in required_fields if field not in config]
    
    if missing_fields:
        error_msg = f"Missing required batch configuration fields: {', '.join(missing_fields)}"
        logger.error(error_msg)
        return False, error_msg
    
    # Validate batch type
    valid_types = ['blast', 'hhsearch', 'full', 'demo']
    if config['type'] not in valid_types:
        error_msg = f"Invalid batch type: {config['type']}. Must be one of {valid_types}"
        logger.error(error_msg)
        return False, error_msg
    
    # Validate output directory exists or can be created
    try:
        import os
        os.makedirs(config['output_dir'], exist_ok=True)
    except OSError as e:
        error_msg = f"Invalid output directory: {config['output_dir']}. Error: {str(e)}"
        logger.error(error_msg)
        return False, error_msg
    
    return True, None

def validate_file_path(file_path: str, must_exist: bool = False, 
                     file_type: Optional[str] = None) -> Tuple[bool, Optional[str]]:
    """Validate a file path
    
    Args:
        file_path: Path to validate
        must_exist: Whether the file must exist
        file_type: Expected file type/extension
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    import os
    
    # Check if path is valid
    if not file_path:
        error_msg = "Empty file path provided"
        logger.error(error_msg)
        return False, error_msg
    
    # Check if file exists if required
    if must_exist and not os.path.exists(file_path):
        error_msg = f"File does not exist: {file_path}"
        logger.error(error_msg)
        return False, error_msg
    
    # Check file type if specified
    if file_type and not file_path.lower().endswith(f".{file_type.lower()}"):
        error_msg = f"Invalid file type for {file_path}. Expected .{file_type}"
        logger.error(error_msg)
        return False, error_msg
    
    return True, None

def validate_blast_output(file_path: str) -> Tuple[bool, Optional[str]]:
    """Validate BLAST output file
    
    Args:
        file_path: Path to BLAST output file
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    # First check if file exists and has content
    import os
    
    if not os.path.exists(file_path):
        error_msg = f"BLAST output file does not exist: {file_path}"
        logger.error(error_msg)
        return False, error_msg
    
    if os.path.getsize(file_path) == 0:
        error_msg = f"BLAST output file is empty: {file_path}"
        logger.error(error_msg)
        return False, error_msg
    
    # Check if file is valid XML
    try:
        import xml.etree.ElementTree as ET
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Check for BLAST specific elements
        if root.tag != 'BlastOutput':
            error_msg = f"Not a valid BLAST XML file: {file_path}"
            logger.error(error_msg)
            return False, error_msg
        
        return True, None
    except ET.ParseError as e:
        error_msg = f"Invalid XML in BLAST output file {file_path}: {str(e)}"
        logger.error(error_msg)
        return False, error_msg
    except Exception as e:
        error_msg = f"Error validating BLAST output file {file_path}: {str(e)}"
        logger.error(error_msg)
        return False, error_msg

def validate_hhsearch_output(file_path: str) -> Tuple[bool, Optional[str]]:
    """Validate HHSearch output file
    
    Args:
        file_path: Path to HHSearch output file
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    # HHSearch output is in a custom format, not XML
    import os
    
    if not os.path.exists(file_path):
        error_msg = f"HHSearch output file does not exist: {file_path}"
        logger.error(error_msg)
        return False, error_msg
    
    if os.path.getsize(file_path) == 0:
        error_msg = f"HHSearch output file is empty: {file_path}"
        logger.error(error_msg)
        return False, error_msg
    
    # Check for basic HHSearch header signature
    try:
        with open(file_path, 'r') as f:
            header_lines = [f.readline() for _ in range(5)]
            header_text = ''.join(header_lines)
            
            if 'HHsearch' not in header_text and 'Query' not in header_text:
                error_msg = f"Not a valid HHSearch output file: {file_path}"
                logger.error(error_msg)
                return False, error_msg
        
        return True, None
    except Exception as e:
        error_msg = f"Error validating HHSearch output file {file_path}: {str(e)}"
        logger.error(error_msg)
        return False, error_msg