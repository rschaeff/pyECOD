#!/usr/bin/env python3
"""
ECOD Pipeline Utilities Module
"""
from .file import (
    ensure_dir, safe_open, atomic_write, calculate_md5,
    calculate_md5_string, check_file_exists, get_file_info,
    read_text_file, write_text_file, read_fasta, write_fasta,
    clean_directory
)
from .sequence import (
    validate_sequence, calculate_md5, parse_range, format_range,
    get_positions_from_range, extract_sequence_by_range, calculate_coverage,
    merge_ranges, calculate_sequence_identity, get_sequence_composition,
    translate_alignment_to_sequence_range
)
from .validation import (
    validate_pdb_id, validate_chain_id, validate_domain_id,
    validate_range_format, validate_fasta_format, validate_path_safety,
    validate_email, validate_url, validate_dict_keys,
    validate_date_format, validate_float_range, validate_int_range
)

__all__ = [
    # File utilities
    'ensure_dir', 'safe_open', 'atomic_write', 'calculate_md5',
    'calculate_md5_string', 'check_file_exists', 'get_file_info',
    'read_text_file', 'write_text_file', 'read_fasta', 'write_fasta',
    'clean_directory',
    
    # Sequence utilities
    'validate_sequence', 'calculate_md5', 'parse_range', 'format_range',
    'get_positions_from_range', 'extract_sequence_by_range', 'calculate_coverage',
    'merge_ranges', 'calculate_sequence_identity', 'get_sequence_composition',
    'translate_alignment_to_sequence_range',
    
    # Validation utilities
    'validate_pdb_id', 'validate_chain_id', 'validate_domain_id',
    'validate_range_format', 'validate_fasta_format', 'validate_path_safety',
    'validate_email', 'validate_url', 'validate_dict_keys',
    'validate_date_format', 'validate_float_range', 'validate_int_range'
]