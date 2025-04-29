# ecod/utils/path_utils.py
"""
Standardized file path utilities for pyECOD

This module provides consistent path construction for all file types
used in the pyECOD processing pipeline, ensuring standardization across
different components.
"""

import os
import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

def get_standardized_paths(batch_path: str, pdb_id: str, chain_id: str, 
                          ref_version: str, create_dirs: bool = True) -> Dict[str, str]:
    """Get standardized file paths for a protein chain
    
    Args:
        batch_path: Path to batch directory
        pdb_id: PDB ID
        chain_id: Chain ID
        ref_version: Reference version
        create_dirs: Whether to create directories if they don't exist
        
    Returns:
        Dictionary of standardized paths
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    # Define all subdirectories following the established structure
    dirs = {
        'fastas': os.path.join(batch_path, "fastas"),
        'hhsearch': os.path.join(batch_path, "hhsearch"),
        'hhsearch_profiles': os.path.join(batch_path, "hhsearch", "profiles"),
        'blast_chain': os.path.join(batch_path, "blast", "chain"),
        'blast_domain': os.path.join(batch_path, "blast", "domain"),
        'domains': os.path.join(batch_path, "domains"),
    }
    
    # Create directories if requested
    if create_dirs:
        for directory in dirs.values():
            try:
                os.makedirs(directory, exist_ok=True)
                logger.debug(f"Ensured directory exists: {directory}")
            except Exception as e:
                logger.warning(f"Failed to create directory {directory}: {str(e)}")
    
    # Define paths dictionary with standardized naming
    paths = {
        # FASTA files
        'fasta': os.path.join(dirs['fastas'], f"{pdb_chain}.fa"),
        
        # HHblits profile files
        'a3m': os.path.join(dirs['hhsearch_profiles'], f"{pdb_chain}.a3m"),
        'hhm': os.path.join(dirs['hhsearch_profiles'], f"{pdb_chain}.hhm"),
        
        # HHSearch result files
        'hhr': os.path.join(dirs['hhsearch'], f"{pdb_chain}.{ref_version}.hhr"),
        'hh_xml': os.path.join(dirs['hhsearch'], f"{pdb_chain}.{ref_version}.xml"),
        
        # BLAST result files
        'chain_blast': os.path.join(dirs['blast_chain'], f"{pdb_chain}.{ref_version}.xml"),
        'domain_blast': os.path.join(dirs['blast_domain'], f"{pdb_chain}.{ref_version}.xml"),
        
        # Domain files - clear distinction between summary and partition
        'domain_summary': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.domain_summary.xml"),
        'domain_partition': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.domains.xml"),

        # Blast-only variants
        'blast_only_summary': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.blast_only.domain_summary.xml"),
        'blast_only_partition': os.path.join(dirs['domains'], f"{pdb_chain}.{ref_version}.blast_only.domains.xml"),
    }

    return paths

def get_file_type_from_path(file_path: str) -> str:
    """Determine file type from file path

    Args:
        file_path: Path to file

    Returns:
        File type identifier or None if unknown
    """
    basename = os.path.basename(file_path).lower()

    # Determine type by extension and naming pattern
    if basename.endswith('.fa'):
        return 'fasta'
    elif basename.endswith('.a3m'):
        return 'a3m'
    elif basename.endswith('.hhm'):
        return 'hhm'
    elif basename.endswith('.hhr'):
        return 'hhr'
    elif basename.endswith('.domains.xml'):
        if 'blast_only' in basename:
            return 'blast_only_partition'
        else:
            return 'domain_partition'
    elif basename.endswith('.domain_summary.xml'):
        if 'blast_only' in basename:
            return 'blast_only_summary'
        else:
            return 'domain_summary'
    elif basename.endswith('.hhsearch.xml') or basename.endswith('.hh.xml'):
        return 'hh_xml'
    elif basename.endswith('.xml') and '/blast/chain/' in file_path.replace('\\', '/'):
        return 'chain_blast'
    elif basename.endswith('.xml') and '/blast/domain/' in file_path.replace('\\', '/'):
        return 'domain_blast'
    elif basename.endswith('.chainwise_blast.xml'):
        return 'chain_blast'
    elif basename.endswith('.blast.xml'):
        return 'domain_blast'

    # Default case - unknown file type
    return None

def find_files_with_legacy_paths(batch_path: str, pdb_id: str, chain_id: str,
                              ref_version: str) -> Dict[str, Dict[str, str]]:
    """Find files using both standard and legacy path patterns

    Args:
        batch_path: Path to batch directory
        pdb_id: PDB ID
        chain_id: Chain ID
        ref_version: Reference version

    Returns:
        Dictionary with file types and paths for both standard and legacy formats
    """
    pdb_chain = f"{pdb_id}_{chain_id}"
    results = {}

    # Define legacy patterns with clear distinction between domain_summary and domain_partition
    legacy_patterns = {
        'hhr': [
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.hhsearch.{ref_version}.hhr"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}_{ref_version}.hhr"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.{ref_version}.hhr")
        ],
        'hh_xml': [
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.hhsearch.{ref_version}.xml"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}.{ref_version}.hhsearch.xml"),
            os.path.join(batch_path, "hhsearch", f"{pdb_chain}_{ref_version}.xml")
        ],
        'chain_blast': [
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.chainwise_blast.xml"),
            os.path.join(batch_path, "blast", "chain", f"{pdb_chain}.chainwise.{ref_version}.xml"),
        ],
        'domain_blast': [
            os.path.join(batch_path, "blast", f"{pdb_chain}.{ref_version}.blast.xml"),
            os.path.join(batch_path, "blast", "domain", f"{pdb_chain}.{ref_version}.blast.xml"),
        ],
        'domain_summary': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.domain_summary.xml"),
        ],
        'domain_partition': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.domains.xml"),
        ],
        'blast_only_summary': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.blast_only.domain_summary.xml"),
        ],
        'blast_only_partition': [
            os.path.join(batch_path, "domains", f"{pdb_chain}.{ref_version}.blast_only.domains.xml"),
        ]
    }

    # Get standard paths
    standard_paths = get_standardized_paths(batch_path, pdb_id, chain_id, ref_version, create_dirs=False)

    # Check both standard and legacy paths
    for file_type, standard_path in standard_paths.items():
        results[file_type] = {
            'standard': standard_path,
            'legacy': None,
            'exists_at': None
        }

        # Check standard path
        if os.path.exists(standard_path):
            results[file_type]['exists_at'] = standard_path
            continue

        # Check legacy paths if file type has legacy patterns
        if file_type in legacy_patterns:
            for legacy_path in legacy_patterns[file_type]:
                if os.path.exists(legacy_path):
                    results[file_type]['legacy'] = legacy_path
                    results[file_type]['exists_at'] = legacy_path
                    break

    return results
