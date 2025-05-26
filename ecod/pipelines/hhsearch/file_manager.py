#!/usr/bin/env python3
"""
File management for HHSearch results
"""
import os
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

from ecod.utils.path_utils import (
    get_standardized_paths, 
    find_files_with_legacy_paths,
    migrate_file_to_standard_path,
    get_file_db_path
)
from .models import HHSearchFile


class HHSearchFileManager:
    """Manages HHSearch result file discovery and migration"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
    
    def find_hhsearch_files(self, base_path: str, pdb_id: str, chain_id: str, 
                           ref_version: str) -> Dict[str, HHSearchFile]:
        """
        Find HHSearch result files for a chain
        
        Returns:
            Dictionary mapping file types to HHSearchFile objects
        """
        results = {}
        
        # Get standard paths
        standard_paths = get_standardized_paths(base_path, pdb_id, chain_id, ref_version)
        
        # Check for files at legacy locations
        legacy_files = find_files_with_legacy_paths(base_path, pdb_id, chain_id, ref_version)
        
        # Check HHR file
        hhr_file = self._find_best_file_location(
            standard_paths['hhr'],
            legacy_files.get('hhr', {}).get('exists_at'),
            additional_patterns=self._get_hhsearch_patterns(base_path, pdb_id, chain_id, ref_version)
        )
        
        if hhr_file:
            results['hhr'] = hhr_file
        
        # Check XML file
        xml_path = standard_paths['hh_xml']
        if os.path.exists(xml_path):
            results['hh_xml'] = self._create_file_info(
                process_id=0,  # Will be set later
                pdb_id=pdb_id,
                chain_id=chain_id,
                file_type='hh_xml',
                file_path=xml_path
            )
        
        return results
    
    def _get_hhsearch_patterns(self, base_path: str, pdb_id: str, 
                              chain_id: str, ref_version: str) -> List[str]:
        """Get additional search patterns specific to HHSearch"""
        patterns = [
            # Common HHSearch output patterns
            os.path.join(base_path, "hhsearch", f"{pdb_id}_{chain_id}.hhsearch.{ref_version}.hhr"),
            os.path.join(base_path, "hhsearch", f"{pdb_id}_{chain_id}.hhr"),
            os.path.join(base_path, f"{pdb_id}_{chain_id}", "hhsearch", f"{pdb_id}_{chain_id}.hhr"),
            
            # Absolute path patterns sometimes seen
            f"/data/ecod/hhsearch/{pdb_id}_{chain_id}.hhsearch.{ref_version}.hhr",
        ]
        
        return [p for p in patterns if p]  # Filter out empty patterns
    
    def _find_best_file_location(self, standard_path: str, 
                                legacy_path: Optional[str],
                                additional_patterns: List[str]) -> Optional[HHSearchFile]:
        """Find the best location for a file"""
        # Check standard path first
        if os.path.exists(standard_path) and os.path.getsize(standard_path) > 0:
            return self._create_file_info(
                process_id=0,
                pdb_id="",  # Will be set by caller
                chain_id="",
                file_type='hhr',
                file_path=standard_path
            )
        
        # Check legacy path
        if legacy_path and os.path.exists(legacy_path) and os.path.getsize(legacy_path) > 0:
            return self._create_file_info(
                process_id=0,
                pdb_id="",
                chain_id="",
                file_type='hhr',
                file_path=legacy_path
            )
        
        # Check additional patterns
        for pattern in additional_patterns:
            if os.path.exists(pattern) and os.path.getsize(pattern) > 0:
                if self.validate_hhsearch_file(pattern):
                    return self._create_file_info(
                        process_id=0,
                        pdb_id="",
                        chain_id="",
                        file_type='hhr',
                        file_path=pattern
                    )
        
        return None
    
    def validate_hhsearch_file(self, file_path: str) -> bool:
        """Validate that a file is HHSearch output (not HHblits)"""
        try:
            with open(file_path, 'r') as f:
                # Read first 20 lines
                header_lines = []
                for _ in range(20):
                    line = f.readline()
                    if not line:
                        break
                    header_lines.append(line)
                
                header_text = ''.join(header_lines)
                
                # Check for HHSearch indicators
                if 'Command' in header_text and 'hhsearch' in header_text.lower():
                    return True
                
                # Check for ECOD domain hits (e.g., "e6igz91")
                import re
                if re.search(r'e\d\w+\d+', header_text):
                    return True
                
                # Check if it's HHblits (which we don't want)
                if 'hhblits' in header_text.lower():
                    self.logger.debug(f"File {file_path} is HHblits output, not HHsearch")
                    return False
                
                return False
                
        except Exception as e:
            self.logger.warning(f"Error validating file {file_path}: {str(e)}")
            return False
    
    def migrate_to_standard_location(self, source_file: HHSearchFile, 
                                    target_path: str) -> bool:
        """Migrate a file to the standard location"""
        try:
            # Ensure target directory exists
            os.makedirs(os.path.dirname(target_path), exist_ok=True)
            
            # Copy file
            shutil.copy2(str(source_file.file_path), target_path)
            
            self.logger.info(f"Migrated {source_file.file_path} to {target_path}")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to migrate file: {str(e)}")
            return False
    
    def _create_file_info(self, process_id: int, pdb_id: str, chain_id: str,
                         file_type: str, file_path: str) -> HHSearchFile:
        """Create HHSearchFile object with file stats"""
        path = Path(file_path)
        
        file_info = HHSearchFile(
            process_id=process_id,
            pdb_id=pdb_id,
            chain_id=chain_id,
            file_type=file_type,
            file_path=path,
            exists=path.exists()
        )
        
        if file_info.exists:
            try:
                stat = path.stat()
                file_info.size = stat.st_size
                file_info.last_modified = datetime.fromtimestamp(stat.st_mtime)
            except Exception:
                pass
        
        return file_info
    
    def get_relative_path(self, base_path: str, file_path: str) -> str:
        """Get relative path for database storage"""
        return get_file_db_path(base_path, file_path)
