#!/usr/bin/env python3
"""
HHSearch result registration module for ECOD pipeline
Handles HHR file processing, XML conversion, and database registration
"""

import os
import logging
import re
from datetime import datetime
from typing import Dict, Any, List, Optional, Tuple

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.exceptions import PipelineError, FileOperationError
# Import consolidated utility classes
from ecod.utils.hhsearch_utils import HHRParser, HHRToXMLConverter


class HHResultRegistrar:
    """Register HHSearch results and convert to XML for domain analysis"""
    
    def __init__(self, context=None):
        """Initialize with application context"""
        self.context = context or ApplicationContext()
        self.db = self.context.db
        self.config = self.context.config_manager.config
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.hhresult_registrar")
        
        # Initialize using the consolidated utility classes
        self.parser = HHRParser(self.logger)
        self.converter = HHRToXMLConverter(self.logger)
    
    def register_batch_results(self, batch_id: int, force_regenerate: bool = False) -> int:
        """
        Find, convert, and register HHSearch results for a batch
        
        Args:
            batch_id: Batch ID to process
            force_regenerate: Force regeneration of XML files
            
        Returns:
            Number of registered files
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 0
        
        # Get chains with HHR files
        chains = self._find_chains_with_hhr(batch_id)
        if not chains:
            self.logger.warning(f"No chains with HHR files found in batch {batch_id}")
            return 0
        
        # Register files
        registered_count = 0
        for chain in chains:
            if self._register_chain_results(
                chain['process_id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_info['ref_version'],
                batch_info['base_path'],
                force_regenerate
            ):
                registered_count += 1
                
        self.logger.info(f"Registered HHR results for {registered_count} out of {len(chains)} chains")
        return registered_count
    
    def register_specific_chains(self, batch_id: int, chain_ids: List[str], 
                               force_regenerate: bool = False) -> int:
        """
        Find, convert, and register HHSearch results for specific chains
        
        Args:
            batch_id: Batch ID
            chain_ids: List of chain IDs to process (format: "pdbid_chainid")
            force_regenerate: Force regeneration of XML files
            
        Returns:
            Number of registered files
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 0
        
        # Get specific chains
        chains = []
        for chain_id in chain_ids:
            try:
                pdb_id, chain_letter = chain_id.split('_')
                result = self._find_specific_chain(batch_id, pdb_id, chain_letter)
                if result:
                    chains.append(result)
            except ValueError:
                self.logger.warning(f"Invalid chain ID format: {chain_id}, expected pdbid_chainid")
        
        if not chains:
            self.logger.warning(f"No specified chains found in batch {batch_id}")
            return 0
        
        # Register files
        registered_count = 0
        for chain in chains:
            if self._register_chain_results(
                chain['process_id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_info['ref_version'],
                batch_info['base_path'],
                force_regenerate
            ):
                registered_count += 1
                
        self.logger.info(f"Registered HHR results for {registered_count} out of {len(chains)} chains")
        return registered_count
    
    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information"""
        query = """
        SELECT 
            id, 
            batch_name, 
            base_path, 
            ref_version
        FROM 
            ecod_schema.batch
        WHERE 
            id = %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            if rows:
                return rows[0]
        except Exception as e:
            self.logger.error(f"Error retrieving batch information: {str(e)}")
            
        return None
    
    def _find_chains_with_hhr(self, batch_id: int) -> List[Dict[str, Any]]:
        """Find chains with HHR files in a batch"""
        # First attempt to find chains with registered HHR files
        query = """
        SELECT 
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            pf.file_path as hhr_path,
            b.base_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.process_file pf ON ps.id = pf.process_id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            ps.batch_id = %s
            AND pf.file_type = 'hhr'
            AND pf.file_exists = TRUE
        """
        
        chains = self.db.execute_dict_query(query, (batch_id,))
        if chains:
            self.logger.info(f"Found {len(chains)} chains with registered HHR files")
            return chains
        
        # If no registered HHR files, search on filesystem using process entries
        self.logger.info("No registered HHR files found, searching filesystem...")
        
        query = """
        SELECT 
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            b.base_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            ps.batch_id = %s
            AND ps.status IN ('success', 'processing')
        """
        
        result = []
        chains = self.db.execute_dict_query(query, (batch_id,))
        
        if not chains:
            return []
            
        # Check if HHR files exist in filesystem
        ref_version = self._get_batch_ref_version(batch_id)
        for chain in chains:
            # Make sure required values are not None
            if None in [chain['pdb_id'], chain['chain_id'], ref_version, chain['base_path']]:
                self.logger.warning(f"Missing required path component for chain {chain.get('pdb_id', 'unknown')}_{chain.get('chain_id', 'unknown')}")
                continue
                
            # Check standard locations
            potential_hhr_paths = [
                # New standard naming convention for HHsearch results
                os.path.join(chain['base_path'], "hhsearch", 
                             f"{chain['pdb_id']}_{chain['chain_id']}.hhsearch.{ref_version}.hhr"),
                
                # Older naming patterns as fallbacks
                os.path.join(chain['base_path'], "hhsearch", 
                             f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr"),
            ]
            
            # Only add path with relative_path if it's not None
            if chain.get('relative_path'):
                potential_hhr_paths.append(
                    # Old chain-specific directory structure
                    os.path.join(chain['base_path'], chain['relative_path'], f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr")
                )
            
            # Add other potential locations
            potential_hhr_paths.extend([
                os.path.join(chain['base_path'], f"{chain['pdb_id']}_{chain['chain_id']}", f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr"),
                os.path.join(chain['base_path'], "ecod_dump", f"{chain['pdb_id']}_{chain['chain_id']}", f"{chain['pdb_id']}_{chain['chain_id']}.{ref_version}.hhr")
            ])
            
            for path in potential_hhr_paths:
                if os.path.exists(path) and os.path.getsize(path) > 0:
                    chain['hhr_path'] = path
                    result.append(chain)
                    break
        
        self.logger.info(f"Found {len(result)} chains with unregistered HHR files")
        return result
    
    def _find_specific_chain(self, batch_id: int, pdb_id: str, chain_id: str) -> Optional[Dict[str, Any]]:
        """Find a specific chain in a batch"""
        query = """
        SELECT 
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            b.base_path
        FROM 
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        JOIN
            ecod_schema.batch b ON ps.batch_id = b.id
        WHERE 
            ps.batch_id = %s
            AND p.pdb_id = %s
            AND p.chain_id = %s
        LIMIT 1
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id, pdb_id, chain_id))
            if rows:
                chain = rows[0]
                
                # Check for HHR file
                ref_version = self._get_batch_ref_version(batch_id)
                
                # Validate required values
                if None in [pdb_id, chain_id, ref_version, chain['base_path']]:
                    self.logger.error(f"Missing required path component for chain {pdb_id}_{chain_id}")
                    return None
                
                # Handle relative_path being None
                relative_path = chain.get('relative_path', '')
                
                # Check standard locations
                potential_hhr_paths = [
                    # New standard location (flat structure in hhsearch dir)
                    os.path.join(chain['base_path'], "hhsearch", f"{pdb_id}_{chain_id}.{ref_version}.hhr"),
                ]
                
                # Only add path with relative_path if it's not None or empty
                if relative_path:
                    potential_hhr_paths.append(
                        # Old chain-specific directory structure
                        os.path.join(chain['base_path'], relative_path, f"{pdb_id}_{chain_id}.{ref_version}.hhr")
                    )
                
                # Add other potential locations
                potential_hhr_paths.extend([
                    os.path.join(chain['base_path'], f"{pdb_id}_{chain_id}", f"{pdb_id}_{chain_id}.{ref_version}.hhr"),
                    os.path.join(chain['base_path'], "ecod_dump", f"{pdb_id}_{chain_id}", f"{pdb_id}_{chain_id}.{ref_version}.hhr")
                ])
                
                for path in potential_hhr_paths:
                    if os.path.exists(path) and os.path.getsize(path) > 0:
                        chain['hhr_path'] = path
                        return chain
                
                self.logger.warning(f"No HHR file found for {pdb_id}_{chain_id}")
                return None
        except Exception as e:
            self.logger.error(f"Error finding chain {pdb_id}_{chain_id}: {str(e)}")
            
        return None
    
    def _get_batch_ref_version(self, batch_id: int) -> str:
        """Get reference version for a batch"""
        query = "SELECT ref_version FROM ecod_schema.batch WHERE id = %s"
        
        try:
            rows = self.db.execute_query(query, (batch_id,))
            if rows:
                return rows[0][0]
        except Exception as e:
            self.logger.error(f"Error getting reference version: {str(e)}")
            
        # Default to configured reference version
        return self.config.get('reference', {}).get('current_version', 'develop291')
    
    def _register_chain_results(self, process_id: int, pdb_id: str, chain_id: str, 
                              ref_version: str, base_path: str, 
                              force_regenerate: bool = False) -> bool:
        """Register HHSearch results for a chain"""
        try:
            # Define standard paths
            pdb_chain = f"{pdb_id}_{chain_id}"
            hhsearch_dir = os.path.join(base_path, "hhsearch")
            os.makedirs(hhsearch_dir, exist_ok=True)
            
            hhr_file = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhr")
            xml_file = os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml")
            
            # Check if XML file already exists (unless force regenerate)
            if os.path.exists(xml_file) and os.path.getsize(xml_file) > 0 and not force_regenerate:
                self.logger.info(f"XML file already exists for {pdb_chain}, using existing file")
                
                # Register XML file if needed
                self._check_and_register_file(process_id, 'hhsearch_xml', xml_file, base_path)
                
                # Update process status
                self._update_process_status(process_id, "hhsearch_complete")
                
                return True
            
            # Check if HHR file exists
            if not os.path.exists(hhr_file) or os.path.getsize(hhr_file) == 0:
                self.logger.warning(f"HHR file missing or empty: {hhr_file}")
                
                # Try to find HHR file in alternative locations
                found_hhr = self._find_and_move_hhr(pdb_id, chain_id, ref_version, base_path, hhr_file)
                if not found_hhr:
                    return False
            
            # Register HHR file
            self._check_and_register_file(process_id, 'hhr', hhr_file, base_path)
            
            # Convert HHR to XML using the consolidated utility classes
            self.logger.info(f"Converting {hhr_file} to XML...")
            hhr_data = self.parser.parse(hhr_file)
            if not hhr_data:
                self.logger.error(f"Failed to parse HHR file: {hhr_file}")
                return False
            
            # Convert to XML
            xml_string = self.converter.convert(hhr_data, pdb_id, chain_id, ref_version)
            if not xml_string:
                self.logger.error(f"Failed to convert HHR data to XML for {pdb_chain}")
                return False
            
            # Save XML file
            if not self.converter.save(xml_string, xml_file):
                self.logger.error(f"Failed to save XML file: {xml_file}")
                return False
            
            # Register XML file
            self._check_and_register_file(process_id, 'hhsearch_xml', xml_file, base_path)
            
            # Update process status
            self._update_process_status(process_id, "hhsearch_complete")
            
            self.logger.info(f"Successfully registered HHSearch results for {pdb_chain}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error registering results for {pdb_id}_{chain_id}: {str(e)}")
            return False
    
    def _find_and_move_hhr(self, pdb_id: str, chain_id: str, ref_version: str, 
                         base_path: str, target_path: str) -> bool:
        """Find HHR file in alternative locations and move to standard location"""
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Check alternative locations
        alternative_paths = [
            os.path.join(base_path, pdb_chain, f"{pdb_chain}.{ref_version}.hhr"),
            os.path.join(base_path, "ecod_dump", pdb_chain, f"{pdb_chain}.{ref_version}.hhr"),
            os.path.join(base_path, "ecod_dump", f"{pdb_chain}.{ref_version}.hhr")
        ]
        
        for path in alternative_paths:
            if os.path.exists(path) and os.path.getsize(path) > 0:
                self.logger.info(f"Found HHR file at alternative location: {path}")
                
                # Copy file to standard location
                try:
                    os.makedirs(os.path.dirname(target_path), exist_ok=True)
                    
                    # Use shutil.copy2 to preserve metadata
                    import shutil
                    shutil.copy2(path, target_path)
                    
                    self.logger.info(f"Copied HHR file to standard location: {target_path}")
                    return True
                except Exception as e:
                    self.logger.error(f"Error copying HHR file: {str(e)}")
                    return False
        
        self.logger.error(f"Could not find HHR file for {pdb_chain} in any location")
        return False
    
    def _check_and_register_file(self, process_id: int, file_type: str, 
                                file_path: str, base_path: str) -> bool:
        """Check if file is registered in database and register if not"""
        try:
            # Check if file already registered
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """
            
            existing = self.db.execute_query(query, (process_id, file_type))
            
            # Get relative path (for database storage)
            rel_path = os.path.relpath(file_path, base_path)
            
            if existing:
                # Update existing record
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path)
                    },
                    "id = %s",
                    (existing[0][0],)
                )
                self.logger.info(f"Updated {file_type} record for process {process_id}")
            else:
                # Create new record
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": process_id,
                        "file_type": file_type,
                        "file_path": rel_path,
                        "file_exists": True,
                        "file_size": os.path.getsize(file_path)
                    }
                )
                self.logger.info(f"Created new {file_type} record for process {process_id}")
            
            return True
        except Exception as e:
            self.logger.error(f"Error registering file {file_path}: {str(e)}")
            return False
    
    def _update_process_status(self, process_id: int, stage: str, error_message: str = None) -> bool:
        """Update process status in database"""
        try:
            status = "error" if error_message else "success"
            
            update_data = {
                "current_stage": stage,
                "status": status
            }
            
            if error_message:
                update_data["error_message"] = error_message
                
            self.db.update(
                "ecod_schema.process_status",
                update_data,
                "id = %s",
                (process_id,)
            )
            
            return True
        except Exception as e:
            self.logger.error(f"Error updating process status for {process_id}: {str(e)}")
            return False


def register_hhsearch_results(context, batch_id, force_regenerate=False):
    """
    Find HHR files and register them in the database, converting to XML if needed
    
    Args:
        context: Application context
        batch_id: Batch ID to process
        force_regenerate: Force regeneration of XML files
        
    Returns:
        Number of registered files
    """
    registrar = HHResultRegistrar(context)
    return registrar.register_batch_results(batch_id, force_regenerate)