#!/usr/bin/env python3
"""
HHSearch Registration Service

High-level service for registering HHSearch results and converting them to
domain evidence XML format for downstream analysis.
"""
import os
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.exceptions import PipelineError, FileOperationError

from .processor import HHRParser, HHRToXMLConverter
from .file_manager import HHSearchFileManager
from .models import (
    HHSearchFile, RegistrationResult, BatchRegistrationResult,
    RegistrationStatus, ServiceConfig
)


class HHSearchRegistrationService:
    """
    Service for registering and converting HHSearch results.
    
    This service handles:
    1. Finding HHR files in standard and legacy locations
    2. Validating that files are HHSearch (not HHblits) output
    3. Converting HHR to standardized XML format
    4. Registering files in the database
    5. Updating process status
    """
    
    def __init__(self, context: ApplicationContext, 
                 service_config: Optional[Dict[str, Any]] = None):
        """
        Initialize the HHSearch registration service.
        
        Args:
            context: Application context
            service_config: Optional service-specific configuration
        """
        self.context = context
        self.config = context.config_manager.config
        self.db = context.db
        self.logger = logging.getLogger(__name__)
        
        # Initialize service configuration
        config_dict = service_config or {}
        self.service_config = ServiceConfig(**config_dict)
        
        # Initialize components
        self.parser = HHRParser(self.logger)
        self.converter = HHRToXMLConverter(self.logger)
        self.file_manager = HHSearchFileManager(self.logger)
        
        # Service statistics
        self.stats = {
            'files_found': 0,
            'files_converted': 0,
            'files_registered': 0,
            'files_skipped': 0,
            'errors': 0,
            'start_time': datetime.now()
        }
        
        self.logger.info("HHSearchRegistrationService initialized")
    
    def register_batch(self, batch_id: int, limit: Optional[int] = None,
                      chain_ids: Optional[List[str]] = None) -> BatchRegistrationResult:
        """
        Register HHSearch results for a batch.
        
        Args:
            batch_id: Batch ID to process
            limit: Optional limit on number of chains to process
            chain_ids: Optional list of specific chain IDs to process
            
        Returns:
            BatchRegistrationResult with processing statistics
        """
        self.logger.info(f"Starting HHSearch registration for batch {batch_id}")
        
        # Initialize results
        results = BatchRegistrationResult(batch_id=batch_id)
        
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            error_msg = f"Batch {batch_id} not found"
            self.logger.error(error_msg)
            results.finalize()
            return results
        
        # Get chains to process
        if chain_ids:
            chains = self._get_specific_chains(batch_id, chain_ids)
        else:
            chains = self._get_chains_with_hhr(batch_id, limit)
        
        if not chains:
            self.logger.warning(f"No chains with HHR files found in batch {batch_id}")
            results.finalize()
            return results
        
        self.logger.info(f"Processing {len(chains)} chains")
        
        # Process chains
        if self.service_config.parallel_processing and len(chains) > 1:
            self._process_batch_parallel(chains, batch_info, results)
        else:
            self._process_batch_sequential(chains, batch_info, results)
        
        # Finalize results
        results.finalize()
        
        # Update statistics
        self.stats['files_registered'] += results.registered
        self.stats['files_skipped'] += results.skipped
        self.stats['errors'] += results.failed
        
        self.logger.info(
            f"Batch {batch_id} registration complete in {results.processing_time:.2f}s. "
            f"Registered: {results.registered}, Skipped: {results.skipped}, "
            f"Failed: {results.failed} ({results.success_rate:.1f}% success rate)"
        )
        
        return results
    
    def register_chain(self, process_id: int, pdb_id: str, chain_id: str,
                      batch_info: Dict[str, Any]) -> RegistrationResult:
        """
        Register HHSearch results for a single chain.
        
        Args:
            process_id: Process ID
            pdb_id: PDB identifier
            chain_id: Chain identifier
            batch_info: Batch information dictionary
            
        Returns:
            RegistrationResult
        """
        start_time = datetime.now()
        
        result = RegistrationResult(
            process_id=process_id,
            pdb_id=pdb_id,
            chain_id=chain_id,
            status=RegistrationStatus.PENDING
        )
        
        try:
            self.logger.debug(f"Processing {pdb_id}_{chain_id}")
            
            # Find HHSearch files
            files = self.file_manager.find_hhsearch_files(
                batch_info['base_path'],
                pdb_id,
                chain_id,
                batch_info['ref_version']
            )
            
            if 'hhr' not in files:
                result.status = RegistrationStatus.FAILED
                result.error = "HHR file not found"
                return result
            
            hhr_file = files['hhr']
            hhr_file.process_id = process_id
            
            # Check if XML already exists and we're not forcing regeneration
            if 'hh_xml' in files and not self.service_config.force_regenerate:
                self.logger.debug(f"XML file already exists for {pdb_id}_{chain_id}")
                
                # Register existing files
                self._register_file(hhr_file, batch_info['base_path'])
                result.hhr_registered = True
                
                xml_file = files['hh_xml']
                xml_file.process_id = process_id
                self._register_file(xml_file, batch_info['base_path'])
                result.xml_registered = True
                
                # Update process status
                self._update_process_status(process_id, "hhsearch_complete")
                
                result.status = RegistrationStatus.REGISTERED
                return result
            
            # Validate HHR file if configured
            if self.service_config.validate_hhsearch:
                if not self.file_manager.validate_hhsearch_file(str(hhr_file.file_path)):
                    result.status = RegistrationStatus.FAILED
                    result.error = "Invalid HHR file (not HHSearch output)"
                    return result
            
            # Parse HHR file
            self.logger.debug(f"Parsing HHR file: {hhr_file.file_path}")
            hhr_data = self.parser.parse(str(hhr_file.file_path))
            if not hhr_data:
                result.status = RegistrationStatus.FAILED
                result.error = "Failed to parse HHR file"
                return result
            
            self.stats['files_found'] += 1
            
            # Convert to XML
            self.logger.debug(f"Converting to XML with min_probability={self.service_config.min_probability}")
            xml_string = self.converter.convert(
                hhr_data,
                pdb_id,
                chain_id,
                batch_info['ref_version'],
                min_probability=self.service_config.min_probability
            )
            
            if not xml_string:
                result.status = RegistrationStatus.FAILED
                result.error = "Failed to convert HHR to XML"
                return result
            
            # Save XML file
            from ecod.utils.path_utils import get_standardized_paths
            paths = get_standardized_paths(
                batch_info['base_path'],
                pdb_id,
                chain_id,
                batch_info['ref_version']
            )
            
            xml_path = paths['hh_xml']
            if self.converter.save(xml_string, xml_path):
                result.xml_generated = True
                self.stats['files_converted'] += 1
            else:
                result.status = RegistrationStatus.FAILED
                result.error = "Failed to save XML file"
                return result
            
            # Register files in database
            self._register_file(hhr_file, batch_info['base_path'])
            result.hhr_registered = True
            
            xml_file = HHSearchFile(
                process_id=process_id,
                pdb_id=pdb_id,
                chain_id=chain_id,
                file_type='hh_xml',
                file_path=Path(xml_path),
                exists=True,
                size=os.path.getsize(xml_path)
            )
            self._register_file(xml_file, batch_info['base_path'])
            result.xml_registered = True
            
            # Update process status
            self._update_process_status(process_id, "hhsearch_complete")
            
            result.status = RegistrationStatus.REGISTERED
            
        except Exception as e:
            self.logger.error(f"Error processing {pdb_id}_{chain_id}: {str(e)}")
            result.status = RegistrationStatus.FAILED
            result.error = str(e)
            self.stats['errors'] += 1
        
        finally:
            result.processing_time = (datetime.now() - start_time).total_seconds()
        
        return result
    
    def _process_batch_sequential(self, chains: List[Dict[str, Any]], 
                                 batch_info: Dict[str, Any],
                                 results: BatchRegistrationResult) -> None:
        """Process chains sequentially"""
        for i, chain in enumerate(chains):
            result = self.register_chain(
                chain['process_id'],
                chain['pdb_id'],
                chain['chain_id'],
                batch_info
            )
            
            results.add_result(result)
            
            # Log progress
            if (i + 1) % 10 == 0:
                self.logger.info(f"Processed {i + 1}/{len(chains)} chains")
    
    def _process_batch_parallel(self, chains: List[Dict[str, Any]], 
                               batch_info: Dict[str, Any],
                               results: BatchRegistrationResult) -> None:
        """Process chains in parallel"""
        max_workers = min(self.service_config.max_workers, len(chains))
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_chain = {
                executor.submit(
                    self.register_chain,
                    chain['process_id'],
                    chain['pdb_id'],
                    chain['chain_id'],
                    batch_info
                ): chain
                for chain in chains
            }
            
            # Collect results
            completed = 0
            for future in as_completed(future_to_chain):
                chain = future_to_chain[future]
                completed += 1
                
                try:
                    result = future.result(timeout=self.service_config.process_timeout)
                    results.add_result(result)
                except Exception as e:
                    self.logger.error(
                        f"Error processing {chain['pdb_id']}_{chain['chain_id']}: {e}"
                    )
                    
                    # Create error result
                    error_result = RegistrationResult(
                        process_id=chain['process_id'],
                        pdb_id=chain['pdb_id'],
                        chain_id=chain['chain_id'],
                        status=RegistrationStatus.FAILED,
                        error=str(e)
                    )
                    results.add_result(error_result)
                
                # Log progress
                if completed % 10 == 0:
                    self.logger.info(f"Completed {completed}/{len(chains)} chains")
    
    def _get_batch_info(self, batch_id: int) -> Optional[Dict[str, Any]]:
        """Get batch information from database"""
        query = """
        SELECT id, batch_name, base_path, ref_version
        FROM ecod_schema.batch
        WHERE id = %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
            return rows[0] if rows else None
        except Exception as e:
            self.logger.error(f"Error getting batch info: {str(e)}")
            return None
    
    def _get_chains_with_hhr(self, batch_id: int, 
                            limit: Optional[int] = None) -> List[Dict[str, Any]]:
        """Get chains with HHR files that need registration"""
        query = """
        SELECT
            p.id as protein_id,
            p.pdb_id,
            p.chain_id,
            ps.id as process_id,
            ps.relative_path,
            pf_hhr.file_path as hhr_path
        FROM
            ecod_schema.process_status ps
        JOIN
            ecod_schema.protein p ON ps.protein_id = p.id
        LEFT JOIN
            ecod_schema.process_file pf_hhr ON (
                ps.id = pf_hhr.process_id AND 
                pf_hhr.file_type = 'hhr'
            )
        LEFT JOIN
            ecod_schema.process_file pf_xml ON (
                ps.id = pf_xml.process_id AND 
                pf_xml.file_type = 'hh_xml'
            )
        WHERE
            ps.batch_id = %s
            AND ps.current_stage IN ('hhsearch_complete', 'hhsearch_running', 'profile_complete')
            AND (
                -- Has HHR file registered but no XML
                (pf_hhr.id IS NOT NULL AND pf_xml.id IS NULL)
                OR
                -- Or has completed HHSearch but nothing registered yet
                (ps.current_stage = 'hhsearch_complete' AND pf_hhr.id IS NULL)
            )
        ORDER BY p.id
        """
        
        if limit:
            query += f" LIMIT {limit}"
        
        try:
            return list(self.db.execute_dict_query(query, (batch_id,)))
        except Exception as e:
            self.logger.error(f"Error getting chains with HHR: {str(e)}")
            return []
    
    def _get_specific_chains(self, batch_id: int, 
                           chain_ids: List[str]) -> List[Dict[str, Any]]:
        """Get specific chains by ID"""
        results = []
        
        for chain_id in chain_ids:
            try:
                pdb_id, chain_letter = chain_id.split('_')
                
                query = """
                SELECT
                    p.id as protein_id,
                    p.pdb_id,
                    p.chain_id,
                    ps.id as process_id,
                    ps.relative_path
                FROM
                    ecod_schema.process_status ps
                JOIN
                    ecod_schema.protein p ON ps.protein_id = p.id
                WHERE
                    ps.batch_id = %s
                    AND p.pdb_id = %s
                    AND p.chain_id = %s
                LIMIT 1
                """
                
                rows = self.db.execute_dict_query(query, (batch_id, pdb_id, chain_letter))
                if rows:
                    results.append(rows[0])
                    
            except ValueError:
                self.logger.warning(f"Invalid chain ID format: {chain_id}")
        
        return results
    
    def _register_file(self, file: HHSearchFile, base_path: str) -> bool:
        """Register file in database"""
        try:
            # Get relative path
            rel_path = self.file_manager.get_relative_path(base_path, str(file.file_path))
            
            # Check if already registered
            query = """
            SELECT id FROM ecod_schema.process_file
            WHERE process_id = %s AND file_type = %s
            """
            
            existing = self.db.execute_query(query, (file.process_id, file.file_type))
            
            if existing:
                # Update existing
                self.db.update(
                    "ecod_schema.process_file",
                    {
                        "file_path": rel_path,
                        "file_exists": file.exists,
                        "file_size": file.size
                    },
                    "id = %s",
                    (existing[0][0],)
                )
            else:
                # Insert new
                self.db.insert(
                    "ecod_schema.process_file",
                    {
                        "process_id": file.process_id,
                        "file_type": file.file_type,
                        "file_path": rel_path,
                        "file_exists": file.exists,
                        "file_size": file.size
                    }
                )
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error registering file: {str(e)}")
            return False
    
    def _update_process_status(self, process_id: int, stage: str,
                             error_message: Optional[str] = None) -> bool:
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
            self.logger.error(f"Error updating process status: {str(e)}")
            return False
    
    def get_service_statistics(self) -> Dict[str, Any]:
        """Get comprehensive service statistics"""
        runtime = (datetime.now() - self.stats['start_time']).total_seconds()
        
        return {
            'service': {
                **self.stats,
                'runtime_seconds': runtime,
                'files_per_minute': (self.stats['files_converted'] / runtime * 60) if runtime > 0 else 0
            },
            'configuration': {
                'force_regenerate': self.service_config.force_regenerate,
                'min_probability': self.service_config.min_probability,
                'validate_hhsearch': self.service_config.validate_hhsearch,
                'parallel_processing': self.service_config.parallel_processing,
                'max_workers': self.service_config.max_workers
            }
        }
    
    def validate_setup(self) -> Dict[str, bool]:
        """Validate service setup"""
        validations = {}
        
        # Check database
        try:
            self.db.test_connection()
            validations['database'] = True
        except Exception:
            validations['database'] = False
        
        # Check components
        validations['parser_available'] = self.parser is not None
        validations['converter_available'] = self.converter is not None
        validations['file_manager_available'] = self.file_manager is not None
        
        return validations


# Convenience functions

def create_service(config_path: Optional[str] = None,
                  service_config: Optional[Dict[str, Any]] = None) -> HHSearchRegistrationService:
    """
    Create an HHSearch registration service.
    
    Args:
        config_path: Optional configuration file path
        service_config: Optional service-specific configuration
        
    Returns:
        Configured HHSearchRegistrationService
    """
    if config_path:
        context = ApplicationContext(config_path)
    else:
        config_path = os.environ.get('ECOD_CONFIG_PATH', 'config/config.yml')
        context = ApplicationContext(config_path)
    
    return HHSearchRegistrationService(context, service_config)


def register_batch_results(batch_id: int, config_path: Optional[str] = None,
                          **kwargs) -> BatchRegistrationResult:
    """
    Convenience function to register HHSearch results for a batch.
    
    Args:
        batch_id: Batch ID to process
        config_path: Optional configuration file path
        **kwargs: Additional options
        
    Returns:
        BatchRegistrationResult
    """
    service = create_service(config_path)
    return service.register_batch(batch_id, **kwargs)
