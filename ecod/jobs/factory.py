#!/usr/bin/env python3
"""
Job manager factory for the ECOD pipeline.
"""
import logging
from typing import Dict, Any, Optional

from ecod.config import ConfigManager
from .base import JobManager
from .local import LocalJobManager
from .slurm import SlurmJobManager

logger = logging.getLogger("ecod.jobs.factory")

def create_job_manager(config: Dict[str, Any], manager_type: Optional[str] = None) -> JobManager:
    """Create a job manager based on configuration
    
    Args:
        config: Configuration dictionary
        manager_type: Optional manager type override ('local' or 'slurm')
        
    Returns:
        JobManager instance
    """
    # Determine job manager type
    if manager_type is None:
        # Get from config
        manager_type = config.get('job_manager', {}).get('type', 'local')
    
    logger.debug(f"Creating job manager of type: {manager_type}")
    
    # Create appropriate manager
    if manager_type.lower() == 'slurm':
        # Try to check if SLURM is available
        try:
            manager = SlurmJobManager(config)
            logger.info("Using SLURM job manager")
            return manager
        except Exception as e:
            logger.warning(f"Failed to create SLURM job manager: {str(e)}")
            logger.warning("Falling back to local job manager")
            return LocalJobManager(config)
    else:
        # Use local job manager
        logger.info("Using local job manager")
        return LocalJobManager(config)