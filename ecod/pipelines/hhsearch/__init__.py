# ecod/pipelines/hhsearch/__init__.py
"""
HHSearch processing pipeline for the ECOD system.
"""

from .processor import HHRToXMLConverter, HHSearchProcessor
from .service import (
    HHSearchRegistrationService,
    create_service,
    register_batch_results
)
from .models import (
    RegistrationStatus,
    HHSearchFile,
    RegistrationResult,
    BatchRegistrationResult,
    ServiceConfig
)

__all__ = [
    # Processors
    'HHRToXMLConverter',
    'HHSearchProcessor',

    # Service
    'HHSearchRegistrationService',
    'create_service',
    'register_batch_results',

    # Models
    'RegistrationStatus',
    'HHSearchFile',
    'RegistrationResult',
    'BatchRegistrationResult',
    'ServiceConfig'
]
