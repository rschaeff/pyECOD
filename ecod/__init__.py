#!/usr/bin/env python3
"""
ECOD (Evolutionary Classification of Protein Domains) Pipeline

A Python framework for analyzing and classifying protein domains.
"""

__version__ = '0.1.0'
__author__ = 'ECOD Team'
__email__ = 'example@example.org'
__license__ = 'MIT'

# Import core modules for easier access
from .core.context import ApplicationContext
from .exceptions import ECODError
from .error_handlers import handle_exceptions

# Make key classes available at package level
__all__ = ['ApplicationContext', 'ECODError', 'handle_exceptions']