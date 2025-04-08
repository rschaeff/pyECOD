#!/usr/bin/env python3
"""
Exception hierarchy for the ECOD pipeline.
All custom exceptions should inherit from ECODError.
"""
from typing import Dict, Any, Optional


class ECODError(Exception):
    """Base exception for all ECOD-related errors"""
    
    def __init__(self, message: str, details: Optional[Dict[str, Any]] = None):
        """Initialize with error message and optional details
        
        Args:
            message: Error message
            details: Optional details dictionary with context
        """
        self.message = message
        self.details = details or {}
        super().__init__(message)


class ConfigurationError(ECODError):
    """Error related to configuration issues"""
    pass


class DatabaseError(ECODError):
    """Base class for database-related errors"""
    pass


class ConnectionError(DatabaseError):
    """Error connecting to a database"""
    pass


class QueryError(DatabaseError):
    """Error executing a database query"""
    pass


class FileOperationError(ECODError):
    """Error during file operations"""
    pass


class ValidationError(ECODError):
    """Data validation error"""
    pass


class PipelineError(ECODError):
    """Error in pipeline processing"""
    pass


class JobError(ECODError):
    """Base class for job-related errors"""
    pass


class JobSubmissionError(JobError):
    """Error submitting a job to a cluster"""
    pass


class JobExecutionError(JobError):
    """Error during job execution"""
    pass