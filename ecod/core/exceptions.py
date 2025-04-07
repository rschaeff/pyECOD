# ecod/core/exceptions.py
class ECODError(Exception):
    """Base exception for all ECOD pipeline errors"""
    pass

class ConfigurationError(ECODError):
    """Error in configuration settings"""
    pass

class DatabaseError(ECODError):
    """Base class for database-related errors"""
    pass

class ConnectionError(DatabaseError):
    """Error connecting to the database"""
    pass

class QueryError(DatabaseError):
    """Error executing database query"""
    pass

class PipelineError(ECODError):
    """Base class for pipeline process errors"""
    pass

class JobSubmissionError(PipelineError):
    """Error submitting a job to the cluster"""
    pass

class JobExecutionError(PipelineError):
    """Error during job execution"""
    pass

class FileOperationError(ECODError):
    """Error during file operations"""
    pass

class DataValidationError(ECODError):
    """Error validating input data"""
    pass