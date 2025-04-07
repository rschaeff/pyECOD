# ecod/core/logging_config.py
import logging
import os
from datetime import datetime
from typing import Optional, Dict, Any

class LoggingManager:
    """Centralized logging configuration for the ECOD pipeline"""
    
    DEFAULT_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    DEFAULT_DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
    
    @staticmethod
    def configure(
        verbose: bool = False, 
        log_file: Optional[str] = None,
        component: str = "ecod",
        log_dir: Optional[str] = None,
        config: Optional[Dict[str, Any]] = None
    ) -> logging.Logger:
        """Configure logging for the application
        
        Args:
            verbose: Enable debug logging if True
            log_file: Specific log file path (overrides automatic naming)
            component: Component name for logger and automatic log file naming
            log_dir: Directory for log files (defaults to logs/ in current directory)
            config: Configuration dictionary that may contain logging settings
            
        Returns:
            Configured logger instance
        """
        # Determine log level
        log_level = logging.DEBUG if verbose else logging.INFO
        
        # Create handlers
        handlers = [logging.StreamHandler()]
        
        # Setup log file if specified or if log_dir is provided
        if log_file or log_dir:
            if not log_file and log_dir:
                # Auto-generate log filename based on component and timestamp
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                if not log_dir:
                    log_dir = "logs"
                os.makedirs(log_dir, exist_ok=True)
                log_file = os.path.join(log_dir, f"{component}_{timestamp}.log")
            
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(logging.Formatter(LoggingManager.DEFAULT_FORMAT))
            handlers.append(file_handler)
        
        # Get format from config if provided
        log_format = config.get('logging', {}).get('format', LoggingManager.DEFAULT_FORMAT) if config else LoggingManager.DEFAULT_FORMAT
        
        # Configure root logger
        logging.basicConfig(
            level=log_level,
            format=log_format,
            handlers=handlers
        )
        
        # Create and return specific logger
        logger = logging.getLogger(component)
        
        # Add contextual information
        logger.info(f"Logging initialized for {component} at level {logging.getLevelName(log_level)}")
        if log_file:
            logger.info(f"Log file: {log_file}")
        
        return logger
    
    @staticmethod
    def get_logger(name: str) -> logging.Logger:
        """Get a named logger with the configured settings
        
        Args:
            name: Logger name, typically module.class format
            
        Returns:
            Logger instance
        """
        return logging.getLogger(name)