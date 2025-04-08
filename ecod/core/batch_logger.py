# ecod/core/batch_logger.py
import os
import json
import logging
import time
from datetime import datetime
from typing import Dict, Any, List, Optional
from ecod.utils.file import atomic_write, ensure_dir

class BatchLogger:
    """Logger for batch operations with structured event logging"""
    
    def __init__(self, batch_id: int, batch_dir: str, log_dir: Optional[str] = None):
        """Initialize batch logger
        
        Args:
            batch_id: Batch ID
            batch_dir: Batch directory
            log_dir: Log directory (defaults to batch_dir/logs)
        """
        self.batch_id = batch_id
        self.batch_dir = batch_dir
        self.log_dir = log_dir or os.path.join(batch_dir, "logs")
        self.logger = logging.getLogger(f"ecod.batch.{batch_id}")
        
        # Create log directory
        ensure_dir(self.log_dir)
        
        # Setup batch log file
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_file = os.path.join(self.log_dir, f"batch_{batch_id}_{timestamp}.log")
        self.events_file = os.path.join(self.log_dir, f"batch_{batch_id}_events.jsonl")
        
        # Add file handler for this batch
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))
        self.logger.addHandler(file_handler)
        
        self.logger.info(f"Initialized batch logger for batch {batch_id}")
        self.logger.info(f"Log file: {self.log_file}")
        self.logger.info(f"Events file: {self.events_file}")
        
        # Log initial batch event
        self.log_event("batch_start", {
            "batch_id": batch_id,
            "start_time": datetime.now().isoformat(),
            "log_file": self.log_file
        })
    
    def log_event(self, event_type: str, data: Dict[str, Any]) -> None:
        """Log a structured event to the events file
        
        Args:
            event_type: Type of event
            data: Event data
        """
        event = {
            "timestamp": datetime.now().isoformat(),
            "batch_id": self.batch_id,
            "event_type": event_type,
            "data": data
        }
        
        # Log to events file
        try:
            with atomic_write(self.events_file, mode='a', encoding='utf-8') as f:
                json.dump(event, f)
                f.write('\n')
                
            # Also log to standard logger
            self.logger.info(f"Event: {event_type} - {json.dumps(data)[:100]}...")
            
        except Exception as e:
            self.logger.error(f"Error logging event {event_type}: {str(e)}")
    
    def log_job_start(self, job_id: str, job_type: str, items: List[Dict[str, Any]]) -> None:
        """Log job start event
        
        Args:
            job_id: Job ID
            job_type: Type of job
            items: Items in the job
        """
        self.log_event("job_start", {
            "job_id": job_id,
            "job_type": job_type,
            "items_count": len(items),
            "items": [item.get('id') for item in items[:10]]  # Log only first 10 item IDs
        })
    
    def log_job_completion(self, job_id: str, status: str, 
                         duration: float, items_processed: int) -> None:
        """Log job completion event
        
        Args:
            job_id: Job ID
            status: Job status
            duration: Job duration in seconds
            items_processed: Number of items processed
        """
        self.log_event("job_completion", {
            "job_id": job_id,
            "status": status,
            "duration_seconds": duration,
            "items_processed": items_processed
        })
    
    def log_item_processing(self, item_id: int, stage: str, status: str, 
                          details: Optional[Dict[str, Any]] = None) -> None:
        """Log item processing event
        
        Args:
            item_id: Item ID
            stage: Processing stage
            status: Processing status
            details: Additional details
        """
        event_data = {
            "item_id": item_id,
            "stage": stage,
            "status": status
        }
        
        if details:
            event_data["details"] = details
            
        self.log_event("item_processing", event_data)
    
    def log_batch_completion(self, status: str, items_processed: int, 
                           items_failed: int, duration: float) -> None:
        """Log batch completion event
        
        Args:
            status: Batch status
            items_processed: Number of items processed
            items_failed: Number of items failed
            duration: Batch duration in seconds
        """
        self.log_event("batch_completion", {
            "status": status,
            "items_processed": items_processed,
            "items_failed": items_failed,
            "duration_seconds": duration,
            "end_time": datetime.now().isoformat()
        })
    
    def log_error(self, error_type: str, message: str, 
                context: Optional[Dict[str, Any]] = None) -> None:
        """Log error event
        
        Args:
            error_type: Type of error
            message: Error message
            context: Error context
        """
        event_data = {
            "error_type": error_type,
            "message": message
        }
        
        if context:
            event_data["context"] = context
            
        self.log_event("error", event_data)
        
        # Also log to standard logger as error
        self.logger.error(f"Error ({error_type}): {message}")