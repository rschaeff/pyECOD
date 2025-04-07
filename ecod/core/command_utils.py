# ecod/core/command_utils.py
import subprocess
import logging
import time
import random
from typing import List, Dict, Any, Optional, Tuple, Union

logger = logging.getLogger("ecod.command_utils")

def run_command_with_retry(
    cmd: Union[List[str], str],
    max_retries: int = 3,
    retry_delay: float = 2.0,
    retry_backoff: float = 2.0,
    retry_jitter: float = 0.5,
    timeout: Optional[float] = None,
    check: bool = True,
    **kwargs
) -> subprocess.CompletedProcess:
    """Run a command with retry logic
    
    Args:
        cmd: Command to run (list of strings or string)
        max_retries: Maximum number of retries
        retry_delay: Initial delay between retries in seconds
        retry_backoff: Backoff multiplier for retry delay
        retry_jitter: Random jitter added to retry delay
        timeout: Command timeout in seconds
        check: Whether to check for non-zero return code
        **kwargs: Additional arguments for subprocess.run
        
    Returns:
        CompletedProcess object
        
    Raises:
        subprocess.CalledProcessError: If command fails after all retries
    """
    cmd_str = cmd if isinstance(cmd, str) else ' '.join(cmd)
    logger.debug(f"Running command with retry: {cmd_str}")
    
    error = None
    retry_count = 0
    current_delay = retry_delay
    
    while retry_count <= max_retries:
        try:
            logger.debug(f"Attempt {retry_count + 1}/{max_retries + 1}: {cmd_str}")
            result = subprocess.run(
                cmd,
                text=True,
                capture_output=True,
                timeout=timeout,
                check=check,
                **kwargs
            )
            return result
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            error = e
            retry_count += 1
            
            if retry_count > max_retries:
                logger.error(f"Command failed after {max_retries + 1} attempts: {cmd_str}")
                logger.error(f"Last error: {str(error)}")
                if hasattr(error, 'stdout'):
                    logger.error(f"Stdout: {error.stdout}")
                if hasattr(error, 'stderr'):
                    logger.error(f"Stderr: {error.stderr}")
                raise
            
            # Calculate delay with jitter
            actual_delay = current_delay + random.uniform(0, retry_jitter)
            logger.warning(f"Command failed (attempt {retry_count}/{max_retries + 1}), retrying in {actual_delay:.2f}s: {cmd_str}")
            if hasattr(error, 'stderr'):
                logger.warning(f"Error: {error.stderr}")
                
            time.sleep(actual_delay)
            current_delay *= retry_backoff

def check_command_availability(command: str) -> bool:
    """Check if a command is available in the system path
    
    Args:
        command: Command to check
        
    Returns:
        True if command is available
    """
    try:
        result = subprocess.run(
            ["which", command],
            text=True,
            capture_output=True,
            check=False
        )
        return result.returncode == 0
    except Exception:
        return False

def check_tool_requirements(required_tools: List[str]) -> Tuple[bool, List[str]]:
    """Check if all required tools are available
    
    Args:
        required_tools: List of required tool commands
        
    Returns:
        Tuple of (all_available, missing_tools)
    """
    missing_tools = []
    for tool in required_tools:
        if not check_command_availability(tool):
            missing_tools.append(tool)
    
    return len(missing_tools) == 0, missing_tools