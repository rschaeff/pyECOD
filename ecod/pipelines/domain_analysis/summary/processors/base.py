#!/usr/bin/env python3
"""
Base classes and utilities for evidence processors.

This module defines the abstract interface that all evidence processors must
implement, along with common utilities and helper classes.
"""

import logging
import xml.etree.ElementTree as ET
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Set, Union
import hashlib
import json

from ecod.models.pipeline import Evidence
from ecod.exceptions import FileOperationError, ValidationError


class ProcessingStatus(Enum):
    """Status of evidence processing"""
    PENDING = auto()
    PROCESSING = auto()
    SUCCESS = auto()
    FAILED = auto()
    SKIPPED = auto()


@dataclass
class ProcessingResult:
    """Result of processing a single file"""
    file_path: Path
    status: ProcessingStatus
    evidence_count: int = 0
    evidence: List[Evidence] = field(default_factory=list)
    error: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def success(self) -> bool:
        """Check if processing was successful"""
        return self.status == ProcessingStatus.SUCCESS

    def add_warning(self, warning: str) -> None:
        """Add a warning message"""
        self.warnings.append(warning)

    def set_error(self, error: str) -> None:
        """Set error and mark as failed"""
        self.error = error
        self.status = ProcessingStatus.FAILED


@dataclass
class ValidationResult:
    """Result of file validation"""
    valid: bool
    file_exists: bool = False
    file_size: int = 0
    format_valid: bool = False
    has_content: bool = False
    error: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


class EvidenceProcessor(ABC):
    """
    Abstract base class for all evidence processors.

    Each processor is responsible for:
    1. Validating input files
    2. Parsing file content
    3. Extracting evidence
    4. Converting to standardized Evidence objects
    """

    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize processor.

        Args:
            logger: Logger instance (creates one if not provided)
        """
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        self._file_cache: Dict[str, Any] = {}
        self._validation_cache: Dict[str, ValidationResult] = {}

    @property
    @abstractmethod
    def supported_formats(self) -> Set[str]:
        """
        Return set of supported file extensions.

        Returns:
            Set of file extensions (e.g., {'.xml', '.hhr'})
        """
        pass

    @property
    @abstractmethod
    def evidence_type(self) -> str:
        """
        Return the type of evidence this processor produces.

        Returns:
            Evidence type string (e.g., 'blast', 'hhsearch')
        """
        pass

    @abstractmethod
    def validate_content(self, content: Any) -> ValidationResult:
        """
        Validate parsed file content.

        Args:
            content: Parsed file content (format depends on processor)

        Returns:
            ValidationResult with detailed validation information
        """
        pass

    @abstractmethod
    def extract_evidence(self, content: Any, file_path: Path) -> List[Evidence]:
        """
        Extract evidence from validated content.

        Args:
            content: Validated file content
            file_path: Path to source file (for metadata)

        Returns:
            List of Evidence objects
        """
        pass

    def validate_file(self, file_path: Path) -> ValidationResult:
        """
        Validate file before processing.

        Args:
            file_path: Path to file to validate

        Returns:
            ValidationResult with validation details
        """
        # Check cache first
        cache_key = str(file_path)
        if cache_key in self._validation_cache:
            return self._validation_cache[cache_key]

        result = ValidationResult(valid=False)

        # Check file exists
        if not file_path.exists():
            result.error = f"File not found: {file_path}"
            self._validation_cache[cache_key] = result
            return result

        result.file_exists = True

        # Check file size
        result.file_size = file_path.stat().st_size
        if result.file_size == 0:
            result.error = "File is empty"
            self._validation_cache[cache_key] = result
            return result

        # Check file extension
        if file_path.suffix.lower() not in self.supported_formats:
            result.error = f"Unsupported format: {file_path.suffix}"
            result.warnings.append(f"Expected one of: {self.supported_formats}")
            self._validation_cache[cache_key] = result
            return result

        # Try to parse file
        try:
            content = self._parse_file(file_path)
            if content is None:
                result.error = "Failed to parse file"
                self._validation_cache[cache_key] = result
                return result

            # Validate content
            content_result = self.validate_content(content)
            result.format_valid = content_result.format_valid
            result.has_content = content_result.has_content
            result.warnings.extend(content_result.warnings)
            result.metadata.update(content_result.metadata)

            if not content_result.valid:
                result.error = content_result.error or "Content validation failed"
            else:
                result.valid = True

        except Exception as e:
            result.error = f"Validation error: {str(e)}"
            self.logger.debug(f"Validation failed for {file_path}", exc_info=True)

        # Cache result
        self._validation_cache[cache_key] = result
        return result

    def process(self, file_path: Path) -> ProcessingResult:
        """
        Process file and return Evidence objects.

        Args:
            file_path: Path to file to process

        Returns:
            ProcessingResult with evidence and status
        """
        result = ProcessingResult(
            file_path=file_path,
            status=ProcessingStatus.PENDING
        )

        # Validate file first
        self.logger.debug(f"Processing {file_path}")
        validation = self.validate_file(file_path)

        if not validation.valid:
            result.set_error(validation.error or "Validation failed")
            result.warnings.extend(validation.warnings)
            return result

        result.status = ProcessingStatus.PROCESSING

        try:
            # Parse file content
            content = self._parse_file(file_path)

            # Extract evidence
            evidence_list = self.extract_evidence(content, file_path)

            # Filter and validate evidence
            valid_evidence = []
            for evidence in evidence_list:
                if self._validate_evidence(evidence):
                    valid_evidence.append(evidence)
                else:
                    result.add_warning(f"Invalid evidence: {evidence.source_id}")

            # Set results
            result.evidence = valid_evidence
            result.evidence_count = len(valid_evidence)
            result.status = ProcessingStatus.SUCCESS

            # Add metadata
            result.metadata['file_size'] = validation.file_size
            result.metadata['evidence_type'] = self.evidence_type

            self.logger.info(
                f"Processed {file_path}: {result.evidence_count} evidence items"
            )

        except Exception as e:
            result.set_error(f"Processing error: {str(e)}")
            self.logger.error(f"Error processing {file_path}", exc_info=True)

        return result

    def process_batch(self, file_paths: List[Path]) -> Dict[Path, ProcessingResult]:
        """
        Process multiple files.

        Args:
            file_paths: List of file paths to process

        Returns:
            Dictionary mapping file paths to results
        """
        results = {}

        for file_path in file_paths:
            results[file_path] = self.process(file_path)

        # Log summary
        success_count = sum(1 for r in results.values() if r.success)
        total_evidence = sum(r.evidence_count for r in results.values())

        self.logger.info(
            f"Batch processing complete: {success_count}/{len(file_paths)} files, "
            f"{total_evidence} total evidence items"
        )

        return results

    def _parse_file(self, file_path: Path) -> Any:
        """
        Parse file content (can be overridden by subclasses).

        Args:
            file_path: Path to file

        Returns:
            Parsed content (format depends on file type)
        """
        # Check cache
        cache_key = str(file_path)
        if cache_key in self._file_cache:
            return self._file_cache[cache_key]

        content = None

        # Default XML parsing
        if file_path.suffix.lower() == '.xml':
            content = self._parse_xml(file_path)
        elif file_path.suffix.lower() == '.json':
            content = self._parse_json(file_path)
        else:
            # Subclasses should override for other formats
            raise NotImplementedError(
                f"Parser not implemented for {file_path.suffix}"
            )

        # Cache parsed content
        if content is not None:
            self._file_cache[cache_key] = content

        return content

    def _parse_xml(self, file_path: Path) -> Optional[ET.Element]:
        """Parse XML file safely."""
        try:
            tree = ET.parse(file_path)
            return tree.getroot()
        except ET.ParseError as e:
            self.logger.error(f"XML parse error in {file_path}: {e}")
            return None
        except Exception as e:
            self.logger.error(f"Error reading {file_path}: {e}")
            return None

    def _parse_json(self, file_path: Path) -> Optional[Dict]:
        """Parse JSON file safely."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            self.logger.error(f"JSON parse error in {file_path}: {e}")
            return None
        except Exception as e:
            self.logger.error(f"Error reading {file_path}: {e}")
            return None

    def _validate_evidence(self, evidence: Evidence) -> bool:
        """
        Validate individual evidence item.

        Args:
            evidence: Evidence to validate

        Returns:
            True if valid
        """
        # Basic validation
        if not evidence.source_id:
            return False

        # Must have some identifying information
        if not any([evidence.domain_id, evidence.query_range]):
            return False

        # Type must match processor
        if evidence.type != self.evidence_type:
            self.logger.warning(
                f"Evidence type mismatch: expected {self.evidence_type}, "
                f"got {evidence.type}"
            )
            return False

        return True

    def clear_cache(self) -> None:
        """Clear all internal caches."""
        self._file_cache.clear()
        self._validation_cache.clear()
        self.logger.debug("Cleared processor caches")

    def get_cache_info(self) -> Dict[str, Any]:
        """Get information about cache usage."""
        return {
            'file_cache_size': len(self._file_cache),
            'validation_cache_size': len(self._validation_cache),
            'cached_files': list(self._file_cache.keys()),
            'memory_usage_estimate': sum(
                len(str(v)) for v in self._file_cache.values()
            )
        }


class XMLProcessor(EvidenceProcessor):
    """Base class for XML-based evidence processors."""

    @property
    def supported_formats(self) -> Set[str]:
        """XML processors support .xml files."""
        return {'.xml'}

    def validate_content(self, content: Any) -> ValidationResult:
        """Validate XML content."""
        result = ValidationResult(valid=False)

        if not isinstance(content, ET.Element):
            result.error = "Content is not an XML Element"
            return result

        result.format_valid = True

        # Check if root tag is expected
        if hasattr(self, 'expected_root_tag'):
            if content.tag != self.expected_root_tag:
                result.error = f"Unexpected root tag: {content.tag}"
                result.warnings.append(f"Expected: {self.expected_root_tag}")
                return result

        # Check for required elements
        if hasattr(self, 'required_elements'):
            for element in self.required_elements:
                if content.find(element) is None:
                    result.warnings.append(f"Missing element: {element}")

        # Check if has meaningful content
        result.has_content = len(content) > 0 or bool(content.text)

        if result.has_content:
            result.valid = True
        else:
            result.error = "No content found in XML"

        return result


class CompositeProcessor(EvidenceProcessor):
    """Processor that combines multiple processors."""

    def __init__(self, processors: List[EvidenceProcessor],
                 logger: Optional[logging.Logger] = None):
        """
        Initialize composite processor.

        Args:
            processors: List of processors to combine
            logger: Logger instance
        """
        super().__init__(logger)
        self.processors = processors

    @property
    def supported_formats(self) -> Set[str]:
        """Combine supported formats from all processors."""
        formats = set()
        for processor in self.processors:
            formats.update(processor.supported_formats)
        return formats

    @property
    def evidence_type(self) -> str:
        """Return composite type."""
        return "composite"

    def validate_content(self, content: Any) -> ValidationResult:
        """Validate using first applicable processor."""
        # This would be implemented based on content type detection
        raise NotImplementedError("Subclass must implement")

    def extract_evidence(self, content: Any, file_path: Path) -> List[Evidence]:
        """Extract evidence using all applicable processors."""
        all_evidence = []

        for processor in self.processors:
            if file_path.suffix.lower() in processor.supported_formats:
                try:
                    result = processor.process(file_path)
                    if result.success:
                        all_evidence.extend(result.evidence)
                except Exception as e:
                    self.logger.error(
                        f"Error in {processor.__class__.__name__}: {e}"
                    )

        return all_evidence


# Utility functions

def create_processor(processor_type: str, **kwargs) -> EvidenceProcessor:
    """
    Factory function to create processors.

    Args:
        processor_type: Type of processor to create
        **kwargs: Additional arguments for processor

    Returns:
        EvidenceProcessor instance

    Raises:
        ValueError: If processor type is unknown
    """
    # This would be implemented with actual processor imports
    processors = {
        # 'blast': BlastEvidenceProcessor,
        # 'hhsearch': HHSearchEvidenceProcessor,
        # etc.
    }

    if processor_type not in processors:
        raise ValueError(f"Unknown processor type: {processor_type}")

    return processors[processor_type](**kwargs)
