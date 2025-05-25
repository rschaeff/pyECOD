#!/usr/bin/env python3
"""
File location and path resolution for evidence files.

This module handles finding and validating evidence files across different
locations, with database integration and fallback strategies.
"""

import os
import re
import logging
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Set, Union
from datetime import datetime, timedelta

from ecod.db import DBManager
from ecod.exceptions import FileOperationError


class FileType(Enum):
    """Evidence file types"""
    FASTA = "fasta"
    CHAIN_BLAST = "chain_blast_result"
    DOMAIN_BLAST = "domain_blast_result"
    HHSEARCH = "hhsearch_result"
    SELF_COMPARISON = "self_comparison"
    DOMAIN_SUMMARY = "domain_summary"
    DOMAIN_PARTITION = "domain_partition"
    HHR = "hhr"
    HHSEARCH_XML = "hhsearch_xml"
    PROFILE = "hhm_profile"


@dataclass
class FileLocation:
    """Information about a located file"""
    path: Path
    exists: bool
    size: int = 0
    modified_time: Optional[datetime] = None
    source: str = "unknown"  # 'database', 'standard', 'alternative'
    file_type: Optional[FileType] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_valid(self) -> bool:
        """Check if file is valid (exists and has content)"""
        return self.exists and self.size > 0

    @property
    def age_days(self) -> Optional[float]:
        """Get file age in days"""
        if self.modified_time:
            return (datetime.now() - self.modified_time).total_seconds() / 86400
        return None


@dataclass
class SearchResult:
    """Result of file search operation"""
    requested_type: FileType
    protein_id: str
    locations: List[FileLocation] = field(default_factory=list)
    search_paths: List[Path] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    @property
    def found(self) -> bool:
        """Check if any valid files were found"""
        return any(loc.is_valid for loc in self.locations)

    @property
    def best_location(self) -> Optional[FileLocation]:
        """Get best (most recent, largest) valid location"""
        valid_locations = [loc for loc in self.locations if loc.is_valid]
        if not valid_locations:
            return None

        # Prefer database source, then standard, then alternative
        for source in ['database', 'standard', 'alternative']:
            source_locs = [loc for loc in valid_locations if loc.source == source]
            if source_locs:
                # Within source, prefer most recent
                return max(source_locs, key=lambda x: x.modified_time or datetime.min)

        return valid_locations[0]


class EvidenceFileLocator:
    """Handles file location and path resolution for evidence files"""

    def __init__(self, db_manager: DBManager, job_dump_dir: Union[str, Path],
                 reference: str = "develop291", cache_ttl: int = 3600):
        """
        Initialize file locator.

        Args:
            db_manager: Database manager instance
            job_dump_dir: Base directory for job files
            reference: Reference version for naming
            cache_ttl: Cache time-to-live in seconds
        """
        self.db = db_manager
        self.job_dump_dir = Path(job_dump_dir)
        self.reference = reference
        self.cache_ttl = cache_ttl
        self.logger = logging.getLogger(__name__)

        # Initialize cache
        self._cache: Dict[str, Tuple[SearchResult, datetime]] = {}

        # Define standard directory structure
        self.dir_structure = {
            FileType.FASTA: ["sequences", "fasta", "."],
            FileType.CHAIN_BLAST: ["blast", "chain_blast", "."],
            FileType.DOMAIN_BLAST: ["blast", "domain_blast", "."],
            FileType.HHSEARCH: ["hhsearch", "."],
            FileType.SELF_COMPARISON: ["self_comparisons", "."],
            FileType.DOMAIN_SUMMARY: ["domains", "."],
            FileType.DOMAIN_PARTITION: ["domains", "."],
            FileType.PROFILE: ["profiles", "hhm", "."]
        }

        # Define naming patterns
        self.naming_patterns = self._initialize_naming_patterns()

    def _initialize_naming_patterns(self) -> Dict[FileType, List[str]]:
        """Initialize file naming patterns"""
        return {
            FileType.FASTA: [
                "{pdb_id}_{chain_id}.fasta",
                "{pdb_id}_{chain_id}.fa",
                "{pdb_id}{chain_id}.fasta"
            ],
            FileType.CHAIN_BLAST: [
                "{pdb_id}_{chain_id}.chain_blast.xml",
                "{pdb_id}_{chain_id}.blast_chain.xml",
                "{pdb_id}_{chain_id}_chain.blast.xml"
            ],
            FileType.DOMAIN_BLAST: [
                "{pdb_id}_{chain_id}.domain_blast.xml",
                "{pdb_id}_{chain_id}.blast_domain.xml",
                "{pdb_id}_{chain_id}_domain.blast.xml"
            ],
            FileType.HHSEARCH: [
                "{pdb_id}_{chain_id}.{reference}.hhsearch.xml",
                "{pdb_id}_{chain_id}.hhsearch.{reference}.xml",
                "{pdb_id}_{chain_id}.{reference}.hh_summ.xml"
            ],
            FileType.HHR: [
                "{pdb_id}_{chain_id}.{reference}.hhr",
                "{pdb_id}_{chain_id}.hhr"
            ],
            FileType.SELF_COMPARISON: [
                "{pdb_id}_{chain_id}.self_comp.xml",
                "{pdb_id}_{chain_id}.self_comparison.xml",
                "{pdb_id}_{chain_id}_self_comp.xml"
            ],
            FileType.DOMAIN_SUMMARY: [
                "{pdb_id}_{chain_id}.{reference}.domain_summary.xml",
                "{pdb_id}_{chain_id}.{reference}.domain_summary.blast_only.xml",
                "{pdb_id}_{chain_id}.domain_summ.xml"
            ],
            FileType.DOMAIN_PARTITION: [
                "{pdb_id}_{chain_id}.{reference}.domains.xml",
                "{pdb_id}_{chain_id}.{reference}.domains_v14.xml",
                "{pdb_id}_{chain_id}.domains.xml"
            ],
            FileType.PROFILE: [
                "{pdb_id}_{chain_id}.{reference}.hhm",
                "{pdb_id}_{chain_id}.hhm"
            ]
        }

    def find_file(self, pdb_id: str, chain_id: str, file_type: FileType,
                  use_cache: bool = True) -> SearchResult:
        """
        Find a specific file type for a protein.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            file_type: Type of file to find
            use_cache: Whether to use cached results

        Returns:
            SearchResult with found locations
        """
        protein_id = f"{pdb_id}_{chain_id}"
        cache_key = f"{protein_id}:{file_type.value}"

        # Check cache
        if use_cache and cache_key in self._cache:
            result, timestamp = self._cache[cache_key]
            if (datetime.now() - timestamp).total_seconds() < self.cache_ttl:
                self.logger.debug(f"Using cached result for {cache_key}")
                return result

        # Create search result
        result = SearchResult(
            requested_type=file_type,
            protein_id=protein_id
        )

        # Search strategies in order
        strategies = [
            self._search_database,
            self._search_standard_locations,
            self._search_alternative_locations
        ]

        for strategy in strategies:
            try:
                locations = strategy(pdb_id, chain_id, file_type)
                result.locations.extend(locations)

                # Stop if found valid files
                if any(loc.is_valid for loc in locations):
                    break

            except Exception as e:
                error_msg = f"Error in {strategy.__name__}: {str(e)}"
                result.errors.append(error_msg)
                self.logger.warning(error_msg)

        # Cache result
        self._cache[cache_key] = (result, datetime.now())

        # Log result
        if result.found:
            best = result.best_location
            self.logger.info(
                f"Found {file_type.value} for {protein_id} at {best.path} "
                f"(source: {best.source})"
            )
        else:
            self.logger.warning(
                f"No {file_type.value} found for {protein_id}. "
                f"Searched {len(result.search_paths)} locations"
            )

        return result

    def _search_database(self, pdb_id: str, chain_id: str,
                        file_type: FileType) -> List[FileLocation]:
        """Search database for file paths"""
        locations = []

        query = """
        SELECT pf.file_path, pf.file_exists, pf.file_size, pf.last_checked
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
        AND pf.file_type = %s
        AND pf.file_exists = TRUE
        ORDER BY pf.id DESC
        """

        try:
            rows = self.db.execute_query(query, (pdb_id, chain_id, file_type.value))

            for row in rows:
                db_path = row[0]

                # Resolve path
                if os.path.isabs(db_path):
                    full_path = Path(db_path)
                else:
                    full_path = self.job_dump_dir / db_path

                full_path = full_path.resolve()

                # Validate existence
                exists = full_path.exists()
                size = full_path.stat().st_size if exists else row[2] or 0

                location = FileLocation(
                    path=full_path,
                    exists=exists,
                    size=size,
                    modified_time=row[3],
                    source="database",
                    file_type=file_type,
                    metadata={"db_exists": row[1], "db_size": row[2]}
                )

                locations.append(location)

                # Only keep first valid result from database
                if location.is_valid:
                    break

        except Exception as e:
            self.logger.error(f"Database search error: {e}")

        return locations

    def _search_standard_locations(self, pdb_id: str, chain_id: str,
                                  file_type: FileType) -> List[FileLocation]:
        """Search standard directory locations"""
        locations = []

        # Get directories to search
        search_dirs = self.dir_structure.get(file_type, ["."])

        # Get naming patterns
        patterns = self.naming_patterns.get(file_type, [])

        for dir_name in search_dirs:
            search_dir = self.job_dump_dir / dir_name

            if not search_dir.exists():
                continue

            for pattern in patterns:
                # Format pattern
                filename = pattern.format(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    reference=self.reference
                )

                file_path = search_dir / filename

                # Check if exists
                if file_path.exists():
                    stat = file_path.stat()
                    location = FileLocation(
                        path=file_path,
                        exists=True,
                        size=stat.st_size,
                        modified_time=datetime.fromtimestamp(stat.st_mtime),
                        source="standard",
                        file_type=file_type
                    )
                    locations.append(location)

                    # Log non-standard naming
                    if pattern != patterns[0]:
                        self.logger.debug(
                            f"Found {file_type.value} with non-standard "
                            f"naming: {filename}"
                        )

        return locations

    def _search_alternative_locations(self, pdb_id: str, chain_id: str,
                                     file_type: FileType) -> List[FileLocation]:
        """Search alternative/legacy locations"""
        locations = []
        protein_id = f"{pdb_id}_{chain_id}"

        # Alternative search paths
        alt_paths = [
            self.job_dump_dir / protein_id,  # Protein-specific directory
            self.job_dump_dir / "ecod_dump" / protein_id,  # Legacy location
            self.job_dump_dir / pdb_id / chain_id,  # Split by PDB/chain
        ]

        # Get naming patterns
        patterns = self.naming_patterns.get(file_type, [])

        for alt_dir in alt_paths:
            if not alt_dir.exists():
                continue

            for pattern in patterns:
                filename = pattern.format(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    reference=self.reference
                )

                file_path = alt_dir / filename

                if file_path.exists():
                    stat = file_path.stat()
                    location = FileLocation(
                        path=file_path,
                        exists=True,
                        size=stat.st_size,
                        modified_time=datetime.fromtimestamp(stat.st_mtime),
                        source="alternative",
                        file_type=file_type,
                        metadata={"alt_dir": str(alt_dir)}
                    )
                    locations.append(location)

        return locations

    def find_blast_file(self, pdb_id: str, chain_id: str,
                       blast_type: str = "domain") -> Optional[Path]:
        """
        Find BLAST result file with validation.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            blast_type: Type of BLAST ('domain' or 'chain')

        Returns:
            Path to valid BLAST file or None
        """
        file_type = FileType.DOMAIN_BLAST if blast_type == "domain" else FileType.CHAIN_BLAST
        result = self.find_file(pdb_id, chain_id, file_type)

        if result.best_location:
            # Validate BLAST file has hits
            if self._validate_blast_file(result.best_location.path):
                return result.best_location.path
            else:
                self.logger.warning(
                    f"BLAST file exists but has no valid hits: "
                    f"{result.best_location.path}"
                )

        return None

    def find_hhsearch_file(self, pdb_id: str, chain_id: str,
                          strict: bool = True) -> Optional[Path]:
        """
        Find HHSearch file with optional strict naming enforcement.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            strict: Whether to enforce standard naming only

        Returns:
            Path to HHSearch file or None
        """
        result = self.find_file(pdb_id, chain_id, FileType.HHSEARCH)

        if not result.found:
            # Try HHR format as fallback
            result = self.find_file(pdb_id, chain_id, FileType.HHR)

        if result.best_location:
            # If strict mode, only accept standard naming
            if strict:
                expected_name = f"{pdb_id}_{chain_id}.{self.reference}.hhsearch.xml"
                if result.best_location.path.name != expected_name:
                    self.logger.warning(
                        f"Found non-standard HHSearch file: "
                        f"{result.best_location.path.name}. "
                        f"Expected: {expected_name}"
                    )
                    return None

            return result.best_location.path

        return None

    def find_sequence_file(self, pdb_id: str, chain_id: str) -> Optional[Path]:
        """Find FASTA sequence file"""
        result = self.find_file(pdb_id, chain_id, FileType.FASTA)

        if result.best_location:
            return result.best_location.path

        return None

    def find_all_evidence_files(self, pdb_id: str, chain_id: str,
                               include_types: Optional[List[FileType]] = None) -> Dict[FileType, Path]:
        """
        Find all available evidence files for a protein.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            include_types: Specific file types to include (default: all)

        Returns:
            Dictionary mapping file types to paths
        """
        if include_types is None:
            include_types = [
                FileType.FASTA,
                FileType.CHAIN_BLAST,
                FileType.DOMAIN_BLAST,
                FileType.HHSEARCH,
                FileType.SELF_COMPARISON
            ]

        found_files = {}

        for file_type in include_types:
            result = self.find_file(pdb_id, chain_id, file_type)
            if result.best_location:
                found_files[file_type] = result.best_location.path

        self.logger.info(
            f"Found {len(found_files)}/{len(include_types)} evidence files "
            f"for {pdb_id}_{chain_id}"
        )

        return found_files

    def _validate_blast_file(self, blast_path: Path) -> bool:
        """Check if BLAST file contains valid hits"""
        try:
            import xml.etree.ElementTree as ET
            tree = ET.parse(blast_path)
            root = tree.getroot()

            # Check for hits
            hits = root.findall(".//Hit")
            return len(hits) > 0

        except Exception as e:
            self.logger.error(f"Error validating BLAST file {blast_path}: {e}")
            return False

    def check_file_ages(self, pdb_id: str, chain_id: str,
                       max_age_days: int = 30) -> Dict[FileType, bool]:
        """
        Check if evidence files are recent enough.

        Args:
            pdb_id: PDB identifier
            chain_id: Chain identifier
            max_age_days: Maximum acceptable age in days

        Returns:
            Dictionary mapping file types to freshness status
        """
        freshness = {}

        for file_type in FileType:
            result = self.find_file(pdb_id, chain_id, file_type)

            if result.best_location:
                age = result.best_location.age_days
                freshness[file_type] = age is not None and age <= max_age_days
            else:
                freshness[file_type] = False

        return freshness

    def clear_cache(self, pattern: Optional[str] = None):
        """
        Clear cache entries.

        Args:
            pattern: Optional pattern to match cache keys
        """
        if pattern:
            keys_to_remove = [
                key for key in self._cache.keys()
                if pattern in key
            ]
            for key in keys_to_remove:
                del self._cache[key]

            self.logger.info(f"Cleared {len(keys_to_remove)} cache entries")
        else:
            self._cache.clear()
            self.logger.info("Cleared all cache entries")

    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        total_entries = len(self._cache)

        # Count by file type
        type_counts = {}
        for key in self._cache:
            file_type = key.split(":")[-1]
            type_counts[file_type] = type_counts.get(file_type, 0) + 1

        # Calculate age distribution
        ages = []
        for _, (_, timestamp) in self._cache.items():
            age = (datetime.now() - timestamp).total_seconds()
            ages.append(age)

        return {
            "total_entries": total_entries,
            "type_counts": type_counts,
            "oldest_entry_seconds": max(ages) if ages else 0,
            "newest_entry_seconds": min(ages) if ages else 0,
            "average_age_seconds": sum(ages) / len(ages) if ages else 0
        }

    def register_custom_location(self, file_type: FileType,
                               directory: Union[str, Path]):
        """
        Register a custom directory for a file type.

        Args:
            file_type: File type to register directory for
            directory: Directory path to add
        """
        if file_type not in self.dir_structure:
            self.dir_structure[file_type] = []

        dir_path = str(Path(directory).relative_to(self.job_dump_dir))
        if dir_path not in self.dir_structure[file_type]:
            self.dir_structure[file_type].insert(0, dir_path)
            self.logger.info(
                f"Registered custom directory for {file_type.value}: {dir_path}"
            )

    def register_custom_pattern(self, file_type: FileType, pattern: str):
        """
        Register a custom naming pattern for a file type.

        Args:
            file_type: File type to register pattern for
            pattern: Naming pattern with placeholders
        """
        if file_type not in self.naming_patterns:
            self.naming_patterns[file_type] = []

        if pattern not in self.naming_patterns[file_type]:
            self.naming_patterns[file_type].insert(0, pattern)
            self.logger.info(
                f"Registered custom pattern for {file_type.value}: {pattern}"
            )
