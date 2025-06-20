# mini/core/writer.py (integrated version)
"""Write domain partition results with comprehensive provenance tracking"""

import xml.etree.ElementTree as ET
import subprocess
import hashlib
import os
from datetime import datetime
from typing import List, Optional
from .models import Domain, PartitionMetadata, DomainLayout

def get_git_commit_hash() -> str:
    """Get current git commit hash for version tracking"""
    try:
        result = subprocess.run(['git', 'rev-parse', 'HEAD'],
                              capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"

def get_git_version() -> str:
    """Get semantic version + commit hash"""
    try:
        # Try to get latest git tag
        result = subprocess.run(['git', 'describe', '--tags', '--always'],
                              capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return f"mini_pyecod_unknown_{get_git_commit_hash()[:8]}"

def calculate_file_hash(file_path: str) -> Optional[str]:
    """Calculate SHA256 hash of source file"""
    if not os.path.exists(file_path):
        return None

    hash_sha256 = hashlib.sha256()
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_sha256.update(chunk)
        return hash_sha256.hexdigest()
    except IOError:
        return None

def write_domain_partition(domains: List[Domain],
                          metadata: PartitionMetadata,
                          output_path: str,
                          reference: str = "mini_pyecod") -> None:
    """Write domains to XML file with comprehensive provenance tracking"""

    # Ensure metadata has versioning information
    if metadata.algorithm_version is None:
        metadata.algorithm_version = get_git_version()
    if metadata.git_commit_hash is None:
        metadata.git_commit_hash = get_git_commit_hash()
    if metadata.processing_timestamp is None:
        metadata.processing_timestamp = datetime.now()

    # Calculate source file hash if path provided
    if metadata.source_domain_summary_path and not metadata.source_domain_summary_hash:
        metadata.source_domain_summary_hash = calculate_file_hash(metadata.source_domain_summary_path)

    # Create root element
    root = ET.Element("domain_partition")
    root.set("pdb_id", metadata.pdb_id)
    root.set("chain_id", metadata.chain_id)
    root.set("reference", reference)
    root.set("is_classified", "true" if domains else "false")

    # Add comprehensive metadata
    metadata_elem = ET.SubElement(root, "metadata")

    # Version tracking
    version_elem = ET.SubElement(metadata_elem, "version")
    version_elem.set("algorithm", metadata.algorithm_version)
    version_elem.set("git_commit", metadata.git_commit_hash)
    version_elem.set("timestamp", metadata.processing_timestamp.isoformat())

    # Source provenance
    if metadata.source_domain_summary_path:
        source_elem = ET.SubElement(metadata_elem, "source")
        source_elem.set("domain_summary_path", metadata.source_domain_summary_path)
        if metadata.source_domain_summary_hash:
            source_elem.set("file_hash", metadata.source_domain_summary_hash)
        if metadata.batch_id:
            source_elem.set("batch_id", metadata.batch_id)

    # Processing parameters
    if metadata.process_parameters:
        params_elem = ET.SubElement(metadata_elem, "parameters")
        for key, value in metadata.process_parameters.items():
            param_elem = ET.SubElement(params_elem, "parameter")
            param_elem.set("name", key)
            param_elem.set("value", str(value))

    # Sequence statistics
    if metadata.sequence_length is not None:
        stats_elem = ET.SubElement(metadata_elem, "statistics")
        stats_elem.set("sequence_length", str(metadata.sequence_length))
        stats_elem.set("domain_count", str(len(domains)))

        # Calculate coverage statistics
        if domains:
            total_assigned = sum(d.length for d in domains)
            coverage = total_assigned / metadata.sequence_length
            stats_elem.set("total_coverage", f"{coverage:.4f}")
            stats_elem.set("residues_assigned", str(total_assigned))

            # Optimization statistics
            optimized_domains = [d for d in domains if d.was_optimized()]
            stats_elem.set("domains_optimized", str(len(optimized_domains)))
            if optimized_domains:
                total_position_change = sum(
                    d.length - len(set(d.original_range.to_positions_simple()))
                    for d in optimized_domains
                )
                stats_elem.set("optimization_position_change", str(total_position_change))

    # Domain definitions
    domains_elem = ET.SubElement(root, "domains")

    for domain in domains:
        d_elem = ET.SubElement(domains_elem, "domain")
        d_elem.set("id", domain.id)
        d_elem.set("range", str(domain.range))
        d_elem.set("family", domain.family)
        d_elem.set("source", domain.source)
        d_elem.set("evidence_count", str(domain.evidence_count))
        d_elem.set("is_discontinuous", str(domain.range.is_discontinuous).lower())

        # Classification hierarchy
        if domain.t_group:
            d_elem.set("t_group", domain.t_group)
        if domain.h_group:
            d_elem.set("h_group", domain.h_group)
        if domain.x_group:
            d_elem.set("x_group", domain.x_group)

        # Domain-level provenance
        if domain.confidence_score is not None:
            d_elem.set("confidence", f"{domain.confidence_score:.4f}")

        if domain.reference_ecod_domain_id:
            d_elem.set("reference_ecod_domain_id", domain.reference_ecod_domain_id)

        # Primary evidence details
        if domain.primary_evidence:
            evidence_elem = ET.SubElement(d_elem, "primary_evidence")
            evidence_elem.set("source_type", domain.primary_evidence.type)

            # Source identification
            source_id = domain.primary_evidence.source_pdb
            if domain.primary_evidence.source_chain_id:
                source_id += f"_{domain.primary_evidence.source_chain_id}"
            evidence_elem.set("source_id", source_id)

            if domain.primary_evidence.domain_id:
                evidence_elem.set("domain_id", domain.primary_evidence.domain_id)

            if domain.primary_evidence.evalue is not None:
                evidence_elem.set("evalue", str(domain.primary_evidence.evalue))

            evidence_elem.set("evidence_range", str(domain.primary_evidence.query_range))

            if domain.primary_evidence.hit_range:
                evidence_elem.set("hit_range", str(domain.primary_evidence.hit_range))

            if domain.primary_evidence.hsp_count is not None:
                evidence_elem.set("hsp_count", str(domain.primary_evidence.hsp_count))

            evidence_elem.set("discontinuous", str(domain.primary_evidence.discontinuous).lower())

        # Boundary optimization tracking
        if domain.was_optimized():
            optimization_elem = ET.SubElement(d_elem, "boundary_optimization")
            optimization_elem.set("original_range", str(domain.original_range))
            optimization_elem.set("optimized_range", str(domain.range))

            if domain.optimization_actions:
                optimization_elem.set("actions", ",".join(domain.optimization_actions))

            # Calculate position change
            original_positions = len(set(domain.original_range.to_positions_simple()))
            position_change = domain.length - original_positions
            optimization_elem.set("position_change", str(position_change))

    # Write and calculate output file hash
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    tree.write(output_path, encoding="utf-8", xml_declaration=True)

    # Update metadata with output file information
    metadata.output_xml_path = output_path
    metadata.output_xml_hash = calculate_file_hash(output_path)

def write_domain_partition_from_layout(layout: DomainLayout,
                                     metadata: PartitionMetadata,
                                     output_path: str) -> None:
    """Convenience wrapper to write from DomainLayout with automatic statistics"""

    # Update metadata with layout statistics
    if metadata.sequence_length is None:
        metadata.sequence_length = layout.sequence_length

    # Add coverage statistics to process parameters
    coverage_stats = layout.get_coverage_stats()
    metadata.process_parameters.update({
        'boundary_optimization_enabled': True,
        'final_coverage_percent': coverage_stats['coverage_percent'],
        'num_gaps_remaining': coverage_stats['num_gaps'],
        'small_fragments_merged': coverage_stats['small_fragments'],
        'large_gaps_remaining': coverage_stats['large_gaps']
    })

    write_domain_partition(
        domains=layout.domains,
        metadata=metadata,
        output_path=output_path
    )

def create_metadata_from_batch(pdb_id: str, chain_id: str,
                             batch_path: str,
                             batch_id: Optional[str] = None) -> PartitionMetadata:
    """Create PartitionMetadata from batch processing context"""

    # Construct domain summary path
    domain_summary_path = os.path.join(batch_path, "domains",
                                     f"{pdb_id}_{chain_id}.develop291.domain_summary.xml")

    # Extract batch_id from path if not provided
    if batch_id is None:
        batch_id = os.path.basename(batch_path)

    return PartitionMetadata(
        pdb_id=pdb_id,
        chain_id=chain_id,
        source_domain_summary_path=domain_summary_path,
        batch_id=batch_id
    )

# Legacy compatibility function
def write_domain_partition_with_provenance(domains: List[Domain],
                                         pdb_id: str, chain_id: str,
                                         output_path: str, batch_path: str):
    """Legacy wrapper for backward compatibility"""

    metadata = create_metadata_from_batch(pdb_id, chain_id, batch_path)
    write_domain_partition(domains, metadata, output_path)
