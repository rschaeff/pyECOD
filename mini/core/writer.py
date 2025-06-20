# mini/core/writer.py (enhanced with comprehensive evidence reporting)
"""Write domain partition results with comprehensive provenance tracking and evidence metrics"""

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

def calculate_original_hhsearch_probability(evidence) -> Optional[float]:
    """Calculate original HHsearch probability from confidence score"""
    if evidence.type != 'hhsearch':
        return None

    # Reverse the confidence calculation from evidence_utils.py
    # confidence = base_confidence * 0.95 (type multiplier)
    # base_confidence = probability (if probability was 0-1 scale)

    base_confidence = evidence.confidence / 0.95  # Reverse type multiplier

    # Convert back to percentage if it was originally 0-100 scale
    # Since we normalized to 0-1 in evidence_utils, multiply by 100
    original_probability = base_confidence * 100

    # Clamp to reasonable range
    return max(0.0, min(100.0, original_probability))

def calculate_reference_coverage(evidence) -> Optional[float]:
    """Calculate reference coverage from evidence"""
    if not evidence.hit_range or not evidence.reference_length:
        return None

    if evidence.reference_length <= 0:
        return None

    return evidence.hit_range.total_length / evidence.reference_length

def format_evalue(evalue: Optional[float]) -> str:
    """Format E-value for display"""
    if evalue is None:
        return "N/A"

    if evalue == 0.0:
        return "0.0"
    elif evalue < 1e-10:
        return f"{evalue:.1e}"
    elif evalue < 1e-3:
        return f"{evalue:.6f}"
    elif evalue < 1.0:
        return f"{evalue:.4f}"
    else:
        return f"{evalue:.1f}"

def write_domain_partition(domains: List[Domain],
                          metadata: PartitionMetadata,
                          output_path: str,
                          reference: str = "mini_pyecod") -> None:
    """Write domains to XML file with comprehensive provenance tracking and evidence metrics"""

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

        # ENHANCED: Primary evidence details with comprehensive metrics
        if domain.primary_evidence:
            evidence_elem = ET.SubElement(d_elem, "primary_evidence")
            evidence = domain.primary_evidence

            evidence_elem.set("source_type", evidence.type)

            # Source identification
            source_id = evidence.source_pdb
            if evidence.source_chain_id:
                source_id += f"_{evidence.source_chain_id}"
            evidence_elem.set("source_id", source_id)

            if evidence.domain_id:
                evidence_elem.set("domain_id", evidence.domain_id)

            # ENHANCED: Evidence type-specific metrics
            if evidence.type == 'hhsearch':
                # Report original HHsearch probability
                original_prob = calculate_original_hhsearch_probability(evidence)
                if original_prob is not None:
                    evidence_elem.set("hh_probability", f"{original_prob:.1f}")

                # Report E-value if available
                if evidence.evalue is not None:
                    evidence_elem.set("evalue", format_evalue(evidence.evalue))

            elif evidence.type in ['domain_blast', 'chain_blast', 'chain_blast_decomposed']:
                # Report E-value for BLAST evidence
                if evidence.evalue is not None:
                    evidence_elem.set("evalue", format_evalue(evidence.evalue))

            # Confidence score (processed value)
            evidence_elem.set("confidence", f"{evidence.confidence:.4f}")

            # Sequence ranges
            evidence_elem.set("evidence_range", str(evidence.query_range))

            if evidence.hit_range:
                evidence_elem.set("hit_range", str(evidence.hit_range))

                # ENHANCED: Reference coverage calculation and reporting
                ref_coverage = calculate_reference_coverage(evidence)
                if ref_coverage is not None:
                    evidence_elem.set("reference_coverage", f"{ref_coverage:.3f}")
                    evidence_elem.set("reference_coverage_percent", f"{ref_coverage:.1%}")

            # Reference length for context
            if evidence.reference_length:
                evidence_elem.set("reference_length", str(evidence.reference_length))

            # Additional evidence metadata
            if evidence.hsp_count is not None:
                evidence_elem.set("hsp_count", str(evidence.hsp_count))

            evidence_elem.set("discontinuous", str(evidence.discontinuous).lower())

            # ENHANCED: Quality assessment flags
            quality_flags = []

            # Check confidence level
            if evidence.confidence >= 0.8:
                quality_flags.append("high_confidence")
            elif evidence.confidence < 0.5:
                quality_flags.append("low_confidence")

            # Check reference coverage
            if ref_coverage is not None:
                if ref_coverage >= 0.8:
                    quality_flags.append("high_coverage")
                elif ref_coverage < 0.5:
                    quality_flags.append("low_coverage")

            # Check E-value quality
            if evidence.evalue is not None:
                if evidence.evalue < 1e-10:
                    quality_flags.append("excellent_evalue")
                elif evidence.evalue > 1.0:
                    quality_flags.append("poor_evalue")

            if quality_flags:
                evidence_elem.set("quality_flags", ",".join(quality_flags))

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

        # ENHANCED: All evidence summary (if multiple evidence items)
        if len(domain.evidence_items) > 1:
            all_evidence_elem = ET.SubElement(d_elem, "supporting_evidence")
            all_evidence_elem.set("count", str(len(domain.evidence_items)))

            evidence_types = [e.type for e in domain.evidence_items]
            all_evidence_elem.set("types", ",".join(set(evidence_types)))

            # Average confidence across all evidence
            avg_confidence = sum(e.confidence for e in domain.evidence_items) / len(domain.evidence_items)
            all_evidence_elem.set("average_confidence", f"{avg_confidence:.4f}")

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

# ENHANCED: Evidence quality reporting
def generate_evidence_quality_report(domains: List[Domain]) -> str:
    """Generate a human-readable evidence quality report"""

    if not domains:
        return "No domains found."

    report_lines = []
    report_lines.append("EVIDENCE QUALITY REPORT")
    report_lines.append("=" * 50)

    for i, domain in enumerate(domains, 1):
        evidence = domain.primary_evidence
        if not evidence:
            continue

        report_lines.append(f"\nDomain {i}: {domain.family} ({domain.range})")
        report_lines.append(f"  Source: {evidence.type}")
        report_lines.append(f"  Confidence: {evidence.confidence:.3f}")

        # Type-specific metrics
        if evidence.type == 'hhsearch':
            original_prob = calculate_original_hhsearch_probability(evidence)
            if original_prob is not None:
                report_lines.append(f"  HH Probability: {original_prob:.1f}%")

        if evidence.evalue is not None:
            report_lines.append(f"  E-value: {format_evalue(evidence.evalue)}")

        # Reference coverage
        ref_coverage = calculate_reference_coverage(evidence)
        if ref_coverage is not None:
            report_lines.append(f"  Reference Coverage: {ref_coverage:.1%} ({evidence.hit_range.total_length}/{evidence.reference_length})")
        else:
            report_lines.append(f"  Reference Coverage: Not available")

        # Quality assessment
        quality_issues = []
        if evidence.confidence < 0.5:
            quality_issues.append("LOW CONFIDENCE")
        if ref_coverage is not None and ref_coverage < 0.5:
            quality_issues.append("POOR REFERENCE COVERAGE")
        if evidence.evalue is not None and evidence.evalue > 1.0:
            quality_issues.append("POOR E-VALUE")

        if quality_issues:
            report_lines.append(f"  ⚠️  Quality Issues: {', '.join(quality_issues)}")
        else:
            report_lines.append(f"  ✓ Quality: Good")

    return "\n".join(report_lines)

# Legacy compatibility function
def write_domain_partition_with_provenance(domains: List[Domain],
                                         pdb_id: str, chain_id: str,
                                         output_path: str, batch_path: str):
    """Legacy wrapper for backward compatibility"""

    metadata = create_metadata_from_batch(pdb_id, chain_id, batch_path)
    write_domain_partition(domains, metadata, output_path)
