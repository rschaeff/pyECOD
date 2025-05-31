# ecod/pipelines/domain_analysis/decomposition/service.py
"""
Chain Blast Decomposition Service

Handles chain-level BLAST hits and decomposes them into domain-level evidence.
This is a separate service due to complexity and need for careful monitoring.
"""

import logging
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass
from enum import Enum
import statistics

class DecompositionStatus(Enum):
    SUCCESS = "success"
    FAILED_SHORT_DOMAINS = "failed_short_domains"
    FAILED_POOR_COVERAGE = "failed_poor_coverage"
    FAILED_NO_ARCHITECTURE = "failed_no_architecture"
    FAILED_ALIGNMENT_ISSUES = "failed_alignment_issues"
    FAILED_UNKNOWN = "failed_unknown"

@dataclass
class DecompositionResult:
    """Result of chain blast decomposition"""
    status: DecompositionStatus
    success: bool
    domains: List[Dict[str, Any]] = None
    quality_metrics: Dict[str, float] = None
    failure_reason: str = ""
    template_info: Dict[str, Any] = None
    quality_warnings: List[str] = None
    
    def __post_init__(self):
        if self.domains is None:
            self.domains = []
        if self.quality_metrics is None:
            self.quality_metrics = {}
        if self.quality_warnings is None:
            self.quality_warnings = []

@dataclass 
class DecompositionConfig:
    """Configuration for chain blast decomposition"""
    
    # Quality thresholds for overriding chain blast precedence
    min_domain_size: int = 30
    min_reference_coverage: float = 0.70  # 70% coverage of reference domain
    min_alignment_identity: float = 0.25   # 25% sequence identity
    min_template_quality: float = 3.0      # Resolution threshold
    
    # Decomposition parameters
    allow_domain_gaps: bool = True
    max_gap_size: int = 20
    min_domain_separation: int = 10
    
    # Fallback thresholds - when to defer to individual evidence
    max_short_domains_ratio: float = 0.3   # If >30% domains are too short, fail
    min_average_coverage: float = 0.50     # Average coverage across all domains
    
    # Monitoring
    log_all_attempts: bool = True
    log_failures_detail: bool = True


class ChainBlastDecompositionService:
    """
    Service for decomposing chain blast hits into domain-level evidence.
    
    This service handles the complex task of taking protein-protein alignments
    and converting them into individual domain boundaries and classifications.
    """
    
    def __init__(self, config: DecompositionConfig, context=None, perl_script_path: Optional[str] = None):
        """FIXED: Get database manager from context"""

        self.config = config
        self.context = context
        self.perl_script_path = perl_script_path
        self.logger = logging.getLogger(__name__)

        # Get database manager from context
        if context and hasattr(context, 'db_manager'):
            self.db_manager = context.db_manager
            self.logger.info("Database manager obtained from context.db_manager")
        elif context and hasattr(context, 'db'):
            self.db_manager = context.db
            self.logger.info("Database manager obtained from context.db")
        else:
            self.db_manager = None
            self.logger.warning("No database manager found in context")
    
    def decompose_chain_blast_hits(self, chain_hits: List[Dict[str, Any]], 
                                 sequence_length: int,
                                 protein_id: str = "") -> List[DecompositionResult]:
        """
        Main entry point for decomposing chain blast hits.
        
        Args:
            chain_hits: List of chain blast hit dictionaries
            sequence_length: Query protein sequence length
            protein_id: Query protein identifier for logging
            
        Returns:
            List of DecompositionResult objects, one per chain hit
        """
        
        results = []
        self.logger.info(f"Starting decomposition of {len(chain_hits)} chain blast hits for {protein_id}")
        
        for i, hit in enumerate(chain_hits):
            self.stats['total_attempts'] += 1
            
            try:
                result = self._decompose_single_hit(hit, sequence_length, f"{protein_id}_hit_{i}")
                results.append(result)
                
                if result.success:
                    self.stats['successful_decompositions'] += 1
                    self.logger.debug(f"Successful decomposition: {hit.get('pdb_id')}_{hit.get('chain_id')} "
                                    f"-> {len(result.domains)} domains")
                else:
                    self._update_failure_stats(result.status)
                    self.logger.debug(f"Failed decomposition: {hit.get('pdb_id')}_{hit.get('chain_id')} "
                                    f"({result.failure_reason})")
                    
            except Exception as e:
                error_result = DecompositionResult(
                    status=DecompositionStatus.FAILED_UNKNOWN,
                    success=False,
                    failure_reason=f"Unexpected error: {str(e)}"
                )
                results.append(error_result)
                self.logger.error(f"Error decomposing chain blast hit: {e}")
        
        self._log_decomposition_summary(results, protein_id)
        return results
    
    def _decompose_single_hit(self, hit: Dict[str, Any], sequence_length: int, 
                            hit_id: str) -> DecompositionResult:
        """Decompose a single chain blast hit"""
        
        template_pdb = hit.get('pdb_id', '')
        template_chain = hit.get('chain_id', '')
        
        # Step 1: Get template domain architecture
        template_domains = self._get_template_architecture(template_pdb, template_chain)
        if not template_domains:
            return DecompositionResult(
                status=DecompositionStatus.FAILED_NO_ARCHITECTURE,
                success=False,
                failure_reason=f"No domain architecture found for {template_pdb}_{template_chain}"
            )
        
        # Step 2: Project domains onto query sequence
        projected_domains = self._project_domains_to_query(hit, template_domains, sequence_length)
        if not projected_domains:
            return DecompositionResult(
                status=DecompositionStatus.FAILED_ALIGNMENT_ISSUES,
                success=False,
                failure_reason="Failed to project template domains onto query sequence"
            )
        
        # Step 3: Quality assessment - check for short domains
        short_domains = [d for d in projected_domains if d['end'] - d['start'] + 1 < self.config.min_domain_size]
        if len(short_domains) / len(projected_domains) > self.config.max_short_domains_ratio:
            return DecompositionResult(
                status=DecompositionStatus.FAILED_SHORT_DOMAINS,
                success=False,
                failure_reason=f"{len(short_domains)}/{len(projected_domains)} domains too short "
                              f"(< {self.config.min_domain_size} residues)",
                domains=projected_domains,  # Include for debugging
                quality_warnings=[f"Short domain: {d['start']}-{d['end']}" for d in short_domains]
            )
        
        # Step 4: Quality assessment - check reference coverage
        coverage_issues = []
        for domain in projected_domains:
            coverage = domain.get('reference_coverage', 0.0)
            if coverage < self.config.min_reference_coverage:
                coverage_issues.append(f"Domain {domain.get('domain_id', 'unknown')}: "
                                     f"{coverage:.2f} < {self.config.min_reference_coverage}")
        
        if len(coverage_issues) / len(projected_domains) > 0.5:  # >50% have poor coverage
            return DecompositionResult(
                status=DecompositionStatus.FAILED_POOR_COVERAGE,
                success=False,
                failure_reason=f"Poor reference coverage in {len(coverage_issues)} domains",
                domains=projected_domains,
                quality_warnings=coverage_issues
            )
        
        # Step 5: Calculate quality metrics
        quality_metrics = self._calculate_decomposition_quality(projected_domains, hit, template_domains)
        
        # Success!
        return DecompositionResult(
            status=DecompositionStatus.SUCCESS,
            success=True,
            domains=projected_domains,
            quality_metrics=quality_metrics,
            template_info={
                'pdb_id': template_pdb,
                'chain_id': template_chain,
                'domain_count': len(template_domains),
                'template_domains': template_domains
            },
            quality_warnings=coverage_issues if coverage_issues else []
        )
    
    def _get_template_architecture(self, pdb_id: str, chain_id: str) -> List[Dict[str, Any]]:
        """
        Get domain architecture for template structure - FIXED for actual database schema
        """
        if not self.db_manager:
            self.logger.warning("No database manager available for template architecture lookup")
            return []

        try:
            # FIXED: Updated query for actual database schema
            query = """
                SELECT d.domain_id, d.ecod_domain_id, d.range,
                       d.t_group, d.h_group, d.x_group, d.a_group,
                       d.is_manual_rep, d.is_f70, d.is_f40, d.is_f99
                FROM pdb_analysis.domain d
                JOIN pdb_analysis.protein p ON d.protein_id = p.id
                WHERE p.pdb_id = %s AND p.chain_id = %s
                ORDER BY d.range
            """

            if hasattr(self.db_manager, 'execute_dict_query'):
                results = self.db_manager.execute_dict_query(query, (pdb_id, chain_id))
            else:
                cursor = self.db_manager.execute(query, (pdb_id, chain_id))
                rows = cursor.fetchall()
                columns = [desc[0] for desc in cursor.description]
                results = [dict(zip(columns, row)) for row in rows]

            domains = []
            for row in results:
                # FIXED: Parse range format like "A:225-384" or "d2:7-68"
                range_str = row['range']
                start, end = self._parse_domain_range_with_chain(range_str, chain_id)

                if start == 0 and end == 0:
                    self.logger.warning(f"Could not parse range '{range_str}' for domain {row['domain_id']}")
                    continue

                domain = {
                    'domain_id': row['domain_id'],
                    'ecod_domain_id': row['ecod_domain_id'],
                    'start': start,
                    'end': end,
                    'range': range_str,
                    't_group': row['t_group'],
                    'h_group': row['h_group'],
                    'x_group': row['x_group'],
                    'a_group': row['a_group'],
                    'is_manual_rep': row.get('is_manual_rep', False),
                    'is_f70': row.get('is_f70', False),
                    'is_f40': row.get('is_f40', False),
                    'is_f99': row.get('is_f99', False)
                }
                domains.append(domain)
            
            self.logger.debug(f"Found {len(domains)} template domains for {pdb_id}_{chain_id}")
            return domains
            
        except Exception as e:
            self.logger.error(f"Error getting template architecture for {pdb_id}_{chain_id}: {e}")
            return []

    def _parse_domain_range_with_chain(self, range_str: str, expected_chain: str) -> Tuple[int, int]:
        """
        Parse domain range string with chain information

        Examples:
        - "A:225-384" -> (225, 384)
        - "d2:7-68" -> (7, 68)
        - "BI:81-266" -> (81, 266)
        """
        try:
            if ':' in range_str:
                # Format: "chain:start-end"
                chain_part, range_part = range_str.split(':', 1)

                # For now, we'll process any chain (could add chain validation later)
                if '-' in range_part:
                    start, end = range_part.split('-')
                    return int(start), int(end)
                else:
                    # Single position
                    pos = int(range_part)
                    return pos, pos
            else:
                # Format: "start-end" (no chain)
                if '-' in range_str:
                    start, end = range_str.split('-')
                    return int(start), int(end)
                else:
                    pos = int(range_str)
                    return pos, pos

        except ValueError as e:
            self.logger.warning(f"Could not parse range '{range_str}': {e}")
            return 0, 0
    
    def _project_domains_to_query(self, hit: Dict[str, Any], template_domains: List[Dict[str, Any]],
                                sequence_length: int) -> List[Dict[str, Any]]:
        """
        Project template domains onto query sequence using alignment coordinates.
        
        This is the core algorithmic challenge of chain blast decomposition.
        """
        
        # Parse alignment coordinates
        query_ranges = self._parse_alignment_ranges(hit.get('query_range', ''))
        hit_ranges = self._parse_alignment_ranges(hit.get('hit_range', ''))
        
        if not query_ranges or not hit_ranges or len(query_ranges) != len(hit_ranges):
            self.logger.warning(f"Invalid alignment ranges: query={query_ranges}, hit={hit_ranges}")
            return []
        
        projected_domains = []
        
        for template_domain in template_domains:
            template_start = template_domain['start']
            template_end = template_domain['end']
            
            # Find alignment segments that overlap with this template domain
            overlapping_segments = []
            
            for (q_start, q_end), (h_start, h_end) in zip(query_ranges, hit_ranges):
                # Check for overlap between template domain and alignment segment
                overlap_start = max(template_start, h_start)
                overlap_end = min(template_end, h_end)
                
                if overlap_start <= overlap_end:
                    # Calculate how this maps to query coordinates
                    template_offset = overlap_start - h_start
                    overlap_length = overlap_end - overlap_start + 1
                    
                    query_proj_start = q_start + template_offset
                    query_proj_end = query_proj_start + overlap_length - 1
                    
                    overlapping_segments.append({
                        'query_start': query_proj_start,
                        'query_end': query_proj_end,
                        'template_start': overlap_start,
                        'template_end': overlap_end,
                        'length': overlap_length
                    })
            
            if not overlapping_segments:
                continue
            
            # Merge overlapping segments and create projected domain
            merged_segments = self._merge_alignment_segments(overlapping_segments)
            
            if merged_segments:
                # Calculate reference coverage
                total_projected_length = sum(s['length'] for s in merged_segments)
                template_length = template_end - template_start + 1
                reference_coverage = total_projected_length / template_length
                
                # Create projected domain (use first/last coordinates if multiple segments)
                projected_start = min(s['query_start'] for s in merged_segments)
                projected_end = max(s['query_end'] for s in merged_segments)
                
                # Validate coordinates
                if (projected_start > 0 and projected_end <= sequence_length and 
                    projected_start <= projected_end):
                    
                    projected_domain = {
                        'start': projected_start,
                        'end': projected_end,
                        'domain_id': template_domain['domain_id'],
                        't_group': template_domain['t_group'],
                        'h_group': template_domain['h_group'],
                        'x_group': template_domain['x_group'],
                        'a_group': template_domain['a_group'],
                        'reference_coverage': reference_coverage,
                        'template_range': f"{template_start}-{template_end}",
                        'alignment_segments': merged_segments,
                        'projection_confidence': self._calculate_projection_confidence(
                            reference_coverage, len(merged_segments)
                        )
                    }
                    projected_domains.append(projected_domain)
        
        return projected_domains
    
    def _parse_alignment_ranges(self, range_str: str) -> List[Tuple[int, int]]:
        """Parse alignment range string into list of (start, end) tuples"""
        ranges = []
        if not range_str:
            return ranges
        
        for segment in range_str.split(','):
            segment = segment.strip()
            if '-' in segment:
                try:
                    start, end = segment.split('-')
                    ranges.append((int(start), int(end)))
                except ValueError:
                    continue
        
        return ranges
    
    def _parse_domain_range(self, range_str: str) -> Tuple[int, int]:
        """Parse domain range string into start, end positions"""
        try:
            if '-' in range_str:
                start, end = range_str.split('-')
                return int(start), int(end)
            else:
                pos = int(range_str)
                return pos, pos
        except ValueError:
            return 0, 0
    
    def _merge_alignment_segments(self, segments: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Merge overlapping alignment segments"""
        if not segments:
            return []
        
        # Sort by query start position
        sorted_segments = sorted(segments, key=lambda x: x['query_start'])
        merged = [sorted_segments[0]]
        
        for current in sorted_segments[1:]:
            last = merged[-1]
            
            # Check if segments can be merged (allow small gaps)
            if current['query_start'] <= last['query_end'] + self.config.max_gap_size:
                # Merge segments
                merged[-1] = {
                    'query_start': last['query_start'],
                    'query_end': max(last['query_end'], current['query_end']),
                    'template_start': min(last['template_start'], current['template_start']),
                    'template_end': max(last['template_end'], current['template_end']),
                    'length': max(last['query_end'], current['query_end']) - last['query_start'] + 1
                }
            else:
                merged.append(current)
        
        return merged
    
    def _calculate_projection_confidence(self, reference_coverage: float, 
                                       segment_count: int) -> float:
        """Calculate confidence in domain projection"""
        
        # Base confidence from reference coverage
        coverage_factor = min(1.0, reference_coverage / 0.8)  # Full confidence at 80% coverage
        
        # Penalty for fragmented alignments
        fragmentation_penalty = 1.0 / (1.0 + (segment_count - 1) * 0.2)  # Penalty for multiple segments
        
        return coverage_factor * fragmentation_penalty
    
    def _calculate_decomposition_quality(self, projected_domains: List[Dict[str, Any]], 
                                       hit: Dict[str, Any], 
                                       template_domains: List[Dict[str, Any]]) -> Dict[str, float]:
        """Calculate overall quality metrics for decomposition"""
        
        if not projected_domains:
            return {}
        
        # Coverage metrics
        coverages = [d['reference_coverage'] for d in projected_domains]
        avg_coverage = statistics.mean(coverages)
        min_coverage = min(coverages)
        
        # Size metrics  
        sizes = [d['end'] - d['start'] + 1 for d in projected_domains]
        avg_size = statistics.mean(sizes)
        min_size = min(sizes)
        
        # Projection confidence
        confidences = [d['projection_confidence'] for d in projected_domains]
        avg_confidence = statistics.mean(confidences)
        
        return {
            'average_reference_coverage': avg_coverage,
            'minimum_reference_coverage': min_coverage,
            'average_domain_size': avg_size,
            'minimum_domain_size': min_size,
            'average_projection_confidence': avg_confidence,
            'domain_count': len(projected_domains),
            'template_domain_count': len(template_domains),
            'decomposition_completeness': len(projected_domains) / len(template_domains),
            'alignment_evalue': hit.get('evalue', 1.0),
            'alignment_identity': hit.get('identity', 0.0)
        }
    
    def _update_failure_stats(self, status: DecompositionStatus):
        """Update failure statistics"""
        if status == DecompositionStatus.FAILED_SHORT_DOMAINS:
            self.stats['failed_short_domains'] += 1
        elif status == DecompositionStatus.FAILED_POOR_COVERAGE:
            self.stats['failed_poor_coverage'] += 1
        elif status == DecompositionStatus.FAILED_NO_ARCHITECTURE:
            self.stats['failed_no_architecture'] += 1
        elif status == DecompositionStatus.FAILED_ALIGNMENT_ISSUES:
            self.stats['failed_alignment_issues'] += 1
    
    def _log_decomposition_summary(self, results: List[DecompositionResult], protein_id: str):
        """Log summary of decomposition results"""
        
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        
        if successful:
            total_domains = sum(len(r.domains) for r in successful)
            avg_domains = total_domains / len(successful)
            self.stats['average_domains_per_decomposition'] = avg_domains
        
        self.logger.info(f"Decomposition summary for {protein_id}: "
                        f"{len(successful)} successful, {len(failed)} failed")
        
        if self.config.log_failures_detail and failed:
            failure_reasons = {}
            for result in failed:
                reason = result.status.value
                failure_reasons[reason] = failure_reasons.get(reason, 0) + 1
            
            self.logger.info(f"Failure breakdown: {failure_reasons}")
    
    def get_service_statistics(self) -> Dict[str, Any]:
        """Get service performance statistics"""
        
        total = self.stats['total_attempts']
        success_rate = (self.stats['successful_decompositions'] / total * 100) if total > 0 else 0
        
        return {
            'total_attempts': total,
            'success_rate_percent': success_rate,
            'successful_decompositions': self.stats['successful_decompositions'],
            'failure_breakdown': {
                'short_domains': self.stats['failed_short_domains'],
                'poor_coverage': self.stats['failed_poor_coverage'],
                'no_architecture': self.stats['failed_no_architecture'],
                'alignment_issues': self.stats['failed_alignment_issues']
            },
            'average_domains_per_success': self.stats['average_domains_per_decomposition']
        }
    
    def reset_statistics(self):
        """Reset service statistics"""
        for key in self.stats:
            if isinstance(self.stats[key], (int, float)):
                self.stats[key] = 0
        
        self.logger.info("Service statistics reset")
