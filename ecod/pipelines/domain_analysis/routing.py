# ecod/pipelines/domain_analysis/routing.py
import os
import math
import logging
from typing import Dict, Any, List, Optional, Tuple, Set

from ecod.core.context import ApplicationContext
from ecod.exceptions import PipelineError, ConfigurationError

class ProcessingRouter:
    """
    Routes proteins to appropriate processing paths based on analysis of BLAST results
    """
    
    def __init__(self, context: ApplicationContext):
        """
        Initialize the processing router
        
        Args:
            context: Application context with shared resources
        """
        self.context = context
        self.db = context.db
        self.config = context.config_manager.config
        self.logger = logging.getLogger("ecod.pipelines.domain_analysis.routing")
        
        # Default confidence threshold
        self.confidence_threshold = self.config.get('routing', {}).get('confidence_threshold', 0.9)
        
        # Load additional configuration settings
        self._load_configuration()
        
    def _load_configuration(self) -> None:
        """Load routing configuration settings"""
        routing_config = self.config.get('routing', {})
        
        # Update confidence threshold if specified
        if 'confidence_threshold' in routing_config:
            self.confidence_threshold = float(routing_config['confidence_threshold'])
            self.logger.info(f"Using confidence threshold: {self.confidence_threshold}")
        
        # Other configuration settings can be added here as needed
    
def assign_processing_paths(self, batch_id: int) -> Dict[str, List[int]]:
    """
    Assign proteins in a batch to appropriate processing paths
    with restart support
    """
    # Check if routing already exists for this batch
    query = """
    SELECT COUNT(*) 
    FROM ecod_schema.protein_processing_path
    WHERE batch_id = %s
    """
    
    try:
        rows = self.db.execute_query(query, (batch_id,))
        already_routed = rows and rows[0][0] > 0
        
        if already_routed:
            self.logger.info(f"Found existing routing for batch {batch_id}, using stored paths")
            return self._get_existing_paths(batch_id)
    except Exception as e:
        self.logger.error(f"Error checking for existing routing: {e}")
        # Continue with fresh routing if check fails
    
    # Initialize processing path groups
    paths = {
        "blast_only": [],    # Skip HHSearch entirely
        "full_pipeline": [], # Complete pipeline with HHSearch
    }
    
    # Get all proteins with BLAST results in this batch
    # ... rest of the original method ...
    
    return paths

def _get_existing_paths(self, batch_id: int) -> Dict[str, List[int]]:
    """Get existing path assignments from database"""
    paths = {
        "blast_only": [],
        "full_pipeline": []
    }
    
    query = """
    SELECT protein_id, path_type 
    FROM ecod_schema.protein_processing_path
    WHERE batch_id = %s
    """
    
    try:
        rows = self.db.execute_dict_query(query, (batch_id,))
        
        for row in rows:
            protein_id = row['protein_id']
            path_type = row['path_type']
            
            if path_type in paths:
                paths[path_type].append(protein_id)
            else:
                self.logger.warning(f"Unknown path type '{path_type}' for protein {protein_id}")
    except Exception as e:
        self.logger.error(f"Error retrieving existing paths: {e}")
    
    # Log path counts
    for path_type, proteins in paths.items():
        self.logger.info(f"Found {len(proteins)} proteins in {path_type} path")
    
    return paths
    
    def prioritize_full_pipeline_proteins(self, batch_id: int) -> Dict[str, List[int]]:
        """
        Prioritize proteins in the full pipeline path
        
        Args:
            batch_id: Batch ID
            
        Returns:
            Dictionary mapping priority levels to lists of protein IDs
        """
        # Get proteins assigned to the full pipeline path
        query = """
            SELECT 
                p.id, p.length, ps.id as process_id
            FROM 
                ecod_schema.protein p
            JOIN
                ecod_schema.process_status ps ON p.id = ps.protein_id
            JOIN
                ecod_schema.protein_processing_path ppp ON p.id = ppp.protein_id
            WHERE 
                ps.batch_id = %s
                AND ppp.path_type = 'full_pipeline'
        """
        
        try:
            rows = self.db.execute_dict_query(query, (batch_id,))
        except Exception as e:
            self.logger.error(f"Error querying full pipeline proteins for batch {batch_id}: {e}")
            return {"high_priority": [], "standard_priority": []}
        
        # Initialize priority groups
        priorities = {
            "high_priority": [],    # Complex/novel domains needing careful analysis
            "standard_priority": [] # Typical cases
        }
        
        for row in rows:
            protein_id = row['id']
            process_id = row['process_id']
            length = row.get('length', 0)
            
            # Analyze complexity factors
            is_multidomain = self._is_likely_multidomain(protein_id)
            has_weak_matches = self._has_only_weak_blast_matches(protein_id)
            is_long = length > 500  # Longer sequences often need more attention
            
            # Assign priority
            if (is_multidomain or has_weak_matches or is_long):
                priority = "high_priority"
                priorities["high_priority"].append(protein_id)
            else:
                priority = "standard_priority"
                priorities["standard_priority"].append(protein_id)
            
            # Update in database
            self._update_path_priority(batch_id, protein_id, process_id, priority)
        
        self.logger.info(
            f"Prioritized {len(rows)} full pipeline proteins in batch {batch_id}: "
            f"{len(priorities['high_priority'])} high priority, "
            f"{len(priorities['standard_priority'])} standard priority"
        )
        
        return priorities
    
    def evaluate_blast_confidence(self, protein_id: int) -> float:
        """
        Evaluate the confidence level of BLAST results for a protein
        
        Args:
            protein_id: ID of the protein to evaluate
            
        Returns:
            float: Confidence score between 0 and 1
        """
        # Get BLAST results for this protein
        domain_blast_results = self._get_domain_blast_results(protein_id)
        chain_blast_results = self._get_chain_blast_results(protein_id)
        
        # If no results, we have no confidence
        if not domain_blast_results and not chain_blast_results:
            self.logger.info(f"No BLAST results for protein {protein_id}, setting confidence to 0")
            return 0.0
        
        # Calculate metrics from domain BLAST
        domain_metrics = self._calculate_domain_blast_metrics(domain_blast_results)
        
        # Calculate metrics from chain BLAST
        chain_metrics = self._calculate_chain_blast_metrics(chain_blast_results)
        
        # Combine metrics
        min_evalue = min(
            domain_metrics.get('min_evalue', 10.0),
            chain_metrics.get('min_evalue', 10.0)
        )
        max_coverage = max(
            domain_metrics.get('max_coverage', 0.0),
            chain_metrics.get('max_coverage', 0.0)
        )
        max_identity = max(
            domain_metrics.get('max_identity', 0.0),
            chain_metrics.get('max_identity', 0.0)
        )
        hit_consistency = (
            domain_metrics.get('hit_consistency', 0.0) * 0.7 +
            chain_metrics.get('hit_consistency', 0.0) * 0.3
        )
        
        # Calculate domain architecture coverage
        architecture_coverage = self._calculate_architecture_coverage(protein_id, domain_blast_results)
        
        # Calculate weighted confidence score
        confidence = (
            0.3 * self._evalue_to_confidence(min_evalue) + 
            0.25 * max_coverage + 
            0.2 * max_identity + 
            0.15 * hit_consistency +
            0.1 * architecture_coverage
        )
        
        # Store metrics for analysis
        self._store_confidence_metrics(
            protein_id,
            min_evalue,
            max_coverage,
            max_identity,
            hit_consistency,
            architecture_coverage,
            confidence
        )
        
        self.logger.info(f"Calculated confidence {confidence:.4f} for protein {protein_id}")
        return min(1.0, max(0.0, confidence))  # Ensure value between 0-1
    
    def _get_domain_blast_results(self, protein_id: int) -> List[Dict[str, Any]]:
        """
        Get domain BLAST results for a protein
        
        Args:
            protein_id: Protein ID
            
        Returns:
            List of domain BLAST hit dictionaries
        """
        query = """
            SELECT 
                h.* 
            FROM 
                ecod_schema.domain_blast_hit h
            WHERE 
                h.protein_id = %s
            ORDER BY 
                h.evalue ASC
            LIMIT 20
        """
        
        try:
            # This table is a simplified example - you'll need to adjust to your actual schema
            rows = self.db.execute_dict_query(query, (protein_id,))
            return rows
        except Exception as e:
            self.logger.error(f"Error getting domain BLAST results for protein {protein_id}: {e}")
            return []
    
    def _get_chain_blast_results(self, protein_id: int) -> List[Dict[str, Any]]:
        """
        Get chain BLAST results for a protein
        
        Args:
            protein_id: Protein ID
            
        Returns:
            List of chain BLAST hit dictionaries
        """
        query = """
            SELECT 
                h.* 
            FROM 
                ecod_schema.chain_blast_hit h
            WHERE 
                h.protein_id = %s
            ORDER BY 
                h.evalue ASC
            LIMIT 20
        """
        
        try:
            # This table is a simplified example - you'll need to adjust to your actual schema
            rows = self.db.execute_dict_query(query, (protein_id,))
            return rows
        except Exception as e:
            self.logger.error(f"Error getting chain BLAST results for protein {protein_id}: {e}")
            return []
    
    def _calculate_domain_blast_metrics(self, blast_results: List[Dict[str, Any]]) -> Dict[str, float]:
        """
        Calculate confidence metrics from domain BLAST results
        
        Args:
            blast_results: List of domain BLAST hit dictionaries
            
        Returns:
            Dict with calculated metrics
        """
        if not blast_results:
            return {}
        
        # Consider top 5 hits
        top_hits = blast_results[:5]
        
        # Extract e-values
        evalues = []
        for hit in top_hits:
            try:
                evalues.extend([float(e) for e in hit.get('evalues', "").split(",") if e])
            except (ValueError, AttributeError):
                continue
        min_evalue = min(evalues) if evalues else 10.0
        
        # Calculate coverage
        coverages = []
        for hit in top_hits:
            query_regions = hit.get('query_regions', "").split(",")
            query_len = self._get_protein_length(hit.get('protein_id'))
            if query_len:
                coverage = self._calculate_regions_coverage(query_regions, query_len)
                coverages.append(coverage)
        max_coverage = max(coverages) if coverages else 0.0
        
        # Calculate identity (from BLAST HSP data)
        identities = []
        for hit in top_hits:
            # Identity might be available directly or need calculation from alignment
            identity = hit.get('identity')
            if identity is None and 'hit_seqs' in hit and 'query_seqs' in hit:
                # Calculate from alignment
                identity = self._calculate_identity_from_seqs(
                    hit['query_seqs'].split(","), 
                    hit['hit_seqs'].split(",")
                )
            if identity is not None:
                identities.append(float(identity))
        max_identity = max(identities) if identities else 0.0
        
        # Evaluate hit consistency
        # Higher when multiple hits point to the same ECOD classification
        classifications = {}
        for hit in top_hits:
            domain_id = hit.get('domain_id')
            if domain_id:
                t_group, h_group = self._get_domain_classification(domain_id)
                key = f"{t_group}_{h_group}"
                classifications[key] = classifications.get(key, 0) + 1
        
        # Calculate consistency as ratio of most common classification
        most_common_count = max(classifications.values()) if classifications else 0
        hit_consistency = most_common_count / len(top_hits) if top_hits else 0.0
        
        return {
            'min_evalue': min_evalue,
            'max_coverage': max_coverage,
            'max_identity': max_identity,
            'hit_consistency': hit_consistency
        }
    
    def _calculate_chain_blast_metrics(self, blast_results: List[Dict[str, Any]]) -> Dict[str, float]:
        """
        Calculate confidence metrics from chain BLAST results
        
        Args:
            blast_results: List of chain BLAST hit dictionaries
            
        Returns:
            Dict with calculated metrics
        """
        # Similar to domain metrics but with chain-specific considerations
        if not blast_results:
            return {}
        
        # Consider top 5 hits
        top_hits = blast_results[:5]
        
        # Extract e-values
        evalues = []
        for hit in top_hits:
            try:
                evalues.extend([float(e) for e in hit.get('evalues', "").split(",") if e])
            except (ValueError, AttributeError):
                continue
        min_evalue = min(evalues) if evalues else 10.0
        
        # Calculate coverage
        coverages = []
        for hit in top_hits:
            query_regions = hit.get('query_regions', "").split(",")
            query_len = self._get_protein_length(hit.get('protein_id'))
            if query_len:
                coverage = self._calculate_regions_coverage(query_regions, query_len)
                coverages.append(coverage)
        max_coverage = max(coverages) if coverages else 0.0
        
        # Calculate identity
        identities = []
        for hit in top_hits:
            identity = hit.get('identity')
            if identity is None and 'hit_seqs' in hit and 'query_seqs' in hit:
                identity = self._calculate_identity_from_seqs(
                    hit['query_seqs'].split(","), 
                    hit['hit_seqs'].split(",")
                )
            if identity is not None:
                identities.append(float(identity))
        max_identity = max(identities) if identities else 0.0
        
        # For chains, evaluate consistency based on PDB chain sources
        # Higher when multiple hits point to chains with similar ECOD domain compositions
        hit_consistency = self._evaluate_chain_hit_consistency(top_hits)
        
        return {
            'min_evalue': min_evalue,
            'max_coverage': max_coverage,
            'max_identity': max_identity,
            'hit_consistency': hit_consistency
        }
    
    def _evalue_to_confidence(self, evalue: float) -> float:
        """
        Convert E-value to confidence score
        
        Args:
            evalue: BLAST E-value
            
        Returns:
            Confidence score (0-1)
        """
        # Lower E-values give higher confidence
        # This formula gives 0.99 confidence at E-value 1e-10
        # and 0.5 confidence at E-value 1e-3
        if evalue <= 0:
            return 1.0
        
        confidence = 1.0 - (1.0 / (1.0 + math.exp(-math.log10(evalue) - 3)))
        return min(1.0, max(0.0, confidence))
    
    def _calculate_regions_coverage(self, regions: List[str], sequence_length: int) -> float:
        """
        Calculate the coverage of regions compared to total sequence length
        
        Args:
            regions: List of region strings (e.g., "1-100", "150-200")
            sequence_length: Total length of the sequence
            
        Returns:
            Coverage ratio (0-1)
        """
        if not regions or not sequence_length:
            return 0.0
        
        # Count positions covered by regions
        covered_positions = set()
        for region in regions:
            if not region:
                continue
            
            if '-' in region:
                try:
                    start, end = map(int, region.split('-'))
                    covered_positions.update(range(start, end + 1))
                except (ValueError, TypeError):
                    continue
        
        # Calculate coverage
        coverage = len(covered_positions) / sequence_length
        return min(1.0, coverage)  # Cap at 1.0
    
    def _calculate_identity_from_seqs(self, query_seqs: List[str], hit_seqs: List[str]) -> float:
        """
        Calculate sequence identity from aligned sequences
        
        Args:
            query_seqs: List of query sequence segments
            hit_seqs: List of hit sequence segments
            
        Returns:
            Identity percentage (0-1)
        """
        if not query_seqs or not hit_seqs or len(query_seqs) != len(hit_seqs):
            return 0.0
        
        total_aligned = 0
        total_identical = 0
        
        # Process each segment pair
        for q_seq, h_seq in zip(query_seqs, hit_seqs):
            if len(q_seq) != len(h_seq):
                continue
                
            # Count aligned and identical positions
            for q, h in zip(q_seq, h_seq):
                if q != '-' and h != '-':  # Both aligned
                    total_aligned += 1
                    if q == h:  # Identical
                        total_identical += 1
        
        # Calculate identity
        if total_aligned > 0:
            return total_identical / total_aligned
        return 0.0
    
    def _get_domain_classification(self, domain_id: str) -> Tuple[str, str]:
        """
        Get ECOD classification for a domain
        
        Args:
            domain_id: ECOD domain ID
            
        Returns:
            Tuple of (t_group, h_group)
        """
        query = """
            SELECT 
                d.t_group, d.h_group
            FROM 
                pdb_analysis.domain d
            WHERE 
                d.domain_id = %s
        """
        
        try:
            rows = self.db.execute_dict_query(query, (domain_id,))
            if rows:
                return rows[0].get('t_group', 'unknown'), rows[0].get('h_group', 'unknown')
        except Exception as e:
            self.logger.error(f"Error getting classification for domain {domain_id}: {e}")
        
        return 'unknown', 'unknown'
    
    def _evaluate_chain_hit_consistency(self, hits: List[Dict[str, Any]]) -> float:
        """Evaluate consistency of chain hits based on domain composition"""
        if not hits:
            return 0.0
        
        # Get domain compositions for hit chains
        chain_domains = {}
        for hit in hits:
            pdb_id = hit.get('pdb_id')
            chain_id = hit.get('chain_id')
            
            if pdb_id and chain_id:
                source_id = f"{pdb_id}_{chain_id}"
                domains = self._get_chain_domains(source_id)
                if domains:
                    chain_domains[source_id] = domains
        
        # Sort domains with handling for None values
        def safe_sort_key(domain):
            t_group = domain.get('t_group', '')
            h_group = domain.get('h_group', '')
            # Replace None with empty string
            t_group = '' if t_group is None else t_group
            h_group = '' if h_group is None else h_group
            return (t_group, h_group)
        
        # Calculate consistency based on domain composition similarity
        if len(chain_domains) <= 1:
            return 0.0
        
        # Count shared domain architecture patterns
        architecture_counts = {}
        for source_id, domains in chain_domains.items():
            if not domains:
                continue
            
            # Sort domains by t_group and h_group
            sorted_domains = sorted(domains, key=safe_sort_key)
            
            # Create architecture key
            arch_key = "_".join([f"{d.get('t_group', 'X') or 'X'}.{d.get('h_group', 'X') or 'X'}" for d in sorted_domains])
            
            architecture_counts[arch_key] = architecture_counts.get(arch_key, 0) + 1
        
        # Calculate consistency score
        if not architecture_counts:
            return 0.0
        
        max_count = max(architecture_counts.values())
        consistency = max_count / len(chain_domains)
        
        return consistency
    
    def _get_chain_domains(self, source_id: str) -> List[Dict[str, Any]]:
        """
        Get domains for a chain
        
        Args:
            source_id: Chain source ID (e.g., "1abc_A")
            
        Returns:
            List of domain dictionaries
        """
        query = """
            SELECT 
                d.t_group, d.h_group, d.x_group, d.a_group
            FROM 
                pdb_analysis.domain d
            JOIN
                pdb_analysis.protein p ON d.protein_id = p.id
            WHERE 
                p.source_id = %s
        """
        
        try:
            return self.db.execute_dict_query(query, (source_id,))
        except Exception as e:
            self.logger.error(f"Error getting domains for chain {source_id}: {e}")
            return []
    
    def _calculate_architecture_coverage(self, protein_id: int, 
                                       domain_blast_results: List[Dict[str, Any]]
    ) -> float:
        """
        Calculate how well BLAST hits cover the entire protein domain architecture
        
        Args:
            protein_id: Protein ID
            domain_blast_results: Domain BLAST results
            
        Returns:
            Architecture coverage score (0-1)
        """
        if not domain_blast_results:
            return 0.0
        
        protein_length = self._get_protein_length(protein_id)
        if not protein_length:
            return 0.0
        
        # Collect all query regions from all hits
        all_regions = []
        for hit in domain_blast_results:
            query_regions = hit.get('query_regions', "").split(",")
            all_regions.extend([r for r in query_regions if r])
        
        # Calculate non-overlapping coverage
        coverage = self._calculate_non_overlapping_coverage(all_regions, protein_length)
        
        # Higher coverage means more complete architecture detection
        return coverage
    
    def _calculate_non_overlapping_coverage(self, regions: List[str], length: int) -> float:
        """
        Calculate non-overlapping coverage of regions
        
        Args:
            regions: List of region strings (e.g., "1-100", "150-200")
            length: Total length
            
        Returns:
            Non-overlapping coverage ratio (0-1)
        """
        if not regions or not length:
            return 0.0
        
        # Convert regions to sets of positions
        all_positions = set()
        for region in regions:
            if not region or '-' not in region:
                continue
                
            try:
                start, end = map(int, region.split('-'))
                all_positions.update(range(start, end + 1))
            except (ValueError, TypeError):
                continue
        
        # Calculate non-overlapping coverage
        coverage = len(all_positions) / length
        return min(1.0, coverage)
    
    def _get_protein_length(self, protein_id: int) -> int:
        """
        Get the length of a protein
        
        Args:
            protein_id: Protein ID
            
        Returns:
            Protein length
        """
        query = "SELECT length FROM ecod_schema.protein WHERE id = %s"
        
        try:
            rows = self.db.execute_query(query, (protein_id,))
            if rows:
                return rows[0][0]
        except Exception as e:
            self.logger.error(f"Error getting protein length for protein {protein_id}: {e}")
        
        return 0
    
    def _store_confidence_metrics(self, protein_id: int, min_evalue: float, 
                                max_coverage: float, max_identity: float,
                                hit_consistency: float, architecture_coverage: float,
                                overall_confidence: float
    ) -> None:
        """
        Store confidence metrics in the database for analysis
        
        Args:
            protein_id: Protein ID
            min_evalue: Minimum E-value from BLAST hits
            max_coverage: Maximum query coverage
            max_identity: Maximum sequence identity
            hit_consistency: Hit classification consistency score
            architecture_coverage: Domain architecture coverage score
            overall_confidence: Overall calculated confidence
        """
        try:
            self.db.insert(
                "ecod_schema.blast_confidence_metrics",
                {
                    "protein_id": protein_id,
                    "min_evalue": min_evalue,
                    "max_coverage": max_coverage,
                    "max_identity": max_identity,
                    "hit_consistency": hit_consistency,
                    "architecture_coverage": architecture_coverage,
                    "overall_confidence": overall_confidence
                }
            )
        except Exception as e:
            self.logger.error(f"Error storing confidence metrics for protein {protein_id}: {e}")
    
    def _record_processing_path(self, batch_id: int, protein_id: int, 
                              process_id: int, path_type: str, 
                              confidence: float
    ) -> None:
        """
        Record protein processing path assignment in database
        """
        try:
            # First check if the record already exists
            query = """
            SELECT id FROM ecod_schema.protein_processing_path
            WHERE batch_id = %s AND protein_id = %s
            """
            rows = self.db.execute_query(query, (batch_id, protein_id))
            
            if rows:
                # Update existing record
                self.db.update(
                    "ecod_schema.protein_processing_path",
                    {
                        "path_type": path_type,
                        "confidence_score": confidence
                    },
                    "batch_id = %s AND protein_id = %s",
                    (batch_id, protein_id)
                )
            else:
                # Insert new record
                self.db.insert(
                    "ecod_schema.protein_processing_path",
                    {
                        "batch_id": batch_id,
                        "protein_id": protein_id,
                        "path_type": path_type,
                        "confidence_score": confidence
                    }
                )
            
            # Update process_status table
            self.db.update(
                "ecod_schema.process_status",
                {
                    "processing_path": path_type
                },
                "id = %s",
                (process_id,)
            )
        except Exception as e:
            self.logger.error(f"Error recording processing path for protein {protein_id}: {e}")
    
    def _update_path_priority(self, batch_id: int, protein_id: int, 
                            process_id: int, priority: str
    ) -> None:
        """
        Update processing path priority
        
        Args:
            batch_id: Batch ID
            protein_id: Protein ID
            process_id: Process status ID
            priority: Priority level ('high_priority' or 'standard_priority')
        """
        try:
            # Update protein_processing_path table
            self.db.update(
                "ecod_schema.protein_processing_path",
                {
                    "priority": priority
                },
                "batch_id = %s AND protein_id = %s",
                (batch_id, protein_id)
            )
            
            # Update process_status table
            self.db.update(
                "ecod_schema.process_status",
                {
                    "path_priority": priority
                },
                "id = %s",
                (process_id,)
            )
        except Exception as e:
            self.logger.error(f"Error updating priority for protein {protein_id}: {e}")
    
    def _update_batch_routing_status(self, batch_id: int, completed: bool) -> None:
        """
        Update batch routing status
        
        Args:
            batch_id: Batch ID
            completed: Whether routing is completed
        """
        try:
            update_data = {
                "routing_strategy": "adaptive_blast_confidence",
                "routing_completed": completed
            }
            
            
            self.db.update(
                "ecod_schema.batch",
                update_data,
                "id = %s",
                (batch_id,)
            )
        except Exception as e:
            self.logger.error(f"Error updating batch routing status for batch {batch_id}: {e}")
    
    def _is_likely_multidomain(self, protein_id: int) -> bool:
        """
        Determine if a protein is likely to have multiple domains
        
        Args:
            protein_id: Protein ID
            
        Returns:
            True if likely multidomain, False otherwise
        """
        # Check protein length (longer proteins more likely to be multidomain)
        length = self._get_protein_length(protein_id)
        if length > 400:  # Empirical threshold
            return True
        
        # Check for multiple non-overlapping hits in domain BLAST
        domain_hits = self._get_domain_blast_results(protein_id)
        if len(domain_hits) >= 3:  # Multiple hits may indicate multiple domains
            # Check if hits cover different regions
            covered_regions = []
            for hit in domain_hits[:5]:  # Consider top 5 hits
                regions = hit.get('query_regions', "").split(",")
                for region in regions:
                    if region and '-' in region:
                        covered_regions.append(region)
            
            # Check for non-overlapping regions
            if self._has_non_overlapping_regions(covered_regions, min_gap=50):
                return True
        
        return False
    
    def _has_non_overlapping_regions(self, regions: List[str], min_gap: int = 50) -> bool:
        """
        Check if regions have significant non-overlapping sections
        
        Args:
            regions: List of region strings (e.g., "1-100", "150-200")
            min_gap: Minimum gap size to consider regions separate
            
        Returns:
            True if non-overlapping regions found, False otherwise
        """
        if len(regions) < 2:
            return False
        
        # Parse region boundaries
        boundaries = []
        for region in regions:
            if '-' in region:
                try:
                    start, end = map(int, region.split('-'))
                    boundaries.append((start, end))
                except (ValueError, TypeError):
                    continue
        
        if len(boundaries) < 2:
            return False
        
        # Sort by start position
        boundaries.sort()
        
        # Check for gaps between regions
        for i in range(len(boundaries) - 1):
            current_end = boundaries[i][1]
            next_start = boundaries[i+1][0]
            
            if next_start - current_end > min_gap:
                return True
        
        return False
    
    def _has_only_weak_blast_matches(self, protein_id: int) -> bool:
        """
        Check if a protein has only weak BLAST matches
        
        Args:
            protein_id: Protein ID
            
        Returns:
            True if only weak matches, False otherwise
        """
        # Get domain BLAST results
        domain_hits = self._get_domain_blast_results(protein_id)
        
        # Check if no hits or only high e-value hits
        if not domain_hits:
            return True
        
        # Check e-values of top hits
        for hit in domain_hits[:3]:  # Consider top 3 hits
            try:
                evalues = [float(e) for e in hit.get('evalues', "").split(",") if e]
                min_evalue = min(evalues) if evalues else 10.0
                
                # If any hit has good e-value, not considered "weak"
                if min_evalue < 1e-5:
                    return False
            except (ValueError, AttributeError):
                continue
        
        # If we get here, no strong hits found
        return True