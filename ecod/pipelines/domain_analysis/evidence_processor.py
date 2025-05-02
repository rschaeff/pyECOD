# ecod/pipelines/domain_analysis/evidence_processor.py

class EvidenceProcessor:
    """Processes raw hits into standardized evidence"""
    
    def __init__(self, context=None):
        self.context = context or ApplicationContext()
        self.logger = logging.getLogger("ecod.evidence_processor")
        
    def process_domain_summary(self, summary_path: str) -> DomainSummaryModel:
        """Process domain summary file
        
        Replaces _process_domain_summary in partition.py
        """
        # Implementation using proper model instantiation
        pass
    
    def convert_blast_to_evidence(self, hit: BlastHit) -> Evidence:
        """Convert BLAST hit to standardized evidence
        
        Used by boundary detection to create consistent evidence objects
        """
        evidence = Evidence(
            type=hit.hit_type or "blast",
            source_id=hit.domain_id or "",
            query_range=hit.range,
            hit_range=hit.hit_range,
            confidence=self._calculate_blast_confidence(hit)
        )
        
        # Add source-specific attributes
        evidence.attributes["evalue"] = hit.evalue
        evidence.attributes["hit_id"] = hit.hit_id
        evidence.attributes["pdb_id"] = hit.pdb_id
        evidence.attributes["chain_id"] = hit.chain_id
        
        return evidence
    
    def convert_hhsearch_to_evidence(self, hit: HHSearchHit) -> Evidence:
        """Convert HHSearch hit to standardized evidence"""
        evidence = Evidence(
            type="hhsearch",
            source_id=hit.domain_id or hit.hit_id,
            query_range=hit.range,
            hit_range=hit.hit_range,
            confidence=hit.probability / 100.0  # Normalize to 0-1 range
        )
        
        # Add source-specific attributes
        evidence.attributes["probability"] = hit.probability
        evidence.attributes["evalue"] = hit.evalue
        evidence.attributes["score"] = hit.score
        
        return evidence
    
    def _calculate_blast_confidence(self, hit: BlastHit) -> float:
        """Calculate confidence score for BLAST hit
        
        Returns value between 0 and 1
        """
        # Convert e-value to confidence score
        if hit.evalue >= 10:
            return 0.0
        elif hit.evalue <= 1e-20:
            return 1.0
        else:
            # Log scale conversion
            import math
            confidence = 1.0 - (math.log10(hit.evalue) + 20) / 20
            return max(0.0, min(1.0, confidence))
