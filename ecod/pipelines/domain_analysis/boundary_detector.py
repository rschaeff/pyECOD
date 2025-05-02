# ecod/pipelines/domain_analysis/boundary_detector.py

class BoundaryDetector:
    """Detects domain boundaries from evidence"""
    
    def __init__(self, context=None):
        self.context = context or ApplicationContext()
        self.logger = logging.getLogger("ecod.boundary_detector")
        
    def detect_boundaries(self, summary: DomainSummaryModel) -> List[DomainCandidate]:
        """Main method to detect domain boundaries from summary data
        
        Replaces _determine_domain_boundaries in partition.py
        """
        domain_candidates = []
        
        # Get candidates from HHSearch hits (high confidence)
        hhsearch_candidates = self._analyze_hhsearch_hits(summary.hhsearch_hits)
        if hhsearch_candidates:
            self.logger.info(f"Found {len(hhsearch_candidates)} candidates from HHSearch")
            domain_candidates.extend(hhsearch_candidates)
            
        # Get candidates from chain BLAST hits
        chain_blast_candidates = self._analyze_chain_blast_hits(summary.chain_blast_hits)
        if chain_blast_candidates:
            self.logger.info(f"Found {len(chain_blast_candidates)} candidates from chain BLAST")
            domain_candidates.extend(chain_blast_candidates)
            
        # Get candidates from domain BLAST hits
        domain_blast_candidates = self._analyze_domain_blast_hits(summary.domain_blast_hits)
        if domain_blast_candidates:
            self.logger.info(f"Found {len(domain_blast_candidates)} candidates from domain BLAST")
            domain_candidates.extend(domain_blast_candidates)
            
        # Resolve overlapping boundaries
        final_candidates = self._resolve_overlaps(domain_candidates, summary.sequence_length)
        
        return final_candidates
    
    def _analyze_hhsearch_hits(self, hits: List[HHSearchHit]) -> List[DomainCandidate]:
        # Implementation moved from _identify_domains_from_hhsearch
        pass
    
    def _analyze_chain_blast_hits(self, hits: List[BlastHit]) -> List[DomainCandidate]:
        # Implementation moved from _analyze_chainwise_hits_for_domains
        pass
    
    def _analyze_domain_blast_hits(self, hits: List[BlastHit]) -> List[DomainCandidate]:
        # Implementation moved from _analyze_domain_blast_hits
        pass
    
    def _resolve_overlaps(self, candidates: List[DomainCandidate], 
                         sequence_length: int) -> List[DomainCandidate]:
        # Implementation moved from _resolve_domain_boundaries
        pass
