# ecod/pipelines/domain_analysis/classifier.py

class DomainClassifier:
    """Assigns ECOD classifications to domain candidates"""
    
    def __init__(self, context=None):
        self.context = context or ApplicationContext()
        self.logger = logging.getLogger("ecod.domain_classifier")
        
        # Initialize classification caches
        self.domain_classification_cache = {}
        
    def classify_domains(self, candidates: List[DomainCandidate]) -> List[Domain]:
        """Assign classifications to domain candidates
        
        Replaces _assign_domain_classifications in partition.py
        """
        domains = []
        
        for candidate in candidates:
            # Create domain from candidate
            domain = Domain(
                domain_id=candidate.id or f"domain_{len(domains)+1}",
                range=candidate.range,
                t_group=candidate.t_group,
                h_group=candidate.h_group,
                x_group=candidate.x_group,
                a_group=candidate.a_group
            )
            
            # If classification is missing, try to infer from evidence
            if not all([domain.t_group, domain.h_group]):
                self._infer_classification(domain, candidate.evidence)
                
            domains.append(domain)
            
        return domains
    
    def _infer_classification(self, domain: Domain, evidence: List[Evidence]) -> None:
        """Infer classification from evidence
        
        Updates domain in place
        """
        # Find best evidence for classification
        best_evidence = None
        best_confidence = 0.0
        
        for e in evidence:
            if e.confidence > best_confidence and e.source_id:
                best_evidence = e
                best_confidence = e.confidence
                
        if best_evidence and best_evidence.source_id:
            # Get classification from database
            classification = self._get_classification(best_evidence.source_id)
            if classification:
                domain.t_group = classification.get("t_group") or domain.t_group
                domain.h_group = classification.get("h_group") or domain.h_group
                domain.x_group = classification.get("x_group") or domain.x_group
                domain.a_group = classification.get("a_group") or domain.a_group
    
    def _get_classification(self, domain_id: str) -> Optional[Dict[str, Any]]:
        """Get classification for domain ID with caching
        
        Replaces _get_domain_classification_by_id in partition.py
        """
        # Check cache first
        if domain_id in self.domain_classification_cache:
            return self.domain_classification_cache[domain_id]
            
        # Query database
        query = """
        SELECT 
            d.t_group, d.h_group, d.x_group, d.a_group,
            d.is_manual_rep, d.is_f70, d.is_f40, d.is_f99
        FROM 
            pdb_analysis.domain d
        WHERE 
            d.domain_id = %s
        """
        
        try:
            rows = self.context.db.execute_dict_query(query, (domain_id,))
            if rows:
                classification = {
                    "t_group": rows[0].get("t_group"),
                    "h_group": rows[0].get("h_group"),
                    "x_group": rows[0].get("x_group"),
                    "a_group": rows[0].get("a_group"),
                    "is_manual_rep": False,  # New domains are not manual reps
                    "is_f70": False,
                    "is_f40": False,
                    "is_f99": False
                }
                
                # Cache the result
                self.domain_classification_cache[domain_id] = classification
                return classification
        except Exception as e:
            self.logger.error(f"Error getting classification for {domain_id}: {e}")
            
        return None
