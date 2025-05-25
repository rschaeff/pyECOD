class EvidenceAnalyzer:
    """Evidence validation and analysis"""
    
    def validate_evidence(self, evidence: Evidence, context: str) -> ValidationResult:
        """Validate evidence with detailed error reporting"""
        
    def extract_evidence_from_summary(self, summary_data: Dict) -> List[Evidence]:
        """Extract and validate evidence from summary"""
        
    def map_chain_evidence_to_domains(self, chain_evidence: Evidence, 
                                     ref_domains: List[Dict]) -> List[Evidence]:
        """Map chain-level evidence to domain-level"""
