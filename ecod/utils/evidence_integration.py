from typing import List, Dict, Any, Optional, Tuple, Union
import xml.etree.ElementTree as ET
import logging

from ecod.models.evidence import DomainEvidence
from ecod.models.domain_analysis.domain_model import DomainModel

def identify_domains_from_hhsearch(pdb_id: str, chain_id: str, hhsearch_hits: List[Any]) -> List[Dict[str, Any]]:
    """
    Identify domains from HHSearch hits using DomainEvidence model
    
    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        hhsearch_hits: List of HHSearch hits from summary
        
    Returns:
        List of domain dictionaries from HHSearch
    """
    logger = logging.getLogger(__name__)
    
    # Filter hits by probability threshold (only high confidence hits)
    high_prob_hits = [hit for hit in hhsearch_hits
                     if hasattr(hit, 'probability') and hit.probability >= 90.0]
    
    if not high_prob_hits:
        return []
    
    # Create domain candidates from high probability hits
    domains = []
    for hit in high_prob_hits:
        # Parse range
        ranges = []
        if hasattr(hit, 'range') and hit.range:
            ranges = parse_range(hit.range)
        else:
            continue
        
        if not ranges:
            continue
            
        # Create evidence using DomainEvidence model
        evidence = DomainEvidence.from_hhsearch_hit(hit)
        
        # Create domain for each range segment
        for start, end in ranges:
            domain = {
                'start': start,
                'end': end,
                'range': f"{start}-{end}",
                'source': 'hhsearch',
                'confidence': hit.probability / 100.0 if hasattr(hit, 'probability') else 0.0,
                'source_id': hit.domain_id if hasattr(hit, 'domain_id') else
                             (hit.hit_id if hasattr(hit, 'hit_id') else ''),
                't_group': None,  # Will be set during classification
                'h_group': None,
                'x_group': None,
                'a_group': None,
                'evidence': [evidence]  # Store DomainEvidence object in list
            }
            domains.append(domain)
    
    return domains

def identify_domains_from_chain_blast(pdb_id: str, chain_id: str, chain_blast_hits: List[Any], 
                                    get_reference_domains_func, map_domains_to_query_func) -> List[Dict[str, Any]]:
    """
    Identify domains from chain BLAST hits using DomainEvidence model
    
    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        chain_blast_hits: List of chain BLAST hits from summary
        get_reference_domains_func: Function to get reference domains
        map_domains_to_query_func: Function to map domains to query
        
    Returns:
        List of domain dictionaries from chain BLAST
    """
    logger = logging.getLogger(__name__)
    
    if not chain_blast_hits:
        logger.warning("No chain BLAST hits provided to analyze")
        return []
    
    # Process each hit to find domains
    domains = []
    for i, hit in enumerate(chain_blast_hits[:30]):  # Limit to first 30 hits for performance
        # Skip hits without required fields
        if not (hasattr(hit, 'pdb_id') and hasattr(hit, 'chain_id')):
            continue
        
        # Parse range
        ranges = []
        if hasattr(hit, 'range') and hit.range:
            ranges = parse_range(hit.range)
        elif hasattr(hit, 'range_parsed') and hit.range_parsed:
            ranges = hit.range_parsed
        else:
            continue
        
        if not ranges:
            continue
        
        # Get reference domains for this hit
        hit_domains = get_reference_domains_func(
            hit.pdb_id if hasattr(hit, 'pdb_id') else '',
            hit.chain_id if hasattr(hit, 'chain_id') else ''
        )
        if not hit_domains:
            continue
        
        # Map reference domains to query sequence
        mapped_domains = map_domains_to_query_func(ranges, hit_domains)
        
        if mapped_domains:
            # Create evidence using DomainEvidence model
            evidence = DomainEvidence.from_blast_hit(hit, "chain_blast")
            
            # Update evidence in mapped domains
            for domain in mapped_domains:
                domain['evidence'] = [evidence]  # Store DomainEvidence object in list
            
            logger.info(f"Mapped {len(mapped_domains)} domains from hit #{i+1}")
            domains.extend(mapped_domains)
    
    return domains

def identify_domains_from_domain_blast(pdb_id: str, chain_id: str, domain_blast_hits: List[Any], 
                                    get_domain_classification_func, parse_range_func) -> List[Dict[str, Any]]:
    """
    Identify domains from domain BLAST hits using DomainEvidence model
    
    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        domain_blast_hits: List of domain BLAST hits from summary
        get_domain_classification_func: Function to get domain classification
        parse_range_func: Function to parse range strings
        
    Returns:
        List of domain dictionaries from domain BLAST
    """
    logger = logging.getLogger(__name__)
    
    if not domain_blast_hits:
        return []
    
    # Group hits by query region
    region_hits = {}
    for hit in domain_blast_hits:
        # Parse range
        ranges = []
        if hasattr(hit, 'range') and hit.range:
            ranges = parse_range_func(hit.range)
        elif hasattr(hit, 'range_parsed') and hit.range_parsed:
            ranges = hit.range_parsed
        else:
            continue
        
        for start, end in ranges:
            region_key = f"{start}-{end}"
            if region_key not in region_hits:
                region_hits[region_key] = []
            region_hits[region_key].append(hit)
    
    # Create domain candidates from regions with sufficient hits
    domains = []
    for region, hits in region_hits.items():
        if len(hits) < 3:  # Require at least 3 hits for confidence
            continue
        
        start, end = map(int, region.split('-'))
        
        # Find most common classification among hits
        t_groups = {}
        for hit in hits:
            if hasattr(hit, 'domain_id') and hit.domain_id:
                domain_info = get_domain_classification_func(hit.domain_id)
                if domain_info and 't_group' in domain_info:
                    t_group = domain_info['t_group']
                    t_groups[t_group] = t_groups.get(t_group, 0) + 1
        
        # Use most common t_group if available
        t_group = None
        if t_groups:
            t_group = max(t_groups.items(), key=lambda x: x[1])[0]
        
        # Convert top hits to DomainEvidence objects
        evidence_list = []
        for hit in hits[:5]:  # Include top 5 hits as evidence
            evidence = DomainEvidence.from_blast_hit(hit, "domain_blast")
            evidence_list.append(evidence)
        
        domain = {
            'start': start,
            'end': end,
            'range': region,
            'source': 'domain_blast',
            'confidence': 0.7,  # Medium confidence level
            'source_id': '',  # Multiple hits, no single source
            't_group': t_group,
            'h_group': None,  # Will be set during classification
            'x_group': None,
            'a_group': None,
            'evidence': evidence_list  # Store list of DomainEvidence objects
        }
        domains.append(domain)
    
    return domains

def assign_domain_classifications(domains: List[Dict[str, Any]], reference_classifications: Dict[str, Dict[str, Any]],
                               get_domain_classification_func) -> Dict[str, int]:
    """
    Assign ECOD classifications to domains with enhanced evidence tracking
    
    Args:
        domains: List of domain dictionaries to classify
        reference_classifications: Dictionary of reference classifications by domain_id
        get_domain_classification_func: Function to get domain classification by ID
        
    Returns:
        Dictionary with classification statistics
    """
    logger = logging.getLogger(__name__)
    
    # Track classification statistics
    classification_stats = {
        "total_domains": len(domains),
        "with_evidence": 0,
        "with_reference": 0,
        "classified_from_evidence": 0,
        "classified_from_reference": 0,
        "unclassified": 0
    }
    
    # Assign classifications to domains
    for i, domain in enumerate(domains):
        domain_location = f"Domain {i+1} ({domain.get('start', '?')}-{domain.get('end', '?')})"
        logger.info(f"Assigning classification to {domain_location}")
        
        # Track source of classification
        classification_source = None
        
        # Case 1: If domain has reference, use it directly
        if "reference" in domain and domain["reference"]:
            domain_id = domain.get("domain_id", "")
            if domain_id in reference_classifications:
                domain.update(reference_classifications[domain_id])
                classification_source = {
                    "source": "reference",
                    "domain_id": domain_id
                }
                classification_stats["classified_from_reference"] += 1
                classification_stats["with_reference"] += 1
                logger.info(f"  Applied reference classification for {domain_location}")
                continue
        
        # Case 2: Check evidence for classification
        if "evidence" not in domain or not domain["evidence"]:
            logger.warning(f"  No evidence found for {domain_location}")
            classification_stats["unclassified"] += 1
            continue
        
        # Log evidence count
        evidence_count = len(domain["evidence"])
        classification_stats["with_evidence"] += 1
        logger.info(f"  Found {evidence_count} evidence items for {domain_location}")
        
        # Find the best evidence (highest probability/lowest e-value)
        best_evidence = None
        best_score = 0
        best_source = ""
        
        for j, evidence in enumerate(domain["evidence"]):
            # Ensure we're using DomainEvidence objects
            if not isinstance(evidence, DomainEvidence):
                try:
                    evidence = DomainEvidence.from_dict(evidence) if isinstance(evidence, dict) else evidence
                except Exception as e:
                    logger.warning(f"  Failed to convert evidence to DomainEvidence: {e}")
                    continue
            
            # Log each piece of evidence
            evidence_type = evidence.type
            domain_id = evidence.domain_id or evidence.source_id
            logger.debug(f"  Evidence {j+1}: type={evidence_type}, domain_id={domain_id}")
            
            if not domain_id:
                logger.debug(f"    Skipping evidence {j+1} - no valid domain_id")
                continue
            
            # Calculate score based on evidence type
            score = 0
            if evidence_type == "hhsearch":
                probability = evidence.attributes.get("probability", 0)
                score = probability
                best_source = f"HHSearch prob={probability:.1f}"
                logger.debug(f"    HHSearch evidence - probability: {probability}")
            elif evidence_type in ["blast", "chain_blast", "domain_blast"]:
                e_value = evidence.attributes.get("evalue", 999)
                score = 100.0 / (1.0 + e_value) if e_value < 10 else 0
                best_source = f"{evidence_type} e-value={e_value:.2e}"
                logger.debug(f"    BLAST evidence - evalue: {e_value}, score: {score}")
            else:
                logger.debug(f"    Unknown evidence type: {evidence_type}")
            
            logger.debug(f"    Evidence score for {domain_id}: {score}")
            
            if score > best_score:
                best_score = score
                best_evidence = evidence
                logger.debug(f"    New best evidence: {domain_id} with score {score}")
        
        if best_evidence:
            domain_id = best_evidence.domain_id or best_evidence.source_id
            logger.info(f"  Best evidence found: {domain_id} from {best_source}")
            
            # Get classification for this domain
            classification = get_domain_classification_func(domain_id)
            if classification:
                logger.info(f"  Classification for {domain_id}: {classification}")
                
                # Store classification details
                domain.update(classification)
                classification_stats["classified_from_evidence"] += 1
                
                # Save classification source information
                classification_source = {
                    "source": best_evidence.type,
                    "domain_id": domain_id,
                    "score": best_score,
                    "source_detail": best_source
                }
                
                # Copy classification to evidence for reference
                for cls_attr in ["t_group", "h_group", "x_group", "a_group"]:
                    if cls_attr in classification and classification[cls_attr]:
                        setattr(best_evidence, cls_attr, classification[cls_attr])
            else:
                logger.warning(f"  No classification found for domain_id {domain_id}")
                classification_stats["unclassified"] += 1
                
                # Add empty classification to prevent None values
                domain.update({
                    "t_group": "",
                    "h_group": "",
                    "x_group": "",
                    "a_group": ""
                })
        else:
            logger.warning(f"  No best evidence found for {domain_location}")
            classification_stats["unclassified"] += 1
            
            # Add empty classification to prevent None values
            domain.update({
                "t_group": "",
                "h_group": "",
                "x_group": "",
                "a_group": ""
            })
        
        # Store classification source if available
        if classification_source:
            domain["classification_source"] = classification_source
    
    return classification_stats

def convert_dict_domains_to_models(domains: List[Dict[str, Any]], pdb_id: str, chain_id: str) -> List[DomainModel]:
    """
    Convert dictionary domain representations to DomainModel objects

    Args:
        domains: List of domain dictionaries
        pdb_id: PDB identifier
        chain_id: Chain identifier

    Returns:
        List of DomainModel objects
    """
    from ecod.models.domain_analysis.domain_model import DomainModel
    from ecod.models.evidence import DomainEvidence
    import logging

    logger = logging.getLogger(__name__)

    domain_models = []
    for i, domain_dict in enumerate(domains):
        try:
            # Create domain ID if not present
            domain_id = domain_dict.get('id', f"{pdb_id}_{chain_id}_d{i+1}")

            # Get basic fields
            start = domain_dict.get('start', 0)
            end = domain_dict.get('end', 0)
            range_text = domain_dict.get('range', f"{start}-{end}")

            # Fix the evidence access - always use .get() for dictionaries
            evidence_list = []
            # NOTE: Fixed this line to use .get() instead of direct attribute access
            for evidence in domain_dict.get('evidence', []):
                if isinstance(evidence, DomainEvidence):
                    evidence_list.append(evidence)
                elif isinstance(evidence, dict):
                    try:
                        evidence_list.append(DomainEvidence.from_dict(evidence))
                    except Exception as e:
                        logger.warning(f"Failed to convert evidence dict to DomainEvidence: {str(e)}")
                else:
                    logger.warning(f"Unknown evidence type: {type(evidence)}, trying to convert anyway")
                    try:
                        # Last attempt - try to make a basic evidence object
                        if hasattr(evidence, 'type'):
                            # It's some kind of object with a type attribute
                            evidence_dict = {
                                'type': getattr(evidence, 'type', 'unknown'),
                                'source_id': getattr(evidence, 'source_id', ''),
                                'domain_id': getattr(evidence, 'domain_id', '')
                            }
                            evidence_list.append(DomainEvidence.from_dict(evidence_dict))
                    except Exception as e:
                        logger.warning(f"Could not convert evidence: {str(e)}")

            # Create DomainModel
            domain_model = DomainModel(
                id=domain_id,
                start=start,
                end=end,
                range=range_text,
                t_group=domain_dict.get('t_group'),
                h_group=domain_dict.get('h_group'),
                x_group=domain_dict.get('x_group'),
                a_group=domain_dict.get('a_group'),
                source=domain_dict.get('source', ''),
                confidence=domain_dict.get('confidence', 0.0),
                source_id=domain_dict.get('source_id', ''),
                is_manual_rep=domain_dict.get('is_manual_rep', False),
                is_f70=domain_dict.get('is_f70', False),
                is_f40=domain_dict.get('is_f40', False),
                is_f99=domain_dict.get('is_f99', False),
                evidence=evidence_list
            )
            domain_models.append(domain_model)
        except Exception as e:
            logger.error(f"Error creating DomainModel: {str(e)}")

    return domain_models
