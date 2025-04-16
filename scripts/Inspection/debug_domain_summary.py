#!/usr/bin/env python3
"""
debug_domain_summary.py - Debug tool for analyzing domain summary processing issues

This script focuses on debugging the domain summary XML processing and domain boundary
determination in the pyECOD pipeline.
"""

import os
import sys
import argparse
import logging
import xml.etree.ElementTree as ET
import json
from typing import Dict, Any, List, Optional, Tuple
import re

# Set up logging
def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Configure logging with detailed formatting"""
    log_level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler()]
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    
    return logging.getLogger("ecod.debug_domain_summary")

def inspect_xml_file(xml_file_path: str) -> Dict[str, Any]:
    """
    Inspect a domain summary XML file and return detailed information about its structure
    
    Args:
        xml_file_path: Path to the domain summary XML file
        
    Returns:
        Dict with parsed data and structural information
    """
    logger = logging.getLogger("ecod.debug_domain_summary")
    
    # Make sure file exists
    if not os.path.exists(xml_file_path):
        logger.error(f"XML file does not exist: {xml_file_path}")
        return {"error": f"File does not exist: {xml_file_path}"}
    
    try:
        # Parse the XML
        tree = ET.parse(xml_file_path)
        root = tree.getroot()
        
        # Basic file info
        file_info = {
            "file_path": xml_file_path,
            "file_size": os.path.getsize(xml_file_path),
            "root_tag": root.tag,
            "root_attributes": dict(root.attrib),
            "child_elements": []
        }
        
        # Analyze structure - first level
        for child in root:
            child_info = {
                "tag": child.tag,
                "attributes": dict(child.attrib),
                "has_text": bool(child.text and child.text.strip()),
                "num_children": len(child)
            }
            file_info["child_elements"].append(child_info)
        
        # Get blast_summ element and its attributes
        blast_summ = root.find(".//blast_summ")
        if blast_summ is not None:
            file_info["blast_summ_attributes"] = dict(blast_summ.attrib)
            if blast_summ.get("pdb") and blast_summ.get("chain"):
                file_info["pdb_chain"] = f"{blast_summ.get('pdb')}_{blast_summ.get('chain')}"
        
        # Extract query length
        query_len_elem = root.find(".//query_len")
        if query_len_elem is not None and query_len_elem.text:
            try:
                file_info["sequence_length"] = int(query_len_elem.text.strip())
                logger.debug(f"Found sequence length: {file_info['sequence_length']}")
            except ValueError:
                logger.warning(f"Invalid sequence length: {query_len_elem.text}")
        
        # Look for sequence content
        query_seq_elems = root.findall(".//query_seq")
        if query_seq_elems:
            # Use the first one we find
            seq_text = query_seq_elems[0].text.strip() if query_seq_elems[0].text else ""
            file_info["sequence_sample"] = (seq_text[:50] + "...") if len(seq_text) > 50 else seq_text
            if "sequence_length" not in file_info and seq_text:
                file_info["sequence_length"] = len(seq_text)
        
        # Extract hit information
        hit_counts = {}
        
        # Chain BLAST hits
        chain_blast_hits = root.findall(".//chain_blast_run/hits/hit")
        hit_counts["chain_blast"] = len(chain_blast_hits)
        
        # Track hit statistics for chain BLAST
        if chain_blast_hits:
            blast_stats = {
                "e_values": [],
                "ranges": [],
                "pdb_ids": [],
                "chain_ids": []
            }
            
            for hit in chain_blast_hits[:min(10, len(chain_blast_hits))]:  # Sample first 10
                try:
                    # Extract attributes
                    evalue_str = hit.get("evalues", "1")
                    e_value = float(evalue_str) if evalue_str != "0.0" else 1e-200
                    pdb_id = hit.get("pdb_id", "")
                    chain_id = hit.get("chain_id", "")
                    
                    # Extract query region
                    query_reg_elem = hit.find("query_reg")
                    hit_range = query_reg_elem.text.strip() if query_reg_elem is not None and query_reg_elem.text else ""
                    
                    blast_stats["e_values"].append(e_value)
                    blast_stats["ranges"].append(hit_range)
                    blast_stats["pdb_ids"].append(pdb_id)
                    blast_stats["chain_ids"].append(chain_id)
                except (ValueError, TypeError, AttributeError) as e:
                    logger.warning(f"Could not parse hit values: {str(e)}")
            
            file_info["chain_blast_stats"] = blast_stats
        
        # Domain BLAST hits
        domain_blast_hits = root.findall(".//blast_run/hits/hit")
        hit_counts["domain_blast"] = len(domain_blast_hits)
        
        # Track hit statistics for domain BLAST
        if domain_blast_hits:
            domain_stats = {
                "e_values": [],
                "ranges": [],
                "domain_ids": []
            }
            
            for hit in domain_blast_hits[:min(10, len(domain_blast_hits))]:  # Sample first 10
                try:
                    # Extract attributes
                    evalue_str = hit.get("evalues", "1")
                    e_value = float(evalue_str) if evalue_str != "0.0" else 1e-200
                    domain_id = hit.get("domain_id", "")
                    
                    # Extract query region
                    query_reg_elem = hit.find("query_reg")
                    hit_range = query_reg_elem.text.strip() if query_reg_elem is not None and query_reg_elem.text else ""
                    
                    domain_stats["e_values"].append(e_value)
                    domain_stats["ranges"].append(hit_range)
                    domain_stats["domain_ids"].append(domain_id)
                except (ValueError, TypeError, AttributeError) as e:
                    logger.warning(f"Could not parse hit values: {str(e)}")
            
            file_info["domain_blast_stats"] = domain_stats
        
        file_info["hit_counts"] = hit_counts
        
        # Check for additional hit types (HHSearch, etc.)
        for hit_type in ["hhsearch_hit", "self_hit"]:
            hits = root.findall(f".//{hit_type}")
            if hits:
                hit_counts[hit_type] = len(hits)
        
        # Sample raw XML for first 1000 chars
        with open(xml_file_path, 'r') as f:
            xml_content = f.read(1000)
            file_info["xml_sample"] = xml_content + ("..." if len(xml_content) >= 1000 else "")
        
        return file_info
        
    except ET.ParseError as e:
        logger.error(f"XML parsing error: {str(e)}")
        
        # Try to extract partial XML content
        try:
            with open(xml_file_path, 'r') as f:
                xml_content = f.read(2000)  # Read more to see the error context
                
            return {
                "error": f"XML parsing error: {str(e)}",
                "xml_sample": xml_content + "..."
            }
        except Exception as inner_e:
            return {"error": f"XML parsing error: {str(e)}, couldn't read file: {str(inner_e)}"}
            
    except Exception as e:
        logger.error(f"Error inspecting XML: {str(e)}", exc_info=True)
        return {"error": f"Error inspecting XML: {str(e)}"}

def analyze_hit_quality(xml_file_path: str, thresholds: Dict[str, float] = None) -> Dict[str, Any]:
    """
    Analyze hit quality in a domain summary file with different thresholds
    
    Args:
        xml_file_path: Path to the domain summary XML file
        thresholds: Dict of threshold values for different parameters
        
    Returns:
        Dict with hit quality analysis
    """
    logger = logging.getLogger("ecod.debug_domain_summary")
    
    # Default thresholds
    if thresholds is None:
        thresholds = {
            "blast_evalue": 1e-3,
            "blast_identity": 30.0,
            "blast_coverage": 50.0,
            "hhsearch_probability": 90.0,
            "hhsearch_evalue": 1e-3
        }
    
    try:
        # Parse the XML
        tree = ET.parse(xml_file_path)
        root = tree.getroot()
        
        # Get sequence length
        query_len_elem = root.find(".//query_len")
        sequence_length = int(query_len_elem.text.strip()) if query_len_elem is not None and query_len_elem.text else 0
        
        # If no query_len, try to get from a query_seq
        if sequence_length == 0:
            query_seq_elems = root.findall(".//query_seq")
            if query_seq_elems and query_seq_elems[0].text:
                sequence_length = len(query_seq_elems[0].text.strip())
                
        # Analysis results
        results = {
            "file_path": xml_file_path,
            "sequence_length": sequence_length,
            "thresholds_used": thresholds,
            "chain_blast_hits": {
                "total": 0,
                "significant": 0,
                "coverage": 0.0,
                "hits_by_significance": []
            },
            "domain_blast_hits": {
                "total": 0,
                "significant": 0,
                "coverage": 0.0,
                "hits_by_significance": []
            }
        }
        
        # Get PDB chain info
        blast_summ = root.find(".//blast_summ")
        if blast_summ is not None:
            pdb_id = blast_summ.get("pdb", "")
            chain_id = blast_summ.get("chain", "")
            if pdb_id and chain_id:
                results["pdb_chain"] = f"{pdb_id}_{chain_id}"
        
        # Process chain BLAST hits
        chain_blast_hits = []
        sequence_coverage = [0] * (sequence_length + 1) if sequence_length > 0 else []
        
        for hit_elem in root.findall(".//chain_blast_run/hits/hit"):
            try:
                hit_info = dict(hit_elem.attrib)
                
                # Extract query region
                query_reg_elem = hit_elem.find("query_reg")
                if query_reg_elem is not None and query_reg_elem.text:
                    hit_info["range"] = query_reg_elem.text.strip()
                    hit_info["range_parsed"] = parse_range(hit_info["range"])
                
                # Convert evalue
                if "evalues" in hit_info:
                    evalue_str = hit_info["evalues"]
                    hit_info["evalue"] = float(evalue_str) if evalue_str != "0.0" else 1e-200
                else:
                    hit_info["evalue"] = 1.0  # Default to a poor E-value
                
                # Add identity and coverage (not directly available, estimate from sequences)
                query_seq_elem = hit_elem.find("query_seq")
                hit_seq_elem = hit_elem.find("hit_seq")
                
                if query_seq_elem is not None and hit_seq_elem is not None:
                    query_seq = query_seq_elem.text.strip() if query_seq_elem.text else ""
                    hit_seq = hit_seq_elem.text.strip() if hit_seq_elem.text else ""
                    
                    if query_seq and hit_seq and len(query_seq) == len(hit_seq):
                        # Calculate identity percentage
                        matches = sum(1 for q, h in zip(query_seq, hit_seq) if q == h)
                        identity = (matches / len(query_seq)) * 100
                        hit_info["identity"] = identity
                    else:
                        # Default identity
                        hit_info["identity"] = 90.0  # High identity for chain hits
                else:
                    hit_info["identity"] = 90.0
                
                # Estimate coverage from region
                if "range_parsed" in hit_info and hit_info["range_parsed"] and sequence_length > 0:
                    total_covered = sum(end - start + 1 for start, end in hit_info["range_parsed"])
                    hit_info["query_coverage"] = (total_covered / sequence_length) * 100
                else:
                    hit_info["query_coverage"] = 90.0  # Default high coverage
                
                # Check significance
                is_significant = (
                    hit_info["evalue"] <= thresholds["blast_evalue"] and
                    hit_info["identity"] >= thresholds["blast_identity"] and
                    hit_info["query_coverage"] >= thresholds["blast_coverage"]
                )
                
                hit_info["is_significant"] = is_significant
                
                # Update coverage if significant
                if is_significant and "range_parsed" in hit_info and sequence_length > 0:
                    for range_start, range_end in hit_info["range_parsed"]:
                        for i in range(max(1, range_start), min(sequence_length + 1, range_end + 1)):
                            sequence_coverage[i-1] = 1
                
                chain_blast_hits.append(hit_info)
                
            except Exception as e:
                logger.warning(f"Error processing chain BLAST hit: {str(e)}")
        
        # Sort by E-value
        chain_blast_hits.sort(key=lambda x: x.get("evalue", 1.0))
        
        # Update results
        results["chain_blast_hits"]["total"] = len(chain_blast_hits)
        results["chain_blast_hits"]["significant"] = sum(1 for hit in chain_blast_hits if hit.get("is_significant", False))
        
        # Calculate coverage
        if sequence_length > 0:
            covered_positions = sum(sequence_coverage)
            results["chain_blast_hits"]["coverage"] = (covered_positions / sequence_length) * 100.0
        
        # Store top hits for review
        results["chain_blast_hits"]["hits_by_significance"] = chain_blast_hits[:10]  # Top 10
        
        # Process domain BLAST hits
        domain_blast_hits = []
        domain_coverage = [0] * (sequence_length + 1) if sequence_length > 0 else []
        
        for hit_elem in root.findall(".//blast_run/hits/hit"):
            try:
                hit_info = dict(hit_elem.attrib)
                
                # Extract query region
                query_reg_elem = hit_elem.find("query_reg")
                if query_reg_elem is not None and query_reg_elem.text:
                    hit_info["range"] = query_reg_elem.text.strip()
                    hit_info["range_parsed"] = parse_range(hit_info["range"])
                
                # Convert evalue
                if "evalues" in hit_info:
                    evalue_str = hit_info["evalues"]
                    hit_info["evalue"] = float(evalue_str) if evalue_str != "0.0" else 1e-200
                else:
                    hit_info["evalue"] = 1.0  # Default to a poor E-value
                
                # Add domain ID if present
                if "domain_id" in hit_info:
                    hit_info["target_id"] = hit_info["domain_id"]
                
                # Add identity and coverage (estimate from region)
                hit_info["identity"] = 90.0  # Use a reasonable default for domain hits
                
                # Estimate coverage from region
                if "range_parsed" in hit_info and hit_info["range_parsed"] and sequence_length > 0:
                    total_covered = sum(end - start + 1 for start, end in hit_info["range_parsed"])
                    hit_info["query_coverage"] = (total_covered / sequence_length) * 100
                else:
                    hit_info["query_coverage"] = 0.0
                
                # Check significance
                is_significant = (
                    hit_info["evalue"] <= thresholds["blast_evalue"] and
                    hit_info["identity"] >= thresholds["blast_identity"] and
                    hit_info["query_coverage"] >= thresholds["blast_coverage"]
                )
                
                hit_info["is_significant"] = is_significant
                
                # Update coverage if significant
                if is_significant and "range_parsed" in hit_info and sequence_length > 0:
                    for range_start, range_end in hit_info["range_parsed"]:
                        for i in range(max(1, range_start), min(sequence_length + 1, range_end + 1)):
                            domain_coverage[i-1] = 1
                
                domain_blast_hits.append(hit_info)
                
            except Exception as e:
                logger.warning(f"Error processing domain BLAST hit: {str(e)}")
        
        # Sort by E-value
        domain_blast_hits.sort(key=lambda x: x.get("evalue", 1.0))
        
        # Update results
        results["domain_blast_hits"]["total"] = len(domain_blast_hits)
        results["domain_blast_hits"]["significant"] = sum(1 for hit in domain_blast_hits if hit.get("is_significant", False))
        
        # Calculate coverage
        if sequence_length > 0:
            covered_positions = sum(domain_coverage)
            results["domain_blast_hits"]["coverage"] = (covered_positions / sequence_length) * 100.0
        
        # Store top hits for review
        results["domain_blast_hits"]["hits_by_significance"] = domain_blast_hits[:10]  # Top 10
        
        # Check for HHSearch hits (if any)
        hhsearch_hits = root.findall(".//hhsearch_hit")
        if hhsearch_hits:
            results["hhsearch_hits"] = {
                "total": len(hhsearch_hits),
                "significant": 0,  # Would need to process these
                "coverage": 0.0,
                "hits_by_significance": []
            }
        
        # Overall assessment
        results["overall_assessment"] = {
            "has_significant_chain_blast_hits": results["chain_blast_hits"]["significant"] > 0,
            "has_significant_domain_blast_hits": results["domain_blast_hits"]["significant"] > 0,
            "has_any_significant_hits": (
                results["chain_blast_hits"]["significant"] > 0 or 
                results["domain_blast_hits"]["significant"] > 0
            )
        }
        
        # Suggest domains based on hits
        results["domain_suggestion"] = suggest_domains_revised(
            chain_blast_hits, domain_blast_hits, sequence_length,
            thresholds
        )
        
        return results
        
    except Exception as e:
        logger.error(f"Error analyzing hit quality: {str(e)}", exc_info=True)
        return {"error": f"Error analyzing hit quality: {str(e)}"}

def parse_range(range_str: str) -> List[Tuple[int, int]]:
    """
    Parse a range string like "1-100,150-200" into a list of tuples [(1,100), (150,200)]
    
    Args:
        range_str: Range string to parse
        
    Returns:
        List of (start, end) tuples
    """
    if not range_str or range_str.strip() == "":
        return []
    
    ranges = []
    parts = range_str.split(',')
    
    for part in parts:
        if '-' in part:
            try:
                start, end = part.split('-')
                start_num = int(re.sub(r'[^0-9]', '', start))
                end_num = int(re.sub(r'[^0-9]', '', end))
                ranges.append((start_num, end_num))
            except (ValueError, IndexError):
                pass
    
    return ranges

def suggest_domains_revised(chain_blast_hits: List[Dict[str, Any]], 
                        domain_blast_hits: List[Dict[str, Any]], 
                        sequence_length: int,
                        thresholds: Dict[str, float]) -> Dict[str, Any]:
    """
    Suggest domain boundaries based on hit information from the actual XML format
    
    This is updated to work with the actual domain summary XML structure
    
    Args:
        chain_blast_hits: List of chain BLAST hits with parsed attributes
        domain_blast_hits: List of domain BLAST hits with parsed attributes
        sequence_length: Length of the protein sequence
        thresholds: Threshold values for significance
        
    Returns:
        Dict with domain suggestions and reasoning
    """
    if sequence_length <= 0:
        return {"error": "Invalid sequence length"}
    
    # Initialize coverage arrays
    position_coverage = [0] * (sequence_length + 1)  # 1-indexed
    
    # Track hits contributing to each position
    position_evidence = [[] for _ in range(sequence_length + 1)]
    
    # Process significant chain BLAST hits
    for hit in chain_blast_hits:
        if not hit.get("is_significant", False):
            continue
            
        if "range_parsed" not in hit or not hit["range_parsed"]:
            continue
            
        hit_id = f"{hit.get('pdb_id', '')}_{hit.get('chain_id', '')} (Chain)"
        
        for start, end in hit["range_parsed"]:
            for i in range(max(1, start), min(sequence_length + 1, end + 1)):
                position_coverage[i-1] += 1
                position_evidence[i-1].append(hit_id)
    
    # Process significant domain BLAST hits
    for hit in domain_blast_hits:
        if not hit.get("is_significant", False):
            continue
            
        if "range_parsed" not in hit or not hit["range_parsed"]:
            continue
            
        hit_id = hit.get("domain_id", "unknown") + " (Domain)"
        
        for start, end in hit["range_parsed"]:
            for i in range(max(1, start), min(sequence_length + 1, end + 1)):
                position_coverage[i-1] += 1
                position_evidence[i-1].append(hit_id)
    
    # Find contiguous regions with coverage
    regions = []
    current_region = None
    
    for i in range(sequence_length):
        if position_coverage[i] > 0:
            if current_region is None:
                current_region = {
                    "start": i + 1,
                    "evidence": set(position_evidence[i])
                }
        else:
            if current_region is not None:
                current_region["end"] = i
                regions.append(current_region)
                current_region = None
    
    # Add final region if exists
    if current_region is not None:
        current_region["end"] = sequence_length
        regions.append(current_region)
    
    # Merge overlapping or close regions (within 30 residues)
    merged_regions = []
    if regions:
        merged_regions.append(regions[0])
        
        for region in regions[1:]:
            prev_region = merged_regions[-1]
            
            # Check if regions are close enough to merge
            if region["start"] - prev_region["end"] <= 30:
                # Merge regions
                prev_region["end"] = region["end"]
                prev_region["evidence"] = prev_region["evidence"].union(region["evidence"])
            else:
                # Add as separate region
                merged_regions.append(region)
    
    # Calculate region sizes and coverage
    for region in merged_regions:
        region["size"] = region["end"] - region["start"] + 1
        region["evidence"] = list(region["evidence"])
        region["coverage"] = region["size"] / sequence_length * 100.0
    
    # Make domain recommendations
    domains = []
    
    # If no regions or very small regions, recommend whole chain
    if not merged_regions or all(r["size"] < 30 for r in merged_regions):
        domains.append({
            "start": 1,
            "end": sequence_length,
            "confidence": "low",
            "reason": "No significant hits found, using whole chain"
        })
    else:
        # Add regions as domains with confidence
        for region in merged_regions:
            confidence = "high" if region["size"] >= 50 else "medium"
            
            domains.append({
                "start": region["start"],
                "end": region["end"],
                "size": region["size"],
                "coverage": region["coverage"],
                "evidence_count": len(region["evidence"]),
                "evidence": region["evidence"][:5],  # Limit to first 5 pieces of evidence
                "confidence": confidence,
                "reason": f"Domain determined by {len(region['evidence'])} significant hits"
            })
    
    return {
        "domains": domains,
        "num_domains": len(domains),
        "sequence_length": sequence_length,
        "covered_positions": sum(1 for c in position_coverage if c > 0),
        "coverage_percent": sum(1 for c in position_coverage if c > 0) / sequence_length * 100.0 if sequence_length > 0 else 0.0
    }

def suggest_domains(blast_hits: List[Dict[str, Any]], 
                   hhsearch_hits: List[Dict[str, Any]], 
                   sequence_length: int,
                   thresholds: Dict[str, float]) -> Dict[str, Any]:
    """
    Legacy function to maintain compatibility
    """
    # Convert to new format
    logger = logging.getLogger("ecod.debug_domain_summary")
    logger.warning("Using legacy suggest_domains function - consider updating to suggest_domains_revised")
    
    # Forward to new function
    return suggest_domains_revised(blast_hits, hhsearch_hits, sequence_length, thresholds)

def output_report(results: Dict[str, Any], output_file: Optional[str] = None) -> None:
    """
    Output a report of the domain analysis
    
    Args:
        results: Analysis results
        output_file: Optional output file for JSON results
    """
    logger = logging.getLogger("ecod.debug_domain_summary")
    
    # Print summary to console
    print(f"\n=== Domain Summary Analysis Report ===\n")
    
    # Basic file info
    if "error" in results:
        print(f"ERROR: {results['error']}")
        if "xml_sample" in results:
            print("\nPartial XML content:")
            print(results["xml_sample"][:500])
        return
    
    # File information
    print(f"File: {results.get('file_path', 'Unknown')}")
    print(f"Sequence length: {results.get('sequence_length', 'Unknown')}")
    
    # Hit information
    if "hit_counts" in results:
        print("\nHit counts:")
        for hit_type, count in results["hit_counts"].items():
            print(f"  - {hit_type}: {count}")
    
    # BLAST statistics summary
    if "blast_stats" in results:
        stats = results["blast_stats"]
        if stats.get("e_values"):
            min_evalue = min(stats["e_values"])
            max_identity = max(stats["identities"]) if stats.get("identities") else 0
            max_coverage = max(stats["query_coverage"]) if stats.get("query_coverage") else 0
            
            print("\nBLAST statistics:")
            print(f"  - Best E-value: {min_evalue:.2e}")
            print(f"  - Best identity: {max_identity:.1f}%")
            print(f"  - Best query coverage: {max_coverage:.1f}%")
    
    # Hit quality analysis
    if "blast_hits" in results and "blast_hits" not in results.get("hit_counts", {}):
        print("\nBLAST hit analysis:")
        print(f"  - Total hits: {results['blast_hits']['total']}")
        print(f"  - Significant hits: {results['blast_hits']['significant']}")
        print(f"  - Query coverage: {results['blast_hits']['coverage']:.1f}%")
        
        print("\nHHSearch hit analysis:")
        print(f"  - Total hits: {results['hhsearch_hits']['total']}")
        print(f"  - Significant hits: {results['hhsearch_hits']['significant']}")
        print(f"  - Query coverage: {results['hhsearch_hits']['coverage']:.1f}%")
    
    # Domain suggestions
    if "domain_suggestion" in results:
        ds = results["domain_suggestion"]
        print("\nDomain suggestions:")
        print(f"  - Number of domains: {ds.get('num_domains', 0)}")
        print(f"  - Sequence coverage: {ds.get('coverage_percent', 0):.1f}%")
        
        if "domains" in ds:
            print("\nProposed domains:")
            for i, domain in enumerate(ds["domains"]):
                print(f"  Domain {i+1}:")
                print(f"    - Range: {domain.get('start', '?')}-{domain.get('end', '?')}")
                print(f"    - Size: {domain.get('size', '?')} residues")
                print(f"    - Confidence: {domain.get('confidence', 'unknown')}")
                print(f"    - Reason: {domain.get('reason', 'unknown')}")
                
                if "evidence" in domain:
                    evidence_str = ", ".join(domain["evidence"])
                    if len(domain.get("evidence", [])) < domain.get("evidence_count", 0):
                        evidence_str += f" and {domain.get('evidence_count', 0) - len(domain.get('evidence', []))} more"
                    print(f"    - Evidence: {evidence_str}")
                    
                print("")
    
    # Save full results to file if requested
    if output_file:
        try:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            print(f"\nDetailed results saved to: {output_file}")
        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")
            print(f"\nError saving results: {str(e)}")

def test_new_thresholds(xml_file_path: str) -> Dict[str, Any]:
    """
    Test multiple threshold combinations to find optimal values
    
    Args:
        xml_file_path: Path to the domain summary XML file
        
    Returns:
        Dict with results from different threshold combinations
    """
    logger = logging.getLogger("ecod.debug_domain_summary")
    
    # Define threshold combinations to test
    threshold_sets = [
        {
            "name": "Default",
            "thresholds": {
                "blast_evalue": 1e-3,
                "blast_identity": 30.0,
                "blast_coverage": 50.0,
                "hhsearch_probability": 90.0,
                "hhsearch_evalue": 1e-3
            }
        },
        {
            "name": "Relaxed BLAST",
            "thresholds": {
                "blast_evalue": 1e-1,
                "blast_identity": 20.0,
                "blast_coverage": 30.0,
                "hhsearch_probability": 90.0,
                "hhsearch_evalue": 1e-3
            }
        },
        {
            "name": "Stringent BLAST",
            "thresholds": {
                "blast_evalue": 1e-5,
                "blast_identity": 40.0,
                "blast_coverage": 70.0,
                "hhsearch_probability": 90.0,
                "hhsearch_evalue": 1e-3
            }
        },
        {
            "name": "Relaxed HHSearch",
            "thresholds": {
                "blast_evalue": 1e-3,
                "blast_identity": 30.0,
                "blast_coverage": 50.0,
                "hhsearch_probability": 70.0,
                "hhsearch_evalue": 1e-1
            }
        }
    ]
    
    results = {}
    
    for threshold_set in threshold_sets:
        logger.info(f"Testing threshold set: {threshold_set['name']}")
        try:
            analysis = analyze_hit_quality(xml_file_path, threshold_set["thresholds"])
            
            # Extract relevant summary data
            summary = {
                "thresholds": threshold_set["thresholds"],
                "blast_significant_hits": analysis.get("blast_hits", {}).get("significant", 0),
                "hhsearch_significant_hits": analysis.get("hhsearch_hits", {}).get("significant", 0),
                "num_domains": analysis.get("domain_suggestion", {}).get("num_domains", 0),
                "domains": analysis.get("domain_suggestion", {}).get("domains", [])
            }
            
            # Simplify domain data for comparison
            domain_summary = []
            for domain in summary.get("domains", []):
                domain_summary.append({
                    "range": f"{domain.get('start', 0)}-{domain.get('end', 0)}",
                    "size": domain.get("size", 0),
                    "confidence": domain.get("confidence", "unknown")
                })
            
            summary["domain_summary"] = domain_summary
            results[threshold_set["name"]] = summary
            
        except Exception as e:
            logger.error(f"Error with threshold set {threshold_set['name']}: {str(e)}")
            results[threshold_set["name"]] = {"error": str(e)}
    
    return results

def main():
    """Main entry point for the script"""
    parser = argparse.ArgumentParser(description='Debug domain summary processing')
    parser.add_argument('--file', type=str, required=True,
                      help='Domain summary XML file to analyze')
    parser.add_argument('--output', type=str, default=None,
                      help='Output JSON file for detailed results')
    parser.add_argument('--test-thresholds', action='store_true',
                      help='Test multiple threshold combinations')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    logger.info(f"Analyzing domain summary file: {args.file}")
    
    # Basic structure check
    file_info = inspect_xml_file(args.file)
    
    if "error" in file_info:
        logger.error(f"Error inspecting file: {file_info['error']}")
        output_report(file_info, args.output)
        return 1
    
    # Analyze hit quality
    hit_analysis = analyze_hit_quality(args.file)
    if "error" in hit_analysis:
        logger.error(f"Error analyzing hits: {hit_analysis['error']}")
        file_info.update({"hit_analysis_error": hit_analysis["error"]})
        output_report(file_info, args.output)
        return 1
    
    # Combine results
    combined_results = {**file_info, **hit_analysis}
    
    # Test thresholds if requested
    if args.test_thresholds:
        logger.info("Testing multiple threshold combinations")
        threshold_results = test_new_thresholds(args.file)
        combined_results["threshold_tests"] = threshold_results
    
    # Output report
    output_report(combined_results, args.output)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())