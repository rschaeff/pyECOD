#!/usr/bin/env python3
"""
analyze_domain_summaries.py - Comprehensive analysis of domain summary files
"""

import os
import sys
import json
import logging
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Any, Optional
from collections import Counter, defaultdict

# Add parent directories to path for imports
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, project_root)

from ecod.core.context import ApplicationContext

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Configure logging"""
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

def analyze_xml_file(file_path: str, detailed: bool = False) -> Dict[str, Any]:
    """Analyze a domain summary XML file with comprehensive evidence checking"""
    logger = logging.getLogger("ecod.analyze.file")
    
    result = {
        "file_path": file_path,
        "valid_xml": False,
        "error": None,
        "evidence_types": [],
        "evidence_counts": {}
    }
    
    if not os.path.exists(file_path):
        result["error"] = "File not found"
        return result
    
    try:
        # Parse XML
        tree = ET.parse(file_path)
        root = tree.getroot()
        
        # Mark as valid XML
        result["valid_xml"] = True
        
        # Get basic information
        if root.tag == 'domain_summ_doc':
            # Extract metadata
            metadata = root.find("metadata")
            if metadata is not None:
                result["pdb_id"] = metadata.findtext("pdb_id", "unknown")
                result["chain_id"] = metadata.findtext("chain_id", "unknown")
                result["reference"] = metadata.findtext("reference", "unknown")
                result["creation_date"] = metadata.findtext("creation_date", "unknown")
            else:
                # Try alternative locations for PDB/chain info
                blast_summ = root.find('.//blast_summ')
                if blast_summ is not None:
                    result["pdb_id"] = blast_summ.get('pdb', "unknown")
                    result["chain_id"] = blast_summ.get('chain', "unknown")
            
            # Check evidence types
            # 1. Chain BLAST evidence
            chain_blast = root.find("chain_blast_evidence")
            if chain_blast is not None:
                chain_hits = chain_blast.findall(".//hit")
                result["chain_blast_hits"] = len(chain_hits)
                if len(chain_hits) > 0:
                    result["evidence_types"].append("chain_blast")
                    result["evidence_counts"]["chain_blast"] = len(chain_hits)
                    
                    # Top hits details if requested
                    if detailed and chain_hits:
                        chain_hit_details = []
                        for hit in chain_hits[:5]:  # Top 5 hits
                            hit_detail = {
                                "pdb_id": hit.get("pdb_id", "unknown"),
                                "chain_id": hit.get("chain_id", "unknown"),
                                "evalues": hit.get("evalues", "unknown")
                            }
                            
                            # Get regions
                            query_reg = hit.find("query_reg")
                            if query_reg is not None and query_reg.text:
                                hit_detail["query_region"] = query_reg.text
                                
                            hit_reg = hit.find("hit_reg")
                            if hit_reg is not None and hit_reg.text:
                                hit_detail["hit_region"] = hit_reg.text
                                
                            chain_hit_details.append(hit_detail)
                            
                        result["chain_blast_top_hits"] = chain_hit_details
            
            # 2. Domain BLAST evidence
            domain_blast = root.find("domain_blast_evidence")
            if domain_blast is not None:
                domain_hits = domain_blast.findall(".//hit")
                result["domain_blast_hits"] = len(domain_hits)
                if len(domain_hits) > 0:
                    result["evidence_types"].append("domain_blast")
                    result["evidence_counts"]["domain_blast"] = len(domain_hits)
                    
                    # Hit details if requested
                    if detailed and domain_hits:
                        domain_hit_details = []
                        for hit in domain_hits[:5]:  # Top 5 hits
                            hit_detail = {
                                "domain_id": hit.get("domain_id", "unknown"),
                                "pdb_id": hit.get("pdb_id", "unknown"),
                                "chain_id": hit.get("chain_id", "unknown"),
                                "evalues": hit.get("evalues", "unknown")
                            }
                            
                            # Get regions
                            query_reg = hit.find("query_reg")
                            if query_reg is not None and query_reg.text:
                                hit_detail["query_region"] = query_reg.text
                                
                            hit_reg = hit.find("hit_reg")
                            if hit_reg is not None and hit_reg.text:
                                hit_detail["hit_region"] = hit_reg.text
                                
                            domain_hit_details.append(hit_detail)
                            
                        result["domain_blast_top_hits"] = domain_hit_details
            
            # 3. HHSearch evidence
            hhsearch = root.find("hhsearch_evidence")
            if hhsearch is not None:
                hh_hits = hhsearch.findall(".//hh_hit")
                result["hhsearch_hits"] = len(hh_hits)
                if len(hh_hits) > 0:
                    result["evidence_types"].append("hhsearch")
                    result["evidence_counts"]["hhsearch"] = len(hh_hits)
                    
                    # Hit details if requested
                    if detailed and hh_hits:
                        hh_hit_details = []
                        for hit in hh_hits[:5]:  # Top 5 hits
                            hit_detail = {
                                "hit_num": hit.get("hit_num", "unknown"),
                                "domain_id": hit.get("domain_id", "unknown"),
                                "probability": hit.get("hh_prob", "unknown"),
                                "score": hit.get("hh_score", "unknown")
                            }
                            
                            # Get regions
                            query_reg = hit.find("query_reg")
                            if query_reg is not None and query_reg.text:
                                hit_detail["query_region"] = query_reg.text
                                
                            hit_reg = hit.find("hit_reg")
                            if hit_reg is not None and hit_reg.text:
                                hit_detail["hit_region"] = hit_reg.text
                                
                            hh_hit_details.append(hit_detail)
                            
                        result["hhsearch_top_hits"] = hh_hit_details
            
            # 4. Domain definitions
            domains = root.findall(".//domain")
            result["domain_count"] = len(domains)
            if len(domains) > 0:
                result["has_domains"] = True
                
                # Domain details if requested
                if detailed and domains:
                    domain_details = []
                    for domain in domains:
                        domain_detail = {
                            "id": domain.get("id", "unknown"),
                            "range": domain.get("range", "unknown")
                        }
                        
                        # Get classification if available
                        classification = domain.find("classification")
                        if classification is not None:
                            domain_detail["t_group"] = classification.get("t_group", "")
                            domain_detail["h_group"] = classification.get("h_group", "")
                            domain_detail["x_group"] = classification.get("x_group", "")
                            domain_detail["a_group"] = classification.get("a_group", "")
                            
                        domain_details.append(domain_detail)
                        
                    result["domain_details"] = domain_details
                    
                    # Count classifications
                    t_groups = Counter()
                    h_groups = Counter()
                    
                    for domain in domain_details:
                        if "t_group" in domain and domain["t_group"]:
                            t_groups[domain["t_group"]] += 1
                        if "h_group" in domain and domain["h_group"]:
                            h_groups[domain["h_group"]] += 1
                    
                    result["t_group_counts"] = dict(t_groups)
                    result["h_group_counts"] = dict(h_groups)
            else:
                result["has_domains"] = False
                
            # Determine if no evidence whatsoever
            if not result["evidence_types"]:
                result["no_evidence"] = True
                logger.warning(f"NO EVIDENCE FOUND in {file_path}")
            else:
                result["no_evidence"] = False
                
                # Create combined evidence key
                result["evidence_key"] = "+".join(sorted(result["evidence_types"]))
                
        # Alternative format (blast_summ_doc)
        elif root.tag == 'blast_summ_doc':
            blast_summ = root.find('blast_summ')
            if blast_summ is not None:
                result["pdb_id"] = blast_summ.get('pdb')
                result["chain_id"] = blast_summ.get('chain')
                
                # Check for peptide flag
                if blast_summ.get('is_peptide') == 'true':
                    result["is_peptide"] = True
                    result["sequence_length"] = int(blast_summ.get('sequence_length', '0'))
                
                # Get chain BLAST run info
                chain_blast = blast_summ.find('chain_blast_run')
                if chain_blast is not None:
                    result["blast_program"] = chain_blast.get('program')
                    result["blast_version"] = chain_blast.get('version')
                    
                    # Get blast DB
                    blast_db = chain_blast.find('blast_db')
                    if blast_db is not None and blast_db.text:
                        result["blast_db"] = blast_db.text
                    
                    # Get query length
                    query_len = chain_blast.find('query_len')
                    if query_len is not None and query_len.text:
                        result["query_length"] = int(query_len.text)
                    
                    # Get hits
                    hits_elem = chain_blast.find('hits')
                    if hits_elem is not None:
                        chain_hits = hits_elem.findall('hit')
                        result["chain_blast_hits"] = len(chain_hits)
                        if len(chain_hits) > 0:
                            result["evidence_types"].append("chain_blast")
                            result["evidence_counts"]["chain_blast"] = len(chain_hits)
                
                # Get domain BLAST run info
                domain_blast = blast_summ.find('blast_run')
                if domain_blast is not None:
                    # Get hits
                    hits_elem = domain_blast.find('hits')
                    if hits_elem is not None:
                        domain_hits = hits_elem.findall('hit')
                        result["domain_blast_hits"] = len(domain_hits)
                        if len(domain_hits) > 0:
                            result["evidence_types"].append("domain_blast")
                            result["evidence_counts"]["domain_blast"] = len(domain_hits)
                
                # Check for HH search evidence
                hh_run = blast_summ.find('hh_run')
                if hh_run is not None:
                    hits_elem = hh_run.find('hits')
                    if hits_elem is not None:
                        hh_hits = hits_elem.findall('hit')
                        result["hhsearch_hits"] = len(hh_hits)
                        if len(hh_hits) > 0:
                            result["evidence_types"].append("hhsearch")
                            result["evidence_counts"]["hhsearch"] = len(hh_hits)
                
                # Determine if no evidence whatsoever
                if not result["evidence_types"]:
                    result["no_evidence"] = True
                    logger.warning(f"NO EVIDENCE FOUND in {file_path}")
                else:
                    result["no_evidence"] = False
                    
                    # Create combined evidence key
                    result["evidence_key"] = "+".join(sorted(result["evidence_types"]))
                    
        else:
            result["error"] = f"Unknown root tag: {root.tag}"
            
        return result
    
    except ET.ParseError as e:
        result["error"] = f"XML parsing error: {str(e)}"
        return result
    except Exception as e:
        result["error"] = f"Error analyzing file: {str(e)}"
        logger.error(f"Error analyzing {file_path}: {str(e)}", exc_info=True)
        return result

class DomainSummaryAnalyzer:
    """Analyzer for domain summary files"""
    
    def __init__(self, context):
        """Initialize with application context"""
        self.context = context
        self.logger = logging.getLogger("ecod.domain_analyzer")
        
        # Results tracking
        self.total_summaries = 0
        self.evidence_stats = Counter()
        self.no_evidence_proteins = []
        self.evidence_counts = defaultdict(list)
        self.parsing_errors = []
        self.batch_stats = {}
    
    def analyze_batch(self, batch_id: int, sample_size: int = None, detailed: bool = False):
        """Analyze domain summaries for a batch, prioritizing full pipeline summaries"""
        self.logger.info(f"Analyzing domain summaries for batch {batch_id}")
        
        try:
            # Get batch info
            batch_query = """
            SELECT id, batch_name, base_path, ref_version, total_items 
            FROM ecod_schema.batch 
            WHERE id = %s
            """
            
            batch_info = self.context.db.execute_dict_query(batch_query, (batch_id,))
            if not batch_info:
                self.logger.error(f"Batch {batch_id} not found")
                return False
            
            base_path = batch_info[0]['base_path']
            ref_version = batch_info[0]['ref_version']
            batch_name = batch_info[0]['batch_name']
            total_items = batch_info[0]['total_items']
            
            self.logger.info(f"Processing batch: {batch_name} with reference {ref_version}")
            
            # Use a simpler query first to confirm database connectivity
            test_query = "SELECT COUNT(*) FROM ecod_schema.process_status WHERE batch_id = %s"
            count_result = self.context.db.execute_query(test_query, (batch_id,))
            self.logger.info(f"Found {count_result[0][0]} process records for batch {batch_id}")
            
            # Get representative proteins with domain summaries
            if sample_size is None:
                # Include representative filter
                protein_query = """
                SELECT 
                    p.id as protein_id, 
                    p.pdb_id, 
                    p.chain_id, 
                    ps.id as process_id, 
                    pf.file_path,
                    CASE 
                        WHEN pf.file_path LIKE '%%blast_only%%' THEN 'blast_only'
                        ELSE 'full'
                    END as summary_type
                FROM 
                    ecod_schema.process_status ps
                JOIN
                    ecod_schema.protein p ON ps.protein_id = p.id
                JOIN
                    ecod_schema.process_file pf ON ps.id = pf.process_id 
                WHERE 
                    ps.batch_id = %s
                    AND ps.is_representative = TRUE
                    AND pf.file_type = 'domain_summary'
                    AND pf.file_exists = TRUE
                ORDER BY p.pdb_id, p.chain_id
                """
            else:
                # Skip representative filter
                protein_query = """
                SELECT 
                    p.id as protein_id, 
                    p.pdb_id, 
                    p.chain_id, 
                    ps.id as process_id, 
                    pf.file_path,
                    CASE 
                        WHEN pf.file_path LIKE '%%blast_only%%' THEN 'blast_only'
                        ELSE 'full'
                    END as summary_type
                FROM 
                    ecod_schema.process_status ps
                JOIN
                    ecod_schema.protein p ON ps.protein_id = p.id
                JOIN
                    ecod_schema.process_file pf ON ps.id = pf.process_id 
                WHERE 
                    ps.batch_id = %s
                    AND pf.file_type = 'domain_summary'
                    AND pf.file_exists = TRUE
                ORDER BY p.pdb_id, p.chain_id
                """
            
            # Execute the query
            proteins = self.context.db.execute_dict_query(protein_query, (batch_id,))
            
            # Create batch stats entry
            self.batch_stats[batch_id] = {
                "id": batch_id,
                "name": batch_name,
                "base_path": base_path,
                "ref_version": ref_version,
                "total_items": total_items,
                "analyzed_items": len(proteins),
                "evidence_types": Counter(),
                "no_evidence_count": 0,
                "evidence_hit_counts": defaultdict(list)
            }

            # Analyze each summary file
            for protein in proteins:
                pdb_id = protein['pdb_id']
                chain_id = protein['chain_id']
                file_path = protein['file_path']
                
                # Handle both relative and absolute paths
                full_path = file_path
                if not os.path.isabs(file_path):
                    full_path = os.path.join(base_path, file_path)
                
                # Analyze summary file
                try:
                    self.total_summaries += 1
                    result = analyze_xml_file(full_path, detailed)
                    
                    if result["valid_xml"]:
                        # Track evidence types
                        if "evidence_types" in result:
                            evidence_types = result["evidence_types"]
                            
                            if not evidence_types:
                                # No evidence case
                                self.no_evidence_proteins.append((pdb_id, chain_id))
                                self.evidence_stats["no_evidence"] += 1
                                self.batch_stats[batch_id]["no_evidence_count"] += 1
                            else:
                                # Track individual evidence types
                                for evidence_type in evidence_types:
                                    self.evidence_stats[evidence_type] += 1
                                    self.batch_stats[batch_id]["evidence_types"][evidence_type] += 1
                                    
                                    # Track hit counts
                                    hit_count = result["evidence_counts"].get(evidence_type, 0)
                                    self.evidence_counts[evidence_type].append((pdb_id, chain_id, hit_count))
                                    self.batch_stats[batch_id]["evidence_hit_counts"][evidence_type].append(hit_count)
                                
                                # Track combined evidence types
                                evidence_key = "+".join(sorted(evidence_types))
                                self.evidence_stats[evidence_key] += 1
                                self.batch_stats[batch_id]["evidence_types"][evidence_key] += 1
                    else:
                        # Track parse errors
                        self.parsing_errors.append((pdb_id, chain_id, result.get("error", "Unknown error")))
                        
                except Exception as e:
                    self.logger.error(f"Error processing {pdb_id}_{chain_id}: {str(e)}")
                    self.parsing_errors.append((pdb_id, chain_id, str(e)))

            # Calculate batch statistics
            batch_stats = self.batch_stats[batch_id]

            # Calculate overall evidence percentages
            if batch_stats["analyzed_items"] > 0:
                no_evidence_pct = (batch_stats["no_evidence_count"] / batch_stats["analyzed_items"]) * 100
                batch_stats["no_evidence_percentage"] = no_evidence_pct
                
                # Calculate evidence type percentages
                evidence_percentages = {}
                for evidence_type, count in batch_stats["evidence_types"].items():
                    if "+" not in evidence_type:  # Individual types only
                        percentage = (count / batch_stats["analyzed_items"]) * 100
                        evidence_percentages[evidence_type] = percentage
                
                batch_stats["evidence_percentages"] = evidence_percentages
                
                # Calculate hit count statistics
                hit_stats = {}
                for evidence_type, counts in batch_stats["evidence_hit_counts"].items():
                    if counts:
                        hit_stats[evidence_type] = {
                            "avg": sum(counts) / len(counts),
                            "min": min(counts),
                            "max": max(counts),
                            "median": sorted(counts)[len(counts)//2]
                        }
                
                batch_stats["hit_statistics"] = hit_stats
            # Log statistics
            self.logger.info(f"Completed analysis for batch {batch_id}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error in analyze_batch: {str(e)}", exc_info=True)
            return False
        
    def print_summary(self):
        """Print summary of analysis results"""
        self.logger.info(f"\n{'='*80}\nDomain Summary Analysis Results\n{'='*80}")
        self.logger.info(f"Total summaries analyzed: {self.total_summaries}")
        
        # Evidence type distributions
        self.logger.info("\nEvidence Type Distribution:")
        for evidence_type, count in sorted(self.evidence_stats.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / self.total_summaries) * 100
            self.logger.info(f"  {evidence_type}: {count} ({percentage:.1f}%)")
        
        # Evidence hit count statistics
        self.logger.info("\nHit Count Statistics:")
        for evidence_type in ["chain_blast", "domain_blast", "hhsearch"]:
            if evidence_type in self.evidence_counts and self.evidence_counts[evidence_type]:
                hit_counts = [count for _, _, count in self.evidence_counts[evidence_type]]
                avg_hits = sum(hit_counts) / len(hit_counts)
                max_hits = max(hit_counts)
                min_hits = min(hit_counts)
                self.logger.info(f"  {evidence_type}: avg={avg_hits:.1f}, min={min_hits}, max={max_hits}")
        
        # No evidence proteins
        if self.no_evidence_proteins:
            self.logger.info("\nProteins with NO evidence:")
            for pdb_id, chain_id in self.no_evidence_proteins[:10]:  # Show first 10
                self.logger.info(f"  {pdb_id}_{chain_id}")
            
            if len(self.no_evidence_proteins) > 10:
                self.logger.info(f"  ... and {len(self.no_evidence_proteins) - 10} more")
                
        # Parse errors
        if self.parsing_errors:
            self.logger.info("\nParse errors:")
            for pdb_id, chain_id, error in self.parsing_errors[:10]:  # Show first 10
                self.logger.info(f"  {pdb_id}_{chain_id}: {error}")
            
            if len(self.parsing_errors) > 10:
                self.logger.info(f"  ... and {len(self.parsing_errors) - 10} more")
    
    def write_report(self, output_file):
        """Write detailed analysis report to file"""
        try:
            with open(output_file, 'w') as f:
                f.write(f"Domain Summary Analysis Report\n")
                f.write(f"============================\n\n")
                f.write(f"Total summaries analyzed: {self.total_summaries}\n\n")
                
                # Per-batch statistics
                f.write("Batch Statistics:\n")
                f.write("================\n\n")
                
                for batch_id, stats in self.batch_stats.items():
                    f.write(f"Batch {batch_id} ({stats['name']}):\n")
                    f.write(f"  Analyzed: {stats['analyzed_items']}/{stats['total_items']} proteins\n")
                    
                    if "no_evidence_percentage" in stats:
                        f.write(f"  No evidence: {stats['no_evidence_count']} proteins ({stats['no_evidence_percentage']:.1f}%)\n")
                    
                    if "evidence_percentages" in stats:
                        f.write("  Evidence types:\n")
                        for evidence_type, percentage in sorted(stats["evidence_percentages"].items(), key=lambda x: x[1], reverse=True):
                            count = stats["evidence_types"][evidence_type]
                            f.write(f"    {evidence_type}: {count} proteins ({percentage:.1f}%)\n")
                    
                    if "hit_statistics" in stats:
                        f.write("  Hit counts:\n")
                        for evidence_type, hit_stats in stats["hit_statistics"].items():
                            f.write(f"    {evidence_type}: avg={hit_stats['avg']:.1f}, min={hit_stats['min']}, max={hit_stats['max']}, median={hit_stats['median']}\n")
                    
                    f.write("\n")
                
                # Overall statistics
                f.write("Overall Evidence Type Distribution:\n")
                for evidence_type, count in sorted(self.evidence_stats.items(), key=lambda x: x[1], reverse=True):
                    percentage = (count / self.total_summaries) * 100
                    f.write(f"  {evidence_type}: {count} ({percentage:.1f}%)\n")
                
                # No evidence proteins - full list
                if self.no_evidence_proteins:
                    f.write("\nProteins with NO evidence:\n")
                    for pdb_id, chain_id in self.no_evidence_proteins:
                        f.write(f"  {pdb_id}_{chain_id}\n")
                
                # Parse errors - full list
                if self.parsing_errors:
                    f.write("\nParse errors:\n")
                    for pdb_id, chain_id, error in self.parsing_errors:
                        f.write(f"  {pdb_id}_{chain_id}: {error}\n")
            
            self.logger.info(f"Report written to {output_file}")
            
        except Exception as e:
            self.logger.error(f"Error writing report: {str(e)}")

def main():
    """Main entry point"""
    # Get script directory for proper path handling when run from anywhere
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(script_dir))
    
    parser = argparse.ArgumentParser(description='Analyze ECOD Domain Summaries')
    parser.add_argument('--config', type=str, default=os.path.join(project_root, 'config', 'config.yml'),
                      help='Path to configuration file')
    parser.add_argument('--batch-id', type=int, default=31,
                      help='Batch ID to analyze (default: 31)')
    parser.add_argument('--sample-size', type=int,
                      help='Number of proteins to sample from batch (random sample)')
    parser.add_argument('--output', type=str,
                      help='Output file for detailed report')
    parser.add_argument('--json-output', type=str,
                      help='Output JSON file for structured results')
    parser.add_argument('--file', type=str,
                      help='Analyze a specific file')
    parser.add_argument('--log-file', type=str,
                      help='Log file path')
    parser.add_argument('--detailed', action='store_true',
                      help='Include detailed hit information in analysis')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')

    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)
    
    logger = logging.getLogger("main")
    
    # Ensure all paths are absolute
    if args.config and not os.path.isabs(args.config):
        args.config = os.path.abspath(args.config)
    
    if args.output and not os.path.isabs(args.output):
        args.output = os.path.abspath(args.output)
        
    if args.json_output and not os.path.isabs(args.json_output):
        args.json_output = os.path.abspath(args.json_output)
        
    if args.log_file and not os.path.isabs(args.log_file):
        args.log_file = os.path.abspath(args.log_file)
    
    # Analyze a specific file if requested
    if args.file:
        # Ensure file path is absolute
        if not os.path.isabs(args.file):
            args.file = os.path.abspath(args.file)
            
        logger.info(f"Analyzing file: {args.file}")
        result = analyze_xml_file(args.file, detailed=True)
        
        print(f"\nFile: {os.path.basename(args.file)}")
        print(f"Valid XML: {result['valid_xml']}")
        
        if result["valid_xml"]:
            print(f"PDB ID: {result.get('pdb_id', 'Unknown')}")
            print(f"Chain ID: {result.get('chain_id', 'Unknown')}")
            
            # Evidence summary
            print(f"\nEvidence types: {', '.join(result.get('evidence_types', []))}") 
            
            # Individual evidence types
            if "chain_blast_hits" in result:
                print(f"Chain BLAST hits: {result['chain_blast_hits']}")
            if "domain_blast_hits" in result:
                print(f"Domain BLAST hits: {result['domain_blast_hits']}")
            if "hhsearch_hits" in result:
                print(f"HHSearch hits: {result['hhsearch_hits']}")
            if "domain_count" in result:
                print(f"Domain count: {result['domain_count']}")
            
            # Detailed hit info if available
            for ev_type in ["chain_blast_top_hits", "domain_blast_top_hits", "hhsearch_top_hits"]:
                if ev_type in result and result[ev_type]:
                    print(f"\nTop {ev_type.replace('_top_hits', '')} hits:")
                    for i, hit in enumerate(result[ev_type][:5]):
                        hit_id = hit.get("pdb_id", "") + "_" + hit.get("chain_id", "")
                        if "domain_id" in hit:
                            hit_id = hit["domain_id"]
                        print(f"  Hit {i+1}: {hit_id}")
                        for k, v in hit.items():
                            if k not in ["pdb_id", "chain_id", "domain_id"]:
                                print(f"    {k}: {v}")
            
            # Domain details if available
            if "domain_details" in result and result["domain_details"]:
                print("\nDomains:")
                for i, domain in enumerate(result["domain_details"]):
                    print(f"  Domain {i+1}: {domain.get('id', 'Unknown')}")
                    print(f"    Range: {domain.get('range', 'Unknown')}")
                    if "t_group" in domain and domain["t_group"]:
                        print(f"    Classification: {domain.get('t_group', '')}, {domain.get('h_group', '')}")
        
        if result.get("error"):
            print(f"Error: {result['error']}")
        
        # Also output as JSON if requested
        if args.json_output:
            with open(args.json_output, 'w') as f:
                json.dump(result, f, indent=2)
            logger.info(f"Results written to {args.json_output}")
        
        return 0
    
    # Initialize application context
    logger.info(f"Starting analysis of domain summaries for batch {args.batch_id}")
    context = ApplicationContext(args.config)
    
    # Create analyzer
    analyzer = DomainSummaryAnalyzer(context)
    
    # Run analysis
    success = analyzer.analyze_batch(args.batch_id, args.sample_size, args.detailed)
    
    if success:
        # Print summary
        analyzer.print_summary()
        
        # Write report if requested
        if args.output:
            analyzer.write_report(args.output)
        
        # Write JSON output if requested
        if args.json_output:
            try:
                with open(args.json_output, 'w') as f:
                    json_data = {
                        "summary": {
                            "total_summaries": analyzer.total_summaries,
                            "evidence_stats": dict(analyzer.evidence_stats),
                            "no_evidence_count": len(analyzer.no_evidence_proteins),
                            "parsing_errors_count": len(analyzer.parsing_errors)
                        },
                        "batch_stats": analyzer.batch_stats,
                        "no_evidence_proteins": [{"pdb_id": p[0], "chain_id": p[1]} for p in analyzer.no_evidence_proteins],
                        "parsing_errors": [{"pdb_id": p[0], "chain_id": p[1], "error": p[2]} for p in analyzer.parsing_errors]
                    }
                    json.dump(json_data, f, indent=2)
                logger.info(f"JSON results written to {args.json_output}")
            except Exception as e:
                logger.error(f"Error writing JSON output: {str(e)}")
        
        logger.info(f"Successfully analyzed domain summaries for batch {args.batch_id}")
        return 0
    else:
        logger.error(f"Failed to analyze domain summaries for batch {args.batch_id}")
        return 1

if __name__ == "__main__":
    sys.exit(main())