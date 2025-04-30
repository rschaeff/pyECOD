#!/usr/bin/env python3
"""
batch_analyzer.py - Comprehensive Batch Analysis Tool for pyECOD

This script provides a unified interface for analyzing ECOD batches:
- structure: Analyze PDB structure metadata (resolution, dates, methods)
- coverage: Assess domain coverage and completeness
- quality: Evaluate domain assignments and partitioning quality
- summary: Generate comprehensive batch summary report
- validation: Validate XML structure and content

Each mode provides specific metrics and can output JSON reports for further processing.
"""

import os
import sys
import logging
import argparse
import json
import re
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Set
from collections import Counter, defaultdict
import statistics

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Import core modules
from ecod.core.context import ApplicationContext
from ecod.utils.path_utils import (
    get_standardized_paths,
    resolve_file_path,
    find_files_with_legacy_paths,
    get_file_type_from_path
)
from ecod.utils.xml_utils import element_to_dict, ensure_dict, ensure_list_of_dicts

def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Configure logging with appropriate handlers and format"""
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

class BatchAnalyzer:
    """Comprehensive batch analyzer with multiple analysis capabilities"""
    
    def __init__(self, context: ApplicationContext, batch_id: int = None, batch_path: str = None):
        """Initialize with application context and batch information"""
        self.context = context
        self.logger = logging.getLogger("ecod.batch_analyzer")
        self.batch_id = batch_id
        self.batch_path = batch_path
        self.ref_version = None
        self.batch_name = None
        
        # Initialize batch information if batch ID is provided
        if batch_id:
            self._initialize_from_db()
        elif batch_path:
            self._initialize_from_fs()
        else:
            self.logger.warning("Neither batch_id nor batch_path provided")
    
    def _initialize_from_db(self):
        """Initialize batch information from database using batch ID"""
        if not self.batch_id:
            self.logger.error("No batch ID provided")
            return False
        
        # Query for batch info
        query = """
        SELECT id, batch_name, base_path, ref_version, total_items
        FROM ecod_schema.batch
        WHERE id = %s
        """
        
        batch_info = self.context.db.execute_dict_query(query, (self.batch_id,))
        if not batch_info:
            self.logger.error(f"Batch {self.batch_id} not found in database")
            return False
        
        self.batch_path = batch_info[0]['base_path']
        self.ref_version = batch_info[0]['ref_version']
        self.batch_name = batch_info[0]['batch_name']
        self.total_items = batch_info[0]['total_items']
        
        self.logger.info(f"Initialized batch {self.batch_id} ({self.batch_name}) with ref {self.ref_version}")
        return True
    
    def _initialize_from_fs(self):
        """Initialize batch information from filesystem using batch path"""
        if not self.batch_path:
            self.logger.error("No batch path provided")
            return False
        
        if not os.path.exists(self.batch_path):
            self.logger.error(f"Batch path doesn't exist: {self.batch_path}")
            return False
        
        # Try to determine reference version from paths
        domains_dir = os.path.join(self.batch_path, "domains")
        if os.path.exists(domains_dir):
            # Look for domain files to extract reference version
            for file in os.listdir(domains_dir):
                if file.endswith(".domain_summary.xml") or file.endswith(".domains.xml"):
                    parts = file.split(".")
                    if len(parts) >= 3:
                        # Reference version is typically the second-to-last part before file type
                        potential_ref = parts[-3] if parts[-2] in ["domain_summary", "domains"] else parts[-2]
                        if re.match(r'develop\d+|release\d+', potential_ref):
                            self.ref_version = potential_ref
                            break
        
        # If we couldn't find reference version, assume default
        if not self.ref_version:
            self.ref_version = "develop291"  # Default reference version
            self.logger.warning(f"Could not determine reference version, using default: {self.ref_version}")
        
        # Extract batch name from path
        self.batch_name = os.path.basename(self.batch_path)
        
        # Try to estimate total items from directories
        fastas_dir = os.path.join(self.batch_path, "fastas")
        if os.path.exists(fastas_dir):
            fasta_files = [f for f in os.listdir(fastas_dir) if f.endswith('.fa')]
            self.total_items = len(fasta_files)
        else:
            self.total_items = 0
        
        self.logger.info(f"Initialized batch {self.batch_name} with ref {self.ref_version}")
        return True
    
    def analyze_structure_metadata(self) -> Dict[str, Any]:
        """Analyze structure metadata (resolution, dates, methods)"""
        self.logger.info(f"Analyzing structure metadata for batch {self.batch_id or self.batch_name}")
        
        result = {
            "batch_id": self.batch_id,
            "batch_name": self.batch_name,
            "total_structures": 0,
            "resolution_stats": {},
            "deposition_dates": {},
            "experimental_methods": {},
            "pdb_id_prefixes": {}
        }
        
        if self.batch_id:
            # Use database to get structure metadata
            query = """
            SELECT 
                p.pdb_id, 
                ps.resolution, 
                ps.experimental_method, 
                ps.deposition_date
            FROM 
                ecod_schema.process_status eps
            JOIN 
                ecod_schema.protein ep ON eps.protein_id = ep.id
            LEFT JOIN 
                pdb_analysis.protein p ON ep.source_id = p.source_id
            LEFT JOIN 
                pdb_analysis.protein_structure ps ON p.id = ps.protein_id
            WHERE 
                eps.batch_id = %s
            """
            
            structures = self.context.db.execute_dict_query(query, (self.batch_id,))
            result["total_structures"] = len(structures)
            
            # Process structure data
            resolutions = []
            deposition_years = Counter()
            methods = Counter()
            id_prefixes = Counter()
            
            for structure in structures:
                # Resolution statistics
                if structure.get('resolution') is not None:
                    resolutions.append(float(structure['resolution']))
                
                # Deposition dates
                if structure.get('deposition_date'):
                    year = structure['deposition_date'].year
                    deposition_years[year] += 1
                
                # Experimental methods
                if structure.get('experimental_method'):
                    methods[structure['experimental_method']] += 1
                
                # PDB ID prefixes (1xxx, 2xxx, etc.)
                if structure.get('pdb_id'):
                    prefix = structure['pdb_id'][0] + 'xxx'
                    id_prefixes[prefix] += 1
            
            # Calculate resolution statistics
            if resolutions:
                result["resolution_stats"] = {
                    "count": len(resolutions),
                    "min": min(resolutions),
                    "max": max(resolutions),
                    "mean": statistics.mean(resolutions),
                    "median": statistics.median(resolutions),
                    "bins": self._bin_resolutions(resolutions)
                }
            
            # Sort and convert counters to dictionaries
            result["deposition_dates"] = dict(sorted(deposition_years.items()))
            result["experimental_methods"] = dict(methods.most_common())
            result["pdb_id_prefixes"] = dict(id_prefixes.most_common())
            
            # Calculate date ranges
            if result["deposition_dates"]:
                years = list(result["deposition_dates"].keys())
                result["date_ranges"] = {
                    "min_year": min(years),
                    "max_year": max(years),
                    "span_years": max(years) - min(years) + 1
                }
                
                # Group into ranges for reporting
                result["year_groups"] = {
                    "pre_2010": sum(result["deposition_dates"].get(year, 0) for year in range(1900, 2010)),
                    "2010_2015": sum(result["deposition_dates"].get(year, 0) for year in range(2010, 2016)),
                    "2016_2020": sum(result["deposition_dates"].get(year, 0) for year in range(2016, 2021)),
                    "2021_present": sum(result["deposition_dates"].get(year, 0) for year in range(2021, 2026))
                }
        else:
            # For filesystem mode, we need to infer structure info from filenames
            # This is less accurate but can still provide some information
            fastas_dir = os.path.join(self.batch_path, "fastas")
            if os.path.exists(fastas_dir):
                pdb_ids = []
                for file in os.listdir(fastas_dir):
                    if file.endswith('.fa'):
                        parts = file.split('.')
                        pdb_chain = parts[0]  # Format: pdbid_chain
                        if '_' in pdb_chain:
                            pdb_id = pdb_chain.split('_')[0]
                            pdb_ids.append(pdb_id)
                
                result["total_structures"] = len(set(pdb_ids))
                
                # Count PDB ID prefixes
                id_prefixes = Counter()
                for pdb_id in pdb_ids:
                    if len(pdb_id) >= 1:
                        prefix = pdb_id[0] + 'xxx'
                        id_prefixes[prefix] += 1
                
                result["pdb_id_prefixes"] = dict(id_prefixes.most_common())
        
        # Log summary information
        self.logger.info(f"Found {result['total_structures']} structures")
        if "resolution_stats" in result and result["resolution_stats"]:
            self.logger.info(f"Resolution range: {result['resolution_stats']['min']:.2f} - {result['resolution_stats']['max']:.2f} Å")
            self.logger.info(f"Median resolution: {result['resolution_stats']['median']:.2f} Å")
        
        if "date_ranges" in result:
            self.logger.info(f"Deposition years: {result['date_ranges']['min_year']} - {result['date_ranges']['max_year']}")
            if "year_groups" in result:
                self.logger.info(f"Recent structures (2021+): {result['year_groups']['2021_present']} structures")
        
        return result
    
    def _bin_resolutions(self, resolutions: List[float]) -> Dict[str, int]:
        """Bin resolutions into categories for reporting"""
        bins = {
            "ultrahigh_<1A": 0,
            "excellent_1-1.5A": 0,
            "high_1.5-2A": 0,
            "medium_2-2.5A": 0,
            "low_2.5-3A": 0,
            "poor_>3A": 0
        }
        
        for res in resolutions:
            if res < 1.0:
                bins["ultrahigh_<1A"] += 1
            elif res < 1.5:
                bins["excellent_1-1.5A"] += 1
            elif res < 2.0:
                bins["high_1.5-2A"] += 1
            elif res < 2.5:
                bins["medium_2-2.5A"] += 1
            elif res < 3.0:
                bins["low_2.5-3A"] += 1
            else:
                bins["poor_>3A"] += 1
        
        return bins
    
    def analyze_domain_coverage(self) -> Dict[str, Any]:
        """Analyze domain coverage and completeness"""
        self.logger.info(f"Analyzing domain coverage for batch {self.batch_id or self.batch_name}")
        
        result = {
            "batch_id": self.batch_id,
            "batch_name": self.batch_name,
            "total_proteins": 0,
            "domain_coverage": {},
            "file_types": {},
            "process_stages": {},
            "domain_statistics": {}
        }
        
        if self.batch_id:
            # Use database to get coverage information
            
            # Get process status stats
            process_query = """
            SELECT 
                current_stage, status, COUNT(*) as count
            FROM 
                ecod_schema.process_status
            WHERE 
                batch_id = %s
            GROUP BY 
                current_stage, status
            ORDER BY 
                current_stage, status
            """
            
            process_stats = self.context.db.execute_dict_query(process_query, (self.batch_id,))
            
            # Process the stats into a structured format
            stages = defaultdict(dict)
            for row in process_stats:
                stage = row['current_stage']
                status = row['status']
                count = row['count']
                stages[stage][status] = count
            
            result["process_stages"] = dict(stages)
            
            # Get file type counts
            file_query = """
            SELECT 
                file_type, 
                COUNT(*) as total_count,
                SUM(CASE WHEN file_exists = TRUE THEN 1 ELSE 0 END) as exists_count
            FROM 
                ecod_schema.process_file pf
            JOIN 
                ecod_schema.process_status ps ON pf.process_id = ps.id
            WHERE 
                ps.batch_id = %s
            GROUP BY 
                file_type
            ORDER BY 
                file_type
            """
            
            file_stats = self.context.db.execute_dict_query(file_query, (self.batch_id,))
            
            # Process file stats
            file_types = {}
            for row in file_stats:
                file_type = row['file_type']
                total = row['total_count']
                exists = row['exists_count']
                file_types[file_type] = {
                    "total": total,
                    "exists": exists,
                    "completion": (exists / total * 100) if total > 0 else 0
                }
            
            result["file_types"] = file_types
            
            # Calculate domain statistics if domain files exist
            if "domain_partition" in file_types and file_types["domain_partition"]["exists"] > 0:
                domain_query = """
                SELECT
                    COUNT(d.id) as total_domains,
                    COUNT(DISTINCT d.protein_id) as proteins_with_domains,
                    COUNT(DISTINCT CASE WHEN t_group IS NOT NULL THEN protein_id END) as classified_proteins
                FROM
                    ecod_schema.process_status ps
                JOIN
                    ecod_schema.protein p ON ps.protein_id = p.id
                LEFT JOIN
                    ecod_schema.domain d ON p.id = d.protein_id
                WHERE
                    ps.batch_id = %s
                """
                
                domain_stats = self.context.db.execute_dict_query(domain_query, (self.batch_id,))
                
                if domain_stats:
                    result["domain_statistics"] = {
                        "total_domains": domain_stats[0]["total_domains"],
                        "proteins_with_domains": domain_stats[0]["proteins_with_domains"],
                        "classified_proteins": domain_stats[0]["classified_proteins"]
                    }
            
            # Get total protein count
            count_query = "SELECT COUNT(*) FROM ecod_schema.process_status WHERE batch_id = %s"
            count_result = self.context.db.execute_query(count_query, (self.batch_id,))
            result["total_proteins"] = count_result[0][0] if count_result else 0
            
            # Calculate domain summary completion percentage
            if "domain_summary" in file_types:
                result["domain_coverage"] = {
                    "completed": file_types["domain_summary"]["exists"],
                    "total": result["total_proteins"],
                    "percentage": (file_types["domain_summary"]["exists"] / result["total_proteins"] * 100) 
                        if result["total_proteins"] > 0 else 0
                }
        else:
            # For filesystem mode, scan the directories
            domains_dir = os.path.join(self.batch_path, "domains")
            fastas_dir = os.path.join(self.batch_path, "fastas")
            
            # Get total proteins from FASTA files
            total_proteins = 0
            if os.path.exists(fastas_dir):
                fasta_files = [f for f in os.listdir(fastas_dir) if f.endswith('.fa')]
                total_proteins = len(fasta_files)
            
            result["total_proteins"] = total_proteins
            
            # Count domain summary and partition files
            if os.path.exists(domains_dir):
                summary_files = [f for f in os.listdir(domains_dir) if "domain_summary" in f]
                partition_files = [f for f in os.listdir(domains_dir) if f.endswith(".domains.xml")]
                
                result["file_types"] = {
                    "domain_summary": {
                        "total": total_proteins,
                        "exists": len(summary_files),
                        "completion": (len(summary_files) / total_proteins * 100) if total_proteins > 0 else 0
                    },
                    "domain_partition": {
                        "total": total_proteins,
                        "exists": len(partition_files),
                        "completion": (len(partition_files) / total_proteins * 100) if total_proteins > 0 else 0
                    }
                }
                
                result["domain_coverage"] = {
                    "completed": len(summary_files),
                    "total": total_proteins,
                    "percentage": (len(summary_files) / total_proteins * 100) if total_proteins > 0 else 0
                }
        
        # Log summary information
        self.logger.info(f"Total proteins: {result['total_proteins']}")
        
        if "domain_coverage" in result:
            coverage = result["domain_coverage"]
            self.logger.info(f"Domain coverage: {coverage['completed']}/{coverage['total']} proteins ({coverage['percentage']:.1f}%)")
        
        if "file_types" in result:
            for file_type, stats in result["file_types"].items():
                if "completion" in stats:
                    self.logger.info(f"{file_type}: {stats['exists']}/{stats['total']} files ({stats['completion']:.1f}%)")
        
        return result
    
    def analyze_quality(self, sample_size: int = 50) -> Dict[str, Any]:
        """Evaluate domain assignment quality and file structure validity"""
        self.logger.info(f"Analyzing quality for batch {self.batch_id or self.batch_name}")
        
        result = {
            "batch_id": self.batch_id,
            "batch_name": self.batch_name,
            "total_analyzed": 0,
            "valid_xml": 0,
            "summary_stats": {},
            "validation_issues": [],
            "evidence_types": {},
            "domain_counts": {}
        }
        
        # Get a sample of domain summary files to analyze
        summary_files = []
        
        if self.batch_id:
            # Use database to get domain summary files
            query = """
            SELECT 
                p.pdb_id, p.chain_id, pf.file_path
            FROM 
                ecod_schema.process_file pf
            JOIN 
                ecod_schema.process_status ps ON pf.process_id = ps.id
            JOIN 
                ecod_schema.protein p ON ps.protein_id = p.id
            WHERE 
                ps.batch_id = %s
                AND pf.file_type = 'domain_summary'
                AND pf.file_exists = TRUE
            ORDER BY 
                RANDOM()
            LIMIT %s
            """
            
            rows = self.context.db.execute_dict_query(query, (self.batch_id, sample_size))
            
            for row in rows:
                file_path = row['file_path']
                if not os.path.isabs(file_path):
                    file_path = os.path.join(self.batch_path, file_path)
                
                summary_files.append({
                    'pdb_id': row['pdb_id'],
                    'chain_id': row['chain_id'],
                    'path': file_path
                })
        else:
            # For filesystem mode, scan the domains directory
            domains_dir = os.path.join(self.batch_path, "domains")
            if os.path.exists(domains_dir):
                # Get all summary files
                all_files = []
                for file in os.listdir(domains_dir):
                    if "domain_summary" in file:
                        file_path = os.path.join(domains_dir, file)
                        # Extract PDB ID and chain ID from filename
                        parts = file.split('.')
                        pdb_chain = parts[0]  # Format: pdbid_chain
                        if '_' in pdb_chain:
                            pdb_id, chain_id = pdb_chain.split('_')
                            all_files.append({
                                'pdb_id': pdb_id,
                                'chain_id': chain_id,
                                'path': file_path
                            })
                
                # Select random sample
                import random
                if len(all_files) > sample_size:
                    summary_files = random.sample(all_files, sample_size)
                else:
                    summary_files = all_files
        
        result["total_analyzed"] = len(summary_files)
        
        # Analyze each summary file
        valid_count = 0
        evidence_types = Counter()
        domain_counts = Counter()
        issues = []
        
        for file_info in summary_files:
            file_path = file_info['path']
            pdb_id = file_info['pdb_id']
            chain_id = file_info['chain_id']
            
            try:
                # Parse XML
                tree = ET.parse(file_path)
                root = tree.getroot()
                
                # Check if valid structure
                if root.tag != "domain_summ_doc":
                    issues.append(f"{pdb_id}_{chain_id}: Invalid root element '{root.tag}'")
                    continue
                
                valid_count += 1
                
                # Check for evidence types
                has_chain_blast = False
                has_domain_blast = False
                has_hhsearch = False
                
                # Check chain blast evidence
                chain_blast = root.find("chain_blast_evidence")
                if chain_blast is not None:
                    hits = chain_blast.findall(".//hit")
                    if hits:
                        has_chain_blast = True
                        evidence_types["chain_blast"] += 1
                
                # Check domain blast evidence
                domain_blast = root.find("domain_blast_evidence")
                if domain_blast is not None:
                    hits = domain_blast.findall(".//hit")
                    if hits:
                        has_domain_blast = True
                        evidence_types["domain_blast"] += 1
                
                # Check HHSearch evidence
                hhsearch = root.find("hhsearch_evidence")
                if hhsearch is not None:
                    hits = hhsearch.findall(".//hh_hit")
                    if not hits:
                        hits = hhsearch.findall(".//hit")
                    if hits:
                        has_hhsearch = True
                        evidence_types["hhsearch"] += 1
                
                # Track combined evidence types
                evidence_key = []
                if has_chain_blast:
                    evidence_key.append("chain_blast")
                if has_domain_blast:
                    evidence_key.append("domain_blast")
                if has_hhsearch:
                    evidence_key.append("hhsearch")
                
                if evidence_key:
                    evidence_types["+".join(sorted(evidence_key))] += 1
                else:
                    evidence_types["no_evidence"] += 1
                
                # Check domain count
                domain_elems = root.findall(".//domain")
                domain_count = len(domain_elems)
                domain_counts[domain_count] += 1
                
                # Check for missing elements
                if chain_blast is None:
                    issues.append(f"{pdb_id}_{chain_id}: Missing chain_blast_evidence section")
                
                if domain_blast is None:
                    issues.append(f"{pdb_id}_{chain_id}: Missing domain_blast_evidence section")
                
                if hhsearch is None and not "blast_only" in file_path:
                    issues.append(f"{pdb_id}_{chain_id}: Missing hhsearch_evidence section in full pipeline file")
                
                domain_suggestions = root.find("domain_suggestions")
                if domain_suggestions is None:
                    issues.append(f"{pdb_id}_{chain_id}: Missing domain_suggestions section")
                
            except ET.ParseError as e:
                issues.append(f"{pdb_id}_{chain_id}: XML parsing error - {str(e)}")
            except Exception as e:
                issues.append(f"{pdb_id}_{chain_id}: Error - {str(e)}")
        
        # Compile results
        result["valid_xml"] = valid_count
        
        if result["total_analyzed"] > 0:
            result["summary_stats"] = {
                "valid_percentage": (valid_count / result["total_analyzed"] * 100),
                "issues_percentage": (len(issues) / result["total_analyzed"] * 100)
            }
        
        result["validation_issues"] = issues[:10]  # Limit to first 10 issues
        result["validation_issue_count"] = len(issues)
        
        # Convert counters to dictionaries
        result["evidence_types"] = dict(evidence_types.most_common())
        result["domain_counts"] = dict(sorted(domain_counts.items()))
        
        # Calculate additional statistics
        if domain_counts:
            domain_values = []
            for count, frequency in domain_counts.items():
                domain_values.extend([count] * frequency)
            
            if domain_values:
                result["domain_stats"] = {
                    "mean": statistics.mean(domain_values),
                    "median": statistics.median(domain_values),
                    "min": min(domain_values),
                    "max": max(domain_values)
                }
        
        # Log summary information
        self.logger.info(f"Analyzed {result['total_analyzed']} summary files")
        self.logger.info(f"Valid XML: {valid_count}/{result['total_analyzed']} ({result['summary_stats'].get('valid_percentage', 0):.1f}%)")
        
        if "evidence_types" in result:
            self.logger.info("Evidence types:")
            for evidence_type, count in result["evidence_types"].items():
                percentage = (count / result["total_analyzed"] * 100) if result["total_analyzed"] > 0 else 0
                self.logger.info(f"  {evidence_type}: {count} ({percentage:.1f}%)")
        
        if issues:
            self.logger.info(f"Found {len(issues)} validation issues")
            for i, issue in enumerate(result["validation_issues"]):
                self.logger.info(f"  {i+1}. {issue}")
            
            if len(issues) > 10:
                self.logger.info(f"  ... and {len(issues) - 10} more issues")
        
        return result
    
    def generate_comprehensive_report(self) -> Dict[str, Any]:
        """Generate a comprehensive report by combining all analysis types"""
        self.logger.info(f"Generating comprehensive report for batch {self.batch_id or self.batch_name}")
        
        # Run all analysis types
        structure_data = self.analyze_structure_metadata()
        coverage_data = self.analyze_domain_coverage()
        quality_data = self.analyze_quality(sample_size=50)
        
        # Combine into a comprehensive report
        report = {
            "batch_id": self.batch_id,
            "batch_name": self.batch_name,
            "ref_version": self.ref_version,
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "structure_metadata": structure_data,
            "domain_coverage": coverage_data,
            "quality_assessment": quality_data,
            "summary": {}
        }
        
        # Create executive summary
        summary = {}
        
        # Structure summary
        if "total_structures" in structure_data:
            summary["total_structures"] = structure_data["total_structures"]
        
        # Add median resolution if available
        if "resolution_stats" in structure_data and "median" in structure_data["resolution_stats"]:
            summary["median_resolution"] = structure_data["resolution_stats"]["median"]
        
        # Add recent structures count if available
        if "year_groups" in structure_data and "2021_present" in structure_data["year_groups"]:
            summary["recent_structures"] = structure_data["year_groups"]["2021_present"]
        
        # Coverage summary
        if "total_proteins" in coverage_data:
            summary["total_proteins"] = coverage_data["total_proteins"]
        
        if "domain_coverage" in coverage_data and "percentage" in coverage_data["domain_coverage"]:
            summary["domain_coverage_pct"] = coverage_data["domain_coverage"]["percentage"]
        
        # Quality summary
        if "valid_xml" in quality_data and "total_analyzed" in quality_data and quality_data["total_analyzed"] > 0:
            summary["xml_validity_pct"] = (quality_data["valid_xml"] / quality_data["total_analyzed"] * 100)
        
        if "domain_stats" in quality_data and "mean" in quality_data["domain_stats"]:
            summary["avg_domains_per_chain"] = quality_data["domain_stats"]["mean"]
        
        report["summary"] = summary
        
        # Log executive summary
        self.logger.info("Executive Summary:")
        for key, value in summary.items():
            if isinstance(value, float):
                self.logger.info(f"  {key}: {value:.2f}")
            else:
                self.logger.info(f"  {key}: {value}")
        
        return report


def structure_mode(args: argparse.Namespace, context: ApplicationContext) -> int:
    """Run structure metadata analysis mode"""
    logger = logging.getLogger("ecod.batch_analyzer.structure")
    
    # Initialize analyzer
    analyzer = BatchAnalyzer(context, args.batch_id, args.batch_path)
    
    # Run analysis
    results = analyzer.analyze_structure_metadata()
    
    # Save results if output file specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results written to {args.output}")
        except Exception as e:
            logger.error(f"Error writing results: {str(e)}")
    
    return 0

def coverage_mode(args: argparse.Namespace, context: ApplicationContext) -> int:
    """Run domain coverage analysis mode"""
    logger = logging.getLogger("ecod.batch_analyzer.coverage")
    
    # Initialize analyzer
    analyzer = BatchAnalyzer(context, args.batch_id, args.batch_path)
    
    # Run analysis
    results = analyzer.analyze_domain_coverage()
    
    # Save results if output file specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results written to {args.output}")
        except Exception as e:
            logger.error(f"Error writing results: {str(e)}")
    
    return 0

def quality_mode(args: argparse.Namespace, context: ApplicationContext) -> int:
    """Run quality assessment mode"""
    logger = logging.getLogger("ecod.batch_analyzer.quality")
    
    # Initialize analyzer
    analyzer = BatchAnalyzer(context, args.batch_id, args.batch_path)
    
    # Run analysis
    results = analyzer.analyze_quality(args.sample_size)
    
    # Save results if output file specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results written to {args.output}")
        except Exception as e:
            logger.error(f"Error writing results: {str(e)}")
    
    return 0

def summary_mode(args: argparse.Namespace, context: ApplicationContext) -> int:
    """Run comprehensive summary mode"""
    logger = logging.getLogger("ecod.batch_analyzer.summary")
    
    # Initialize analyzer
    analyzer = BatchAnalyzer(context, args.batch_id, args.batch_path)
    
    # Run analysis
    results = analyzer.generate_comprehensive_report()
    
    # Save results if output file specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Comprehensive report written to {args.output}")
        except Exception as e:
            logger.error(f"Error writing results: {str(e)}")
    
    return 0

def validation_mode(args: argparse.Namespace, context: ApplicationContext) -> int:
    """Run XML structure validation mode"""
    logger = logging.getLogger("ecod.batch_analyzer.validation")
    
    # Initialize analyzer with larger sample size for validation
    analyzer = BatchAnalyzer(context, args.batch_id, args.batch_path)
    
    # Run analysis with larger sample size
    results = analyzer.analyze_quality(args.sample_size or 100)
    
    # Filter to just validation-related fields
    validation_results = {
        "batch_id": results["batch_id"],
        "batch_name": results["batch_name"],
        "total_analyzed": results["total_analyzed"],
        "valid_xml": results["valid_xml"],
        "valid_percentage": results["summary_stats"].get("valid_percentage", 0),
        "validation_issues": results["validation_issues"],
        "validation_issue_count": results["validation_issue_count"],
        "evidence_types": results["evidence_types"]
    }
    
    # Save results if output file specified
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(validation_results, f, indent=2)
            logger.info(f"Validation results written to {args.output}")
        except Exception as e:
            logger.error(f"Error writing results: {str(e)}")
    
    return 0

def main():
    """Main entry point"""
    # Create top-level parser
    parser = argparse.ArgumentParser(description='Batch Analyzer - Comprehensive Batch Analysis Tool for pyECOD')
    parser.add_argument('--config', type=str, default='config/config.yml', help='Path to configuration file')
    parser.add_argument('--log-file', type=str, help='Log file path')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    
    # Create subparsers for different modes
    subparsers = parser.add_subparsers(dest='mode', help='Analysis mode')
    
    # Structure metadata mode
    structure_parser = subparsers.add_parser('structure', help='Analyze PDB structure metadata')
    structure_parser.add_argument('--batch-id', type=int, help='Batch ID (for database mode)')
    structure_parser.add_argument('--batch-path', type=str, help='Batch path (for filesystem mode)')
    structure_parser.add_argument('--output', type=str, help='Output JSON file')
    
    # Domain coverage mode
    coverage_parser = subparsers.add_parser('coverage', help='Analyze domain coverage and completeness')
    coverage_parser.add_argument('--batch-id', type=int, help='Batch ID (for database mode)')
    coverage_parser.add_argument('--batch-path', type=str, help='Batch path (for filesystem mode)')
    coverage_parser.add_argument('--output', type=str, help='Output JSON file')
    
    # Quality assessment mode
    quality_parser = subparsers.add_parser('quality', help='Evaluate domain assignment quality')
    quality_parser.add_argument('--batch-id', type=int, help='Batch ID (for database mode)')
    quality_parser.add_argument('--batch-path', type=str, help='Batch path (for filesystem mode)')
    quality_parser.add_argument('--sample-size', type=int, default=50, help='Number of files to analyze')
    quality_parser.add_argument('--output', type=str, help='Output JSON file')
    
    # Summary mode
    summary_parser = subparsers.add_parser('summary', help='Generate comprehensive batch summary')
    summary_parser.add_argument('--batch-id', type=int, help='Batch ID (for database mode)')
    summary_parser.add_argument('--batch-path', type=str, help='Batch path (for filesystem mode)')
    summary_parser.add_argument('--output', type=str, help='Output JSON file')
    
    # Validation mode
    validation_parser = subparsers.add_parser('validation', help='Validate XML structure and content')
    validation_parser.add_argument('--batch-id', type=int, help='Batch ID (for database mode)')
    validation_parser.add_argument('--batch-path', type=str, help='Batch path (for filesystem mode)')
    validation_parser.add_argument('--sample-size', type=int, help='Number of files to analyze')
    validation_parser.add_argument('--output', type=str, help='Output JSON file')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger("ecod.batch_analyzer")
    
    # Validate arguments - need either batch_id or batch_path
    if not hasattr(args, 'batch_id') and not hasattr(args, 'batch_path'):
        logger.error("Either --batch-id or --batch-path must be specified")
        return 1
    
    if args.batch_id is None and args.batch_path is None:
        logger.error("Either --batch-id or --batch-path must be provided")
        return 1
    
    # Initialize application context
    try:
        context = ApplicationContext(args.config)
    except Exception as e:
        logger.error(f"Error initializing application context: {str(e)}")
        return 1
    
    # Run appropriate mode
    if args.mode == 'structure':
        return structure_mode(args, context)
    elif args.mode == 'coverage':
        return coverage_mode(args, context)
    elif args.mode == 'quality':
        return quality_mode(args, context)
    elif args.mode == 'summary':
        return summary_mode(args, context)
    elif args.mode == 'validation':
        return validation_mode(args, context)
    else:
        logger.error(f"Unknown mode: {args.mode}")
        parser.print_help()
        return 1

if __name__ == "__main__":
    sys.exit(main())
