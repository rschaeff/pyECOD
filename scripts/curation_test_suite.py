#!/usr/bin/env python3
"""
Curation Test Suite for Domain Partition Algorithm Evaluation

This script creates test sets from manually curated proteins and evaluates
improvements in the domain partition algorithm. Specifically designed to test
improvements in:

1. DISCONTINUOUS DOMAIN HANDLING - Better treatment of multi-segment domains
2. REFERENCE DOMAIN COVERAGE CUTOFFS - More explicit coverage thresholds  
3. BOUNDARY ACCURACY - Critical for domain definition correctness
4. FRAGMENT/PEPTIDE DETECTION - Distinguishing short peptides from domains

Key workflow:
- Extract ~110 curated proteins as test set
- QUICK TEST on sample (10-15 proteins) to verify algorithm improvements work
- If improvements confirmed, snapshot current partition_domains as baseline
- Re-run improved algorithm on test proteins  
- Compare boundary accuracy and fragment detection
- Validate improvements before production deployment

Usage Examples:
    # Safe development workflow
    python curation_test_suite.py extract --config config.yml --min-confidence 4
    python curation_test_suite.py quicktest --test-set-id 1 --sample-size 10
    # Only proceed if quicktest shows improvements!
    python curation_test_suite.py snapshot --test-set-id 1
    python curation_test_suite.py rerun --test-set-id 1 --algorithm-version improved_v1
    python curation_test_suite.py evaluate --test-set-id 1
    
    # After validation, promote to production
    python curation_test_suite.py promote --test-set-id 1 --algorithm-version improved_v1 --force
"""

import os
import sys
import logging
import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, asdict
from collections import defaultdict

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from ecod.config import ConfigManager
from ecod.core.context import ApplicationContext
from ecod.db import DBManager
from ecod.pipelines.domain_analysis.partition import DomainPartitionService
from ecod.models.pipeline.partition import DomainPartitionResult

class AlgorithmIterationType(Enum):
    """Types of algorithm iterations"""
    BASELINE = "baseline"
    COVERAGE_FOCUSED = "coverage_focused"
    EVIDENCE_WEIGHTING = "evidence_weighting"
    ARCHITECTURAL_TRANSFER = "architectural_transfer"
    HYBRID = "hybrid"
    EXPERIMENTAL = "experimental"

@dataclass
class AlgorithmVersion:
    """Comprehensive algorithm version specification"""

    # Identity
    version_id: str  # e.g., "v2.1_chain_blast_priority"
    iteration_type: AlgorithmIterationType
    parent_version: Optional[str] = None

    # Metadata
    name: str = ""
    description: str = ""
    created_by: str = ""
    created_at: datetime = field(default_factory=datetime.now)

    # Core algorithm configuration
    partition_options: Dict[str, Any] = field(default_factory=dict)

    # Evidence processing configuration
    evidence_weights: Dict[str, float] = field(default_factory=dict)
    evidence_filters: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    # Coverage and boundary settings
    coverage_thresholds: Dict[str, float] = field(default_factory=dict)
    boundary_settings: Dict[str, Any] = field(default_factory=dict)

    # Performance and behavioral flags
    performance_settings: Dict[str, Any] = field(default_factory=dict)
    behavioral_flags: Dict[str, bool] = field(default_factory=dict)

    # Testing and deployment
    test_results: List[Dict[str, Any]] = field(default_factory=list)
    deployment_status: str = "development"  # development, testing, production, deprecated

    # Change tracking
    changes_from_parent: List[str] = field(default_factory=list)
    improvement_targets: List[str] = field(default_factory=list)
    known_limitations: List[str] = field(default_factory=list)

def create_baseline_algorithm() -> AlgorithmVersion:
    """Original baseline algorithm (pre-improvements)"""
    return AlgorithmVersion(
        version_id="v1.0_baseline",
        iteration_type=AlgorithmIterationType.BASELINE,
        name="Original Baseline",
        description="Original domain partition algorithm before systematic improvements",

        partition_options={
            "validation_level": "normal",
            "min_domain_size": 20,
            "overlap_threshold": 0.3,
            "merge_gap_tolerance": 20
        },

        evidence_weights={
            "hhsearch": 3.0,
            "domain_blast": 2.5,
            "chain_blast": 2.0,
            "blast": 1.5
        },

        coverage_thresholds={
            "min_reference_coverage": 0.5,  # Very lenient
            "strict_reference_coverage": 0.8,
            "partial_coverage_threshold": 0.2
        },

        behavioral_flags={
            "prefer_hhsearch_classification": True,
            "extend_to_reference_size": False
        },

        improvement_targets=[],
        known_limitations=[
            "Poor fragment detection",
            "Over-segments discontinuous domains",
            "Inconsistent boundary determination",
            "High false positive rate for short peptides"
        ]
    )

def create_coverage_focused_algorithm() -> AlgorithmVersion:
    """First iteration: Reference coverage focused"""
    return AlgorithmVersion(
        version_id="v1.1_coverage_focused",
        iteration_type=AlgorithmIterationType.COVERAGE_FOCUSED,
        parent_version="v1.0_baseline",
        name="Coverage-Focused Algorithm",
        description="Improved reference domain alignment coverage requirements to reduce fragment assignments",

        partition_options={
            "validation_level": "normal",
            "min_domain_size": 20,
            "overlap_threshold": 0.3,
            "merge_gap_tolerance": 20
        },

        evidence_weights={
            "hhsearch": 3.0,
            "domain_blast": 2.5,
            "chain_blast": 2.0,
            "blast": 1.5
        },

        coverage_thresholds={
            "min_reference_coverage": 0.7,  # INCREASED - key improvement
            "strict_reference_coverage": 0.9,  # INCREASED
            "partial_coverage_threshold": 0.3,  # INCREASED
            "extend_to_reference_size": True,  # NEW
            "reference_size_tolerance": 0.15,  # NEW
            "use_ungapped_coverage": True  # NEW
        },

        behavioral_flags={
            "prefer_hhsearch_classification": True,
            "extend_to_reference_size": True  # NEW
        },

        changes_from_parent=[
            "Increased min_reference_coverage from 0.5 to 0.7",
            "Increased strict_reference_coverage from 0.8 to 0.9",
            "Added reference size extension capability",
            "Added ungapped coverage calculation",
            "Improved coverage validation in analyzer"
        ],

        improvement_targets=[
            "Reduce false positive domain assignments in short peptides",
            "Improve fragment vs domain classification accuracy",
            "Ensure adequate reference coverage for domain assignments"
        ],

        known_limitations=[
            "Still struggles with discontinuous domains",
            "Coverage alone cannot determine gap inclusion",
            "May be too restrictive for some valid domains",
            "Doesn't address evidence precedence issues"
        ]
    )

def create_chain_blast_priority_algorithm() -> AlgorithmVersion:
    """Second iteration: Chain BLAST architectural transfer"""
    return AlgorithmVersion(
        version_id="v2.0_chain_blast_priority",
        iteration_type=AlgorithmIterationType.ARCHITECTURAL_TRANSFER,
        parent_version="v1.1_coverage_focused",
        name="Chain BLAST Architectural Transfer",
        description="Prioritize chain BLAST evidence to leverage validated domain architectures and improve discontinuous domain handling",

        partition_options={
            "validation_level": "normal",
            "min_domain_size": 20,
            "overlap_threshold": 0.4,  # INCREASED - more lenient with chain BLAST
            "merge_gap_tolerance": 30   # INCREASED - better gap handling
        },

        evidence_weights={
            "chain_blast": 3.5,      # HIGHEST PRIORITY - key change
            "hhsearch": 2.0,          # REDUCED from 3.0
            "domain_blast": 2.2,      # SLIGHTLY REDUCED
            "blast": 2.0              # INCREASED from 1.5
        },

        coverage_thresholds={
            # Keep coverage improvements from v1.1
            "min_reference_coverage": 0.7,
            "strict_reference_coverage": 0.9,
            "partial_coverage_threshold": 0.3,
            "extend_to_reference_size": True,
            "reference_size_tolerance": 0.15,
            "use_ungapped_coverage": True
        },

        behavioral_flags={
            "prefer_hhsearch_classification": False,  # CHANGED - remove HHSearch bias
            "prefer_chain_blast_architecture": True,  # NEW
            "extend_to_reference_size": True,
            "use_architectural_context": True  # NEW
        },

        changes_from_parent=[
            "Increased chain_blast weight from 2.0 to 3.5 (highest priority)",
            "Reduced hhsearch weight from 3.0 to 2.0",
            "Disabled prefer_hhsearch_classification",
            "Added chain BLAST boundary protection",
            "Increased overlap and gap tolerance",
            "Added architectural context usage"
        ],

        improvement_targets=[
            "Improve discontinuous domain detection and boundaries",
            "Leverage architectural knowledge from reference proteins",
            "Reduce over-fragmentation of natural domain units",
            "Better handle domain architecture variations"
        ],

        known_limitations=[
            "Depends on quality of reference protein architectures",
            "May miss novel domain arrangements not in reference set",
            "Could be less precise than HHSearch for boundary details",
            "Needs comprehensive chain BLAST database"
        ]
    )

def create_chain_blast_only_algorithm() -> AlgorithmVersion:
    """Experimental: Chain BLAST only"""
    return AlgorithmVersion(
        version_id="v2.1_chain_blast_only",
        iteration_type=AlgorithmIterationType.EXPERIMENTAL,
        parent_version="v2.0_chain_blast_priority",
        name="Chain BLAST Only (Experimental)",
        description="Test version using only chain BLAST evidence to isolate architectural transfer effects",

        partition_options={
            "validation_level": "lenient",  # More lenient for testing
            "min_domain_size": 20,
            "overlap_threshold": 0.4,
            "merge_gap_tolerance": 30
        },

        evidence_weights={
            "chain_blast": 3.0,      # Only evidence type
            "hhsearch": 0.0,          # DISABLED
            "domain_blast": 0.0,      # DISABLED
            "blast": 2.0              # Generic fallback
        },

        coverage_thresholds={
            "min_reference_coverage": 0.5,  # REDUCED - more lenient
            "strict_reference_coverage": 0.8,
            "partial_coverage_threshold": 0.2,  # REDUCED
            "extend_to_reference_size": True,
            "use_ungapped_coverage": True
        },

        behavioral_flags={
            "use_chain_blast": True,
            "use_domain_blast": False,  # DISABLED
            "use_hhsearch": False,      # DISABLED
            "prefer_chain_blast_architecture": True,
            "extend_to_reference_size": True
        },

        changes_from_parent=[
            "Disabled all evidence except chain BLAST",
            "Reduced coverage requirements for testing",
            "More lenient validation for chain-blast-only scenarios"
        ],

        improvement_targets=[
            "Isolate chain BLAST architectural transfer effects",
            "Test discontinuous domain improvements without other evidence interference",
            "Validate architectural knowledge transfer hypothesis"
        ],

        known_limitations=[
            "May miss domains with poor chain BLAST coverage",
            "Reduced classification precision without HHSearch",
            "Experimental only - not suitable for production",
            "Limited by chain BLAST database completeness"
        ]
    )

class AlgorithmVersionManager:
    """Simple version of algorithm version manager"""

    def __init__(self, db_manager):
        self.db = db_manager
        self.logger = logging.getLogger(__name__)

        # Registry of available algorithms
        self.algorithms = {
            "v1.0_baseline": create_baseline_algorithm(),
            "v1.1_coverage_focused": create_coverage_focused_algorithm(),
            "v2.0_chain_blast_priority": create_chain_blast_priority_algorithm(),
            "v2.1_chain_blast_only": create_chain_blast_only_algorithm()
        }

    def get_algorithm(self, version_id: str) -> Optional[AlgorithmVersion]:
        """Get algorithm version by ID"""
        return self.algorithms.get(version_id)

    def list_algorithms(self, deployment_status: Optional[str] = None) -> List[AlgorithmVersion]:
        """List available algorithms"""
        algorithms = list(self.algorithms.values())

        if deployment_status:
            algorithms = [a for a in algorithms if a.deployment_status == deployment_status]

        return sorted(algorithms, key=lambda a: a.created_at)

@dataclass
class CurationDecision:
    """Manual curation decision for a protein"""
    protein_id: int
    source_id: str
    pdb_id: str
    chain_id: str
    has_domain: bool
    domain_assigned_correctly: Optional[bool]
    boundaries_correct: Optional[bool]
    is_fragment: bool
    is_repeat_protein: bool
    confidence_level: int
    primary_evidence_type: Optional[str]
    reference_domain_id: Optional[str] 
    reference_pdb_id: Optional[str]
    reference_chain_id: Optional[str]
    evidence_confidence: Optional[float]
    evidence_evalue: Optional[float]
    curator_name: str
    session_id: int
    created_at: datetime


@dataclass 
class PartitionResult:
    """Automated domain partition result"""
    protein_id: int
    source_id: str
    pdb_id: str
    chain_id: str
    is_classified: bool
    is_peptide: bool
    domain_count: int
    domains: List[Dict[str, Any]]
    coverage: float
    sequence_length: int
    confidence_scores: List[float]
    algorithm_version: str
    processing_timestamp: datetime


@dataclass
class TestSet:
    """A test set of curated proteins"""
    test_set_id: int
    name: str
    description: str
    created_at: datetime
    protein_count: int
    curator_breakdown: Dict[str, int]
    decision_breakdown: Dict[str, int]
    proteins: List[CurationDecision]


@dataclass
class ComparisonMetrics:
    """Metrics comparing automated results to manual curation"""
    test_set_id: int
    algorithm_version: str
    
    # Primary metrics (focus areas)
    domain_presence_accuracy: float  # Correctly predicted has_domain vs not
    fragment_detection_accuracy: float  # Correctly identified fragments/peptides
    boundary_agreement_rate: float  # % with correct boundaries 
    
    # Discontinuous domain metrics (new improvement focus)
    discontinuous_domain_detection: float  # How well we detect multi-segment domains
    discontinuous_boundary_accuracy: float  # Boundary accuracy for discontinuous domains
    coverage_cutoff_effectiveness: float  # How well coverage cutoffs work
    
    # Boundary-specific metrics (major focus area)
    exact_boundary_matches: float  # Exact start/end position matches
    boundary_tolerance_5: float   # Within 5 residues
    boundary_tolerance_10: float  # Within 10 residues
    boundary_over_segmentation: float  # Rate of splitting domains incorrectly
    boundary_under_segmentation: float # Rate of merging domains incorrectly
    
    # Fragment/peptide detection (major focus area)
    peptide_precision: float      # Precision for peptide calls
    peptide_recall: float         # Recall for peptide calls
    fragment_vs_domain_accuracy: float  # Correctly distinguishing fragments from domains
    
    # Secondary metrics
    domain_count_accuracy: float  # How often we get the right number of domains
    classification_agreement_rate: float  # % with correct classification (when available)
    
    # Detailed breakdowns
    confusion_matrix: Dict[str, int]
    boundary_error_distribution: Dict[str, List[int]]  # Error sizes for boundary mismatches
    improvement_cases: List[Dict[str, Any]]
    regression_cases: List[Dict[str, Any]]
    discontinuous_cases: List[Dict[str, Any]]  # Specific discontinuous domain cases


class CurationTestManager:
    """Manager for curation-based testing of domain partition improvements"""
    
    def __init__(self, context: ApplicationContext):
        self.context = context
        self.db = context.db
        self.logger = logging.getLogger("curation_test")
        
        # Initialize test set tracking table
        self._ensure_test_tables()
    
    def _ensure_test_tables(self):
        """Create tables for tracking test sets and results"""
        
        # Test sets table
        self.db.execute_query("""
            CREATE TABLE IF NOT EXISTS pdb_analysis.curation_test_sets (
                id SERIAL PRIMARY KEY,
                name VARCHAR(255) NOT NULL,
                description TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                protein_count INTEGER,
                curator_breakdown JSONB,
                decision_breakdown JSONB,
                extraction_criteria JSONB
            )
        """)
        
        # Test proteins table
        self.db.execute_query("""
            CREATE TABLE IF NOT EXISTS pdb_analysis.curation_test_proteins (
                id SERIAL PRIMARY KEY,
                test_set_id INTEGER REFERENCES pdb_analysis.curation_test_sets(id),
                protein_id INTEGER REFERENCES pdb_analysis.protein(id),
                curation_decision_id INTEGER REFERENCES pdb_analysis.curation_decision(id),
                source_id VARCHAR(20) NOT NULL,
                pdb_id VARCHAR(4) NOT NULL,
                chain_id VARCHAR(10) NOT NULL,
                has_domain BOOLEAN,
                is_fragment BOOLEAN,
                confidence_level INTEGER,
                curator_name VARCHAR(100),
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        # Algorithm results table  
        self.db.execute_query("""
            CREATE TABLE IF NOT EXISTS pdb_analysis.curation_test_results (
                id SERIAL PRIMARY KEY,
                test_set_id INTEGER REFERENCES pdb_analysis.curation_test_sets(id),
                protein_id INTEGER REFERENCES pdb_analysis.protein(id),
                algorithm_version VARCHAR(50) NOT NULL,
                is_classified BOOLEAN,
                is_peptide BOOLEAN,
                domain_count INTEGER,
                coverage FLOAT,
                sequence_length INTEGER,
                domains JSONB,
                confidence_scores JSONB,
                processing_time FLOAT,
                result_file_path TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        # Comparison metrics table
        self.db.execute_query("""
            CREATE TABLE IF NOT EXISTS pdb_analysis.curation_test_metrics (
                id SERIAL PRIMARY KEY,
                test_set_id INTEGER REFERENCES pdb_analysis.curation_test_sets(id),
                algorithm_version VARCHAR(50) NOT NULL,
                metric_name VARCHAR(100) NOT NULL,
                metric_value FLOAT,
                metric_data JSONB,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
    def extract_curated_proteins(self, 
                                min_confidence: int = 3,
                                exclude_fragments: bool = True,
                                min_curator_count: int = 1,
                                max_proteins: Optional[int] = None) -> TestSet:
        """Extract proteins with manual curation decisions as test set"""
        
        self.logger.info("Extracting curated proteins for test set")
        
        # Build query for curated proteins
        query = """
        SELECT DISTINCT
            cd.protein_id,
            cd.source_id,
            p.pdb_id,
            p.chain_id,
            cd.has_domain,
            cd.domain_assigned_correctly,
            cd.boundaries_correct,
            cd.is_fragment,
            cd.is_repeat_protein,
            cd.confidence_level,
            cd.primary_evidence_type,
            cd.reference_domain_id,
            cd.reference_pdb_id,
            cd.reference_chain_id,
            cd.evidence_confidence,
            cd.evidence_evalue,
            cs.curator_name,
            cd.session_id,
            cd.created_at
        FROM pdb_analysis.curation_decision cd
        JOIN pdb_analysis.curation_session cs ON cd.session_id = cs.id
        JOIN pdb_analysis.protein p ON cd.protein_id = p.id
        WHERE cd.confidence_level >= %s
        """
        
        params = [min_confidence]
        
        if exclude_fragments:
            query += " AND cd.is_fragment = FALSE"
            
        # Only include committed sessions
        query += " AND cs.status = 'committed'"
        
        # Order by confidence and creation date
        query += " ORDER BY cd.confidence_level DESC, cd.created_at DESC"
        
        if max_proteins:
            query += f" LIMIT {max_proteins}"
            
        results = self.db.execute_dict_query(query, params)
        
        if not results:
            raise ValueError("No curated proteins found matching criteria")
            
        # Convert to CurationDecision objects
        proteins = []
        for row in results:
            decision = CurationDecision(
                protein_id=row['protein_id'],
                source_id=row['source_id'],
                pdb_id=row['pdb_id'],
                chain_id=row['chain_id'],
                has_domain=row['has_domain'],
                domain_assigned_correctly=row['domain_assigned_correctly'],
                boundaries_correct=row['boundaries_correct'],
                is_fragment=row['is_fragment'],
                is_repeat_protein=row['is_repeat_protein'],
                confidence_level=row['confidence_level'],
                primary_evidence_type=row['primary_evidence_type'],
                reference_domain_id=row['reference_domain_id'],
                reference_pdb_id=row['reference_pdb_id'],
                reference_chain_id=row['reference_chain_id'],
                evidence_confidence=row['evidence_confidence'],
                evidence_evalue=row['evidence_evalue'],
                curator_name=row['curator_name'],
                session_id=row['session_id'],
                created_at=row['created_at']
            )
            proteins.append(decision)
            
        # Calculate breakdown statistics
        curator_breakdown = defaultdict(int)
        decision_breakdown = defaultdict(int)
        
        for protein in proteins:
            curator_breakdown[protein.curator_name] += 1
            
            if protein.has_domain:
                if protein.domain_assigned_correctly:
                    decision_breakdown['correct_domain'] += 1
                else:
                    decision_breakdown['incorrect_domain'] += 1
            else:
                decision_breakdown['no_domain'] += 1
                
            if protein.is_fragment:
                decision_breakdown['fragment'] += 1
                
        # Create test set
        test_set = TestSet(
            test_set_id=0,  # Will be set when saved
            name=f"curated_proteins_{datetime.now().strftime('%Y%m%d_%H%M')}",
            description=f"Extracted {len(proteins)} curated proteins (min_confidence={min_confidence})",
            created_at=datetime.now(),
            protein_count=len(proteins),
            curator_breakdown=dict(curator_breakdown),
            decision_breakdown=dict(decision_breakdown),
            proteins=proteins
        )
        
        # Save to database
        test_set_id = self._save_test_set(test_set)
        test_set.test_set_id = test_set_id
        
        self.logger.info(f"Created test set {test_set_id} with {len(proteins)} proteins")
        self.logger.info(f"Curator breakdown: {dict(curator_breakdown)}")
        self.logger.info(f"Decision breakdown: {dict(decision_breakdown)}")
        
        return test_set
    
    def _save_test_set(self, test_set: TestSet) -> int:
        """Save test set to database"""
        
        # Insert test set record
        query = """
        INSERT INTO pdb_analysis.curation_test_sets 
        (name, description, protein_count, curator_breakdown, decision_breakdown, extraction_criteria)
        VALUES (%s, %s, %s, %s, %s, %s)
        RETURNING id
        """
        
        extraction_criteria = {
            'min_confidence': 3,  # Would need to pass these as parameters
            'exclude_fragments': True,
            'max_proteins': len(test_set.proteins)
        }
        
        result = self.db.execute_query(query, (
            test_set.name,
            test_set.description, 
            test_set.protein_count,
            json.dumps(test_set.curator_breakdown),
            json.dumps(test_set.decision_breakdown),
            json.dumps(extraction_criteria)
        ))
        
        test_set_id = result[0][0]
        
        # Insert protein records
        for protein in test_set.proteins:
            protein_query = """
            INSERT INTO pdb_analysis.curation_test_proteins
            (test_set_id, protein_id, source_id, pdb_id, chain_id, 
             has_domain, is_fragment, confidence_level, curator_name)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
            """
            
            self.db.execute_query(protein_query, (
                test_set_id,
                protein.protein_id,
                protein.source_id,
                protein.pdb_id,
                protein.chain_id,
                protein.has_domain,
                protein.is_fragment,
                protein.confidence_level,
                protein.curator_name
            ))
            
        return test_set_id
    
    def get_test_set(self, test_set_id: int) -> TestSet:
        """Retrieve test set from database"""
        
        # Get test set metadata
        query = """
        SELECT name, description, created_at, protein_count, 
               curator_breakdown, decision_breakdown
        FROM pdb_analysis.curation_test_sets
        WHERE id = %s
        """
        
        results = self.db.execute_dict_query(query, (test_set_id,))
        if not results:
            raise ValueError(f"Test set {test_set_id} not found")
            
        metadata = results[0]
        
        # Get proteins in test set
        protein_query = """
        SELECT ttp.protein_id, ttp.source_id, ttp.pdb_id, ttp.chain_id,
               ttp.has_domain, ttp.is_fragment, ttp.confidence_level, 
               ttp.curator_name, cd.*
        FROM pdb_analysis.curation_test_proteins ttp
        LEFT JOIN pdb_analysis.curation_decision cd ON ttp.protein_id = cd.protein_id
        WHERE ttp.test_set_id = %s
        ORDER BY ttp.pdb_id, ttp.chain_id
        """
        
        protein_results = self.db.execute_dict_query(protein_query, (test_set_id,))
        
        proteins = []
        for row in protein_results:
            decision = CurationDecision(
                protein_id=row['protein_id'],
                source_id=row['source_id'],
                pdb_id=row['pdb_id'],
                chain_id=row['chain_id'],
                has_domain=row['has_domain'],
                domain_assigned_correctly=row.get('domain_assigned_correctly'),
                boundaries_correct=row.get('boundaries_correct'),
                is_fragment=row['is_fragment'],
                is_repeat_protein=row.get('is_repeat_protein', False),
                confidence_level=row['confidence_level'],
                primary_evidence_type=row.get('primary_evidence_type'),
                reference_domain_id=row.get('reference_domain_id'),
                reference_pdb_id=row.get('reference_pdb_id'),
                reference_chain_id=row.get('reference_chain_id'),
                evidence_confidence=row.get('evidence_confidence'),
                evidence_evalue=row.get('evidence_evalue'),
                curator_name=row['curator_name'],
                session_id=row.get('session_id', 0),
                created_at=row.get('created_at', datetime.now())
            )
            proteins.append(decision)
            
        return TestSet(
            test_set_id=test_set_id,
            name=metadata['name'],
            description=metadata['description'],
            created_at=metadata['created_at'],
            protein_count=metadata['protein_count'],
            curator_breakdown=metadata['curator_breakdown'],
            decision_breakdown=metadata['decision_breakdown'],
            proteins=proteins
        )
    
    def rerun_algorithm(self, test_set_id: int, algorithm_version: str = "improved_v1") -> int:
        """Re-run domain partition algorithm on test set proteins"""
        
        self.logger.info(f"Re-running algorithm on test set {test_set_id}")
        
        # Get test set
        test_set = self.get_test_set(test_set_id)
        
        # Create domain partition service
        service = DomainPartitionService(self.context)
        
        # Process each protein
        results_saved = 0
        
        for protein in test_set.proteins:
            try:
                # Find domain summary files for this protein
                summary_paths = self._find_domain_summaries(protein.pdb_id, protein.chain_id)
                
                if not summary_paths:
                    self.logger.warning(f"No domain summaries found for {protein.pdb_id}_{protein.chain_id}")
                    continue
                    
                # Run partition algorithm
                result = self._run_partition_for_protein(service, protein, summary_paths)
                
                if result:
                    # Save results to database
                    self._save_algorithm_result(test_set_id, protein, result, algorithm_version)
                    results_saved += 1
                    
                    self.logger.info(f"Processed {protein.pdb_id}_{protein.chain_id}: "
                                   f"classified={result.is_classified}, domains={result.domain_count}")
                
            except Exception as e:
                self.logger.error(f"Error processing {protein.pdb_id}_{protein.chain_id}: {e}")
                
        self.logger.info(f"Saved {results_saved}/{len(test_set.proteins)} algorithm results")
        return results_saved
    
    def _find_domain_summaries(self, pdb_id: str, chain_id: str) -> List[str]:
        """Find domain summary files for a protein using the proper ecod_schema"""

        # Now that we have ecod_schema, use the proper database query
        query = """
        SELECT DISTINCT pf.file_path, b.base_path, ps.batch_id
        FROM ecod_schema.process_file pf
        JOIN ecod_schema.process_status ps ON pf.process_id = ps.id
        JOIN ecod_schema.protein p ON ps.protein_id = p.id
        JOIN ecod_schema.batch b ON ps.batch_id = b.id
        WHERE p.pdb_id = %s AND p.chain_id = %s
          AND pf.file_type IN ('domain_summary', 'blast_only_summary')
          AND pf.file_exists = TRUE
        ORDER BY ps.batch_id DESC
        LIMIT 5
        """
        
        results = self.db.execute_dict_query(query, (pdb_id, chain_id))
        
        summary_paths = []
        for row in results:
            full_path = os.path.join(row['base_path'], row['file_path'])
            if os.path.exists(full_path):
                summary_paths.append(full_path)
                
        return summary_paths
    
    def _run_partition_for_protein(self, service: DomainPartitionService, 
                                  protein: CurationDecision,
                                  summary_paths: List[str]) -> Optional[PartitionResult]:
        """Run domain partition algorithm for a single protein"""

        # Use the most recent summary file
        summary_path = summary_paths[0] if summary_paths else None

        # Create a temporary output directory
        temp_output = f"/tmp/curation_test_{protein.pdb_id}_{protein.chain_id}"
        os.makedirs(temp_output, exist_ok=True)

        try:
            # Run partition algorithm
            result = service.partition_protein(
                pdb_id=protein.pdb_id,
                chain_id=protein.chain_id,
                summary_path=summary_path,
                output_dir=temp_output
            )

            if not result.success:
                self.logger.warning(f"Partition failed for {protein.pdb_id}_{protein.chain_id}: {result.error}")
                return None

            # Extract domain information with CORRECT attribute names
            domains = []
            confidence_scores = []

            for domain in result.domains:
                # Use the correct attribute names for DomainModel
                domain_dict = {
                    'range': domain.range,
                    'start_pos': domain.start,      # FIXED: use 'start' not 'start_pos'
                    'end_pos': domain.end,          # FIXED: use 'end' not 'end_pos'
                    'source': domain.source,
                    'source_id': domain.source_id,
                    'confidence': domain.confidence,
                    't_group': domain.t_group,
                    'h_group': domain.h_group,
                    'x_group': domain.x_group,
                    'a_group': domain.a_group
                }
                domains.append(domain_dict)
                
                if domain.confidence:
                    confidence_scores.append(domain.confidence)
                    
            return PartitionResult(
                protein_id=protein.protein_id,
                source_id=protein.source_id,
                pdb_id=protein.pdb_id,
                chain_id=protein.chain_id,
                is_classified=result.is_classified,
                is_peptide=result.is_peptide,
                domain_count=len(result.domains),
                domains=domains,
                coverage=result.coverage,
                sequence_length=result.sequence_length,
                confidence_scores=confidence_scores,
                algorithm_version="improved_v1",
                processing_timestamp=datetime.now()
            )
            
        finally:
            # Clean up temporary files
            import shutil
            if os.path.exists(temp_output):
                shutil.rmtree(temp_output)
    
    def _save_algorithm_result(self, test_set_id: int, protein: CurationDecision,
                              result: PartitionResult, algorithm_version: str):
        """Save algorithm result to database"""
        
        query = """
        INSERT INTO pdb_analysis.curation_test_results
        (test_set_id, protein_id, algorithm_version, is_classified, is_peptide,
         domain_count, coverage, sequence_length, domains, confidence_scores)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        
        self.db.execute_query(query, (
            test_set_id,
            result.protein_id,
            algorithm_version,
            result.is_classified,
            result.is_peptide,
            result.domain_count,
            result.coverage,
            result.sequence_length,
            json.dumps(result.domains),
            json.dumps(result.confidence_scores)
        ))
    
    def evaluate_improvements(self, test_set_id: int, 
                            baseline_version: str = "original",
                            improved_version: str = "improved_v1") -> ComparisonMetrics:
        """Evaluate improvements by comparing algorithm versions against manual curation"""
        
        self.logger.info(f"Evaluating improvements for test set {test_set_id}")
        
        # Get test set and results
        test_set = self.get_test_set(test_set_id)
        baseline_results = self._get_algorithm_results(test_set_id, baseline_version)
        improved_results = self._get_algorithm_results(test_set_id, improved_version)
        
        # Calculate metrics
        metrics = self._calculate_comparison_metrics(
            test_set, baseline_results, improved_results, improved_version
        )
        
        # Save metrics to database
        self._save_comparison_metrics(test_set_id, metrics)
        
        return metrics
    
    def create_baseline_snapshot(self, test_set_id: int, 
                                snapshot_name: str = None) -> str:
        """Create a baseline snapshot of current partition_domains for test proteins"""
        
        if not snapshot_name:
            snapshot_name = f"baseline_{datetime.now().strftime('%Y%m%d_%H%M')}"
            
        self.logger.info(f"Creating baseline snapshot: {snapshot_name}")
        
        # Get test set proteins
        test_set = self.get_test_set(test_set_id)
        protein_ids = [p.protein_id for p in test_set.proteins]
        
        if not protein_ids:
            raise ValueError("No proteins in test set")
            
        # Create snapshot table if it doesn't exist
        self.db.execute_query(f"""
            CREATE TABLE IF NOT EXISTS pdb_analysis.partition_domains_baseline (
                LIKE pdb_analysis.partition_domains INCLUDING ALL
            )
        """)
        
        # Add metadata columns if they don't exist
        self.db.execute_query("""
            ALTER TABLE pdb_analysis.partition_domains_baseline 
            ADD COLUMN IF NOT EXISTS snapshot_name VARCHAR(100),
            ADD COLUMN IF NOT EXISTS snapshot_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        """)
        
        # Copy current results for test proteins
        placeholders = ','.join(['%s'] * len(protein_ids))
        copy_query = f"""
        INSERT INTO pdb_analysis.partition_domains_baseline 
        SELECT pd.*, %s, CURRENT_TIMESTAMP
        FROM pdb_analysis.partition_domains pd
        JOIN pdb_analysis.partition_proteins pp ON pd.protein_id = pp.id
        WHERE pp.id IN ({placeholders})
        """
        
        self.db.execute_query(copy_query, [snapshot_name] + protein_ids)
        
        # Record snapshot metadata
        self.db.execute_query("""
            INSERT INTO pdb_analysis.curation_test_metrics
            (test_set_id, algorithm_version, metric_name, metric_value, metric_data)
            VALUES (%s, %s, %s, %s, %s)
        """, (
            test_set_id,
            "baseline_snapshot",
            "snapshot_created", 
            len(protein_ids),
            json.dumps({
                "snapshot_name": snapshot_name,
                "protein_count": len(protein_ids),
                "created_at": datetime.now().isoformat()
            })
        ))
        
        self.logger.info(f"Created baseline snapshot for {len(protein_ids)} proteins")
        return snapshot_name
    
    def quick_test_improvements(self, test_set_id: int, 
                              algorithm_version: str = "improved_v1",
                              sample_size: int = 10) -> Dict[str, Any]:
        """Quick test of improvements without formal snapshotting"""
        
        self.logger.info(f"Quick testing improvements on {sample_size} proteins from test set {test_set_id}")
        
        # Get test set
        test_set = self.get_test_set(test_set_id)
        
        # Take a sample for quick testing
        import random
        sample_proteins = random.sample(test_set.proteins, min(sample_size, len(test_set.proteins)))
        
        # Create domain partition service
        service = DomainPartitionService(self.context)
        
        results = {
            'total_tested': len(sample_proteins),
            'successful_runs': 0,
            'algorithm_errors': 0,
            'missing_summaries': 0,
            'comparison_results': [],
            'sample_improvements': 0,
            'sample_regressions': 0,
            'fragment_detection_changes': 0,
            'boundary_changes': 0
        }
        
        for protein in sample_proteins:
            try:
                # Find domain summary
                summary_paths = self._find_domain_summaries(protein.pdb_id, protein.chain_id)
                
                if not summary_paths:
                    self.logger.warning(f"No domain summaries found for {protein.pdb_id}_{protein.chain_id}")
                    results['missing_summaries'] += 1
                    continue
                
                # Get current baseline from partition_domains (without formal snapshot)
                current_baseline = self._get_current_partition_result(protein)
                
                # Run improved algorithm
                improved_result = self._run_partition_for_protein(service, protein, summary_paths)
                
                if not improved_result:
                    results['algorithm_errors'] += 1
                    continue
                    
                results['successful_runs'] += 1
                
                # Quick comparison
                comparison = self._quick_compare_results(current_baseline, improved_result, protein)
                results['comparison_results'].append(comparison)
                
                # Track improvement types
                if comparison['overall_improvement'] > 0.1:
                    results['sample_improvements'] += 1
                elif comparison['overall_improvement'] < -0.1:
                    results['sample_regressions'] += 1
                    
                if comparison['fragment_detection_changed']:
                    results['fragment_detection_changes'] += 1
                    
                if comparison['boundary_changed']:
                    results['boundary_changes'] += 1
                    
                self.logger.info(f"Tested {protein.pdb_id}_{protein.chain_id}: "
                               f"improvement={comparison['overall_improvement']:.3f}")
                
            except Exception as e:
                self.logger.error(f"Error testing {protein.pdb_id}_{protein.chain_id}: {e}")
                results['algorithm_errors'] += 1
        
        # Calculate summary statistics
        if results['comparison_results']:
            improvements = [r['overall_improvement'] for r in results['comparison_results']]
            results['avg_improvement'] = sum(improvements) / len(improvements)
            results['positive_improvements'] = len([i for i in improvements if i > 0])
            results['negative_improvements'] = len([i for i in improvements if i < 0])
        
        return results
    
    def _get_current_partition_result(self, protein: CurationDecision) -> Optional[PartitionResult]:
        """Get current partition result from partition_domains table"""
        
        # Find the partition_proteins record
        pp_query = """
        SELECT pp.id, pp.is_classified, pp.coverage, pp.sequence_length
        FROM pdb_analysis.partition_proteins pp
        JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
        WHERE p.id = %s
        ORDER BY pp."timestamp" DESC
        LIMIT 1
        """
        
        pp_result = self.db.execute_dict_query(pp_query, (protein.protein_id,))
        if not pp_result:
            return None
            
        pp_record = pp_result[0]
        partition_protein_id = pp_record['id']
        
        # Get domains
        domain_query = """
        SELECT * FROM pdb_analysis.partition_domains
        WHERE protein_id = %s
        ORDER BY domain_number
        """
        
        domain_results = self.db.execute_dict_query(domain_query, (partition_protein_id,))
        
        # Convert to PartitionResult format
        domains = []
        confidence_scores = []
        
        for domain in domain_results:
            domain_dict = {
                'range': domain['range'],
                'start_pos': domain['start_pos'],
                'end_pos': domain['end_pos'],
                'source': domain['source'],
                'source_id': domain['source_id'],
                'confidence': domain['confidence'],
                't_group': domain['t_group'],
                'h_group': domain['h_group'],
                'x_group': domain['x_group'],
                'a_group': domain['a_group'],
                'pdb_range': domain.get('pdb_range'),
                'length': domain.get('length', 0)
            }
            domains.append(domain_dict)
            
            if domain['confidence']:
                confidence_scores.append(domain['confidence'])
        
        return PartitionResult(
            protein_id=protein.protein_id,
            source_id=protein.source_id,
            pdb_id=protein.pdb_id,
            chain_id=protein.chain_id,
            is_classified=pp_record['is_classified'] or False,
            is_peptide=len(domains) == 0 and pp_record['sequence_length'] and pp_record['sequence_length'] < 50,
            domain_count=len(domains),
            domains=domains,
            coverage=pp_record['coverage'] or 0.0,
            sequence_length=pp_record['sequence_length'] or 0,
            confidence_scores=confidence_scores,
            algorithm_version="current",
            processing_timestamp=datetime.now()
        )
    
    def _quick_compare_results(self, baseline: Optional[PartitionResult], 
                             improved: PartitionResult, 
                             manual: CurationDecision) -> Dict[str, Any]:
        """Quick comparison between baseline and improved results"""
        
        comparison = {
            'protein_id': manual.protein_id,
            'pdb_id': manual.pdb_id,
            'chain_id': manual.chain_id,
            'baseline_available': baseline is not None,
            'overall_improvement': 0.0,
            'fragment_detection_changed': False,
            'boundary_changed': False,
            'domain_count_changed': False,
            'details': {}
        }
        
        if not baseline:
            # No baseline to compare against, just evaluate improved vs manual
            improved_score = self._calculate_boundary_fragment_score(improved, manual)
            comparison['overall_improvement'] = improved_score  # Assume baseline was 0
            comparison['details']['improved_score'] = improved_score
            comparison['details']['note'] = 'no_baseline_comparison'
            return comparison
        
        # Calculate scores
        baseline_score = self._calculate_boundary_fragment_score(baseline, manual)
        improved_score = self._calculate_boundary_fragment_score(improved, manual)
        
        comparison['overall_improvement'] = improved_score - baseline_score
        comparison['details']['baseline_score'] = baseline_score
        comparison['details']['improved_score'] = improved_score
        
        # Check specific changes
        if baseline.is_peptide != improved.is_peptide:
            comparison['fragment_detection_changed'] = True
            comparison['details']['fragment_change'] = {
                'baseline_peptide': baseline.is_peptide,
                'improved_peptide': improved.is_peptide,
                'manual_fragment': manual.is_fragment
            }
        
        if baseline.domain_count != improved.domain_count:
            comparison['domain_count_changed'] = True
            comparison['details']['domain_count_change'] = {
                'baseline_count': baseline.domain_count,
                'improved_count': improved.domain_count
            }
        
        # Check for boundary changes (simplified)
        if baseline.domains and improved.domains:
            baseline_ranges = {d.get('range', '') for d in baseline.domains}
            improved_ranges = {d.get('range', '') for d in improved.domains}
            
            if baseline_ranges != improved_ranges:
                comparison['boundary_changed'] = True
                comparison['details']['boundary_change'] = {
                    'baseline_ranges': list(baseline_ranges),
                    'improved_ranges': list(improved_ranges)
                }
        
        return comparison 

    def _get_baseline_results(self, test_set_id: int, snapshot_name: str = None
    ) -> Dict[int, PartitionResult]:
        """Get baseline results from snapshot or current partition_domains"""

        if snapshot_name:
            # Get from snapshot
            query = """
            SELECT pp.id as protein_id, pd.*
            FROM pdb_analysis.partition_domains_baseline pd
            JOIN pdb_analysis.partition_proteins pp ON pd.protein_id = pp.id
            WHERE pd.snapshot_name = %s
            """
            params = [snapshot_name]
        else:
            # Get current results for test proteins
            test_set = self.get_test_set(test_set_id)
            protein_ids = [p.protein_id for p in test_set.proteins]
            placeholders = ','.join(['%s'] * len(protein_ids))

            query = f"""
            SELECT pp.id as protein_id, pd.*
            FROM pdb_analysis.partition_domains pd
            JOIN pdb_analysis.partition_proteins pp ON pd.protein_id = pp.id
            JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
            WHERE p.id IN ({placeholders})
            """
            params = protein_ids

        results = self.db.execute_dict_query(query, params)

        # Group by protein and convert to PartitionResult
        baseline_results = {}
        protein_domains = defaultdict(list)

        for row in results:
            protein_id = row['protein_id']
            domain_info = {
                'range': row['range'],
                'start_pos': row['start_pos'],
                'end_pos': row['end_pos'],
                'source': row['source'],
                'source_id': row['source_id'],
                'confidence': row['confidence'],
                't_group': row['t_group'],
                'h_group': row['h_group'],
                'x_group': row['x_group'],
                'a_group': row['a_group'],
                'pdb_range': row.get('pdb_range'),
                'length': row.get('length', 0)
            }
            protein_domains[protein_id].append(domain_info)

        # Convert to PartitionResult objects
        for protein_id, domains in protein_domains.items():
            # Calculate metrics
            is_classified = len(domains) > 0
            is_peptide = False  # Would need to determine from data
            domain_count = len(domains)

            confidence_scores = [d['confidence'] for d in domains if d['confidence']]
            total_length = sum(d['length'] for d in domains if d['length'])

            result = PartitionResult(
                protein_id=protein_id,
                source_id="",  # Will be filled from protein lookup if needed
                pdb_id="",
                chain_id="",
                is_classified=is_classified,
                is_peptide=is_peptide,
                domain_count=domain_count,
                domains=domains,
                coverage=0.0,  # Would calculate if sequence length available
                sequence_length=0,
                confidence_scores=confidence_scores,
                algorithm_version="baseline",
                processing_timestamp=datetime.now()
            )
            baseline_results[protein_id] = result

        return baseline_results

    def _get_algorithm_results(self, test_set_id: int, algorithm_version: str) -> Dict[int, PartitionResult]:
        """Get algorithm results for a specific version"""

        query = """
        SELECT protein_id, is_classified, is_peptide, domain_count,
               coverage, sequence_length, domains, confidence_scores
        FROM pdb_analysis.curation_test_results
        WHERE test_set_id = %s AND algorithm_version = %s
        """

        results = self.db.execute_dict_query(query, (test_set_id, algorithm_version))

        algorithm_results = {}
        for row in results:
            result = PartitionResult(
                protein_id=row['protein_id'],
                source_id="",  # Not needed for comparison
                pdb_id="",     # Not needed for comparison
                chain_id="",   # Not needed for comparison
                is_classified=row['is_classified'],
                is_peptide=row['is_peptide'],
                domain_count=row['domain_count'],
                domains=row['domains'] or [],
                coverage=row['coverage'] or 0.0,
                sequence_length=row['sequence_length'] or 0,
                confidence_scores=row['confidence_scores'] or [],
                algorithm_version=algorithm_version,
                processing_timestamp=datetime.now()
            )
            algorithm_results[row['protein_id']] = result

        return algorithm_results

    def _calculate_comparison_metrics(self, test_set: TestSet,
                                    baseline_results: Dict[int, PartitionResult],
                                    improved_results: Dict[int, PartitionResult],
                                    algorithm_version: str
    ) -> ComparisonMetrics:
        """Calculate detailed comparison metrics focused on boundary accuracy and fragment detection"""
        
        # Initialize counters for primary metrics
        domain_presence_correct = 0
        fragment_detection_correct = 0
        exact_boundary_matches = 0
        boundary_tolerance_5 = 0
        boundary_tolerance_10 = 0
        
        # Discontinuous domain tracking
        discontinuous_detected = 0
        discontinuous_total = 0
        discontinuous_boundaries_correct = 0
        
        # Fragment/peptide tracking
        true_peptides = 0
        predicted_peptides = 0
        correct_peptide_predictions = 0
        
        # Boundary error tracking
        boundary_errors = []
        confusion_matrix = defaultdict(int)
        improvement_cases = []
        regression_cases = []
        discontinuous_cases = []
        
        total_proteins = len(test_set.proteins)
        proteins_with_domains = 0
        boundary_comparisons = 0
        
        for protein in test_set.proteins:
            protein_id = protein.protein_id
            
            # Get results for this protein
            baseline = baseline_results.get(protein_id)
            improved = improved_results.get(protein_id)
            
            if not improved:
                continue  # Skip if no improved results
                
            # Domain presence accuracy
            manual_has_domain = protein.has_domain
            predicted_has_domain = improved.is_classified and not improved.is_peptide
            
            if manual_has_domain == predicted_has_domain:
                domain_presence_correct += 1
                
            # Fragment/peptide detection accuracy
            manual_is_fragment = protein.is_fragment
            predicted_is_peptide = improved.is_peptide
            
            if manual_is_fragment:
                true_peptides += 1
            if predicted_is_peptide:
                predicted_peptides += 1
            if manual_is_fragment == predicted_is_peptide:
                fragment_detection_correct += 1
                if manual_is_fragment:  # Both agree it's a fragment
                    correct_peptide_predictions += 1
                    
            # Confusion matrix
            if manual_has_domain and predicted_has_domain:
                confusion_matrix['domain_domain'] += 1
            elif manual_has_domain and predicted_is_peptide:
                confusion_matrix['domain_peptide'] += 1
            elif manual_is_fragment and predicted_has_domain:
                confusion_matrix['peptide_domain'] += 1
            elif manual_is_fragment and predicted_is_peptide:
                confusion_matrix['peptide_peptide'] += 1
            else:
                confusion_matrix['other'] += 1
                
            # Detailed boundary analysis for proteins with domains
            if manual_has_domain and predicted_has_domain:
                proteins_with_domains += 1
                
                # Analyze boundary agreement (this is simplified - in reality you'd need
                # to parse reference domain ranges from manual curation)
                if protein.boundaries_correct is True:
                    # For now, assume single domain case
                    if improved.domain_count == 1 and improved.domains:
                        domain = improved.domains[0]
                        
                        # Check if domain has discontinuous ranges (comma-separated)
                        domain_range = domain.get('range', '')
                        is_discontinuous = ',' in domain_range
                        
                        if is_discontinuous:
                            discontinuous_total += 1
                            discontinuous_cases.append({
                                'protein_id': protein_id,
                                'pdb_id': protein.pdb_id,
                                'chain_id': protein.chain_id,
                                'range': domain_range,
                                'boundaries_correct': protein.boundaries_correct
                            })
                            
                            if protein.boundaries_correct:
                                discontinuous_detected += 1
                                discontinuous_boundaries_correct += 1
                        
                        # For boundary tolerance analysis, we'd need reference ranges
                        # This is a placeholder for when you have explicit boundary data
                        boundary_comparisons += 1
                        if protein.boundaries_correct:
                            exact_boundary_matches += 1
                            boundary_tolerance_5 += 1
                            boundary_tolerance_10 += 1
                        
            # Compare baseline vs improved (focus on boundary and fragment improvements)
            if baseline and improved:
                baseline_score = self._calculate_boundary_fragment_score(baseline, protein)
                improved_score = self._calculate_boundary_fragment_score(improved, protein)
                
                if improved_score > baseline_score + 0.1:  # Significant improvement
                    improvement_cases.append({
                        'protein_id': protein_id,
                        'pdb_id': protein.pdb_id,
                        'chain_id': protein.chain_id,
                        'baseline_score': baseline_score,
                        'improved_score': improved_score,
                        'improvement': improved_score - baseline_score,
                        'issue_type': self._classify_improvement_type(baseline, improved, protein)
                    })
                elif improved_score < baseline_score - 0.1:  # Significant regression
                    regression_cases.append({
                        'protein_id': protein_id,
                        'pdb_id': protein.pdb_id,
                        'chain_id': protein.chain_id,
                        'baseline_score': baseline_score,
                        'improved_score': improved_score,
                        'regression': baseline_score - improved_score,
                        'issue_type': self._classify_improvement_type(baseline, improved, protein)
                    })
        
        # Calculate final metrics
        domain_presence_accuracy = domain_presence_correct / total_proteins if total_proteins > 0 else 0
        fragment_detection_accuracy = fragment_detection_correct / total_proteins if total_proteins > 0 else 0
        
        # Boundary metrics
        boundary_agreement_rate = exact_boundary_matches / boundary_comparisons if boundary_comparisons > 0 else 0
        exact_boundary_rate = exact_boundary_matches / boundary_comparisons if boundary_comparisons > 0 else 0
        tolerance_5_rate = boundary_tolerance_5 / boundary_comparisons if boundary_comparisons > 0 else 0
        tolerance_10_rate = boundary_tolerance_10 / boundary_comparisons if boundary_comparisons > 0 else 0
        
        # Discontinuous domain metrics
        discontinuous_detection_rate = discontinuous_detected / discontinuous_total if discontinuous_total > 0 else 0
        discontinuous_boundary_rate = discontinuous_boundaries_correct / discontinuous_total if discontinuous_total > 0 else 0
        
        # Fragment/peptide metrics  
        peptide_precision = correct_peptide_predictions / predicted_peptides if predicted_peptides > 0 else 0
        peptide_recall = correct_peptide_predictions / true_peptides if true_peptides > 0 else 0
        
        return ComparisonMetrics(
            test_set_id=test_set.test_set_id,
            algorithm_version=algorithm_version,
            
            # Primary metrics
            domain_presence_accuracy=domain_presence_accuracy,
            fragment_detection_accuracy=fragment_detection_accuracy,
            boundary_agreement_rate=boundary_agreement_rate,
            
            # Discontinuous domain metrics
            discontinuous_domain_detection=discontinuous_detection_rate,
            discontinuous_boundary_accuracy=discontinuous_boundary_rate,
            coverage_cutoff_effectiveness=0.0,  # Would calculate from coverage analysis
            
            # Boundary-specific metrics
            exact_boundary_matches=exact_boundary_rate,
            boundary_tolerance_5=tolerance_5_rate,
            boundary_tolerance_10=tolerance_10_rate,
            boundary_over_segmentation=0.0,  # Would calculate from domain count analysis
            boundary_under_segmentation=0.0,  # Would calculate from domain count analysis
            
            # Fragment/peptide metrics
            peptide_precision=peptide_precision,
            peptide_recall=peptide_recall,
            fragment_vs_domain_accuracy=fragment_detection_accuracy,
            
            # Secondary metrics
            domain_count_accuracy=0.0,  # Would calculate from domain count comparison
            classification_agreement_rate=0.0,  # Would calculate from classification comparison
            
            # Detailed data
            confusion_matrix=dict(confusion_matrix),
            boundary_error_distribution={'errors': boundary_errors},
            improvement_cases=improvement_cases,
            regression_cases=regression_cases,
            discontinuous_cases=discontinuous_cases
        )
    
    def _calculate_boundary_fragment_score(self, result: PartitionResult, manual: CurationDecision) -> float:
        """Calculate score focused on boundary accuracy and fragment detection"""
        
        score = 0.0
        
        # Fragment detection (50% weight - major focus)
        manual_is_fragment = manual.is_fragment
        predicted_is_peptide = result.is_peptide
        
        if manual_is_fragment == predicted_is_peptide:
            score += 0.5
            
        # Domain presence (30% weight)
        manual_has_domain = manual.has_domain
        predicted_has_domain = result.is_classified and not result.is_peptide
        
        if manual_has_domain == predicted_has_domain:
            score += 0.3
            
        # Boundary accuracy (20% weight - major focus)
        if manual.boundaries_correct is True and manual_has_domain and predicted_has_domain:
            score += 0.2
        elif manual.boundaries_correct is False and manual_has_domain and predicted_has_domain:
            # Penalty for incorrect boundaries
            score -= 0.1
            
        return max(0.0, score)
    
    def _classify_improvement_type(self, baseline: PartitionResult, 
                                 improved: PartitionResult, 
                                 manual: CurationDecision
    ) -> str:
        """Classify the type of improvement/regression"""
        
        baseline_is_peptide = baseline.is_peptide
        improved_is_peptide = improved.is_peptide
        manual_is_fragment = manual.is_fragment
        
        # Fragment detection improvements
        if manual_is_fragment:
            if not baseline_is_peptide and improved_is_peptide:
                return "fragment_detection_improvement"
            elif baseline_is_peptide and not improved_is_peptide:
                return "fragment_detection_regression"
                
        # Domain boundary improvements (simplified)
        if manual.has_domain:
            if baseline.domain_count != improved.domain_count:
                if manual.boundaries_correct:
                    return "boundary_segmentation_improvement"
                else:
                    return "boundary_segmentation_regression"
                    
        return "other"
    
    def _calculate_protein_score(self, result: PartitionResult, manual: CurationDecision) -> float:
        """Calculate a composite score for how well automated result matches manual curation"""
        
        score = 0.0
        
        # Domain presence agreement (0.4 weight)
        manual_has_domain = manual.has_domain
        predicted_has_domain = result.is_classified and not result.is_peptide
        
        if manual_has_domain == predicted_has_domain:
            score += 0.4
            
        # Fragment detection agreement (0.2 weight)
        if manual.is_fragment == result.is_peptide:
            score += 0.2
            
        # Domain count similarity (0.2 weight) - only if both have domains
        if manual_has_domain and predicted_has_domain:
            # Assume manual has 1 domain if not specified (could be improved)
            manual_domains = 1  # Could extract from reference_domain_id
            predicted_domains = result.domain_count
            
            if manual_domains == predicted_domains:
                score += 0.2
            elif abs(manual_domains - predicted_domains) == 1:
                score += 0.1  # Close but not exact
                
        # Confidence calibration (0.2 weight)
        if result.confidence_scores:
            avg_confidence = sum(result.confidence_scores) / len(result.confidence_scores)
            manual_confidence = manual.confidence_level / 5.0  # Normalize to 0-1
            
            # Reward well-calibrated confidence
            confidence_diff = abs(avg_confidence - manual_confidence)
            score += 0.2 * (1.0 - confidence_diff)
            
        return score
    
    def _save_comparison_metrics(self, test_set_id: int, metrics: ComparisonMetrics):
        """Save comparison metrics to database"""

        # Save individual metrics - FIXED: only use attributes that exist
        metric_items = [
            ('domain_presence_accuracy', metrics.domain_presence_accuracy),
            ('fragment_detection_accuracy', metrics.fragment_detection_accuracy),
            ('boundary_agreement_rate', metrics.boundary_agreement_rate),
            ('discontinuous_domain_detection', metrics.discontinuous_domain_detection),
            ('discontinuous_boundary_accuracy', metrics.discontinuous_boundary_accuracy),
            ('coverage_cutoff_effectiveness', metrics.coverage_cutoff_effectiveness),
            ('exact_boundary_matches', metrics.exact_boundary_matches),
            ('boundary_tolerance_5', metrics.boundary_tolerance_5),
            ('boundary_tolerance_10', metrics.boundary_tolerance_10),
            ('boundary_over_segmentation', metrics.boundary_over_segmentation),
            ('boundary_under_segmentation', metrics.boundary_under_segmentation),
            ('peptide_precision', metrics.peptide_precision),
            ('peptide_recall', metrics.peptide_recall),
            ('fragment_vs_domain_accuracy', metrics.fragment_vs_domain_accuracy),
            ('domain_count_accuracy', metrics.domain_count_accuracy),
            ('classification_agreement_rate', metrics.classification_agreement_rate)
        ]

        for metric_name, metric_value in metric_items:
            query = """
            INSERT INTO pdb_analysis.curation_test_metrics
            (test_set_id, algorithm_version, metric_name, metric_value)
            VALUES (%s, %s, %s, %s)
            """

            self.db.execute_query(query, (
                test_set_id, metrics.algorithm_version, metric_name, metric_value
            ))

        # Save complex metrics as JSON
        complex_metrics = [
            ('confusion_matrix', metrics.confusion_matrix),
            ('boundary_error_distribution', metrics.boundary_error_distribution),
            ('improvement_cases', metrics.improvement_cases),
            ('regression_cases', metrics.regression_cases),
            ('discontinuous_cases', metrics.discontinuous_cases)
        ]

        for metric_name, metric_data in complex_metrics:
            query = """
            INSERT INTO pdb_analysis.curation_test_metrics
            (test_set_id, algorithm_version, metric_name, metric_data)
            VALUES (%s, %s, %s, %s)
            """

            self.db.execute_query(query, (
                test_set_id, metrics.algorithm_version, metric_name, json.dumps(metric_data)
            ))

    def list_test_sets(self) -> List[Dict[str, Any]]:
        """List all available test sets"""
        
        query = """
        SELECT id, name, description, created_at, protein_count,
               curator_breakdown, decision_breakdown
        FROM pdb_analysis.curation_test_sets
        ORDER BY created_at DESC
        """
        
        return self.db.execute_dict_query(query)
    
    def promote_improved_results(self, test_set_id: int, algorithm_version: str,
                               dry_run: bool = True) -> Dict[str, int]:
        """Promote improved results to production partition_domains table"""
        
        self.logger.info(f"Promoting improved results for test set {test_set_id} (dry_run={dry_run})")
        
        # Get test proteins
        test_set = self.get_test_set(test_set_id)
        
        # Get improved results
        improved_results = self._get_algorithm_results(test_set_id, algorithm_version)
        
        if not improved_results:
            raise ValueError(f"No improved results found for algorithm version {algorithm_version}")
            
        # Track operations
        operations = {
            'proteins_updated': 0,
            'domains_deleted': 0,
            'domains_inserted': 0,
            'errors': 0
        }
        
        for protein in test_set.proteins:
            try:
                protein_id = protein.protein_id
                improved_result = improved_results.get(protein_id)
                
                if not improved_result:
                    continue
                    
                # Find corresponding partition_proteins record
                pp_query = """
                SELECT id FROM pdb_analysis.partition_proteins pp
                JOIN pdb_analysis.protein p ON pp.pdb_id = p.pdb_id AND pp.chain_id = p.chain_id
                WHERE p.id = %s
                ORDER BY pp."timestamp" DESC
                LIMIT 1
                """
                
                pp_result = self.db.execute_query(pp_query, (protein_id,))
                if not pp_result:
                    self.logger.warning(f"No partition_proteins record found for protein {protein_id}")
                    continue
                    
                partition_protein_id = pp_result[0][0]
                
                if not dry_run:
                    # Delete existing domains for this protein
                    delete_query = """
                    DELETE FROM pdb_analysis.partition_domains 
                    WHERE protein_id = %s
                    """
                    delete_result = self.db.execute_query(delete_query, (partition_protein_id,))
                    operations['domains_deleted'] += len(delete_result) if delete_result else 0
                    
                    # Insert new domains
                    for domain in improved_result.domains:
                        insert_query = """
                        INSERT INTO pdb_analysis.partition_domains
                        (protein_id, domain_number, start_pos, end_pos, range, source, 
                         source_id, confidence, t_group, h_group, x_group, a_group,
                         pdb_range, length)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        """
                        
                        self.db.execute_query(insert_query, (
                            partition_protein_id,
                            1,  # Would need proper domain numbering
                            domain.get('start_pos'),
                            domain.get('end_pos'),
                            domain.get('range'),
                            domain.get('source'),
                            domain.get('source_id'),
                            domain.get('confidence'),
                            domain.get('t_group'),
                            domain.get('h_group'),
                            domain.get('x_group'),
                            domain.get('a_group'),
                            domain.get('pdb_range'),
                            domain.get('length')
                        ))
                        operations['domains_inserted'] += 1
                    
                    # Update partition_proteins metadata
                    update_pp_query = """
                    UPDATE pdb_analysis.partition_proteins
                    SET is_classified = %s,
                        coverage = %s,
                        residues_assigned = %s,
                        domains_with_evidence = %s,
                        fully_classified_domains = %s,
                        "timestamp" = CURRENT_TIMESTAMP
                    WHERE id = %s
                    """
                    
                    self.db.execute_query(update_pp_query, (
                        improved_result.is_classified,
                        improved_result.coverage,
                        0,  # Would calculate from domains
                        len(improved_result.domains),
                        len(improved_result.domains),
                        partition_protein_id
                    ))
                    
                operations['proteins_updated'] += 1
                
            except Exception as e:
                self.logger.error(f"Error promoting results for protein {protein_id}: {e}")
                operations['errors'] += 1
                
        if dry_run:
            self.logger.info(f"DRY RUN: Would update {operations['proteins_updated']} proteins")
        else:
            self.logger.info(f"Updated {operations['proteins_updated']} proteins in production")
            
        return operations
    
    def validate_test_set_integrity(self, test_set_id: int) -> Dict[str, Any]:
        """Validate that test set proteins have all necessary data"""
        
        test_set = self.get_test_set(test_set_id)
        validation = {
            'total_proteins': len(test_set.proteins),
            'proteins_with_summaries': 0,
            'proteins_with_baseline': 0,
            'proteins_with_improved': 0,
            'missing_summaries': [],
            'missing_baseline': [],
            'missing_improved': []
        }
        
        for protein in test_set.proteins:
            # Check for domain summaries
            summaries = self._find_domain_summaries(protein.pdb_id, protein.chain_id)
            if summaries:
                validation['proteins_with_summaries'] += 1
            else:
                validation['missing_summaries'].append(f"{protein.pdb_id}_{protein.chain_id}")
                
            # Check for baseline results (would implement)
            # Check for improved results (would implement)
            
        return validation
        """Generate a comprehensive evaluation report"""
        
        test_set = self.get_test_set(test_set_id)
        
        report = f"""
CURATION TEST EVALUATION REPORT
===============================

Test Set: {test_set.name}
Description: {test_set.description}
Created: {test_set.created_at}
Proteins: {test_set.protein_count}

CURATOR BREAKDOWN:
{json.dumps(test_set.curator_breakdown, indent=2)}

DECISION BREAKDOWN:
{json.dumps(test_set.decision_breakdown, indent=2)}

ALGORITHM RESULTS:
        """
        
        # Get metrics for available algorithm versions
        metrics_query = """
        SELECT DISTINCT algorithm_version, metric_name, metric_value
        FROM pdb_analysis.curation_test_metrics
        WHERE test_set_id = %s
        ORDER BY algorithm_version, metric_name
        """
        
        metrics_data = self.db.execute_dict_query(metrics_query, (test_set_id,))
        
        if metrics_data:
            current_version = None
            for row in metrics_data:
                if row['algorithm_version'] != current_version:
                    current_version = row['algorithm_version']
                    report += f"\n\n{current_version.upper()}:\n"
                    report += "-" * (len(current_version) + 1) + "\n"
                    
                report += f"{row['metric_name']}: {row['metric_value']:.3f}\n"
        else:
            report += "\nNo algorithm results found.\n"
            
        return report

# Enhanced CurationTestManager with algorithm versioning
class EnhancedCurationTestManager(CurationTestManager):
    """Extended curation test manager with algorithm versioning support"""

    def __init__(self, context: ApplicationContext):
        super().__init__(context)

        # Add algorithm version manager
        self.version_manager = AlgorithmVersionManager(context.db)

        # Register your algorithm progression
        self._register_algorithm_lineage()

    def _register_algorithm_lineage(self):
        """Register the complete algorithm development lineage"""

        # Your development history
        algorithms = [
            create_baseline_algorithm(),
            create_coverage_focused_algorithm(),  # Your first iteration
            create_chain_blast_priority_algorithm(),  # Your next iteration
            create_chain_blast_only_algorithm()  # Experimental test
        ]

        for algorithm in algorithms:
            self.version_manager.register_algorithm(algorithm)
            self.logger.info(f"Registered {algorithm.version_id}: {algorithm.name}")

    def test_chain_blast_iteration(self, test_set_id: int) -> Dict[str, Any]:
        """
        Complete workflow to test chain BLAST iteration against your coverage iteration.

        This leverages your existing infrastructure while testing the new approach.
        """

        results = {
            'test_set_id': test_set_id,
            'algorithms_tested': [],
            'quick_test_results': {},
            'full_evaluation_results': {},
            'recommendations': []
        }

        # Step 1: Quick test chain BLAST priority algorithm
        self.logger.info("Step 1: Quick testing chain BLAST priority algorithm")

        quick_results = self.quick_test_algorithm_version(
            test_set_id=test_set_id,
            version_id="v2.0_chain_blast_priority",
            sample_size=15
        )

        results['quick_test_results']['chain_blast_priority'] = quick_results

        if quick_results['avg_improvement'] > 0.1:
            self.logger.info("✅ Chain BLAST priority shows improvements - proceeding to full test")

            # Step 2: Quick test chain BLAST only (experimental)
            self.logger.info("Step 2: Testing chain BLAST only (experimental)")

            quick_only_results = self.quick_test_algorithm_version(
                test_set_id=test_set_id,
                version_id="v2.1_chain_blast_only",
                sample_size=10
            )

            results['quick_test_results']['chain_blast_only'] = quick_only_results

            # Step 3: Full evaluation if quick tests are promising
            if quick_results['positive_improvements'] >= quick_results['successful_runs'] * 0.6:

                self.logger.info("Step 3: Running full evaluation")

                # Create baseline snapshot from coverage-focused algorithm
                snapshot_name = self.create_algorithm_baseline_snapshot(
                    test_set_id, "v1.1_coverage_focused"
                )

                # Run chain BLAST priority on full test set
                results_saved = self.rerun_algorithm_version(test_set_id, "v2.0_chain_blast_priority")

                # Evaluate improvements
                full_metrics = self.evaluate_algorithm_improvements(
                    test_set_id=test_set_id,
                    baseline_version="v1.1_coverage_focused",
                    improved_version="v2.0_chain_blast_priority"
                )

                results['full_evaluation_results'] = full_metrics
                results['algorithms_tested'] = ["v1.1_coverage_focused", "v2.0_chain_blast_priority"]

                # Generate recommendations
                results['recommendations'] = self._analyze_chain_blast_results(full_metrics)

            else:
                results['recommendations'] = [
                    "Chain BLAST priority shows mixed results in quick test",
                    "Review algorithm integration before full evaluation",
                    f"Positive improvements: {quick_results.get('positive_improvements', 0)}/{quick_results['successful_runs']}"
                ]

        else:
            results['recommendations'] = [
                "Chain BLAST priority algorithm may not be working as expected",
                "Review evidence weighting and integration",
                f"Average improvement: {quick_results.get('avg_improvement', 0):.3f}"
            ]

        return results

    def quick_test_algorithm_version(self, test_set_id: int, version_id: str,
                                   sample_size: int = 10) -> Dict[str, Any]:
        """Quick test a specific algorithm version"""

        algorithm = self.version_manager.get_algorithm(version_id)
        if not algorithm:
            raise ValueError(f"Algorithm version {version_id} not found")

        self.logger.info(f"Quick testing {version_id}: {algorithm.name}")

        # Get test set and sample
        test_set = self.get_test_set(test_set_id)
        import random
        sample_proteins = random.sample(test_set.proteins, min(sample_size, len(test_set.proteins)))

        # Create service with algorithm-specific configuration
        options, analyzer, processor = algorithm.create_service_components(self.context)

        service = DomainPartitionService(self.context)
        service.default_options = options
        service.analyzer = analyzer
        service.processor = processor

        # Track results with algorithm version context
        results = {
            'algorithm_version': version_id,
            'algorithm_name': algorithm.name,
            'total_tested': len(sample_proteins),
            'successful_runs': 0,
            'algorithm_errors': 0,
            'missing_summaries': 0,
            'comparison_results': [],
            'discontinuous_improvements': 0,  # NEW: Track discontinuous domain improvements
            'fragment_improvements': 0,       # Track fragment detection improvements
            'boundary_improvements': 0,        # Track boundary accuracy improvements
            'architectural_transfers': 0       # Track successful architectural transfers
        }

        for protein in sample_proteins:
            try:
                # Find domain summary
                summary_paths = self._find_domain_summaries(protein.pdb_id, protein.chain_id)

                if not summary_paths:
                    results['missing_summaries'] += 1
                    continue

                # Get baseline from coverage-focused algorithm (if available)
                baseline = self._get_current_partition_result(protein)

                # Run new algorithm
                improved_result = self._run_partition_for_protein_with_version(
                    service, protein, summary_paths, algorithm
                )

                if not improved_result:
                    results['algorithm_errors'] += 1
                    continue

                results['successful_runs'] += 1

                # Enhanced comparison focused on chain BLAST improvements
                comparison = self._compare_with_chain_blast_focus(baseline, improved_result, protein)
                results['comparison_results'].append(comparison)

                # Track specific improvement types
                if comparison.get('discontinuous_improvement'):
                    results['discontinuous_improvements'] += 1
                if comparison.get('fragment_improvement'):
                    results['fragment_improvements'] += 1
                if comparison.get('boundary_improvement'):
                    results['boundary_improvements'] += 1
                if comparison.get('architectural_transfer'):
                    results['architectural_transfers'] += 1

            except Exception as e:
                self.logger.error(f"Error testing {protein.pdb_id}_{protein.chain_id}: {e}")
                results['algorithm_errors'] += 1

        # Calculate chain BLAST specific metrics
        if results['comparison_results']:
            improvements = [r['overall_improvement'] for r in results['comparison_results']]
            results['avg_improvement'] = sum(improvements) / len(improvements)
            results['positive_improvements'] = len([i for i in improvements if i > 0])
            results['negative_improvements'] = len([i for i in improvements if i < 0])

            # Chain BLAST specific success rates
            results['discontinuous_improvement_rate'] = (
                results['discontinuous_improvements'] / results['successful_runs']
            )
            results['architectural_transfer_rate'] = (
                results['architectural_transfers'] / results['successful_runs']
            )

        return results

    def rerun_algorithm_version(self, test_set_id: int, version_id: str) -> int:
        """Re-run a specific algorithm version on test set"""

        algorithm = self.version_manager.get_algorithm(version_id)
        if not algorithm:
            raise ValueError(f"Algorithm version {version_id} not found")

        self.logger.info(f"Running {version_id} on test set {test_set_id}")

        # Get test set
        test_set = self.get_test_set(test_set_id)

        # Create service with algorithm configuration
        options, analyzer, processor = algorithm.create_service_components(self.context)

        service = DomainPartitionService(self.context)
        service.default_options = options
        service.analyzer = analyzer
        service.processor = processor

        # Track algorithm run in database
        run_id = self._start_algorithm_test_run(test_set_id, version_id, "full_evaluation")

        # Process proteins
        results_saved = 0

        for protein in test_set.proteins:
            try:
                summary_paths = self._find_domain_summaries(protein.pdb_id, protein.chain_id)

                if not summary_paths:
                    continue

                result = self._run_partition_for_protein_with_version(
                    service, protein, summary_paths, algorithm
                )

                if result:
                    # Save with algorithm version tracking
                    self._save_algorithm_result_versioned(
                        test_set_id, protein, result, version_id, run_id
                    )
                    results_saved += 1

            except Exception as e:
                self.logger.error(f"Error processing {protein.pdb_id}_{protein.chain_id}: {e}")

        # Complete algorithm run tracking
        self._complete_algorithm_test_run(run_id, results_saved, len(test_set.proteins))

        return results_saved

    def evaluate_algorithm_improvements(self, test_set_id: int,
                                      baseline_version: str,
                                      improved_version: str) -> ComparisonMetrics:
        """Evaluate algorithm improvements with version tracking"""

        baseline_algorithm = self.version_manager.get_algorithm(baseline_version)
        improved_algorithm = self.version_manager.get_algorithm(improved_version)

        if not baseline_algorithm or not improved_algorithm:
            raise ValueError("Algorithm versions not found")

        # Get test results
        test_set = self.get_test_set(test_set_id)
        baseline_results = self._get_algorithm_results(test_set_id, baseline_version)
        improved_results = self._get_algorithm_results(test_set_id, improved_version)

        # Calculate metrics with algorithm context
        metrics = self._calculate_comparison_metrics_with_context(
            test_set, baseline_results, improved_results,
            baseline_algorithm, improved_algorithm
        )

        # Save with version lineage
        self._save_comparison_metrics_versioned(test_set_id, metrics, baseline_version, improved_version)

        return metrics

    def _compare_with_chain_blast_focus(self, baseline: Optional[PartitionResult],
                                      improved: PartitionResult,
                                      manual: CurationDecision) -> Dict[str, Any]:
        """Enhanced comparison focused on chain BLAST improvements"""

        comparison = self._quick_compare_results(baseline, improved, manual)

        # Add chain BLAST specific analysis
        comparison.update({
            'discontinuous_improvement': False,
            'fragment_improvement': False,
            'boundary_improvement': False,
            'architectural_transfer': False
        })

        # Check for discontinuous domain improvements
        if improved.domains:
            for domain in improved.domains:
                domain_range = domain.get('range', '')
                if ',' in domain_range:  # Discontinuous domain detected
                    comparison['discontinuous_improvement'] = True
                    comparison['details']['discontinuous_domains'] = comparison['details'].get('discontinuous_domains', 0) + 1

        # Check for architectural transfer evidence
        if improved.domains:
            chain_blast_sources = sum(1 for d in improved.domains
                                     if d.get('source', '').startswith('chain'))
            if chain_blast_sources > 0:
                comparison['architectural_transfer'] = True
                comparison['details']['chain_blast_domains'] = chain_blast_sources

        # Enhanced fragment improvement detection
        if baseline and improved:
            baseline_peptide = baseline.is_peptide
            improved_peptide = improved.is_peptide
            manual_fragment = manual.is_fragment

            if manual_fragment and not baseline_peptide and improved_peptide:
                comparison['fragment_improvement'] = True
            elif not manual_fragment and baseline_peptide and not improved_peptide:
                comparison['fragment_improvement'] = True

        return comparison

    def _analyze_chain_blast_results(self, metrics: ComparisonMetrics) -> List[str]:
        """Generate recommendations based on chain BLAST algorithm results"""

        recommendations = []

        # Discontinuous domain analysis
        if metrics.discontinuous_domain_detection > 0.7:
            recommendations.append("✅ Excellent discontinuous domain detection improvement")
        elif metrics.discontinuous_domain_detection > 0.4:
            recommendations.append("⚠️ Moderate discontinuous domain improvement")
        else:
            recommendations.append("❌ Limited discontinuous domain improvement")

        # Fragment detection analysis
        if metrics.fragment_detection_accuracy > 0.8:
            recommendations.append("✅ Strong fragment detection (coverage algorithm success retained)")
        elif metrics.fragment_detection_accuracy > 0.6:
            recommendations.append("⚠️ Fragment detection acceptable but could be better")
        else:
            recommendations.append("❌ Fragment detection degraded - may need to retain coverage thresholds")

        # Boundary accuracy analysis
        if metrics.boundary_agreement_rate > 0.7:
            recommendations.append("✅ Good boundary accuracy via architectural transfer")
        else:
            recommendations.append("⚠️ Boundary accuracy needs improvement")

        # Overall assessment
        improvement_count = len([r for r in recommendations if r.startswith("✅")])

        if improvement_count >= 2:
            recommendations.append("🎯 RECOMMENDATION: Promote chain BLAST algorithm to production")
            recommendations.append("   - Shows improvements in target areas (discontinuous domains)")
            recommendations.append("   - Maintains fragment detection accuracy")
        elif improvement_count >= 1:
            recommendations.append("🔧 RECOMMENDATION: Refine chain BLAST algorithm")
            recommendations.append("   - Shows promise but needs optimization")
            recommendations.append("   - Consider hybrid approach with coverage algorithm")
        else:
            recommendations.append("🔄 RECOMMENDATION: Rethink chain BLAST approach")
            recommendations.append("   - May need different evidence weighting")
            recommendations.append("   - Consider alternative architectural transfer methods")

        return recommendations

    def create_algorithm_baseline_snapshot(self, test_set_id: int, algorithm_version: str) -> str:
        """Create baseline snapshot from specific algorithm version"""

        algorithm = self.version_manager.get_algorithm(algorithm_version)
        if not algorithm:
            raise ValueError(f"Algorithm version {algorithm_version} not found")

        snapshot_name = f"baseline_{algorithm_version}_{datetime.now().strftime('%Y%m%d_%H%M')}"

        self.logger.info(f"Creating baseline snapshot from {algorithm_version}: {snapshot_name}")

        # Use existing snapshot functionality but with version tracking
        baseline_snapshot = self.create_baseline_snapshot(test_set_id, snapshot_name)

        # Add version metadata
        self.db.execute_query("""
            INSERT INTO pdb_analysis.curation_test_metrics
            (test_set_id, algorithm_version, metric_name, metric_data)
            VALUES (%s, %s, %s, %s)
        """, (
            test_set_id,
            algorithm_version,
            "baseline_snapshot_metadata",
            json.dumps({
                "snapshot_name": snapshot_name,
                "source_algorithm": algorithm_version,
                "source_algorithm_name": algorithm.name,
                "created_for_comparison": True
            })
        ))

        return snapshot_name

    # Additional helper methods for algorithm version tracking
    def _start_algorithm_test_run(self, test_set_id: int, version_id: str, run_type: str) -> int:
        """Start tracking an algorithm test run"""

        query = """
        INSERT INTO pdb_analysis.algorithm_test_runs
        (version_id, test_set_id, run_type, test_parameters)
        VALUES (%s, %s, %s, %s)
        RETURNING id
        """

        result = self.db.execute_query(query, (
            version_id, test_set_id, run_type,
            json.dumps({"started_at": datetime.now().isoformat()})
        ))

        return result[0][0]

    def _complete_algorithm_test_run(self, run_id: int, successful: int, total: int):
        """Complete algorithm test run tracking"""

        query = """
        UPDATE pdb_analysis.algorithm_test_runs
        SET completed_at = CURRENT_TIMESTAMP,
            proteins_tested = %s,
            successful_runs = %s,
            status = 'completed'
        WHERE id = %s
        """

        self.db.execute_query(query, (total, successful, run_id))

# Usage example: Testing your chain BLAST iteration
def test_chain_blast_iteration_example():
    """
    Example workflow for testing your chain BLAST iteration
    using the enhanced curation test suite.
    """

    # Initialize enhanced test manager
    context = ApplicationContext("config/config.yml")
    manager = EnhancedCurationTestManager(context)

    # Test the chain BLAST iteration
    test_set_id = 1  # Your existing test set

    results = manager.test_chain_blast_iteration(test_set_id)

    print("CHAIN BLAST ITERATION TEST RESULTS")
    print("=" * 50)

    # Quick test results
    quick_results = results['quick_test_results'].get('chain_blast_priority', {})
    print(f"Quick Test Results:")
    print(f"  Discontinuous improvements: {quick_results.get('discontinuous_improvements', 0)}")
    print(f"  Architectural transfers: {quick_results.get('architectural_transfers', 0)}")
    print(f"  Fragment improvements: {quick_results.get('fragment_improvements', 0)}")
    print(f"  Average improvement: {quick_results.get('avg_improvement', 0):.3f}")

    # Full evaluation (if run)
    if 'full_evaluation_results' in results:
        full_metrics = results['full_evaluation_results']
        print(f"\nFull Evaluation Results:")
        print(f"  Discontinuous domain detection: {full_metrics.discontinuous_domain_detection:.3f}")
        print(f"  Fragment detection accuracy: {full_metrics.fragment_detection_accuracy:.3f}")
        print(f"  Boundary agreement rate: {full_metrics.boundary_agreement_rate:.3f}")

    # Recommendations
    print(f"\nRecommendations:")
    for rec in results['recommendations']:
        print(f"  {rec}")

    return results

def setup_logging(verbose: bool = False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )



def main():
    """Main CLI interface"""
    parser = argparse.ArgumentParser(description="Curation Test Suite for Domain Partition Evaluation")
    parser.add_argument('--config', type=str, default='config/config.yml', help='Configuration file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Extract command
    extract_parser = subparsers.add_parser('extract', help='Extract curated proteins as test set')
    extract_parser.add_argument('--min-confidence', type=int, default=3, 
                               help='Minimum curation confidence level (1-5)')
    extract_parser.add_argument('--include-fragments', action='store_true',
                               help='Include proteins marked as fragments')
    extract_parser.add_argument('--max-proteins', type=int,
                               help='Maximum number of proteins to extract')
    
    # Quick test subparser - for testing improvements before formal evaluation
    quicktest_parser = subparsers.add_parser('quicktest', help='Quick test improvements on sample proteins')
    quicktest_parser.add_argument('--test-set-id', type=int, required=True,
                                help='Test set ID to sample from')
    quicktest_parser.add_argument('--algorithm-version', type=str, default='improved_v1',
                                help='Algorithm version to test')
    quicktest_parser.add_argument('--sample-size', type=int, default=10,
                                help='Number of proteins to test (default: 10)')
    
    # Rerun command
    rerun_parser = subparsers.add_parser('rerun', help='Re-run algorithm on test set')
    rerun_parser.add_argument('--test-set-id', type=int, required=True,
                             help='Test set ID to process')
    rerun_parser.add_argument('--algorithm-version', type=str, default='improved_v1',
                             help='Algorithm version identifier')
    
    # Evaluate command
    eval_parser = subparsers.add_parser('evaluate', help='Evaluate algorithm improvements')
    eval_parser.add_argument('--test-set-id', type=int, required=True,
                            help='Test set ID to evaluate')
    eval_parser.add_argument('--baseline-version', type=str, default='original',
                            help='Baseline algorithm version')
    eval_parser.add_argument('--improved-version', type=str, default='improved_v1',
                            help='Improved algorithm version')
    
    # List command
    list_parser = subparsers.add_parser('list', help='List available test sets')
    
    # Create baseline snapshot subparser
    snapshot_parser = subparsers.add_parser('snapshot', help='Create baseline snapshot')
    snapshot_parser.add_argument('--test-set-id', type=int, required=True,
                                help='Test set ID for snapshot')
    snapshot_parser.add_argument('--snapshot-name', type=str,
                                help='Custom snapshot name')
    
    # Validate subparser
    validate_parser = subparsers.add_parser('validate', help='Validate test set integrity')
    validate_parser.add_argument('--test-set-id', type=int, required=True,
                               help='Test set ID to validate')
    
    # Promote subparser
    promote_parser = subparsers.add_parser('promote', help='Promote improved results to production')
    promote_parser.add_argument('--test-set-id', type=int, required=True,
                              help='Test set ID to promote')
    promote_parser.add_argument('--algorithm-version', type=str, default='improved_v1',
                              help='Algorithm version to promote')
    promote_parser.add_argument('--dry-run', action='store_true', default=True,
                              help='Perform dry run (default: True)')
    promote_parser.add_argument('--force', action='store_true',
                              help='Actually perform the promotion (overrides dry-run)')
    
    # Report command
    report_parser = subparsers.add_parser('report', help='Generate evaluation report')
    report_parser.add_argument('--test-set-id', type=int, required=True,
                              help='Test set ID for report')
    report_parser.add_argument('--include-details', action='store_true',
                              help='Include detailed improvement/regression cases')

    # Test algorithm version
    test_version_parser = subparsers.add_parser('test-version', help='Test specific algorithm version')
    test_version_parser.add_argument('--test-set-id', type=int, required=True)
    test_version_parser.add_argument('--version-id', type=str, required=True)
    test_version_parser.add_argument('--sample-size', type=int, default=10)

    # Test chain BLAST iteration
    test_chain_parser = subparsers.add_parser('test-chain-blast', help='Test chain BLAST iteration')
    test_chain_parser.add_argument('--test-set-id', type=int, required=True)

    # List algorithm versions
    list_versions_parser = subparsers.add_parser('list-versions', help='List algorithm versions')
    list_versions_parser.add_argument('--status', type=str, choices=['development', 'testing', 'production'])

    # Register new algorithm version
    register_parser = subparsers.add_parser('register-algorithm', help='Register new algorithm version')
    register_parser.add_argument('--version-file', type=str, required=True, help='JSON file with algorithm definition')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
        
    setup_logging(args.verbose)
    
    # Initialize context
    context = ApplicationContext(args.config)
    manager = CurationTestManager(context)
    
    try:
        if args.command == 'extract':
            test_set = manager.extract_curated_proteins(
                min_confidence=args.min_confidence,
                exclude_fragments=not args.include_fragments,
                max_proteins=args.max_proteins
            )
            print(f"Created test set {test_set.test_set_id} with {test_set.protein_count} proteins")
            print(f"Use --test-set-id {test_set.test_set_id} for subsequent commands")

        elif args.command == 'quicktest':
            results = manager.quick_test_improvements(
                args.test_set_id, args.algorithm_version, args.sample_size
            )

            print(f"QUICK TEST RESULTS for algorithm version {args.algorithm_version}")
            print("=" * 60)
            print(f"Total proteins tested: {results['total_tested']}")
            print(f"Successful runs: {results['successful_runs']}")
            print(f"Algorithm errors: {results['algorithm_errors']}")
            print(f"Missing summaries: {results['missing_summaries']}")

            if results['successful_runs'] > 0:
                print(f"\nIMPROVEMENT ANALYSIS:")
                print(f"  Average improvement: {results.get('avg_improvement', 0):.3f}")
                print(f"  Proteins improved: {results.get('positive_improvements', 0)}")
                print(f"  Proteins regressed: {results.get('negative_improvements', 0)}")
                print(f"  Fragment detection changes: {results['fragment_detection_changes']}")
                print(f"  Boundary changes: {results['boundary_changes']}")

                # Show some specific examples
                if results['comparison_results']:
                    print(f"\nTOP EXAMPLES:")
                    sorted_results = sorted(results['comparison_results'],
                                          key=lambda x: x['overall_improvement'], reverse=True)

                    print("  Best improvements:")
                    for result in sorted_results[:3]:
                        improvement = result['overall_improvement']
                        pdb_chain = f"{result['pdb_id']}_{result['chain_id']}"
                        print(f"    {pdb_chain}: +{improvement:.3f}")
                        if result['fragment_detection_changed']:
                            print(f"      -> Fragment detection changed")
                        if result['boundary_changed']:
                            print(f"      -> Boundaries changed")

                    print("  Potential regressions:")
                    for result in sorted_results[-2:]:
                        if result['overall_improvement'] < 0:
                            improvement = result['overall_improvement']
                            pdb_chain = f"{result['pdb_id']}_{result['chain_id']}"
                            print(f"    {pdb_chain}: {improvement:.3f}")

                # Decision guidance
                positive_rate = results.get('positive_improvements', 0) / results['successful_runs']
                avg_improvement = results.get('avg_improvement', 0)

                print(f"\nRECOMMENDATION:")
                if positive_rate > 0.6 and avg_improvement > 0.1:
                    print("  ✅ Algorithm shows good improvements!")
                    print("  ➡️  Proceed with formal evaluation:")
                    print(f"     python curation_test_suite.py snapshot --test-set-id {args.test_set_id}")
                    print(f"     python curation_test_suite.py rerun --test-set-id {args.test_set_id}")
                    print(f"     python curation_test_suite.py evaluate --test-set-id {args.test_set_id}")
                elif positive_rate > 0.4:
                    print("  ⚠️  Mixed results - some improvements, some regressions")
                    print("  ➡️  Consider reviewing regression cases before proceeding")
                else:
                    print("  ❌ Algorithm may not be working as expected")
                    print("  ➡️  Review algorithm integration and test cases")
            else:
                print("\n❌ No successful runs - check algorithm integration and data availability")

        elif args.command == 'rerun':
            results_saved = manager.rerun_algorithm(args.test_set_id, args.algorithm_version)
            print(f"Processed {results_saved} proteins with algorithm version {args.algorithm_version}")

        elif args.command == 'snapshot':
            snapshot_name = manager.create_baseline_snapshot(args.test_set_id, args.snapshot_name)
            print(f"Created baseline snapshot: {snapshot_name}")

        elif args.command == 'validate':
            validation = manager.validate_test_set_integrity(args.test_set_id)
            print(f"Test set validation for {args.test_set_id}:")
            print(f"  Total proteins: {validation['total_proteins']}")
            print(f"  With summaries: {validation['proteins_with_summaries']}")
            print(f"  With baseline: {validation['proteins_with_baseline']}")
            print(f"  With improved: {validation['proteins_with_improved']}")

            if validation['missing_summaries']:
                print(f"  Missing summaries: {len(validation['missing_summaries'])}")
                if len(validation['missing_summaries']) <= 5:
                    print(f"    {', '.join(validation['missing_summaries'])}")

        elif args.command == 'promote':
            dry_run = args.dry_run and not args.force
            operations = manager.promote_improved_results(
                args.test_set_id, args.algorithm_version, dry_run=dry_run
            )

            if dry_run:
                print(f"DRY RUN - Would promote {operations['proteins_updated']} proteins")
            else:
                print(f"Promoted {operations['proteins_updated']} proteins to production")

            print(f"  Domains deleted: {operations['domains_deleted']}")
            print(f"  Domains inserted: {operations['domains_inserted']}")
            if operations['errors'] > 0:
                print(f"  Errors: {operations['errors']}")

        elif args.command == 'evaluate':
            # FIXED: This section had syntax errors in the original
            baseline_results = manager._get_baseline_results(args.test_set_id)
            improved_results = manager._get_algorithm_results(args.test_set_id, args.improved_version)

            if not baseline_results:
                print("No baseline results found. Create a snapshot first with 'snapshot' command.")
                return 1

            if not improved_results:
                print(f"No improved results found for version {args.improved_version}. Run 'rerun' command first.")
                return 1

            test_set = manager.get_test_set(args.test_set_id)
            metrics = manager._calculate_comparison_metrics(
                test_set, baseline_results, improved_results, args.improved_version
            )
            manager._save_comparison_metrics(args.test_set_id, metrics)

            print(f"Evaluation complete for test set {args.test_set_id}")
            print("\nPRIMARY METRICS (Focus Areas):")
            print(f"  Domain presence accuracy: {metrics.domain_presence_accuracy:.3f}")
            print(f"  Fragment detection accuracy: {metrics.fragment_detection_accuracy:.3f}")
            print(f"  Boundary agreement rate: {metrics.boundary_agreement_rate:.3f}")

            print("\nDISCONTINUOUS DOMAIN METRICS (New Improvements):")
            print(f"  Discontinuous detection: {metrics.discontinuous_domain_detection:.3f}")
            print(f"  Discontinuous boundary accuracy: {metrics.discontinuous_boundary_accuracy:.3f}")

            print("\nBOUNDARY METRICS (Major Focus):")
            print(f"  Exact boundary matches: {metrics.exact_boundary_matches:.3f}")
            print(f"  Within 5 residues: {metrics.boundary_tolerance_5:.3f}")
            print(f"  Within 10 residues: {metrics.boundary_tolerance_10:.3f}")

            print("\nFRAGMENT/PEPTIDE METRICS (Major Focus):")
            print(f"  Peptide precision: {metrics.peptide_precision:.3f}")
            print(f"  Peptide recall: {metrics.peptide_recall:.3f}")

            print(f"\nIMPROVEMENT SUMMARY:")
            print(f"  Improvements: {len(metrics.improvement_cases)}")
            print(f"  Regressions: {len(metrics.regression_cases)}")
            print(f"  Discontinuous cases: {len(metrics.discontinuous_cases)}")

            # Show top improvements/regressions
            if metrics.improvement_cases:
                print(f"\nTop improvements:")
                for case in sorted(metrics.improvement_cases,
                                 key=lambda x: x['improvement'], reverse=True)[:3]:
                    print(f"  {case['pdb_id']}_{case['chain_id']}: +{case['improvement']:.3f} ({case['issue_type']})")

            if metrics.regression_cases:
                print(f"\nTop regressions:")
                for case in sorted(metrics.regression_cases,
                                 key=lambda x: x['regression'], reverse=True)[:3]:
                    print(f"  {case['pdb_id']}_{case['chain_id']}: -{case['regression']:.3f} ({case['issue_type']})")

        elif args.command == 'report':
            report = manager.generate_report(args.test_set_id)
            print(report)

            if hasattr(args, 'include_details') and args.include_details:
                # Add detailed cases
                metrics_query = """
                SELECT metric_data FROM pdb_analysis.curation_test_metrics
                WHERE test_set_id = %s AND metric_name IN ('improvement_cases', 'regression_cases')
                ORDER BY metric_name
                """

                details = manager.db.execute_dict_query(metrics_query, (args.test_set_id,))
                for detail in details:
                    if detail['metric_data']:
                        cases = json.loads(detail['metric_data'])
                        if cases:
                            print(f"\nDetailed cases: {len(cases)} found")
                            for case in cases[:5]:  # Show first 5
                                print(f"  {case}")

        elif args.command == 'list':
            test_sets = manager.list_test_sets()
            print("Available test sets:")
            print(f"{'ID':>3} {'Name':30} {'Proteins':>8} {'Created':20}")
            print("-" * 65)
            for ts in test_sets:
                print(f"{ts['id']:>3} {ts['name'][:30]:30} {ts['protein_count']:>8} {ts['created_at'].strftime('%Y-%m-%d %H:%M'):20}")

        elif args.command == 'test-version':
            manager = EnhancedCurationTestManager(context)
            results = manager.quick_test_algorithm_version(
                args.test_set_id, args.version_id, args.sample_size
            )

            print(f"Algorithm Version Test: {args.version_id}")
            print(f"  Success rate: {results['successful_runs']}/{results['total_tested']}")
            print(f"  Average improvement: {results.get('avg_improvement', 0):.3f}")
            print(f"  Discontinuous improvements: {results.get('discontinuous_improvements', 0)}")
            print(f"  Architectural transfers: {results.get('architectural_transfers', 0)}")

        elif args.command == 'test-chain-blast':
            manager = EnhancedCurationTestManager(context)
            results = manager.test_chain_blast_iteration(args.test_set_id)

            # Display comprehensive results
            # [Implementation as shown in example above]

        elif args.command == 'list-versions':
            manager = EnhancedCurationTestManager(context)
            algorithms = manager.version_manager.list_algorithms(args.status)

            print("Algorithm Versions:")
            print(f"{'Version ID':30} {'Name':40} {'Type':20} {'Status':15}")
            print("-" * 105)

            for alg in algorithms:
                print(f"{alg.version_id:30} {alg.name[:40]:40} {alg.iteration_type.value:20} {alg.deployment_status:15}")

        elif args.command == 'register-algorithm':
            # Load algorithm definition from JSON file and register
            with open(args.version_file) as f:
                alg_data = json.load(f)

            # Convert to AlgorithmVersion object and register
            # [Implementation details]

            return 0

    except Exception as e:
        logging.getLogger().error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())

# TYPICAL WORKFLOW EXAMPLE:
"""
# 1. Extract curated proteins as test set
python curation_test_suite.py extract --config config.yml --min-confidence 4
# Output: Created test set 1 with 87 proteins

# 2. QUICK TEST - Test algorithm on sample before committing to evaluation
python curation_test_suite.py quicktest --test-set-id 1 --sample-size 10
# Output: Quick test results showing if algorithm improvements are working
# ✅ If improvements look good, proceed to formal evaluation
# ❌ If not working, fix algorithm integration first

# 3. Create baseline snapshot of current results (after confirming improvements work)
python curation_test_suite.py snapshot --test-set-id 1
# Output: Created baseline snapshot: baseline_20250529_1430

# 4. Re-run improved algorithm on ALL test proteins
python curation_test_suite.py rerun --test-set-id 1 --algorithm-version improved_v1
# Output: Processed 85 proteins with algorithm version improved_v1

# 5. Evaluate improvements vs baseline
python curation_test_suite.py evaluate --test-set-id 1 --improved-version improved_v1
# Output: Evaluation results focusing on boundary and fragment improvements

# 6. Generate comprehensive report
python curation_test_suite.py report --test-set-id 1 --include-details
# Output: Detailed report with recommendations

# 7. If results look good, promote to production (dry run first)
python curation_test_suite.py promote --test-set-id 1 --algorithm-version improved_v1 --dry-run
# Output: DRY RUN - Would promote 85 proteins

# 8. Actually promote if satisfied
python curation_test_suite.py promote --test-set-id 1 --algorithm-version improved_v1 --force
# Output: Promoted 85 proteins to production

# Key focus areas for your improvements:
# - Fragment/peptide detection accuracy (sigmoid confidence behavior)
# - Domain boundary accuracy (may need structural features)  
# - Discontinuous domain handling (coverage cutoffs)
# - Reference domain coverage thresholds

# NEW: Quick test workflow for safe development
# python curation_test_suite.py extract --config config.yml --min-confidence 4
# python curation_test_suite.py quicktest --test-set-id 1 --sample-size 15
# -> Only proceed with snapshot/evaluation if quicktest shows improvements!
"""
