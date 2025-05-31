#!/usr/bin/env python3
"""
Example: Integrating Algorithm Versioning with Domain Analysis Pipeline
"""

from ecod.core.context import ApplicationContext
from ecod.evaluation.algorithm_versions.manager import AlgorithmVersionManager

def run_domain_partition_with_versioning(batch_id: int, algorithm_version_id: str = None):
    """
    Run domain partition with algorithm version tracking
    """
    
    # Initialize context and algorithm manager
    context = ApplicationContext("config/config.yml")
    version_manager = AlgorithmVersionManager(context)
    
    # Get algorithm version to use
    if algorithm_version_id:
        algorithm = version_manager.get_version(algorithm_version_id)
        if not algorithm:
            raise ValueError(f"Algorithm version {algorithm_version_id} not found")
    else:
        # Use current production algorithm
        algorithm = version_manager.get_production_version()
        if not algorithm:
            raise ValueError("No production algorithm found")
    
    print(f"Using algorithm: {algorithm.version_id} - {algorithm.name}")
    
    # Start tracking this algorithm run
    run_id = version_manager.start_algorithm_run(
        algorithm.version_id, 
        batch_id=batch_id,
        run_type="production"
    )
    
    try:
        # Get algorithm configuration for domain analysis
        config_dict = algorithm.to_config_dict()
        domain_config = config_dict.get('domain_analysis', {})
        
        print(f"Algorithm configuration:")
        print(f"  - Evidence weights: {domain_config.get('evidence_weights', {})}")
        print(f"  - Coverage thresholds: {domain_config.get('coverage_thresholds', {})}")
        print(f"  - Behavioral flags: {domain_config.get('behavioral_flags', {})}")
        
        # TODO: Use this configuration in your DomainPartitionService
        # partition_service = DomainPartitionService(context, algorithm_config=domain_config)
        # results = partition_service.process_batch(batch_id)
        
        # Simulate processing results
        results = {
            "proteins_processed": 150,
            "successful": 145,
            "failed": 5,
            "success_rate": 0.967,
            "avg_domains_per_protein": 2.3
        }
        
        print(f"Processing results: {results}")
        
        # Complete the algorithm run tracking
        version_manager.complete_algorithm_run(run_id, results)
        
        print(f"✅ Algorithm run {run_id} completed successfully")
        
        return results
        
    except Exception as e:
        # Mark run as failed
        version_manager.complete_algorithm_run(run_id, {"error": str(e), "status": "failed"})
        print(f"❌ Algorithm run {run_id} failed: {e}")
        raise

def compare_algorithms_on_test_set():
    """
    Example: Compare different algorithms on a test set
    """
    context = ApplicationContext("config/config.yml")
    version_manager = AlgorithmVersionManager(context)
    
    # Get algorithms to compare
    baseline = version_manager.get_version("v1.0_baseline")
    improved = version_manager.get_version("v2.0_improved_coverage")
    experimental = version_manager.get_version("v2.1_chain_blast_priority")
    
    algorithms_to_test = [baseline, improved, experimental]
    
    for algorithm in algorithms_to_test:
        if not algorithm:
            continue
            
        print(f"\n🧪 Testing algorithm: {algorithm.version_id}")
        
        # Start evaluation run
        run_id = version_manager.start_algorithm_run(
            algorithm.version_id,
            run_type="evaluation"
        )
        
        # Simulate running algorithm on test set
        # results = run_algorithm_on_test_set(algorithm)
        
        # Mock results for example
        results = {
            "test_proteins": 50,
            "boundary_accuracy": 0.85 + (0.1 if "improved" in algorithm.version_id else 0),
            "fragment_detection": 0.78 + (0.15 if "chain_blast" in algorithm.version_id else 0),
            "processing_time": 120.5
        }
        
        print(f"  Results: {results}")
        
        version_manager.complete_algorithm_run(run_id, results)

def manage_algorithm_lifecycle():
    """
    Example: Complete algorithm lifecycle management
    """
    context = ApplicationContext("config/config.yml")
    version_manager = AlgorithmVersionManager(context)
    
    print("📊 Current Algorithm Status:")
    print("-" * 50)
    
    algorithms = version_manager.list_versions()
    for algo in algorithms:
        print(f"{algo.version_id:30} {algo.status.value:12} {algo.name}")
    
    print(f"\n🏆 Production Algorithm:")
    production = version_manager.get_production_version()
    if production:
        print(f"  {production.version_id} - {production.name}")
        
        # Show configuration summary
        config = production.to_config_dict()
        evidence_weights = config.get('domain_analysis', {}).get('evidence_weights', {})
        print(f"  Evidence weights: {evidence_weights}")
    else:
        print("  No production algorithm found")

if __name__ == "__main__":
    print("🧬 Algorithm Versioning Integration Examples")
    print("=" * 60)
    
    try:
        # Example 1: Run domain partition with specific algorithm
        print("\n1️⃣ Running Domain Partition with Version Tracking")
        run_domain_partition_with_versioning(
            batch_id=999, 
            algorithm_version_id="v2.0_improved_coverage"
        )
        
        # Example 2: Algorithm comparison
        print("\n2️⃣ Comparing Algorithms on Test Set")
        compare_algorithms_on_test_set()
        
        # Example 3: Lifecycle management
        print("\n3️⃣ Algorithm Lifecycle Management")
        manage_algorithm_lifecycle()
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
