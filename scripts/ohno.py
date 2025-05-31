# COMPLETE FIX: Database access using correct method names

# ==========================================
# FIX 1: Update EvidenceAnalyzer.__init__ to get database from context
# ==========================================

def __init__(self, options: PartitionOptions, context=None):
    """FIXED: Get database manager from context"""
    
    self.options = options
    self.context = context
    self.logger = logging.getLogger(__name__)
    
    # Get database manager from context
    if context and hasattr(context, 'db_manager'):
        self.db_manager = context.db_manager
        self.logger.info("Database manager obtained from context.db_manager")
    elif context and hasattr(context, 'db'):
        self.db_manager = context.db
        self.logger.info("Database manager obtained from context.db")
    else:
        self.db_manager = None
        self.logger.warning("No database manager found in context")
    
    # ... existing initialization code ...
    
    # Initialize chain blast decomposition service
    decomp_config = DecompositionConfig(
        min_domain_size=options.min_domain_size,
        min_reference_coverage=getattr(options, 'min_reference_coverage', 0.70),
        log_all_attempts=True,
        log_failures_detail=True
    )
    
    self.decomposition_service = ChainBlastDecompositionService(
        config=decomp_config,
        context=context  # Pass context to decomposition service
    )
    
    # Debug database connectivity
    if self.db_manager:
        self.logger.info("Decomposition service initialized WITH database manager from context")
    else:
        self.logger.warning("Decomposition service initialized WITHOUT database manager")

# ==========================================
# FIX 2: Update DecompositionService.__init__ to get database from context
# ==========================================

def __init__(self, config: DecompositionConfig, context=None, perl_script_path: Optional[str] = None):
    """FIXED: Get database manager from context"""
    
    self.config = config
    self.context = context
    self.perl_script_path = perl_script_path
    self.logger = logging.getLogger(__name__)
    
    # Get database manager from context
    if context and hasattr(context, 'db_manager'):
        self.db_manager = context.db_manager
        self.logger.info("Database manager obtained from context.db_manager")
    elif context and hasattr(context, 'db'):
        self.db_manager = context.db
        self.logger.info("Database manager obtained from context.db")
    else:
        self.db_manager = None
        self.logger.warning("No database manager found in context")
    
    # ... existing stats initialization ...

# ==========================================
# FIX 3: Update _get_template_architecture to use correct method names
# ==========================================

def _get_template_architecture(self, pdb_id: str, chain_id: str) -> List[Dict[str, Any]]:
    """FIXED: Use correct DBManager method names"""
    
    if not self.db_manager:
        self.logger.warning("No database manager available for template architecture lookup")
        return []

    try:
        # FIXED: Query for actual database schema
        query = """
            SELECT d.domain_id, d.ecod_domain_id, d.range, 
                   d.t_group, d.h_group, d.x_group, d.a_group,
                   d.is_manual_rep, d.is_f70, d.is_f40, d.is_f99
            FROM pdb_analysis.domain d
            JOIN pdb_analysis.protein p ON d.protein_id = p.id
            WHERE p.pdb_id = %s AND p.chain_id = %s
            ORDER BY d.range
        """
        
        # FIXED: Use execute_dict_query instead of execute
        results = self.db_manager.execute_dict_query(query, (pdb_id, chain_id))
        
        domains = []
        for row in results:
            # FIXED: Parse range format like "A:225-384"
            range_str = row['range']
            start, end = self._parse_domain_range_with_chain(range_str)
            
            if start == 0 and end == 0:
                self.logger.warning(f"Could not parse range '{range_str}' for domain {row['domain_id']}")
                continue
            
            domain = {
                'domain_id': row['domain_id'],
                'ecod_domain_id': row['ecod_domain_id'],
                'start': start,
                'end': end,
                'range': range_str,
                't_group': row['t_group'],
                'h_group': row['h_group'],
                'x_group': row['x_group'],
                'a_group': row['a_group'],
                'is_manual_rep': row.get('is_manual_rep', False),
                'is_f70': row.get('is_f70', False),
                'is_f40': row.get('is_f40', False),
                'is_f99': row.get('is_f99', False)
            }
            domains.append(domain)
        
        self.logger.debug(f"Found {len(domains)} template domains for {pdb_id}_{chain_id}")
        return domains
        
    except Exception as e:
        self.logger.error(f"Error getting template architecture for {pdb_id}_{chain_id}: {e}")
        return []

def _parse_domain_range_with_chain(self, range_str: str) -> Tuple[int, int]:
    """FIXED: Parse range format like 'A:225-384' or 'd2:7-68'"""
    try:
        if ':' in range_str:
            # Format: "chain:start-end"
            chain_part, range_part = range_str.split(':', 1)
            
            if '-' in range_part:
                start, end = range_part.split('-', 1)
                return int(start), int(end)
            else:
                # Single position
                pos = int(range_part)
                return pos, pos
        else:
            # Format: "start-end" (no chain prefix)
            if '-' in range_str:
                start, end = range_str.split('-', 1)
                return int(start), int(end)
            else:
                pos = int(range_str)
                return pos, pos
                
    except ValueError as e:
        self.logger.warning(f"Could not parse range '{range_str}': {e}")
        return 0, 0

# ==========================================
# FIX 4: Update DomainPartitionService to pass context
# ==========================================

# In DomainPartitionService.__init__, change:
def __init__(self, context):
    """Initialize with application context"""
    self.context = context
    self.logger = logging.getLogger(__name__)
    
    # Create partition options
    partition_options = PartitionOptions(
        # ... your existing options ...
    )
    
    # FIXED: Pass context to analyzer
    self.analyzer = EvidenceAnalyzer(
        options=partition_options,
        context=context  # Pass context instead of db_manager
    )
    
    self.logger.info("DomainPartitionService initialized with context")

# ==========================================
# COMPLETE TEST: Verify all fixes work together
# ==========================================

def test_complete_database_fix():
    """Test that all database fixes work together"""
    
    print("=== COMPLETE DATABASE FIX TEST ===\n")
    
    from ecod.core.context import ApplicationContext
    from ecod.pipelines.domain_analysis.partition.service import DomainPartitionService
    
    # Test 1: Context database access
    print("1. Testing context database access:")
    context = ApplicationContext('config/config.yml')
    
    if hasattr(context, 'db_manager') and context.db_manager:
        print(f"   ✓ context.db_manager available: {type(context.db_manager)}")
        
        # Test correct method call
        try:
            results = context.db_manager.execute_dict_query("SELECT 1 as test")
            print(f"   ✓ execute_dict_query works: {results}")
        except Exception as e:
            print(f"   ✗ execute_dict_query failed: {e}")
            return False
    else:
        print(f"   ✗ No db_manager in context")
        return False
    
    # Test 2: Service creation with context
    print(f"\n2. Testing service creation:")
    try:
        service = DomainPartitionService(context)
        
        # Check analyzer database
        analyzer_has_db = hasattr(service.analyzer, 'db_manager') and service.analyzer.db_manager is not None
        print(f"   Analyzer has database: {analyzer_has_db}")
        
        # Check decomposition service database
        decomp_has_db = hasattr(service.analyzer.decomposition_service, 'db_manager') and service.analyzer.decomposition_service.db_manager is not None
        print(f"   Decomposition service has database: {decomp_has_db}")
        
        if not (analyzer_has_db and decomp_has_db):
            print(f"   ✗ Database not reaching services")
            return False
            
    except Exception as e:
        print(f"   ✗ Service creation failed: {e}")
        return False
    
    # Test 3: Template architecture lookup
    print(f"\n3. Testing template architecture lookup:")
    
    # Find test structures
    try:
        test_query = """
            SELECT p.pdb_id, p.chain_id, COUNT(d.id) as domain_count
            FROM pdb_analysis.protein p
            JOIN pdb_analysis.domain d ON p.id = d.protein_id
            WHERE p.pdb_id IN ('1gfl', '8ovp', '4kpt')
            GROUP BY p.pdb_id, p.chain_id
            HAVING COUNT(d.id) > 0
            LIMIT 3
        """
        test_structures = context.db_manager.execute_dict_query(test_query)
        
        if test_structures:
            print(f"   Found {len(test_structures)} test structures:")
            for struct in test_structures:
                print(f"     - {struct['pdb_id']}_{struct['chain_id']}: {struct['domain_count']} domains")
            
            # Test template lookup on first structure
            test_struct = test_structures[0]
            decomp_service = service.analyzer.decomposition_service
            
            template_domains = decomp_service._get_template_architecture(
                test_struct['pdb_id'], test_struct['chain_id']
            )
            
            print(f"   Template lookup test:")
            print(f"     Expected: {test_struct['domain_count']} domains")
            print(f"     Retrieved: {len(template_domains)} domains")
            
            if len(template_domains) > 0:
                print(f"   ✓ Template architecture lookup WORKING!")
                
                # Show sample domains
                for i, domain in enumerate(template_domains[:2]):
                    print(f"     {i+1}. {domain['domain_id']}: {domain['range']} ({domain['t_group']})")
                    
                return True
            else:
                print(f"   ✗ No domains retrieved")
                return False
        else:
            print(f"   ✗ No test structures found")
            return False
            
    except Exception as e:
        print(f"   ✗ Template lookup test failed: {e}")
        return False

if __name__ == "__main__":
    success = test_complete_database_fix()
    
    if success:
        print(f"\n🎉 ALL DATABASE FIXES WORKING!")
        print(f"\nNow test the full pipeline:")
        print(f"python full_precedence_test.py")
        print(f"\nExpected results:")
        print(f"- Successful decompositions > 0")
        print(f"- Architectural evidence generated")
        print(f"- Better GFP domain predictions")
    else:
        print(f"\n🔧 APPLY THE FIXES ABOVE")
        print(f"1. Update EvidenceAnalyzer.__init__ to get db from context")
        print(f"2. Update DecompositionService.__init__ to get db from context")
        print(f"3. Fix _get_template_architecture to use execute_dict_query")
        print(f"4. Update DomainPartitionService to pass context")
