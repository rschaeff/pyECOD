## processor.py
class HHSearchProcessor:
    """Process HHSearch results and integrate with BLAST evidence"""
    
    def __init__(self, context):
        """Initialize with application context"""
        self.context = context
        self.config = context.config_manager.config
        self.db = context.db
        self.logger = logging.getLogger("ecod.hhsearch_processor")
        
        self.parser = HHRParser(self.logger)
        self.converter = HHRToXMLConverter(self.logger)
        self.collator = DomainEvidenceCollator(self.logger)
    
    def _get_file_paths(self, batch_info, pdb_id, chain_id, ref_version):
        """Get standardized file paths for a protein chain
        
        Args:
            batch_info: Batch information dictionary
            pdb_id: PDB ID
            chain_id: Chain ID
            ref_version: Reference version
            
        Returns:
            Dictionary of file paths
        """
        base_path = batch_info['base_path']
        pdb_chain = f"{pdb_id}_{chain_id}"
        
        # Define standardized directories
        hhsearch_dir = os.path.join(base_path, "hhsearch")
        domains_dir = os.path.join(base_path, "domains")
        
        # Create directories if they don't exist
        os.makedirs(hhsearch_dir, exist_ok=True)
        os.makedirs(domains_dir, exist_ok=True)
        
        # Define file paths
        paths = {
            # HHSearch files
            'a3m': os.path.join(hhsearch_dir, f"{pdb_chain}.a3m"),
            'hhm': os.path.join(hhsearch_dir, f"{pdb_chain}.hhm"),
            'hhr': os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhr"),
            'hh_xml': os.path.join(hhsearch_dir, f"{pdb_chain}.{ref_version}.hhsearch.xml"),
            
            # BLAST result files
            'chain_blast': os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.chainwise_blast.xml"),
            'domain_blast': os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.blast.xml"),
            
            # Output domain summary file
            'domain_summary': os.path.join(domains_dir, f"{pdb_chain}.{ref_version}.domains.xml")
        }
        
        return paths
    
    def process_batch(self, batch_id):
        """Process HHSearch results for a batch
        
        Args:
            batch_id: Batch ID to process
            
        Returns:
            Number of successfully processed chains
        """
        # Get batch information
        batch_info = self._get_batch_info(batch_id)
        if not batch_info:
            self.logger.error(f"Batch {batch_id} not found")
            return 0
        
        # Get chains with completed HHSearch results
        chains = self._get_chains_with_hhsearch(batch_id)
        if not chains:
            self.logger.warning(f"No chains with completed HHSearch results found in batch {batch_id}")
            return 0
        
        self.logger.info(f"Found {len(chains)} chains with HHSearch results to process")
        
        processed_count = 0
        for chain in chains:
            success = self._process_chain(
                chain['pdb_id'], 
                chain['chain_id'], 
                chain['process_id'],
                batch_info,
                batch_info['ref_version']
            )
            
            if success:
                processed_count += 1
                
        self.logger.info(f"Successfully processed {processed_count} out of {len(chains)} chains")
        return processed_count
    
    def _process_chain(self, pdb_id, chain_id, process_id, batch_info, ref_version):
        """Process HHSearch results for a chain
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            process_id: Process ID
            batch_info: Batch information dictionary
            ref_version: Reference version
            
        Returns:
            True if successful
        """
        try:
            # Get file paths
            paths = self._get_file_paths(batch_info, pdb_id, chain_id, ref_version)
            
            # Ensure HHR file exists
            if not os.path.exists(paths['hhr']):
                self.logger.error(f"HHR file not found: {paths['hhr']}")
                return False
            
            # Parse HHR file
            hhr_data = self.parser.parse(paths['hhr'])
            if not hhr_data:
                self.logger.error(f"Failed to parse HHR file: {paths['hhr']}")
                return False
            
            # Convert to XML
            xml_string = self.converter.convert(hhr_data, pdb_id, chain_id, ref_version)
            if not xml_string:
                self.logger.error(f"Failed to convert HHR data to XML for {pdb_id}_{chain_id}")
                return False
            
            # Save HHSearch XML
            if not self.converter.save(xml_string, paths['hh_xml']):
                self.logger.error(f"Failed to save HHSearch XML: {paths['hh_xml']}")
                return False
            
            # Register HHSearch XML in database
            self._register_file(process_id, "hhsearch_xml", paths['hh_xml'], os.path.getsize(paths['hh_xml']))
            
            # Collate evidence and create domain summary
            summary_path = self._collate_evidence(pdb_id, chain_id, paths, ref_version)
            if not summary_path:
                self.logger.error(f"Failed to collate evidence for {pdb_id}_{chain_id}")
                return False
            
            # Register domain summary in database
            self._register_file(process_id, "domain_summary", summary_path, os.path.getsize(summary_path))
            
            # Update process status
            self._update_process_status(process_id, "domain_summary_complete")
            
            self.logger.info(f"Successfully processed HHSearch results for {pdb_id}_{chain_id}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error processing chain {pdb_id}_{chain_id}: {str(e)}")
            self._update_process_status(process_id, "error", str(e))
            return False
    
    def _collate_evidence(self, pdb_id, chain_id, paths, ref_version):
        """Collate domain evidence from different sources
        
        Args:
            pdb_id: PDB ID
            chain_id: Chain ID
            paths: Dictionary of file paths
            ref_version: Reference version
            
        Returns:
            Path to domain summary file or None if failed
        """
        try:
            # Check if required files exist
            for file_type in ['chain_blast', 'domain_blast', 'hh_xml']:
                if not os.path.exists(paths[file_type]):
                    self.logger.warning(f"Evidence file not found: {paths[file_type]}")
                    return None
            
            # Parse evidence files
            chain_blast_data = self._parse_xml(paths['chain_blast'])
            domain_blast_data = self._parse_xml(paths['domain_blast'])
            hhsearch_data = self._parse_xml(paths['hh_xml'])
            
            if not chain_blast_data or not domain_blast_data or not hhsearch_data:
                self.logger.error(f"Failed to parse one or more evidence files for {pdb_id}_{chain_id}")
                return None
            
            # Create collated XML
            root = ET.Element("domain_summ_doc")
            
            # Add metadata
            metadata = ET.SubElement(root, "metadata")
            ET.SubElement(metadata, "pdb_id").text = pdb_id
            ET.SubElement(metadata, "chain_id").text = chain_id
            ET.SubElement(metadata, "reference").text = ref_version
            ET.SubElement(metadata, "creation_date").text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            # Add chain blast evidence
            chain_blast_elem = ET.SubElement(root, "chain_blast_evidence")
            self._append_child_elements(chain_blast_elem, chain_blast_data.find("blast_summ_doc"))
            
            # Add domain blast evidence
            domain_blast_elem = ET.SubElement(root, "domain_blast_evidence")
            self._append_child_elements(domain_blast_elem, domain_blast_data.find("blast_summ_doc"))
            
            # Add HHSearch evidence
            hhsearch_elem = ET.SubElement(root, "hhsearch_evidence")
            self._append_child_elements(hhsearch_elem, hhsearch_data.find("hh_summ_doc"))
            
            # Add domain suggestions based on evidence
            domains_elem = ET.SubElement(root, "domain_suggestions")
            self._generate_domain_suggestions(domains_elem, chain_blast_data, domain_blast_data, hhsearch_data)
            
            # Write to file
            rough_string = ET.tostring(root, 'utf-8')
            reparsed = minidom.parseString(rough_string)
            pretty_xml = reparsed.toprettyxml(indent="  ")
            
            with open(paths['domain_summary'], 'w', encoding='utf-8') as f:
                f.write(pretty_xml)
                
            return paths['domain_summary']
            
        except Exception as e:
            self.logger.error(f"Error collating evidence: {str(e)}")
            return None