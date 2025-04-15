def create_summary(self, pdb_id: str, chain_id: str, reference: str, 
                 job_dump_dir: str, blast_only: bool = False) -> str:
    """Create domain summary for a protein chain"""
    # Define paths and check for existing files
    pdb_chain = f"{pdb_id}_{chain_id}"
    
    # Define output directory and filename once
    domains_dir = os.path.join(job_dump_dir, "domains")
    os.makedirs(domains_dir, exist_ok=True)
    
    # Define the output filename with proper context information
    suffix = ".blast_only" if blast_only else ""
    output_filename = f"{pdb_chain}.{reference}.blast_summ{suffix}.xml"
    output_path = os.path.join(domains_dir, output_filename)
    
    # Check for existing file
    if os.path.exists(output_path) and not self.config.get('force_overwrite', False):
        self.logger.warning(f"Output file {output_path} already exists, skipping...")
        return output_path
    
    # Rest of the method including reading sequence and checking if it's a peptide
    
    # When creating a special summary for peptides, use the already defined output_path:
    if sequence and len(sequence) < 30:
        # Create peptide summary XML
        peptide_summary = ET.Element("blast_summ_doc")
        # ...
        
        # Write to the previously defined output_path
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        tree = ET.ElementTree(peptide_summary)
        tree.write(output_path, encoding='utf-8', xml_declaration=True)
        
        self.logger.info(f"Created peptide summary: {output_path}")
        return output_path
    
    # Similarly, for regular summaries at the end of the method:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    tree = ET.ElementTree(root)
    tree.write(output_path, encoding='utf-8', xml_declaration=True)
    
    self.logger.info(f"Created domain summary: {output_path}")
    return output_path