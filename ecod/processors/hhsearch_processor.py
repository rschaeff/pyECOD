# ecod/processors/hhsearch_processor.py
class HHSearchProcessor:
    """Processor for HHSearch files"""
    
    def __init__(self, context=None):
        self.context = context
        self.logger = logging.getLogger("ecod.processors.hhsearch")
        self.parser = HHRParser(self.logger)
        self.converter = HHRToXMLConverter(self.logger)
    
    def process_hhr_file(self, hhr_path, output_path, metadata):
        """Process HHR file to XML"""
        # Implementation...
    
    def extract_hits(self, xml_path):
        """Extract hits from HHSearch XML"""
        # Implementation...
