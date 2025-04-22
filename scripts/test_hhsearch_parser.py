# test_hhsearch_parser.py
import os, sys
import logging

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from ecod.pipelines.hhsearch import HHRParser, HHRToXMLConverter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("test")

# Test parsing
hhr_file = "./8ca5_t.develop291.hhr"
output_file = "./test_output.xml"
if os.path.exists(hhr_file):
    parser = HHRParser(logger)
    converter = HHRToXMLConverter(logger)
    
    # Parse HHR file
    hhr_data = parser.parse(hhr_file)
    if hhr_data:
        logger.info(f"Successfully parsed HHR file with {len(hhr_data['hits'])} hits")
        
        # Convert to XML
        xml_string = converter.convert(hhr_data, "test", "A", "develop")
        if xml_string:
            logger.info("Successfully converted HHR data to XML")
            
            # Save XML file
            if converter.save(xml_string, output_file):
                logger.info("Successfully saved XML file")
            else:
                logger.error("Failed to save XML file")
        else:
            logger.error("Failed to convert HHR data to XML")
    else:
        logger.error("Failed to parse HHR file")
else:
    logger.error(f"HHR file not found: {hhr_file}")