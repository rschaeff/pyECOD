# test_hhsearch_parser.py
import os
import logging
from ecod.pipelines.hhsearch import HHRParser, HHRToXMLConverter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("test")

# Test parsing
hhr_file = "path/to/test.hhr"
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
            if converter.save(xml_string, "test_output.xml"):
                logger.info("Successfully saved XML file")
            else:
                logger.error("Failed to save XML file")
        else:
            logger.error("Failed to convert HHR data to XML")
    else:
        logger.error("Failed to parse HHR file")
else:
    logger.error(f"HHR file not found: {hhr_file}")