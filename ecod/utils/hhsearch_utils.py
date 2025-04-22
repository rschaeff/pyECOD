	#!/usr/bin/env python3
"""
HHRParser - A robust parser for HHSearch result files (HHR format)
"""

import re
import logging
from typing import Dict, List, Any, Optional, Tuple


class HHRParser:
    """
    A comprehensive parser for HHSearch result files (HHR format).
    
    This parser handles the structured nature of HHR files with proper state
    transitions between the header section, hit summary table, and detailed
    alignment sections.
    """
    
    def __init__(self, logger=None):
        """Initialize parser with logger"""
        self.logger = logger or logging.getLogger("hhr_parser")
    
    def parse(self, hhr_file_path: str) -> Optional[Dict[str, Any]]:
        """
        Parse HHR file and extract structured data
        
        Args:
            hhr_file_path: Path to HHR file
            
        Returns:
            Dictionary with parsed data or None if parsing failed
        """
        try:
            with open(hhr_file_path, 'r') as f:
                content = f.read()
                
            # Parse file in distinct sections
            header = self._parse_header(content)
            hit_summary, alignment_start_line = self._parse_hit_summary_table(content)
            alignments = self._parse_alignments(content, alignment_start_line, hit_summary)
            
            # Combine hit summary with alignment details
            hits = self._merge_hit_data(hit_summary, alignments)
            
            return {
                'header': header,
                'hits': hits
            }
        except Exception as e:
            self.logger.error(f"Error parsing HHR file {hhr_file_path}: {str(e)}")
            self.logger.exception("Stack trace:")
            return None
    
    def _parse_header(self, content: str) -> Dict[str, Any]:
        """
        Parse header section of HHR file
        
        Args:
            content: Content of HHR file
            
        Returns:
            Dictionary with header information
        """
        header = {}
        lines = content.split('\n')
        
        # Find the end of header (the hit table header line)
        header_end = 0
        for i, line in enumerate(lines):
            if line.startswith(' No Hit'):
                header_end = i
                break
        
        # Process header lines
        for i in range(header_end):
            line = lines[i].strip()
            if not line:
                continue
                
            # Extract key-value pairs
            if ':' in line:
                key, value = line.split(':', 1)
                header[key.strip()] = value.strip()
            elif ' ' in line:
                # Handle space-separated key-value pairs
                parts = line.split(' ', 1)
                if len(parts) == 2:
                    key, value = parts
                    header[key.strip()] = value.strip()
        
        # Extract specific values using more specific logic
        for line in lines[:header_end]:
            if line.startswith('Query'):
                parts = line.split(None, 1)  # Split on first whitespace
                if len(parts) > 1:
                    header['query_id'] = parts[1].strip()
            elif line.startswith('Match_columns'):
                parts = line.split(None, 1)  # Split on first whitespace
                if len(parts) > 1:
                    try:
                        header['match_columns'] = int(parts[1].strip())
                    except ValueError:
                        header['match_columns'] = parts[1].strip()
            elif line.startswith('No_of_seqs'):
                parts = line.split(None, 1)  # Split on first whitespace
                if len(parts) > 1:
                    # Handle format like "150 out of 1487"
                    value_parts = parts[1].split('out of')
                    if len(value_parts) > 1:
                        try:
                            header['no_of_seqs'] = int(value_parts[0].strip())
                            header['total_seqs'] = int(value_parts[1].strip())
                        except ValueError:
                            header['no_of_seqs'] = parts[1].strip()
                    else:
                        try:
                            header['no_of_seqs'] = int(parts[1].strip())
                        except ValueError:
                            header['no_of_seqs'] = parts[1].strip()
            elif line.startswith('Neff'):
                parts = line.split(None, 1)  # Split on first whitespace
                if len(parts) > 1:
                    try:
                        header['neff'] = float(parts[1].strip())
                    except ValueError:
                        header['neff'] = parts[1].strip()
            elif line.startswith('Searched_HMMs'):
                parts = line.split(None, 1)  # Split on first whitespace
                if len(parts) > 1:
                    try:
                        header['searched_hmms'] = int(parts[1].strip())
                    except ValueError:
                        header['searched_hmms'] = parts[1].strip()
            elif line.startswith('Date'):
                parts = line.split(None, 1)  # Split on first whitespace
                if len(parts) > 1:
                    header['date'] = parts[1].strip()
            elif line.startswith('Command'):
                parts = line.split(None, 1)  # Split on first whitespace
                if len(parts) > 1:
                    header['command'] = parts[1].strip()
        
        return header
    
    def _parse_hit_summary_table(self, content: str) -> Tuple[List[Dict[str, Any]], int]:
        """
        Parse the hit summary table section of HHR file
        
        Args:
            content: Content of HHR file
            
        Returns:
            Tuple containing:
            - List of hit summary dictionaries
            - Line number where alignments start
        """
        lines = content.split('\n')
        hits = []
        
        # Find the hit table header
        table_start = None
        for i, line in enumerate(lines):
            if line.startswith(' No Hit'):
                table_start = i + 1
                break
        
        if table_start is None:
            self.logger.warning("Hit table header not found in HHR file")
            return [], len(lines)
        
        # Find the end of the table (first empty line after table_start)
        table_end = None
        for i in range(table_start, len(lines)):
            if not lines[i].strip():
                table_end = i
                break
        
        if table_end is None:
            table_end = len(lines)
        
        # Parse hit lines
        for i in range(table_start, table_end):
            line = lines[i].strip()
            if not line:
                continue
            
            # Parse hit line using regex to handle variable whitespace
            match = re.match(
                r'^\s*(\d+)\s+(\S+)\s+(\S+)\s+([0-9.-]+[Ee]?[+-]?\d*)\s+([0-9.-]+[Ee]?[+-]?\d*)\s+([0-9.-]+)\s+([0-9.-]+)\s+(\S+)?(.*)$',
                line
            )
            
            if match:
                hit_num = int(match.group(1))
                hit_id = match.group(2)
                description = match.group(9).strip() if match.group(9) else ""
                
                # Extract probability and other values
                try:
                    probability = float(match.group(3))
                except ValueError:
                    probability = match.group(3)
                
                try:
                    e_value = self._parse_scientific(match.group(4))
                except ValueError:
                    e_value = match.group(4)
                
                try:
                    p_value = self._parse_scientific(match.group(5))
                except ValueError:
                    p_value = match.group(5)
                
                try:
                    score = float(match.group(6))
                except ValueError:
                    score = match.group(6)
                
                try:
                    ss_score = float(match.group(7))
                except ValueError:
                    ss_score = match.group(7)
                
                # Extract alignment range info if present
                cols = match.group(8) if match.group(8) else ""
                
                # Process description to extract range information
                template_range = ""
                if description:
                    # Extract range info from description if available
                    range_match = re.search(r'([a-zA-Z]):([-\d]+)-([-\d]+)', description)
                    if range_match:
                        chain = range_match.group(1)
                        start = range_match.group(2)
                        end = range_match.group(3)
                        template_range = f"{chain}:{start}-{end}"
                
                # Create hit object
                hit = {
                    'hit_num': hit_num,
                    'hit_id': hit_id,
                    'probability': probability,
                    'e_value': e_value,
                    'p_value': p_value,
                    'score': score,
                    'ss_score': ss_score,
                    'cols': cols,
                    'description': description,
                    'template_range': template_range
                }
                
                hits.append(hit)
            else:
                self.logger.warning(f"Could not parse hit line: {line}")
        
        # Return hits and the line where alignments begin
        return hits, table_end + 1
    
    def _parse_alignments(self, content: str, alignment_start_line: int, hit_summary: List[Dict[str, Any]]) -> Dict[int, Dict[str, Any]]:
        """
        Parse the detailed alignment sections
        
        Args:
            content: Content of HHR file
            alignment_start_line: Line number where alignments start
            hit_summary: List of hit summaries to match alignments to
            
        Returns:
            Dictionary mapping hit numbers to alignment details
        """
        lines = content.split('\n')
        alignments = {}
        
        # Create lookup for hit IDs
        hit_id_to_num = {hit['hit_id']: hit['hit_num'] for hit in hit_summary}
        
        # State variables
        current_hit_num = None
        current_hit_id = None
        in_alignment = False
        alignment_blocks = []
        current_block = {
            'query_id': None, 'query_start': None, 'query_seq': "", 'query_end': None,
            'template_id': None, 'template_start': None, 'template_seq': "", 'template_end': None,
            'match_seq': "", 'pp_seq': ""  # pp = posterior probabilities
        }
        
        # Process alignment lines
        for i in range(alignment_start_line, len(lines)):
            line = lines[i].strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Check for hit header (starts with ">")
            if line.startswith('>'):
                # Store previous hit data if exists
                if current_hit_num is not None and in_alignment:
                    # Store the current alignment block if it has data
                    if current_block['query_seq'] and current_block['template_seq']:
                        alignment_blocks.append(current_block)
                    
                    # Store the alignment blocks for this hit
                    alignments[current_hit_num] = {
                        'hit_id': current_hit_id,
                        'alignment_blocks': alignment_blocks
                    }
                
                # Extract hit details
                match = re.search(r'>(\S+)', line)
                if match:
                    current_hit_id = match.group(1)
                    
                    # Find hit number from summary
                    current_hit_num = hit_id_to_num.get(current_hit_id)
                    if current_hit_num is None:
                        self.logger.warning(f"Hit ID {current_hit_id} not found in summary table")
                    
                    # Reset alignment state for new hit
                    in_alignment = True
                    alignment_blocks = []
                    current_block = {
                        'query_id': None, 'query_start': None, 'query_seq': "", 'query_end': None,
                        'template_id': None, 'template_start': None, 'template_seq': "", 'template_end': None,
                        'match_seq': "", 'pp_seq': ""
                    }
            
            # Process alignment lines
            elif in_alignment:
                # Identify line type using regex patterns
                q_match = re.match(r'^Q\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)', line)
                t_match = re.match(r'^T\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)', line)
                consensus_match = re.match(r'^Q\s+Consensus\s+\d+\s+(\S+)\s+\d+', line)
                template_consensus_match = re.match(r'^T\s+Consensus\s+\d+\s+(\S+)\s+\d+', line)
                confidence_match = re.match(r'^Confidence\s+(\S+)', line)
                match_line_match = re.match(r'^\s+(\S+)', line)
                
                # Check if this is the start of a new alignment block
                if q_match and not current_block['query_seq']:
                    # Start of a new alignment block
                    current_block = {
                        'query_id': q_match.group(1),
                        'query_start': int(q_match.group(2)),
                        'query_seq': q_match.group(3),
                        'query_end': int(q_match.group(4)),
                        'template_id': None,
                        'template_start': None,
                        'template_seq': "",
                        'template_end': None,
                        'match_seq': "",
                        'pp_seq': ""
                    }
                elif q_match and current_block['query_seq'] and not current_block['template_seq']:
                    # Continuing the query sequence
                    current_block['query_seq'] += q_match.group(3)
                    current_block['query_end'] = int(q_match.group(4))
                elif q_match and current_block['query_seq'] and current_block['template_seq']:
                    # Save previous block and start a new one
                    alignment_blocks.append(current_block)
                    current_block = {
                        'query_id': q_match.group(1),
                        'query_start': int(q_match.group(2)),
                        'query_seq': q_match.group(3),
                        'query_end': int(q_match.group(4)),
                        'template_id': None,
                        'template_start': None,
                        'template_seq': "",
                        'template_end': None,
                        'match_seq': "",
                        'pp_seq': ""
                    }
                elif t_match:
                    # Template alignment line
                    if not current_block['template_id']:
                        current_block['template_id'] = t_match.group(1)
                        current_block['template_start'] = int(t_match.group(2))
                    current_block['template_seq'] += t_match.group(3)
                    current_block['template_end'] = int(t_match.group(4))
                elif consensus_match:
                    # Query consensus line (optional)
                    pass
                elif template_consensus_match:
                    # Template consensus line (optional)
                    pass
                elif confidence_match:
                    # Confidence line
                    current_block['pp_seq'] = confidence_match.group(1)
                elif match_line_match and current_block['query_seq'] and not current_block['pp_seq']:
                    # Match line (residue quality indicators)
                    current_block['match_seq'] = match_line_match.group(1)
        
        # Store the last hit
        if current_hit_num is not None and in_alignment:
            # Store the current alignment block if it has data
            if current_block['query_seq'] and current_block['template_seq']:
                alignment_blocks.append(current_block)
            
            # Store the alignment blocks for this hit
            alignments[current_hit_num] = {
                'hit_id': current_hit_id,
                'alignment_blocks': alignment_blocks
            }
        
        return alignments
    
    def _merge_hit_data(self, hit_summary: List[Dict[str, Any]], alignments: Dict[int, Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Merge hit summary data with detailed alignment information
        
        Args:
            hit_summary: List of hit summaries
            alignments: Dictionary of alignment details keyed by hit number
            
        Returns:
            List of complete hit data dictionaries
        """
        merged_hits = []
        
        for hit in hit_summary:
            hit_num = hit['hit_num']
            
            # Create a copy of the hit data
            merged_hit = hit.copy()
            
            # Add alignment data if available
            if hit_num in alignments:
                alignment_data = alignments[hit_num]
                merged_hit['alignment_blocks'] = alignment_data['alignment_blocks']
                
                # Calculate query and template ranges from alignments
                if alignment_data['alignment_blocks']:
                    query_ranges = []
                    template_ranges = []
                    
                    for block in alignment_data['alignment_blocks']:
                        # Calculate query range with gapped positions removed
                        query_seq = block['query_seq']
                        query_start = block['query_start']
                        
                        ranges = self._calculate_range_from_alignment(query_seq, query_start)
                        query_ranges.extend(ranges)
                        
                        # Calculate template range with gapped positions removed
                        template_seq = block['template_seq']
                        template_start = block['template_start']
                        
                        ranges = self._calculate_range_from_alignment(template_seq, template_start)
                        template_ranges.extend(ranges)
                    
                    # Format ranges as comma-separated list of start-end pairs
                    merged_hit['query_range'] = ','.join([f"{start}-{end}" for start, end in query_ranges])
                    merged_hit['template_range'] = ','.join([f"{start}-{end}" for start, end in template_ranges])
                    
                    # Calculate alignment statistics
                    if alignment_data['alignment_blocks']:
                        stats = self._calculate_alignment_statistics(alignment_data['alignment_blocks'])
                        merged_hit.update(stats)
            
            merged_hits.append(merged_hit)
        
        return merged_hits
    
    def _calculate_range_from_alignment(self, seq: str, start_pos: int) -> List[Tuple[int, int]]:
        """
        Calculate sequence ranges from an alignment string, handling gaps
        
        Args:
            seq: Alignment string (with gaps)
            start_pos: Starting position
            
        Returns:
            List of (start, end) tuples representing ranges
        """
        ranges = []
        current_range_start = None
        current_pos = start_pos
        
        for char in seq:
            if char != '-':  # Not a gap
                if current_range_start is None:
                    current_range_start = current_pos
                current_pos += 1
            else:  # Gap
                if current_range_start is not None:
                    ranges.append((current_range_start, current_pos - 1))
                    current_range_start = None
        
        # Add the last range if exists
        if current_range_start is not None:
            ranges.append((current_range_start, current_pos - 1))
        
        return ranges
    
    def _calculate_alignment_statistics(self, alignment_blocks: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Calculate statistics from alignment blocks
        
        Args:
            alignment_blocks: List of alignment block dictionaries
            
        Returns:
            Dictionary of calculated statistics
        """
        stats = {
            'aligned_cols': 0,
            'identical_positions': 0,
            'similar_positions': 0,
            'identity_percentage': 0.0,
            'similarity_percentage': 0.0
        }
        
        total_aligned_cols = 0
        total_identical = 0
        total_similar = 0
        
        for block in alignment_blocks:
            query_seq = block['query_seq']
            template_seq = block['template_seq']
            match_seq = block.get('match_seq', '')
            
            # Count aligned columns (non-gap positions in both sequences)
            aligned_cols = sum(1 for q, t in zip(query_seq, template_seq) if q != '-' and t != '-')
            total_aligned_cols += aligned_cols
            
            # Count identical residues
            identical = sum(1 for q, t in zip(query_seq, template_seq) if q != '-' and t != '-' and q == t)
            total_identical += identical
            
            # Count similar residues (if match_seq has '+' symbols)
            if match_seq:
                similar = sum(1 for m, q, t in zip(match_seq, query_seq, template_seq) 
                             if q != '-' and t != '-' and m == '+')
                total_similar += similar
        
        stats['aligned_cols'] = total_aligned_cols
        stats['identical_positions'] = total_identical
        
        # Calculate percentages if there are aligned columns
        if total_aligned_cols > 0:
            stats['identity_percentage'] = (total_identical / total_aligned_cols) * 100.0
            stats['similarity_percentage'] = ((total_identical + total_similar) / total_aligned_cols) * 100.0
        
        return stats
    
    def _parse_scientific(self, value: str) -> float:
        """
        Parse scientific notation like 5.3E-16
        
        Args:
            value: String representation of a number
            
        Returns:
            Float value
        """
        # Handle different formats of scientific notation
        if 'E' in value or 'e' in value:
            parts = re.split('[Ee]', value)
            if len(parts) == 2:
                mantissa = float(parts[0])
                exponent = int(parts[1])
                return mantissa * (10 ** exponent)
        return float(value)


# Test code
if __name__ == "__main__":
    import sys
    
    # Configure logging
    logging.basicConfig(level=logging.INFO,
                      format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger("test")
    
    # Check if file path is provided
    if len(sys.argv) < 2:
        logger.error("Usage: python hhr_parser.py <hhr_file_path>")
        sys.exit(1)
    
    hhr_file = sys.argv[1]
    parser = HHRParser(logger)
    
    # Parse HHR file
    result = parser.parse(hhr_file)
    
    if result:
        # Print some basic info
        logger.info(f"Successfully parsed HHR file")
        logger.info(f"Query: {result['header'].get('query_id', 'Unknown')}")
        logger.info(f"Number of hits: {len(result['hits'])}")
        
        # Print first few hits
        for hit in result['hits'][:5]:
            logger.info(f"Hit {hit['hit_num']}: {hit['hit_id']} - Prob: {hit['probability']}, E-value: {hit['e_value']}")
            if 'alignment_blocks' in hit:
                for i, block in enumerate(hit['alignment_blocks']):
                    logger.info(f"  Block {i+1}: Query {block['query_start']}-{block['query_end']}, "
                               f"Template {block['template_start']}-{block['template_end']}")
    else:
        logger.error("Failed to parse HHR file")