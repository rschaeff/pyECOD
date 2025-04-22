    #!/usr/bin/env python3
"""
HHRParser - A robust parser for HHSearch result files (HHR format)
"""

import os
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
            # Check if file exists
            if not os.path.exists(hhr_file_path):
                self.logger.error(f"File not found: {hhr_file_path}")
                return None
                
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
            
            # Handle the specific format of HHR hit lines
            try:
                # First extract the hit number and ID at the beginning
                hit_match = re.match(r'^\s*(\d+)\s+(\S+)', line)
                if hit_match:
                    hit_num = int(hit_match.group(1))
                    hit_id = hit_match.group(2)
                    
                    # Extract the description, which is the remainder after hit_id
                    description_parts = line.split(hit_id, 1)
                    description = description_parts[1].strip() if len(description_parts) > 1 else ""
                    
                    # Now extract probability, e-value, etc. which appear in a regular pattern
                    # with multiple whitespace characters between them
                    values_match = re.search(r'\s+([\d\.]+)\s+([\dE\.\+\-]+)\s+([\dE\.\+\-]+)\s+([\d\.]+)\s+([\d\.]+)\s+(\d+)\s+([\d\-]+)\s+([\d\-]+)', description)
                    
                    if values_match:
                        try:
                            probability = float(values_match.group(1))
                        except ValueError:
                            probability = values_match.group(1)
                        
                        try:
                            e_value = self._parse_scientific(values_match.group(2))
                        except ValueError:
                            e_value = values_match.group(2)
                        
                        try:
                            p_value = self._parse_scientific(values_match.group(3))
                        except ValueError:
                            p_value = values_match.group(3)
                        
                        try:
                            score = float(values_match.group(4))
                        except ValueError:
                            score = values_match.group(4)
                        
                        try:
                            ss_score = float(values_match.group(5))
                        except ValueError:
                            ss_score = values_match.group(5)
                        
                        cols = values_match.group(6)
                        query_range = values_match.group(7)
                        template_range = values_match.group(8)
                        
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
                            'query_range': query_range,
                            'template_range': template_range,
                            'description': description
                        }
                        
                        hits.append(hit)
                    else:
                        # Fallback method if the regex doesn't match
                        # Split by whitespace and try to extract values by position
                        parts = re.split(r'\s+', line.strip())
                        if len(parts) >= 9:
                            try:
                                probability = float(parts[2])
                            except ValueError:
                                probability = parts[2]
                            
                            try:
                                e_value = self._parse_scientific(parts[3])
                            except ValueError:
                                e_value = parts[3]
                            
                            try:
                                p_value = self._parse_scientific(parts[4])
                            except ValueError:
                                p_value = parts[4]
                            
                            try:
                                score = float(parts[5])
                            except ValueError:
                                score = parts[5]
                            
                            try:
                                ss_score = float(parts[6])
                            except ValueError:
                                ss_score = parts[6]
                            
                            cols = parts[7]
                            remaining = ' '.join(parts[8:])
                            
                            hit = {
                                'hit_num': hit_num,
                                'hit_id': hit_id,
                                'probability': probability,
                                'e_value': e_value,
                                'p_value': p_value,
                                'score': score,
                                'ss_score': ss_score,
                                'cols': cols,
                                'description': description
                            }
                            
                            hits.append(hit)
                        else:
                            self.logger.warning(f"Could not parse hit line: {line}")
                else:
                    self.logger.warning(f"Could not parse hit line: {line}")
            except Exception as e:
                self.logger.warning(f"Error parsing hit line: {line}, reason: {str(e)}")
        
        # Return hits and the line where alignments begin
        return hits, table_end + 1
    
	def _parse_alignments(self, content: str, alignment_start_line: int, hit_summary: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
	    """
	    Parse the detailed alignment sections with explicit distinction between sequence types
	    
	    Args:
	        content: Content of HHR file
	        alignment_start_line: Line number where alignments start
	        hit_summary: List of hit summaries to match alignments to
	            
	    Returns:
	        Dictionary mapping hit IDs to alignment details
	    """
	    lines = content.split('\n')
	    alignments = {}
	    
	    # State variables
	    current_hit_id = None
	    current_hit_data = {
	        'query_seq': "",
	        'template_seq': "",
	        'match_seq': "",
	        'query_consensus': "",
	        'template_consensus': "",
	        'stats': {}
	    }
	    in_alignment_block = False
	    
	    # Process alignment lines
	    i = alignment_start_line
	    while i < len(lines):
	        line = lines[i].strip()
	        
	        # Check for hit header (starts with ">")
	        if line.startswith('>'):
	            # Save previous hit data if it exists
	            if current_hit_id is not None and (current_hit_data['query_seq'] or current_hit_data['template_seq']):
	                alignments[current_hit_id] = {
	                    'alignment_blocks': [{
	                        'query_seq': current_hit_data['query_seq'],
	                        'template_seq': current_hit_data['template_seq'],
	                        'match_seq': current_hit_data['match_seq'],
	                        'consensus_q': current_hit_data['query_consensus'],
	                        'consensus_t': current_hit_data['template_consensus'],
	                        'stats': current_hit_data['stats']
	                    }]
	                }
	            
	            # Extract hit ID from the header
	            header_parts = line[1:].split(None, 1)
	            if header_parts:
	                current_hit_id = header_parts[0]
	                # Reset hit data for the new hit
	                current_hit_data = {
	                    'query_seq': "",
	                    'template_seq': "",
	                    'match_seq': "",
	                    'query_consensus': "",
	                    'template_consensus': "",
	                    'stats': {}
	                }
	                in_alignment_block = False
	            else:
	                self.logger.warning(f"Malformed alignment header: {line}")
	                current_hit_id = None
	        
	        # Check for alignment probability line (starts with "Probab=")
	        elif line.startswith('Probab=') and current_hit_id is not None:
	            # Parse stats from probability line
	            stats = {}
	            for stat_pair in line.split():
	                if '=' in stat_pair:
	                    key, value = stat_pair.split('=', 1)
	                    try:
	                        value = float(value.replace('%', '')) / 100.0 if '%' in value else float(value)
	                    except ValueError:
	                        pass
	                    stats[key] = value
	            
	            current_hit_data['stats'] = stats
	            in_alignment_block = True
	        
	        # Process alignment lines within a block
	        elif in_alignment_block and current_hit_id is not None:
	            # QUERY SEQUENCE LINE - starts with "Q" followed by a sequence ID, contains uppercase letters
	            # Pattern: Q followed by an identifier that is NOT "Consensus"
	            if line.startswith('Q ') and 'Consensus' not in line and re.search(r'[A-Z]', line):
	                parts = line.split()
	                if len(parts) >= 5:  # Q id start seq end
	                    query_id = parts[1]
	                    query_start = int(parts[2])
	                    query_seq = parts[3]
	                    query_end = int(parts[4])
	                    
	                    # Check if this is an actual sequence (typically uppercase) vs consensus (has ~)
	                    if '~' not in query_seq and query_id != "Consensus":
	                        current_hit_data['query_seq'] += query_seq
	            
	            # QUERY CONSENSUS LINE - explicitly contains "Consensus" and typically has tilde (~)
	            elif line.startswith('Q Consensus'):
	                parts = line.split()
	                if len(parts) >= 5:  # Q Consensus start seq end
	                    consensus_start = int(parts[2])
	                    consensus_seq = parts[3]
	                    consensus_end = int(parts[4])
	                    
	                    # Verify this is a consensus line (should contain ~)
	                    if '~' in consensus_seq:
	                        current_hit_data['query_consensus'] += consensus_seq
	            
	            # MATCH LINE - line of spaces and symbols between query and template
	            elif line.strip() and line[0] == ' ' and not line.startswith(' Q') and not line.startswith(' T'):
	                match_seq = line.strip()
	                current_hit_data['match_seq'] += match_seq
	            
	            # TEMPLATE SEQUENCE LINE - starts with "T" followed by template ID
	            elif line.startswith('T ') and 'Consensus' not in line and re.search(r'[A-Z]', line):
	                parts = line.split()
	                if len(parts) >= 5:  # T id start seq end
	                    template_id = parts[1]
	                    template_start = int(parts[2])
	                    template_seq = parts[3]
	                    template_end = int(parts[4])
	                    
	                    # Check if this is an actual sequence (typically uppercase) vs consensus (has ~)
	                    if '~' not in template_seq and template_id != "Consensus":
	                        current_hit_data['template_seq'] += template_seq
	            
	            # TEMPLATE CONSENSUS LINE - explicitly contains "Consensus" and typically has tilde (~)
	            elif line.startswith('T Consensus'):
	                parts = line.split()
	                if len(parts) >= 5:  # T Consensus start seq end
	                    consensus_start = int(parts[2])
	                    consensus_seq = parts[3]
	                    consensus_end = int(parts[4])
	                    
	                    # Verify this is a consensus line (should contain ~)
	                    if '~' in consensus_seq:
	                        current_hit_data['template_consensus'] += consensus_seq
	        
	        i += 1
	    
	    # Add the last hit if it exists
	    if current_hit_id is not None and (current_hit_data['query_seq'] or current_hit_data['template_seq']):
	        alignments[current_hit_id] = {
	            'alignment_blocks': [{
	                'query_seq': current_hit_data['query_seq'],
	                'template_seq': current_hit_data['template_seq'],
	                'match_seq': current_hit_data['match_seq'],
	                'consensus_q': current_hit_data['query_consensus'],
	                'consensus_t': current_hit_data['template_consensus'],
	                'stats': current_hit_data['stats']
	            }]
	        }
	    
	    return alignments  
    def _merge_hit_data(self, hit_summary: List[Dict[str, Any]], alignments: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Merge hit summary data with detailed alignment information
        
        Args:
            hit_summary: List of hit summaries
            alignments: Dictionary of alignment details keyed by hit ID
            
        Returns:
            List of complete hit data dictionaries
        """
        merged_hits = []
        
        for hit in hit_summary:
            hit_id = hit['hit_id']
            
            # Create a copy of the hit data
            merged_hit = hit.copy()
            
            # Add alignment data if available
            if hit_id in alignments:
                alignment_data = alignments[hit_id]
                merged_hit['alignment_blocks'] = alignment_data['alignment_blocks']
                
                # Calculate query and template ranges from alignments
                if alignment_data['alignment_blocks']:
                    # Extract statistics from first alignment block
                    if 'stats' in alignment_data['alignment_blocks'][0]:
                        merged_hit.update(alignment_data['alignment_blocks'][0]['stats'])
                    
                    # Calculate query and template ranges
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
                    if query_ranges:
                        merged_hit['query_range_calculated'] = ','.join([f"{start}-{end}" for start, end in query_ranges])
                    
                    if template_ranges:
                        merged_hit['template_range_calculated'] = ','.join([f"{start}-{end}" for start, end in template_ranges])
                    
                    # Calculate alignment statistics
                    stats = self._calculate_alignment_statistics(alignment_data['alignment_blocks'])
                    for key, value in stats.items():
                        if key not in merged_hit:
                            merged_hit[key] = value
            
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

    def _merge_hit_data(self, hit_summary: List[Dict[str, Any]], alignments: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Merge hit summary data with detailed alignment information
        
        Args:
            hit_summary: List of hit summaries
            alignments: Dictionary of alignment details keyed by hit ID
            
        Returns:
            List of complete hit data dictionaries
        """
        merged_hits = []
        
        for hit in hit_summary:
            hit_id = hit['hit_id']
            
            # Create a copy of the hit data
            merged_hit = hit.copy()
            
            # Add alignment data if available
            if hit_id in alignments:
                alignment_data = alignments[hit_id]
                merged_hit['alignment_blocks'] = alignment_data['alignment_blocks']
                
                # Extract statistics from first alignment block if available
                if alignment_data['alignment_blocks'] and 'stats' in alignment_data['alignment_blocks'][0]:
                    stats = alignment_data['alignment_blocks'][0]['stats']
                    for key, value in stats.items():
                        if key not in merged_hit:
                            merged_hit[key] = value
                
                # Prioritize range information - use a single source of truth
                # First check if we already have valid ranges from hit summary
                has_query_range = 'query_range' in merged_hit and self._is_valid_range(merged_hit['query_range'])
                has_template_range = 'template_range' in merged_hit and self._is_valid_range(merged_hit['template_range'])
                
                # Calculate missing ranges from alignment blocks if needed
                if not has_query_range or not has_template_range:
                    query_ranges = []
                    template_ranges = []
                    
                    # Extract ranges from alignment blocks
                    for block in alignment_data['alignment_blocks']:
                        # Calculate query range
                        if not has_query_range and 'query_seq' in block and 'query_start' in block:
                            ranges = self._calculate_range_from_alignment(
                                block['query_seq'], 
                                block['query_start']
                            )
                            query_ranges.extend(ranges)
                        
                        # Calculate template range
                        if not has_template_range and 'template_seq' in block and 'template_start' in block:
                            ranges = self._calculate_range_from_alignment(
                                block['template_seq'], 
                                block['template_start']
                            )
                            template_ranges.extend(ranges)
                    
                    # Process and store calculated ranges
                    if query_ranges and not has_query_range:
                        # Sort and merge overlapping ranges
                        merged_ranges = self._merge_overlapping_ranges(query_ranges)
                        # Store as the primary range - don't use "calculated" suffix
                        merged_hit['query_range'] = ','.join([f"{start}-{end}" for start, end in merged_ranges])
                    
                    if template_ranges and not has_template_range:
                        # Sort and merge overlapping ranges
                        merged_ranges = self._merge_overlapping_ranges(template_ranges)
                        # Store as the primary range - don't use "calculated" suffix
                        merged_hit['template_range'] = ','.join([f"{start}-{end}" for start, end in merged_ranges])
                
                # Calculate additional alignment statistics
                additional_stats = self._calculate_alignment_statistics(alignment_data['alignment_blocks'])
                for key, value in additional_stats.items():
                    if key not in merged_hit:
                        merged_hit[key] = value
                
                # Remove any redundant calculated ranges
                if 'query_range_calculated' in merged_hit:
                    merged_hit.pop('query_range_calculated')
                if 'template_range_calculated' in merged_hit:
                    merged_hit.pop('template_range_calculated')
            
            merged_hits.append(merged_hit)
        
        return merged_hits

    def _is_valid_range(self, range_str: str) -> bool:
        """
        Check if a range string is properly formatted
        
        Args:
            range_str: Range string to check
            
        Returns:
            True if valid, False otherwise
        """
        if not range_str:
            return False
        
        # Check format: "start-end,start-end,..."
        range_pattern = r'^\d+-\d+(?:,\d+-\d+)*$'
        return bool(re.match(range_pattern, range_str))

    def _merge_overlapping_ranges(self, ranges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """
        Merge overlapping ranges
        
        Args:
            ranges: List of (start, end) tuples
            
        Returns:
            List of merged (start, end) tuples
        """
        if not ranges:
            return []
        
        # Sort ranges by start position
        sorted_ranges = sorted(ranges, key=lambda x: x[0])
        
        merged = []
        current = sorted_ranges[0]
        
        for next_range in sorted_ranges[1:]:
            # Check if ranges overlap or are adjacent
            if next_range[0] <= current[1] + 1:
                # Merge ranges
                current = (current[0], max(current[1], next_range[1]))
            else:
                # No overlap, add current to result and move to next
                merged.append(current)
                current = next_range
        
        # Add the last range
        merged.append(current)
        
        return merged