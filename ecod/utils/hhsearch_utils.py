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
            if line.strip().startswith('No Hit'):
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

        # Parse hit lines using a more flexible approach
        for i in range(table_start, table_end):
            line = lines[i].strip()
            if not line:
                continue

            try:
                # Split the line by whitespace and extract fields
                parts = line.split()

                if len(parts) < 10:
                    self.logger.warning(f"Hit line has too few fields: {line}")
                    continue

                # Extract hit number and ID
                hit_num = int(parts[0])
                hit_id = parts[1]

                # Find where the numeric values start (probability)
                # This handles cases where the hit description might contain spaces
                numeric_start_idx = 2
                while numeric_start_idx < len(parts):
                    try:
                        # Try to parse as float to find the probability field
                        float(parts[numeric_start_idx])
                        break
                    except ValueError:
                        numeric_start_idx += 1

                # If we couldn't find the numeric fields, skip this line
                if numeric_start_idx >= len(parts) - 8:  # Need at least 8 more fields
                    self.logger.warning(f"Could not find numeric fields in hit line: {line}")
                    continue

                # Extract the remaining hit description before the probability
                description = " ".join(parts[2:numeric_start_idx])

                # Extract numeric fields - ALL AS FLOATS
                try:
                    probability = float(parts[numeric_start_idx])
                    e_value = self._parse_scientific(parts[numeric_start_idx + 1])
                    p_value = self._parse_scientific(parts[numeric_start_idx + 2])
                    score = float(parts[numeric_start_idx + 3])
                    ss_score = float(parts[numeric_start_idx + 4])  # This is a float, not an int!
                    cols = int(parts[numeric_start_idx + 5])  # This is the only int
                    query_range = parts[numeric_start_idx + 6]
                    template_range = parts[numeric_start_idx + 7]

                    # Extract template info in parentheses if available
                    template_length = None
                    remaining_parts = parts[(numeric_start_idx + 8):]
                    if remaining_parts and remaining_parts[-1].startswith('(') and remaining_parts[-1].endswith(')'):
                        try:
                            template_length = int(remaining_parts[-1].strip('()'))
                        except ValueError:
                            pass

                    # Create hit object
                    hit = {
                        'hit_num': hit_num,
                        'hit_id': hit_id,
                        'description': description,
                        'probability': probability,
                        'e_value': e_value,
                        'p_value': p_value,
                        'score': score,
                        'ss_score': ss_score,
                        'cols': cols,
                        'query_range': query_range,
                        'template_range': template_range
                    }

                    if template_length is not None:
                        hit['template_length'] = template_length

                    hits.append(hit)
                except (ValueError, IndexError) as e:
                    self.logger.warning(f"Error parsing numeric fields in hit line: {line}, error: {str(e)}, numeric_start_idx: {numeric_start_idx}")
                    continue

            except Exception as e:
                self.logger.warning(f"Could not parse hit line: {line}, reason: {str(e)}")

        return hits, table_end + 1
    
    # Improved method for parsing alignments
    def _parse_alignments(self, content: str, alignment_start_line: int, hit_summary: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """Parse the detailed alignment sections"""
        lines = content.split('\n')
        alignments = {}

        # Map hit numbers to IDs for easier lookup
        hit_id_map = {hit['hit_num']: hit['hit_id'] for hit in hit_summary}

        # State variables
        current_hit_id = None
        current_hit_num = None
        in_alignment_section = False
        current_alignment_block = None

        # Process lines
        i = alignment_start_line
        while i < len(lines):
            line = lines[i].strip()

            # Detect alignment section start - look for "No X" where X is the hit number
            match = re.match(r'^No\s+(\d+)$', line)
            if match:
                hit_num = int(match.group(1))
                if hit_num in hit_id_map:
                    # Save previous alignment block if exists
                    if current_hit_id and current_alignment_block:
                        if current_hit_id not in alignments:
                            alignments[current_hit_id] = {'alignment_blocks': []}
                        alignments[current_hit_id]['alignment_blocks'].append(current_alignment_block)

                    # Start new alignment block
                    current_hit_num = hit_num
                    current_hit_id = hit_id_map[hit_num]
                    in_alignment_section = True
                    current_alignment_block = {
                        'query_seq': "",
                        'template_seq': "",
                        'match_seq': "",
                        'query_consensus': "",
                        'template_consensus': "",
                        'stats': {},
                        'query_start': None,
                        'template_start': None
                    }

                    i += 1  # Skip the hit ID line
                    continue

            # Check for the probability line
            if in_alignment_section and line.startswith('Probab='):
                stats = {}
                for stat_pair in line.split():
                    if '=' in stat_pair:
                        key, value = stat_pair.split('=', 1)
                        try:
                            value = float(value.replace('%', '')) / 100.0 if '%' in value else float(value)
                        except ValueError:
                            pass
                        stats[key] = value

                if current_alignment_block:
                    current_alignment_block['stats'] = stats

            # Process alignment lines
            elif in_alignment_section and current_alignment_block:
                # Query sequence line
                if line.startswith('Q ') and ' ' in line:
                    parts = line.split()
                    if len(parts) >= 5:  # Q id start seq end
                        if parts[1] == "Consensus":
                            # This is a consensus line
                            current_alignment_block['query_consensus'] += parts[3]
                        else:
                            # This is a sequence line
                            current_alignment_block['query_seq'] += parts[3]
                            if current_alignment_block['query_start'] is None:
                                try:
                                    current_alignment_block['query_start'] = int(parts[2])
                                except ValueError:
                                    pass

                # Match line
                elif in_alignment_section and current_alignment_block and line and line[0] == ' ' and not line.startswith(' Q') and not line.startswith(' T'):
                    match_seq = line.strip()
                    current_alignment_block['match_seq'] += match_seq

                # Template sequence line
                elif line.startswith('T ') and ' ' in line:
                    parts = line.split()
                    if len(parts) >= 5:  # T id start seq end
                        if parts[1] == "Consensus":
                            # This is a consensus line
                            current_alignment_block['template_consensus'] += parts[3]
                        else:
                            # This is a sequence line
                            current_alignment_block['template_seq'] += parts[3]
                            if current_alignment_block['template_start'] is None:
                                try:
                                    current_alignment_block['template_start'] = int(parts[2])
                                except ValueError:
                                    pass

            # Check for end of alignment section (empty line followed by a new No X or end of file)
            if in_alignment_section and not line:
                # Look ahead to see if this is the end of the section
                if i + 1 < len(lines) and (lines[i+1].startswith('No ') or not lines[i+1].strip()):
                    # Save current alignment block
                    if current_hit_id and current_alignment_block:
                        if current_hit_id not in alignments:
                            alignments[current_hit_id] = {'alignment_blocks': []}
                        alignments[current_hit_id]['alignment_blocks'].append(current_alignment_block)

                        # Reset for next section
                        in_alignment_section = False
                        current_alignment_block = None

            i += 1

        # Add the last alignment block if it exists
        if current_hit_id and current_alignment_block:
            if current_hit_id not in alignments:
                alignments[current_hit_id] = {'alignment_blocks': []}
            alignments[current_hit_id]['alignment_blocks'].append(current_alignment_block)

        return alignments
    
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
