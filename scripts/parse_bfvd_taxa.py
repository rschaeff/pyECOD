import csv
import json
import re
import psycopg2
from datetime import datetime

def parse_bfvd_tsv(input_file, output_file):
    """
    Parse the BFVD TSV file containing viral taxonomies and prepare for DB upload
    """
    proteins = []
    
    with open(input_file, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        
        for row in reader:
            if len(row) < 5:
                continue
            
            pdb_file = row[0]
            taxid = int(row[1])
            rank = row[2]
            scientific_name = row[3]
            lineage_string = row[4]
            
            # Extract UniProt accession from filename
            unp_acc = pdb_file.split('_')[0]
            
            # Parse lineage to extract family, genus, species
            parsed_lineage = parse_lineage(lineage_string)
            standard_ranks = parsed_lineage["standard"]
            
            # Extract taxonomy info for the existing table
            taxonomy_record = {
                "unp_acc": unp_acc,
                "virus_family": standard_ranks.get("family", ""),
                "virus_genus": standard_ranks.get("genus", ""),
                "virus_species": standard_ranks.get("species", ""),
                "host_organism": "",  # Not available in the TSV
                "host_tax_id": None,  # Not available in the TSV
                "taxid": taxid,  # Additional field not in current schema
                "rank": rank,    # Additional field not in current schema
                "scientific_name": scientific_name  # Additional field not in current schema
            }
            
            # Create a detailed lineage record for the new table
            lineage_record = {
                "unp_acc": unp_acc,
                "taxid": taxid,
                "rank": rank,
                "scientific_name": scientific_name,
                "full_lineage": lineage_string,
                "domain": standard_ranks.get("domain", ""),
                "kingdom": standard_ranks.get("kingdom", ""),
                "phylum": standard_ranks.get("phylum", ""),
                "class": standard_ranks.get("class", ""),
                "order": standard_ranks.get("order", ""),
                "family": standard_ranks.get("family", ""),
                "genus": standard_ranks.get("genus", ""),
                "species": standard_ranks.get("species", ""),
                "informal_ranks": parsed_lineage["informal"]
            }
            
            # Create protein model record
            model_info = extract_model_info(pdb_file)
            model_record = {
                "unp_acc": unp_acc,
                "source_id": pdb_file,  # Using PDB filename as source_id
                "avg_plddt": None,  # Not available in the TSV
                "ptm_score": None,  # Not available in the TSV
                "version": model_info.get("model_version", ""),
                "file_path": pdb_file,
                "file_size": None,  # Not available in the TSV
                "file_date": datetime.now().isoformat(),
                "is_current": True
            }
            
            proteins.append({
                "protein_data": {
                    "unp_acc": unp_acc,
                    "model_count": 1,
                    "is_split": False,
                    "sequence_length": None,  # Not available in TSV
                    "sequence": None,         # Not available in TSV
                    "sequence_md5": None      # Not available in TSV
                },
                "model_data": model_record,
                "taxonomy_data": taxonomy_record,
                "lineage_data": lineage_record
            })
    
    # Write to JSON file
    with open(output_file, 'w') as json_file:
        json.dump({"proteins": proteins}, json_file, indent=2)
    
    return len(proteins)

def extract_model_info(filename):
    """Extract model information from the AlphaFold PDB filename"""
    # Example: A0A514CX81_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb
    parts = filename.split('_')
    
    # The issue was with this part - we need to handle the seed more carefully
    seed_value = 0
    if len(parts) > 9:
        # Find the seed part and extract the number correctly
        seed_part = parts[9].split('.')[0]
        # Only try to convert to int if it's actually digits
        if seed_part.isdigit():
            seed_value = int(seed_part)
    
    model_info = {
        "uniprot_id": parts[0],
        "relaxation_status": parts[1] if len(parts) > 1 else "unknown",
        "rank": int(parts[3]) if len(parts) > 3 and parts[3].isdigit() else 0,
        "model_version": parts[4] if len(parts) > 4 else "unknown",
        "ptm_status": "ptm" if len(parts) > 5 and parts[5] == "ptm" else "non_ptm",
        "model_number": int(parts[7]) if len(parts) > 7 and parts[7].isdigit() else 0,
        "seed": seed_value
    }
    
    return model_info

def parse_lineage(lineage_string):
    """
    Parse the lineage string into a structured format,
    separating standard and non-standard ranks
    """
    # Split the lineage by semicolons
    lineage_parts = lineage_string.split(';')
    
    # Standard taxonomic ranks
    standard_ranks = {
        "d_": "domain",
        "k_": "kingdom",
        "p_": "phylum",
        "c_": "class",
        "o_": "order",
        "f_": "family",
        "g_": "genus",
        "s_": "species"
    }
    
    # Initialize results
    parsed = {
        "standard": {},
        "informal": []
    }
    
    # Process each part of the lineage
    for part in lineage_parts:
        part = part.strip()
        if not part:
            continue
            
        # Check if this is a standard rank
        is_standard = False
        for prefix, rank_name in standard_ranks.items():
            if part.startswith(prefix):
                parsed["standard"][rank_name] = part[len(prefix):]
                is_standard = True
                break
        
        # If not a standard rank, it's an informal designation
        if not is_standard:
            if part.startswith('-_'):
                # This is explicitly marked as informal
                informal_name = part[2:]  # Remove the '-_' prefix
                category = "subfamily" if "virinae" in informal_name else "unknown"
                parsed["informal"].append({
                    "value": informal_name,
                    "type": category
                })
            else:
                # Unknown format - store as is
                parsed["informal"].append({
                    "value": part,
                    "type": "unknown"
                })
    
    return parsed