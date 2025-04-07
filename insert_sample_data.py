#!/usr/bin/env python3
"""
insert_sample_data.py - Script to insert sample protein data into the PyECOD database
"""

import os
import sys
import argparse
import logging
from pathlib import Path
import hashlib

# Add parent directory to path if needed
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from ecod.core.config import ConfigManager
from ecod.core.db_manager import DBManager

def setup_logging(verbose: bool = False):
    """Configure logging"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    )

def calculate_md5(sequence: str) -> str:
    """Calculate MD5 hash of a sequence"""
    return hashlib.md5(sequence.encode()).hexdigest()

def insert_sample_protein(db: DBManager, pdb_id: str, chain_id: str, sequence: str) -> int:
    """Insert a sample protein and its sequence"""
    logger = logging.getLogger("ecod.insert_sample")
    
    # Create the source_id (format used in the code is pdb_id_chain_id)
    source_id = f"{pdb_id}_{chain_id}"
    
    # Check if protein already exists
    check_query = "SELECT id FROM ecod_schema.protein WHERE source_id = %s"
    existing = db.execute_query(check_query, (source_id,))
    
    if existing:
        logger.info(f"Protein {source_id} already exists with ID {existing[0][0]}")
        return existing[0][0]
    
    # Insert the protein
    logger.info(f"Inserting new protein {source_id}")
    protein_id = db.insert(
        "ecod_schema.protein",
        {
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "source_id": source_id,
            "length": len(sequence)
        },
        "id"
    )
    
    # Calculate MD5 hash of the sequence
    md5_hash = calculate_md5(sequence)
    
    # Insert the sequence
    db.insert(
        "ecod_schema.protein_sequence",
        {
            "protein_id": protein_id,
            "sequence": sequence,
            "md5_hash": md5_hash
        }
    )
    
    logger.info(f"Inserted protein {source_id} with ID {protein_id} and its sequence")
    return protein_id

def main():
    parser = argparse.ArgumentParser(description='Insert Sample PyECOD Protein Data')
    parser.add_argument('--config', type=str, default='config/config.yml', 
                      help='Path to configuration file')
    parser.add_argument('--pdb-id', type=str, default='1abc',
                      help='PDB ID of the sample protein')
    parser.add_argument('--chain-id', type=str, default='A',
                      help='Chain ID of the sample protein')
    parser.add_argument('--sequence', type=str,
                      help='Amino acid sequence. If not provided, a sample sequence will be used.')
    parser.add_argument('-v', '--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    setup_logging(args.verbose)
    
    # Default sample sequence if not provided
    if not args.sequence:
        args.sequence = "MKKRLTESQFQEGTARIMSSFGLSKKGVANITVAVTNGYSGKQLISQQEDITKSLLVTELGSSRTPETVRMVLSNMEKLSSADFVFLTADQVEEKILSKHQGVQVLESTLMTIPDTPASMLLKLTDVPAKLKEALDATLNDKEFLERKAQEKEAFLEALKTMTKDLPRAQKAAFEEQEKQWLEAAQVAPATVAPQFPTVEAAQELKDLGKEKLSAVTFADLNEALEHHDRKLLSSVQSRAFTTPTNTEMSALDIEMDKVREADLEKATKEAQSIVNMVSKQAIIRQRDLNGFLEKLNKQEVQKSVQEILERKQQVLDAMEAREAALQSEAVTATKQALINKIKQQQKAYKEALEDMVAVKLKKAEAENAQKQLAEQTAKSLQMASEEGIEQLEQQLGMKVSELQASQEELTNLSAKIAEAQGDVKGAASILGHLNDQKQALAEEQKNAEELKAEVATIVADLKKAFEKLEHEYEALKEELTKLKEQLEVQAKDVERALSSAEETVKMLEQKLNVKAEQEKLELEKLRLQLEEAKRQAEEVRVQVEQEKRAQLQKLRAELDKTAQDLTSVKEEAEKAAEQLNKTLSEKEAALSTAQSQLEGLLNTKAELEAELAEAQAKLEEQYAVLSKRAADAEAEANTKALEAVKGLERDVNQSQKSLEAQVGTLKGKVGELASAIESLKQRLKEVEATRQAVQQAEDALRELQSSLDQLKAAKTELEQTVSQLESQAKEQADSLQEVQQARDQAETEKYKTEIDQMKRQLEETNTSLEAELIKIRNELQEEKKKLDEETQKRIEELEAELAKAQEKLQAEAEAKVDLENKTKQELQAAIKDVESQKQSLAAVNDAKTRLEQVNLEREEVGRLKSQLEEARIRLLEISEEKLAQVASGKDELLKHIQSLNLEKDAALENRLRELAKVKGELEAATEELEKKIAQLEQEKRALQSQLEEAKRSLATALSEAEARADQLKQVRTEFRRQLEEAQNAADASRQLQEEMEDAKRQLQEAEAAVAALKAKLDAIQGEITSLTAQIAELEAELRAVQSEIDRKIAEKDEEIQVLRAQVEELMTAKAEAERALKQAEELRQQLEIVRNQLDEAKRQLQQAEEHVGTLQALQAELATAKAEAEVALRAMENEIKQLKQQVEEMKTEAEQTQKALAEANGLLDQLRAQLEQADMARDTLLKSASVIAEAQEKLQAELESVRTALATAVSEAEARAAQLRQEVLQLRAQLEEAERAATAANELQKSLQRTSAELDRTKQTISQLQAEVEELRATVEQAEADKQQLDEKLKQVTAQLSMEQALAAAKDAISTAEALFQTKAELESKLAQAEEQLRDAETVAQRVLRAQEEAIQEKDLMQQKLQDVMALQGEVEDLRTKLAETNAALQAERDRADQALAELQKEIAATEDALRQARAALEAQEAKLGRLEEERGLLEREAAQLQTKAESLAASKSQLEQAVSALKAEAEAALIKTKDEIAVQTKAFDQVRSELENTLVATNAAIDEQKELWKQLREASEKAQAVLTQLEEAEEALQKAQAAAAELESLKGKFSDASEALKEQLLDVQHRLDEAEQNAKKAEDQLKKAQEEIDAAQAQLQEAEAKVDVSIKQIEAEKKEIEVLKASKESIAVTKTRIDLDESLKAATDRLDEAESRAASEKQELTESLKSAVDQTDRLIDELEKVRQEKDEAEAQAKQLRQDLTELNDRLAEAESRAASEKQELTESLKSAVDQTDRLIDELEKVRQEKDEAEAQAKQLRQDLTELNDRLAEAESRAASEKQELTESLKSAVDQTDRLIDELEKVRQEKDEAEAQAKQLRQDLTELNDRLAEAESRAASEKQELTESLKSAVDQTDRLIDELEKVRQEKDEAEAQAKQLRQDLTELNDRLAEAESRAAEEKQELTDSLKSA"
    
    config_manager = ConfigManager(args.config)
    db_config = config_manager.get_db_config()
    db = DBManager(db_config)
    
    try:
        protein_id = insert_sample_protein(db, args.pdb_id, args.chain_id, args.sequence)
        print(f"Successfully inserted sample protein with ID: {protein_id}")
    except Exception as e:
        logging.error(f"Error inserting sample protein: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()