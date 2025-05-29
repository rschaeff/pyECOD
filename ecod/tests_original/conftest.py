#!/usr/bin/env python3
"""
Test fixtures for domain analysis module
"""
import os
import pytest
import tempfile
import shutil
from unittest.mock import MagicMock

from ecod.config import ConfigManager
from ecod.db.manager import DBManager


@pytest.fixture
def temp_test_dir():
    """Create a temporary directory for test files"""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def mock_config():
    """Mock configuration for testing"""
    return {
        'database': {
            'host': 'localhost',
            'port': 5432,
            'user': 'test_user',
            'password': 'test_password',
            'database': 'test_db'
        },
        'paths': {
            'output_dir': '/tmp/ecod_test'
        },
        'reference': {
            'current_version': 'develop291',
            'domain_db': '/data/ref/domain_db',
            'chain_db': '/data/ref/chain_db'
        },
        'tools': {
            'blast_path': '/usr/bin/blast',
            'hhblits_path': '/usr/bin/hhblits',
            'hhsearch_path': '/usr/bin/hhsearch',
            'hhmake_path': '/usr/bin/hhmake'
        },
        'force_overwrite': True
    }


@pytest.fixture
def mock_config_manager(mock_config):
    """Mock ConfigManager instance"""
    config_manager = MagicMock(spec=ConfigManager)
    config_manager.config = mock_config
    config_manager.get_db_config.return_value = mock_config['database']
    return config_manager


@pytest.fixture
def mock_db_manager():
    """Mock DBManager instance"""
    db_manager = MagicMock(spec=DBManager)
    
    # Setup common database responses
    db_manager.execute_dict_query.return_value = []
    
    return db_manager


@pytest.fixture
def sample_blast_xml(temp_test_dir):
    """Create a sample BLAST XML file for testing"""
    xml_content = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.10.0+</BlastOutput_version>
  <BlastOutput_db>/data/ecod/blast_db/ecod.chain.fa</BlastOutput_db>
  <BlastOutput_query-def>1abc_A</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
          <Hit_def>2xyz B</Hit_def>
          <Hit_len>95</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>180.6</Hsp_bit-score>
              <Hsp_score>456</Hsp_score>
              <Hsp_evalue>1e-50</Hsp_evalue>
              <Hsp_query-from>5</Hsp_query-from>
              <Hsp_query-to>95</Hsp_query-to>
              <Hsp_hit-from>3</Hsp_hit-from>
              <Hsp_hit-to>93</Hsp_hit-to>
              <Hsp_align-len>90</Hsp_align-len>
              <Hsp_qseq>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ</Hsp_qseq>
              <Hsp_hseq>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ</Hsp_hseq>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
    
    file_path = os.path.join(temp_test_dir, 'blast_result.xml')
    with open(file_path, 'w') as f:
        f.write(xml_content)
        
    return file_path


@pytest.fixture
def sample_domain_blast_xml(temp_test_dir):
    """Create a sample domain BLAST XML file for testing"""
    xml_content = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.10.0+</BlastOutput_version>
  <BlastOutput_db>/data/ecod/blast_db/ecod.domain.fa</BlastOutput_db>
  <BlastOutput_query-def>1abc_A</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
          <Hit_def>d2xyzB1 2xyz:B</Hit_def>
          <Hit_len>50</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>95.5</Hsp_bit-score>
              <Hsp_score>237</Hsp_score>
              <Hsp_evalue>1e-25</Hsp_evalue>
              <Hsp_query-from>5</Hsp_query-from>
              <Hsp_query-to>50</Hsp_query-to>
              <Hsp_hit-from>3</Hsp_hit-from>
              <Hsp_hit-to>48</Hsp_hit-to>
              <Hsp_align-len>45</Hsp_align-len>
              <Hsp_qseq>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRST</Hsp_qseq>
              <Hsp_hseq>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRST</Hsp_hseq>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|1</Hit_id>
          <Hit_def>d2xyzB2 2xyz:B</Hit_def>
          <Hit_len>45</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>85.5</Hsp_bit-score>
              <Hsp_score>210</Hsp_score>
              <Hsp_evalue>1e-20</Hsp_evalue>
              <Hsp_query-from>60</Hsp_query-from>
              <Hsp_query-to>95</Hsp_query-to>
              <Hsp_hit-from>5</Hsp_hit-from>
              <Hsp_hit-to>40</Hsp_hit-to>
              <Hsp_align-len>35</Hsp_align-len>
              <Hsp_qseq>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHI</Hsp_qseq>
              <Hsp_hseq>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHI</Hsp_hseq>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""
    
    file_path = os.path.join(temp_test_dir, 'domain_blast_result.xml')
    with open(file_path, 'w') as f:
        f.write(xml_content)
        
    return file_path


@pytest.fixture
def sample_blast_summary_xml(temp_test_dir):
    """Create a sample BLAST summary XML file for testing"""
    xml_content = """<?xml version="1.0" encoding="utf-8"?>
<blast_summ_doc>
  <blast_summ pdb="1abc" chain="A">
    <chain_blast_run program="blastp" version="BLASTP 2.10.0+">
      <blast_db>/data/ecod/blast_db/ecod.chain.fa</blast_db>
      <blast_query>1abc_A</blast_query>
      <query_len>100</query_len>
      <hits>
        <hit num="1" pdb_id="2xyz" chain_id="B" hsp_count="1" evalues="1e-50">
          <query_reg>5-95</query_reg>
          <hit_reg>3-93</hit_reg>
          <query_seq>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ</query_seq>
          <hit_seq>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ</hit_seq>
        </hit>
      </hits>
    </chain_blast_run>
    <blast_run program="blastp" version="BLASTP 2.10.0+">
      <blast_db>/data/ecod/blast_db/ecod.domain.fa</blast_db>
      <blast_query>1abc_A</blast_query>
      <query_len>100</query_len>
      <hits>
        <hit num="1" domain_id="d2xyzB1" pdb_id="2xyz" chain_id="B" hsp_count="1" evalues="1e-25">
          <query_reg>5-50</query_reg>
          <hit_reg>3-48</hit_reg>
        </hit>
        <hit num="2" domain_id="d2xyzB2" pdb_id="2xyz" chain_id="B" hsp_count="1" evalues="1e-20">
          <query_reg>60-95</query_reg>
          <hit_reg>5-40</hit_reg>
        </hit>
      </hits>
    </blast_run>
    <hh_run program="hhsearch" db="hora_full">
      <hits>
        <hit num="1" domain_id="d2xyzB1" hh_prob="95.5" hh_score="150.2" hit_cover="0.9">
          <query_reg>5-50</query_reg>
          <hit_reg>3-48</hit_reg>
        </hit>
        <hit num="2" domain_id="d2xyzB2" hh_prob="92.3" hh_score="145.8" hit_cover="0.85">
          <query_reg>60-95</query_reg>
          <hit_reg>5-40</hit_reg>
        </hit>
      </hits>
    </hh_run>
    <self_comp_run programs="dali">
      <hits>
        <hit aligner="dali" z_score="8.5">
          <query_reg>5-50</query_reg>
          <hit_reg>60-95</hit_reg>
        </hit>
      </hits>
    </self_comp_run>
  </blast_summ>
</blast_summ_doc>"""
    
    file_path = os.path.join(temp_test_dir, 'blast_summary.xml')
    with open(file_path, 'w') as f:
        f.write(xml_content)
        
    return file_path


@pytest.fixture
def sample_domain_xml(temp_test_dir):
    """Create a sample domain XML file for testing"""
    xml_content = """<?xml version="1.0" encoding="utf-8"?>
<domain_doc pdb="1abc" chain="A" reference="develop291">
  <domain_list>
    <domain pdb="1abc" chain="A" range="5-50" t_group="a" h_group="b" x_group="c" a_group="d">
      <range>5-50</range>
      <evidence>
        <match domain_id="d2xyzB1" type="blast" evalue="1e-25">
          <query_range>5-50</query_range>
          <hit_range>3-48</hit_range>
        </match>
        <match domain_id="d2xyzB1" type="hhsearch" probability="95.5">
          <query_range>5-50</query_range>
          <hit_range>3-48</hit_range>
        </match>
      </evidence>
    </domain>
    <domain pdb="1abc" chain="A" range="60-95" t_group="e" h_group="f" x_group="g" a_group="h">
      <range>60-95</range>
      <evidence>
        <match domain_id="d2xyzB2" type="blast" evalue="1e-20">
          <query_range>60-95</query_range>
          <hit_range>5-40</hit_range>
        </match>
        <match domain_id="d2xyzB2" type="hhsearch" probability="92.3">
          <query_range>60-95</query_range>
          <hit_range>5-40</hit_range>
        </match>
      </evidence>
    </domain>
  </domain_list>
</domain_doc>"""
    
    file_path = os.path.join(temp_test_dir, 'domains.xml')
    with open(file_path, 'w') as f:
        f.write(xml_content)
        
    return file_path


@pytest.fixture
def sample_fasta_file(temp_test_dir):
    """Create a sample FASTA file for testing"""
    fasta_content = """>1abc_A
ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ"""
    
    file_path = os.path.join(temp_test_dir, '1abc_A.fa')
    with open(file_path, 'w') as f:
        f.write(fasta_content)
        
    return file_path


@pytest.fixture
def sample_self_comparison_xml(temp_test_dir):
    """Create a sample self-comparison XML file for testing"""
    xml_content = """<?xml version="1.0" encoding="utf-8"?>
<self_comparison>
  <structural_repeat aligner="dali" zscore="8.5">
    <ref_range>5-50</ref_range>
    <mob_range>60-95</mob_range>
  </structural_repeat>
  <sequence_repeat_set aligner="hhrepid" type="alpha">
    <sequence_repeat>
      <range>10-30</range>
    </sequence_repeat>
    <sequence_repeat>
      <range>40-60</range>
    </sequence_repeat>
  </sequence_repeat_set>
</self_comparison>"""
    
    file_path = os.path.join(temp_test_dir, '1abc_A.self_comp.xml')
    with open(file_path, 'w') as f:
        f.write(xml_content)
        
    return file_path


@pytest.fixture
def sample_hhsearch_xml(temp_test_dir):
    """Create a sample HHSearch summary XML file for testing"""
    xml_content = """<?xml version="1.0" encoding="utf-8"?>
<hh_summ_doc>
  <metadata>
    <pdb_id>1abc</pdb_id>
    <chain_id>A</chain_id>
    <reference>develop291</reference>
    <creation_date>2023-01-01 12:00:00</creation_date>
  </metadata>
  <hh_hit_list>
    <hh_hit hit_num="1" hh_prob="95.5" hh_score="150.2" hh_evalue="1e-20" ecod_domain_id="d2xyzB1">
      <query_range>5-50</query_range>
      <template_seqid_range coverage="0.9">3-48</template_seqid_range>
      <alignment>
        <query_ali>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRST</query_ali>
        <template_ali>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRST</template_ali>
      </alignment>
    </hh_hit>
    <hh_hit hit_num="2" hh_prob="92.3" hh_score="145.8" hh_evalue="1e-18" ecod_domain_id="d2xyzB2">
      <query_range>60-95</query_range>
      <template_seqid_range coverage="0.85">5-40</template_seqid_range>
      <alignment>
        <query_ali>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHI</query_ali>
        <template_ali>ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHI</template_ali>
      </alignment>
    </hh_hit>
  </hh_hit_list>
</hh_summ_doc>"""
    
    file_path = os.path.join(temp_test_dir, '1abc_A.develop291.hh_summ.xml')
    with open(file_path, 'w') as f:
        f.write(xml_content)
        
    return file_path


@pytest.fixture
def sample_protein_chains():
    """Create sample protein chain data"""
    return [
        {
            'id': 1,
            'pdb_id': '1abc',
            'chain_id': 'A',
            'source_id': '1abc_A',
            'length': 100,
            'sequence': 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ'
        },
        {
            'id': 2,
            'pdb_id': '2xyz',
            'chain_id': 'B',
            'source_id': '2xyz_B',
            'length': 95,
            'sequence': 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUV'
        }
    ]


@pytest.fixture
def sample_domain_classifications():
    """Create sample domain classification data"""
    return [
        {
            'domain_id': 'd2xyzB1',
            'ecod_uid': 1001,
            't_group': 'a',
            'h_group': 'b',
            'x_group': 'c',
            'a_group': 'd',
            'is_manual_rep': True,
            'is_f70': True,
            'is_f40': False,
            'is_f99': True
        },
        {
            'domain_id': 'd2xyzB2',
            'ecod_uid': 1002,
            't_group': 'e',
            'h_group': 'f',
            'x_group': 'g',
            'a_group': 'h',
            'is_manual_rep': False,
            'is_f70': True,
            'is_f40': False,
            'is_f99': True
        }
    ]


@pytest.fixture
def mock_batch_info():
    """Create mock batch information"""
    return {
        'id': 1,
        'batch_name': 'test_batch_20230101_1200',
        'base_path': '/tmp/ecod_test/batches/test_batch_20230101_1200',
        'type': 'blast',
        'ref_version': 'develop291',
        'total_items': 10,
        'completed_items': 5,
        'status': 'processing',
        'created_at': '2023-01-01 12:00:00',
        'completed_at': None
    }


@pytest.fixture
def setup_test_environment(temp_test_dir, sample_fasta_file, sample_blast_xml, sample_domain_blast_xml, 
                         sample_blast_summary_xml, sample_domain_xml, sample_self_comparison_xml, 
                         sample_hhsearch_xml):
    """Setup a complete test environment with all necessary files"""
    # Create directory structure
    pdb_chain_dir = os.path.join(temp_test_dir, '1abc_A')
    os.makedirs(pdb_chain_dir, exist_ok=True)
    
    chain_blast_dir = os.path.join(temp_test_dir, 'chain_blast_results')
    os.makedirs(chain_blast_dir, exist_ok=True)
    
    domain_blast_dir = os.path.join(temp_test_dir, 'domain_blast_results')
    os.makedirs(domain_blast_dir, exist_ok=True)
    
    # Copy sample files to appropriate locations
    shutil.copy(sample_fasta_file, os.path.join(pdb_chain_dir, '1abc_A.fa'))
    shutil.copy(sample_blast_xml, os.path.join(chain_blast_dir, '1abc_A.chainwise_blast.xml'))
    shutil.copy(sample_domain_blast_xml, os.path.join(domain_blast_dir, '1abc_A.domain_blast.xml'))
    shutil.copy(sample_blast_summary_xml, os.path.join(pdb_chain_dir, '1abc_A.develop291.blast_summ.xml'))
    shutil.copy(sample_domain_xml, os.path.join(pdb_chain_dir, 'domains_v12.1abc_A.develop291.xml'))
    shutil.copy(sample_self_comparison_xml, os.path.join(pdb_chain_dir, '1abc_A.self_comp.xml'))
    shutil.copy(sample_hhsearch_xml, os.path.join(pdb_chain_dir, '1abc_A.develop291.hh_summ.xml'))
    
    return temp_test_dir
="A" reference="develop291">
  <domain_list>
    <domain pdb="1abc" chain