import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import pandas as pd
import numpy as np
import ms_module as ms
#
#
pep_fname = "../peptides.xls"
spec_fname = "../specs.xls"
fasta_fname = "dec23_1701.fasta"
#
pep_info = pd.read_csv(pep_fname,sep='\t')
# cols = ['Protein accession numbers','Assigned','Other Proteins']
spec_info = pd.read_csv(spec_fname,sep='\t')
#
fasta = SeqIO.to_dict(SeqIO.parse(fasta_fname,"fasta"),key_function=lambda _: _.id.split('|')[1])
#

























###########################################
# 'Protein accession numbers'
# 'Exclusive unique peptide count'
# 'Exclusive unique spectrum count'
# 'Total spectrum count'
# 'Spectrum name'
# 'Peptide sequence'
# 'Mascot Ion score'
# 'Mascot Identity score'
# 'Mascot Delta Ion Score'
# 'Number of enzymatic termini'
# 'Fixed modifications identified by spectrum'
# 'Variable modifications identified by spectrum'
# 'Observed m/z'
# 'Actual peptide mass (AMU)'
# 'Calculated +1H Peptide Mass (AMU)'
# 'Spectrum charge'
# 'Actual minus calculated peptide mass (AMU)'
# 'Actual minus calculated peptide mass (PPM)'
# 'Retention Time (sec)'
# 'Precursor Intensity'
# 'Total Ion Current'
###########################################
# cc = ['Protein accession numbers',
# # 'Exclusive unique peptide count',
# # 'Exclusive unique spectrum count',
# # 'Total spectrum count',
# # 'Spectrum name',
# 'Peptide sequence',
# # 'Mascot Ion score',
# # 'Mascot Identity score',
# # 'Mascot Delta Ion Score',
# # 'Number of enzymatic termini',
# 'Fixed modifications identified by spectrum',
# 'Variable modifications identified by spectrum',
# # 'Observed m/z',
# # 'Actual peptide mass (AMU)',
# # 'Calculated +1H Peptide Mass (AMU)',
# # 'Spectrum charge',
# # 'Actual minus calculated peptide mass (AMU)',
# # 'Actual minus calculated peptide mass (PPM)',
# # 'Retention Time (sec)',
# # 'Precursor Intensity',
# # 'Total Ion Current'
# ]
###########################################

















