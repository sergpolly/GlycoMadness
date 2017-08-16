import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import pandas as pd
import numpy as np
import ms_module as ms
import re
############################
from Bio import Entrez
from Bio import SeqIO
from StringIO import StringIO
import time
from urllib2 import HTTPError  # for Python 2
#
import argparse

import warnings
from Bio import BiopythonWarning, BiopythonParserWarning


full_sorted_output_with_criteria = False

#
# do some arguments parsing to make the script looks civilized ...
parser = argparse.ArgumentParser()
parser.add_argument("-f","--pept_with_fetch",
    help="specify file name of peptide summary with fetchids (with/without path)",required=True)
parser.add_argument("-g","--genbank",
    help="specify file name of genbank records with pulled proteins (with/without path)",required=True)
# we don't need spectrum file for downloading proteins, it is too redundant for that purpose ...
# parser.add_argument("--verbose", help="verbose output", action="store_true")
parser.add_argument("--prefix", help="specify common part of the path for peptide and spectrum files")
parser.add_argument("--threshold", type=int, default=130, help="Threshold for protein per peptide choice, defaults at 130")
args = parser.parse_args()
#
###############################################
if args.prefix is not None:
    pep_fetch_fname = os.path.join( args.prefix, args.pept_with_fetch )
    gb_fname = os.path.join( args.prefix, args.genbank )
else:
    pep_fetch_fname = args.pept_with_fetch
    gb_fname = args.genbank
# get the common path for later use ...
common_path = os.path.commonprefix([pep_fetch_fname,gb_fname])
common_path = os.path.dirname(common_path)
#
# # don'r forget to provide you email
# Entrez.email = args.email if args.email else "your_email@mail_server.com"
crit_threshold = args.threshold
# Reading genbank mindfully next ...
# gbrecs = ms.genebank_fix_n_read(gb_fname)
gbrecs = ms.genebank_fix_n_read(gb_fname,key_func_type='id')
######################################
# assign some module internal stuff ...
ms.gbrecs = gbrecs


############################
# READING file containing GeneName(and/or locus) and FetchID association ...
print "Reading %s with the updated spectrum that includes fetchid column ..."%pep_fetch_fname
pep_fetch = pd.read_csv(pep_fetch_fname)


# here is the NEW plan!:
# first, we try to assign a single protein to each peptide
# we collect peptide-protein pairs that failed to match, declare them BAD and send them to manuall processing ...



#
############################################
#  columns that needs to be delivered ...  #
############################################
# A gsites, 1 per line
# B pept, 1 per line
# B1 enzyme, G or T, derive from 'Biological sample category', like this: {'TrypsinSample1':'T','GluC_Sample2':'G'}
# C peptide_start, 1 per line accordingly
# D all_uids, REPLACE WITH col:H
# E prot_seq, try to get those from NCBI, not from UniProt ...
# F protein, ??? sequence, name or what???
# G uid_max, UID for major form instead or something like that ...
# H prot_name, parsed out human-readable name from 'Protein name'
# H1 gene_name, parsed out GN=xxx from 'Protein name'
# I uniq_peptide_count, discrad that column ...
# J pept_probability, output number not the string - this would be the criteria 
# K gsites_predicted, OK
# L gsites_predicted_number, OK
# M gsite_start, beware of 0 or 1 type of indexing ...
# N,O,P - gsites AAs in separate columns
# M1, NOP combined, gsite sequence basically!
# Q signal, from GeneBank record on the protein, simply Y,N on whether there is a 'Signal' in gb.
# R signal_location, location of the signal from Q
# S tm_span, Y,N just for the fact of having TM span as a protein feature.
#
#
#
#
#####################################################################################################

# SIMPLIFY THING BY RENAMING SOME COLUMNS ...
col_rename = {
 'Peptide start index':'pept_start',
 'Peptide stop index':'pept_stop',
 'Previous amino acid':'prev_aa',
 'Next amino acid':'next_aa',
 'Peptide sequence':'pept',
}
# rename ...
pep_fetch.rename(columns=col_rename,inplace=True)
# simplified column names to be further used in the processing ...
cols_simple = [
 'enzyme',
 'prev_aa',
 'next_aa',
 'pept_start',
 'pept_stop',
 'pept',
 'GN',
 'OS',
 'locus',
 'prot_name',
 'uid',
 'fetchid',
 'fetchacc',
 'Best Mascot Ion Score',
 'Best Mascot Identity Score',
 'Best Mascot Delta Ion Score']
cols_simple += [
'uid_fetched',
'GN_fetched',
'signal',
'signal_loc',
'tm_span']
#
# FILL IN SOME COLUMNS ...
pep_fetch['enzyme']      = pep_fetch['Biological sample category'].apply(ms.get_enzyme)
# pep_fetch['prot_ident_probab'] = pep_fetch['Protein identification probability'].str.strip('%').apply(float)
# pep_fetch['pept_ident_probab'] = pep_fetch['Peptide identification probability'].str.strip('%').apply(float)
pep_fetch['uid_fetched'] = pep_fetch['fetchacc'].apply(lambda fidx: gbrecs[fidx].id if pd.notnull(fidx) else None)
# pep_fetch['uid_fetched'] = pep_fetch['fetchid'].apply(lambda fidx: gbrecs[str(int(fidx))].id if pd.notnull(fidx) else None)
pep_fetch['GN_fetched']  = pep_fetch['fetchacc'].apply( ms.get_genename )
#
pep_fetch['signal']      = pep_fetch['fetchacc'].apply( ms.get_signal )
pep_fetch['signal_loc']  = pep_fetch['fetchacc'].apply( ms.get_signal_loc )
pep_fetch['tm_span']     = pep_fetch['fetchacc'].apply( ms.get_tm_span )
#
# NOW SIMPLIFY THE DATAFRAME TO CONSIDER ONLY IMPORTNAT COLUMNS AND REMOVE DUPLICATES ...
pep_simple = pep_fetch[cols_simple].drop_duplicates()
#
#
#
#
#
############################################
#  CRITERIA FOR PERFECT PROTEIN MATCH...   #
############################################
# GN_fetched == GN (?)
# uid_fetched == uid (as split into tuples) #Q9UHG3-2
# uid_fetched == uid (major part only)
# pept is in protein sequence
# pept start and pept stop indeces MATCH 
# previous AA and next AA MATCH
#
pep_simple['crit_GN'] = (pep_simple['GN']==pep_simple['GN_fetched'])
# uid splitter - or . AMD uid major taker ...
uid_split = lambda uid: tuple(re.split('[-\.]',uid)) if pd.notnull(uid) else None
uid_major = lambda uid: re.split('[-\.]',uid)[0] if pd.notnull(uid) else None
pep_simple['crit_uid_full'] = (pep_simple['uid'].apply(uid_split)==pep_simple['uid_fetched'].apply(uid_split))
pep_simple['crit_uid_maj'] = (pep_simple['uid'].apply(uid_major)==pep_simple['uid_fetched'].apply(uid_major))
# pept_isin = lambda fidx: gbrecs[str(int(fidx))].id if pd.notnull(fidx) else None)


pep_simple['crit_pept_in'] = pep_simple[['pept','fetchacc']].apply(ms.pept_isin,axis=1)
# extract pept_info fetched first ...
###############################################################
###############################################################
pep_simple = pep_simple.merge( pep_simple[['pept','fetchacc']].apply(ms.pept_info,axis=1),left_index=True,right_index=True )
#
pep_simple['crit_start'] = pep_simple['pept_start'] == pep_simple['start_fetched']
pep_simple['crit_stop'] = pep_simple['pept_stop'] == pep_simple['stop_fetched']
pep_simple['crit_prev_aa'] = pep_simple['prev_aa'] == pep_simple['prev_aa_fetched']
pep_simple['crit_next_aa'] = pep_simple['next_aa'] == pep_simple['next_aa_fetched']
#
#
#
#
crit_cols = [cn for cn in pep_simple.columns if 'crit_' in cn]
# criteria weights ...
crit_weight = pd.Series({'crit_GN':10,
 'crit_uid_full':1,
 'crit_uid_maj':10,
 'crit_pept_in':100,
 'crit_start':1,
 'crit_stop':1,
 'crit_prev_aa':10,
 'crit_next_aa':10})
# get a weighted sum of all criteria ...
pep_simple['SCORE'] = pep_simple[crit_cols].fillna(False).mul(crit_weight).sum(axis=1)
#
# ###################################################################
# #  EVALUATE THESE CRITERIA AND GET A COLUMN WITH SUM(AXIS=1)...   #
# ###################################################################
# # THEN DECIDE HOW MANY GENEBANK PROTEIN RECORDS QUALIFY CRITERIA,
# # IF IT IS JUST 1, THEN PROCEED WITH THE ONE, ELSE PRINT ALL THE INFO FOR
# # FURTHER INVESTIGATION ...
cols = [
 'pept',
 'fetchid',
 'fetchacc',
 'GN',
 'GN_fetched',
 'prev_aa',
 'next_aa',
 'prev_aa_fetched',
 'next_aa_fetched',
 'pept_start',
 'pept_stop',
 'start_fetched',
 'stop_fetched',
 'enzyme',
 'OS',
 'locus',
 'prot_name',
 'uid',
 'uid_fetched',
 'signal',
 'signal_loc',
 'tm_span',
 'crit_GN',
 'crit_uid_full',
 'crit_uid_maj',
 'crit_pept_in',
 'crit_start',
 'crit_stop',
 'crit_prev_aa',
 'crit_next_aa',
 'Best Mascot Ion Score',
 'Best Mascot Identity Score',
 'Best Mascot Delta Ion Score']


cols_short = [
 'pept',
 'fetchid',
 'fetchacc',
 'enzyme',
 'GN',
 'GN_fetched',
 'prev_aa',
 'next_aa',
 'prev_aa_fetched',
 'next_aa_fetched',
 'pept_start',
 'pept_stop',
 'start_fetched',
 'stop_fetched',
 'locus',
 'prot_name',
 'uid',
 'uid_fetched',
 'signal',
 'signal_loc',
 'tm_span',
 'Best Mascot Ion Score',
 'Best Mascot Identity Score',
 'Best Mascot Delta Ion Score']


# output stuff ...
pep_simple_sorted = pep_simple.sort_values(by=['pept','SCORE'],inplace=False)[cols+['SCORE']]
if full_sorted_output_with_criteria:
    pep_simple_sorted.to_csv('FULL_SORTED_pepts_and_scores.csv',index=False)

# group data by peptides ...
pep_grouped = pep_simple_sorted.groupby(by='pept')

# let's choose a single protein per peptide:
# idea: of max SCORE in a peptide-group is above Threshold AND is the only one in a group, THEN
# we assign that protein to peptide.
# OTHERWISE, (multiple max_SCORE proteins, no proteins that qualify criteria by Threshold) 
crit_threshold = 110
num_qualify_prots = lambda pep_grp: (pep_grp==pep_grp.max()).sum() if pep_grp.max()>=crit_threshold else 0
idx_qualify_prots = lambda pep_grp:  pep_grp.idxmax()  if num_qualify_prots(pep_grp)==1 else None

# now we have to retrieve those peptide-protein combinations that didn't work out
# (ambiguous or lacking qualified protein)
print
print "Pept-protein map relies on the uniqness of the max_SCORE value ..."
print "Here are max_SCORE occurences in the sample:"
print str(pep_grouped['SCORE'].apply(num_qualify_prots).value_counts())
print "Those with number of occurences != 1  will go to BAD_PEPTS file ..."
print

# indexes of qualified pept-protein pairs, including peptides with no certail match as 'None' ...
qual_prot_idxs = pep_grouped['SCORE'].apply(idx_qualify_prots)

# separate BAD from qualified pept-protein pairs ...
pept_prot_idxs = qual_prot_idxs[qual_prot_idxs.notnull()].map(int)
BAD_PEPTS = qual_prot_idxs[qual_prot_idxs.isnull()].index

# GOOD ONE ARE OUT FIRST ...
# qualified pept-protein pairs goes straight to output for further use ...
pept_prot_map = pep_simple_sorted.loc[pept_prot_idxs][cols_short + ['SCORE','crit_pept_in']]
pept_prot_map.to_csv(os.path.join(common_path,'pept_prot_map.csv'),index=False)
# common_path


# FINISH IT UP WITH BAD ONES ...
bad_pept_prot_output = pd.concat( (pep_grouped.get_group(peptide) for peptide in BAD_PEPTS) ).reset_index(drop=True)[cols+['SCORE',]]
bad_pept_prot_output.to_csv(os.path.join(common_path,'BAD_pept_prot.csv'),index=False)


# STAGE 2 APPEARS TO BE WORKING PROPERLY ...











