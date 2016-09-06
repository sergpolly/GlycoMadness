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
parser.add_argument("-f","--raw_fetch", help="speicfy input data fname (with fetchid column(!), with/without path)",required=True)
# parser.add_argument("-f","--pept_with_fetch",
#     help="specify file name of peptide summary with fetchids (with/without path)",required=True)
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
    raw_fetch_fname = os.path.join( args.prefix, args.raw_fetch )
    gb_fname = os.path.join( args.prefix, args.genbank )
else:
    raw_fetch_fname = args.raw_fetch
    gb_fname = args.genbank
# get the common path for later use ...
common_path = os.path.commonprefix([raw_fetch_fname,gb_fname])
common_path = os.path.dirname(common_path)
#
# # don'r forget to provide you email
# Entrez.email = args.email if args.email else "your_email@mail_server.com"
crit_threshold = args.threshold
# Reading genbank mindfully next ...
gbrecs = ms.genebank_fix_n_read(gb_fname)
######################################
# assign some module internal stuff ...
ms.gbrecs = gbrecs


############################
# READING file containing GeneName(and/or locus) and FetchID association ...
print "Reading %s with the updated spectrum that includes fetchid column ..."%raw_fetch_fname
raw_fetch = pd.read_csv(raw_fetch_fname)


# here is the NEW plan!:
# first, we try to assign a single protein to each peptide
# we collect peptide-protein pairs that failed to match, declare them BAD and send them to manuall processing ...
#####################################################################################################

# SIMPLIFY THING BY RENAMING SOME COLUMNS ...
col_rename = {
 # 'Peptide start index':'pept_start',
 # 'Peptide stop index':'pept_stop',
 # 'Previous amino acid':'prev_aa',
 # 'Next amino acid':'next_aa',
 'MS Sample':'exp_name',
 'Spectrum Name':'spec_name',
 'Peptide Sequence':'pept',
  # Channel 1  Channel 2   Normalized Channel 1    Normalized Channel 2
 'Channel 1':'ch1',
 'Channel 2':'ch2',
 'Normalized Channel 1':'norm_ch1',
 'Normalized Channel 2':'norm_ch2'}
# rename ...
raw_fetch.rename(columns=col_rename,inplace=True)
# simplified column names to be further used in the processing ...
cols_simple = [
 'enzyme',
 # 'prev_aa',
 # 'next_aa',
 # 'pept_start',
 # 'pept_stop',
 'spec_name',
 'exp_name',
 'pept',
 'GN',
 'OS',
 'locus',
 'prot_name',
 'uid',
 'fetchid',
 'ch1',
 'ch2',
 'norm_ch1',
 'norm_ch2']
cols_simple += [
'uid_fetched',
'GN_fetched',
'signal',
'signal_loc',
'tm_span']
#
# FILL IN SOME COLUMNS ...
raw_fetch['enzyme']      = 'T' # all SILAC are trypsin ...
# raw_fetch['enzyme']      = raw_fetch['Biological sample category'].apply(ms.get_enzyme)
raw_fetch['uid_fetched'] = raw_fetch['fetchid'].apply(lambda fidx: gbrecs[str(int(fidx))].id if pd.notnull(fidx) else None)
raw_fetch['GN_fetched']  = raw_fetch['fetchid'].apply( ms.get_genename )
#
raw_fetch['signal']      = raw_fetch['fetchid'].apply( ms.get_signal )
raw_fetch['signal_loc']  = raw_fetch['fetchid'].apply( ms.get_signal_loc )
raw_fetch['tm_span']     = raw_fetch['fetchid'].apply( ms.get_tm_span )
#
# NOW SIMPLIFY THE DATAFRAME TO CONSIDER ONLY IMPORTNAT COLUMNS AND REMOVE DUPLICATES ...
raw_simple = raw_fetch[cols_simple].drop_duplicates()
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
raw_simple['crit_GN'] = (raw_simple['GN']==raw_simple['GN_fetched'])
# uid splitter - or . AMD uid major taker ...
uid_split = lambda uid: tuple(re.split('[-\.]',uid)) if pd.notnull(uid) else None
uid_major = lambda uid: re.split('[-\.]',uid)[0] if pd.notnull(uid) else None
raw_simple['crit_uid_full'] = (raw_simple['uid'].apply(uid_split)==raw_simple['uid_fetched'].apply(uid_split))
raw_simple['crit_uid_maj'] = (raw_simple['uid'].apply(uid_major)==raw_simple['uid_fetched'].apply(uid_major))
# pept_isin = lambda fidx: gbrecs[str(int(fidx))].id if pd.notnull(fidx) else None)


raw_simple['crit_pept_in'] = raw_simple[['pept','fetchid']].apply(ms.pept_isin,axis=1)
# extract pept_info fetched first ...
###############################################################
###############################################################
raw_simple = raw_simple.merge( raw_simple[['pept','fetchid']].apply(ms.pept_info,axis=1),left_index=True,right_index=True )
#
raw_simple['crit_start']   = False #raw_simple['pept_start'] == raw_simple['start_fetched']
raw_simple['crit_stop']    = False #raw_simple['pept_stop'] == raw_simple['stop_fetched']
raw_simple['crit_prev_aa'] = False #raw_simple['prev_aa'] == raw_simple['prev_aa_fetched']
raw_simple['crit_next_aa'] = False #raw_simple['next_aa'] == raw_simple['next_aa_fetched']
#
#
#
#
crit_cols = [cn for cn in raw_simple.columns if 'crit_' in cn]
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
raw_simple['SCORE'] = raw_simple[crit_cols].fillna(False).mul(crit_weight).sum(axis=1)
#
# ###################################################################
# #  EVALUATE THESE CRITERIA AND GET A COLUMN WITH SUM(AXIS=1)...   #
# ###################################################################
# # THEN DECIDE HOW MANY GENEBANK PROTEIN RECORDS QUALIFY CRITERIA,
# # IF IT IS JUST 1, THEN PROCEED WITH THE ONE, ELSE PRINT ALL THE INFO FOR
# # FURTHER INVESTIGATION ...
cols = [
 'spec_name',
 'exp_name',
 'pept',
 'fetchid',
 'GN',
 'GN_fetched',
 # 'prev_aa',
 # 'next_aa',
 'prev_aa_fetched',
 'next_aa_fetched',
 # 'pept_start',
 # 'pept_stop',
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
 'ch1',
 'ch2',
 'norm_ch1',
 'norm_ch2']


cols_short = [
 'spec_name',
 'exp_name',
 'pept',
 'fetchid',
 'enzyme',
 'GN',
 'GN_fetched',
 # 'prev_aa',
 # 'next_aa',
 'prev_aa_fetched',
 'next_aa_fetched',
 # 'pept_start',
 # 'pept_stop',
 'start_fetched',
 'stop_fetched',
 'locus',
 'prot_name',
 'uid',
 'uid_fetched',
 'signal',
 'signal_loc',
 'tm_span',
 'ch1',
 'ch2',
 'norm_ch1',
 'norm_ch2']


# output stuff ...
raw_simple_sorted = raw_simple.sort(columns=['pept','SCORE'],inplace=False)[cols+['SCORE']]
if full_sorted_output_with_criteria:
    raw_simple_sorted.to_csv('FULL_SORTED_pepts_and_scores.csv',index=False)

# group data by peptides ...
raw_grouped = raw_simple_sorted.groupby(by='pept')

# let's choose a single protein per peptide:
# idea: of max SCORE in a peptide-group is above Threshold AND is the only one in a group, THEN
# we assign that protein to peptide.
# OTHERWISE, (multiple max_SCORE proteins, no proteins that qualify criteria by Threshold) 
crit_threshold = 100
print "\nAttention!!! crit_threshold is redefined manually to 100 ...\n"
num_qualify_prots = lambda pep_grp: (pep_grp==pep_grp.max()).sum() if pep_grp.max()>=crit_threshold else 0
idx_qualify_prots = lambda pep_grp:  pep_grp.idxmax()  if num_qualify_prots(pep_grp)==1 else None
idx_ambig_prots   = lambda pep_grp:  pep_grp.idxmax()  if num_qualify_prots(pep_grp)>=2 else None
idx_bad_prots     = lambda pep_grp:  pep_grp.idxmax()  if num_qualify_prots(pep_grp)==0 else None

# now we have to retrieve those peptide-protein combinations that didn't work out
# (ambiguous or lacking qualified protein)
print
print "Pept-protein map relies on the uniqness of the max_SCORE value ..."
print "Here are max_SCORE occurences in the sample:"
print str(raw_grouped['SCORE'].apply(num_qualify_prots).value_counts())
print "Those with number of occurences != 1  will go to BAD_PEPTS file ..."
print

# indexes of qualified pept-protein pairs, including peptides with no certail match as 'None' ...
qual_prot_idxs  = raw_grouped['SCORE'].apply(idx_qualify_prots)
ambig_prot_idxs = raw_grouped['SCORE'].apply(idx_ambig_prots)
bad_prot_idxs   = raw_grouped['SCORE'].apply(idx_bad_prots)


# separate BAD from qualified pept-protein pairs ...
pept_prot_idxs = qual_prot_idxs[qual_prot_idxs.notnull()].map(int)
AMBIG_PEPTS = ambig_prot_idxs[ambig_prot_idxs.notnull()].index
BAD_PEPTS = bad_prot_idxs[bad_prot_idxs.notnull()].index


# GOOD ONE ARE OUT FIRST ...
# qualified pept-protein pairs goes straight to output for further use ...
pept_prot_map = raw_simple_sorted.loc[pept_prot_idxs][cols_short + ['SCORE','crit_pept_in']]
pept_prot_map.to_csv(os.path.join(common_path,'pept_prot_map.csv'),index=False)
# common_path


# FINISH IT UP WITH BAD ONES ...
###############################################################################
# modify BAD pepts, so that if a peptide is 'explained' by multiple proteins/fetchids
# and there both with pept in and not in, show those with in - only!!!
###############################################################################
filter_ambig = lambda df: df[df['SCORE']>=crit_threshold].reset_index(drop=True)
bad_pept_prot_output = pd.concat(
                    [raw_grouped.get_group(peptide) for peptide in BAD_PEPTS]+
                    [filter_ambig(raw_grouped.get_group(peptide)) for peptide in AMBIG_PEPTS]
                    ).reset_index(drop=True)[cols+['SCORE',]]
bad_pept_prot_output.to_csv(os.path.join(common_path,'BAD_pept_prot.csv'),index=False)


# STAGE 2 APPEARS TO BE WORKING PROPERLY ...











