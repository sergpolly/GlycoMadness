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
from StringIO import StringIO
#
import argparse
#
#
# HOW TO LAUNCH THIS THING ...
# %run stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/011216\ glycocapture\ 90-90 -m pept_prot_map.csv -g pulled_proteins.gb -s specs.xls
#
# do some arguments parsing to make the script looks civilized ...
parser = argparse.ArgumentParser()
parser.add_argument("-e","--exp_num",
    help="specify number of an experiment",required=True)
parser.add_argument("-m","--pept_prot_map",
    help="specify file name of peptide summary with unique of fetchids of matching proteins (with/without path)",required=True)
parser.add_argument("-g","--genbank",
    help="specify file name of genbank records with pulled proteins (with/without path)",required=True)
parser.add_argument("-s","--spec_summary", help="speicfy spectrum file name (with/without path)",required=True)
parser.add_argument("--prefix", help="specify common part of the path for peptide and spectrum files")
parser.add_argument("--separator", help="speicfy separator type in the input data",default='tab')
args = parser.parse_args()
#
###############################################
if args.prefix is not None:
    pep_map_fname = os.path.join( args.prefix, args.pept_prot_map )
    spec_fname    = os.path.join( args.prefix, args.spec_summary )
    gb_fname      = os.path.join( args.prefix, args.genbank )
else:
    pep_map_fname = args.pept_prot_map
    spec_fname    = args.spec_summary
    gb_fname      = args.genbank
# get the common path for later use ...
common_path = os.path.commonprefix([pep_map_fname,spec_fname,gb_fname])
common_path = os.path.dirname(common_path)
#
# Reading genbank mindfully next ...
gbrecs = ms.genebank_fix_n_read(gb_fname,key_func_type='id')
######################################
# assign some module internal stuff ...
ms.gbrecs = gbrecs
#
# separator type choice is needed only for the ORIGINAL input files ...
if args.separator == "tab":
    separator = '\t'
elif args.separator == "comma":
    separator = ','
else:
    separator = '\t'
#################
pep_info = pd.read_csv(pep_map_fname)
spec_info = pd.read_csv(spec_fname,sep=separator)
# fix their peptide sequence thing right away ...
spec_info['pept'] = spec_info['Peptide sequence'].str.upper()
pep_info['fetchid'] = pep_info['fetchid'].apply(int)
# this is an UGLY fix that we'd have to implement here just to save everything ...
if args.exp_num=='1':
    pep_info['enzyme'] = 'T'
#
# fasta = SeqIO.to_dict(SeqIO.parse(fasta_fname,"fasta"),key_function=lambda _: _.id.split('|')[1])
# 1-BASED NOTATION FOR PROTEINS INDEXING ENFORCED ...
# pep_df = pd.read_csv(uniq_pept_fname)

# connection between peptide info and spectrum info to be established ...
##########################################################################
# unroll that spec table to have 1 deamid per row ...
spec_info_unrolled = ms.unroll_by_mfunc(spec_info,'Variable modifications identified by spectrum',ms.extract_deamids,'deamid_info')
if spec_info_unrolled['Protein identification probability'].dtype != 'float':
    spec_info_unrolled['prot_ident_probab'] = spec_info_unrolled['Protein identification probability'].str.strip('%').apply(float)
else:
    spec_info_unrolled['prot_ident_probab'] = spec_info_unrolled['Protein identification probability']
####################################################################################################################
if spec_info_unrolled['Peptide identification probability'].dtype != 'float':
    spec_info_unrolled['pept_ident_probab'] = spec_info_unrolled['Peptide identification probability'].str.strip('%').apply(float)
else:
    spec_info_unrolled['pept_ident_probab'] = spec_info_unrolled['Peptide identification probability']
####################################################################################################################

# so far the following merge seems to be 100% sufficient for the desired final output ...
# we could add on extra features if needed ...
spec_n_pep = spec_info_unrolled[['pept',
                    'deamid_info',
                    'prot_ident_probab',
                    'pept_ident_probab']].merge(pep_info,how='right',on='pept',suffixes=('','_x'))

# Now, extract those gsites ...
dg_func = lambda x: pd.Series( ms.deamid_to_gsite(x['deamid_info'], x['start_fetched'], str(gbrecs[x['fetchacc']].seq)) )
# and add them back to the main table ...
gs_res = spec_n_pep[['deamid_info','start_fetched','fetchacc']].apply( dg_func, axis=1 )
spec_n_pep = spec_n_pep.merge(gs_res,left_index=True,right_index=True)





print 
print "Now we'd need to add theoretical glycosilation sites as a separate column ..."
print "full protein sequence and its length is added as well ..."
# this analysis must be done, once for each 'fetchid', and then merged back to the main table ...

get_theor_sites_fid = lambda facc: ms.get_theor_sites(str(gbrecs[str(facc)].seq))
get_theor_sites_number_fid = lambda facc: ms.get_theor_sites_number(str(gbrecs[str(facc)].seq))
theor_sites_info = lambda facc: pd.Series(
                                    {'fetchacc':facc,
                                     'gsites_predicted':get_theor_sites_fid(facc),
                                     'gsites_predicted_number':get_theor_sites_number_fid(facc),
                                     'prot_seq':str(gbrecs[str(facc)].seq),
                                     'prot_len':len(str(gbrecs[str(facc)].seq))} )
###################################################
predicted_gsite_info = spec_n_pep['fetchacc'].drop_duplicates().apply(theor_sites_info)
# add back to the main table ...
spec_n_pep = spec_n_pep.merge(predicted_gsite_info,on='fetchacc',how='right')

print "done ..."
print "numbering appears to be 1-based and overall correct!"
print
# print " 'gsites_predicted' column uses 1-based numbering. Enforced and checked."

# SOME FINAL PREPARATIONS TO COMPLY WITH THE REQUESTED FORMAT ...

# extract gsite AAs as separate columns ...
spec_n_pep['gsite_AA1'] = spec_n_pep['gsite_seq'].str.get(0)
spec_n_pep['gsite_AA2'] = spec_n_pep['gsite_seq'].str.get(1)
spec_n_pep['gsite_AA3'] = spec_n_pep['gsite_seq'].str.get(2)



# locus  protein_name  uid  Protein_ID_PERCENT    peptides    best_peptide    Peptide_probability protease    Expt_NUMBER  prev_aa next_aa pept_start  pept_stop   Location match  g_site  gsite_start gsites_AA1_N    gsites_AA2_XbutP    gsites_AA3_TS   Best Mascot Ion score   Best Mascot Identity score  Best Mascot Delta Ion score Prot_seq    signal  signal_loc  tm_span protein length


requested_cols = ['gsite',
'pept',
'enzyme',
'start_fetched',
'prot_name',
'fetchid',
'fetchacc',
'uid_fetched',
'GN_fetched',
'pept_ident_probab',
'gsites_predicted',
'gsites_predicted_number',
'gsite_seq',
'gstart',
'gsite_AA1',
'gsite_AA2',
'gsite_AA3',
'signal',
'signal_loc',
'tm_span']


requested_cols = ['locus',
'prot_name',
'uid_fetched',
'prot_ident_probab',
# 'peptides',  # THIS IS NEW STUFF ...
'pept',
'fetchid',
'fetchacc',
# 'best_pept', # THIS IS NEW STUFF ...
'pept_ident_probab', # BEWARE, pept ID % of the BEST PEPTIDE ...
'enzyme',
# 'experiment_num', # THIS IS NEW STUFF ...
###########################
'prev_aa',
'next_aa',
'prev_aa_fetched',
'next_aa_fetched',
'pept_start',
'pept_stop',
'start_fetched',
'stop_fetched',
###########################
'gsite_seq',
'gstart',
'gsite_AA1',
'gsite_AA2',
'gsite_AA3',
'Best Mascot Ion Score',
'Best Mascot Identity Score',
'Best Mascot Delta Ion Score',
'prot_seq', # PROTEIN SEQUENCE TO BE ADDED ...
'signal',
'signal_loc',
'tm_span',
'prot_len', # PROTEIN LENGTH TO BE ADDED ...
'SCORE',
'crit_pept_in']







# ###################################################
# # TO BE CONTINUED ...
THE_MOST_FINAL_DF = spec_n_pep[requested_cols].drop_duplicates().reset_index(drop=True)

# choose peptide with highest Pept_ident_probab 
# Let's collpase (gsite,pept,fetchid) using the highest pept_ident_probab ...
THE_MOST_FINAL_DF_max_prob = THE_MOST_FINAL_DF.loc[THE_MOST_FINAL_DF.groupby(['gsite_seq','gstart','pept','fetchid','fetchacc','enzyme'],sort=False)['pept_ident_probab'].idxmax() ].reset_index(drop=True)



##################################################
#   TO BE CONTINUED ....
#################################################

# BEST PEPTIDE CHOSING IS WORKING, FINISH ALL PEPTIDES AND OTHER STUFF ...
def choose_best_pept(df,dp_margin=1.0):
    # pandas 0.18.1 broke something here: "TypeError: Series.name must be a hashable type" ...
    # we're fixing it, introducing one more layer of control here ...
    assert 'pept' in df.columns
    assert 'pept_ident_probab' in df.columns
    # this should work on any df with the 2 specified columns ...
    if df.shape[0]>1:
        shortes_idx  = df['pept'].str.len().idxmin()
        max_prob_idx = df['pept_ident_probab'].idxmax()
        if df['pept_ident_probab'][shortes_idx] < (df['pept_ident_probab'][max_prob_idx]-dp_margin):
            # return df['pept'][max_prob_idx]
            # print "returning max_prob_idx",max_prob_idx
            return max_prob_idx
        else:
            # return df['pept'][shortes_idx]
            # print "returning shortes_idx",shortes_idx
            return shortes_idx
    else:
        # print "returning one thing", df['pept_ident_probab'].idxmax()
        return df['pept_ident_probab'].idxmax()


def unstack_pept(series):
    return ','.join(series)


# fff=THE_MOST_FINAL_DF_max_prob.groupby(['gsite_seq','gstart','locus','enzyme'],sort=False).size()
# fff[fff>1]


THE_MOST_FINAL_DF_max_prob.groupby(['gsite_seq', 'gstart', 'locus', 'enzyme'],sort=False)['pept'].apply(unstack_pept).reset_index(drop=True)




# choosing BEST PEPTIDE containing given gsite 
# according to the Reid's rule of best peptide ...
THE_MOST_FINAL_DF_uniq_pept = \
        THE_MOST_FINAL_DF_max_prob.loc[
                THE_MOST_FINAL_DF_max_prob.groupby(['gsite_seq',
                                                    'gstart',
                                                    'locus',
                                                    'enzyme'],sort=False).apply(choose_best_pept)
                                         ].reset_index(drop=True)

# given the best peptides are choosen, unstacking multiple peptides corresponding to the same gsite
# and forming a column 'peptides' ...
THE_MOST_FINAL_DF_uniq_pept['peptides'] = \
        THE_MOST_FINAL_DF_max_prob.groupby(['gsite_seq',
                                              'gstart',
                                              'locus',
                                              'enzyme'],sort=False)['pept'].apply(unstack_pept).reset_index(drop=True)


# rename pept to best_pept AND enzyme to protease ...
THE_MOST_FINAL_DF_uniq_pept = THE_MOST_FINAL_DF_uniq_pept.rename(columns={'pept':'best_pept','enzyme':'protease',})

# add experiment number, something new ...
THE_MOST_FINAL_DF_uniq_pept['exp_num'] = int(args.exp_num)
# location match instead of fetched/Scaffold-based resutls ...
THE_MOST_FINAL_DF_uniq_pept['loc_match'] = THE_MOST_FINAL_DF_uniq_pept['pept_start'] == THE_MOST_FINAL_DF_uniq_pept['start_fetched']
THE_MOST_FINAL_DF_uniq_pept['loc_match'] = THE_MOST_FINAL_DF_uniq_pept['loc_match'].map({True:'Y',False:'N'})


requested_cols = ['locus',
'prot_name',
'uid_fetched',
'prot_ident_probab',
'peptides',  # THIS IS NEW STUFF ...
# 'pept',
# 'fetchid',
'fetchacc',
'best_pept', # THIS IS NEW STUFF ...
'pept_ident_probab', # BEWARE, pept ID % of the BEST PEPTIDE ...
'protease',
'exp_num', # THIS IS NEW STUFF ...
###########################
'prev_aa',
'next_aa',
# 'prev_aa_fetched',
# 'next_aa_fetched',
'pept_start',
'pept_stop',
'loc_match',
# 'start_fetched',
# 'stop_fetched',
###########################
'gsite_seq',
'gstart',
'gsite_AA1',
'gsite_AA2',
'gsite_AA3',
'Best Mascot Ion Score',
'Best Mascot Identity Score',
'Best Mascot Delta Ion Score',
'prot_seq', # PROTEIN SEQUENCE TO BE ADDED ...
'signal',
'signal_loc',
'tm_span',
'prot_len', # PROTEIN LENGTH TO BE ADDED ...
'SCORE',
'crit_pept_in']




THE_MOST_FINAL_DF_uniq_pept[requested_cols].to_csv(os.path.join(common_path,'FINAL_gsite_anthology.csv'),index=False)



# THE_MOST_FINAL_DF_uniq.to_csv(os.path.join(common_path,'FINAL_gsite_anthology.csv'),index=False)



# # DESIREd COLUMNS ...
# # ############################################
# # #  columns that needs to be delivered ...  #
# # ############################################
# # # A gsites, 1 per line
# # # B pept, 1 per line
# # # B1 enzyme, G or T, derive from 'Biological sample category', like this: {'TrypsinSample1':'T','GluC_Sample2':'G'}
# # # C peptide_start, 1 per line accordingly
# # # D all_uids, REPLACE WITH col:H
# # # E prot_seq, try to get those from NCBI, not from UniProt ...
# # # F protein, ??? sequence, name or what???
# # # G uid_max, UID for major form instead or something like that ...
# # # H prot_name, parsed out human-readable name from 'Protein name'
# # # H1 gene_name, parsed out GN=xxx from 'Protein name'
# # # I uniq_peptide_count, discrad that column ...
# # # J pept_probability, output number not the string - this would be the criteria 
# # # K gsites_predicted, OK
# # # L gsites_predicted_number, OK
# # # M gsite_start, beware of 0 or 1 type of indexing ...
# # # N,O,P - gsites AAs in separate columns
# # # M1, NOP combined, gsite sequence basically!
# # # Q signal, from GeneBank record on the protein, simply Y,N on whether there is a 'Signal' in gb.
# # # R signal_location, location of the signal from Q
# # # S tm_span, Y,N just for the fact of having TM span as a protein feature.


# # THINGS TO BE DONE:

# # 1) MERGE PEP_INFO WITH SPEC_INFO, SO THAT PEPT-PROT RELATIONSHIP WOULD BE SET UP ...
# #     probably getting read of many-many columns in the spec_info along the way ...

# # 2) MODIFY AND TEST THE 'deamid_to_gsite' FUNCTION (PAY ATTENTION TO 1-BASED AND 0-BASED NUMBERING OF AAs)

# # 3) USE 'deamid_to_gsite' TO FINALLY EXTRACT GSITES ...

# # 4) ANALYSE(?) GSITES: GROUPBY UNIQUE GSITES (GSITE[seq],GSITE_START,PROT_IDENTIFICATOR) TO SEE HOW MANY PEPTS YIELD THE SAME GSITES ...

# # 5) SELECT A SINGLE PEPT PER GSITE??????? USING REAID'S CRITERIA, OR LEAVE ALL GSITE-PEPTS PAIRS ???????????
#                      #######################################
# # 6) COMPLY WITH THE #columns that needs to be delivered...# FOR THE VERY VERY FINAL OUTPUT ...
#                      #######################################









