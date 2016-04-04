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

import warnings
from Bio import BiopythonWarning, BiopythonParserWarning
#
#
#
########################################################################################################################
def fix_str(input_str,known_errors,known_fixes):
    """function that would be replacing known error strings in the input by the known fix replacement"""
    fixed_input = str(input_str)
    for err,fix in zip(known_errors,known_fixes):
        fixed_input = fixed_input.replace(err,fix)
    return fixed_input
########################################################################################################################
def fix_genbank(error_msg,genbank_fname):
    error_pattern = "(\d+\^\d+)"
    known_errors = re.findall(error_pattern,str(error_msg))
    known_fixes  = [err.replace('^','..') for err in known_errors]
    #
    if known_errors:
        with open(genbank_fname,'r') as fp:
            file_string_io = [ fix_str(line,known_errors,known_fixes) for line in fp ]
        file_string_io = StringIO(''.join(file_string_io))
        print "Error if fixed locally, proceeding ..."
        return file_string_io
    else:
        print "Error in the genbank could not be resolved. Termination"
        sys.exit(1)
########################################################################################################################
g_site = re.compile(r'(?=(N[ACDEFGHIKLMNQRSTVWY][TS]))')
#########################################################################################
def get_theor_sites_number(prot_seq):
    # find all sites ...
    all_sites = [(site.start(),site.groups()[0]) for site in g_site.finditer(prot_seq)]
    N_sites = len(all_sites)
    return N_sites
#########################################################################################
def get_theor_sites(prot_seq):
    # BEWARE ZERO-BASED INDEXING TO 1-BASED INDEXING TRANSLATION ...
    # find all sites ...
    all_sites = [(site.start(),site.groups()[0]) for site in g_site.finditer(prot_seq)]
    # N_sites = len(all_sites)
    return ';'.join( (gsite_seq+'(%d)'%(pos+1)) for pos,gsite_seq in all_sites) # pos+1 - indexing transition ...
##################################################################################################################
#
##################################################################################################################
# looking for Deamidated sites only ...
deamid = re.compile('[D,d]eamidat')
def extract_deamids(mod_str):
    mod_list = [ mod.strip() for mod in mod_str.split(',') ]
    # return pd.Series( ms.parse_spectrum_modifications(mod) for mod in mod_list if bool(deamid.search(mod)) )
    # collecting results to return ...
    to_return = []
    for mod in mod_list:
        if bool(deamid.search(mod)):
            type_aa, gpos_pept, value = ms.parse_spectrum_modifications(mod)
            # let's filter the aspartic ones with the value 3 right away ...
            if (type_aa in ['n','N']) and (np.abs(value-3)<0.01):
                to_return.append( (type_aa, gpos_pept, value) )
    return pd.Series( to_return )
##################################################################################################################
# unroll/exapnd spec table to account for multiple deamid-sites/gsites per peptide ...
def unroll_by_mfunc(df,col_name,mfunc,unrolled_colname='unrolled'):
    # get the column with arguments of the 'mfunc' ...
    thecolumn = df[col_name]
    # it appears as multi-column thing right after mfunc application ...
    multicol_mfunc_result = thecolumn.apply(mfunc)
    # stacked - multiindex is in use, level=1 index is the column name from 'multicol_mfunc_result' ...
    unrolled_mindex = multicol_mfunc_result.stack()
    # index is no longer uniq after the following operation ... we dropped inner indexing part.
    # IMPORTANT: indexes correspond to those from the original df (use it for merging later) ...
    unrolled_origindex = pd.DataFrame(unrolled_mindex.reset_index(level=1,drop=True),columns=[unrolled_colname,])
    # merge unrolled_origindex (a single column with ambiguous index) with the original df ...
    # 'unrolled_origindex' must be DataFrame to merge: Series are not mergeable for some reason ...
    # the whole df is to be unrolled after the following operation.
    unrolled_df = df.merge(unrolled_origindex,left_index=True,right_index=True)
    #
    return unrolled_df.reset_index(drop=True)
##################################################################################################################
def deamid_to_gsite(deamid_mod, pept_pos, prot_seq):
    type_aa, gpos_pept, value = deamid_mod
    gpos_pept = int(gpos_pept)
    pept_pos = int(pept_pos)
    assert type_aa in ['n','N']
    assert np.abs(value-3)<0.01
    # 'pept_pos' - is 1-based absolute poisition of the peptide in the protein ...
    # 'gpos_pept' - is 1-based relative position of gsite_start_N in the peptide ...
    gsite_start = pept_pos + gpos_pept-1 # 1-based coordinate ...
    gsite_stop  = pept_pos + gpos_pept-1 + 3-1 # 1-based coordinate ...
    # Due to slicing rules, we need [start-1:stop], no have position 'stop' included ...
    gsite_seq = prot_seq[gsite_start-1:gsite_stop]
    ############################################################
    # gstart must be 1-based for output ...
    return {'gsite':"%s(%d)"%(gsite_seq,gsite_start), 'gsite_seq':gsite_seq, 'gstart':gsite_start}
#########################################################################################
#
#
# HOW TO LAUNCH THIS THING ...
# %run stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/011216\ glycocapture\ 90-90 -m pept_prot_map.csv -g pulled_proteins.gb -s specs.xls
#
# do some arguments parsing to make the script looks civilized ...
parser = argparse.ArgumentParser()
parser.add_argument("-m","--pept_prot_map",
    help="specify file name of peptide summary with unique of fetchids of matching proteins (with/without path)",required=True)
parser.add_argument("-g","--genbank",
    help="specify file name of genbank records with pulled proteins (with/without path)",required=True)
parser.add_argument("-s","--spec_summary", help="speicfy spectrum file name (with/without path)",required=True)
parser.add_argument("--prefix", help="specify common part of the path for peptide and spectrum files")
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
#
#
############################
# READING SEQUENCES FROM THE FILE ...
print "Reading %s with genebank records from the NCBI fetch ..."%gb_fname
with warnings.catch_warnings():
    # e.g. BiopythonParserWarning: Dropping bond qualifier in feature location
    #
    warnings.simplefilter("ignore", BiopythonParserWarning)
    #
    gb_recs_iter = SeqIO.parse(gb_fname,'gb')
    try:
        gbrecs = SeqIO.to_dict( gb_recs_iter, key_function=lambda rec: rec.annotations['gi'] )
    except ValueError, er_msg:
        print "Catched ValueError: %s"%str(er_msg)
        # #
        # Invalid between location '1^593'
        # Try to fix that thing, by replacing '^' with '..'
        file_string_io = fix_genbank(er_msg,gb_fname)
        gb_recs_iter = SeqIO.parse(file_string_io,'gb')
        gbrecs = SeqIO.to_dict( gb_recs_iter, key_function=lambda rec: rec.annotations['gi'] )
#
#
#
#
pep_info = pd.read_csv(pep_map_fname,sep=',')
spec_info = pd.read_csv(spec_fname,sep='\t')
# fix their peptide sequence thing right away ...
spec_info['pept'] = spec_info['Peptide sequence'].str.upper()
pep_info['fetchid'] = pep_info['fetchid'].apply(int)
# fasta = SeqIO.to_dict(SeqIO.parse(fasta_fname,"fasta"),key_function=lambda _: _.id.split('|')[1])
# 1-BASED NOTATION FOR PROTEINS INDEXING ENFORCED ...
# pep_df = pd.read_csv(uniq_pept_fname)


# connection between peptide info and spectrum info to be established ...
##########################################################################
# unroll that spec table to have 1 deamid per row ...
spec_info_unrolled = unroll_by_mfunc(spec_info,'Variable modifications identified by spectrum',extract_deamids,'deamid_info')
spec_info_unrolled['prot_ident_probab'] = spec_info_unrolled['Protein identification probability'].str.strip('%').apply(float)
spec_info_unrolled['pept_ident_probab'] = spec_info_unrolled['Peptide identification probability'].str.strip('%').apply(float)
##########################################################

# so far the following merge seems to be 100% sufficient for the desired final output ...
# we could add on extra features if needed ...
spec_n_pep = spec_info_unrolled[['pept',
                    'deamid_info',
                    'prot_ident_probab',
                    'pept_ident_probab',
                    # 'Protein name'
                    ]].merge(pep_info,how='right',on='pept')

# Now, extract those gsites ...
dg_func = lambda x: pd.Series( deamid_to_gsite(x['deamid_info'], x['start_fetched'], str(gbrecs[str(int(x['fetchid']))].seq)) )
# and add them back to the main table ...
gs_res = spec_n_pep[['deamid_info','start_fetched','fetchid']].apply( dg_func, axis=1 )
spec_n_pep = spec_n_pep.merge(gs_res,left_index=True,right_index=True)





print 
print "Now we'd need to add theoretical glycosilation sites as a separate column ..."
# this analysis must be done, once for each 'fetchid', and then merged back to the main table ...

get_theor_sites_fid = lambda fid: get_theor_sites(str(gbrecs[str(fid)].seq))
get_theor_sites_number_fid = lambda fid: get_theor_sites_number(str(gbrecs[str(fid)].seq))
theor_sites_info = lambda fid: pd.Series(
                                    {'fetchid':fid,
                                     'gsites_predicted':get_theor_sites_fid(fid),
                                     'gsites_predicted_number':get_theor_sites_number_fid(fid)}
                                        )
###################################################
predicted_gsite_info = spec_n_pep['fetchid'].drop_duplicates().apply(theor_sites_info)
# add back to the main table ...
spec_n_pep = spec_n_pep.merge(predicted_gsite_info,on='fetchid',how='right')

print "done ..."
print "numbering appears to be 1-based and overall correct!"
print
# print " 'gsites_predicted' column uses 1-based numbering. Enforced and checked."

# SOME FINAL PREPARATIONS TO COMPLY WITH THE REQUESTED FORMAT ...

# extract gsite AAs as separate columns ...
spec_n_pep['gsite_AA1'] = spec_n_pep['gsite_seq'].str.get(0)
spec_n_pep['gsite_AA2'] = spec_n_pep['gsite_seq'].str.get(1)
spec_n_pep['gsite_AA3'] = spec_n_pep['gsite_seq'].str.get(2)



requested_cols = ['gsite',
'pept',
'enzyme',
'start_fetched',
'prot_name',
'fetchid',
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
'tm_span',]

###################################################
# TO BE CONTINUED ...
THE_MOST_FINAL_DF = spec_n_pep[requested_cols]


# DESIREd COLUMNS ...
# ############################################
# #  columns that needs to be delivered ...  #
# ############################################
# # A gsites, 1 per line
# # B pept, 1 per line
# # B1 enzyme, G or T, derive from 'Biological sample category', like this: {'TrypsinSample1':'T','GluC_Sample2':'G'}
# # C peptide_start, 1 per line accordingly
# # D all_uids, REPLACE WITH col:H
# # E prot_seq, try to get those from NCBI, not from UniProt ...
# # F protein, ??? sequence, name or what???
# # G uid_max, UID for major form instead or something like that ...
# # H prot_name, parsed out human-readable name from 'Protein name'
# # H1 gene_name, parsed out GN=xxx from 'Protein name'
# # I uniq_peptide_count, discrad that column ...
# # J pept_probability, output number not the string - this would be the criteria 
# # K gsites_predicted, OK
# # L gsites_predicted_number, OK
# # M gsite_start, beware of 0 or 1 type of indexing ...
# # N,O,P - gsites AAs in separate columns
# # M1, NOP combined, gsite sequence basically!
# # Q signal, from GeneBank record on the protein, simply Y,N on whether there is a 'Signal' in gb.
# # R signal_location, location of the signal from Q
# # S tm_span, Y,N just for the fact of having TM span as a protein feature.


# THINGS TO BE DONE:

# 1) MERGE PEP_INFO WITH SPEC_INFO, SO THAT PEPT-PROT RELATIONSHIP WOULD BE SET UP ...
#     probably getting read of many-many columns in the spec_info along the way ...

# 2) MODIFY AND TEST THE 'deamid_to_gsite' FUNCTION (PAY ATTENTION TO 1-BASED AND 0-BASED NUMBERING OF AAs)

# 3) USE 'deamid_to_gsite' TO FINALLY EXTRACT GSITES ...

# 4) ANALYSE(?) GSITES: GROUPBY UNIQUE GSITES (GSITE[seq],GSITE_START,PROT_IDENTIFICATOR) TO SEE HOW MANY PEPTS YIELD THE SAME GSITES ...

# 5) SELECT A SINGLE PEPT PER GSITE??????? USING REAID'S CRITERIA, OR LEAVE ALL GSITE-PEPTS PAIRS ???????????
                     #######################################
# 6) COMPLY WITH THE #columns that needs to be delivered...# FOR THE VERY VERY FINAL OUTPUT ...
                     #######################################









