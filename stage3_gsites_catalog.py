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
#
#
#
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
    pep_map_fname = os.path.join( args.prefix, args.pept_with_fetch )
    spec_fname    = os.path.join( args.prefix, args.spec_summary )
    gb_fname      = os.path.join( args.prefix, args.genbank )
else:
    pep_map_fname = args.pept_with_fetch
    spec_fname    = args.spec_summary
    gb_fname      = args.genbank
# get the common path for later use ...
common_path = os.path.commonprefix([pep_map_fname,spec_fname,gb_fname])
common_path = os.path.dirname(common_path)
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
pep_info = pd.read_csv(pep_map_fname,sep='\t')
spec_info = pd.read_csv(spec_fname,sep='\t')
#
# fasta = SeqIO.to_dict(SeqIO.parse(fasta_fname,"fasta"),key_function=lambda _: _.id.split('|')[1])
# 1-BASED NOTATION FOR PROTEINS INDEXING ENFORCED ...
# pep_df = pd.read_csv(uniq_pept_fname)


# connection between peptide info and spectrum info to be established ...
################################
#
# looking for Deamidated sites only ...
deamid = re.compile('[D,d]eamidat')
spec_info['pept'] = spec_info['Peptide sequence'].str.upper()
#

##########################################################################
#  MADE IT THUS FAR ...
#########################################################################
# CONTINUE LATER ON ............................
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#


# parse_pname = lambda pn: pd.Series(ms.parse_prot_name(pn))
# spec_info = spec_info.merge(spec_info['Protein name'].apply(parse_pname), left_index=True, right_index=True)
# pep_info = pep_info.merge(pep_info['Protein name'].apply(parse_pname), left_index=True, right_index=True)
#
def extract_deamids(mod_str):
    mod_list = [ mod.strip() for mod in mod_str.split(',') ]
    return pd.Series( ms.parse_spectrum_modifications(mod) for mod in mod_list if bool(deamid.search(mod)) )
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
##########################################################################
# unroll that spec table to have 1 deamid per row ...
spec_info_unrolled = unroll_by_mfunc(spec_info,'Variable modifications identified by spectrum',extract_deamids,'deamid_info')





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








# enumerating all tractable peptides from pep_df ...
glyco_mod = []
for uniq_pept, pept_pos, prot_sr in pep_df[['pept','peptide_start','prot_seqrec']].itertuples(index=False):
    prot_seq = str(prot_sr)
    pept_spectrum = spec_info[ spec_info['pept']==uniq_pept ]
    # let's check all present modifications in the flatten-out list of lists ...
    modifs = [ mod for cmod in pept_spectrum['Variable modifications identified by spectrum'] for mod in cmod.strip().split(',') ]
    # and get those that are unique ...
    modifs = np.unique(modifs)
    # now extract type,position and value for each of them ...
    modifs = [ ms.parse_spectrum_modifications(mod) for mod in modifs if bool(deamid.search(mod)) ]
    # now extrating meaningfull glycosilation sites ...
    glyco_sites = []
    glyco_start = []
    # looks like the inner loop here is the only place where we do need 0-based indexing switching ...
    for type_aa,gpos_pept,value in modifs:
        if (type_aa in ['n','N']) and (np.abs(value-3)<0.1):
            # 'pept_pos' - is 1-based absolute poisition of the peptide in the protein ...
            # 'gpos_pept' - is 1-based relative position of gsite_start_N in the peptide ...
            gsite_start = pept_pos + gpos_pept-1 # 1-based coordinate ...
            gsite_stop  = pept_pos + gpos_pept-1 + 3-1 # 1-based coordinate ...
            glyco_start.append(gsite_start)
            # Due to slicing rules, we need [start-1:stop], no have position 'stop' included ...
            glyco_sites.append(prot_seq[gsite_start-1:gsite_stop])
    ############################################################
    # gstart must be 1-based for output ...
    glyco_mod_str = ';'.join([ gsite+("(%d)"%gstart) for gsite,gstart in zip(glyco_sites,glyco_start)])
    glyco_mod.append(glyco_mod_str)
############################
#

#
############################
pep_df['gsites'] = glyco_mod
############################
print
print " We need to unroll glyco-sites and expand the DataFrame ..."
# we need to unroll those rows whith more than 1 items of gsites ...
pep_df_cols = ['all_uids', 'pept', 'peptide_start', 'prot_seqrec', 'protlen', 'uid_max', 'gsites','aa_before','aa_after',"prot_name", "uniq_pept_count", "pept_probab"]
pep_df_gsite_extracted = []
for all_uids,pept,pept_start,prot_sr,protlen,uid_max,gsites,aa_bef,aa_aft,prot_name,upept_count,pept_probab in pep_df[pep_df_cols].itertuples(index=False):
    # unrolling ...
    for gsite in gsites.split(';'):
        pep_df_gsite_extracted.append((all_uids,pept,pept_start,prot_sr,protlen,uid_max,gsite,aa_bef,aa_aft,prot_name,upept_count,pept_probab))
########################
pep_df_ext = pd.DataFrame(pep_df_gsite_extracted,columns = pep_df_cols)
print "mission accomplished ..."
# so far we didn't use gsite index for adressing, so it was 1-based safely ...

print
print "Now we'd need to collapse those identical glyco sites that orginiate from overlaing peptides ..."
gsite_uid_combined = pep_df_ext['gsites'] +':'+ pep_df_ext['all_uids']
gsites_uniq = gsite_uid_combined.unique()
print " There are %d unique (glyco-site, protein id) pairs that includes empty sites, however."%gsites_uniq.size
# now collapsing on those uniq sites ...
final_dataframe = []
for site_uniq in gsites_uniq:
    # indexes of a given unque gsite ...
    sites_index = (gsite_uid_combined == site_uniq)
    pep_site = pep_df_ext[sites_index]
    #
    pepts = ';'.join( '('+pep_site['aa_before']+')'+pep_site['pept']+'('+pep_site['aa_after']+')' )
    all_uids, = pep_site['all_uids'].unique()
    # peptide_start must be 1-based, we won't use it as index anymore ...
    pepts_start = ';'.join( str(pos_value) for pos_value in pep_site['peptide_start'] )
    prot_seq = str(pep_site['prot_seqrec'].unique()[0])
    prot_len, = pep_site['protlen'].unique()
    uid_max, = pep_site['uid_max'].unique()
    #
    prot_name, = pep_site["prot_name"].unique()
    upept_count = pep_site["uniq_pept_count"].unique()[0] # temporary solution, due to inconsistencies ...
    pept_probab = ';'.join(str(_) for _ in pep_site["pept_probab"])
    #
    # make sure there is just a single gsite ...
    pep_site_gsite_uniq, = pep_site['gsites'].unique()
    # appending to final dataframe ...
    final_dataframe.append((pep_site_gsite_uniq,pepts,pepts_start,all_uids,prot_seq,prot_len,uid_max, prot_name, upept_count, pept_probab))
# turn to dataframe ...
final_dataframe = pd.DataFrame(final_dataframe,columns=['gsites', 'pept', 'peptide_start','all_uids', 'prot_seq', 'protlen', 'uid_max','prot_name', 'uniq_pept_count', 'pept_probab'])



#
print 
print "Now we'd need to add theoretical glycosilation sites as a separate column ..."
g_site = re.compile(r'(?=(N[ACDEFGHIKLMNQRSTVWY][TS]))')


def get_theor_sites_number(prot_seq):
    # find all sites ...
    all_sites = [(site.start(),site.groups()[0]) for site in g_site.finditer(prot_seq)]
    N_sites = len(all_sites)
    return N_sites

def get_theor_sites(prot_seq):
    # BEWARE ZERO-BASED INDEXING TO 1-BASED INDEXING TRANSLATION ...
    # find all sites ...
    all_sites = [(site.start(),site.groups()[0]) for site in g_site.finditer(prot_seq)]
    # N_sites = len(all_sites)
    return ';'.join( (gsite_seq+'(%d)'%(pos+1)) for pos,gsite_seq in all_sites) # pos+1 - indexing transition ...


##################################################################################################################
print " 'gsites_predicted' column uses 1-based numbering. Enforced and checked."
final_dataframe['gsites_predicted'] = final_dataframe['prot_seq'].apply(get_theor_sites)
final_dataframe['gsites_predicted_number'] = final_dataframe['prot_seq'].apply(get_theor_sites_number)
# gsite_start must be 1-based ...
final_dataframe['gsite_start'] = final_dataframe['gsites'].apply( lambda x: int(x[3:].strip('()')) if x else None )
final_dataframe['gsites_AA1_N'] = final_dataframe['gsites'].apply( lambda x: x[0] if x else None )
final_dataframe['gsites_AA2_XbutP'] = final_dataframe['gsites'].apply( lambda x: x[1] if x else None )
final_dataframe['gsites_AA3_TS'] = final_dataframe['gsites'].apply( lambda x: x[2] if x else None )
# final_dataframe['gsites']NFT(1098)

# #
# # let's change within a cell separators to ';' instead of ',', which is used for column separation ...
# for col,dtype in final_dataframe.dtypes.iteritems():
#     if dtype=='object':
#         final_dataframe[col].str.replace(';',' ').str.replace(',',';')


final_dataframe.to_csv(out_fname,index=False)


# 1    NIS(126)                                     (R)ALSNISLR(F)
# 2    NTT(794)  (R)CVYEALCNTTSECPPPVITR(Q),(E)ALCNTTSECPPPVITR...
# 3   NDT(1048)  (R)EAESLQPMTVVGTDYVFHNDTK(V),(E)SLQPMTVVGTDYVF...
# 4    NET(732)                             (K)LSHDANETLPLHLYVK(Y)
# 5    NMS(527)                              (K)SCVAVTSAQPQNMSR(A)
# 6   NVT(1001)                               (R)SINVTGQGFSLIQR(A)
# 7   NFT(1098)  (R)TEAGAFEYVPDPTFENFTGGVK(Q),(R)TEAGAFEYVPDPTF...
# 8   NLT(1067)                  (K)VVFLSPAVPEEPEAYNLTVLIEMDGHR(L)











