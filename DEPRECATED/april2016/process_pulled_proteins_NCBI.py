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
Entrez.email = "sergey.venev@umassmed.edu"
#
#
# # for these pair: there are peptides present in spec, but missing from pep isn't it odd?
# # suggesting, that different cutoffs were used for Scaffold program to output pep and spec ...
# pep_fname = "../raw_data/New_files_to_analyze/original_input_before_dec25/peptides.xls"
# spec_fname = "../raw_data/New_files_to_analyze/original_input_before_dec25/specs.xls"
#
# pept from spec are all in pep and all Peptide sequences are in
pep_fname = "../raw_data/New_files_to_analyze/011616 glycocapture 90-80/peptides.xls"
spec_fname = "../raw_data/New_files_to_analyze/011616 glycocapture 90-80/specs.xls"
#
# files generated after the NCBI genebank records search/fetch ...
spec_info_with_fetch_fname = 'spec_info_with_fetch.csv'
pulled_gb_recs_fname = "pulled_proteins.gb"
#
#
pep_info = pd.read_csv(pep_fname,sep='\t')
spec_info = pd.read_csv(spec_fname,sep='\t')
#
# fasta = SeqIO.to_dict(SeqIO.parse(fasta_fname,"fasta"),key_function=lambda _: _.id.split('|')[1])
# 1-BASED NOTATION FOR PROTEINS INDEXING ENFORCED ...
# pep_df = pd.read_csv(uniq_pept_fname)
#
#
############################
# READING SEQUENCES FROM THE FILE ...
print "Reading %s with genebank records from the NCBI fetch ..."%pulled_gb_recs_fname
gb_recs_iter = SeqIO.parse(pulled_gb_recs_fname,'gb')
gbrecs = SeqIO.to_dict( gb_recs_iter, key_function=lambda rec: rec.annotations['gi'] )
############################
# # TOUGH CHOISE ON WHICH PROTEIN TO FINALLY ASSIGN TO THE GN ...
# # yet to be made ...
#
############################
# READING file containing GeneName and FetchID association ...
print "Reading %s with the updated spectrum that includes fetchid column ..."%spec_info_with_fetch_fname
spec_info_with_fetch = pd.read_csv(spec_info_with_fetch_fname,sep=',')
############################
# # TOUGH CHOISE ON WHICH PROTEIN TO FINALLY ASSIGN TO THE GN ...
# # yet to be made ...
#
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
#############################################################################################################
# """Simply looking into DSC2 and DSC3 case tell us that there have to be sophisticated algorithm to choose
# the best matching fetchid for each gene name. Simply taking the first answer wouldn't work the best."""
# # cccf = ['pept', 'GN', 'OS', 'locus', 'prot_name', 'uid', 'deamid_info', 'fetchid']
# # ccc = ['pept', 'GN', 'OS', 'locus', 'prot_name', 'uid', 'deamid_info']
# #########################
# # LOOKS LIKE UID IS A GREAT CANDIDATE TO DISCERN
# # BETWEEN RIGHT AND WRONG SEQUENCES RETURNED BY NCBI FOR A GENENAME REQUEST ...
# #################################
spec_info_with_fetch['uid_fetched'] = spec_info_with_fetch['fetchid'].apply(lambda fidx: gbrecs[str(int(fidx))].id if pd.notnull(fidx) else None)
spec_info_with_fetch['GN_fetched']  = spec_info_with_fetch['fetchid'].apply(lambda fidx: gbrecs[str(int(fidx))].features[1].qualifiers['gene'][0] if pd.notnull(fidx) else None)
# # cccf = ['pept', 'GN','GN_pull', 'locus', 'prot_name', 'uid','uid_pull', 'deamid_info', 'fetchid']
# # # TOUGH CHOISE ON WHICH PROTEIN TO FINALLY ASSIGN TO THE GN ...
# yet to be masde ...
#


######################################################################################################
yn_map = {True:'Y',False:'N'}
######################################################################################################
def get_tm_span(fidx):
    if pd.notnull(fidx):
        fidx_str = str(int(fidx))
        feats_descr = []
        for feat in gbrecs[fidx_str].features:
            quals = feat.qualifiers
            feats_descr.append( yn_map["Transmembrane" in ''.join(quals['region_name'])] if ('region_name' in quals) else None )
        return 'Y' if ('Y' in feats_descr) else 'N'
    else:
        return None
######################################################################################################
def get_signal(fidx):
    if pd.notnull(fidx):
        fidx_str = str(int(fidx))
        feats_descr = []
        for feat in gbrecs[fidx_str].features:
            quals = feat.qualifiers
            feats_descr.append( yn_map['Signal' in quals['region_name']] if ('region_name' in quals) else None )
        return 'Y' if ('Y' in feats_descr) else 'N'
    else:
        return None
######################################################################################################
# to be edited to trun into feature locator ...
def get_signal_loc(fidx):
    if pd.notnull(fidx):
        fidx_str = str(int(fidx))
        for feat in gbrecs[fidx_str].features:
            quals = feat.qualifiers
            if ('region_name' in quals):
                if 'Signal' in quals['region_name']:
                    # start,end = (feat.location.start.position, feat.location.end.position)
                    return "%d..%d"%(feat.location.start.position, feat.location.end.position)
        return None
    else:
        return None
######################################################################################################


spec_info_with_fetch['signal']  = spec_info_with_fetch['fetchid'].apply( get_signal )
spec_info_with_fetch['signal_loc']  = spec_info_with_fetch['fetchid'].apply( get_signal_loc )
spec_info_with_fetch['tm_span']  = spec_info_with_fetch['fetchid'].apply( get_tm_span )

# ['signal', 'signal_loc', 'tm_span']

# #
# ############################
# results = []
# for idx in record['IdList']:
#     print "fetching",idx
#     handle = Entrez.efetch(db='protein',id=idx,rettype='gb',retmode='text')
#     filelike = StringIO(handle.read())
#     seqrec = SeqIO.read(filelike,format='gb')
#     results.append( seqrec )
#########################################################################################
# # The features that Reid is asking about are accessible through:
# r0 = results[0]
# f4 = r0.features[4]
# f4.qualifiers['region_name']







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
spec_info_with_fetch['crit_GN'] = (spec_info_with_fetch['GN']==spec_info_with_fetch['GN_fetched'])
# uid splitter - or . AMD uid major taker ...
uid_split = lambda uid: tuple(re.split('[-\.]',uid)) if pd.notnull(uid) else None
uid_major = lambda uid: re.split('[-\.]',uid)[0] if pd.notnull(uid) else None
spec_info_with_fetch['crit_uid_full'] = (spec_info_with_fetch['uid'].apply(uid_split)==spec_info_with_fetch['uid_fetched'].apply(uid_split))
spec_info_with_fetch['crit_uid_maj'] = (spec_info_with_fetch['uid'].apply(uid_major)==spec_info_with_fetch['uid_fetched'].apply(uid_major))
# pept_isin = lambda fidx: gbrecs[str(int(fidx))].id if pd.notnull(fidx) else None)
def pept_isin(row):
    pept,fetchid = row
    if pd.notnull(fetchid):
        fidx = str(int(fetchid))
        prot_seq = gbrecs[fidx].seq
        return (pept in prot_seq)
    else:
        None
spec_info_with_fetch['crit_pept_in'] = spec_info_with_fetch[['pept','fetchid']].apply(pept_isin,axis=1)
# extract pept_info fetched first ...
###############################################################
def pept_info(row):
    pept,fetchid = row
    if pd.notnull(fetchid):
        fidx = str(int(fetchid))
        prot_seq = gbrecs[fidx].seq
        # find pept in prot:
        pept_found = prot_seq.find(pept)
        if pept_found > -1:
            # 1-based indexing right away ...
            start = prot_seq.find(pept) + 1
            stop  = prot_seq.find(pept) + len(pept)
            # because of 1-based indexing ...
            if stop >= len(prot_seq):
                print "peptide is at the end: ",pept,fidx
                next_aa = None
            else:
                next_aa = prot_seq[stop]
            ################################
            if start <= 1:
                print "peptide is at the start: ",pept,fidx
                prev_aa = None
            else:
                prev_aa = prot_seq[start-2]
            # return 4 columns ...
            return pd.Series( {'start_fetched': start,
                'stop_fetched': stop,
                'prev_aa_fetched': prev_aa,
                'next_aa_fetched': next_aa} )
        else:
            print "(!!!) peptide not found: ",pept,fidx
            None
    else:
        print "fetchid is None for pept",pept
        None
###############################################################
spec_info_with_fetch = spec_info_with_fetch.merge( spec_info_with_fetch[['pept','fetchid']].apply(pept_info,axis=1),left_index=True,right_index=True )
#
spec_info_with_fetch['crit_start'] = spec_info_with_fetch['Peptide start index'] == spec_info_with_fetch['start_fetched']
spec_info_with_fetch['crit_stop'] = spec_info_with_fetch['Peptide stop index'] == spec_info_with_fetch['stop_fetched']
spec_info_with_fetch['crit_prev_aa'] = spec_info_with_fetch['Previous amino acid'] == spec_info_with_fetch['prev_aa_fetched']
spec_info_with_fetch['crit_next_aa'] = spec_info_with_fetch['Next amino acid'] == spec_info_with_fetch['next_aa_fetched']
###################################################################
#  EVALUATE THESE CRITERIA AND GET A COLUMN WITH SUM(AXIS=1)...   #
###################################################################
# THEN DECIDE HOW MANY GENEBANK PROTEIN RECORDS QUALIFY CRITERIA,
# IF IT IS JUST 1, THEN PROCEED WITH THE ONE, ELSE PRINT ALL THE INFO FOR
# FURTHER INVESTIGATION ...
# THERE ARE EXAMPLES WHERE 
cols = [#'Biological sample category',
 # 'Protein accession numbers',
 # 'Assigned',
 # 'Spectrum name',
 # 'Protein identification probability',
 # 'Peptide identification probability',
 'Previous amino acid',
 'Next amino acid',
 'prev_aa_fetched',
 'next_aa_fetched',
 # 'Number of enzymatic termini',
 # 'Fixed modifications identified by spectrum',
 # 'Variable modifications identified by spectrum',
 'Peptide start index',
 'Peptide stop index',
 'start_fetched',
 'stop_fetched',
 # 'Other Proteins',
 'pept',
 'GN',
 'GN_fetched',
 'OS',
 'locus',
 'prot_name',
 'uid',
 'uid_fetched',
 'fetchid',
 'crit_GN',
 'crit_uid_full',
 'crit_uid_maj',
 'crit_pept_in',
 'crit_start',
 'crit_stop',
 'crit_prev_aa',
 'crit_next_aa']
# spec_info_with_fetch.drop_duplicates(subset=cols)[cols].head(40)


# how many crits are True ...


crt_cols = [_ for _ in spec_info_with_fetch.columns if ('crit_' in _)]
spec_info_with_fetch['comb_crt'] = spec_info_with_fetch[ crt_cols ].sum(axis=1)

spec_by_prot = spec_info_with_fetch.drop_duplicates(subset=cols).groupby(by='prot_name',sort=False)
rrr=spec_by_prot.apply( lambda x: x.groupby(by='fetchid',sort=False)['comb_crt'].mean().max() )
# would be interesting to look at these ...
print rrr[rrr<4]
#[cols].head(30)
print spec_info_with_fetch[cols][spec_info_with_fetch['prot_name']=='HLA class I histocompatibility antigen, alpha chain E'].iloc[0]
# spec_by_prot.apply(lambda x: x.drop_duplicates(subset=cols)).groupby(by='fetchid',sorted=False)


#########################################################################################################
#                                grouped by pept                                                        #
#########################################################################################################
spec_by_pept = spec_info_with_fetch.drop_duplicates(subset=cols).groupby(by='pept',sort=False)
uuu=spec_by_pept.apply( lambda x: x.groupby(by='fetchid',sort=False)['comb_crt'].mean().max() )
yyy=spec_by_pept.apply( lambda x: x.groupby(by='fetchid',sort=False)['comb_crt'].max() )
# # would be interesting to look at these ...
# print uuu[uuu<4]
# #[cols].head(30)
# print spec_info_with_fetch[cols][spec_info_with_fetch['prot_name']=='HLA class I histocompatibility antigen, alpha chain E'].iloc[0]




# #
# ############################
# results = []
# for idx in record['IdList']:
#     print "fetching",idx
#     handle = Entrez.efetch(db='protein',id=idx,rettype='gb',retmode='text')
#     filelike = StringIO(handle.read())
#     seqrec = SeqIO.read(filelike,format='gb')
#     results.append( seqrec )
#########################################################################################
# # The features that Reid is asking about are accessible through:
# r0 = results[0]
# f4 = r0.features[4]
# f4.qualifiers['region_name']























