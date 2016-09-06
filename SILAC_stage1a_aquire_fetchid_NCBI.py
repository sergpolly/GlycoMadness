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

import argparse

# do some arguments parsing to make the script looks civilized ...
parser = argparse.ArgumentParser()
parser.add_argument("-r",
                    "--raw_silac",
                    help = "speicfy raw data SILAC file name (with/without path)\n \
                    In reality, any file with 'Accession Number' and 'Protein Name' columns would go.",
                    required=True)
# we don't need spectrum file for downloading proteins, it is too redundant for that purpose ...
# parser.add_argument("--search", help="Specify NCBI search term: locus or GN (Gene Name)?", default="locus")
parser.add_argument("--verbose", help="verbose output", action="store_true")
parser.add_argument("--swiss", help="Swiss prot input (glyco4 and beyond ...)", action="store_true")
parser.add_argument("--prefix", help="specify common part of the path for peptide and spectrum files")
parser.add_argument("--email", help="Provide your email for NCBI servers abuse-feedback")
parser.add_argument("--separator", help="speicfy separator type in the input data",default='tab')
args = parser.parse_args()

# print args
###############################################
if args.verbose:
    print "Verbose output is to follow ...\n\n"
###############################################
if args.prefix is not None:
    raw_fname = os.path.join( args.prefix, args.raw_silac )
else:
    raw_fname = args.raw_silac
# get the common path for later use ...
raw_path = os.path.dirname(raw_fname)
#
# don'r forget to provide you email
Entrez.email = args.email if args.email else "your_email@mail_server.com"

#
if args.separator == "tab":
    separator = '\t'
elif args.separator == "comma":
    separator = ','
else:
    separator = '\t'
#
raw_info = pd.read_csv(raw_fname,sep=separator)
# parse protein names to get GeneName Organisms etc ...
parse_pname = lambda pn: pd.Series(ms.parse_prot_name(pn))
# add columns with the parsed information on the proteins ...
raw_info = raw_info.merge(raw_info['Protein Name'].apply(parse_pname), left_index=True, right_index=True)
# there is a little fix, for datasets Glyco4 and beyond ...
if args.swiss:
    raw_info['locus'] = raw_info['Accession Numbers'].str.strip('[]')
    raw_info['uid']   = None






# "RecName[Title] AND UID[accession] AND LOCUS[locus] AND GN[Gene Name] AND OS[Organism]"
def term_func(row):
    locus,GN,uid,OS = row['locus'],row['GN'],row['uid'],row['OS']
    if pd.isnull(locus) and pd.isnull(GN) and pd.isnull(uid):
        return None
    if pd.notnull(locus):
        term_l = "RecName[Title] AND \"%s\"[locus] AND \"%s\"[Organism]"%(locus,OS) if pd.notnull(OS) else "RecName[Title] AND \"%s\"[locus]"%locus
        term_u = ""
        term_g = ""
    if pd.notnull(uid):
        # trim version off of the accession numbers (UniprotID) ...
        uid = uid.split('-')[0]
        uid = uid.split('.')[0]
        term_u = "\"%s\"[accession] AND \"%s\"[Organism]"%(uid,OS) if pd.notnull(OS) else "\"%s\"[accession]"%uid
        term_g = ""
    if pd.notnull(GN):
        term_g = "RecName[Title] AND \"%s\"[Gene Name] AND \"%s\"[Organism]"%(GN,OS) if pd.notnull(OS) else "RecName[Title] AND \"%s\"[Gene Name]"%GN
    # returning ....
    return ','.join([term_l,term_u,term_g])


def get_term_type(term):
    if "[Gene Name]" in term:
        return "gene"
    if "[locus]" in term:
        return "locus"
    if "[accession]" in term:
        return "uid"



def get_best_fetchids(df):
    fid_hist = df.groupby('fetchid').size()
    if (fid_hist==3).any():
        return df[ df['fetchid'].isin(fid_hist[fid_hist==3].index) ]
    else:
        return df



# # GNs of some proteins evaulted to NaN, skip them ...
# gn_notnull = raw_info[raw_info['GN'].notnull()]
search_used_features = ['GN','locus','OS','uid']
search_term_arguments = raw_info[search_used_features].drop_duplicates().reset_index(drop=True)
terms_array = search_term_arguments.apply(term_func,axis=1).dropna() # dropna just in case ...
# dropna, keeps index, so terms_array is still mergeable with search_term_argumentsby index ...
#
def search_term(terms,speedup=True,wait_time=10.0):
    if pd.isnull(terms):
        return pd.Series([])
    else:
        terms = terms.strip(',')
        print "Sending request for: %s ..."%terms
        fetchid_list = []
        prev_list = None
        for term in terms.split(','):
            if not term:
                continue
            # infer term type ...
            term_type = get_term_type(term)
            # once all set, try actually searching something, using several attempts ...
            for attempt in range(10):
                try:
                    if attempt >= 1:
                        print "Sending next request %d for: %s ..."%(attempt+1,terms)
                    #
                    handle = Entrez.esearch(db="protein", term=term, usehistory="y")
                    # handle=Entrez.esearch(db="protein",terms="RecName[TITL] AND PPT1[GENE] AND \"Homo sapiens\"[ORGN]")
                    record = Entrez.read(handle)
                except:
                    # TRY SETTING A DELAY OR SOMETHING, IN CASE NCBI WOULD START COMPLAINING ...
                    # NCBI might not like the ever frequent requests to its servers - that's why!
                    print "wait for %.1f sec ..."%wait_time
                    time.sleep(wait_time*1000./1000.) # try to sleep for 100 miliseconds to make it easier on NCBI ...
                    continue
                # break only when you succeed ...
                break
            # #
            # #
            handle.close()
            print "Results:",record['IdList']
            # fetchid list will containt information on the source where it has been obtained ...
            new_id_list = list(record['IdList'])
            # for speedup, we are assuming that if search returned 2 identical (notnull) resutls,
            # we might skip the search and proceed to the next bunch of terms ...
            if speedup and (prev_list==new_id_list) and bool(prev_list):
                break
            elif speedup:
                prev_list = new_id_list
            # otherwise , simply proceed ...
            fetchid_list += map(lambda x:':'.join((term_type,x)), new_id_list)
            # fetch_idlist.append(record['IdList'])
            #
            #
        # TRY SETTING A DELAY OR SOMETHING, IN CASE NCBI WOULD START COMPLAINING ...
        # NCBI might not like the ever frequent requests to its servers - that's why!
        time.sleep(100./1000.) # try to sleep for 100 miliseconds to make it easier on NCBI ...
        return pd.Series(fetchid_list)
    ########################################
multicol_fetch_res = terms_array.apply(search_term)
# stacked - multiindex is in use, level=1 index is the column name from 'multicol_mfunc_result' ...
unrolled_fetch_res = multicol_fetch_res.stack()
# index is no longer uniq after the following operation ... we dropped inner indexing part.
# IMPORTANT: indexes correspond to those from the original df (use it for merging later) ...
unrolled_origfetch = unrolled_fetch_res.reset_index(level=1,drop=True).str.split(':',expand=True)
unrolled_origfetch = unrolled_origfetch.rename(columns={0:'source',1:'fetchid'})
# merge unrolled_origindex (a single column with ambiguous index) with the original df ...
# 'unrolled_origindex' must be DataFrame to merge: Series are not mergeable for some reason ...
# the whole df is to be unrolled after the following operation.
unrolled_fetchdf = search_term_arguments.merge(unrolled_origfetch,left_index=True,right_index=True).drop_duplicates().reset_index(drop=True)
#
#
# Now do some analysis using the source of data ...
# sss.groupby('o').apply(get_best_fetchids).drop_duplicates(['b','o']).reset_index(drop=True).drop('a',axis=1)
unrolled_fetchdf = unrolled_fetchdf.fillna('').groupby(search_used_features).apply(get_best_fetchids).drop('source',axis=1).drop_duplicates().reset_index(drop=True)
# by now unrolled_fetchdf must be ried of majority of redundancies ...
#
# """Simply looking into DSC2 and DSC3 case tell us that there have to be sophisticated algorithm to choose
# the best matching fetchid for each gene name. Simply taking the first answer wouldn't work the best."""
# TO BE CONTINUED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pep_info_with_fetch = raw_info.fillna('').merge(unrolled_fetchdf,how='outer', on=search_used_features).reset_index(drop=True)
# #########################
# # LOOKS LIKE UID IS A GREAT CANDIDATE TO DISCERN
# # BETWEEN RIGHT AND WRONG SEQUENCES RETURNED BY NCBI FOR A GENENAME REQUEST ...
# #################################
#############################
# # # TOUGH CHOISE ON WHICH PROTEIN TO FINALLY ASSIGN TO THE GN ...
# yet to be made ...
######################################
# OUTPUT THE UPDATED SPEC_INFO (WITH FETCHIDS) TO BE USED ON THE NEXT STAGE ...
print
print "FetchIDs are ready."
print "Storing updated raw data file to include FetchIDs for every GeneName ..."
raw_info_with_fetch_fname = os.path.join( raw_path, 'SILAC_plus_fetch.csv' )
pep_info_with_fetch.to_csv(raw_info_with_fetch_fname,index=False)
print "file is stored as %s"%raw_info_with_fetch_fname
print
#













