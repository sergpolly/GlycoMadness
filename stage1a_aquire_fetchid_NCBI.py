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
parser.add_argument("-p","--pept_summary", help="speicfy peptide summary file name (with/without path)",required=True)
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
    pep_fname = os.path.join( args.prefix, args.pept_summary )
else:
    pep_fname = args.pept_summary
# get the common path for later use ...
pep_path = os.path.dirname(pep_fname)
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
pep_info = pd.read_csv(pep_fname,sep=separator)
# parse protein names to get GeneName Organisms etc ...
parse_pname = lambda pn: pd.Series(ms.parse_prot_name(pn))
# add columns with the parsed information on the proteins ...
pep_info = pep_info.merge(pep_info['Protein name'].apply(parse_pname), left_index=True, right_index=True)
# there is a little fix, for datasets Glyco4 and beyond ...
if args.swiss:
    pep_info['locus'] = pep_info['Protein accession numbers']
    pep_info['uid']   = None



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
# # this would be all possible fields that you can search in the protein database ...
# print "These are the fields one can search in db='protein':"
# handle = Entrez.einfo(db='protein')
# prot_db_fields = Entrez.read(handle)
# for idx,field in enumerate(prot_db_fields['DbInfo']['FieldList']):
#     print idx+1, field['Name'], field['FullName'], field['Description']
# handle.close()
#


# def term_func_GN(gn,os):
#     if pd.isnull(gn) and pd.isnull(os):
#         return None
#     else:
#         term = ["RecName[Title]",]
#         if not pd.isnull(gn):
#             term.append("\"%s\"[Gene Name]"%gn)
#         if not pd.isnull(os):
#             term.append("\"%s\"[Organism]"%os)
#         return ' AND '.join(term)
# def term_func_LOCUS(loc,os):
#     if pd.isnull(loc) and pd.isnull(os):
#         return None
#     else:
#         term = ["RecName[Title]",]
#         if not pd.isnull(loc):
#             term.append("\"%s\"[locus]"%loc)
#         if not pd.isnull(os):
#             term.append("\"%s\"[Organism]"%os)
#         return ' AND '.join(term)
# # GENERATE AN ARRAY OF ACCESS TERMS:
# # locus search term ...
# # choose 
# if args.search == 'locus':
#     term_func = lambda row: term_func_LOCUS(row['locus'],row['OS'])
# elif args.search == 'GN':
#     term_func = lambda row: term_func_GN(row['GN'],row['OS'])
# else:
#     print
#     print "search argument is wrong, read help:"
#     print
#     parser.print_help()
#     sys.exit(1)

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
# gn_notnull = pep_info[pep_info['GN'].notnull()]
search_used_features = ['GN','locus','OS','uid']
search_term_arguments = pep_info[search_used_features].drop_duplicates().reset_index(drop=True)
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
pep_info_with_fetch = pep_info.fillna('').merge(unrolled_fetchdf,how='outer', on=search_used_features).reset_index(drop=True)
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
print "Storing updated spectrum file to include FetchIDs for every GeneName ..."
pep_info_with_fetch_fname = os.path.join( pep_path, 'peptides_with_fetch.csv' )
pep_info_with_fetch.to_csv(pep_info_with_fetch_fname,index=False)
print "file is stored as %s"%pep_info_with_fetch_fname
print
#












