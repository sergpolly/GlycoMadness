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
parser.add_argument("--verbose", help="verbose output", action="store_true")
parser.add_argument("--prefix", help="specify common part of the path for peptide and spectrum files")
parser.add_argument("--email", help="Provide your email for NCBI servers abuse-feedback")
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
#
pep_info = pd.read_csv(pep_fname,sep='\t')
# parse protein names to get GeneName Organisms etc ...
parse_pname = lambda pn: pd.Series(ms.parse_prot_name(pn))
# add columns with the parsed information on the proteins ...
pep_info = pep_info.merge(pep_info['Protein name'].apply(parse_pname), left_index=True, right_index=True)



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
# GENERATE AN ARRAY OF ACCESS TERMS:
term_func = lambda (gn,os): "RecName[Title] AND \"%s\"[Gene Name] AND \"%s\"[Organism]"%(gn,os)
# GNs of some proteins evaulted to NaN, skip them ...
gn_notnull = pep_info[pep_info['GN'].notnull()]
gn_os_for_terms = gn_notnull[['GN','OS']].drop_duplicates().reset_index(drop=True)
terms_array = gn_os_for_terms.apply(term_func,axis=1)
#
def search_term(term):
    # once all set, try actually searching something ...
    handle = Entrez.esearch(db="protein", term=term)
    # handle=Entrez.esearch(db="protein",term="RecName[TITL] AND PPT1[GENE] AND \"Homo sapiens\"[ORGN]")
    record = Entrez.read(handle)
    handle.close()
    #
    # fetch_idlist.append(record['IdList'])
    #
    print "Sending request for: %s ..."%term
    print "Results:",record['IdList']
    #
    # TRY SETTING A DELAY OR SOMETHING, IN CASE NCBI WOULD START COMPLAINING ...
    # NCBI might not like the ever frequent requests to its servers - that's why!
    time.sleep(500./1000.) # try to sleep for 100 miliseconds to make it easier on NCBI ...
    return pd.Series(list(record['IdList']))
    ########################################
multicol_fetch_res = terms_array.apply(search_term)
# stacked - multiindex is in use, level=1 index is the column name from 'multicol_mfunc_result' ...
unrolled_fetch_res = multicol_fetch_res.stack()
# index is no longer uniq after the following operation ... we dropped inner indexing part.
# IMPORTANT: indexes correspond to those from the original df (use it for merging later) ...
unrolled_origfetch = pd.DataFrame(unrolled_fetch_res.reset_index(level=1,drop=True),columns=['fetchid',])
# merge unrolled_origindex (a single column with ambiguous index) with the original df ...
# 'unrolled_origindex' must be DataFrame to merge: Series are not mergeable for some reason ...
# the whole df is to be unrolled after the following operation.
unrolled_fetchdf = gn_os_for_terms.merge(unrolled_origfetch,left_index=True,right_index=True).reset_index(drop=True)
#
#
#
# """Simply looking into DSC2 and DSC3 case tell us that there have to be sophisticated algorithm to choose
# the best matching fetchid for each gene name. Simply taking the first answer wouldn't work the best."""
# TO BE CONTINUED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pep_info_with_fetch = pep_info.merge(unrolled_fetchdf[['GN','fetchid']],how='outer', on='GN').reset_index(drop=True)
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
#
print
print "Posting and fetching genebank records corresponding to the available FetchIDs from the Protein DB ..."
pulled_gb_recs_fname = os.path.join( pep_path, "pulled_proteins.gb" )
batch_size = 60
attempts_limit = 3
# THEN WE'D NEED TO DO POST AND ONLY AFTER EFETCH ...
search_results = Entrez.read( Entrez.epost("protein", id=",".join( unrolled_fetchdf['fetchid'].unique() )) )
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
# download results in batches using history and coockies ....
count = unrolled_fetchdf.shape[0]
out_handle = open(pulled_gb_recs_fname, "w")
for start in range(0, count, batch_size):
    end = min(count, start+batch_size)
    print("Going to download record %i to %i" % (start+1, end))
    attempt = 0
    while attempt < attempts_limit:
        attempt += 1
        try:
            fetch_handle = Entrez.efetch(db="protein", rettype="gb", retmode="text",
                                         retstart=start, retmax=batch_size,
                                         webenv=webenv, query_key=query_key)
            break # skip subsequent attempts is succeeded ...
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %d of %d"%(attempt,attempts_limit))
                # attempt += 1
                time.sleep(15)
            else:
                print "oh Shut! %d"%attempt
                raise
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()
#
print "Fetched genebank records are stored in %s."%pulled_gb_recs_fname
print "Check for BioPython gb consistency before processing ..."
print "THE END"

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























