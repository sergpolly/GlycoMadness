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
parser.add_argument("-f","--raw_fetch", help="speicfy input data fname (with fetchid column(!), with/without path)",required=True)
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
    raw_info_with_fetch_fname = os.path.join( args.prefix, args.raw_fetch )
else:
    raw_info_with_fetch_fname = args.raw_fetch
# get the common path for later use ...
raw_path = os.path.dirname(raw_info_with_fetch_fname)
#
# don'r forget to provide you email
Entrez.email = args.email if args.email else "your_email@mail_server.com"

#
# peptides_with_fetch.csv
# raw_info_with_fetch_fname
# raw_info_with_fetch
raw_info_with_fetch = pd.read_csv(raw_info_with_fetch_fname)
assert 'fetchid' in raw_info_with_fetch.columns


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
# H prot_name, parsed out human-readable name from 'Protein Name'
# H1 gene_name, parsed out GN=xxx from 'Protein Name'
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
print
print "Posting and fetching genebank records corresponding to the available FetchIDs from the Protein DB ..."
pulled_gb_recs_fname = os.path.join( raw_path, "pulled_proteins.gb" )
batch_size = 60
attempts_limit = 3
# THEN WE'D NEED TO DO POST AND ONLY AFTER EFETCH ...
# there might be some EMPTY fetchids ...
non_empty_fetchids  = raw_info_with_fetch['fetchid'][raw_info_with_fetch['fetchid'].notnull()].apply(int)
with_empty_fetchids = raw_info_with_fetch[raw_info_with_fetch['fetchid'].isnull()]
#
print
print "BEWARE! There are %d empty fetchids ..."%with_empty_fetchids.shape[0]
print with_empty_fetchids[['Protein Name','Peptide Sequence']]
print
#
search_results = Entrez.read( Entrez.epost("protein", id=",".join( non_empty_fetchids.apply(str).unique() )) )
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
# download results in batches using history and coockies ....
count, = non_empty_fetchids.unique().shape
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




























