"""
Script is searching for unique UNIPROT id-s in columns "Protein accession numbers" and "Other Proteins",
and then fetches corresponding protein sequences from uniprot.org to create a local sequence database.
"""
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import pandas as pd
import numpy as np
import requests
import ms_module as ms
import StringIO
#
# requests session to buffer requests to uniprot.org
session = requests.Session()
#
input_fname = "../peptides.xls"
output_fname = "backup_uids_output.fasta"
# interesting column names ...
# cols = ['Protein accession numbers','Assigned','Other Proteins']
pep_info = pd.read_csv(input_fname,sep='\t')

# let's just get all Uids from columns 'Protein accession numbers' and 'Other Proteins' ...
uid_ambig_list = []
for uids in pep_info['Protein accession numbers'].unique():
    # avoid missing data ...
    if pd.notnull(uids):
        # break them up (they are comma-separated) ...
        for uid in uids.strip(',').split(','):
            # avoid Unknowns ...
            if len(uid.split('|'))>1:
                uid_ambig_list.append(uid)
for uids in pep_info['Other Proteins'].unique():
    # avoid missing data ...
    if pd.notnull(uids):
        # break them up (they are comma-separated) ...
        for uid in uids.strip(',').split(','):
            # avoid Unknowns ...
            if len(uid.split('|'))>1:
                uid_ambig_list.append(uid)
# now we got them ...

# once uids are extracted  make sure they are unique ...
extract_uid = lambda _: _.split('|')[1]
unique_uid_list = [extract_uid(_) for _ in np.unique(uid_ambig_list)]

print
print "Start fetching protein sequences from Uniprot database ..."
unique_fasta_list = []
# see http://www.uniprot.org/help/programmatic_access for details
# the way we form Uniprot ID request ...
get_uid_url = lambda _: "http://www.uniprot.org/uniprot/%s.fasta"%_
# web request for a given uid ...
for uid in unique_uid_list:
    uid_url = get_uid_url(uid)
    # make a request ...
    req_res = session.get(uid_url)
    # check request status ...
    if req_res.status_code==200 and bool(req_res.content):
        # to be read by __SeqIO ...
        string_as_handle = StringIO.StringIO(req_res.content)
        seq_rec = SeqIO.read(string_as_handle,'fasta')
        unique_fasta_list.append(seq_rec)
    elif req_res.status_code==200:
        print ms.__web_request_status_collection[req_res.status_code]
        print "... But, the content is empty for accession number %s!"%uid
        sys.exit(1)
    elif req_res.status_code in ms.__web_request_status_collection:
        print ms.__web_request_status_collection[req_res.status_code]
        sys.exit(1)
    else:
        print "Unknown status code returned!"
        sys.exit(1)
print
print "Fetching is complete, writing output to fasta file ..."
#############################################################################
SeqIO.write(unique_fasta_list,output_fname,"fasta")




