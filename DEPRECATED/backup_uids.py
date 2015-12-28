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
#
session = requests.Session()
#
#
pep_info = pd.read_csv("../peptides.xls",sep='\t')
# interesting column names ...
# cols = ['Protein accession numbers','Assigned','Other Proteins']
spec_info = pd.read_csv("../specs.xls",sep='\t')
#
#
# BRIEF PEPTIDE SUMMARY ...
peptides_num = pep_info.shape[0]
peptides_unambiguously_assigned = pep_info['Assigned'].sum()
# getting the number of extractable Uid-s ... 
check_uid = lambda line: bool( line.split('|')[1] if len(line.split('|'))>1 else None )
extractable_uids = pep_info['Protein accession numbers'].apply(check_uid).sum()
# number of peptides matched amboguously ...
check_ambig_uid = lambda row: bool(True if (len(row[0].split('|'))>1)and(len(str(row[1]).split('|'))>1) else False)
ambig_uids_extracted = pep_info[['Protein accession numbers','Other Proteins']].apply(check_ambig_uid,axis='columns').sum()
#
print
print "--- PEPTIDE REPORT ---"
print
print "Total peptides detected %d"%peptides_num
print "There are %d Unknown (not-mapped) peptides out of these %d"%((peptides_num - extractable_uids),peptides_num)
print "Accordingly, %d peptides have at least 1 associated protein in the databsae."%extractable_uids
print "Among them, %d are claimed unassigned (ambiguously assigned), while %d of them have extractable ids"%((peptides_num - peptides_unambiguously_assigned),ambig_uids_extracted)
print "However, there are many peptides assigned to more than 1 proteins (i.e., multiple Uniprot ids)"
print "Yet, these ids are usually reffering to the same exact protein of its isoform."
print
print "--- PEPTIDE REPORT COMPLETE ---"
print
#
##########################################################
print 
print "comparing spectrum data with the peptide data ..."
unique_ids_pept = pep_info['Protein accession numbers'].unique().shape[0]
unique_ids_spec = spec_info['Protein accession numbers'].unique().shape[0]
###############################################################
combine_uids = lambda row: ','.join(map(str,row))# if row[1] is not None else row[0]
###############################################################
unique_ids_pept_ambig = pep_info[['Protein accession numbers','Other Proteins']].apply(combine_uids,axis='columns').unique().shape[0]
unique_ids_spec_ambig = spec_info[['Protein accession numbers','Other Proteins']].apply(combine_uids,axis='columns').unique().shape[0]
###############################################################
print "Number of unique uniprot ids in SPECTRUM: %d and in PEPTIDE: %d"%(unique_ids_spec,unique_ids_pept)
print "Same, but combining original accessions with 'Other Proteins'"
print "Number of unique uniprot ids in SPECTRUM: %d and in PEPTIDE: %d"%(unique_ids_spec_ambig,unique_ids_pept_ambig)
print "Those numbers included Unknowns!"
##########################################################
#


# let's just get all Uids from columns 'Protein accession numbers' and 'Other Proteins' ...
uid_ambig_list = []
for uids in pep_info['Protein accession numbers'].unique():
    # avoid missing data ...
    if pd.notnull(uids):
        # break them up ...
        for uid in uids.strip(',').split(','):
            # avoid Unknowns ...
            if len(uid.split('|'))>1:
                uid_ambig_list.append(uid)
for uids in pep_info['Other Proteins'].unique():
    # avoid missing data ...
    if pd.notnull(uids):
        # break them up ...
        for uid in uids.strip(',').split(','):
            # avoid Unknowns ...
            if len(uid.split('|'))>1:
                uid_ambig_list.append(uid)
# now we got them ...

# once uids are extracted ...
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

SeqIO.write(unique_fasta_list,"backup_uids_output.fasta","fasta")




