import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import pandas as pd
import numpy as np
import requests
import ms_module as ms
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
##########################################################
def combine_uids(row):
    first  = '' if pd.isnull(row[0]) else ('' if len(row[0].split('|'))<=1 else row[0])
    second = row[1] if pd.notnull(row[1]) else ''
    return ','.join([first,second])
# # let's generate exhaustive list of all Uniprot Ids that are present in the pe_info ...
uid_list = []
for uids in pep_info[['Protein accession numbers','Other Proteins']].apply(combine_uids,axis='columns').unique():
    for uid in uids.strip(',').split(','):
        uid_list.append(uid)
# once uids are extracted ...
extract_uid = lambda _: _.split('|')[1]
unique_uid_list = [extract_uid(_) for _ in np.unique(uid_list)]
#   to be continued ...
#   to be continued ...
#   to be continued ...
#   to be continued ...
#   to be continued ...
#   to be continued ...
#   to be continued ...
#   to be continued ...

# extracting UID from protein accession numbers ...
# this way we return None for the Unknown entries ...
extract_uid = lambda line: line.split('|')[1] if len(line.split('|'))>1 else None
# get a single unique Uniprot ID ...
pep_info['uid'] = pep_info['Protein accession numbers'].apply(extract_uid)
# fetch protein sequence for each of the Uid-s ...
fetching = False
if fetching:
    print "fetching from uniprot.org ..."
    pep_info['fasta'] = pep_info['uid'].apply(lambda _: ms.get_uniprot(session,_))
    print "fetching complete"
    # Align peptide sequence to the extracted protein sequence and find the peptide starting position ...
    pep_info['my_start'] = pep_info[ ['Peptide sequence','fasta'] ].apply(lambda _:ms.stupid_aligner(*_),axis='columns')



# c = ['Protein name',
# 'Protein accession numbers',
# 'Database sources',
# 'Exclusive unique peptide count',
# 'Peptide sequence',
# 'Previous amino acid',
# 'Next amino acid',
# 'Peptide start index',
# 'Peptide stop index',
# 'Star Category',
# 'Assigned',
# 'Other Proteins',
# 'uid']
##########################################
# c2 = ['Protein name',
# 'Peptide sequence',
# 'Other Proteins']
##########################################
# ['Protein name','Other Proteins']
##########################################
# sss = ">sp|P15586-2|GNS_HUMAN Isoform 2 of N-acetylglucosamine-6-sulfatase OS=Homo sapiens GN=GNS\n\
# MRLLPLAPGRLRRGSPRHLPSCSPALLLLVLGGCLGVFGVAAGTRRPNVVLLLTDDQDEV\
# LGGMYVPSALCCPSRASILTGKYPHNHHVVNNTLEGNCSSKSWQKIQEPNTFPAILRSMC\
# GYQTFFAGKYLNEYGAPDAGGLEHVPLGWSYWYALEKNSKYYNYTLSINGKARKHGENYS\
# VDYLTDVLANVSLDFLDYKSNFEPFFMMIATPAPHSPWTAAPQYQKAFQNVFAPRNKNFN\
# IHGTNKHWLIRQAKTPMTNSSIQFLDNAFRKRWQTLLSVDDLVEKLVKRLEFTGELNNTY\
# IFYTSDNGYHTGQFSLPIDKRQLYEFDIKVPLLVRGPGIKPNQTSKMLVANIDLGPTILD\
# IAGYDLNKTQMDGMSLLPILRGASNLTWRSDVLVEYQGEGRNVTDPTCPSLSPGVSQCFP\
# DCVCEDAYNNTYACVRTMSALWNLQYCEFDDQEVFVEVYNLTADPDQITNIAKTIDPELL\
# GKMNYRLMMLQSCSGPTCRTPGVFDPGYRFDPRLMFSNRGSVRTRRFSKHLL"







