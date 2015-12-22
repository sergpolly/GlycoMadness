import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import pandas as pd
import numpy as np
import requests
import ms_module as ms


session = requests.Session()


pep_info = pd.read_csv("../peptides.xls",sep='\t')
# interesting column names ...
# cols = ['Protein accession numbers','Assigned','Other Proteins']


# BRIEF PEPTIDE SUMMARY ...
peptides_num = pep_info.shape[0]
peptides_unambiguously_assigned = pep_info['Assigned'].sum()
# getting the number of extractable Uid-s ... 
check_uid = lambda line: bool( line.split('|')[1] if len(line.split('|'))>1 else None )
extractable_uids = pep_info['Protein accession numbers'].apply(check_uid).sum()
# number of peptides matched amboguously ...
check_ambig_uid = lambda row: bool(True if (len(row[0].split('|'))>1)and(len(str(row[1]).split('|'))>1) else False)
ambig_uids_extracted = pep_info[['Protein accession numbers','Other Proteins']].apply(check_ambig_uid,axis='columns').sum()

print
print "Simple report "
print
print
print "Peptides detected %d"%peptides_num
print "Out of these %d are unambiguously assigned (including Unknown)"%peptides_unambiguously_assigned
print 
print "Number of peptide with Uniprot id that can be extracted %d"%extractable_uids
print "Thus number of strict Unknowns is %d"%(peptides_num - extractable_uids)
print
print "Ambigously matching peptides with extractable Uniprot ids %d"%ambig_uids_extracted
print
print
print
print



# extracting UID from protein accession numbers ...
# this way we return None for the Unknown entries ...
extract_uid = lambda line: line.split('|')[1] if len(line.split('|'))>1 else None
# get a single unique Uniprot ID ...
pep_info['uid'] = pep_info['Protein accession numbers'].apply(extract_uid)
# fetch protein sequence for each of the Uid-s ...
print "fetchign from uniprot.org ..."
pep_info['fasta'] = pep_info['uid'].apply(lambda _: ms.get_uniprot(session,_))
print "fetchign complete"
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







