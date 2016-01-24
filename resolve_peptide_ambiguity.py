import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import pandas as pd
import numpy as np
import ms_module as ms
import re
#
# #
# pep_fname = "../peptides.xls"
# fasta_fname = "dec23_1701.fasta"
if len(sys.argv)<4:
    print "Command line arguments required! Launch example:"
    print "%s ../input_peptides.xls required_prot_dec23_1701.fasta uniq_peptides_catalog.sv"%sys.argv[0]
    sys.exit(1)
pep_fname = sys.argv[1]
fasta_fname = sys.argv[2]
out_fname = sys.argv[3]
#
pep_info = pd.read_csv(pep_fname,sep='\t')
#
fasta = SeqIO.to_dict(SeqIO.parse(fasta_fname,"fasta"),key_function=lambda _: _.id.split('|')[1])
#
# at first let's get a list of unique peptide sequences ...
uniq_pept = pep_info['Peptide sequence'].unique()


# now let's find all Uniprot ids associated with each peptide ...
def extract_uids(peptide_seq,pep_dat_info,columns = ["Protein accession numbers","Other Proteins"]):
    """ looking for all uid in specified columns excluding Unknowns."""
    pep_dat = pep_dat_info[pep_dat_info["Peptide sequence"]==peptide_seq]
    uid_ambig_list = []
    for col in columns:
        for uids in pep_dat[col].unique():
            # avoid missing data ...
            if pd.notnull(uids):
                # break them up (they are comma-separated) ...
                for uid in uids.strip(',').split(','):
                    # avoid Unknowns ...
                    if len(uid.split('|'))>1:
                        uid_ambig_list.append(uid)
    #############################################
    uid_uniq_list = [ uid_uniq_item.replace(';',' ') for uid_uniq_item in set(uid_ambig_list) ]
    to_return = ';'.join(uid_uniq_list) if uid_uniq_list else None
    return (pep_dat,to_return)


def get_single_fasta(uids,fa_dict):
    #######################
    uid_list = [uid.split('|')[1] for uid in uids.strip(',').split(',')]
    fasta_list = [fa_dict[uid] for uid in uid_list]
    fasta_len_list = [ len(sr.seq) for sr in fasta_list ]
    max_len_index = np.argmax(fasta_len_list)
    #######################
    return (uid_list[max_len_index],fasta_len_list[max_len_index],fasta_list[max_len_index])


# extract_uids - is working
interesting_peptide = []
uids_list = []
uid_of_maxlen_list = []
prot_len = []
prot_fasta = []
pept_positions = []
prot_name = []
uniq_pept_count = []
pept_probab = []
#
aa_before = []
aa_after = []
for pept in uniq_pept:
    pep_dat_pept,uids = extract_uids(pept,pep_info)
    if uids:
        interesting_peptide.append(pept)
        uids_list.append(uids)
        #################################
        _1,prot_len_fasta,prot_seq_fasta = get_single_fasta(uids,fasta)
        uid_of_maxlen_list.append(_1)
        prot_len.append(prot_len_fasta)
        prot_fasta.append(str(prot_seq_fasta.seq))
        peptide_start_in_protein = ms.stupid_aligner(pept,prot_seq_fasta) # 1-based ...
        peptide_stop_in_protein = peptide_start_in_protein + len(pept) # 1-based ...
        pept_positions.append(peptide_start_in_protein)
        prot_name.append(prot_seq_fasta.description.replace(',',' ')) # long protein name here ...
        # uniq peptide count taken from pep_dat_pept, for definition look up extract_uids...
        uniq_pept_count_val = pep_dat_pept['Exclusive unique peptide count'].unique()[0]
        # uniq_pept_count_val, = pep_dat_pept['Exclusive unique peptide count'][pep_dat_pept['Exclusive unique peptide count']>0].unique()
        uniq_pept_count.append(uniq_pept_count_val)
        # some kind of peptide probability (like a quality score from experimental data)?!
        pept_probab_val, = pep_dat_pept['Best Peptide identification probability'].unique()
        pept_probab.append(pept_probab_val)
        #################################
        # BEWARE: 1-BASED INDEXING ALL THE WAY ACROSS SO FAR...
        # peptide can start right at N-terminus, so there will be no AminoAcid preceding it, call it a START
        aa_before.append(str(prot_seq_fasta.seq)[peptide_start_in_protein-2] if peptide_start_in_protein>1 else 'START')
        # peptide can end right at C-terminus, so there will be no AminoAcid after it, call it an END
        aa_after.append(str(prot_seq_fasta.seq)[peptide_stop_in_protein+0] if peptide_stop_in_protein<prot_len_fasta else 'END')
#########################################
dict_df = {
    "pept":interesting_peptide,
    "all_uids":uids_list,
    "uid_max":uid_of_maxlen_list,
    "protlen":prot_len,
    "peptide_start":pept_positions,
    "prot_seqrec":prot_fasta,
    "aa_before":aa_before,
    "aa_after":aa_after,
    "prot_name":prot_name,
    "uniq_pept_count":uniq_pept_count,
    "pept_probab":pept_probab
}
##############################################################################################################################
pep_df = pd.DataFrame(dict_df)
# within-cell separators are ';' now ...
# # let's change within a cell separators to ';' instead of ',', which is used for column separation ...
# for col,dtype in pep_df.dtypes.iteritems():
#     if dtype=='object':
#         pep_df[col].str.replace(';',' ').str.replace(',',';')
##############################################################################################################################
pep_df.to_csv(out_fname,index=False)


# print "Beware! Some columns have string values with the comma-characters ',' in them, turned out pandas "
# print "deals with the problem graciously, placing such value in quotes, that makes reading such csv files "
# print "an easy task. At least both pandas read_csv and Apple Numbers(like Excel) interpret the situation correctly!"














