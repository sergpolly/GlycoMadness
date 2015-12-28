import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import pandas as pd
import numpy as np
import ms_module as ms
import re
#
#
pep_fname = "../peptides.xls"
spec_fname = "../specs.xls"
fasta_fname = "dec23_1701.fasta"
#
pep_info = pd.read_csv(pep_fname,sep='\t')
#
spec_info = pd.read_csv(spec_fname,sep='\t')
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
    uid_uniq_list = list(set(uid_ambig_list))
    to_return = ','.join(uid_uniq_list) if uid_uniq_list else None
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
        peptide_start_in_protein = ms.stupid_aligner(pept,prot_seq_fasta)
        peptide_stop_in_protein = peptide_start_in_protein + len(pept)
        pept_positions.append(ms.stupid_aligner(pept,prot_seq_fasta))
        prot_name.append(prot_seq_fasta.description.replace(',',' ')) # long protein name here ...
        # uniq peptide count taken from pep_dat_pept, for definition look up extract_uids...
        uniq_pept_count_val = pep_dat_pept['Exclusive unique peptide count'].unique()[0]
        # uniq_pept_count_val, = pep_dat_pept['Exclusive unique peptide count'][pep_dat_pept['Exclusive unique peptide count']>0].unique()
        uniq_pept_count.append(uniq_pept_count_val)
        # some kind of peptide probability (like a quality score from experimental data)?!
        pept_probab_val, = pep_dat_pept['Best Peptide identification probability'].unique()
        pept_probab.append(pept_probab_val)
        #################################
        aa_before.append(str(prot_seq_fasta.seq)[peptide_start_in_protein-1])
        aa_after.append(str(prot_seq_fasta.seq)[peptide_stop_in_protein+1] if peptide_stop_in_protein+1<prot_len_fasta else 'END')
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
##########################################
pep_df = pd.DataFrame(dict_df)


#############################################
# TO  BE CONTINUED ....

# connection between peptide info and spectrum info to be established ...
################################
spec_info['pept'] = spec_info['Peptide sequence'].str.upper()

# now let's check is all of the peptides from pep_df are in spec_info ...
spec_pept_present_index = pep_df['pept'].isin(spec_info['pept'])
print "FYI"
print " There are %d peptides are present in spectrum file out of %d in the peptide info file ..."%(sum(spec_pept_present_index),pep_df.shape[0])
print
spec_pept_uniq_num = spec_info['pept'].unique().size
peptdf_present_num = sum(pd.Series(spec_info['pept'].unique()).isin(pep_df['pept']))
print "And vice a versa"
print "There are %d spectrum peptides present in peptide file out of %d"%(peptdf_present_num,spec_pept_uniq_num)
print "Nevertheless, we move on with whatever we can do here ..."


# looking for Daemidated sites only ...
deamid = re.compile('[D,d]eamidat')
# enumerating all tractable peptides from pep_df ...
glyco_mod = []
for _,uniq_pept,pept_pos,prot_sr in pep_df[['pept','peptide_start','prot_seqrec']].itertuples():
    prot_seq = str(prot_sr.seq)
    pept_spectrum = spec_info[ spec_info['pept']==uniq_pept ]
    # let's check all present modifications ...
    modifs = ','.join(pept_spectrum['Variable modifications identified by spectrum'])
    modifs = [mod.strip() for mod in modifs.strip().split(',')]
    # and get those that are unique ...
    modifs = np.unique(modifs)
    # now extract type,position and value for each of them ...
    modifs = [ ms.parse_spectrum_modifications(mod) for mod in modifs if bool(deamid.search(mod)) ]
    # now extrating meaningfull glycosilation sites ...
    glyco_sites = []
    glyco_start = []
    for type_aa,pos,value in modifs:
        if (type_aa in ['n','N']) and (np.abs(value-3)<0.1):
            # zero-based index of the g-size ...
            gsite_start = pept_pos + pos - 1
            gsite_stop  = pept_pos + pos - 1 + 3
            glyco_start.append(gsite_start)
            glyco_sites.append(prot_seq[gsite_start:gsite_stop])
    ############################################################
    glyco_mod_str = ','.join([ gsite+("(%d)"%gstart) for gsite,gstart in zip(glyco_sites,glyco_start)])
    glyco_mod.append(glyco_mod_str)


pep_df['gsites'] = glyco_mod

print
print " We need to unroll glyco-sites and expand the DataFrame ..."
# we need to unroll those rows whith more than 1 items of gsites ...
pep_df_cols = ['all_uids', 'pept', 'peptide_start', 'prot_seqrec', 'protlen', 'uid_max', 'gsites','aa_before','aa_after']
pep_df_gsite_extracted = []
for _,all_uids,pept,pept_start,prot_sr,protlen,uid_max,gsites,aa_bef,aa_aft in pep_df[pep_df_cols].itertuples():
    # unrolling ...
    for gsite in gsites.split(','):
        pep_df_gsite_extracted.append((all_uids,pept,pept_start,prot_sr,protlen,uid_max,gsite,aa_bef,aa_aft))
########################
pep_df_ext = pd.DataFrame(pep_df_gsite_extracted,columns = pep_df_cols)
print "mission accomplished ..."


print
print "Now we'd need to collapse those identical glyco sites that orginiate from overlaing peptides ..."
gsite_uid_combined = pep_df_ext['gsites'] +':'+ pep_df_ext['all_uids']
gsites_uniq = gsite_uid_combined.unique()
print " There are %d unique (glyco-site, protein id) pairs that includes empty sites, however."%gsites_uniq.size
# now collapsing on those uniq sites ...
final_dataframe = []
for site_uniq in gsites_uniq:
    sites_index = (gsite_uid_combined == site_uniq)
    pep_site = pep_df_ext[sites_index]
    #
    pepts = ','.join( '('+pep_site['aa_before']+')'+pep_site['pept']+'('+pep_site['aa_after']+')' )
    all_uids, = pep_site['all_uids'].unique()
    pepts_start = ','.join(str(_) for _ in pep_site['peptide_start'])
    prot_seq = str(pep_site['prot_seqrec'].unique()[0].seq)
    prot_len, = pep_site['protlen'].unique()
    uid_max, = pep_site['uid_max'].unique()
    # make sure there is just a single gsite ...
    pep_site_gsite_uniq, = pep_site['gsites'].unique()
    # appending to final dataframe ...
    final_dataframe.append((pep_site_gsite_uniq,pepts,pepts_start,all_uids,prot_seq,prot_len,uid_max))
# turn to dataframe ...
final_dataframe = pd.DataFrame(final_dataframe,columns=['gsites', 'pept', 'peptide_start','all_uids', 'prot_seq', 'protlen', 'uid_max'])



# c = ['Protein accession numbers',
# 'Database sources',
# 'Exclusive unique peptide count',
# 'Exclusive unique spectrum count',
# 'Total spectrum count',
# 'Assigned',
# 'Spectrum name',
# 'Peptide sequence',
# 'Variable modifications identified by spectrum',
# 'Spectrum charge',
# 'Peptide start index',
# 'Peptide stop index',
# 'Exclusive',
# 'Other Proteins']
###########################################
# 'Protein accession numbers'
# 'Exclusive unique peptide count'
# 'Exclusive unique spectrum count'
# 'Total spectrum count'
# 'Spectrum name'
# 'Peptide sequence'
# 'Mascot Ion score'
# 'Mascot Identity score'
# 'Mascot Delta Ion Score'
# 'Number of enzymatic termini'
# 'Fixed modifications identified by spectrum'
# 'Variable modifications identified by spectrum'
# 'Observed m/z'
# 'Actual peptide mass (AMU)'
# 'Calculated +1H Peptide Mass (AMU)'
# 'Spectrum charge'
# 'Actual minus calculated peptide mass (AMU)'
# 'Actual minus calculated peptide mass (PPM)'
# 'Retention Time (sec)'
# 'Precursor Intensity'
# 'Total Ion Current'
###########################################
# cc = ['Protein accession numbers',
# # 'Exclusive unique peptide count',
# # 'Exclusive unique spectrum count',
# # 'Total spectrum count',
# # 'Spectrum name',
# 'Peptide sequence',
# # 'Mascot Ion score',
# # 'Mascot Identity score',
# # 'Mascot Delta Ion Score',
# # 'Number of enzymatic termini',
# 'Fixed modifications identified by spectrum',
# 'Variable modifications identified by spectrum',
# # 'Observed m/z',
# # 'Actual peptide mass (AMU)',
# # 'Calculated +1H Peptide Mass (AMU)',
# # 'Spectrum charge',
# # 'Actual minus calculated peptide mass (AMU)',
# # 'Actual minus calculated peptide mass (PPM)',
# # 'Retention Time (sec)',
# # 'Precursor Intensity',
# # 'Total Ion Current'
# ]
###########################################

















