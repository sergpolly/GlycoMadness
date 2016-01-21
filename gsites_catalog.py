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
# spec_fname = "../specs.xls"
# fasta_fname = "dec23_1701.fasta"
# uniq_pept_fname = "uniq_peptides_catalog.csv"
# #
if len(sys.argv)<6:
    print "Command line arguments are required! Launch example:"
    print "%s  ../peptides.xls  ../specs.xls  dec23_1701.fasta  uniq_peptides_catalog.csv  gsites_antology.csv"%sys.argv[0]
pep_fname =         sys.argv[1]
spec_fname =        sys.argv[2]
fasta_fname =       sys.argv[3]
uniq_pept_fname =   sys.argv[4]
out_fname =         sys.argv[5]
#
pep_info = pd.read_csv(pep_fname,sep='\t')
#
spec_info = pd.read_csv(spec_fname,sep='\t')
#
fasta = SeqIO.to_dict(SeqIO.parse(fasta_fname,"fasta"),key_function=lambda _: _.id.split('|')[1])
# ZERO-BASED NOTATION FOR PROTEINS INDEXING ENFORCED ...
pep_df = pd.read_csv(uniq_pept_fname)
print "1-based indexing if expected in the input peptide file, %s. Affected columns expected: %s"%(uniq_pept_fname,"peptide_start")
print "BEWARE: Enforcing 0-based notation back!"
pep_df["peptide_start"] = pep_df["peptide_start"] - 1


# connection between peptide info and spectrum info to be established ...
################################
spec_info['pept'] = spec_info['Peptide sequence'].str.upper()

# now let's check if all of the peptides from pep_df are in spec_info ...
spec_pept_present_index = pep_df['pept'].isin(spec_info['pept'])
print "FYI"
print " There are %d peptides are present in spectrum file out of %d in the peptide info file ..."%(sum(spec_pept_present_index),pep_df.shape[0])
print
spec_pept_uniq_num = spec_info['pept'].unique().size
peptdf_present_num = sum(pd.Series(spec_info['pept'].unique()).isin(pep_df['pept']))
print "And vice a versa"
print "There are %d spectrum peptides present in peptide file out of %d"%(peptdf_present_num,spec_pept_uniq_num)
print "Nevertheless, we move on with whatever we can do here ..."
#
#
# looking for Daemidated sites only ...
deamid = re.compile('[D,d]eamidat')
# enumerating all tractable peptides from pep_df ...
glyco_mod = []
for uniq_pept, pept_pos, prot_sr in pep_df[['pept','peptide_start','prot_seqrec']].itertuples(index=False):
    prot_seq = str(prot_sr)
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
    # looks like the inner loop here is the only place where we do need 0-based indexing switching ...
    for type_aa,gpos_pept,value in modifs:
        if (type_aa in ['n','N']) and (np.abs(value-3)<0.1):
            # 'gpos_pept' - is 1-based relative position of gsite_start_N in the peptide ...
            gsite_start = pept_pos + gpos_pept - 1 # zero-based coordinate ...
            gsite_stop  = pept_pos + gpos_pept - 1 + 3 - 1 # zero-based coordinate ...
            glyco_start.append(gsite_start)
            # Due to slicing rules, we need [start:stop+1], no have position 'stop' included ...
            glyco_sites.append(prot_seq[gsite_start:gsite_stop+1])
    ############################################################
    # gstart must be 1-based for output, so it's important to correct it right here ...
    glyco_mod_str = ','.join([ gsite+("(%d)"%(gstart+1)) for gsite,gstart in zip(glyco_sites,glyco_start)])
    glyco_mod.append(glyco_mod_str)
############################
#
#
#
############################
pep_df['gsites'] = glyco_mod
############################
print
print " We need to unroll glyco-sites and expand the DataFrame ..."
# we need to unroll those rows whith more than 1 items of gsites ...
pep_df_cols = ['all_uids', 'pept', 'peptide_start', 'prot_seqrec', 'protlen', 'uid_max', 'gsites','aa_before','aa_after',"prot_name", "uniq_pept_count", "pept_probab"]
pep_df_gsite_extracted = []
for all_uids,pept,pept_start,prot_sr,protlen,uid_max,gsites,aa_bef,aa_aft,prot_name,upept_count,pept_probab in pep_df[pep_df_cols].itertuples(index=False):
    # unrolling ...
    for gsite in gsites.split(','):
        pep_df_gsite_extracted.append((all_uids,pept,pept_start,prot_sr,protlen,uid_max,gsite,aa_bef,aa_aft,prot_name,upept_count,pept_probab))
########################
pep_df_ext = pd.DataFrame(pep_df_gsite_extracted,columns = pep_df_cols)
print "mission accomplished ..."
# so far we didn't use gsite index for adressing, so it was 1-based safely ...

print
print "Now we'd need to collapse those identical glyco sites that orginiate from overlaing peptides ..."
gsite_uid_combined = pep_df_ext['gsites'] +':'+ pep_df_ext['all_uids']
gsites_uniq = gsite_uid_combined.unique()
print " There are %d unique (glyco-site, protein id) pairs that includes empty sites, however."%gsites_uniq.size
# now collapsing on those uniq sites ...
final_dataframe = []
for site_uniq in gsites_uniq:
    # indexes of a given unque gsite ...
    sites_index = (gsite_uid_combined == site_uniq)
    pep_site = pep_df_ext[sites_index]
    #
    pepts = ','.join( '('+pep_site['aa_before']+')'+pep_site['pept']+'('+pep_site['aa_after']+')' )
    all_uids, = pep_site['all_uids'].unique()
    # peptide_start must be corrected back to 1-based, we won't use it as index anymore ...
    pepts_start = ','.join(str(pos_value+1) for pos_value in pep_site['peptide_start'])
    prot_seq = str(pep_site['prot_seqrec'].unique()[0])
    prot_len, = pep_site['protlen'].unique()
    uid_max, = pep_site['uid_max'].unique()
    #
    prot_name, = pep_site["prot_name"].unique()
    upept_count = pep_site["uniq_pept_count"].unique()[0] # temporary solution, due to inconsistencies ...
    pept_probab = ','.join(str(_) for _ in pep_site["pept_probab"])
    #
    # make sure there is just a single gsite ...
    pep_site_gsite_uniq, = pep_site['gsites'].unique()
    # appending to final dataframe ...
    final_dataframe.append((pep_site_gsite_uniq,pepts,pepts_start,all_uids,prot_seq,prot_len,uid_max, prot_name, upept_count, pept_probab))
# turn to dataframe ...
final_dataframe = pd.DataFrame(final_dataframe,columns=['gsites', 'pept', 'peptide_start','all_uids', 'prot_seq', 'protlen', 'uid_max','prot_name', 'uniq_pept_count', 'pept_probab'])



#
print 
print "Now we'd need to add theoretical glycosilation sites as a separate column ..."
g_site = re.compile(r'(?=(N[ACDEFGHIKLMNQRSTVWY][TS]))')


def get_theor_sites_number(prot_seq):
    # find all sites ...
    all_sites = [(site.start(),site.groups()[0]) for site in g_site.finditer(prot_seq)]
    N_sites = len(all_sites)
    return N_sites

def get_theor_sites(prot_seq):
    # BEWARE ZERO-BASED INDEXING TO 1-BASED INDEXING TRANSLATION ...
    # find all sites ...
    all_sites = [(site.start(),site.groups()[0]) for site in g_site.finditer(prot_seq)]
    # N_sites = len(all_sites)
    return ','.join( (gsite_seq+'(%d)'%(pos+1)) for pos,gsite_seq in all_sites) # pos+1 - indexing transition ...



# # find all contexts ...
# all_contexts = [seq[((pos-6)if pos>6 else 0):((pos+7)if pos<seq_len-7 else seq_len)] for pos,_ in all_sites ]
# nterm_cys_pos = [-1]*N_sites
# cterm_cys_pos = [-1]*N_sites
# for i,(pos,_) in zip(range(N_sites),all_sites):
#     # check including and excluding thing here ...
#     # N terminal cysteine ...
#     nterm_seq = seq[:pos]
#     res = nterm_cys.search(nterm_seq)
#     if res is not None:
#         nterm_cys_pos[i] = res.start()
#     # C terminal cysteine ...
#     cterm_seq = seq[(pos+3):]
#     res = cterm_cys.search(cterm_seq)
#     if res is not None:
#         cterm_cys_pos[i] = pos+res.end()-1
# # find QXT/S sites, according to the new Reid's idea, they must be depleated from Glycoproteome ...
print " 'gsites_predicted' column uses 1-based numbering. Enforced and checked."
final_dataframe['gsites_predicted'] = final_dataframe['prot_seq'].apply(get_theor_sites)
final_dataframe['gsites_predicted_number'] = final_dataframe['prot_seq'].apply(get_theor_sites_number)
# gsite_start must be 1-based ...
final_dataframe['gsite_start'] = final_dataframe['gsites'].apply( lambda x: int(x[3:].strip('()')) if x else None )
final_dataframe['gsites_AA1_N'] = final_dataframe['gsites'].apply( lambda x: x[0] if x else None )
final_dataframe['gsites_AA2_XbutP'] = final_dataframe['gsites'].apply( lambda x: x[1] if x else None )
final_dataframe['gsites_AA3_TS'] = final_dataframe['gsites'].apply( lambda x: x[2] if x else None )
# final_dataframe['gsites']NFT(1098)


final_dataframe.to_csv(out_fname,index=False)


# 1    NIS(126)                                     (R)ALSNISLR(F)
# 2    NTT(794)  (R)CVYEALCNTTSECPPPVITR(Q),(E)ALCNTTSECPPPVITR...
# 3   NDT(1048)  (R)EAESLQPMTVVGTDYVFHNDTK(V),(E)SLQPMTVVGTDYVF...
# 4    NET(732)                             (K)LSHDANETLPLHLYVK(Y)
# 5    NMS(527)                              (K)SCVAVTSAQPQNMSR(A)
# 6   NVT(1001)                               (R)SINVTGQGFSLIQR(A)
# 7   NFT(1098)  (R)TEAGAFEYVPDPTFENFTGGVK(Q),(R)TEAGAFEYVPDPTF...
# 8   NLT(1067)                  (K)VVFLSPAVPEEPEAYNLTVLIEMDGHR(L)











