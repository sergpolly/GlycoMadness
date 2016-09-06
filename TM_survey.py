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
#
import argparse

import warnings
from Bio import BiopythonWarning, BiopythonParserWarning


#############################################################################################################################################
def derive_cterm_topo(row):
    # # All topo locations found during the survey ...
    # 'Cytoplasmic'
    # 'Extracellular'
    # 'Lumenal'
    # 'Mitochondrial'
    # 'Nuclear'
    # 'Perinuclear'
    # 'Vacuolar'
    # our understanding of opposite and same side topology ...
    topo_opposite_dict = {'Cytoplasmic':'Extracellular','Extracellular':'Cytoplasmic','Lumenal':'Cytoplasmic'}
    ###################################################################
    tm_num = int(row['TM_num']) if pd.notnull(row['TM_num']) else None
    nterm_topo = row['N-term'] if pd.notnull(row['N-term']) else None
    if None in [tm_num,nterm_topo]:
        print "Cannot derive C-terminus topoly for %s"%row['uid_fetched']
        return None
    else:
        # even number of TM-spans yields the same side for N-term topology
        if tm_num%2 == 0:
            return nterm_topo
        # odd number of TM-spans yields the opposite side for N-term topology
        else:
            return topo_opposite_dict[nterm_topo] if (nterm_topo in topo_opposite_dict) else 'Opposite'
#############################################################################################################################################
def get_topo_descr(row):
    topo_short_dict = {'Cytoplasmic':'cyt','Extracellular':'lum','Lumenal':'lum'}
    ###################################################################
    tm_num = int(row['TM_num']) if pd.notnull(row['TM_num']) else None
    nterm_topo = row['N-term'] if pd.notnull(row['N-term']) else None
    cterm_topo = row['C-term-derived'] if pd.notnull(row['C-term-derived']) else None
    if None in [tm_num,nterm_topo,cterm_topo]:
        print "Cannot derive description for %s"%row['uid_fetched']
        return None
    else:
        nterm_topo_short = topo_short_dict[nterm_topo] if (nterm_topo in topo_short_dict) else nterm_topo[:3].lower() 
        cterm_topo_short = topo_short_dict[cterm_topo] if (cterm_topo in topo_short_dict) else cterm_topo[:3].lower() 
        return "%dTM N%s-C%s"%(tm_num,nterm_topo_short,cterm_topo_short)
#############################################################################################################################################
def compare_cterm_topo(row):
    similar_topos = ['Extracellular','Lumenal']
    # compare N-term topology from genbank and Derived,
    # assuming Lumenal and Extracellular are the same ...
    if pd.notnull(row['C-term']) and pd.notnull(row['C-term-derived']):
        if row['C-term']==row['C-term-derived']:
            return True
        elif (row['C-term'] in similar_topos) and (row['C-term-derived'] in similar_topos):
            return True
        else:
            return False
    else:
        return False
#############################################################################################################################################
def get_TM_boundaries(row):
    if pd.notnull(row['TM_locs']):
        gstart = row['gstart']
        # extract TM locs in the format appropriate for analysis ...
        str_to_span = lambda span: tuple( map(int,span.split('..')) ) 
        tm_locs = [ str_to_span(span) for span in row['TM_locs'].strip().split(',')]
        ##############
        first_tm_start = tm_locs[0][0]
        last_tm_end    = tm_locs[-1][1]
        if gstart < first_tm_start:
            return pd.Series({'N-TM boundary':None, 'C-TM boundary':first_tm_start})            
        if gstart > last_tm_end:
            return pd.Series({'N-TM boundary':last_tm_end, 'C-TM boundary':None})            
        for tm_span_pair in zip(tm_locs,tm_locs[1:]):
            prev_tm, next_tm = tm_span_pair
            prev_tm_end   = prev_tm[1]
            next_tm_start = next_tm[0]
            if prev_tm_end < gstart < next_tm_start:
                return pd.Series({'N-TM boundary':prev_tm_end, 'C-TM boundary':next_tm_start})
        # given nontrivial row['TM_locs'], gsite nearset TM should have been found!!!!!!!!
        print "Unable to locate nearest TM for gstart %d at %s,\n tm_locs: %s"%(gstart,row['uid_fetched'],row['TM_locs'])
        return pd.Series({'N-TM boundary':None, 'C-TM boundary':None})
    else:
        return pd.Series({'N-TM boundary':None, 'C-TM boundary':None})
#############################################################################################################################################


dest = "../ANALYSIS"
input_fname = "Memb Prot input file.xlsx"
tm_survey_fname = os.path.join(dest, input_fname)


gb_dest = ".."
gb_fname = "total_pulled_prots_GI.gb"
gb_fname = "total_pulled_prots_UID.gb"
gb_fname = os.path.join( gb_dest, gb_fname )



# (???)
# Reading genbank mindfully next ...
gbrecs = ms.genebank_fix_n_read(gb_fname,'id')
######################################
# (?????) assign some module internal stuff ...
ms.gbrecs = gbrecs


############################
# READING file containing GeneName(and/or locus) and FetchID association ...
print "Reading %s with the protein list to make TM-survey on ..."%tm_survey_fname
tm_survey_list = pd.read_excel(tm_survey_fname)



############################################
# COLUMN NAMES IN THE INPUT FILE ...
############################################
# 'locus'
# 'prot_name'
# 'uid_fetched'
# 'protein length'
# 'signal seq'
# 'TM span'
# 'gsite_seq'
# 'gstart'
############################################



############################################
# DESIRED OUTPUT ...
############################################
# 'N-term'
# 'C-term'
# 'TM_num'
# 'TM_locs'
# 'C-term-derived'
# 'C-term-match'
# 'Topology'
# 'N-TM boundary'
# 'C-TM boundary'
############################################




# extract topolgy from genbanks ...
tm_survey_list = tm_survey_list.merge( tm_survey_list['uid_fetched'].apply( ms.get_topo ) , left_index=True, right_index=True)
# extract all TM spans info from genbank ...
tm_survey_list = tm_survey_list.merge( tm_survey_list['uid_fetched'].apply( ms.get_all_tms ) , left_index=True, right_index=True)

tm_survey_list['C-term-derived'] = tm_survey_list.apply(derive_cterm_topo, axis=1)
tm_survey_list['C-term-match'] = tm_survey_list.apply(compare_cterm_topo, axis=1)
tm_survey_list['Topology'] = tm_survey_list.apply(get_topo_descr, axis=1)


tm_survey_list = tm_survey_list.merge( tm_survey_list.apply(get_TM_boundaries, axis=1), left_index=True, right_index=True)



# ['TM_locs', 'gstart', 'N-TM boundary', 'C-TM boundary']
# TODO...
# 1 Define bad entries ...
# BAD, are the ones where C-term-derived and N-term does not match AND/OR  tm-span is not found in a GENBANK ...
bad_index = (~tm_survey_list['C-term-match']) | (tm_survey_list['TM_num'].isnull())
good_index = (~bad_index)



# 2 output ...
with open(os.path.join(dest,'FINAL_TM_SURVEY_OUT.csv'),'w') as fp:
    # output formatted BAD stuff first ...
    tm_survey_list[bad_index].to_csv(fp,index=False)
    #
    fp.write("\n\nBAD and ambigous entries are over. Entries that worked out to follow:\n\n")
    #
    tm_survey_list[good_index].to_csv(fp,index=False)





































