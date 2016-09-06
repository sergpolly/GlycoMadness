import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import pandas as pd
import numpy as np
import ms_module as ms
############################
from StringIO import StringIO
import warnings
from Bio import BiopythonWarning, BiopythonParserWarning
import subprocess as sub



dest = "../PULLED_PROTEINS_TOTAL"

# get file names of all the pulled files from destination ...
pulled_files = sub.check_output(['ls',dest])
pulled_files = pulled_files.strip().split('\n')
pulled_files = ['/'.join([dest,fname]) for fname in pulled_files]

# # # Reading genbank mindfully next ...
gbrecs_id = [ms.genebank_fix_n_read(fn,'id') for fn in pulled_files]
gbrecs_gi = [ms.genebank_fix_n_read(fn,'gi') for fn in pulled_files]
# # # ######################################
# # # # assign some module internal stuff ...
# # # ms.gbrecs = gbrecs


# perform some simple and stupid tests ...
id_keys = [np.asarray(gg.keys()) for gg in gbrecs_id]
gi_keys = [np.asarray(gg.keys()) for gg in gbrecs_gi]
#####################################################################################
print "Proteins in each of id groups are unique, True or False:"
answ1 = np.asarray([f.size==np.unique(f).size for f in id_keys]).all()
print answ1
#####################################################################################
print "Proteins in each of gi groups are unique, True or False:"
answ2 = np.asarray([f.size==np.unique(f).size for f in gi_keys]).all()
print answ2
#####################################################################################
print "Number of proteins in corresponding id and gi groups matches, True or False:"
answ3 = np.asarray([f.size==s.size for f,s in zip(id_keys,gi_keys)]).all()
print answ3
#####################################################################################
print 
print "All preceding answers must be True, otherwise, there is something wrong with the dataset ..."
if (answ3 and answ2 and answ1):
    print "Looks like we're all set!!!"
else:
    print "Check what's going on?!?!"
#####################################################################################



# uniq IDs and GIs ...
uniq_ids = np.unique(np.concatenate(id_keys))
uniq_gis = np.unique(np.concatenate(gi_keys))

print "Number of uniq IDs and GIs, %d and %d are supposed to match..."%(uniq_ids.size,uniq_gis.size)


################################################################################################
# now we'd need to look up all Genebank records
# from both lists 'uniq_ids' 'uniq_gis' and store them as separate files ...
# TO BE CONTINUED ....
################################################################################################

# gbrecs_id
# gbrecs_gi


# compose an array of unique genbank protein-records based on the uniq_ids list ...
uniq_gbrecs_id = {}
for uid in uniq_ids:
    for gbrec in gbrecs_id:
        if uid in gbrec:
            uniq_gbrecs_id[uid] = gbrec[uid]
            break

with open('/'.join(['..','total_pulled_prots_UID.gb']),'w') as fp:
    SeqIO.write(uniq_gbrecs_id.values(),fp,'genbank')


# compose an array of unique genbank protein-records based on the uniq_gis list ...
uniq_gbrecs_gi = {}
for gi in uniq_gis:
    for gbrec in gbrecs_gi:
        if gi in gbrec:
            uniq_gbrecs_gi[gi] = gbrec[gi]
            break


with open('/'.join(['..','total_pulled_prots_GI.gb']),'w') as fp:
    SeqIO.write(uniq_gbrecs_gi.values(),fp,'genbank')



# within the context of all the ground work and based on the 'total_pulled_prots_UID.gb' and 'total_pulled_prots_GI.gb' sizes, 
# we can safely assume, that these files store identical list of protein records. Order is not guaranteed though ...


# THE END ...



















































