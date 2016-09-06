import os
import sys
import argparse
import pandas as pd
import ms_module as ms
############################
#
#
# do some arguments parsing to make the script looks civilized ...
parser = argparse.ArgumentParser()
parser.add_argument("-p","--pept_summary", help="speicfy peptide summary file name (with/without path)",required=True)
parser.add_argument("-s","--spec_summary", help="speicfy spectrum file name (with/without path)",required=True)
parser.add_argument("--separator", help="speicfy separator type in the input data",default='tab')
parser.add_argument("--verbose", help="verbose output", action="store_true")
parser.add_argument("--prefix", help="specify common part of the path for peptide and spectrum files")
# parser.add_argument("--pipeline", help="Act as part of the pipeline: print NCBI search terms for proteins", action="store_true")
args = parser.parse_args()

# print args
###############################################
if args.verbose:
    print "Verbose output is to follow ..."
    print
###############################################
if args.prefix is not None:
    pep_fname = os.path.join( args.prefix, args.pept_summary )
    spec_fname = os.path.join( args.prefix, args.spec_summary )
else:
    pep_fname = args.pept_summary
    spec_fname = args.spec_summary
#
#
#
#
if args.separator == "tab":
    separator = '\t'
elif args.separator == "comma":
    separator = ','
else:
    separator = '\t'
#
# peptide summary and spectrum file names must be specified as command line arguments beforehand ...
pep_info = pd.read_csv(pep_fname,sep=separator)
spec_info = pd.read_csv(spec_fname,sep=separator)
#####################################################################
# REQUIRED PRE PROCESSING OF THE PEP AND SPEC TABLES ...
# connection between peptide info and spectrum info to be established ...
spec_info['pept'] = spec_info['Peptide sequence'].str.upper()
# parse_pname = lambda pn: pd.Series(ms.parse_prot_name(pn,verbose=args.verbose))
# # add columns with the parsed information on the proteins ...
# spec_info = spec_info.merge(spec_info['Protein name'].apply(parse_pname), left_index=True, right_index=True)
# pep_info = pep_info.merge(pep_info['Protein name'].apply(parse_pname), left_index=True, right_index=True)
# #####################################################################
#
#
verbose_msg = "\n\n"
quite_msg = "\n\n"
success_status = True
#
#
# SOME VALIDATIONS ...


###############################
# nice form of validation in a convenient table format ...
spec_uniq_pepts = spec_info['pept'].unique().size # # of unique peptides in the spectrum file ...
pep_uniq_pepts = pep_info['Peptide sequence'].unique().size # # of unique peptides in the peptide summary file ...
spec_peps_in_peps = spec_info['pept'].isin(pep_info['Peptide sequence']) # spectrum pepts that are in peptide summary ...
spec_peps_NOT_in_peps = spec_info['pept'][~spec_peps_in_peps].unique().size
pep_peps_in_spec = pep_info['Peptide sequence'].isin(spec_info['pept']) # spectrum pepts that are in peptide summary ...
pep_peps_NOT_in_spec = pep_info['Peptide sequence'][~pep_peps_in_spec].unique().size

print
print 'New style output: # of unique PEPTIDES in both files and their pairwise NOT_INs ...'
arr = {'spectrum':[spec_uniq_pepts,spec_peps_NOT_in_peps],'pept_sum':[pep_peps_NOT_in_spec,pep_uniq_pepts]}
print pd.DataFrame(arr,index=['spectrum','pept_sum']).to_string(columns=['spectrum','pept_sum'])
print


#########################################################
# let's generate some files if there is a miss ...
if spec_peps_NOT_in_peps:
    df_notin = spec_info[~spec_peps_in_peps]
    # let's see where, corresponding Protein names are going ...
    df_inother = pep_info[ pep_info['Protein name'].isin(df_notin['Protein name']) ]
    fname_notin = os.path.join( args.prefix,"Spectrum_pepts_NOT_IN_pept_summary.csv")
    with open(fname_notin,'w') as fp:
        df_notin.to_csv(fp, index=False)
        fp.write('\n')
        fp.write('Where the Protein names corresponding to these peptides are ended up ...')
        fp.write('\n')
        df_inother.to_csv(fp,index=False)
if pep_peps_NOT_in_spec:
    df_notin = pep_info[~pep_peps_in_spec]
    # let's see where, corresponding Protein names are going ...
    df_inother = spec_info[ spec_info['Protein name'].isin(df_notin['Protein name']) ]
    fname_notin = os.path.join( args.prefix,  "Pept_summary_pepts_NOT_IN_Spectrum.csv")
    with open(fname_notin,'w') as fp:
        df_notin.to_csv(fp, index=False)
        fp.write('\n')
        fp.write('Where the Protein names corresponding to these peptides are ended up ...')
        fp.write('\n')
        df_inother.to_csv(fp,index=False)
#########################################################

#THERE IS A MISTAKE ... pepts correspond to spec  and vice versa ...
pept_uniq_prots = spec_info['Protein name'].unique().size
spec_uniq_prots = pep_info['Protein name'].unique().size
spec_prots_in_peps = spec_info['Protein name'].isin(pep_info['Protein name'])
spec_prots_NOT_in_peps = spec_info['Protein name'][~spec_prots_in_peps].unique().size
pep_prots_in_spec = pep_info['Protein name'].isin(spec_info['Protein name'])
pep_prots_NOT_in_spec = pep_info['Protein name'][~pep_prots_in_spec].unique().size

print
print 'New style output: # of unique PROTEINS in both files and their pairwise NOT_INs ...'
arr = {'spectrum':[spec_uniq_prots,spec_prots_NOT_in_peps],'pept_sum':[pep_prots_NOT_in_spec,pept_uniq_prots]}
print pd.DataFrame(arr,index=['spectrum','pept_sum']).to_string(columns=['spectrum','pept_sum'])
print


#########################################################
# let's generate some files if there is a miss ...
if spec_prots_NOT_in_peps:
    df_notin = spec_info[~spec_prots_in_peps]
    # let's see where, corresponding Peptides are going ...
    df_inother = pep_info[ pep_info['Peptide sequence'].isin(df_notin['pept']) ]
    fname_notin = os.path.join( args.prefix,"Spectrum_Proteins_NOT_IN_pept_summary.csv")
    with open(fname_notin,'w') as fp:
        df_notin.to_csv(fp, index=False)
        fp.write('\n')
        fp.write('Where the Peptides corresponding to these Protein names are ended up ...')
        fp.write('\n')
        df_inother.to_csv(fp,index=False)
if pep_prots_NOT_in_spec:
    df_notin = pep_info[~pep_prots_in_spec]
    # let's see where, corresponding Peptides are going ...
    df_inother = spec_info[ spec_info['pept'].isin(df_notin['Peptide sequence']) ]
    fname_notin = os.path.join( args.prefix, "Pept_summary_Proteins_NOT_IN_Spectrum.csv")
    with open(fname_notin,'w') as fp:
        df_notin.to_csv(fp, index=False)
        fp.write('\n')
        fp.write('Where the Peptides corresponding to these Protein names are ended up ...')
        fp.write('\n')
        df_inother.to_csv(fp,index=False)
#########################################################






# spec_peps_in_peps = spec_info['pept'].isin(pep_info['Peptide sequence']) # spectrum pepts that are in peptide summary ...
# spec_uniq_pepts = spec_info['pept'].unique().size # # of unique peptides in the spectrum file ...
# pep_uniq_pepts = pep_info['Peptide sequence'].unique().size # # of unique peptides in the peptide summary file ...
if ( spec_peps_in_peps.all() and (spec_uniq_pepts==pep_uniq_pepts) ):
    quite_msg += "Peptides: OK\n"
    verbose_msg += "All peptides from spectrum file are present in the peptide summary file!\n"
    verbose_msg += "Sounds very logicall.\n"
    verbose_msg += "There are %d of such peptides.\n\n"%pep_uniq_pepts
    success_status = success_status and True
else:
    quite_msg += "Peptides: not OK!\n"
    verbose_msg += "There are some peptides from spectrum file that are not present in the peptide summary file:\n"
    verbose_msg += str(spec_info[~spec_peps_in_peps][['Protein name','pept']])+"\n"
    verbose_msg += "There are %d unique peptides in the summary\n"%pep_uniq_pepts
    verbose_msg += "and %d peptides in spectrum file."%spec_uniq_pepts
    verbose_msg += "\n"
    verbose_msg += """It is a strange situation,
        suggesting that different parameters were used in Scaffold
        to generate peptide summary and spectrum file.\n"""
    verbose_msg += """We proceed dismissing this fact,
        and using all the gsite/peptides pairs present in spectrum,
        thus assuming self-sufficiency of the spectrum file
        and its prevalence over peptide summary.
        In other words, there is nothing in the peptide summary file,
        that cannot be deduced from the spectrum file.(? seem to be true, but is it general?)\n\n"""
    success_status = success_status and False
##################################################################################################
spec_prots_in_peps = spec_info['Protein name'].isin(pep_info['Protein name'])
pept_uniq_prots = spec_info['Protein name'].unique().size
spec_uniq_prots = pep_info['Protein name'].unique().size
if ( spec_prots_in_peps.all() and (pept_uniq_prots==spec_uniq_prots) ):
    quite_msg += "Proteins: OK\n"
    verbose_msg += "Peptide summary and spectrum files are refferring to the same %d protein names.\n"%pept_uniq_prots
    verbose_msg += "It is a good sign!\n\n"
    success_status = success_status and True
else:
    quite_msg += "Proteins: not OK!\n"
    verbose_msg += "Pept.summary file is refferring to %d unique protein names.\n"%pept_uniq_prots
    verbose_msg += "Spectrum file is refferring to %d unique protein names.\n"%spec_uniq_prots
    verbose_msg += str(spec_info[~spec_prots_in_peps][['Protein name','pept']])+"\n"
    verbose_msg += "This is unexpected discrepancy: proceed using data stored in the spectrum file.\n\n"
    success_status = success_status and False

# print that summary message finally ...
if args.verbose:
    print verbose_msg
else:
    print quite_msg


# just so one can stare at the screen for half a minute ...
time.sleep(20.0)

# exit with the proper status code ...
sys.exit(0 if success_status else 1)

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
# H prot_name, parsed out human-readable name from 'Protein name'
# H1 gene_name, parsed out GN=xxx from 'Protein name'
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





















