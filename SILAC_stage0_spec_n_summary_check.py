import os
import sys
import argparse
import pandas as pd
import ms_module as ms
import time
############################
#
#
# do some arguments parsing to make the script looks civilized ...
parser = argparse.ArgumentParser()
parser.add_argument("-r","--raw_silac", help="speicfy raw data SILAC file name (with/without path)",required=True)
parser.add_argument("-q","--quant_silac", help="speicfy quantification SILAC file name (with/without path)",required=True)
parser.add_argument("--separator", help="speicfy separator type in the input data",default='tab')
parser.add_argument("--verbose", help="verbose output", action="store_true")
parser.add_argument("--prefix", help="specify common part of the path for peptide and quantification files")
# parser.add_argument("--pipeline", help="Act as part of the pipeline: print NCBI search terms for proteins", action="store_true")
args = parser.parse_args()

# print args
###############################################
if args.verbose:
    print "Verbose output is to follow ..."
    print
###############################################
if args.prefix is not None:
    raw_fname = os.path.join( args.prefix, args.raw_silac )
    quant_fname = os.path.join( args.prefix, args.quant_silac )
else:
    raw_fname = args.raw_silac
    quant_fname = args.quant_silac
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
# peptide summary and quantification file names must be specified as command line arguments beforehand ...
raw_info = pd.read_csv(raw_fname,sep=separator)
quant_info = pd.read_csv(quant_fname,sep=separator)
#####################################################################
# REQUIRED PRE PROCESSING OF THE PEP AND SPEC TABLES ...
# connection between peptide info and quantification info to be established ...
quant_info['pept'] = quant_info['Sequence'].str.upper()
# parse_pname = lambda pn: pd.Series(ms.parse_prot_name(pn,verbose=args.verbose))
# # add columns with the parsed information on the proteins ...
# quant_info = quant_info.merge(quant_info['Protein Name'].apply(parse_pname), left_index=True, right_index=True)
# raw_info = raw_info.merge(raw_info['Protein Name'].apply(parse_pname), left_index=True, right_index=True)
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
quant_uniq_pepts = quant_info['pept'].unique().size # # of unique peptides in the quant file ...
raw_uniq_pepts = raw_info['Peptide Sequence'].unique().size # # of unique peptides in the peptide summary file ...
quant_peps_in_raw = quant_info['pept'].isin(raw_info['Peptide Sequence']) # quantification pepts that are in peptide summary ...
quant_peps_NOT_in_raw = quant_info['pept'][~quant_peps_in_raw].unique().size
raw_peps_in_quant = raw_info['Peptide Sequence'].isin(quant_info['pept']) # quantification pepts that are in peptide summary ...
raw_peps_NOT_in_quant = raw_info['Peptide Sequence'][~raw_peps_in_quant].unique().size

print
print 'New style output: # of unique PEPTIDES in both files and their pairwise NOT_INs ...'
arr = {'quantification':[quant_uniq_pepts,quant_peps_NOT_in_raw],'raw_data':[raw_peps_NOT_in_quant,raw_uniq_pepts]}
print pd.DataFrame(arr,index=['quantification','raw_data']).to_string(columns=['quantification','raw_data'])
print


#########################################################
# let's generate some files if there is a miss ...
if quant_peps_NOT_in_raw:
    df_notin = quant_info[~quant_peps_in_raw]
    # let's see where, corresponding Protein names are going ...
    df_inother = raw_info[ raw_info['Protein Name'].isin(df_notin['Protein Name']) ]
    fname_notin = os.path.join( args.prefix,"Quant_pepts_NOT_IN_raw.csv")
    with open(fname_notin,'w') as fp:
        df_notin.to_csv(fp, index=False)
        fp.write('\n')
        fp.write('Where the Protein names corresponding to these peptides are ended up ...')
        fp.write('\n')
        df_inother.to_csv(fp,index=False)
if raw_peps_NOT_in_quant:
    df_notin = raw_info[~raw_peps_in_quant]
    # let's see where, corresponding Protein names are going ...
    df_inother = quant_info[ quant_info['Protein Name'].isin(df_notin['Protein Name']) ]
    fname_notin = os.path.join( args.prefix,  "Raw_pepts_NOT_IN_Quant.csv")
    with open(fname_notin,'w') as fp:
        df_notin.to_csv(fp, index=False)
        fp.write('\n')
        fp.write('Where the Protein names corresponding to these peptides are ended up ...')
        fp.write('\n')
        df_inother.to_csv(fp,index=False)
#########################################################


quant_uniq_prots = quant_info['Protein Name'].unique().size
raw_uniq_prots = raw_info['Protein Name'].unique().size
quant_prots_in_raw = quant_info['Protein Name'].isin(raw_info['Protein Name'])
quant_prots_NOT_in_raw = quant_info['Protein Name'][~quant_prots_in_raw].unique().size
raw_prots_in_quant = raw_info['Protein Name'].isin(quant_info['Protein Name'])
raw_prots_NOT_in_quant = raw_info['Protein Name'][~raw_prots_in_quant].unique().size

print
print 'New style output: # of unique PROTEINS in both files and their pairwise NOT_INs ...'
arr = {'quantification':[quant_uniq_prots,quant_prots_NOT_in_raw],'raw_data':[raw_prots_NOT_in_quant,raw_uniq_prots]}
print pd.DataFrame(arr,index=['quantification','raw_data']).to_string(columns=['quantification','raw_data'])
print


#########################################################
# let's generate some files if there is a miss ...
if quant_prots_NOT_in_raw:
    df_notin = quant_info[~quant_prots_in_raw]
    # let's see where, corresponding Peptides are going ...
    df_inother = raw_info[ raw_info['Peptide Sequence'].isin(df_notin['pept']) ]
    fname_notin = os.path.join( args.prefix,"Quant_Proteins_NOT_IN_raw.csv")
    with open(fname_notin,'w') as fp:
        df_notin.to_csv(fp, index=False)
        fp.write('\n')
        fp.write('Where the Peptides corresponding to these Protein names are ended up ...')
        fp.write('\n')
        df_inother.to_csv(fp,index=False)
if raw_prots_NOT_in_quant:
    df_notin = raw_info[~raw_prots_in_quant]
    # let's see where, corresponding Peptides are going ...
    df_inother = quant_info[ quant_info['pept'].isin(df_notin['Peptide Sequence']) ]
    fname_notin = os.path.join( args.prefix, "Raw_Proteins_NOT_IN_Quant.csv")
    with open(fname_notin,'w') as fp:
        df_notin.to_csv(fp, index=False)
        fp.write('\n')
        fp.write('Where the Peptides corresponding to these Protein names are ended up ...')
        fp.write('\n')
        df_inother.to_csv(fp,index=False)
#########################################################


# just so one can stare at the screen for half a minute ...
time.sleep(20.0)


# # quant_peps_in_raw = quant_info['pept'].isin(raw_info['Peptide Sequence']) # quantification pepts that are in peptide summary ...
# # quant_uniq_pepts = quant_info['pept'].unique().size # # of unique peptides in the quantification file ...
# # raw_uniq_pepts = raw_info['Peptide Sequence'].unique().size # # of unique peptides in the peptide summary file ...
# if ( quant_peps_in_raw.all() and (quant_uniq_pepts==raw_uniq_pepts) ):
#     quite_msg += "Peptides: OK\n"
#     verbose_msg += "All peptides from quantification file are present in the peptide summary file!\n"
#     verbose_msg += "Sounds very logicall.\n"
#     verbose_msg += "There are %d of such peptides.\n\n"%raw_uniq_pepts
#     success_status = success_status and True
# else:
#     quite_msg += "Peptides: not OK!\n"
#     verbose_msg += "There are some peptides from quantification file that are not present in the peptide summary file:\n"
#     verbose_msg += str(quant_info[~quant_peps_in_raw][['Protein Name','pept']])+"\n"
#     verbose_msg += "There are %d unique peptides in the summary\n"%raw_uniq_pepts
#     verbose_msg += "and %d peptides in quantification file."%quant_uniq_pepts
#     verbose_msg += "\n"
#     verbose_msg += """It is a strange situation,
#         suggesting that different parameters were used in Scaffold
#         to generate peptide summary and quantification file.\n"""
#     verbose_msg += """We proceed dismissing this fact,
#         and using all the gsite/peptides pairs present in quantification,
#         thus assuming self-sufficiency of the quantification file
#         and its prevalence over peptide summary.
#         In other words, there is nothing in the peptide summary file,
#         that cannot be deduced from the quantification file.(? seem to be true, but is it general?)\n\n"""
#     success_status = success_status and False
# ##################################################################################################
# quant_prots_in_raw = quant_info['Protein Name'].isin(raw_info['Protein Name'])
# quant_uniq_prots = quant_info['Protein Name'].unique().size
# raw_uniq_prots = raw_info['Protein Name'].unique().size
# if ( quant_prots_in_raw.all() and (quant_uniq_prots==raw_uniq_prots) ):
#     quite_msg += "Proteins: OK\n"
#     verbose_msg += "Peptide summary and quantification files are refferring to the same %d protein names.\n"%quant_uniq_prots
#     verbose_msg += "It is a good sign!\n\n"
#     success_status = success_status and True
# else:
#     quite_msg += "Proteins: not OK!\n"
#     verbose_msg += "Pept.summary file is refferring to %d unique protein names.\n"%quant_uniq_prots
#     verbose_msg += "quantification file is refferring to %d unique protein names.\n"%raw_uniq_prots
#     verbose_msg += str(quant_info[~quant_prots_in_raw][['Protein Name','pept']])+"\n"
#     verbose_msg += "This is unexpected discrepancy: proceed using data stored in the quantification file.\n\n"
#     success_status = success_status and False

# # print that summary message finally ...
# if args.verbose:
#     print verbose_msg
# else:
#     print quite_msg

# # exit with the proper status code ...
# sys.exit(0 if success_status else 1)

















