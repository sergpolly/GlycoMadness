import sys as __sys
import requests as __requests
# from Bio import Seq as Seq
from Bio import SeqIO as __SeqIO
from Bio import SeqRecord as __SeqRecord
import StringIO as __StringIO
import re as __re
import pandas as __pd
import numpy  as __np

__web_request_status_collection = {200: "The request was processed successfully.",
400: "Bad request. There is a problem with your input.",
404: "Not found. The resource you requested doesnt exist.",
410: "Gone. The resource you requested was removed.",
500: "Internal server error. Most likely a temporary problem, but if the problem persists please contact us.",
503: "Service not available. The server is being updated, try again later."}


# # typical Uniprot ID (for isoform 2) ...
# uid = "P16870-2"


def get_uniprot(session,uid,seq_format='fasta'):
    """ see http://www.uniprot.org/help/programmatic_access for details """
    # treat missing or inknown data fairly ...
    if uid is None:
        return None
    # the way we form Uniprot ID request ...
    get_uid_url = lambda _: "http://www.uniprot.org/uniprot/%s.fasta"%_
    # web request for a given uid ...
    uid_url = get_uid_url(uid)
    # make a request ...
    req_res = session.get(uid_url)
    # check request status ...
    if req_res.status_code==200 and bool(req_res.content):
        # to be read by __SeqIO ...
        string_as_handle = __StringIO.StringIO(req_res.content)
        seq_rec = __SeqIO.read(string_as_handle,seq_format)
        return seq_rec
    elif req_res.status_code==200:
        print __web_request_status_collection[req_res.status_code]
        print "... But, the content is empty for accession number %s!"%uid
        __sys.exit(1)
    elif req_res.status_code in __web_request_status_collection:
        print __web_request_status_collection[req_res.status_code]
        __sys.exit(1)
    else:
        print "Unknown status code returned!"
        __sys.exit(1)




def stupid_aligner(peptide,protein):
    """ 1-based number of peptide occurance in the protein ..."""
    # small subfunction to get the hamming distance ...
    def hamming_distance(seq1,seq2):
        assert len(seq1)==len(seq2)
        mismatches = sum( (l1!=l2) for l1,l2 in zip(seq1,seq2) )
        return mismatches
    # treat missing data fairly ...
    if protein is None:
        return None
    # sequences to strings, making sure they are SeqRec or Seq entyties ...
    peptide_str = str(peptide.seq) if type(peptide)==__SeqRecord.SeqRecord else str(peptide)
    protein_str = str(protein.seq) if type(protein)==__SeqRecord.SeqRecord else str(protein)
    # lengths ...
    pept_len = len(peptide_str)
    prot_len = len(protein_str)
    # stupid alignment ...
    min_mismatch = prot_len
    align_frame = 0
    for f in range(prot_len - pept_len + 1):
        prot_substring = protein_str[f:f+pept_len]
        delta_hd = hamming_distance(peptide_str,prot_substring)
        # in case perfect alignment is found ...
        if delta_hd == 0:
            align_frame = f
            return align_frame
        # or keep searching minimum mismatch alignment ...
        if delta_hd < min_mismatch:
            align_frame, min_mismatch = f, delta_hd
    # make a verbose report after the whole protein was scanned ...
    print "Beware! Best alignment found has %d mismatches for peptide %s"%(min_mismatch,peptide_str)
    # Here we're enforcing the 1-based indexing ...
    return align_frame + 1





# modifications can come in various forms ...
# n4: Deamidated:18O(1) (+2.99)
# Deamidated:18O(1) (+3)
def parse_spectrum_modifications(modifier):
    # print modifier
    loc_mod = modifier.strip()
    loc_mod = loc_mod.split(' ')
    if len(loc_mod)==3:
        mod_type = loc_mod[0].strip(':')
        # aa modified and peptide position ...
        # POSITIONS ARE MOST LIKELY TO HAVE 1-BASED INDEX ...
        mod_type_aa, mod_type_pos = mod_type[0], int(mod_type[1:])
        # value ...
        mod_val = float(loc_mod[2].strip('()'))
        #
        return (mod_type_aa, mod_type_pos, mod_val)
    elif len(loc_mod)==2:
        mod_val = float(loc_mod[1].strip('()'))
        #
        return (None, None, mod_val)
    else:
        print "Unknown format of modification description: %s"%modifier
        sys.exit(1)



# protein name example:
# sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha chain OS=Homo sapiens GN=HLA-A PE=1 SV=2
# regexp to use often:
# first part of full protein name with Uid and Locus ...
__left_part = __re.compile("[a-z]{2}\|[A-Z0-9\-\_\.\s]+\|[A-Z0-9\_]+")
# features like OS,GN etc, simply to extract what fatures are present ...
__feature_types = __re.compile("[A-Z]{2}=")
def parse_prot_name(prot_name,verbose=True):
    """Functions parses provided protein name and returns a dict with: uid,locus,prot_name,GeneName,OrgSource
    -------------------------
    prot name expected format:
    'sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha chain OS=Homo sapiens GN=HLA-A PE=1 SV=2'"""
    #
    dict_to_return = {}
    #
    # GO THROUGH FEATURES, like OS=Homo sapiens GN=HLA-A PE=1 SV=2 ...
    # what features are present ...
    f_types = [ f.strip('=') for f in  __feature_types.findall(prot_name) ]
    if f_types:
        # generate regexp based on the combination of features:
        # for example: "GN=(.+)PE=(.+)SV=(.+)"
        f_pattern = ''.join(['%s=(.+)'%f for f in f_types])
        # find the whole pattern in the protein name:
        f_values, = __re.findall(f_pattern,prot_name)
        # right features part for stripping ...
        rp_extracted = ''.join( "%s=%s"%(k,v) for k,v in zip(f_types,f_values) )
    else:
        rp_extracted = ''
    # store everything in an f_dict ...
    f_dict = dict( zip(f_types,f_values) ) if f_types else {}
    #
    #
    # extract left most part (if it's extractable)...:
    lp_extracted = __left_part.findall(prot_name)
    if len(lp_extracted)!=1:
        if verbose:
            print "Could not match left part of protein name: %s"%prot_name
            print "Extraction result is: ", lp_extracted
        lp_extracted = ""
    else:
        lp_extracted, = lp_extracted
        _,uid,locus = lp_extracted.split('|')
        dict_to_return['uid'] = uid.strip()
        dict_to_return['locus'] = locus.strip()
    #
    # strip left and right part of the full prot name, to get the human readable
    prot_name_extracted = prot_name.replace(lp_extracted,'').replace(rp_extracted,'').strip()
    dict_to_return['prot_name'] = prot_name_extracted.strip()
    #
    # returning all extracted information ...
    if ('GN' not in f_types)or('OS' not in f_types) :
        if verbose:
            print "There is no GeneName or OrganismSource in the protein name: %s"%prot_name
            print "Feature types extracted are: ",f_types
    else:
        dict_to_return['GN'] = f_dict['GN'].strip()
        dict_to_return['OS'] = f_dict['OS'].strip()
    # returning the result regardless of the internals ...
    return dict_to_return





########################################################################################################################
# READING GENEBANK CAREFULLY ...
#########################################################
import warnings as __warnings 
from Bio import BiopythonWarning as __BiopythonWarning
from Bio import BiopythonParserWarning as __BiopythonParserWarning


def __fix_str(input_str,known_errors,known_fixes):
    """function that would be replacing known error strings in the input by the known fix replacement"""
    fixed_input = str(input_str)
    for err,fix in zip(known_errors,known_fixes):
        fixed_input = fixed_input.replace(err,fix)
    return fixed_input
########################################################################################################################
def __fix_genbank(error_msg,genbank_fname):
    error_pattern = "(\d+\^\d+)"
    known_errors = __re.findall(error_pattern,str(error_msg))
    known_fixes  = [err.replace('^','..') for err in known_errors]
    #
    if known_errors:
        with open(genbank_fname,'r') as fp:
            file_string_io = [ __fix_str(line,known_errors,known_fixes) for line in fp ]
        file_string_io = __StringIO.StringIO(''.join(file_string_io))
        print "Error if fixed locally, proceeding ..."
        return file_string_io
    else:
        print "Error in the genbank could not be resolved. Termination"
        __sys.exit(1)
########################################################################################################################
# READING SEQUENCES FROM THE FILE ...
def genebank_fix_n_read(gb_fname,key_func_type='gi'):
    """locations formatted as  '1^593' cause BioPython error while reading genbanks ...
    We are addressing that by fixing genbank source on the fly ..."""
    print "Reading %s with genebank records from the NCBI fetch ..."%gb_fname
    # choose the key-function based on the 'key_func_type' argument:
    if key_func_type == 'gi':
        key_function=lambda rec: rec.annotations['gi']
    if key_func_type == 'id':
        key_function=lambda rec: rec.id
    print "Using %s as a key."%key_func_type
    #
    with __warnings.catch_warnings():
        # e.g. BiopythonParserWarning: Dropping bond qualifier in feature location
        #
        __warnings.simplefilter("ignore", __BiopythonParserWarning)
        #
        gb_recs_iter = __SeqIO.parse(gb_fname,'gb')
        try:
            gbrecs = __SeqIO.to_dict( gb_recs_iter, key_function=key_function )
        except ValueError, er_msg:
            print "Catched ValueError: %s"%str(er_msg)
            # #
            # Invalid between location '1^593'
            # Try to fix that thing, by replacing '^' with '..'
            file_string_io = __fix_genbank(er_msg,gb_fname)
            gb_recs_iter = __SeqIO.parse(file_string_io,'gb')
            gbrecs = __SeqIO.to_dict( gb_recs_iter, key_function=key_function )
    return gbrecs########################################################################################################################
# READING GENEBANK CAREFULLY ...
#########################################################






###########################################################################
#  STAGE 2 FUNCTIONS AND METHODS ...
###########################################################################

gbrecs = None

########################################################################################################################
def get_enzyme(sample_cat):
    if __pd.notnull(sample_cat):
        if 'Try' in sample_cat:
            return 'T'
        elif 'Glu' in sample_cat:
            return 'G'
        else:
            return None
    else:
        return None
######################################################################################################
__yn_map = {True:'Y',False:'N'}
######################################################################################################
def get_tm_span(fidx):
    if __pd.notnull(fidx):
        try:
            fidx_str = str(int(fidx))
        except:
            fidx_str = str(fidx)
        feats_descr = []
        for feat in gbrecs[fidx_str].features:
            quals = feat.qualifiers
            feats_descr.append( __yn_map["Transmembrane" in ''.join(quals['region_name'])] if ('region_name' in quals) else None )
        return 'Y' if ('Y' in feats_descr) else 'N'
    else:
        return None
######################################################################################################
def get_genename(fidx):
    if __pd.notnull(fidx):
        try:
            fidx_str = str(int(fidx))
        except:
            fidx_str = str(fidx)
        fid_features = gbrecs[fidx_str].features
        if 'gene' in fid_features[1].qualifiers:
            return fid_features[1].qualifiers['gene'][0]
        else:
            for feat in fid_features:
                if 'gene' in feat.qualifiers:
                    return feat.qualifiers['gene'][0]
        # if GN wasn't found in any of the features, return None ...
        return None
    else:
        return None
######################################################################################################
def get_signal(fidx):
    if __pd.notnull(fidx):
        try:
            fidx_str = str(int(fidx))
        except:
            fidx_str = str(fidx)
        feats_descr = []
        for feat in gbrecs[fidx_str].features:
            quals = feat.qualifiers
            feats_descr.append( __yn_map['Signal' in quals['region_name']] if ('region_name' in quals) else None )
        return 'Y' if ('Y' in feats_descr) else 'N'
    else:
        return None
######################################################################################################
# to be edited to trun into feature locator ...
def get_signal_loc(fidx):
    if __pd.notnull(fidx):
        try:
            fidx_str = str(int(fidx))
        except:
            fidx_str = str(fidx)
        for feat in gbrecs[fidx_str].features:
            quals = feat.qualifiers
            if ('region_name' in quals):
                if 'Signal' in quals['region_name']:
                    # start,end = (feat.location.start.position+1, feat.location.end.position)
                    return "%d..%d"%(feat.location.start.position+1, feat.location.end.position)
        return None
    else:
        return None
######################################################################################################
# NEW TO BE TESTED ...
def get_topo(fidx):
    def extract_topo_info(quals):
        if 'note' in quals:
            # Cytoplasmic or Extracellular
            if "Cytoplasmic" in ' '.join(quals['note']):
                return 'Cytoplasmic'
            elif "Extracellular" in ' '.join(quals['note']):
                return 'Extracellular'
            elif "Lumenal" in ' '.join(quals['note']):
                return 'Lumenal'
            elif "Mitochondrial" in ' '.join(quals['note']):
                return 'Mitochondrial'
            elif "Nuclear" in ' '.join(quals['note']):
                return 'Nuclear'
            elif "Perinuclear" in ' '.join(quals['note']):
                return 'Perinuclear'
            elif "Vacuolar" in ' '.join(quals['note']):
                return 'Vacuolar'
            else:
                print "Unidentified localization of topo domain %s"%str(quals)
                return None
        else:
            print "topo domain has no note, quals: %s"%str(quals)
            return None
    # seems to be working ok ...
    if __pd.notnull(fidx):
        try:
            fidx_str = str(int(fidx))
        except:
            fidx_str = str(fidx)
        # at first, simply survey all Regions, and record all topo-domains ...
        topo_domains = {}
        for feat_id,feat in enumerate(gbrecs[fidx_str].features):
            quals = feat.qualifiers
            if 'region_name' in quals:
                if "Topological domain" in ''.join(quals['region_name']):
                    start = feat.location.start.position+1
                    end = feat.location.end.position
                    topo_domains[(start,end)] = feat_id
        #################################################################
        if not topo_domains:
            print "No Topological domains detcted for %s"%fidx_str
            return __pd.Series({'N-term':None,'C-term':None})
        else:
            topo_domains_description = {}
            N_key = min(topo_domains)
            C_key = max(topo_domains)
            N_domain_info = gbrecs[fidx_str].features[ topo_domains[N_key] ].qualifiers
            C_domain_info = gbrecs[fidx_str].features[ topo_domains[C_key] ].qualifiers
            topo_domains_description['N-term'] = extract_topo_info(N_domain_info)
            topo_domains_description['C-term'] = extract_topo_info(C_domain_info)
            #####################################################
            return __pd.Series(topo_domains_description)
    else:
        return None
######################################################################################################
def get_all_tms(fidx):
    #  seems to be working fine ...
    if __pd.notnull(fidx):
        try:
            fidx_str = str(int(fidx))
        except:
            fidx_str = str(fidx)
        tm_locs = []
        for feat in gbrecs[fidx_str].features:
            quals = feat.qualifiers
            if 'region_name' in quals:
                if "Transmembrane" in ''.join(quals['region_name']):
                    start = feat.location.start.position+1
                    end = feat.location.end.position
                    tm_locs.append("%d..%d"%(start,end))
        return __pd.Series({'TM_num':len(tm_locs),'TM_locs':','.join(tm_locs)}) if tm_locs else __pd.Series({'TM_num':None,'TM_locs':None})
    else:
        return __pd.Series({'TM_num':None,'TM_locs':None})
######################################################################################################
# THAT WHAT IS REMAINING ...
# 5.  Closest TM boundaries to g-site. Here, I want the closest TM boundary on both sides (if present).
# For example, entry 5 on the output example, has a g site at 212, and TM spans at 133..153 and 225..245  that flank this site.
# 153 is the closest N-terminal boundary, 225 is the closest C-terminal boundary.
# For proteins with a single TM span (lines 2 and 3 of sample output file, one of these output columns will be empty.
######################################################################################################
#######################################################################################################
########################################################################################################################
def pept_isin(row):
    pept,fetchid = row
    if __pd.notnull(fetchid):
        try:
            fidx_str = str(int(fetchid))
        except:
            fidx_str = str(fetchid)
        prot_seq = gbrecs[fidx_str].seq
        return (pept in prot_seq)
    else:
        None
##############################
########################################################################################################################
def pept_info(row):
    pept,fetchid = row
    if __pd.notnull(fetchid):
        try:
            fidx_str = str(int(fetchid))
        except:
            fidx_str = str(fetchid)
        prot_seq = gbrecs[fidx_str].seq
        # find pept in prot:
        pept_found = prot_seq.find(pept)
        if pept_found > -1:
            # 1-based indexing right away ...
            start = prot_seq.find(pept) + 1
            stop  = prot_seq.find(pept) + len(pept)
            # because of 1-based indexing ...
            if stop >= len(prot_seq):
                print "peptide is at the end: ",pept,fidx_str
                next_aa = None
            else:
                next_aa = prot_seq[stop]
            ################################
            if start <= 1:
                print "peptide is at the start: ",pept,fidx_str
                prev_aa = None
            else:
                prev_aa = prot_seq[start-2]
            # return 4 columns ...
            return __pd.Series( {'start_fetched': start,
                'stop_fetched': stop,
                'prev_aa_fetched': prev_aa,
                'next_aa_fetched': next_aa} )
        else:
            print "(!!!) peptide not found: ",pept,fidx_str
            return __pd.Series({})
    else:
        print "fetchid is None for pept",pept
        return __pd.Series({})
########################################################################################################################
########################################################################################################################





#################################################################
# STAGE 3 STUFF, manuall gsites extraction ...
#################################################################
#
#
########################################################################################################################
__g_site = __re.compile(r'(?=(N[ACDEFGHIKLMNQRSTVWY][TS]))')
#########################################################################################
def get_theor_sites_number(prot_seq):
    # find all sites ...
    all_sites = [(site.start(),site.groups()[0]) for site in __g_site.finditer(prot_seq)]
    N_sites = len(all_sites)
    return N_sites
#########################################################################################
def get_theor_sites(prot_seq):
    # BEWARE ZERO-BASED INDEXING TO 1-BASED INDEXING TRANSLATION ...
    # find all sites ...
    all_sites = [(site.start(),site.groups()[0]) for site in __g_site.finditer(prot_seq)]
    # N_sites = len(all_sites)
    return ';'.join( (gsite_seq+'(%d)'%(pos+1)) for pos,gsite_seq in all_sites) # pos+1 - indexing transition ...
##################################################################################################################
#
##################################################################################################################
# looking for Deamidated sites only ...
__deamid = __re.compile('[D,d]eamidat')
def extract_deamids(mod_str,sequeon=None):
    # "IADDkYnDTFWk" with modifications: "Label:13C(6)15N(2) (+8), Deamidated:18O(1) (+3), Label:13C(6)15N(2) (+8)"
    mod_list = [ mod.strip() for mod in mod_str.split(',') ]
    if sequeon is not None:
        mod_locations = [(aa,idx+1) for idx,aa in enumerate(sequeon) if aa.islower()]
    else:
        mod_locations = [False for _ in mod_list]
    # return pd.Series( ms.parse_spectrum_modifications(mod) for mod in mod_list if bool(deamid.search(mod)) )
    # collecting results to return ...
    to_return = []
    for mod,mloc in zip(mod_list,mod_locations):
        if bool(__deamid.search(mod)):
            type_aa, gpos_pept, value = parse_spectrum_modifications(mod)
            if (type_aa is None) and (gpos_pept is None) and mloc:
                type_aa, gpos_pept    = mloc
            elif (type_aa is not None) and (gpos_pept is not None) and mloc:
                assert type_aa==mloc[0]
                assert gpos_pept==mloc[1]
            elif (type_aa is not None) and (gpos_pept is not None) and (not mloc):
                pass
            else:
                print "unknown case of deamidation extraction!!! mod: %s, seq: %s"%(mod_str,str(sequeon))
                sys.exit(1)
            # let's filter the aspartic ones with the value 3 right away ...
            if (type_aa in ['n','N']) and (__np.abs(value-3)<0.01):
                to_return.append( (type_aa, gpos_pept, value) )
    return __pd.Series( to_return )
##################################################################################################################
# unroll/exapnd spec table to account for multiple deamid-sites/gsites per peptide ...
def unroll_by_mfunc(df,col_names,mfunc,unrolled_colname='unrolled'):
    # get the column with arguments of the 'mfunc' ...
    thecolumns = df[col_names]
    if type(col_names) == list:
        # it appears as multi-column thing right after mfunc application ...
        multicol_mfunc_result = thecolumns.apply(mfunc,axis=1)
    elif type(col_names) == str:
        # it appears as multi-column thing right after mfunc application ...
        multicol_mfunc_result = thecolumns.apply(mfunc)
    else:
        print "unroll_by_mfunc expected col_names as a list OR a string ..."
        print "%s fails to meet these expectations. Exit"%str(col_names)
        __sys.exit(1)
    # stacked - multiindex is in use, level=1 index is the column name from 'multicol_mfunc_result' ...
    unrolled_mindex = multicol_mfunc_result.stack()
    # index is no longer uniq after the following operation ... we dropped inner indexing part.
    # IMPORTANT: indexes correspond to those from the original df (use it for merging later) ...
    unrolled_origindex = __pd.DataFrame(unrolled_mindex.reset_index(level=1,drop=True),columns=[unrolled_colname,])
    # merge unrolled_origindex (a single column with ambiguous index) with the original df ...
    # 'unrolled_origindex' must be DataFrame to merge: Series are not mergeable for some reason ...
    # the whole df is to be unrolled after the following operation.
    unrolled_df = df.merge(unrolled_origindex,left_index=True,right_index=True)
    #
    return unrolled_df.reset_index(drop=True)
##################################################################################################################
def deamid_to_gsite(deamid_mod, pept_pos, prot_seq):
    type_aa, gpos_pept, value = deamid_mod
    gpos_pept = int(gpos_pept)
    pept_pos = int(pept_pos)
    assert type_aa in ['n','N']
    assert __np.abs(value-3)<0.01
    # 'pept_pos' - is 1-based absolute poisition of the peptide in the protein ...
    # 'gpos_pept' - is 1-based relative position of gsite_start_N in the peptide ...
    gsite_start = pept_pos + gpos_pept-1 # 1-based coordinate ...
    gsite_stop  = pept_pos + gpos_pept-1 + 3-1 # 1-based coordinate ...
    # Due to slicing rules, we need [start-1:stop], no have position 'stop' included ...
    gsite_seq = prot_seq[gsite_start-1:gsite_stop]
    ############################################################
    # gstart must be 1-based for output ...
    return {'gsite':"%s(%d)"%(gsite_seq,gsite_start), 'gsite_seq':gsite_seq, 'gstart':gsite_start}
#########################################################################################
#







####################################################################################
## Code  ##   Description
####################################################################################
##  200  ##   The request was processed successfully.
##  400  ##   Bad request. There is a problem with your input.
##  404  ##   Not found. The resource you requested doesnt exist.
##  410  ##   Gone. The resource you requested was removed.
##  500  ##   Internal server error. Most likely a temporary problem, but if the problem persists please contact us.
##  503  ##   Service not available. The server is being updated, try again later.
######################################################################################



if __name__ == "__main__":
    pass





