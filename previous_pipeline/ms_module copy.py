import sys as __sys
import requests as __requests
# from Bio import Seq as Seq
from Bio import SeqIO as __SeqIO
from Bio import SeqRecord as __SeqRecord
import StringIO as __StringIO
import re as __re

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






# n4: Deamidated:18O(1) (+2.99)
def parse_spectrum_modifications(modifier):
    # print modifier
    loc_mod = modifier.strip()
    loc_mod = loc_mod.split(' ')
    mod_type = loc_mod[0].strip(':')
    # aa modified and peptide position ...
    # POSITIONS ARE MOST LIKELY TO HAVE 1-BASED INDEX ...
    mod_type_aa, mod_type_pos = mod_type[0], int(mod_type[1:])
    # value ...
    mod_val = float(loc_mod[2].strip('()'))
    #
    return (mod_type_aa, mod_type_pos, mod_val)




# protein name example:
# sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha chain OS=Homo sapiens GN=HLA-A PE=1 SV=2
# regexp to use often:
# first part of full protein name with Uid and Locus ...
__left_part = __re.compile("[a-z]{2}\|[A-Z0-9\-]+\|[A-Z0-9\_]+")
# features like OS,GN etc, simply to extract what fatures are present ...
__feature_types = __re.compile("[A-Z]{2}=")
def parse_prot_name(prot_name):
    """Functions parses provided protein name and returns a dict with: uid,locus,prot_name,GeneName,OrgSource
    -------------------------
    prot name expected format:
    'sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha chain OS=Homo sapiens GN=HLA-A PE=1 SV=2'"""
    # extract left most part:
    lp_extracted = __left_part.findall(prot_name)
    if len(lp_extracted)!=1:
        print "Could not match left part of protein name: %s"%prot_name
        print "Extraction result is: ", lp_extracted
        # sys.exit(1)
        # return None
        return { 'prot_name':prot_name }
    lp_extracted, = lp_extracted
    #
    # what features are present ...
    f_types = [ f.strip('=') for f in  __feature_types.findall(prot_name) ]
    # generate regexp based on the combination of features:
    # for example: "GN=(.+)PE=(.+)SV=(.+)"
    f_pattern = ''.join(['%s=(.+)'%f for f in f_types])
    # find the whole pattern in the protein name:
    f_values, = __re.findall(f_pattern,prot_name)
    # right features part for stripping ...
    rp_extracted = ''.join( "%s=%s"%(k,v) for k,v in zip(f_types,f_values) )
    #
    # strip left and right part of the full prot name, to get the human readable
    prot_name_extracted = prot_name.replace(lp_extracted,'').replace(rp_extracted,'').strip()
    #
    _,uid,locus = lp_extracted.split('|')
    f_dict = dict( zip(f_types,f_values) )
    # returning all extracted information ...
    if ('GN' not in f_types)or('OS' not in f_types) :
        print "There is no GeneName or OrganismSource in the protein name: %s"%prot_name
        print "Feature types extracted are: ",f_types
        # sys.exit(1)
        return { 'uid':uid,
                'locus':locus,
                'prot_name':prot_name_extracted }
    else:
        return { 'uid':uid,
                'locus':locus,
                'prot_name':prot_name_extracted,
                'GN':f_dict['GN'].strip(),
                'OS':f_dict['OS'].strip() }







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





