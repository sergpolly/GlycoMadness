import sys as __sys
import requests as __requests
# from Bio import Seq as Seq
from Bio import SeqIO as __SeqIO
from Bio import SeqRecord as __SeqRecord
import StringIO as __StringIO


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
    """ Zero-based number of peptide occurance in the protein ..."""
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
    return align_frame






# n4: Deamidated:18O(1) (+2.99)
def parse_spectrum_modifications(modifier):
    # print modifier
    loc_mod = modifier.strip()
    loc_mod = loc_mod.split(' ')
    mod_type = loc_mod[0].strip(':')
    # aa modified and peptide position ...
    mod_type_aa, mod_type_pos = mod_type[0], int(mod_type[1:])
    # value ...
    mod_val = float(loc_mod[2].strip('()'))
    #
    return (mod_type_aa, mod_type_pos, mod_val)




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





