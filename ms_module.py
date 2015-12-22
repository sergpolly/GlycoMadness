import sys
import requests
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import StringIO


####################################################################################
# Code  #   Description
####################################################################################
#  200  #   The request was processed successfully.
#  400  #   Bad request. There is a problem with your input.
#  404  #   Not found. The resource you requested doesn’t exist.
#  410  #   Gone. The resource you requested was removed.
#  500  #   Internal server error. Most likely a temporary problem, but if the problem persists please contact us.
#  503  #   Service not available. The server is being updated, try again later.
######################################################################################
web_request_status_collection = {200:"The request was processed successfully.",
400:"Bad request. There is a problem with your input.",
404:"Not found. The resource you requested doesn’t exist.",
410:"Gone. The resource you requested was removed.",
500:"Internal server error. Most likely a temporary problem, but if the problem persists please contact us.",
503:"Service not available. The server is being updated, try again later."}


# # typical Uniprot ID (for isoform 2) ...
# uid = "P16870-2"


def get_uniprot(uid):
    # the way we form Uniprot ID request ...
    get_uid_url = lambda _: "http://www.uniprot.org/uniprot/%s.fasta"%_
    # web request for a given uid ...
    uid_url = get_uid_url(uid)
    # make a request ...
    req_res = requests.get(uid_url)
    # check request status ...
    if req_res.status_code==200 and bool(req_res.content):
        # to be read by SeqIO ...
        string_as_handle = StringIO.StringIO(req_res.content)
        seq_rec = SeqIO.read(like_a_filehandle)
        return seq_rec
    elif req_res.status_code==200:
        print web_request_status_collection[req_res.status_code]
        print "... But, the content is empty for accession number %s!"%uid
        sys.exit(1)
    elif req_res.status_code in web_request_status_collection:
        print web_request_status_collection[req_res.status_code]
        sys.exit(1)
    else:
        print "Unknown status code returned!"
        sys.exit(1)




def stupid_aligner(peptide,protein):
    # small subfunction to get the hamming distance ...
    def hamming_distance(seq1,seq2):
        assert len(seq1)==len(seq2)
        mismatches = sum( (l1!=l2) for l1,l2 in zip(seq1,seq2) )
        return mismatches
    # sequences to strings, making sure they are SeqRec or Seq entyties ...
    peptide_str = str(peptide.seq) if type(peptide)==SeqRecord.SeqRecord else str(peptide)
    protein_str = str(protein.seq) if type(protein)==SeqRecord.SeqRecord else str(protein)
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
    print "Beware! Best alignment found has %d mismatches"%min_mismatch
    return align_frame










