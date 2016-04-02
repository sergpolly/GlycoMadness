from Bio import Entrez
from Bio import SeqIO
from StringIO import StringIO

Entrez.email = "sergpolly@gmail.com"


# this would be all possible fields that you can search in the protein database ...
print "These are the fields one can search in db='protein':"
handle = Entrez.einfo(db='protein')
prot_db_fields = Entrez.read(handle)
for idx,field in enumerate(prot_db_fields['DbInfo']['FieldList']):
    print idx+1, field['Name'], field['FullName'], field['Description']
handle.close()


# once all set, try actually searching something ...
handle = Entrez.esearch(db="protein", term="RecName[Title] AND PPT1[Gene Name] AND \"Homo sapiens\"[Organism]")
# handle = Entrez.esearch(db="protein", term="RecName[TITL] AND PPT1[GENE] AND \"Homo sapiens\"[ORGN]")
record = Entrez.read(handle)
handle.close()

results = []
for idx in record['IdList']:
    print "fetching",idx
    handle = Entrez.efetch(db='protein',id=idx,rettype='gb',retmode='text')
    filelike = StringIO(handle.read())
    seqrec = SeqIO.read(filelike,format='gb')
    results.append( seqrec )


# The features that Reid is asking about are accessible through:
r0 = results[0]
f4 = r0.features[4]
f4.qualifiers['region_name']










































