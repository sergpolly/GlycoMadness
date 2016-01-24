INPUT="../raw_data/New_files_to_analyze/121515 glycocapture 90-80"
OUTPUT="../Outputs/new_output_large_121515"


# goto the GlycoMadness folder ...
# and execute from there:

echo "Backing up required protein sequences from UniProt ..."
# get needed protein sequences from uniprot and store em locally ...
# kinda like that:
# python  backup_uids.py  ../input_peptides.xls  backup_uids_output.fasta
if [ ! -f "$OUTPUT"/backup_uids_output.fasta ];
then
	echo "No pre-fetched .fasta file detected, fetching ..."
	python backup_uids.py \
	"$INPUT"/peptides.xls \
	"$OUTPUT"/backup_uids_output.fasta
else
	echo "Pre-fetched .fasta file detected. Skip fetching ..."
fi




echo "Unique peptides identification and clean up ..."
# then do some clean up in the 'proteins' file using the same proteins file and proteins backup
# python  resolve_peptide_ambiguity.py  ../input_peptides.xls  backup_uids_output.fasta  uniq_peptides_catalog.sv
python resolve_peptide_ambiguity.py \
	"$INPUT"/peptides.xls \
	"$OUTPUT"/backup_uids_output.fasta \
	"$OUTPUT"/uniq_peptides_catalog.csv



echo "Finally, extracting G-sites ..."
# finally extract all gsites and make final output ...
# python  gsites_catalog.py  ../peptides.xls  ../specs.xls  backup_uids_output.fasta  uniq_peptides_catalog.csv gsites_antology.csv
python  gsites_catalog.py \
	"$INPUT"/peptides.xls \
	"$INPUT"/specs.xls \
	"$OUTPUT"/backup_uids_output.fasta \
	"$OUTPUT"/uniq_peptides_catalog.csv \
	"$OUTPUT"/gsites_antology.csv














