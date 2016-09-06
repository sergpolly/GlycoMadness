INPUT="../raw_data/New_files_to_analyze/SILAC9NGI"
EXPNUM=20
RAW="SilacNGI-1_raw datareport.csv"
QUANT="SILACNGI1-1_quantspectra report.csv"

# goto the GlycoMadness folder ...
# and execute from there:

echo "stage 0"
# stage 0 is simply checking the consistency between spectrum report and peptide summary
# it produces readable outputs: unique peptides in both files, their intersections along with
# unique proteins in both files and their intersection ...
# Proteins and peptides are compared independently, using "Protein name" and "Peptide sequence" columns ...
python SILAC_stage0_spec_n_summary_check.py \
    --prefix "$INPUT" \
    --verbose \
    --separator comma \
    -r "$RAW" \
    -q "$QUANT"
# NOT_IN kinda files are also produced by this script ...



echo "stage 1a and 1b ..."
# stage 1 performs protein records fetching from NCBI, using parsed Protein name from Peptide summary
# it identifies fetchid-s first, and then it downloads corresponding records in genbank format ...
# not all proteins might end up with the fetchids, moreso,
# there might be discrepancies between spectrum and pept_sum in terms of # of proteins these files are reffering to.
python SILAC_stage1a_aquire_fetchid_NCBI.py \
    --prefix "$INPUT" \
    --separator comma \
    -r "$RAW" \
    --email sergey.venev@umassmed.edu \
    --swiss 
# pull em ...
python SILAC_stage1b_pull_proteins_NCBI.py \
    --prefix "$INPUT" \
    -f SILAC_plus_fetch.csv \
    --email sergey.venev@umassmed.edu 



echo "stage 2"
# Stage2 is simply assigning pulled protein sequences to the peptides listed in summary...
# We're trying to choose best peptide-protein combinations according to some rules and criteria.
# Ambiguous and/or underachieving are reported in the BAD file ...
python SILAC_stage2_assign_pulled_proteins.py \
    --prefix "$INPUT" \
    -f SILAC_plus_fetch.csv \
    -g pulled_proteins.gb \
    --threshold 100




echo "stage 3"
# Stage 3 is the place where we are finally touching upon actual Peptide spectrum(!) ...
# we're extracting glycosilation information and MERGING already processed peptide summary with spectrum report.
# Spec and pep are merged using peptides, so if there were any mismatches between 2 files on the "Protein" level,
# than Peptide summary associated nomenclature would end up in the FINAL file ...
python SILAC_stage3_gsites_catalog.py \
    -e "$EXPNUM" \
    --prefix "$INPUT" \
    --separator comma \
    -m pept_prot_map.csv \
    -g pulled_proteins.gb \
    -q "$QUANT"



echo "stage 4"
# combine BAD and GOOD in one excel spreadsheet ...
python SILAC_stage4_eval_bad_combine_out.py \
    --prefix "$INPUT" \
    -e "$EXPNUM" \
    -b BAD_pept_prot.csv \
    -o FINAL_gsite_anthology.csv \
    -g pulled_proteins.gb \
    -q "$QUANT" \
    --separator comma



# TO DO
# adadress 'enzyme' issue ...
# discrepancies between pepts and spec summay - what to do?!
# END OF FILE?




