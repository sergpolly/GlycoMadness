%run SILAC_stage0_spec_n_summary_check.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer --verbose -r SILAC1(STT3A)_rawdatareport.csv -q SILAC1(STT3A)_QuantSpecReport.csv --separator comma
# %run SILAC_stage0_spec_n_summary_check.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer --verbose -p Glycocap4-peptides.csv -s Glycocap4-spectra.csv --separator comma
%run SILAC_stage0_spec_n_summary_check.py --prefix ../raw_data/New_files_to_analyze/SILAC10B --verbose -r SILAC10B_rawdatareport.csv -q SILAC10B_quantspectrareport.csv --separator comma
                                             
        
# next stage simply requires any file with Accession Number column, and Protein Name. So potentially tha can be Quant  either/or Raw_data
%run SILAC_stage1a_aquire_fetchid_NCBI.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -r SILAC1(STT3A)_rawdatareport.csv --email sergey.venev@umassmed.edu --separator comma  --swiss
# %run stage1a_aquire_fetchid_NCBI.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -p Glycocap4-peptides.csv --email sergey.venev@umassmed.edu --separator comma --swiss
%run SILAC_stage1a_aquire_fetchid_NCBI.py --prefix ../raw_data/New_files_to_analyze/SILAC10B -r SILAC10B_rawdatareport.csv --email sergey.venev@umassmed.edu --separator comma  --swiss


%run SILAC_stage1b_pull_proteins_NCBI.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -f SILAC_plus_fetch.csv --email sergey.venev@umassmed.edu
# %run SILAC_stage1b_pull_proteins_NCBI.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -f peptides_with_fetch.csv --email sergey.venev@umassmed.edu
%run SILAC_stage1b_pull_proteins_NCBI.py --prefix ../raw_data/New_files_to_analyze/SILAC10B -f SILAC_plus_fetch.csv --email sergey.venev@umassmed.edu


%run  SILAC_stage2_assign_pulled_proteins.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -f SILAC_plus_fetch.csv -g pulled_proteins.gb --threshold 100
# %run  SILAC_stage2_assign_pulled_proteins.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -f peptides_with_fetch.csv -g pulled_proteins.gb
%run  SILAC_stage2_assign_pulled_proteins.py --prefix ../raw_data/New_files_to_analyze/SILAC10B -f SILAC_plus_fetch.csv -g pulled_proteins.gb --threshold 100


# 	not working!!!!! modifications extraction to be fixed later on ...
%run  SILAC_stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -e 4 -m pept_prot_map.csv -g pulled_proteins.gb -q SILAC1(STT3A)_QuantSpecReport.csv --separator comma
%run  SILAC_stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/SILAC5A -e 12 -m pept_prot_map.csv -g pulled_proteins.gb -q SILAC5A_quant\ spectra.csv --separator comma
# %run  SILAC_stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -e 4 -m pept_prot_map.csv -g pulled_proteins.gb -s Glycocap4-spectra.csv --separator comma
%run  SILAC_stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/SILAC10B -e 21 -m pept_prot_map.csv -g pulled_proteins.gb -q SILAC10B_quantspectrareport.csv --separator comma
%run  SILAC_stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/SILAC10A -e 22 -m pept_prot_map.csv -g pulled_proteins.gb -q Silac10A_quantspectrareport.csv --separator comma



# ATTENTION!!! in case of ambiguous peptide - proteins matching, we need more than 1 column to match quant data and raw_w_fetch ...
%run  SILAC_stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -e 4 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -q SILAC1(STT3A)_QuantSpecReport.csv --separator comma
%run  SILAC_stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/SILAC5A -e 12 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -q SILAC5A_quant\ spectra.csv --separator comma
# %run  SILAC_stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/silac_test_summer -e 4 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -s Glycocap4-spectra.csv --separator comma
%run  SILAC_stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/SILAC10B -e 21 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -q SILAC10B_quantspectrareport.csv --separator comma



# RecName[Title] AND "FAM234A_HUMAN"[locus] AND "Homo sapiens"[Organism],
# "Q9HOX4"[accession] AND "Homo sapiens"[Organism],
# RecName[Title] AND "ITFG3"[Gene Name] AND "Homo sapiens"[Organism]


