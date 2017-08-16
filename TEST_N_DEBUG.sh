%run stage0_spec_n_summary_check.py --prefix ../raw_data/New_files_to_analyze/glycocapture1_by_Reid --verbose -p Glycocap1-peptides.csv -s Glycocap1-spectra.csv --separator comma
%run stage0_spec_n_summary_check.py --prefix ../raw_data/New_files_to_analyze/glycocapture2_by_Reid --verbose -p Glycocap2-peptides.csv -s Glycocap2-spectra.csv --separator comma
%run stage0_spec_n_summary_check.py --prefix ../raw_data/New_files_to_analyze/glycocapture3_by_Reid --verbose -p Glycocap3-peptides.csv -s Glycocap3-spectra.csv --separator comma
%run stage0_spec_n_summary_check.py --prefix ../raw_data/New_files_to_analyze/glycocapture4_by_Reid --verbose -p Glycocap4-peptides.csv -s Glycocap4-spectra.csv --separator comma
                                             
                                             
%run stage1a_aquire_fetchid_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture1_by_Reid -p Glycocap1-peptides.csv --email sergey.venev@umassmed.edu --separator comma
%run stage1a_aquire_fetchid_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture2_by_Reid -p Glycocap2-peptides.csv --email sergey.venev@umassmed.edu --separator comma
%run stage1a_aquire_fetchid_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture3_by_Reid -p Glycocap3-peptides.csv --email sergey.venev@umassmed.edu --separator comma
%run stage1a_aquire_fetchid_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture4_by_Reid -p Glycocap4-peptides.csv --email sergey.venev@umassmed.edu --separator comma --swiss
%run stage1a_aquire_fetchid_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture20_by_Reid -p Glycocapture20_peptide\ report.csv --email sergey.venev@umassmed.edu --separator comma --swiss


%run stage1b_pull_proteins_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture1_by_Reid -f peptides_with_fetch.csv --email sergey.venev@umassmed.edu
%run stage1b_pull_proteins_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture2_by_Reid -f peptides_with_fetch.csv --email sergey.venev@umassmed.edu
%run stage1b_pull_proteins_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture4_by_Reid -f peptides_with_fetch.csv --email sergey.venev@umassmed.edu
%run stage1b_pull_proteins_NCBI.py --prefix ../raw_data/New_files_to_analyze/glycocapture20_by_Reid -f peptides_with_fetch.csv --email sergey.venev@umassmed.edu




%run  stage2_assign_pulled_proteins.py --prefix ../raw_data/New_files_to_analyze/glycocapture1_by_Reid -f peptides_with_fetch.csv -g pulled_proteins.gb
%run  stage2_assign_pulled_proteins.py --prefix ../raw_data/New_files_to_analyze/glycocapture2_by_Reid -f peptides_with_fetch.csv -g pulled_proteins.gb
%run  stage2_assign_pulled_proteins.py --prefix ../raw_data/New_files_to_analyze/glycocapture3_by_Reid -f peptides_with_fetch.csv -g pulled_proteins.gb
%run  stage2_assign_pulled_proteins.py --prefix ../raw_data/New_files_to_analyze/glycocapture4_by_Reid -f peptides_with_fetch.csv -g pulled_proteins.gb
%run  stage2_assign_pulled_proteins.py --prefix ../raw_data/New_files_to_analyze/glycocapture20_by_Reid -f peptides_with_fetch.csv -g pulled_proteins.gb




%run  stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/glycocapture1_by_Reid -e 1 -m pept_prot_map.csv -g pulled_proteins.gb -s Glycocap1-spectra.csv --separator comma
%run  stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/glycocapture2_by_Reid -e 2 -m pept_prot_map.csv -g pulled_proteins.gb -s Glycocap2-spectra.csv --separator comma
%run  stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/glycocapture3_by_Reid -e 3 -m pept_prot_map.csv -g pulled_proteins.gb -s Glycocap3-spectra.csv --separator comma
%run  stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/glycocapture4_by_Reid -e 4 -m pept_prot_map.csv -g pulled_proteins.gb -s Glycocap4-spectra.csv --separator comma
%run  stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/glycocapture12_by_Reid -e 12 -m pept_prot_map.csv -g pulled_proteins.gb -s Glycocap12_spectra\ report.csv --separator comma
%run  stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/glycocapture13_by_Reid -e 13 -m pept_prot_map.csv -g pulled_proteins.gb -s Glycocap13_spectra\ report.csv --separator comma
%run  stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/glycocapture17_by_Reid -e 17 -m pept_prot_map.csv -g pulled_proteins.gb -s Glcocap17_spectra\ report.csv --separator comma
%run  stage3_gsites_catalog.py --prefix ../raw_data/New_files_to_analyze/glycocapture20_by_Reid -e 20 -m pept_prot_map.csv -g pulled_proteins.gb -s Glycocap20_spectrumreport.csv --separator comma


%run  stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/glycocapture1_by_Reid -e 1 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -s Glycocap1-spectra.csv --separator comma
%run  stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/glycocapture2_by_Reid -e 2 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -s Glycocap2-spectra.csv --separator comma
%run  stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/glycocapture4_by_Reid -e 4 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -s Glycocap4-spectra.csv --separator comma
%run  stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/glycocapture12_by_Reid -e 12 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -s Glycocap12_spectra\ report.csv --separator comma
%run  stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/glycocapture13_by_Reid -e 13 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -s Glycocap13_spectra\ report.csv --separator comma
%run  stage4_eval_bad_combine_out.py --prefix ../raw_data/New_files_to_analyze/glycocapture20_by_Reid -e 20 -b BAD_pept_prot.csv -o FINAL_gsite_anthology.csv -g pulled_proteins.gb -s Glycocap20_spectrumreport.csv --separator comma


# RecName[Title] AND "FAM234A_HUMAN"[locus] AND "Homo sapiens"[Organism],
# "Q9HOX4"[accession] AND "Homo sapiens"[Organism],
# RecName[Title] AND "ITFG3"[Gene Name] AND "Homo sapiens"[Organism]


