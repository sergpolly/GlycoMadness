
echo ""
echo "(1) 011616 glycocapture 90-80 ..."
echo ""
sh pipeline_generic.sh \
"../raw_data/New_files_to_analyze/011616 glycocapture 90-80" \
"../Outputs/output_01212016" > \
../Outputs/output_01212016/report.txt


echo ""
echo "(2) before december 25 ..."
echo ""
sh pipeline_generic.sh \
"../raw_data/New_files_to_analyze/original_input_before_dec25" \
"../Outputs/original_output_before_dec25" > \
../Outputs/original_output_before_dec25/report.txt


echo ""
echo "(3) 111815 glycocapture 90-80 ..."
echo ""
sh pipeline_generic.sh \
"../raw_data/New_files_to_analyze/111815 glycocapture 90-80" \
"../Outputs/new_output_small_111815" > \
../Outputs/new_output_small_111815/report.txt


echo ""
echo "(4) 121515 glycocapture 90-80 ..."
echo ""
sh pipeline_generic.sh \
"../raw_data/New_files_to_analyze/121515 glycocapture 90-80" \
"../Outputs/new_output_large_121515" > \
../Outputs/new_output_large_121515/report.txt













