DIR="../raw_data/New_files_to_analyze"
FNAME="pulled_proteins.gb"
DESTINATION="../PULLED_PROTEINS_TOTAL"

# list of the folders with the pulled_proteins.gb files of interest ...
datalist="
    glycocapture10_by_Reid
    glycocapture11_by_Reid
    glycocapture12_by_Reid
    glycocapture13_by_Reid
    glycocapture14_by_Reid
    glycocapture15_by_Reid
    glycocapture1_by_Reid
    glycocapture2_by_Reid
    glycocapture3_by_Reid
    glycocapture4_by_Reid
    glycocapture5_by_Reid
    glycocapture6_by_Reid
    glycocapture7_by_Reid
    glycocapture8_by_Reid
    glycocapture9_by_Reid
    SILAC1A
    SILAC1B
    SILAC2A
    SILAC2B
    SILAC3A
    SILAC3B
    SILAC4A
    SILAC4B
    SILAC5A
    SILAC5B
    SILAC6A
    SILAC6B"

count=0
for d in $datalist; do
    DEST_FNAME="$count-$FNAME"
    echo "copying file number" $count":"  $DIR/$d/$FNAME " to " $DEST_FNAME
    cp $DIR/$d/$FNAME $DESTINATION/$DEST_FNAME
    count=$(( $count + 1 ))
done














