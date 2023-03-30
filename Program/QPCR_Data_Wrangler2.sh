
# This ensures that proper directories are made
Home=$(pwd)
read -p  "How long is your amplicon?" amp

Rscript --vanilla Graph_maker.R $Home/int_files/edit_me.txt $amp

#./Graph_maker.R $Home/int_files/edit_me.txt $amp
mv *.png output/
mv Sample_Data.csv output/
rm *.pdf



