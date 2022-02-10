
# This ensures that proper directories are made
Home=$(pwd)


./Graph_maker.R $Home/int_files/edit_me.txt
mv *.png output/
rm *.pdf



