#!/bin/bash

# objective: This will combine the copy number data and the data for the
# alkaloid analysis into one file.

Home=$(pwd)


read -p  "Input copy number / biomass file " infile
infile2="${infile:1:-1}"
echo
echo "****************************************"
echo "biomass file is" $infile2
echo "****************************************"

read -p  "input alkaloid file " infile3
infile4="${infile3:1:-1}"
echo
echo "****************************************"
echo "alkaloid file is" $infile4
echo "****************************************"

Rscript --vanilla Combo_Maker.R $infile2 $infile4

sed -i 's/Sample/<Trait>/g' comboData.csv
sed -i 's/,/\t/g' comboData.csv
sed -i 's/"//g' comboData.csv
mv comboData.csv comboData.txt


mv comboData.txt Output/

