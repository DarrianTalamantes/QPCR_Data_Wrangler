
# This code  is specifically for my biomass Analysis . Not for just any QPCR Data
Home=$(pwd)
read -p  "Input your Endophyte Data Created from QPCR_Data_Wrangler2.sh" endo
endo2="${endo:1:-1}"
echo
echo
read -p  "Input your Fescue Data Created from QPCR_Data_Wrangler2.sh" fescue
fescue2="${fescue:1:-1}"
echo
echo


Rscript --vanilla Relative_Biomass.R $endo2 $fescue2
mv Biomass_Data.csv output/Data_for_project/Biomass_Data/
echo "Analysis complete"