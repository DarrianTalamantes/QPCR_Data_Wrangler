
# This ensures that proper directories are made
Home=$(pwd)
int_files=$Home/int_files
output=$Home/output
if [ ! -e $int_files ] ; then mkdir $int_files; fi
if [ ! -e $output ] ; then mkdir $output; fi


read -p  "Input file from qPCR machine " infile
infile2="${infile:1:-1}"
echo "Formatting" $infile2
cat $infile2 | sed 1d | cut -f 3,5 | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' | sed s'/Cp/Cp\tTreatment\tEndoPos_Neg_Water\tPrimer_Set\tConcentration/'g > $int_files/edit_me.txt
echo
echo "Done formatting please input qualitative data in designated spots"
echo "For standards label them with stdX, X being a unique number for every standard, under the Treatment section. Input the concentration of each standard under the concentration column."
echo "For water control use water"
echo "everything else can be whatever you want"




