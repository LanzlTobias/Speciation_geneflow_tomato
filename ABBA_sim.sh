#!/bin/bash
## This script simulates populations after the format of the ABBA BABA test
# -m takes the migration rate between H3 and H2
# -d takes a file with the three divergence times separated by tabs per line in the format <H1-H2\tH2-H3\tH3-H4\tNAME> with name indicating a string that will be added to the simulation name
# -u takes a constant mutation rate
# -p takes the number of repetition that should be performed per simulation scenario
# -n takes the Ne 
# -v takes the time at which the migration ends
# -o takes the name of the 
# -c takes a file with the chromosome names and lengths separated by tabs (one chromosome per line)
# -s takes the fraction of Ne that H2 has
# -x takes the fraction of Ne that H1 has


while getopts m:d:u:t:p:n:v:o:c: option
do
case "${option}"
in
m) migRate=${OPTARG};;
d) divTime=${OPTARG};;
u) mu=${OPTARG};;
p) rep=${OPTARG};;
n) Ne=${OPTARG};;
v) migTime=${OPTARG};;
o) output=${OPTARG};;
c) chromosome=${OPTARG};;
s) frac_h2=${OPTARG};;
x) frac_h1=${OPTARG};;
esac
done

n_div=$(cat $divTime | wc -l)


for (( i=1; i<=$rep; i++ ))
do
echo 'Starting repetion '$i
while read j; do	
	div_name=$(echo $j | awk '{print $4}')
	sim_name=$(echo $output'_'$div_name'_'$i)
echo 'Starting '$sim_name
divT1=$(echo -e $Ne"\t"$j | awk '{print ($2/5)/($1*4)}')   ## divergence time H2 and H1
divT2=$(echo -e $Ne"\t"$j | awk '{print ($3/5)/($1*4)}')  ## divergence time H3 and H2
divT3=$(echo -e $Ne"\t"$j | awk '{print ($4/5)/($1*4)}') ## divergence time H4 and H3

migT=$(echo -e $Ne"\t"$migTime | awk '{print ($2/5)/($1*4)}')   ## migration time from H3 to H2
migR=$(echo $Ne"\t"$migRate | awk '{print $2*$1*4}')          ## migration rate from H3 to H2



while read p; do
	name=$(echo $p | awk '{print $1}')
	size=$(echo $p | awk '{print $2}')
        theta=$(echo $Ne $mu $size | awk '{print $1*$2*$3*4}')
        rho=$theta
/data/proj/teaching/NGS_course/Softwares/scrm 24 1 -t $theta -r $rho $size -I 4 2 2 10 10 \
        -m 2 3 0 -em $migT 2 3 $migR \
	-en 0 2 $frac_h2 -en 0 1 $frac_h1 \
        -ej $divT1 1 2 \
        -ej $divT2 2 3 \
        -ej $divT3 3 4


done < $chromosome > $sim_name

echo "###scrm finished" 

gsize=$(awk '{sum+=$2} END {print sum}' $chromosome)

/data/proj/teaching/NGS_course/bin/ms2vcf $sim_name 1 $gsize
mv d0.vcf $sim_name.vcf
rm $sim_name

echo "### ms2vcf finished"
done < $divTime
done
