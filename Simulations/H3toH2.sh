#!/bin/bash
Time0=$(date +%s)

## pop order H1 (1), H2 (2), H3 (3), H4 (4)

for i in {1..100}
do
sim_name=$(echo "H3toH2_$i")

echo "$sim_name"

Ne=2e4
mu=5.1e-9

divT1=$(echo $Ne | awk '{print (200000/5)/($1*4)}')   ## divergence time H2 and H1
divT2=$(echo $Ne | awk '{print (1000000/5)/($1*4)}')  ## divergence time H3 and H2
divT3=$(echo $Ne | awk '{print (2000000/5)/($1*4)}') ## divergence time H4 and H3

migT=$(echo $Ne | awk '{print (50000/5)/($1*4)}')   ## migration time from H3 to H2
migR=$(echo $Ne | awk '{print 1e-4*$1*4}')          ## migration rate from H3 to H2

while read p; do
        name=`echo $p | awk '{print $1}'`
        size=`echo $p | awk '{print $2}'`
        theta=$(echo $Ne $mu $size | awk '{print $1*$2*$3*4}')
        rho=$theta
/data/proj/teaching/NGS_course/Softwares/scrm 24 1 -t $theta -r $rho $size -I 4 2 2 10 10 \
        -m 2 3 $migR -em $migT 2 3 0 \
        -ej $divT1 1 2 -ej $divT2 2 3 -ej $divT3 3 4
done < chroms_lengths.txt > sim_$sim_name

echo "###scrm Finished, run time in min:" 
awk "BEGIN {print ($(date +%s)-$Time0)/60}"

gsize=`awk '{sum+=$2} END {print sum}' chroms_lengths.txt `

/data/proj/teaching/NGS_course/bin/ms2vcf sim_$sim_name 1 $gsize
mv d0.vcf sim_$sim_name.vcf

echo "### ms2vcf finished, run time in min:"
awk "BEGIN {print ($(date +%s)-$Time0)/60}"

python parseVCF.py -i sim_"$sim_name".vcf > sim_"$sim_name".geno

echo "### parseVCFs finished, run time in min:"
awk "BEGIN {print ($(date +%s)-$Time0)/60}"

sed -i 's/|//g' sim_"$sim_name".geno

sed -i '1d' sim_"$sim_name".geno

perl convertDoGeno2egglib.pl sim_"$sim_name".geno > sim_"$sim_name".geno.txt

gzip sim_"$sim_name".geno

cat header.txt sim_"$sim_name".geno.txt > sim_"$sim_name".geno_WH.txt

rm sim_"$sim_name".geno.txt

echo "### Transforming finished, run time in min: "
awk "BEGIN {print ($(date +%s)-$Time0)/60}"

perl ABBA_BABA.v1.pl sim_$sim_name.geno_WH.txt pop.file pop.order > D_fd_Fhom_ANGSD_$sim_name.txt

perl ABBA_out_blocker.pl D_fd_Fhom_ANGSD_$sim_name.txt > D_fd_Fhom_ANGSD_block_$sim_name.txt

Rscript Jackknife_ABBA_pipe.R D_fd_Fhom_ANGSD_block_$sim_name.txt D_fd_Fhom_ANGSD_jackknife_$sim_name.txt

tail -n 3 D_fd_Fhom_ANGSD_$sim_name.txt > tmp1
tail -n 3 D_fd_Fhom_ANGSD_jackknife_$sim_name.txt > tmp2

paste tmp1 tmp2 > D_fd_Fhom_ANGSD_2_$sim_name.txt
rm tmp1
rm tmp2

Rscript ABBA_pvalue.R D_fd_Fhom_ANGSD_2_$sim_name.txt D_fd_Fhom_ANGSD_final_$sim_name.txt

echo "### AbbaBaba finished, run time in min: "
awk "BEGIN {print ($(date +%s)-$Time0)/60}"

done
