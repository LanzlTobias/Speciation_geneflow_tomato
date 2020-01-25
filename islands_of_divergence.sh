#!/bin/bash
#Required scripts are: 
#popgenWindows.py available at https://github.com/simonhmartin/genomics_general
#parseVCF.py available at: https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing
#convertDoGeno2egglib.pl available at: https://github.com/owensgl/abba_baba#
#outlier_fst_dxy.R available at: https://github.com/LanzlTobias/Speciation_geneflow_tomato


FILE=${1%.gz}
FILE=${file%.vcf}

if [ -s $FILE.vcf.gz ]
then 
 # Get a .map and .ped file
 vcftools --gzvcf $FILE.vcf.gz --SNPdensity 10000 --out $FILE"_10kb"
else
 FILE=${FILE%.vcf}
  vcftools --vcf $FILE.vcf --SNPdensity 10000 --out $FILE"_10kb"
fi

python $script_dir/parseVCF.py -i $FILE.vcf.gz > $FILE.geno

sed -i 's/|//g' $FILE.geno

sed -i '1d' $FILE.geno

perl $script_dir/convertDoGeno2egglib.pl $FILE.geno > $FILE.geno.txt

cat header.txt $FILE.geno.txt > $FILE"_WH.geno"

rm $FILE.geno.txt $FILE.geno

python popgenWindows.py -w 10000 -g $FILE"_WH.geno" -o popgenWindows_10Kb.csv -f diplo -T 5 -p LA1963 --popsFile pop.file2 -p LA2931 --popsFile pop.file2 -p LA2932 --popsFile pop.file2 -p LA3111 --popsFile pop.file2 -p LA4107 --popsFile pop.file2 -p LA4330 --popsFile pop.file2 -p LA_Peruvianum --popsFile pop.file2 -p PI_Peruvianum --popsFile pop.file2

python popgenWindows.py -w 100000 -g $FILE"_WH.geno" -o popgenWindows.csv -f diplo -T 5 -p LA1963 --popsFile pop.file2 -p LA2931 --popsFile pop.file2 -p LA2932 --popsFile pop.file2 -p LA3111 --popsFile pop.file2 -p LA4107 --popsFile pop.file2 -p LA4330 --popsFile pop.file2 -p LA_Peruvianum --popsFile pop.file2 -p PI_Peruvianum --popsFile pop.file2

Rscript outlier_fst_dxy.R
