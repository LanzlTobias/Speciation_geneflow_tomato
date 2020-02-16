#!/bin/bash
#Script to prune SNPs for LD and then run a PCA with plink and ADMIXTURE them.
#First argument is the VCF file (gzipped or not)

FILE=${1%.gz}
FILE=${FILE%.vcf}


#Extract sites in high linkage
vcftools --gzvcf $FILE.vcf.gz --plink --out $FILE 2> tmp
rm tmp
sed -i 's/^0\t/1\t/g' $FILE.map
plink --file $FILE --indep-pairwise 100 20 0.1 --out $FILE --noweb --silent
sed -i 's/:/\t/g' $FILE.prune.in

#Prune the extracted sites
vcftools --gzvcf $FILE.vcf.gz --out $FILE.pruning --positions $FILE.prune.in --stdout --recode | gzip > $FILE.LDpruned.vcf.gz

#ADMIXTURE
plink --vcf $FILE.LDpruned.vcf.gz --make-bed --out $FILE --allow-extra-chr

awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim
rm $FILE.bim.tmp

for i in {2..8}
do
 admixture --cv=20 $FILE.bed $i > log$i.out
done
awk '/CV/ {print $3,$4}' log*out | cut -c 4,7-20 > $FILE.cv.error

#PCA
plink --vcf $FILE.LDpruned.vcf.gz --allow-extra-chr --pca --out $FILE

