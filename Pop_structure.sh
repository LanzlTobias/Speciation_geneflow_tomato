#Script to prune SNPs for LD and then run a PCA with plink, ADMIXTURE, and TreeMix on them.
#The transformation of the file for TreeMix was taken from https://speciationgenomics.github.io/Treemix/ 

#First argument is the VCF file (gzipped or not)
#Second argument is the cluster file required for TreeMix (Example file solanum.clust in the repo)
#The python script plink2treemix.py has to be in the working directory, otherwise a path has to be given as a third argument

FILE=${1%.gz}
FILE=${file%.vcf}

if [ -z "$3" ]
  then
    script_dir=$pwd
  else
    script_dir=$3
fi


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

#Treemix
./vcf2treemix.sh $FILE.vcf.gz solanum.clust
## Adjust the file to fit for treemix
vcftools --gzvcf $FILE.LDpruned.vcf.gz --plink --mac 2 --remove-indels --max-alleles 2 --out $FILE
awk -F"\t" '{
        split($2,chr,":")
        $1="1"
        $2="1:"chr[2]
        print $0
}' $FILE.map > better.map
mv better.map $FILE.map
plink --file $FILE --make-bed --out $FILE --allow-no-sex --allow-extra-chr 0
plink --bfile $FILE --freq --missing --within $2 --out $FILE --allow-no-sex --allow-extra-chr 0

gzip $FILE.frq.strat
$script_dir/plink2treemix.py $FILE.frq.strat.gz $FILE.treemix.frq.gz
gunzip $FILE.treemix.frq.gz
gunzip $FILE.frq.strat.gz

awk 'BEGIN{print "scaffold_pos\tscaffold\tpos"}{split($2,pos,":");print $2"\t"pos[1]"\t"pos[2]}' $FILE.map > $FILE.positions
paste $FILE".positions" $FILE".treemix.frq" > $FILE.frequencies

awk '{printf $0
        for(i = 4; i <= NF; i++){
                split($i,values,",")
                if((values[1]+values[2])>0) freq=values[1]/(values[1]+values[2])
                else freq=0
                printf freq"\t"
        }
        printf "\n"}' $FILE.frequencies > $FILE.frequencies2
mv $FILE.frequencies2 $FILE.frequencies

awk 'BEGIN{scaffold="";pos=0;newpos=0}
        {if($2==scaffold){newpos=pos+$3}else{scaffold=$2;pos=newpos};chpos=pos+$3;print $0,chpos}' \
        $FILE.frequencies > $FILE.frequencies.newpos

gzip $FILE.treemix.frq

## Run treemix
for i in {0..5}
do
 treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.$i -root Pennellii -bootstrap -k 500 -noss > treemix_${i}_log
done
