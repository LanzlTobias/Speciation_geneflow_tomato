#!/bin/bash
#############
##SNP calling
#############
for i in {1..12}
do
        chr=$(head -n$i chr.file | tail -n1)
        freebayes -0 -f reference.fasta -L bam.filelist --populations pop.File -r $chr  > $i.vcf &
done
wait

for i in {1..12}
do
        gzip chr$i.vcf
done

####################
##Start of filtering
####################
for i in {1..12} 
do
        #filt1
        vcftools --gzvcf chr$i.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out filt1
        #filt2
        vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" filt1.recode.vcf > filt2.vcf
        #filt3
        vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s filt2.vcf > filt3.vcf
        #filt4
        vcffilter -f "QUAL / DP > 0.25" filt3.vcf > filt4.vcf
        #filt5
                #filt5.step1. create a list of the depth of each locus
                cut -f8 filt4.vcf | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > filt4.DEPTH
                #filt5.step2. create a list of quality scores
                awk '!/#/' filt4.vcf | cut -f1,2,6 > filt4.vcf.loci.qual
                #filt5.step3. calculate the mean depth
                meanDP=`awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' filt4.DEPTH`
                #filt5.step4. mean plus 3X the square of the mean
                distDP=`python -c "print int($meanDP+3*($meanDP**0.5))"`
                #filt5.step5. paste the depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
                paste filt4.vcf.loci.qual filt4.DEPTH | awk -v x=$distDP '$4 > x' | awk '$3 < 2 * $4' > filt4.lowQDloci
                #filt5.step5. remove those sites and recalculate the depth across loci
                /data/proj/teaching/NGS_course/bin/vcftools --vcf filt4.vcf --site-depth --exclude-positions filt4.lowQDloci --out filt4
                #filt5.step6. cut output to only the depth scores
                cut -f3 filt4.ldepth > filt4.site.depth
                #filt5.step7. calculate the average depth by dividing the above file by the number of individuals (x=1)
                awk '!/D/' filt4.site.depth | awk -v x=1 '{print $1/x}' > meandepthpersite
                ##filt5 ## check the Depth histogram
                R --slave -e 'pdf("DPhist.chrom$i.pdf"); hist(read.table(file="meandepthpersite")[,1], xlim=c(0,300), breaks=100000); dev.off()'
                #filt5.step8. combine both filters above to produce new VCF file (remove all loci above a mean depth of 150 and loci that do not have quality scores 2 times the depth)
                vcftools --vcf filt4.vcf --recode --recode-INFO-all --out filt5 --max-meanDP 150 --exclude-positions filt4.lowQDloci
        #filt6
        vcftools --vcf filt5.recode.vcf --max-missing 0.95 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | gzip -c > solanum.$i.filtered.vcf.gz
        #separate phisically linked SNPs
        gunzip -c solanum.$i.filtered.vcf.gz | vcfallelicprimitives -k -g | gzip -c > filt6.vcf.gz
        #filt7
        vcftools --gzvcf filt6.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip -c > solanum.$i.goodVariants.vcf.gz
done

mv solanum.1.goodVariants.vcf.gz solanum.fb.snp.vcf.gz
for i in {2..12}
do
		zcat solanum.$i.goodVariants.vcf.gz | grep -v "^#" >> solanum.fb.snp.vcf.gz
		rm solanum.$i.goodVariants.vcf.gz
done
