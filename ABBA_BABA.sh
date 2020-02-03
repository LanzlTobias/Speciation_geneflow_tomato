#!/bin/bash
## Script for running the ABBA_BABA test
#
#
## Arguments for the script:
#First argument is the VCF file (gzipped or not)
#Second argument is the file containing the paths to all poporder-files. Each path has to be given in a new line.
#Example of a poporder-file is given as pop_order.txt
#Third argument is the populations.txt. Example is given as pop_file.txt
#Fourth argument is optional. It is the number on which poporder it should be started. For multi-threading.
#Fifth argument is optional. It is the number on which poporder it should end. For multi-threading.
#Sixth argument is optional. It is the path to the directory containing the scripts needed for the pipeline.
#parseVCF.py available at: https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing
#ABBA_BABA.v1.pl available at: https://github.com/owensgl/abba_baba
#ABBA_out_blocker_5Mb.pl available at: https://github.com/owensgl/abba_baba
#ABBA_out_blocker.pl available at: https://github.com/owensgl/abba_baba
#convertDoGeno2egglib.pl available at: https://github.com/owensgl/abba_baba
#countbadcolumns.pl available at: https://github.com/owensgl/pop_gen
#Jackknife_ABBA_pipe.R available at: https://github.com/owensgl/abba_baba
#ABBA_pvalue.R available at: https://github.com/owensgl/abba_baba

FILE=${1%.gz}
FILE=${file%.vcf}

if [ -z "$3" ]
  then
    script_dir=$pwd
  else
    script_dir=$3
fi

if test "$5" == "1"
  then

python $script_dir/parseVCF.py -i $FILE.vcf.gz > $FILE.geno

sed -i 's/|//g' $FILE.geno

sed -i '1d' $FILE.geno

perl $script_dir/convertDoGeno2egglib.pl $FILE.geno > $FILE.geno.txt

cat header.txt $FILE.geno.txt > $FILE"_WH.geno"

rm $FILE.geno.txt $FILE.geno
fi

echo "############## Start ABBA BABA############"
END=$(cat $2 | wc -l)

if [ -z "$5" ]
  then
    START=1
  else
    START=$5
fi

if [ -z "$6" ]
  then
   END=$(cat $2 | wc -l)
  else
   END=$6
fi
    

for (( i=$START; i<=$END; i++ ))
do
poporder=$(head -n $i $2 | tail -n 1)
        perl $script_dir/ABBA_BABA.v1.pl $FILE"_WH.geno" $3 $poporder > D_fd_Fhom_ANGSD_$i.txt

perl $script_dir/ABBA_out_blocker_5Mb.pl D_fd_Fhom_ANGSD_$i.txt > D_fd_Fhom_ANGSD_block_$i.txt

Rscript $script_dir/Jackknife_ABBA_pipe.R D_fd_Fhom_ANGSD_block_$i.txt D_fd_Fhom_ANGSD_jackknife_$i.txt

tail -n 3 D_fd_Fhom_ANGSD_$i.txt > tmp1
tail -n 3 D_fd_Fhom_ANGSD_jackknife_$i.txt > tmp2

paste tmp1 tmp2 > D_fd_Fhom_ANGSD_2_$i.txt
rm tmp1
rm tmp2

Rscript $script_dir/ABBA_pvalue.R D_fd_Fhom_ANGSD_2_$i.txt D_fd_Fhom_ANGSD_final_$i.txt
done

echo "############# Finished ABBA BABA #########"

#Summarize the D and Z-score for all the combinations in Zscore.txt
echo -e "Poporder_Nr\tCombination\tMean_D\tStd_D\tZ_score\tp_value" > Zscore.txt
for (( i=$START; i<=$END; i++ ))
do
        echo "$i" > number
        poporder=$(head -n $2 | tail -n)
cut -f1 $poporder | tr '\n' '\t' | sed 's/\t$//' | tr '\t' _ > combination
cut -f3-6 D_fd_Fhom_ANGSD_final_$i.txt | head -n1 > values
paste number combination values >> sample

cat sample >> Zscore_$START"_"$END.txt
done
rm combination
rm values
rm number
rm sample
