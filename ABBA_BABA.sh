

python parseVCF.py -i solanum.fb.vcf.gz > solanum.fb.geno

sed -i 's/|//g' solanum.fb.geno

sed -i '1d' solanum.fb.geno

perl convertDoGeno2egglib.pl solanum.fb.geno > solanum.fb.geno.txt


cat header.txt solanum.fb.geno.txt > $OUTPUT/solanum.fb_WH.geno

mv *pl $OUTPUT

cd $OUTPUT
pwd
ls

echo "############## Start ABBA BABA############"

for i in {1..20}
do
        perl ABBA_BABA.v1.pl solanum.fb_WH.geno ../Population_file2.txt ../poporder/Poporder_$i.txt > D_fd_Fhom_ANGSD_$i.txt

perl ABBA_out_blocker_5Mb.pl D_fd_Fhom_ANGSD_$i.txt > D_fd_Fhom_ANGSD_block_$i.txt

Rscript ../Jackknife_ABBA_pipe2.R D_fd_Fhom_ANGSD_block_$i.txt D_fd_Fhom_ANGSD_jackknife_$i.txt

tail -n 3 D_fd_Fhom_ANGSD_$i.txt > tmp1
tail -n 3 D_fd_Fhom_ANGSD_jackknife_$i.txt > tmp2

paste tmp1 tmp2 > D_fd_Fhom_ANGSD_2_$i.txt
rm tmp1
rm tmp2

Rscript ../ABBA_pvalue.R D_fd_Fhom_ANGSD_2_$i.txt D_fd_Fhom_ANGSD_final_$i.txt
done

echo "############# Finished ABBA BABA #########"

../extract_zscore.sh
mkdir site_windows
cd site_windows
