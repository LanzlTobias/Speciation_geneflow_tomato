#!/bin/bash

./H3toH2.sh
./H3toH2H1.sh
./H3toH2-lm.sh
./H3toH2-lm-early.sh
./H3toH2-admix.sh
./H3toH2-admix-early.sh
./H3toH2-cm.sh
./nomig.sh
./nomig-div.sh

echo -e "Simulation\tMean_D\tStd_D\tZ_score\tp-value" > header
for i in $(ls *final*txt)
do
        name=$(echo $i | cut -f6,7 -d"_"| cut -f1 -d".")
        echo $name > Simulation
cut -f3-6 $i | head -n1 > values
paste Simulation values >> sample
rm Simulation
rm values
done
cat header sample > Zscore_sim.txt
rm header
rm sample

tar -czf D_sim_files.tar.gz D*txt
rm D*txt
