#!/bin/bash

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
