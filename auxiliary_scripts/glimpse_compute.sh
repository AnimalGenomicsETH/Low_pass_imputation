#!/usr/bin/env bash

#Go to the root directory of the project and use as ./glimpse_compute.sh > <output>

echo "panel content replicate coverage CPU memory"

for p in BSW_panel non_BSW_panel multibreed_panel
do
  for i in low_pass_hap_panels/GLIMPSE/${p}/log_folder/impute_phase/*log
  do
     echo $(basename ${i%_*.log}) $(awk '/CPU time|Max Memory/ {print $4}' $i)
  done | awk -v P=${p} '{C[$1]+=$2; if ($3>M[$1]) {M[$1]=$3}} END {for (key in C) {print P,key,C[key],M[key]}}' | sed -E 's/(_shuffle|impute_phase_|_samples|_panel)//g' | sed 's/_/ /g' | sed 's/non /non_/g'
done
