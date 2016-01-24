#!/bin/bash

# "et" \
# "eta" \
# "pfMET" \
for var in \
 "sieieF5x5" 
do
#  "75to100" \
#  "100to125" \
#  "125to145" \
#  "145to155" \
#  "155to165" \
#  "165to175" \
 for ptrange in \
  "175to190" \
  "190to250" \
  "250to400" \
  "400to1000"
 do
#  python ${submitbase}/scripts/plotqcd.py \
#   --version "${version}" \
#   --inpdir "${submitbase}/${version}/analyzed" \
#   --mc_filename "analyzed_GJets_Merged.root" \
#   --data_filename "analyzed_SinglePhoton_2015D.root" \
#   --outdir "${submitbase}/${version}/plots" \
#   --out_filename "plotqcd_${var}_${ptrange}" \
#   --variable "${var}" \
#   --ptrange "${ptrange}" 

  python ${submitbase}/scripts/plotqcd.py \
   --version "${version}" \
   --inpdir "${analyzed}" \
   --mc_filename   "mrg4bins_GJets.root" \
   --data_filename "mrg4bins_SinglePhoton.root" \
   --outdir "${plots}" \
   --out_filename "plotqcd4bins_${var}_${ptrange}" \
   --variable "${var}" \
   --ptrange "${ptrange}" \
   --do_log

 done # for ptrange in "75to100" ..
done # for var in "et" "eta" ..
