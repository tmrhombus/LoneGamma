#!/bin/bash

# "et" \
# "eta" \
# "pfMET" \
for var in \
 "sieieF5x5" 
do
 for ptrange in \
  "175to190" \
  "190to250" \
  "250to400" \
  "400to1000"
 do

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
