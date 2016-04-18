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

 for cut in \
  "a_idnc" \
  "b_idnc_mL30" \
  "c_idnc_trig" \
  "d_idnc_mL30_trig" \
  "e_idnc_t175" \
  "f_idnc_mL30_t175" \
  "g_idnc_t250" \
  "h_idnc_mL30_t250" \
  "i_idnc_mL30_allt" 
 do


  python ${submitbase}/scripts/plotqcd.py \
   --version "${version}" \
   --inpdir "${analyzed}" \
   --mc_filename   "mrg4bins_GJets.root" \
   --qcd_filename  "mrg4bins_QCD.root" \
   --data_filename "mrg4bins_SinglePhoton.root" \
   --data_ele_filename "mrg4bins_DoubleElectron.root" \
   --outdir "${plots}" \
   --out_filename "plotqcd4bins_${var}_${ptrange}" \
   --variable "${var}" \
   --cut "${cut}" \
   --ptrange "${ptrange}" \
   --do_log

  done # for cut in
 done # for ptrange in "75to100" ..
done # for var in "et" "eta" ..
