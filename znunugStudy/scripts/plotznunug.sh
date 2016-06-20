#!/bin/bash

# "et" \
# "eta" \
# "pfMET" \
#  "175to190" \
#  "190to250" \
#  "250to400" \
#  "400to700" \
#  "175to190" \
#  "190to250" \
#  "250to400" \
#  "400to700" \
#  "700to1000" 
 #"dilep_mass" 
 #"dimu_mass" \
 #"diele_mass"
for var in \
 "dilep_mass" 

do
  #"175to1000"
 for ptrange in \
  "175to1000" 
 do

   # "m170_ywnd_ydphi" \
   # "m170_ywnd_ndphi" \
   # "m170_nwnd_ydphi" \
   # "m170_nwnd_ndphi" \
   # "m110_ywnd_ydphi" \
   # "m110_ywnd_ndphi" \
   # "m110_nwnd_ydphi" \
   # "m110_nwnd_ndphi" 
  for sel in \
   "m170_nwnd_ydphi"
  do  

  python -b ${submitbase}/scripts/plotznunug.py \
   --version "${version}" \
   --inpdir "${analyzed}" \
   --mc_filename   "analyzed_ZllGJets.root" \
   --nc_filename   "analyzed_ZnnGJets.root" \
   --lc_filename   "analyzed_ZllJets.root" \
   --data_filename "analyzed_SinglePhoton.root" \
   --outdir "${plots}" \
   --out_filename "plotZnunuG_${var}_${ptrange}_${sel}" \
   --variable "${var}" \
   --ptrange "${ptrange}" \
   --selection "${sel}"
   #--do_log

  done # for sel in 
 done # for ptrange in "75to100" ..
done # for var in "et" "eta" ..
