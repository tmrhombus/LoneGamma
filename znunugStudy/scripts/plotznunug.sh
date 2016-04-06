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
for var in \
 "dilep_mass" 
do
  #"175to1000"
 for ptrange in \
  "175to1000" 
 do

  sf=2320./2240.
  echo $sff

  python ${submitbase}/scripts/plotznunug.py \
   --version "${version}" \
   --inpdir "${analyzed}" \
   --mc_filename   "analyzed_ZLLG.root" \
   --nc_filename   "analyzed_ZnnJ.root" \
   --data_filename "analyzed_SinglePhoton.root" \
   --outdir "${plots}" \
   --scalefactor ${sf} \
   --out_filename "plotZnunuG_${var}_${ptrange}" \
   --variable "${var}" \
   --ptrange "${ptrange}" #\
   #--do_log

 done # for ptrange in "75to100" ..
done # for var in "et" "eta" ..
