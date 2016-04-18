#!/bin/sh

# make directory structure

# webdir

 purpath="/afs/hep.wisc.edu/home/tperry/www/MonoPhoton/qcdPlots/${version}/purity"

for ptrange in \
 "175to190" \
 "190to250" \
 "250to400" \
 "400to700" \
 "700to1000" \
 "175to1000" 
do

 for var in \
 "chiso"     \
 "sieieF5x5" \
 "sieipF5x5" \
 "sipipF5x5" \
 "rho" \
 "nVtx" \
 "phoet" \
 "pfMET" 
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
 
   mkdir -p "${purpath}/pdf/${ptrange}/${var}/${cut}/lines" 
   mkdir -p "${purpath}/png/${ptrange}/${var}/${cut}/lines" 

  done # for cut 
 done # for var 
done # for ptrange





