#!/bin/sh

# Information Printer

 purpath="http://www.hep.wisc.edu/~tperry/MonoPhoton/qcdPlots/${version}/purity"

   printf  " a_idnc           : Photon ID (pT > 175)           \n"  
   printf  " b_idnc_mL30      : ID + MET < 30                  \n"   
   printf  " c_idnc_trig      : ID + HLT_Photon165_HE10_v (12) \n"       
   printf  " d_idnc_mL30_trig : ID + MET + T165                \n"    
   printf  " e_idnc_t175      : ID + HLT_Photon175_v (7)       \n"       
   printf  " f_idnc_mL30_t175 : ID + MET + T175                \n"       
   printf  " g_idnc_t250      : ID + HLT_Photon250_NoHE_v      \n"       
   printf  " h_idnc_mL30_t250 : ID + MET + T250                \n"        
   printf  " i_idnc_mL30_allt : ID + MET + T165 + T175 + T250  \n"         

for plot in \
 "Histo"  \
 "Purity" 
do
 
 for ptrange in \
  "175to190"  \
  "190to250"  \
  "250to400"  \
  "400to700"  \
  "700to1000" \
  "175to1000" 
 
 do
 
  for var in \
  "chiso"     
 # "sieieF5x5" \
 # "sieipF5x5" \
 # "sipipF5x5" \
 # "rho"       \
 # "nVtx"      \
 # "phoet"     \
 # "pfMET"     \
 
  do
 
   for cut in \
    "a_idnc"           
 #   "b_idnc_mL30"      \  
 #   "c_idnc_trig"      \
 #   "d_idnc_mL30_trig" \ 
 #   "e_idnc_t175"      \  
 #   "f_idnc_mL30_t175" \ 
 #   "g_idnc_t250"      \  
 #   "h_idnc_mL30_t250" \ 
 #   "i_idnc_mL30_allt" \
 
   do
 
    printf  "${purpath}/png/${ptrange}/${var}/${cut}/${plot}_${var}_${ptrange}_${cut}.png  \n"   

   done # for plot
  done # for cut 
 done # for var 
done # for ptrange





