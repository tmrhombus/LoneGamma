echo
echo "submitbase: ${submitbase}"
echo "version: ${version}"

mkdir -p ${submitbase}/${version}/eventcomp
printf "${submitbase}/${version}/eventcomp \n"

thedir="${submitbase}/${version}/eventcomp"

## "ZLLG" \
## "ZJets100To200" \
## "ZJets200To400" \
## "ZJets400To600" \
## "ZJets600ToInf" \
for samplename in \
 "SinglePhoton"
do
 for var in \
  n_initial \
  n_trig    \
  n_phokin  \
  n_spike   \
  n_met     \
  n_metfilt \
  n_dphipho \
  n_lepveto \
  n_dphijet
 
 do
   grep -h ${var} /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*/*out > ${thedir}/${var}_${samplename}.txt
 done # for var in ...
done # for samplename in "GJets_*" "SingleP"

#sed -i 's/\s*$//' ${thedir}/events_${samplename}.txt
#
#grep -v -F -x -f ../EventList.txt ${thedir}/events_SinglePhoton.txt > ${thedir}/events_menothey.txt
#grep -F -x -f ../EventList.txt ${thedir}/events_SinglePhoton.txt > ${thedir}/events_common.txt
#grep -v -F -x -f ${thedir}/events_SinglePhoton.txt ../EventList.txt > ${thedir}/events_theynome.txt
#
#
# n_passSpike      3345
# n_passLepRej     3330
# n_passMETfilters 3339
# n_passMET        11  
# n_passdPhiPhoMET 1114
# n_passdPhiJM     2346
# n_initial 
# n_trig    
# n_phokin  
# n_spike   
# n_met     
# n_metfilt 
# n_dphipho 
# n_lepveto 
# n_dphijet 
