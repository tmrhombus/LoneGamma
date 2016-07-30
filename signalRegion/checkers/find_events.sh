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
 for event in \
  "295505278" \
  "285078441" \
  "46661859" \
  "389981922" \
  "2208621600" \
  "2127047219" \
  "1215728566" \
  "85380629" \
  "89128679" \
  "776377515"

 do
   grep  ${event} /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*/*out >> ${thedir}/events_${samplename}_bad.txt
 done # for event in ... 
done # for samplename in "GJets_*" "SingleP"

#sed -i 's/\s*$//' ${thedir}/events_${samplename}.txt
#
#grep -v -F -x -f ../EventList.txt ${thedir}/events_SinglePhoton.txt > ${thedir}/events_menothey.txt
#grep -F -x -f ../EventList.txt ${thedir}/events_SinglePhoton.txt > ${thedir}/events_common.txt
#grep -v -F -x -f ${thedir}/events_SinglePhoton.txt ../EventList.txt > ${thedir}/events_theynome.txt


