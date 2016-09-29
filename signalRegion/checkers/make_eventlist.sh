echo
echo "submitbase: ${submitbase}"
echo "version: ${version}"

mkdir -p ${submitbase}/gitignore/${version}/eventcomp
printf "${submitbase}/gitignore/${version}/eventcomp \n"

thedir="${submitbase}/gitignore/${version}/eventcomp"

## "ZLLG" \
## "ZJets100To200" \
## "ZJets200To400" \
## "ZJets400To600" \
## "ZJets600ToInf" \
for samplename in \
 "SinglePhotonData"
do
 #grep -h 'Total Passing RECO : ' /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*/*out > ${thedir}/events_${samplename}_v.txt

 grep -h '[0-9][0-9][0-9][0-9][0-9][0-9]:[0-9][0-9]*:[0-9]\|Total Passing RECO' /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*/*out > ${thedir}/events_${samplename}_w.txt

 #grep -h -A 2 '[0-9][0-9][0-9][0-9][0-9][0-9]:[0-9][0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*/*out > ${thedir}/events_${samplename}.txt


# sed -i 's/\s*$//' ${thedir}/events_${samplename}.txt
# 
# grep -v -F -x -f eventlist_400.txt ${thedir}/events_${samplename}.txt > ${thedir}/events_menothey.txt
# grep -F -x -f eventlist_400.txt ${thedir}/events_${samplename}.txt > ${thedir}/events_common.txt
# grep -v -F -x -f ${thedir}/events_${samplename}.txt eventlist_400.txt > ${thedir}/events_theynome.txt
   
done # for samplename in "GJets_*" "SingleP"

