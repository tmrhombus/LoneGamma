echo
echo "submitbase: ${submitbase}"
echo "version: ${version}"

countanalyzed="True"

nfstot=0
hdfstot=0

tabs 15

####
# Count 
####

printf " \t\t /nfs \t /hdfs \t stauts \n"
printf "%s " " -----------------------------------------"
printf "%s \n" "-----------------------------------------"

 #"QCDPt15to20" \
 #"QCDPt20to30" \
 #"QCDPt30to50" \
 #"QCDPt50to80" \
 #"DoubleElectron" 
for samplename in \
 "QCDPt80to120" \
 "QCDPt120to170" \
 "QCDPt170to300" \
 "QCDPt300toInf" \
 "GJetsHT40to100" \
 "GJetsHT100to200" \
 "GJetsHT200to400" \
 "GJetsHT400to600" \
 "GJetsHT600toInf" \
 "SinglePhoton" 
do

 # analyzed
 if [ "$countanalyzed" = "True" ]
 then
   nfs_analyzed=$(ls -d /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_QCD/*/ | wc -l)
   hdfs_analyzed=$(ls -1 /hdfs/store/user/tperry/${version}*${samplename}_callpostAnalyzer_QCD/*.root | wc -l) 
   analyzed_diff=$(($nfs_analyzed-$hdfs_analyzed))
   hdfstot=$((${hdfstot}+${hdfs_analyzed}))
   nfstot=$((${nfstot}+${nfs_analyzed}))
  printf "  ${samplename}\t ${nfs_analyzed}\t ${hdfs_analyzed}\t ${analyzed_diff} incomplete analyzers\n"
 fi
done # for samplename in "GJets_*" "SingleP"

printf "%s " " -----------------------------------------"
printf "%s \n" "-----------------------------------------"
diftot=$((${nfstot}-${hdfstot}))
printf "  Total: nfs: ${nfstot}  hdfs:  ${hdfstot}  left:  ${diftot}\n\n"

