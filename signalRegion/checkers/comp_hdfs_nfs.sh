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

for samplename in \
 "SinglePhoton"      \
 "GJetsHT40to100"    \
 "GJetsHT100to200"   \
 "GJetsHT200to400"   \
 "GJetsHT400to600"   \
 "GJetsHT600toInf"   \
 "QCDPt15to20"       \
 "QCDPt20to30"       \
 "QCDPt30to50"       \
 "QCDPt50to80"       \
 "QCDPt80to120"      \
 "QCDPt120to170"     \
 "QCDPt170to300"     \
 "QCDPt300toInf"     \
 "ZllGJets"          \
 "ZnnGJets"          \
 "ZllJetsHT100to200" \
 "ZllJetsHT200to400" \
 "ZllJetsHT400to600" \
 "ZllJetsHT600toInf" \
 "ZnnJetsHT100to200" \
 "ZnnJetsHT200to400" \
 "ZnnJetsHT400to600" \
 "ZnnJetsHT600toInf" \
 "WlnGJets"          \
 "Wen"               \
 "Wmn"               \
 "Wtn"               \
 "TTGJets"           
do

 # analyzed
 if [ "$countanalyzed" = "True" ]
 then
   nfs_analyzed=$(ls -d /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*/ | wc -l)
   hdfs_analyzed=$(ls -1 /hdfs/store/user/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*.root | wc -l) 
   analyzed_diff=$(($nfs_analyzed-$hdfs_analyzed))
   hdfstot=$((${hdfstot}+${hdfs_analyzed}))
   nfstot=$((${nfstot}+${nfs_analyzed}))
  printf "  %20s %6s %6s %s \n" "${samplename}" "${nfs_analyzed}" "${hdfs_analyzed}" "${analyzed_diff} incomplete analyzers"
  #printf "  ${samplename}\t ${nfs_analyzed}\t ${hdfs_analyzed}\t ${analyzed_diff} incomplete analyzers\n"
 fi
done # for samplename in "GJets_*" "SingleP"

printf "%s " " -----------------------------------------"
printf "%s \n" "-----------------------------------------"
diftot=$((${nfstot}-${hdfstot}))
printf "  Total: nfs: ${nfstot}  hdfs:  ${hdfstot}  left:  ${diftot}\n\n"


##### for checking individual files
#### grep -l "exited with status 0" /nfs_scratch/tperry/Earth_DataA_8TeVEle-Ele-PATData/*/*out > dataAEle_good.txt
#### ls /nfs_scratch/tperry/Earth_DataA_8TeVEle-Ele-PATData/*/*out > dataAEle_all.txt
#### grep -Fvf good_dataAEle.txt all_dataAEle.txt > dataAEle_bad.txt


