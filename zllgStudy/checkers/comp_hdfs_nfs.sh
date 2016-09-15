echo
echo "submitbase: ${submitbase}"
echo "version: ${version}"

countanalyzed="True"

nfstot=0
hdfstot=0

tabs 20

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
  "WlnGJets"          \
  "TTGJets" 
do

 # analyzed
 if [ "$countanalyzed" = "True" ]
 then
   nfs_analyzed=$(ls -d /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_ZnunuG/*/ | wc -l)
   hdfs_analyzed=$(ls -1 /hdfs/store/user/tperry/${version}*${samplename}_callpostAnalyzer_ZnunuG/*callpostAnalyzer*.root | wc -l) 
   analyzed_diff=$(($nfs_analyzed-$hdfs_analyzed))
   hdfstot=$((${hdfstot}+${hdfs_analyzed}))
   nfstot=$((${nfstot}+${nfs_analyzed}))
  #printf "  ${samplename}\t ${nfs_analyzed}\t ${hdfs_analyzed}\t ${analyzed_diff} incomplete analyzers\n"
  printf "  %20s %20s %20s %20s incomplete analyzers\n" ${samplename} ${nfs_analyzed} ${hdfs_analyzed} ${analyzed_diff}
 fi
done # for samplename in "GJets_*" "SingleP"

printf "%s " " -----------------------------------------"
printf "%s \n" "-----------------------------------------"
diftot=$((${nfstot}-${hdfstot}))
printf "  Total: nfs: ${nfstot}  hdfs:  ${hdfstot}  left:  ${diftot}\n\n"


##### for checking individual files
# grep -l "exited with status 0" /nfs_scratch/tperry/jaumesse-*/*/*out > dataAEle_good.txt
# ls /nfs_scratch/tperry/jaumesse-*/*/*out > dataAEle_all.txt
# grep -Fvf dataAEle_good.txt dataAEle_all.txt > dataAEle_bad.txt
