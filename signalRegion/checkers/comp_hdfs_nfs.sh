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
  "SinglePhotonData"  \
  "SinglePhotonEle"   \
  "SinglePhotonJet"   \
  "SinglePhotonHalo"  \
  "SinglePhotonSpike" \
  "GJetsHT40to100"    \
  "GJetsHT100to200"   \
  "GJetsHT200to400"   \
  "GJetsHT400to600"   \
  "GJetsHT600toInf"   \
  "ZllGJets"          \
  "ZnnGJets"          \
  "ZllJetsHT100to200" \
  "ZllJetsHT200to400" \
  "ZllJetsHT400to600" \
  "ZllJetsHT600toInf" \
  "WlnGJets"          \
  "Wmn"               \
  "Wtn"               \
  "TTGJets"           \
  "TGJets"            \
  "GGJets"            \
  "WWG"               \
  "WZ"                \
  "ZZ"

do

 # analyzed
 if [ "$countanalyzed" = "True" ]
 then
   nfs_analyzed=$(ls -d /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*/ | wc -l)
   hdfs_analyzed=$(ls -1 /hdfs/store/user/tperry/${version}*${samplename}_callpostAnalyzer_Signal/*callpostAnalyzer*.root | wc -l) 
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

#hadd \
# ${submitbase}/gitignore/${version}/analyzed/analyzed_ZllJets.root \
# ${submitbase}/gitignore/${version}/analyzed/analyzed_ZllJetsHT*.root 
#
#hadd \
# ${submitbase}/gitignore/${version}/analyzed/analyzed_GJets.root \
# ${submitbase}/gitignore/${version}/analyzed/analyzed_GJetsHT*.root 



##### for checking individual files
#### grep -l "exited with status 0" /nfs_scratch/tperry/Earth_DataA_8TeVEle-Ele-PATData/*/*out > dataAEle_good.txt
#### ls /nfs_scratch/tperry/Earth_DataA_8TeVEle-Ele-PATData/*/*out > dataAEle_all.txt
#### grep -Fvf good_dataAEle.txt all_dataAEle.txt > dataAEle_bad.txt
