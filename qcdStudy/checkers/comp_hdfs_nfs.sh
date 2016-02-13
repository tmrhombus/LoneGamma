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
 "QCD_Pt15to20" \
 "QCD_Pt20to30" \
 "QCD_Pt30to50" \
 "QCD_Pt50to80" \
 "QCD_Pt80to120" \
 "QCD_Pt120to170" \
 "QCD_Pt170to300" \
 "QCD_Pt300toInf" \
 "GJets_HT40To100" \
 "GJets_HT100To200" \
 "GJets_HT200To400" \
 "GJets_HT400To600" \
 "GJets_HT600ToInf" \
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


##### for checking individual files
#### grep -l "exited with status 0" /nfs_scratch/tperry/Earth_DataA_8TeVEle-Ele-PATData/*/*out > dataAEle_good.txt
#### ls /nfs_scratch/tperry/Earth_DataA_8TeVEle-Ele-PATData/*/*out > dataAEle_all.txt
#### grep -Fvf good_dataAEle.txt all_dataAEle.txt > dataAEle_bad.txt
