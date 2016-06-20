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
 "ZnnGJets" \
 "ZllGJets" \
 "ZllJetsHT100to200" \
 "ZllJetsHT200to400" \
 "ZllJetsHT400to600" \
 "ZllJetsHT600toInf" \
 "SinglePhoton"
do

 # analyzed
 if [ "$countanalyzed" = "True" ]
 then
   nfs_analyzed=$(ls -d /nfs_scratch/tperry/${version}*${samplename}_callpostAnalyzer_ZnunuG/*/ | wc -l)
   hdfs_analyzed=$(ls -1 /hdfs/store/user/tperry/${version}*${samplename}_callpostAnalyzer_ZnunuG/*.root | wc -l) 
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
#### grep -l "exited with status 0" /nfs_scratch/tperry/Earth_DataA_8TeVEle-Ele-PATData/*/*out > dataAEle_good.txt
#### ls /nfs_scratch/tperry/Earth_DataA_8TeVEle-Ele-PATData/*/*out > dataAEle_all.txt
#### grep -Fvf good_dataAEle.txt all_dataAEle.txt > dataAEle_bad.txt
