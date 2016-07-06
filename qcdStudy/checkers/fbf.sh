echo
echo "submitbase: ${submitbase}"
echo "version: ${version}"

nfstot=0
hdfstot=0

tabs 15

####
# Count 
####

 #"QCDPt15to20" \
 #"QCDPt20to30" \
 #"QCDPt30to50" \
 #"QCDPt50to80" \
 #"DoubleElectron" 
 #"QCDPt170to300" \
 #"GJetsHT40to100" \
 #"GJetsHT100to200" \
 #"GJetsHT200to400" \
 #"GJetsHT400to600" \
 #"GJetsHT600toInf" \
 #"SinglePhoton" 
for samplename in \
 "QCDPt80to120" \
 "QCDPt120to170" \
 "QCDPt300toInf" 
do

##### for checking individual files
 grep -l "exited with status 0" /nfs_scratch/tperry/${version}-${samplename}*/*/*out > ${samplename}_good.txt
 ls /nfs_scratch/tperry/${version}-${samplename}*/*/*out > ${samplename}_all.txt
 grep -Fvf ${samplename}_good.txt ${samplename}_all.txt > ${samplename}_bad.txt
done # for samplename in "GJets_*" "SingleP"

