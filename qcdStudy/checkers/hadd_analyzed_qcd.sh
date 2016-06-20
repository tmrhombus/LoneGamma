echo
echo "version: ${version}"

mkdir -p ${submitbase}/gitignore/${version}/analyzed
mkdir -p ${submitbase}/gitignore/${version}/plots

printf "${submitbase}/gitignore/${version}/analyzed \n"
printf "${submitbase}/gitignore/${version}/plots \n"

## "QCDPt15to20" \
## "QCDPt20to30" \
## "QCDPt30to50" \
## "QCDPt50to80" \
#for samplename in \
# "QCDPt80to120" \
# "QCDPt120to170" \
# "QCDPt170to300" \
# "QCDPt300toInf" \
# "GJetsHT40to100" \
# "GJetsHT100to200" \
# "GJetsHT200to400" \
# "GJetsHT400to600" \
# "GJetsHT600toInf" \
# "SinglePhoton" \
# "DoubleElectron"
#
#do
# hadd \
#  ${submitbase}/gitignore/${version}/analyzed/analyzed_${samplename}.root \
#  ${hdfs}/${version}-${samplename}_callpostAnalyzer_QCD/*root
#done

hadd \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_GJetsMerged.root \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_GJetsHT*root

hadd \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_QCDMerged.root \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_QCDPt*root

