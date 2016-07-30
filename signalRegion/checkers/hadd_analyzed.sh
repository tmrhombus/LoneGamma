echo
echo "version: ${version}"

mkdir -p ${submitbase}/gitignore/${version}/analyzed
mkdir -p ${submitbase}/gitignore/${version}/plots

printf "${submitbase}/gitignore/${version}/analyzed \n"
printf "${submitbase}/gitignore/${version}/plots \n"

for samplename in \
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
 hadd \
  ${submitbase}/gitignore/${version}/analyzed/analyzed_${samplename}.root \
  ${hdfs}/${version}-${samplename}_callpostAnalyzer_Signal/*root
done

hadd \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_GJets.root \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_GJetsHT*root

hadd \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_QCD.root \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_QCDPt*root

hadd \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_ZllJets.root \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_ZllJetsHT*root

hadd \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_ZnnJets.root \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_ZnnJetsHT*root


