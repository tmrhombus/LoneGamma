echo
echo "version: ${version}"

mkdir -p ${submitbase}/gitignore/${version}/analyzed
mkdir -p ${submitbase}/gitignore/${version}/plots

printf "${submitbase}/gitignore/${version}/analyzed \n"
printf "${submitbase}/gitignore/${version}/plots \n"

for samplename in \
 "ZnnGJets" \
 "ZllGJets" \
 "ZllJetsHT100to200" \
 "ZllJetsHT200to400" \
 "ZllJetsHT400to600" \
 "ZllJetsHT600toInf" \
 "SinglePhoton"

do
 hadd \
  ${submitbase}/gitignore/${version}/analyzed/analyzed_${samplename}.root \
  ${hdfs}/${version}-${samplename}_callpostAnalyzer_ZnunuG/*root
done

hadd \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_ZllJets.root \
 ${submitbase}/gitignore/${version}/analyzed/analyzed_ZllJetsHT*.root \
