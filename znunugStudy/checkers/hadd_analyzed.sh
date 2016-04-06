echo
echo "version: ${version}"

mkdir -p ${submitbase}/${version}/analyzed
mkdir -p ${submitbase}/${version}/plots

printf "${submitbase}/${version}/analyzed \n"
printf "${submitbase}/${version}/plots \n"

for samplename in \
 "ZLLG" \
 "ZJetsToNuNu100To200" \
 "ZJetsToNuNu200To400" \
 "ZJetsToNuNu400To600" \
 "ZJetsToNuNu600ToInf" \
 "SinglePhoton"

do
 hadd \
  ${submitbase}/${version}/analyzed/analyzed_${samplename}.root \
  ${hdfs}/${version}-${samplename}_callpostAnalyzer_ZnunuG/*root
done

hadd \
 ${submitbase}/${version}/analyzed/analyzed_ZnnJ.root \
 ${submitbase}/${version}/analyzed/analyzed_ZJetsToNuNu*.root \
