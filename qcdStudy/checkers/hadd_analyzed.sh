echo
echo "version: ${version}"

mkdir -p ${submitbase}/${version}/analyzed
mkdir -p ${submitbase}/${version}/plots

printf "${submitbase}/${version}/analyzed \n"
printf "${submitbase}/${version}/plots \n"

for samplename in \
 "GJets_HT40To100" \
 "GJets_HT100To200" \
 "GJets_HT200To400" \
 "GJets_HT400To600" \
 "GJets_HT600ToInf" \
 "SinglePhoton_2015D"

do
 hadd \
  ${submitbase}/${version}/analyzed/analyzed_${samplename}.root \
  ${hdfs}/${version}-${samplename}_callpostAnalyzer_QCD/*root
done

hadd \
 ${submitbase}/${version}/analyzed/analyzed_GJets_Merged.root \
 ${submitbase}/${version}/analyzed/analyzed_GJets_HT*root

