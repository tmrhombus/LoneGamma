echo
echo "version: ${version}"

mkdir -p ${submitbase}/${version}/analyzed
mkdir -p ${submitbase}/${version}/plots

printf "${submitbase}/${version}/analyzed \n"
printf "${submitbase}/${version}/plots \n"

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
 hadd \
  ${submitbase}/${version}/analyzed/analyzed_${samplename}.root \
  ${hdfs}/${version}-${samplename}_callpostAnalyzer_QCD/*root
done

hadd \
 ${submitbase}/${version}/analyzed/analyzed_GJets_Merged.root \
 ${submitbase}/${version}/analyzed/analyzed_GJets_HT*root

hadd \
 ${submitbase}/${version}/analyzed/analyzed_QCD_Merged.root \
 ${submitbase}/${version}/analyzed/analyzed_QCD_Pt*root

