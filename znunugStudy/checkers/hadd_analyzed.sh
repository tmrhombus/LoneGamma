echo
echo "version: ${version}"

mkdir -p ${submitbase}/${version}/analyzed
mkdir -p ${submitbase}/${version}/plots

printf "${submitbase}/${version}/analyzed \n"
printf "${submitbase}/${version}/plots \n"

for samplename in \
 "SinglePhoton"

do
 hadd \
  ${submitbase}/${version}/analyzed/analyzed_${samplename}.root \
  ${hdfs}/${version}-${samplename}_callpostAnalyzer_ZnunuG/*root
done

