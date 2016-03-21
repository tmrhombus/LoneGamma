#!/bin/bash

#voms-proxy-init --voms cms --valid 100:00

domc=true
dodata=false
dosubmit=false

START=$(date +%s);
printf "Started at `date`\n\n"

mkdir -p "${submitbase}/${version}/lists"
mkdir -p "${submitbase}/${version}/submit"

lumi=2240. # /pb

  #"ZLLG"
if [ ${domc} = true ]
then
 for mc_samplename in \
  "ZJetsToNuNu"
 do

  for htbin in \
   "100To200" \
   "200To400" \
   "400To600" \
   "600ToInf"
  do 

  #submitname="${mc_samplename}"
  submitname="${mc_samplename}${htbin}"

  initevents="${submitbase}/${version}/lists/initialEvents.txt"
  touch ${initevents}

 # count the total number of events, put in a file if it's not there
 if ! grep -F "${submitname} " ${initevents}
 then 
  echo "Need to count events.. making list"

#/hdfs/store/user/jjbuch/ZJetsToNuNu_HT-100To200_13TeV-madgraph/
#/hdfs/store/user/jjbuch/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph/
#/hdfs/store/user/jjbuch/ZJetsToNuNu_HT-200To400_13TeV-madgraph/
#/hdfs/store/user/jjbuch/ZJetsToNuNu_HT-400To600_13TeV-madgraph/

   #find /hdfs/store/user/jjbuch/${mc_samplename}*TuneCUETP8M1_13TeV*/crab_ggNtuplizer_spring15_${mc_samplename}*_try5/151212_*/0000/*root > \
   find /hdfs/store/user/jjbuch/${mc_samplename}*${htbin}*13TeV*/crab_ggNtuplizer_spring15_${mc_samplename}*_try5/151212_*/0000/*root > \
    ${submitbase}/${version}/lists/hdfslist_${submitname}.txt

   # format as xrootd
   cp ${submitbase}/${version}/lists/hdfslist_${submitname}.txt \
      ${submitbase}/${version}/lists/xrdlist_${submitname}.txt 
   xrdlist="${submitbase}/${version}/lists/xrdlist_${submitname}.txt"
   sed -i 's@/hdfs/@root://cmsxrootd.hep.wisc.edu//@g' $xrdlist #

   echo "                    .. counting events"
   python ./eventCounter.py ${submitname} \
      ${submitbase}/${version}/lists/hdfslist_${submitname}.txt \
      ${initevents}
   echo "                    .. counted events"
   
  fi
  xrdlist="${submitbase}/${version}/lists/xrdlist_${submitname}.txt"

  # sample specific parameters..
  treename="ggNtuplizer/EventTree"
  isMC="kTRUE"
  nrE="$(grep -P ${submitname} ${initevents} | sed -n -e "s@${submitname} Events: @@p")"
  xc="$(grep -P ${submitname}  ${initevents} | sed -n -e "s@${submitname} XC: @@p")"
   #printf "\nxc is ${xc}\n"

  # make correct executable xx_callpostAnalyzer_ZnunuG
  cp template_callpostAnalyzer_ZnunuG.cc    "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@SAMPLENAME@${submitname}@g"  "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@TREENAME@${treename}@g"      "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@ISMC@${isMC}@g"              "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@CROSSSEC@${xc}@g"            "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@NREVENTS@${nrE}@g"           "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@LUMI@${lumi}@g"              "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_ZnunuG.cc"

  printf "  done making submit template\n"

  if [ ${dosubmit} = true ]
  then
  printf " submitting\n"

  farmoutAnalysisJobs \
   --infer-cmssw-path \
   --fwklite \
   --input-file-list=${xrdlist} \
   --input-files-per-job=100 \
   --use-hdfs \
   --extra-inputs=${submitbase}/postAnalyzer_ZnunuG.C,${submitbase}/postAnalyzer_ZnunuG.h \
   ${version} \
   "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_ZnunuG.cc"

  fi # ${dosubmit} = true

  done # for htbin in ..
 done # for mc_samplename in ..
fi # ${domc} = true


######## Data ################################

if [ ${dodata} = true ]
then
 for data_samplename in \
  "SinglePhoton"

 do
  printf "\n${data_samplename}\n"
  
  printf " making list of files\n"

  find /hdfs/store/user/gomber/SinglePhoton_Crab_2015D_v3_226fb/SinglePhoton/crab_job_single_photon_13TeV_v3_226fb/160109_082344/0000/*root > \
   ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt
  find /hdfs/store/user/gomber/SinglePhoton_Crab_2015D_v4_226fb/SinglePhoton/crab_job_single_photon_13TeV_v4_226fb/160109_082133/000*/*root >> \
   ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt

  # format as xrootd
  cp ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt \
     ${submitbase}/${version}/lists/xrdlist_${data_samplename}.txt 
  xrdlist="${submitbase}/${version}/lists/xrdlist_${data_samplename}.txt"
  sed -i 's@/hdfs/@root://cmsxrootd.hep.wisc.edu//@g' $xrdlist #

  printf "  done making list of files\n"
  printf " making submit template\n"

  # sample specific parameters..
  treename="ggNtuplizer/EventTree"
  isMC="kFALSE"
  xc="1."
  nrE=${lumi}

  # make correct executable xx_callpostAnalyzer_ZnunuG
  cp template_callpostAnalyzer_ZnunuG.cc         "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@SAMPLENAME@${data_samplename}@g"  "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@TREENAME@${treename}@g"           "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@ISMC@${isMC}@g"                   "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@CROSSSEC@${xc}@g"                 "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@NREVENTS@${nrE}@g"                "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_ZnunuG.cc"
  sed -i "s@LUMI@${lumi}@g"                   "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_ZnunuG.cc"

  printf "  done making submit template\n"

  if [ ${dosubmit} = true ]
  then
  printf " submitting\n"

   farmoutAnalysisJobs \
    --infer-cmssw-path \
    --fwklite \
    --input-file-list=${xrdlist} \
    --input-files-per-job=100 \
    --use-hdfs \
    --extra-inputs=${submitbase}/postAnalyzer_ZnunuG.C,${submitbase}/postAnalyzer_ZnunuG.h \
    ${version} \
    "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_ZnunuG.cc"

  fi # ${dosubmit} = true

 done # for data_samplename in ..
fi # ${dodata}
