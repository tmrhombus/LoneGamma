#!/bin/bash

#voms-proxy-init --voms cms --valid 100:00

domc=true
dodata=true
dosubmit=true

START=$(date +%s);
printf "Started at `date`"

mkdir -p "${submitbase}/${version}/lists"
mkdir -p "${submitbase}/${version}/submit"

# lumi = 

if [ ${domc} = true ]
then
 for mc_samplename in \
  "GJets_HT"  
 do

  for htbin in \
   "40To100" \
   "100To200" \
   "200To400" \
   "400To600" \
   "600ToInf" 
  do

  submitname="${mc_samplename}${htbin}"

  printf "${submitname}\n"
  
  # hdfs/store/user/jjbuch/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-40To100_25ns/151101_093834/0000/*root
  # hdfs/store/user/jjbuch/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-100To200_25ns/151101_092459/0000/*root
  # hdfs/store/user/jjbuch/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-200To400_25ns/151101_092710/0000/*root
  # hdfs/store/user/jjbuch/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-400To600_25ns/151101_093518/0000/*root
  # hdfs/store/user/jjbuch/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-600ToInf_25ns/151101_094056/0000/*root
  printf " making list of files\n"
  find /hdfs/store/user/jjbuch/${mc_samplename}-${htbin}_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-${htbin}_25ns/151101_*/0000/*root > \
   ${submitbase}/${version}/lists/hdfslist_${submitname}.txt

  # format as xrootd
  cp ${submitbase}/${version}/lists/hdfslist_${submitname}.txt \
     ${submitbase}/${version}/lists/xrdlist_${submitname}.txt 
  xrdlist="${submitbase}/${version}/lists/xrdlist_${submitname}.txt"
  sed -i 's@/hdfs/@root://cmsxrootd.hep.wisc.edu//@g' $xrdlist #
  
  printf "  done making list of files\n"
  printf " making submit template\n"

  # sample specific parameters..
  treename="ggNtuplizer/EventTree"
  isMC="kTRUE"
  xc="$(grep -P ${submitname} ${submitbase}/sample_info/xclist.txt | sed -n -e "s@${submitname}  XC : @@p")"
   #printf "\nxc is ${xc}\n"

  # make correct executable xx_callpostAnalyzer_QCD
  cp template_callpostAnalyzer_QCD.cc    "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@SAMPLENAME@${submitname}@g"  "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@TREENAME@${treename}@g"      "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@ISMC@${isMC}@g"              "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@XC@${xc}@g"                  "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"

  printf "  done making submit template\n"

  if [ ${dosubmit} = true ]
  then
  printf " submitting\n"

  farmoutAnalysisJobs \
   --infer-cmssw-path \
   --fwklite \
   --input-file-list=${xrdlist} \
   --input-files-per-job=30 \
   --use-hdfs \
   --extra-inputs=${submitbase}/postAnalyzer_QCD.C,${submitbase}/postAnalyzer_QCD.h \
   ${version} \
   "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"

  fi # ${dosubmit} = true

  done # for htbin in ..
 done # for mc_samplename in ..
fi # ${domc} = true






######## Data ################################

if [ ${dodata} = true ]
then
 for data_samplename in \
  "SinglePhoton_2015D"

 do
  printf "\n${data_samplename}\n"
  
  # hdfs/store/user/gomber/SinglePhoton_2015D_1p2fb1_condor/run_data_2015D_74X-*root
  printf " making list of files\n"

  find /hdfs/store/user/gomber/${data_samplename}_1p2fb1_condor/run_data_2015D_74X-*root > \
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

  # make correct executable xx_callpostAnalyzer_QCD
  cp template_callpostAnalyzer_QCD.cc         "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@SAMPLENAME@${data_samplename}@g"  "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@TREENAME@${treename}@g"           "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@ISMC@${isMC}@g"                   "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@XC@${xc}@g"                       "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"

  printf "  done making submit template\n"

  if [ ${dosubmit} = true ]
  then
  printf " submitting\n"

  farmoutAnalysisJobs \
   --infer-cmssw-path \
   --fwklite \
   --input-file-list=${xrdlist} \
   --input-files-per-job=30 \
   --use-hdfs \
   --extra-inputs=${submitbase}/postAnalyzer_QCD.C,${submitbase}/postAnalyzer_QCD.h \
   ${version} \
   "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"

  fi # ${dosubmit} = true

 done # for data_samplename in ..
fi # ${dodata}
