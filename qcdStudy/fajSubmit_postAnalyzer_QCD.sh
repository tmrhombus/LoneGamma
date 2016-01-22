#!/bin/bash

#voms-proxy-init --voms cms --valid 100:00

domc=false
dodata=true
dosubmit=true

START=$(date +%s);
printf "Started at `date`"

mkdir -p "${submitbase}/${version}/lists"
mkdir -p "${submitbase}/${version}/submit"

# lumi = 

if [ ${domc} = true ]
then
  #"GJets_HT"  
 for mc_samplename in \
  "QCD_Pt"
 do

   #"40To100" \
   #"100To200" \
   #"200To400" \
   #"400To600" \
   #"600ToInf" 
  for htbin in \
   "15to20" \
   "20to30" \
   "30to50" \
   "50to80" \
   "80to120" \
   "120to170" \
   "170to300" \
   "300toInf" 
  do

  submitname="${mc_samplename}${htbin}"

  printf "${submitname}\n"

  # /hdfs/store/user/jjbuch/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-15to20_EMEnriched_try5/151212_075448/0000/
  # /hdfs/store/user/jjbuch/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-20to30_EMEnriched_try5/151212_075603/0000/
  # /hdfs/store/user/jjbuch/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-30to50_EMEnriched_try5/151212_075643/0000/
  # /hdfs/store/user/jjbuch/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-50to80_EMEnriched_try5/151212_075742/0000/
  # /hdfs/store/user/jjbuch/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-80to120_EMEnriched_try5/151212_075831/0000/
  # /hdfs/store/user/jjbuch/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-120to170_EMEnriched_try5/151212_080000/0000/
  # /hdfs/store/user/jjbuch/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-170to300_EMEnriched_try5/151212_080124/0000/
  # /hdfs/store/user/jjbuch/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-300toInf_EMEnriched_try5/151212_080212/0000/

  
  # hdfs/store/user/jjbuch/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-40To100_25ns/151101_093834/0000/*root
  # hdfs/store/user/jjbuch/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-100To200_25ns/151101_092459/0000/*root
  # hdfs/store/user/jjbuch/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-200To400_25ns/151101_092710/0000/*root
  # hdfs/store/user/jjbuch/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-400To600_25ns/151101_093518/0000/*root
  # hdfs/store/user/jjbuch/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-600ToInf_25ns/151101_094056/0000/*root
  printf " making list of files\n"
  #find /hdfs/store/user/jjbuch/${mc_samplename}-${htbin}_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_gjet_HT-${htbin}_25ns/151101_*/0000/*root > \
  find /hdfs/store/user/jjbuch/${mc_samplename}-${htbin}_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_${mc_samplename}-${htbin}_EMEnriched_try5/151212_*/0000/*root > \
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
  "SinglePhoton"

 do
  printf "\n${data_samplename}\n"
  
  printf " making list of files\n"

#  #find /hdfs/store/user/jjbuch/SinglePhoton/crab_ggNtuplizer_spring15_SinglePhoton_Run2015D_PromptReco_v3_try5/151212_080518/0000/*root > \
#  find /hdfs/store/user/gomber/SinglePhoton_Crab_2015D_v3_226fb/SinglePhoton/crab_job_single_photon_13TeV_v3_226fb/160109_082344/0000/*root > \
#   ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt
#  #find /hdfs/store/user/jjbuch/SinglePhoton/crab_ggNtuplizer_spring15_SinglePhoton_Run2015D_PromptReco_v4_try5/151212_080439/0000/*root >> \
#  find /hdfs/store/user/gomber/SinglePhoton_Crab_2015D_v4_226fb/SinglePhoton/crab_job_single_photon_13TeV_v4_226fb/160109_082133//000*/*root >> \
#   ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt
#
#  # format as xrootd
#  cp ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt \
#     ${submitbase}/${version}/lists/xrdlist_${data_samplename}.txt 
#  xrdlist="${submitbase}/${version}/lists/xrdlist_${data_samplename}.txt"
#  sed -i 's@/hdfs/@root://cmsxrootd.hep.wisc.edu//@g' $xrdlist #

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
    --input-files-per-job=100 \
    --use-hdfs \
    --extra-inputs=${submitbase}/postAnalyzer_QCD.C,${submitbase}/postAnalyzer_QCD.h \
    ${version} \
    "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"

  fi # ${dosubmit} = true

 done # for data_samplename in ..
fi # ${dodata}
