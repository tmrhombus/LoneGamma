#!/bin/bash

#voms-proxy-init --voms cms --valid 100:00

domc=false
dodata=true
dosubmit=true

START=$(date +%s);
printf "Started at `date`"

mkdir -p "${submitbase}/${version}/lists"
mkdir -p "${submitbase}/${version}/submit"

lumi=2320. # /pb

if [ ${domc} = true ]
then

 initevents="${submitbase}/${version}/lists/initialEvents.txt"
 echo " " >> ${initevents}
 touch ${initevents}

 for mc_samplename in \
  "GJets_HT" \
  "QCD_Pt"
 do

  if [ ${mc_samplename} = "GJets_HT" ]
   then 
    bins=( "40To100" "100To200" "200To400" "400To600" "600ToInf" )
   fi
  if [ ${mc_samplename} = "QCD_Pt" ]
   then
    bins=( "15to20" "20to30" "30to50" "50to80" "80to120" "120to170" "170to300" "300toInf" )
  fi

  for bin in ${bins[*]}
  do  
   echo $bin

   trynr="5"
   if [ ${bin} = "100To200" ] 
   then
    trynr="4"
   fi

  submitname="${mc_samplename}${bin}"

 # count the total number of events, put in a file if it's not there
 if ! grep -F "${submitname} " ${initevents}
 then 
  echo "Need to count events.. making list"

   # /hdfs/store/user/jjbuch/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-15to20_EMEnriched_try5/151212_075448/0000/
   # /hdfs/store/user/jjbuch/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-20to30_EMEnriched_try5/151212_075603/0000/
   # /hdfs/store/user/jjbuch/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-30to50_EMEnriched_try5/151212_075643/0000/
   # /hdfs/store/user/jjbuch/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-50to80_EMEnriched_try5/151212_075742/0000/
   # /hdfs/store/user/jjbuch/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-80to120_EMEnriched_try5/151212_075831/0000/
   # /hdfs/store/user/jjbuch/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-120to170_EMEnriched_try5/151212_080000/0000/
   # /hdfs/store/user/jjbuch/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-170to300_EMEnriched_try5/151212_080124/0000/
   # /hdfs/store/user/jjbuch/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-300toInf_EMEnriched_try5/151212_080212/0000/

   # /hdfs/store/user/jjbuch/${mc_samplename}-${bin}_*TuneCUETP8M1_13TeV*pythia8/crab_ggNtuplizer_spring15_${mc_samplename}-${${bin}_try${trynr}/151212_*/0000/
   # /hdfs/store/user/jjbuch/${mc_samplename}-${bin}_*TuneCUETP8M1_13TeV*pythia8/crab_ggNtuplizer_spring15_${mc_samplename}-${${bin}_try${trynr}/
   # /hdfs/store/user/jjbuch/${mc_samplename}-${bin}_*TuneCUETP8M1_13TeV*pythia8/crab_ggNtuplizer_spring15_${mc_samplename}-${${bin}_try${trynr}/
   # /hdfs/store/user/jjbuch/${mc_samplename}-${bin}_*TuneCUETP8M1_13TeV*pythia8/crab_ggNtuplizer_spring15_${mc_samplename}-${${bin}_try${trynr}/
   # /hdfs/store/user/jjbuch/${mc_samplename}-${bin}_*TuneCUETP8M1_13TeV*pythia8/crab_ggNtuplizer_spring15_${mc_samplename}-${${bin}_try${trynr}/

   find /hdfs/store/user/jjbuch/${mc_samplename}-${bin}_*TuneCUETP8M1_13TeV*pythia8/crab_ggNtuplizer_spring15_${mc_samplename}-*_try${trynr}/151212_*/0000/*root > \
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

  # make correct executable xx_callpostAnalyzer_QCD
  cp template_callpostAnalyzer_QCD.cc    "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@SAMPLENAME@${submitname}@g"  "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@TREENAME@${treename}@g"      "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@ISMC@${isMC}@g"              "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@CROSSSEC@${xc}@g"            "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@NREVENTS@${nrE}@g"           "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@LUMI@${lumi}@g"              "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"

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
   "${submitbase}/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"

  fi # ${dosubmit} = true

  done # for htbin in ..
 done # for mc_samplename in ..
fi # ${domc} = true


######## Data ################################

  #"SinglePhoton" 
if [ ${dodata} = true ]
then
 for data_samplename in \
  "DoubleElectron"

 do
  printf "\n${data_samplename}\n"
  printf " making list of files\n"

  if [ "$data_samplename" = "SinglePhoton" ]
  then 
   find /hdfs/store/user/gomber/SinglePhoton_Crab_2015D_v3_226fb/SinglePhoton/crab_job_single_photon_13TeV_v3_226fb/160109_082344/0000/*root > \
    ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt
   find /hdfs/store/user/gomber/SinglePhoton_Crab_2015D_v4_226fb/SinglePhoton/crab_job_single_photon_13TeV_v4_226fb/160109_082133/000*/*root >> \
    ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt
  fi

  if [ "$data_samplename" = "DoubleElectron" ]
  then
   find /hdfs/store/user/gomber/DoubleElectron_Crab_2015D_v3_226fb/*/crab_job_*_13TeV_v3_226fb/160216_*/0000/*root > \
    ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt
   find /hdfs/store/user/gomber/DoubleElectron_Crab_2015D_v4_226fb/*/crab_job_*_13TeV_v4_226fb/160216_*/000*/*root >> \
    ${submitbase}/${version}/lists/hdfslist_${data_samplename}.txt
  fi

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
  nrE="2260."

  # make correct executable xx_callpostAnalyzer_QCD
  cp template_callpostAnalyzer_QCD.cc         "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@SAMPLENAME@${data_samplename}@g"  "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@TREENAME@${treename}@g"           "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@ISMC@${isMC}@g"                   "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@CROSSSEC@${xc}@g"                 "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@NREVENTS@${nrE}@g"                "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@LUMI@${lumi}@g"                   "${submitbase}/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"

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
