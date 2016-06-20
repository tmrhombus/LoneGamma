#!/bin/bash

#voms-proxy-init --voms cms --valid 100:00

domc=true
dodata=true
dosubmit=true

START=$(date +%s);
printf "Started at `date`"

mkdir -p "${submitbase}/gitignore/${version}/lists"
mkdir -p "${submitbase}/gitignore/${version}/submit"

lumi=2200. # /pb

if [ ${domc} = true ]
then

 initevents="${CMSSW_BASE}/src/LoneGamma/lists/initialEvents.txt"

 for mc_samplename in \
  "GJetsHT" \
  "QCDPt"
 do

  if [ ${mc_samplename} = "GJetsHT" ]
   then 
    bins=( "40to100" "100to200" "200to400" "400to600" "600toInf" )
   fi
  if [ ${mc_samplename} = "QCDPt" ]
   then
    bins=( "80to120" "120to170" "170to300" "300toInf" )
    #bins=( "15to20" "20to30" "30to50" "50to80" "80to120" "120to170" "170to300" "300toInf" )
  fi

  for bin in ${bins[*]}
  do  

  submitname="${mc_samplename}${bin}"
  printf "${submitname}\n"

   # format as xrootd
   cp ${CMSSW_BASE}/src/LoneGamma/lists/hdfslist_${submitname}.txt \
      ${submitbase}/gitignore/${version}/lists/xrdlist_${submitname}.txt 
   xrdlist="${submitbase}/gitignore/${version}/lists/xrdlist_${submitname}.txt"
   sed -i 's@/hdfs/@root://cmsxrootd.hep.wisc.edu//@g' $xrdlist #

  # sample specific parameters..
  treename="ggNtuplizer/EventTree"
  isMC="kTRUE"
  isELE="kFALSE"
  nrE="$(grep -P ${submitname} ${initevents} | sed -n -e "s@${submitname} Events: @@p")"
  xc="$(grep -P ${submitname}  ${initevents} | sed -n -e "s@${submitname} XC: @@p")"
   #printf "\nxc is ${xc}\n"

  # make correct executable xx_callpostAnalyzer_QCD
  cp template_callpostAnalyzer_QCD.cc    "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@SAMPLENAME@${submitname}@g"  "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@TREENAME@${treename}@g"      "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@ISMC@${isMC}@g"              "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@ISELE@${isELE}@g"            "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@CROSSSEC@${xc}@g"            "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@NREVENTS@${nrE}@g"           "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"
  sed -i "s@LUMI@${lumi}@g"              "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"

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
   "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_QCD.cc"

  fi # ${dosubmit} = true

  done # for htbin in ..
 done # for mc_samplename in ..
fi # ${domc} = true


######## Data ################################

if [ ${dodata} = true ]
then
 for data_samplename in \
  "SinglePhoton" \
  "DoubleElectron"

 do
  printf "\n${data_samplename}\n"
  printf " making list of files\n"

  if [ "$data_samplename" = "SinglePhoton" ]
  then 
   isELE="kFALSE"
  fi

  if [ "$data_samplename" = "DoubleElectron" ]
  then
   isELE="kTRUE"
  fi

  # format as xrootd
   cp ${CMSSW_BASE}/src/LoneGamma/lists/hdfslist_${data_samplename}.txt \
     ${submitbase}/gitignore/${version}/lists/xrdlist_${data_samplename}.txt 
  xrdlist="${submitbase}/gitignore/${version}/lists/xrdlist_${data_samplename}.txt"
  sed -i 's@/hdfs/@root://cmsxrootd.hep.wisc.edu//@g' $xrdlist #

  printf "  done making list of files\n"
  printf " making submit template\n"

  # sample specific parameters..
  treename="ggNtuplizer/EventTree"
  isMC="kFALSE"
  xc="1."
  nrE=${lumi}

  # make correct executable xx_callpostAnalyzer_QCD
  cp template_callpostAnalyzer_QCD.cc         "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@SAMPLENAME@${data_samplename}@g"  "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@TREENAME@${treename}@g"           "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@ISMC@${isMC}@g"                   "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@ISELE@${isELE}@g"                 "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@CROSSSEC@${xc}@g"                 "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@NREVENTS@${nrE}@g"                "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"
  sed -i "s@LUMI@${lumi}@g"                   "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"

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
    "${submitbase}/gitignore/${version}/submit/${data_samplename}_callpostAnalyzer_QCD.cc"

  fi # ${dosubmit} = true

 done # for data_samplename in ..
fi # ${dodata}
