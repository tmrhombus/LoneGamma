#!/bin/bash

#voms-proxy-init --voms cms --valid 100:00

domc=true
dodata=true
dosubmit=true

START=$(date +%s);
printf "Started at `date`\n\n"

mkdir -p "${submitbase}/gitignore/${version}/lists"
mkdir -p "${submitbase}/gitignore/${version}/submit"

lumi=12900. # /pb

if [ ${domc} = true ]
then

 initevents="${CMSSW_BASE}/src/LoneGamma/lists/initialEvents.txt"

  #"QCDPt"     \
#  "Wen"       \
 for mc_samplename in \
  "GJetsHT"   \
  "ZllJetsHT" \
  "ZnnGJets"  \
  "ZllGJets"  \
  "WlnGJets"  \
  "Wmn"       \
  "Wtn"       \
  "WWG"       \
  "TTGJets"   \
  "GGJets"    \
  "TGJets"    \
  "WZ"        \
  "ZZ"

 do

  bins=( "" )
  if [ ${mc_samplename} = "GJetsHT" ]
   then 
    bins=( "40to100" "100to200" "200to400" "400to600" "600toInf" )
   fi
  if [ ${mc_samplename} = "QCDPt" ]
   then
    bins=( "80to120" "120to170" "170to300" "300toInf" )
    #bins=( "15to20" "20to30" "30to50" "50to80" "80to120" "120to170" "170to300" "300toInf" )
  fi
  if [ ${mc_samplename} = "ZllJetsHT" ]
   then
    bins=( "100to200" "200to400" "400to600" "600toInf" )
  fi
  if [ ${mc_samplename} = "ZnnJetsHT" ]
   then
    bins=( "100to200" "200to400" "400to600" "600toInf" )
  fi

  for bin in "${bins[@]}"
  do  

  submitname="${mc_samplename}${bin}"
  printf "\nPreparing  ${submitname}\n"

   # format as xrootd
   cp ${CMSSW_BASE}/src/LoneGamma/lists/hdfslist_${submitname}.txt \
      ${submitbase}/gitignore/${version}/lists/xrdlist_${submitname}.txt 
   xrdlist="${submitbase}/gitignore/${version}/lists/xrdlist_${submitname}.txt"
   sed -i 's@/hdfs/@root://cmsxrootd.hep.wisc.edu//@g' $xrdlist #

  # sample specific parameters..
  treename="ggNtuplizer/EventTree"
  isMC="kTRUE"
  isZnnG="kFALSE"
  ewkWG="kFALSE"
  isEle="kFALSE"
  isJet="kFALSE"
  if [ ${submitname} = "ZnnGJets" ]
   then 
   isZnnG="kTRUE" 
  fi
  if [ ${submitname} = "WlnGJets" ]
   then 
   ewkWG="kTRUE"
  fi
  nrE="$(grep -P ${submitname} ${initevents} | sed -n -e "s@ ${submitname} Events: @@p")"
  xc="$(grep -P ${submitname}  ${initevents} | sed -n -e "s@ ${submitname} XC: @@p")"
   #printf "\nxc is ${xc}\n"

  # make correct executable xx_callpostAnalyzer_WlnG
  cp template_callpostAnalyzer_WlnG.cc    "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@SAMPLENAME@${submitname}@g"     "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@TREENAME@${treename}@g"         "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@ISMC@${isMC}@g"                 "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@ISZNNG@${isZnnG}@g"             "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@ISELE@${isEle}@g"               "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@ISJET@${isJet}@g"               "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@EWKZG@${ewkWG}@g"               "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@CROSSSEC@${xc}@g"               "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@NREVENTS@${nrE}@g"              "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
  sed -i "s@LUMI@${lumi}@g"                 "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"

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
   --extra-inputs=${submitbase}/postAnalyzer_WlnG.C,${submitbase}/postAnalyzer_WlnG.h,${submitbase}/ewk_corr.root \
   ${version} \
   "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"

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

  for datacut in \
   "Data"  \
   "Ele"   \
   "Jet"

  do  

   submitname="${data_samplename}${datacut}"
   printf "\nPreparing  ${submitname}\n"

   # format as xrootd
    cp ${CMSSW_BASE}/src/LoneGamma/lists/hdfslist_${data_samplename}.txt \
      ${submitbase}/gitignore/${version}/lists/xrdlist_${submitname}.txt 
   xrdlist="${submitbase}/gitignore/${version}/lists/xrdlist_${submitname}.txt"
   sed -i 's@/hdfs/@root://cmsxrootd.hep.wisc.edu//@g' $xrdlist #

   printf "  done making list of files\n"
   printf " making submit template\n"

   # sample specific parameters..
   treename="ggNtuplizer/EventTree"
   isMC="kFALSE"
   isZnnG="kFALSE"
   ewkWG="kFALSE"
   xc="1."
   nrE=${lumi}
   isEle="kFALSE"
   isHalo="kFALSE"
   isSpike="kFALSE"
   isJet="kFALSE"
   if [ ${datacut} = "Ele" ]
   then 
    isEle="kTRUE"
   fi  
   if [ ${datacut} = "Halo" ]
   then 
    isHalo="kTRUE"
   fi  
   if [ ${datacut} = "Spike" ]
   then 
    isSpike="kTRUE"
   fi  
   if [ ${datacut} = "Jet" ]
   then 
    isJet="kTRUE"
   fi 

   # make correct executable xx_callpostAnalyzer_WlnG
   cp template_callpostAnalyzer_WlnG.cc      "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@SAMPLENAME@${submitname}@g"       "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@TREENAME@${treename}@g"           "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@ISMC@${isMC}@g"                   "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@ISZNNG@${isZnnG}@g"               "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@EWKZG@${ewkWG}@g"                 "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@ISELE@${isEle}@g"                 "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@ISJET@${isJet}@g"                 "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@CROSSSEC@${xc}@g"                 "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@NREVENTS@${nrE}@g"                "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"
   sed -i "s@LUMI@${lumi}@g"                   "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"

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
     --extra-inputs=${submitbase}/postAnalyzer_WlnG.C,${submitbase}/postAnalyzer_WlnG.h,${submitbase}/ewk_corr.root \
     ${version} \
     "${submitbase}/gitignore/${version}/submit/${submitname}_callpostAnalyzer_WlnG.cc"

   fi # ${dosubmit} = true

  done # for datacut in Data Ele Jet
 done # for submitname in ..
fi # ${dodata}
