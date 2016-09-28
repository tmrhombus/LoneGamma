#!/bin/bash

outdir="${CMSSW_BASE}/src/LoneGamma/lists"
mkdir -p "${outdir}"

doSinglePhoton=false
doDoubleElectron=false
doGJetsMC=false
 # GJetsHT40to100
 # GJetsHT100to200
 # GJetsHT200to400
 # GJetsHT400to600
 # GJetsHT600toInf
doQCDMC=false
 # QCDPt80to120
 # QCDPt120to170
 # QCDPt170to300
 # QCDPt300toInf
doZllJMC=false
 # ZllJetsHT100to200
 # ZllJetsHT200to400
 # ZllJetsHT400to600
 # ZllJetsHT600toInf
doZnnJMC=false
 # ZnnJetsHT100to200
 # ZnnJetsHT200to400
 # ZnnJetsHT400to600
 # ZnnJetsHT600toInf
doZllGJMC=false
 # ZllGJets
doZnnGJMC=true
 # ZnnGJets
doWlnGJMC=false
 # WlnGJets
doWenMC=false
 # Wen
doWmnMC=false
 # Wmn
doWtnMC=false
 # Wtn
doTTGJMC=false
 # TTGJets
doWWGMC=false
 # WWG
doGGJMC=false
 # GGJets
doTGJMC=false
 # TGJets
doWZMC=false
 # WZ
doZZMC=false
 # ZZ


#SinglePhoton
#------------
if [ ${doSinglePhoton} = true ]
then
 printf "Making SinglePhoton\n" 
 ls  /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton*12p9*newbh_1/*/*/*root \
 > "${outdir}/hdfslist_SinglePhoton.txt"
fi


## #DoubleElectron
## #------------
## if [ ${doDoubleElectron} = true ]
## then
##  printf "Making DoubleElectron\n" 
## 
##    find /hdfs/store/user/gomber/DoubleElectron_Crab_2015D_v3_226fb/*/crab_job_*_13TeV_v3_226fb/160216_*/0000/*root \
##      >  "${outdir}/hdfslist_DoubleElectron.txt"
##    find /hdfs/store/user/gomber/DoubleElectron_Crab_2015D_v4_226fb/*/crab_job_*_13TeV_v4_226fb/160216_*/000*/*root \
##      >> "${outdir}/hdfslist_DoubleElectron.txt"
## fi

#GJets MC
#------------
if [ ${doGJetsMC} = true ]
then
 printf "Making GJets " 
 printf "1  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-40To100/160621_133030/0000/*root \
  > "${outdir}/hdfslist_GJetsHT40to100.txt"

 printf "2  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-100To200/160621_133042/0000/*root \
  > "${outdir}/hdfslist_GJetsHT100to200.txt"

 printf "3  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-200To400/160621_133055/0000/*root \
  > "${outdir}/hdfslist_GJetsHT200to400.txt"

 printf "4  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-400To600/160621_133110/0000/*root \
  > "${outdir}/hdfslist_GJetsHT400to600.txt"

 printf "5  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-600ToInf/160621_133123/0000/*root \
  > "${outdir}/hdfslist_GJetsHT600toInf.txt"

 printf "\n" 
fi

#QCD MC
#----------
if [ ${doQCDMC} = true ]
then
 printf "Making QCD " 
# printf "1  " 
# ls /hdfs/store/user/jjbuch/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-15to20_EMEnriched_try5/151212_075448/0000/*root \
#  > "${outdir}/hdfslist_QCDPt15to20.txt"
#
# printf "2  " 
# ls /hdfs/store/user/jjbuch/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-20to30_EMEnriched_try5/151212_075603/0000/*root \
#  > "${outdir}/hdfslist_QCDPt20to30.txt"
#
# printf "3  " 
# ls /hdfs/store/user/jjbuch/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-30to50_EMEnriched_try5/151212_075643/0000/*root \
#  > "${outdir}/hdfslist_QCDPt30to50.txt"
#
# printf "4  " 
# ls /hdfs/store/user/jjbuch/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-50to80_EMEnriched_try5/151212_075742/0000/*root \
#  > "${outdir}/hdfslist_QCDPt50to80.txt"

 printf "5  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-80to120/160621_133148/0000/*root \
  > "${outdir}/hdfslist_QCDPt80to120.txt"

 printf "6  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-120to170/160621_133201/0000/*root \
  > "${outdir}/hdfslist_QCDPt120to170.txt"

 printf "7  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-170to300/160621_133213/0000/*root \
  > "${outdir}/hdfslist_QCDPt170to300.txt"

 printf "8  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-300toInf/160621_133225/0000/*root \
> "${outdir}/hdfslist_QCDPt300toInf.txt"

 printf "\n" 
fi

# Z(ll)Jets MC
#----------
if [ ${doZllJMC} = true ]
then
 printf "Making Z(ll)Jets MC " 
 printf "1  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-100to200/160621_133238/0000/*root \
  > "${outdir}/hdfslist_ZllJetsHT100to200.txt"

 printf "2  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-200to400/160621_133251/0000/*root \
  > "${outdir}/hdfslist_ZllJetsHT200to400.txt"

 printf "3  "
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-400to600/160621_133304/0000/*root \
  > "${outdir}/hdfslist_ZllJetsHT400to600.txt"

 printf "4  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-600toInf/160621_133316/0000/*root \
 > "${outdir}/hdfslist_ZllJetsHT600toInf.txt"

 printf "\n" 
fi

###Z(nunu)Jets MC
###------------
##if [ ${doZnnJMC} = true ]
##then
## printf "Making Z(nunu)Jets MC " 
## printf "1  " 
##  ls /hdfs/store/user/jjbuch/ZJetsToNuNu_HT-100To200_13TeV-madgraph/crab_ggNtuplizer_spring15_ZJetsToNuNu_HT-100to200_try5/151212_074908/0000/ggtree_mc_*root \
##  > "${outdir}/hdfslist_ZnnJetsHT100to200.txt"
##
## printf "2  " 
##  ls /hdfs/store/user/jjbuch/ZJetsToNuNu_HT-200To400_13TeV-madgraph/crab_ggNtuplizer_spring15_ZJetsToNuNu_HT-200to400_try5/151212_074940/0000/ggtree_mc_*root \
##  > "${outdir}/hdfslist_ZnnJetsHT200to400.txt"
##
## printf "3  "
##  ls /hdfs/store/user/jjbuch/ZJetsToNuNu_HT-400To600_13TeV-madgraph/crab_ggNtuplizer_spring15_ZJetsToNuNu_HT-400to600_try5/151212_075020/0000/ggtree_mc_*root \
##  > "${outdir}/hdfslist_ZnnJetsHT400to600.txt"
##
## printf "4  " 
##  ls /hdfs/store/user/jjbuch/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph/crab_ggNtuplizer_spring15_ZJetsToNuNu_HT-600toInf_try5/151212_075101/0000/ggtree_mc_*root \
##  > "${outdir}/hdfslist_ZnnJetsHT600toInf.txt"
##
## printf "\n" 
##
##fi


#Z(ll)Gamma Jets
#------------
if [ ${doZllGJMC} = true ]
then
 printf "Making Z(ll)Gamma Jets MC\n" 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/*root \
  > "${outdir}/hdfslist_ZllGJets.txt"
fi

#Z(nn)Gamma Jets
#------------
if [ ${doZnnGJMC} = true ]
then
 printf "Making Z(nn)Gamma Jets MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets/160621_132818/0000/*root \
  > "${outdir}/hdfslist_ZnnGJets.txt"
fi

#W(ln)Gamma Jets
#------------
if [ ${doWlnGJMC} = true ]
then
 printf "Making Z(ln)Gamma Jets MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/*root \
  > "${outdir}/hdfslist_WlnGJets.txt"
fi

###W(en)
###------------
##if [ ${doWenMC} = true ]
##then
## printf "Making W(en) MC\n" 
##
##  ls /hdfs/store/user/jjbuch/WToENu_M-100_MiniAOD/WToENu_M-100_TuneCUETP8M1_13TeV-pythia8/crab_WToENu_M-100_MiniAOD/160413_214913/0000/ggtree_mc_DiPhoton_MINIAOD_TEST_*root \
##  > "${outdir}/hdfslist_Wen.txt"
##fi

#W(mn)
#------------
if [ ${doWmnMC} = true ]
then
 printf "Making W(mn) MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WToMuNu_M-100_TuneCUETP8M1_13TeV-pythia8/crab_WToMuNu/160621_133001/0000/*root \
  > "${outdir}/hdfslist_Wmn.txt"
fi

#W(tn)
#------------
if [ ${doWtnMC} = true ]
then
 printf "Making W(tn) MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WToTauNu_M-100_TuneCUETP8M1_13TeV-pythia8-tauola/crab_WToTauNu/160621_133014/0000/*root \
  > "${outdir}/hdfslist_Wtn.txt"
fi

#TTGJ
#------------
if [ ${doTTGJMC} = true ]
then
 printf "Making TTGJets MC\n" 

  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/160621_133135/0000/*root \
  > "${outdir}/hdfslist_TTGJets.txt"
fi

#WWG
#------------
if [ ${doWWGMC} = true ]
then
 printf "Making WWG MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_WWG/160630_154410/0000/*.root \
  > "${outdir}/hdfslist_WWG.txt"
fi

#GGJ
#------------
if [ ${doGGJMC} = true ]
then
 printf "Making GGJets MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/crab_DiPhotonJets_MGG-80toInf_amcatnlo/160825_000327/0000/*.root \
  > "${outdir}/hdfslist_GGJets.txt"
fi

#TGJ
#------------
if [ ${doTGJMC} = true ]
then
 printf "Making TGJets MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_TGJets/160823_154537/0000/*.root \
  > "${outdir}/hdfslist_TGJets.txt"
fi

#WZ
#------------
if [ ${doWZMC} = true ]
then
 printf "Making WZ MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/160630_154423/0000/*.root \
  > "${outdir}/hdfslist_WZ.txt"
fi

#ZZ
#------------
if [ ${doZZMC} = true ]
then
 printf "Making ZZ MC\n" 
  ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/160630_154436/0000/*.root \
  > "${outdir}/hdfslist_ZZ.txt"
fi
