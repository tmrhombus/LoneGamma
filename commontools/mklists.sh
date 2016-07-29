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
doZnnGJMC=false
 # ZnnGJets
doWlnGJMC=true
 # WlnGJets
doWenMC=true
 # Wen
doWmnMC=true
 # Wmn
doWtnMC=true
 # Wtn
doTTGJMC=true
 # TTGJets


#SinglePhoton
#------------
if [ ${doSinglePhoton} = true ]
then
 printf "Making SinglePhoton\n" 

  ls /hdfs/store/user/gomber/SinglePhoton_Crab_2015D_v3_226fb/SinglePhoton/crab_job_single_photon_13TeV_v3_226fb/160109_082344/0000/*root \
    >  "${outdir}/hdfslist_SinglePhoton.txt"
  ls /hdfs/store/user/gomber/SinglePhoton_Crab_2015D_v4_226fb/SinglePhoton/crab_job_single_photon_13TeV_v4_226fb/160109_082133/000*/*root \
    >> "${outdir}/hdfslist_SinglePhoton.txt"
fi


#DoubleElectron
#------------
if [ ${doDoubleElectron} = true ]
then
 printf "Making DoubleElectron\n" 

   find /hdfs/store/user/gomber/DoubleElectron_Crab_2015D_v3_226fb/*/crab_job_*_13TeV_v3_226fb/160216_*/0000/*root \
     >  "${outdir}/hdfslist_DoubleElectron.txt"
   find /hdfs/store/user/gomber/DoubleElectron_Crab_2015D_v4_226fb/*/crab_job_*_13TeV_v4_226fb/160216_*/000*/*root \
     >> "${outdir}/hdfslist_DoubleElectron.txt"
fi

#GJets MC
#------------
if [ ${doGJetsMC} = true ]
then
 printf "Making GJets " 
 printf "1  "  
 ls /hdfs/store/user/jjbuch/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_GJets_HT-40to100_try5/151212_074559/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_GJetsHT40to100.txt"

 printf "2  " 
 ls /hdfs/store/user/jjbuch/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_GJets_HT-100to200_try4/151212_042213/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_GJetsHT100to200.txt"

 printf "3  " 
 ls /hdfs/store/user/jjbuch/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_GJets_HT-200to400_try5/151212_074638/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_GJetsHT200to400.txt"

 printf "4  " 
 ls /hdfs/store/user/jjbuch/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_GJets_HT-400to600_try5/151212_074715/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_GJetsHT400to600.txt"

 printf "5  " 
 ls /hdfs/store/user/jjbuch/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ggNtuplizer_spring15_GJets_HT-600toInf_try5/151212_074751/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_GJetsHT600toInf.txt"

 printf "\n" 
fi

#QCD MC
#----------
if [ ${doQCDMC} = true ]
then
 printf "Making QCD " 
 printf "1  " 
 ls /hdfs/store/user/jjbuch/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-15to20_EMEnriched_try5/151212_075448/0000/*root \
  > "${outdir}/hdfslist_QCDPt15to20.txt"

 printf "2  " 
 ls /hdfs/store/user/jjbuch/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-20to30_EMEnriched_try5/151212_075603/0000/*root \
  > "${outdir}/hdfslist_QCDPt20to30.txt"

 printf "3  " 
 ls /hdfs/store/user/jjbuch/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-30to50_EMEnriched_try5/151212_075643/0000/*root \
  > "${outdir}/hdfslist_QCDPt30to50.txt"

 printf "4  " 
 ls /hdfs/store/user/jjbuch/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-50to80_EMEnriched_try5/151212_075742/0000/*root \
  > "${outdir}/hdfslist_QCDPt50to80.txt"

 printf "5  " 
 ls /hdfs/store/user/jjbuch/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-80to120_EMEnriched_try5/151212_075831/0000/*root \
  > "${outdir}/hdfslist_QCDPt80to120.txt"

 printf "6  " 
 ls /hdfs/store/user/jjbuch/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-120to170_EMEnriched_try5/151212_080000/0000/*root \
  > "${outdir}/hdfslist_QCDPt120to170.txt"

 printf "7  " 
 ls /hdfs/store/user/jjbuch/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-170to300_EMEnriched_try5/151212_080124/0000/*root \
  > "${outdir}/hdfslist_QCDPt170to300.txt"

 printf "8  " 
 ls /hdfs/store/user/jjbuch/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_ggNtuplizer_spring15_QCD_Pt-300toInf_EMEnriched_try5/151212_080212/0000/*root \
  > "${outdir}/hdfslist_QCDPt300toInf.txt"

 printf "\n" 
fi

# Z(ll)Jets MC
#----------
if [ ${doZllJMC} = true ]
then
 printf "Making Z(ll)Jets MC " 
 printf "1  " 
  ls /hdfs/store/user/jjbuch/DYJetsToLL_74X/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-100to200/160530_182830/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZllJetsHT100to200.txt"

 printf "2  " 
  ls /hdfs/store/user/jjbuch/DYJetsToLL_74X/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-200to400/160530_182903/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZllJetsHT200to400.txt"

 printf "3  "
  ls /hdfs/store/user/jjbuch/DYJetsToLL_74X/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-400to600/160530_182941/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZllJetsHT400to600.txt"

 printf "4  " 
  ls /hdfs/store/user/jjbuch/DYJetsToLL_74X/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-600toInf/160530_182718/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZllJetsHT600toInf.txt"

 printf "\n" 
fi

#Z(nunu)Jets MC
#------------
if [ ${doZnnJMC} = true ]
then
 printf "Making Z(nunu)Jets MC " 
 printf "1  " 
  ls /hdfs/store/user/jjbuch/ZJetsToNuNu_HT-100To200_13TeV-madgraph/crab_ggNtuplizer_spring15_ZJetsToNuNu_HT-100to200_try5/151212_074908/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZnnJetsHT100to200.txt"

 printf "2  " 
  ls /hdfs/store/user/jjbuch/ZJetsToNuNu_HT-200To400_13TeV-madgraph/crab_ggNtuplizer_spring15_ZJetsToNuNu_HT-200to400_try5/151212_074940/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZnnJetsHT200to400.txt"

 printf "3  "
  ls /hdfs/store/user/jjbuch/ZJetsToNuNu_HT-400To600_13TeV-madgraph/crab_ggNtuplizer_spring15_ZJetsToNuNu_HT-400to600_try5/151212_075020/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZnnJetsHT400to600.txt"

 printf "4  " 
  ls /hdfs/store/user/jjbuch/ZJetsToNuNu_HT-600ToInf_13TeV-madgraph/crab_ggNtuplizer_spring15_ZJetsToNuNu_HT-600toInf_try5/151212_075101/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZnnJetsHT600toInf.txt"

 printf "\n" 

fi


#Z(ll)Gamma Jets
#------------
if [ ${doZllGJMC} = true ]
then
 printf "Making Z(ll)Gamma Jets MC\n" 
 ls /hdfs/store/user/jjbuch/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ggNtuplizer_spring15_ZLLGJets_try5/151212_075328/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZllGJets.txt"
fi

#Z(nn)Gamma Jets
#------------
if [ ${doZnnGJMC} = true ]
then
 printf "Making Z(nn)Gamma Jets MC\n" 
  ls /hdfs/store/user/jjbuch/ZNuNuGJets_MiniAOD/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets_MiniAOD/160315_232734/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_ZnnGJets.txt"
fi

#W(ln)Gamma Jets
#------------
if [ ${doWlnGJMC} = true ]
then
 printf "Making Z(ln)Gamma Jets MC\n" 
  ls /hdfs/store/user/jjbuch/WGJets_MiniAOD/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets_MiniAOD/160315_232658/0000/ggtree_mc_*root \
  > "${outdir}/hdfslist_WlnGJets.txt"
fi

#W(en)
#------------
if [ ${doWenMC} = true ]
then
 printf "Making W(en) MC\n" 

  ls /hdfs/store/user/jjbuch/WToENu_M-100_MiniAOD/WToENu_M-100_TuneCUETP8M1_13TeV-pythia8/crab_WToENu_M-100_MiniAOD/160413_214913/0000/ggtree_mc_DiPhoton_MINIAOD_TEST_*root \
  > "${outdir}/hdfslist_Wen.txt"
fi

#W(mn)
#------------
if [ ${doWmnMC} = true ]
then
 printf "Making W(mn) MC\n" 

  ls /hdfs/store/user/jjbuch/WToMuNu_M-100_MiniAOD/WToMuNu_M-100_TuneCUETP8M1_13TeV-pythia8/crab_WToMuNu_M-100_MiniAOD/160413_214950/0000/ggtree_mc_DiPhoton_MINIAOD_TEST_*root \
  > "${outdir}/hdfslist_Wmn.txt"
fi

#W(tn)
#------------
if [ ${doWtnMC} = true ]
then
 printf "Making W(tn) MC\n" 

  ls /hdfs/store/user/jjbuch/WToTauNu_M-100_MiniAOD/WToTauNu_M-100_TuneCUETP8M1_13TeV-pythia8-tauola/crab_WToTauNu_M-100_MiniAOD/160413_215033/0000/ggtree_mc_DiPhoton_MINIAOD_TEST_*root \
  > "${outdir}/hdfslist_Wtn.txt"
fi

#TTGJ
#------------
if [ ${doTTGJMC} = true ]
then
 printf "Making TTGJets MC\n" 

  ls /hdfs/store/user/jjbuch/TTGJets_MiniAOD/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets_MiniAOD/160413_215119/0000/*root \
  > "${outdir}/hdfslist_TTGJets.txt"
fi
