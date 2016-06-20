#!/bin/bash

outdir="${CMSSW_BASE}/src/LoneGamma/lists"
mkdir -p "${outdir}"

doSinglePhoton=true
doDoubleElectron=true
doGJetsMC=true
 # GJetsHT40to100
 # GJetsHT100to200
 # GJetsHT200to400
 # GJetsHT400to600
 # GJetsHT600toInf
doQCDMC=true
 # QCDPt80to120
 # QCDPt120to170
 # QCDPt170to300
 # QCDPt300toInf
doZllJMC=true
 # ZllJetsHT100to200
 # ZllJetsHT200to400
 # ZllJetsHT400to600
 # ZllJetsHT600toInf
doZnnGJMC=true
 # ZnnGJets
doZllGJMC=true
 # ZllGJets


#SinglePhoton
#------------
if [ ${doSinglePhoton} = true ]
then
 printf "Making SinglePhoton\n" 
 ls /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton_13TeV_2016_2p6fb_allvar/160618_232516/*/*root \
  > "${outdir}/hdfslist_SinglePhoton.txt"
fi


#DoubleElectron
#------------
if [ ${doDoubleElectron} = true ]
then
 printf "Making DoubleElectron\n" 
 ls /hdfs/store/user/gomber/SinglePhoton_2016B/DoubleEG/crab_job_DoubleEG_13TeV_2016_2fb/160610_135127/*/*root \
  > "${outdir}/hdfslist_DoubleElectron.txt"
fi

#GJets MC
#------------
if [ ${doGJetsMC} = true ]
then
 printf "Making GJets " 
 printf "1  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/GJets_HT-40To100_MiniAOD/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-40To100_MiniAOD/160525_113910/0000/*root \
  > "${outdir}/hdfslist_GJetsHT40to100.txt"

 printf "2  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/GJets_HT-100To200_MiniAOD/MiniAODv2/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-100To200_MiniAOD_MiniAODv2/160606_151142/0000/*root \
  > "${outdir}/hdfslist_GJetsHT100to200.txt"

 printf "3  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/GJets_HT-200To400_MiniAOD/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-200To400_MiniAOD/160525_114055/0000/*root \
  > "${outdir}/hdfslist_GJetsHT200to400.txt"

 printf "4  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/GJets_HT-400To600_MiniAOD/MiniAODv2/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-400To600_MiniAOD_MiniAODv2/160610_130056/0000/*root \
  > "${outdir}/hdfslist_GJetsHT400to600.txt"

 printf "5  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/GJets_HT-600ToInf_MiniAOD/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-600ToInf_MiniAOD/160525_114545/0000/*root \
  > "${outdir}/hdfslist_GJetsHT600toInf.txt"

 printf "\n" 
fi

#QCD MC
#----------
if [ ${doQCDMC} = true ]
then
 printf "Making QCD " 
 printf "1  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/QCD_Pt-80to120_MiniAOD/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-80to120_MiniAOD/160603_163503/0000/*root \
  > "${outdir}/hdfslist_QCDPt80to120.txt"

 printf "2  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/QCD_Pt-120to170_MiniAOD/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-120to170_MiniAOD/160527_121408/000*/*root \
  > "${outdir}/hdfslist_QCDPt120to170.txt"

 printf "3  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/QCD_Pt-170to300_MiniAOD/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-170to300_MiniAOD/160525_115107/0000/*root \
  > "${outdir}/hdfslist_QCDPt170to300.txt"

 printf "4  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/QCD_Pt-300toInf_MiniAOD/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-300toInf_MiniAOD/160525_115146/0000/*root \
  > "${outdir}/hdfslist_QCDPt300toInf.txt"

 printf "\n" 
fi

# Z(ll)Jets MC
#----------
if [ ${doZllJMC} = true ]
then
 printf "Counting Z(ll)Jets MC " 
 printf "1  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160612_134238/0000/*root \
  > "${outdir}/hdfslist_ZllJetsHT100to200.txt"

 printf "2  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160612_134311/0000/*root \
  > "${outdir}/hdfslist_ZllJetsHT200to400.txt"

 printf "3  "
 ls /hdfs/store/user/jjbuch/ntuples80X/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160612_134340/0000/*root \
  > "${outdir}/hdfslist_ZllJetsHT400to600.txt"

 printf "4  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MiniAODv2/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/160612_134401/0000/*root \
  > "${outdir}/hdfslist_ZllJetsHT600toInf.txt"

 printf "\n" 
fi

#Z(nunu)Gamma Jets
#------------
if [ ${doZnnGJMC} = true ]
then
 printf "Making Z(nunu)Gamma Jets MC\n" 
 ls /hdfs/store/user/jjbuch/ntuples80X/ZNuNuGJets_MiniAOD/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets_MiniAOD/160525_112045/0000/*root \
  > "${outdir}/hdfslist_ZnnGJets.txt"
fi

#Z(ll)Gamma Jets
#------------
if [ ${doZllGJMC} = true ]
then
 printf "Making Z(ll)Gamma Jets MC\n" 
 ls /hdfs/store/user/jjbuch/ntuples80X/ZLLGJets_MiniAOD/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets_MiniAOD/160525_112640/0000/*root \
  > "${outdir}/hdfslist_ZllGJets.txt"
fi
