#!/bin/bash

outdir="${CMSSW_BASE}/src/LoneGamma/lists"
mkdir -p "${outdir}"

doSinglePhoton=true
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
doZnnGJMC=false
 # ZnnGJets
doZllGJMC=false
 # ZllGJets


#SinglePhoton
#------------
if [ ${doSinglePhoton} = true ]
then
 printf "Making SinglePhoton\n" 
 ls  /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton*12p9*newbh_1/*/*/*root \
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
 printf "1  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-80to120/160621_133148/0000/*root \
  > "${outdir}/hdfslist_QCDPt80to120.txt"

 printf "2  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-120to170/160621_133201/0000/*root \
  > "${outdir}/hdfslist_QCDPt120to170.txt"

 printf "3  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-170to300/160621_133213/0000/*root \
  > "${outdir}/hdfslist_QCDPt170to300.txt"

 printf "4  " 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt-300toInf/160621_133225/0000/*root \
  > "${outdir}/hdfslist_QCDPt300toInf.txt"

 printf "\n" 
fi

# Z(ll)Jets MC
#----------
if [ ${doZllJMC} = true ]
then
 printf "Counting Z(ll)Jets MC " 
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

#Z(nunu)Gamma Jets
#------------
if [ ${doZnnGJMC} = true ]
then
 printf "Making Z(nunu)Gamma Jets MC\n" 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets/160621_132818/0000/*root \
  > "${outdir}/hdfslist_ZnnGJets.txt"
fi


#Z(ll)Gamma Jets
#------------
if [ ${doZllGJMC} = true ]
then
 printf "Making Z(ll)Gamma Jets MC\n" 
 ls /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/*root \
  > "${outdir}/hdfslist_ZllGJets.txt"
fi
