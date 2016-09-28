#!/bin/bash

outdir="${CMSSW_BASE}/src/LoneGamma/lists"
outname="initialEvents.txt"
outfile="${outdir}/${outname}"

doSinglePhoton=false
doDoubleElectron=false
doGJetsMC=false
 # GJetsHT40to100
 # GJetsHT100to200
 # GJetsHT200to400
 # GJetsHT400to600
 # GJetsHT600toInf
doQCDMC=false
 # QCDPt15to20
 # QCDPt20to30
 # QCDPt30to50
 # QCDPt50to80
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
 # WWGets
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
 printf "Counting SinglePhoton\n" 
 python eventCounter.py -i "${outdir}/hdfslist_SinglePhoton.txt" -o "${outfile}" \
 -s "SinglePhoton"
fi


#DoubleElectron
#------------
if [ ${doDoubleElectron} = true ]
then
 printf "Counting DoubleElectron\n" 
 python eventCounter.py -i "${outdir}/hdfslist_DoubleElectron.txt" -o "${outfile}" \
 -s "DoubleElectron"
fi


#GJets MC
#------------
if [ ${doGJetsMC} = true ]
then
 printf "Counting GJets " 
 printf "1  "
 python eventCounter.py -i "${outdir}/hdfslist_GJetsHT40to100.txt" -o "${outfile}" \
 -s "GJetsHT40to100"

 printf "2  "
 python eventCounter.py -i "${outdir}/hdfslist_GJetsHT100to200.txt" -o "${outfile}" \
 -s "GJetsHT100to200"

 printf "3  "
 python eventCounter.py -i "${outdir}/hdfslist_GJetsHT200to400.txt" -o "${outfile}" \
 -s "GJetsHT200to400"

 printf "4  "
 python eventCounter.py -i "${outdir}/hdfslist_GJetsHT400to600.txt" -o "${outfile}" \
 -s "GJetsHT400to600"

 printf "5  "
 python eventCounter.py -i "${outdir}/hdfslist_GJetsHT600toInf.txt" -o "${outfile}" \
 -s "GJetsHT600toInf"

 printf "\n" 
fi


#QCD MC
#----------
if [ ${doQCDMC} = true ]
then
 printf "Counting QCD " 
# printf "1  " 
#  python eventCounter.py -i "${outdir}/hdfslist_QCDPt15to20.txt" -o "${outfile}" \
#  -s "QCDPt15to20"
#
# printf "2  " 
#  python eventCounter.py -i "${outdir}/hdfslist_QCDPt20to30.txt" -o "${outfile}" \
#  -s "QCDPt20to30"
#
# printf "3  " 
#  python eventCounter.py -i "${outdir}/hdfslist_QCDPt30to50.txt" -o "${outfile}" \
#  -s "QCDPt30to50"
#
# printf "4  " 
#  python eventCounter.py -i "${outdir}/hdfslist_QCDPt50to80.txt" -o "${outfile}" \
#  -s "QCDPt50to80"

 printf "5  " 
  python eventCounter.py -i "${outdir}/hdfslist_QCDPt80to120.txt" -o "${outfile}" \
  -s "QCDPt80to120"

 printf "6  " 
  python eventCounter.py -i "${outdir}/hdfslist_QCDPt120to170.txt" -o "${outfile}" \
  -s "QCDPt120to170"

 printf "7  " 
  python eventCounter.py -i "${outdir}/hdfslist_QCDPt170to300.txt" -o "${outfile}" \
  -s "QCDPt170to300"

 printf "8  " 
  python eventCounter.py -i "${outdir}/hdfslist_QCDPt300toInf.txt" -o "${outfile}" \
  -s "QCDPt300toInf"

 printf "\n" 
fi


# Z(ll)Jets MC
#----------
if [ ${doZllJMC} = true ]
then
 printf "Counting Z(ll)Jets MC "
 printf "1  "
 python eventCounter.py -i "${outdir}/hdfslist_ZllJetsHT100to200.txt" -o "${outfile}" \
 -s "ZllJetsHT100to200"

 printf "2  "
 python eventCounter.py -i "${outdir}/hdfslist_ZllJetsHT200to400.txt" -o "${outfile}" \
 -s "ZllJetsHT200to400"

 printf "3  "
 python eventCounter.py -i "${outdir}/hdfslist_ZllJetsHT400to600.txt" -o "${outfile}" \
 -s "ZllJetsHT400to600"

 printf "4  "
 python eventCounter.py -i "${outdir}/hdfslist_ZllJetsHT600toInf.txt" -o "${outfile}" \
 -s "ZllJetsHT600toInf"

 printf "\n"
fi


# Z(ll)Jets MC
#----------
if [ ${doZnnJMC} = true ]
then
 printf "Counting Z(nn)Jets MC "
 printf "1  "
 python eventCounter.py -i "${outdir}/hdfslist_ZnnJetsHT100to200.txt" -o "${outfile}" \
 -s "ZnnJetsHT100to200"

 printf "2  "
 python eventCounter.py -i "${outdir}/hdfslist_ZnnJetsHT200to400.txt" -o "${outfile}" \
 -s "ZnnJetsHT200to400"

 printf "3  "
 python eventCounter.py -i "${outdir}/hdfslist_ZnnJetsHT400to600.txt" -o "${outfile}" \
 -s "ZnnJetsHT400to600"

 printf "4  "
 python eventCounter.py -i "${outdir}/hdfslist_ZnnJetsHT600toInf.txt" -o "${outfile}" \
 -s "ZnnJetsHT600toInf"

 printf "\n"
fi


#Z(nunu)Gamma Jets
#------------
if [ ${doZnnGJMC} = true ]
then
 printf "Counting Z(nunu)Gamma Jets\n"
 python eventCounter.py -i "${outdir}/hdfslist_ZnnGJets.txt" -o "${outfile}" \
 -s "ZnnGJets"
fi


#Z(ll)Gamma Jets
#------------
if [ ${doZllGJMC} = true ]
then
 printf "Counting Z(ll)Gamma Jets\n"
 python eventCounter.py -i "${outdir}/hdfslist_ZllGJets.txt" -o "${outfile}" \
 -s "ZllGJets"
fi


#W(ln)Gamma Jets
#------------
if [ ${doWlnGJMC} = true ]
then
 printf "Counting W(ln)Gamma Jets MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_WlnGJets.txt" -o "${outfile}" \
 -s "WlnGJets"
fi


#W(en)
#------------
if [ ${doWenMC} = true ]
then
 printf "Counting W(en) MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_Wen.txt" -o "${outfile}" \
 -s "Wen"
fi


#W(mn)
#------------
if [ ${doWmnMC} = true ]
then
 printf "Counting W(mn) MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_Wmn.txt" -o "${outfile}" \
 -s "Wmn"
fi


#W(tn)
#------------
if [ ${doWtnMC} = true ]
then
 printf "Counting W(tn) MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_Wtn.txt" -o "${outfile}" \
 -s "Wtn"
fi


#TTGJ
#------------
if [ ${doTTGJMC} = true ]
then
 printf "Counting TTGJets MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_TTGJets.txt" -o "${outfile}" \
 -s "TTGJets"
fi


#WWG
#------------
if [ ${doWWGMC} = true ]
then
 printf "Counting WWG MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_WWG.txt" -o "${outfile}" \
 -s "WWG"
fi


#GGJ
#------------
if [ ${doGGJMC} = true ]
then
 printf "Counting GGJets MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_GGJets.txt" -o "${outfile}" \
 -s "GGJets"
fi


#TGJ
#------------
if [ ${doTGJMC} = true ]
then
 printf "Counting TGJets MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_TGJets.txt" -o "${outfile}" \
 -s "TGJets"
fi


#WZ
#------------
if [ ${doWZMC} = true ]
then
 printf "Counting WZ MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_WZ.txt" -o "${outfile}" \
 -s "WZ"
fi


#ZZ
#------------
if [ ${doZZMC} = true ]
then
 printf "Counting ZZ MC\n" 
 python eventCounter.py -i "${outdir}/hdfslist_ZZ.txt" -o "${outfile}" \
 -s "ZZ"
fi
