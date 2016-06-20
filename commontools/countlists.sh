#!/bin/sh

outdir="${CMSSW_BASE}/src/LoneGamma/lists"
outname="initialEvents.txt"
outfile="${outdir}/${outname}"
#if [ ! -a ${outfile} ] 
#then printf  "" > $outfile fi

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
 printf "1  " 
 python eventCounter.py -i "${outdir}/hdfslist_QCDPt80to120.txt" -o "${outfile}" \
 -s "QCDPt80to120"
 
 printf "2  " 
 python eventCounter.py -i "${outdir}/hdfslist_QCDPt120to170.txt" -o "${outfile}" \
 -s "QCDPt120to170"
 
 printf "3  " 
 python eventCounter.py -i "${outdir}/hdfslist_QCDPt170to300.txt" -o "${outfile}" \
 -s "QCDPt170to300"
 
 printf "4  " 
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

