#!/usr/bin/python
'''
adds up the total number of evnts run over MC in a list of files
pulling info from summary/results

usage:
python eventCounter.py list_of_file_locations.txt

Author: T.M.Perry UW-Madison
'''

from sys import argv, stdout, stderr
import ROOT
import math

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-s","--samplename",  help="name of sample")
parser.add_argument("-i","--inlist",      help="input list of files (full path)")
parser.add_argument("-o","--outfilename", help="output filename (full path)")
args = parser.parse_args()

samplename    = args.samplename  
inlist        = args.inlist
outfilename   = args.outfilename


#Returns the number of lines in the files
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

counter = 0
total_events = 0

# count number of events
numfiles = file_len(inlist)
with open(inlist) as f:
 for filename in f:
  filename = filename.rstrip()
  ntuple_file = ROOT.TFile(filename)
  summary = ntuple_file.Get("ggNtuplizer/hEvents")
  nr_events = summary.GetBinContent(1)
  total_events+=nr_events
  counter = counter+1

# read cross section from list 
# a https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
# b https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
# c mcm
 ## https://cms-pdmv.cern.ch/mcm/requests?produce=%2FZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph%2FRunIIWinter15GS-MCRUN2_71_V1-v1%2FGEN-SIM&page=0&shown=262271
if    samplename == "SinglePhoton"      : xc="1."
elif  samplename == "DoubleElectron"    : xc="1."
elif  samplename == "GJetsHT40to100"    : xc="20790."     # a LO  #Gamma_jets
elif  samplename == "GJetsHT100to200"   : xc="9238."      # a LO  #Gamma_jets  
elif  samplename == "GJetsHT200to400"   : xc="2305."      # a LO  #Gamma_jets  
elif  samplename == "GJetsHT400to600"   : xc="274.4"      # a LO  #Gamma_jets  
elif  samplename == "GJetsHT600toInf"   : xc="93.46"      # a LO  #Gamma_jets
elif  samplename == "QCDPt15to20"       : xc="2302200.0"  # a LO  #QCD 1279000000 * 0.0018
elif  samplename == "QCDPt20to30"       : xc="5352960.0"  # a LO  #QCD  557600000 * 0.0096
elif  samplename == "QCDPt30to50"       : xc="9928000.0"  # a LO  #QCD  136000000 * 0.073
elif  samplename == "QCDPt50to80"       : xc="2890800.0"  # a LO  #QCD  19800000 * 0.146
elif  samplename == "QCDPt80to120"      : xc="350000.0"   # a LO  #QCD  2800000 * 0.125
elif  samplename == "QCDPt120to170"     : xc="62964.0"    # a LO  #QCD  477000 * 0.132
elif  samplename == "QCDPt170to300"     : xc="18810.0"    # a LO  #QCD  114000 * 0.165
elif  samplename == "QCDPt300toInf"     : xc="1350.0"     # a LO  #QCD  9000 * 0.15
elif  samplename == "ZllGJets"          : xc="0.1859"     # c   1.3 * 0.143 
elif  samplename == "ZnnGJets"          : xc="0.1903"     # ?    k-factors
elif  samplename == "ZllJetsHT100to200" : xc="181.302"    # a    147.40*1.23   
elif  samplename == "ZllJetsHT200to400" : xc="50.418"     # a    40.99*1.23   
elif  samplename == "ZllJetsHT400to600" : xc="6.984"      # a    5.678*1.23
elif  samplename == "ZllJetsHT600toInf" : xc="2.704"      # a    2.198*1.23
elif  samplename == "ZnnJetsHT100to200" : xc="344.83"     # a    280.35*1.23   
elif  samplename == "ZnnJetsHT200to400" : xc="95.53"      # a    77.67*1.23   
elif  samplename == "ZnnJetsHT400to600" : xc="13.20"      # a    10.73*1.23
elif  samplename == "ZnnJetsHT600toInf" : xc="5.06"       # a    4.116*1.23
elif  samplename == "WlnGJets"          : xc="0.889"      # "0.6637" x 1.34 
#elif  samplename == "Wen"               : xc="174."       # 20508.9
elif  samplename == "Wmn"               : xc="174."       # 20508.9
elif  samplename == "Wtn"               : xc="165."       # 20508.9
elif  samplename == "TTGJets"           : xc="3.697"
elif  samplename == "WWG"               : xc="0.2147"     # NLO
elif  samplename == "GGJets"            : xc="135.1"      # b aMC@NLO 
elif  samplename == "TGJets"            : xc="2.967"      # b NLO
elif  samplename == "WZ"                : xc="47.13"      # b NLO MCFM
elif  samplename == "ZZ"                : xc="16.523"     # b 


outfile = open(outfilename,'a')
outfile.write(" %s Events: %s\n"%(samplename,int(total_events)))
outfile.write(" %s XC: %s\n"%(samplename,xc))
outfile.close()
