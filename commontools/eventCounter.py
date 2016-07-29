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
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
# https://cms-pdmv.cern.ch/mcm/requests?produce=%2FZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph%2FRunIIWinter15GS-MCRUN2_71_V1-v1%2FGEN-SIM&page=0&shown=262271
if    samplename == "SinglePhoton"      : xc="1."
elif  samplename == "DoubleElectron"    : xc="1."
elif  samplename == "GJetsHT40to100"    : xc="20790."     # LO
elif  samplename == "GJetsHT100to200"   : xc="9238." 
elif  samplename == "GJetsHT200to400"   : xc="2305." 
elif  samplename == "GJetsHT400to600"   : xc="274.4"
elif  samplename == "GJetsHT600toInf"   : xc="93.46"
elif  samplename == "QCDPt15to20"       : xc="2302200.0"  # LO  1279000000 * 0.0018
elif  samplename == "QCDPt20to30"       : xc="5352960.0"  # 557600000 * 0.0096
elif  samplename == "QCDPt30to50"       : xc="9928000.0"  # 136000000 * 0.073
elif  samplename == "QCDPt50to80"       : xc="2890800.0"  # 19800000 * 0.146
elif  samplename == "QCDPt80to120"      : xc="350000.0"   # 2800000 * 0.125
elif  samplename == "QCDPt120to170"     : xc="62964.0"    # 477000 * 0.132
elif  samplename == "QCDPt170to300"     : xc="18810.0"    # 114000 * 0.165
elif  samplename == "QCDPt300toInf"     : xc="1350.0"     # 9000 * 0.15
elif  samplename == "ZllGJets"          : xc="0.1859"     # 1.3 * 0.143 
elif  samplename == "ZnnGJets"          : xc="1."
 # 175 [14.4] 190 [29.7] 250 [17.5] 400 [3.7] 700 [0.3] Inf
elif  samplename == "ZllJetsHT100to200" : xc="181.302"    # 147.40*1.23   
elif  samplename == "ZllJetsHT200to400" : xc="50.418"     # 40.99*1.23   
elif  samplename == "ZllJetsHT400to600" : xc="6.984"      # 5.678*1.23
elif  samplename == "ZllJetsHT600toInf" : xc="2.704"      # 2.198*1.23
elif  samplename == "ZnnJetsHT100to200" : xc="344.83"     # 280.35*1.23   
elif  samplename == "ZnnJetsHT200to400" : xc="95.53"      # 77.67*1.23   
elif  samplename == "ZnnJetsHT400to600" : xc="13.20"      # 10.73*1.23
elif  samplename == "ZnnJetsHT600toInf" : xc="5.06"       # 4.116*1.23
elif  samplename == "WlnGJets"          : xc="243.9"      # 
elif  samplename == "Wen"               : xc="20508.9"
elif  samplename == "Wmn"               : xc="20508.9"
elif  samplename == "Wtn"               : xc="20508.9"
elif  samplename == "TTGJets"           : xc="3.697"



outfile = open(outfilename,'a')
outfile.write("%s Events: %s\n"%(samplename,int(total_events)))
outfile.write("%s XC: %s\n"%(samplename,xc))
outfile.close()
