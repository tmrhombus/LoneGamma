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

#Returns the number of lines in the files
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

counter = 0
total_events = 0

samplename = argv[1]
list_of_files = argv[2]
outfile_name = argv[3]

# count number of events
numfiles = file_len(list_of_files)
with open(list_of_files) as f:
 for filename in f:
  filename = filename.rstrip()
  ntuple_file = ROOT.TFile(filename)
  summary = ntuple_file.Get("ggNtuplizer/hEvents")
  nr_events = summary.GetBinContent(1)
  total_events+=nr_events
  counter = counter+1

# read cross section from list 
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
# https://cms-pdmv.cern.ch/mcm/requests?produce=%2FZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph%2FRunIIWinter15GS-MCRUN2_71_V1-v1%2FGEN-SIM&page=0&shown=262271
if    samplename == "GJets_HT40To100"  : xc="20790."     # LO
elif  samplename == "GJets_HT100To200" : xc="9238." 
elif  samplename == "GJets_HT200To400" : xc="2305." 
elif  samplename == "GJets_HT400To600" : xc="274.4"
elif  samplename == "GJets_HT600ToInf" : xc="93.46"
elif  samplename == "QCD_Pt15to20"     : xc="2302200.0"  # LO  1279000000 * 0.0018
elif  samplename == "QCD_Pt20to30"     : xc="5352960.0"  # 557600000 * 0.0096
elif  samplename == "QCD_Pt30to50"     : xc="9928000.0"  # 136000000 * 0.073
elif  samplename == "QCD_Pt50to80"     : xc="2890800.0"  # 19800000 * 0.146
elif  samplename == "QCD_Pt80to120"    : xc="350000.0"   # 2800000 * 0.125
elif  samplename == "QCD_Pt120to170"   : xc="62964.0"    # 477000 * 0.132
elif  samplename == "QCD_Pt170to300"   : xc="18810.0"    # 114000 * 0.165
elif  samplename == "QCD_Pt300toInf"   : xc="1350.0"     # 9000 * 0.15
elif  samplename == "ZLLG"             : xc="0.1859"     # 1.3 * 0.143 
elif  samplename == "ZNuNuG"           : xc="0.2899"     # 1.3 * 0.223
elif  samplename == "ZJetsToNuNu100To200"    : xc="181.302"    # 147.40*1.23   
elif  samplename == "ZJetsToNuNu200To400"    : xc="50.418"     # 40.99*1.23   
elif  samplename == "ZJetsToNuNu400To600"    : xc="6.984"      # 5.678*1.23
elif  samplename == "ZJetsToNuNu600ToInf"    : xc="2.704"      # 2.198*1.23

outfile = open(outfile_name,'a')
outfile.write("%s Events: %s\n"%(samplename,int(total_events)))
outfile.write("%s XC: %s\n"%(samplename,xc))
outfile.close()
