#!/usr/bin/env python
'''
rescale all histograms (for lumi)
Author: T.M.Perry
'''
import ROOT
from ROOT import TH1F,TFile



samples = [
"ZllGJets",
"ZllJets",
"ZnnGJets"
]

path = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_530_CMSSW_8_0_8/src/LoneGamma/znunugStudy/gitignore/Erste/analyzed"

for sample in samples:
 inname  = path+"/analyzed_"+sample
 outname = path+"/rwt26_"+sample
 
 infile  = TFile(inname+".root")
 outfile = TFile(outname+'.root','RECREATE','rescaled MC')
 
 for key in infile.GetListOfKeys():
  obj = key.ReadObj()
  if(obj.IsA().InheritsFrom("TH1")):
   h_name = obj.GetName()
   h_new = obj.Clone()
   h_new.Scale( 2600./2200. )
   h_new.Write()
 
 infile.Close()
 outfile.Close()
