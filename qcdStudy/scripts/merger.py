#!/usr/bin/env python
'''
Merge Histograms
Author: T.M.Perry
'''
import ROOT
from ROOT import TH1F,TFile,gROOT, TCanvas, TLegend, gStyle
import sys 
import fnmatch as fnm 
import time
import os


infilename =  "mrg4bins_SinglePhoton" # "mrg4bins_GJets" # "mrg4bins_QCD" # "analyzed_QCD_Merged" # "analyzed_GJets_Merged" # "analyzed_SinglePhoton"  # 
outfilename = "mrg3bins_SinglePhoton" # "mrg3bins_GJets" # "mrg3bins_QCD" # "mrg4bins_QCD"        # "mrg4bins_GJets"        # "mrg4bins_SinglePhoton"  # 

rebin = 1
path = os.environ.get('analyzed')

gStyle.SetOptStat('')
gStyle.SetLineWidth(3)
gStyle.SetPadTickY(1)

infile = TFile("%s/%s.root"%(path,infilename))

canx = 900 
cany = 900
c = TCanvas('c','Canvas Named c',canx,cany)

log = open('%s/%s.log'%(path,outfilename),'w')
log.write('\n\n') 
log.write('  infilename %s \n'%( infilename ) ) 
log.write('\n') 

outfile=gROOT.FindObject('%s/%s.root'%(path,outfilename))
if outfile : outFile.Close()
outfile = TFile('%s/%s.root'%(path,outfilename),'RECREATE','QCD from mT<30')

for key in infile.GetListOfKeys():
 obj = key.ReadObj()
 if(obj.IsA().InheritsFrom("TH1")):
  hname = obj.GetName() 
  if fnm.fnmatch(hname,"*250to400*"):    # 250 to 400 
   hname = hname.replace("400","1000")  # 250 to 1000
  #if fnm.fnmatch(hname,"*400to700*"):   # 400 to 700
  # hname = hname.replace("700","1000")  # 400 to 1000
   obj.SetName(hname)
   tbaname = ""
   tbaname = hname.replace("250","400")  # 400 to 1000
   #tbaname = hname.replace("400","700")  # 700 to 1000
   obj.Add(infile.Get(tbaname).Clone())

  if fnm.fnmatch(hname,"*400to1000*"): continue
  #if fnm.fnmatch(hname,"*700to1000*"): continue
  
  if fnm.fnmatch(hname,"h_*_sieieF5x5_*0"): 
   pdfname=outfilename+"_"+hname
   obj_size_inp = obj.Integral(-1,-1)
   obj_max = obj.GetMaximum()
   obj.SetFillColor(0)
   obj.SetLineColor(1)
   obj.SetLineWidth(3)

   #leg=TLegend(0.58,0.7,0.78,0.88)
   #leg.SetFillColor(0)
   #leg.SetBorderSize(0)
   #leg.AddEntry(obj,"Wisconsin")

   obj.GetXaxis().SetLabelSize(0.03)
   obj.GetYaxis().SetLabelSize(0.03)
   obj.GetYaxis().SetTitleOffset(1.5)
   obj.GetYaxis().SetTitle( "Events / %s GeV"%(obj.GetBinWidth(1)) )
   obj.SetMaximum(1.2*obj_max)
   obj.GetXaxis().SetTitle( "#sigma i#eta i#eta" )
   obj.SetTitle( hname )

   obj.Draw("hist,GOFF")
   #leg.Draw('sames')
   time.sleep(1)
   c.Print(path+"/"+pdfname+".pdf")

  obj.Write()

#  log.write("  %s \t & %.1f \t & %.1f\n"%( hname,obj_size_inp,obj_size_scl ) ) 
#
