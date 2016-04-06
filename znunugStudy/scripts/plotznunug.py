#!/usr/bin/env python
'''
Plots Numerator distributions prefit
Author: T.M.Perry UW
'''
import ROOT
from ROOT import THStack,TH1F,TFile
from ROOT import TLegend,TCanvas,TPad,TLatex,TLine,TColor
from ROOT import gROOT,gStyle
import sys
import time
import cmsPrelim as cpr

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v","--version",        help="version name")
parser.add_argument("-id","--inpdir",        help="/full/path/to/input/directory")
parser.add_argument("-mn","--mc_filename",   help="mc filename Z(ll)Gamma (.root)")
parser.add_argument("-nn","--nc_filename",   help="mc filename Z(nunu)Jets (.root)")
parser.add_argument("-dn","--data_filename", help="data filename (.root)")
parser.add_argument("-od","--outdir",        help="/path/to/output/directory")
parser.add_argument("-on","--out_filename",  help="output name (no .root)")
parser.add_argument("-sf","--scalefactor",   help="scale factor for mc")
parser.add_argument("-var","--variable",     help="variable name")
parser.add_argument("-rng","--ptrange",      help="pT range (ex. 250to400")
parser.add_argument("-log","--do_log",  action='store_true', help="Set Log Scale")
args = parser.parse_args()

version        = args.version
inpdir         = args.inpdir        
mc_filename    = args.mc_filename   
nc_filename    = args.nc_filename   
data_filename  = args.data_filename   
outdir         = args.outdir  
out_filename   = args.out_filename  
#sf             = args.scalefactor
sf = 2320./2240.
variable       = args.variable
ptrange        = args.ptrange
do_log         = args.do_log
webdir = "/afs/hep.wisc.edu/user/tperry/www/MonoPhoton/znunugPlots/%s"%(version)

xn = ""
if(do_log): out_filename+="_log"

m_file = TFile("%s/%s"%(inpdir,mc_filename))
n_file = TFile("%s/%s"%(inpdir,nc_filename))
d_file = TFile("%s/%s"%(inpdir,data_filename))

#canvas attributes
canx = 1200
cany = 1100 
gStyle.SetOptStat('')
gStyle.SetLineWidth(3)
gStyle.SetPadTickY(1)

#color scheme
#Int_t ci = 1756; // color index
#TColor *color = new TColor(ci, 0.1, 0.2, 0.3);

c_ntot = ROOT.kBlack
c_nsig = ROOT.kRed
c_zllmc = ROOT.TColor.GetColor("#9CFDAD")
c_znnmc = ROOT.TColor.GetColor("#564DFB")
fs = 1001
# line style
n_ls = 1
d_ls = 7

# TLatex
tex = ROOT.TLatex()
tex.SetTextSize(0.07)
tex.SetTextAlign(13)
tex.SetNDC(True)

rebin = 1

c = TCanvas('c','Canvas Named c',canx,cany)
c.cd()
c.SetLogy(do_log)
c.SetFrameLineWidth(2)
c.SetCanvasSize(canx,cany)

# Data
h_data = d_file.Get("h_sig_%s_%s"%(variable,ptrange))
h_data.SetName("h_data")
h_data.SetTitle("")
h_data.SetLineColor(c_ntot)
h_data.SetLineWidth(4)
h_data.SetFillStyle(0)
h_data.SetLineStyle(n_ls)
h_data.Draw("GOFF")
#h_data.Scale( 1. / max(h_data.Integral(-1,-1),1.) )
themax = h_data.GetMaximum()

# MC (ZG->llG)
h_zllmc = m_file.Get("h_sig_%s_%s"%(variable,ptrange))
h_zllmc.SetName("h_zllmc")
h_zllmc.SetTitle("")
h_zllmc.SetLineColor(1)
h_zllmc.SetLineWidth(2)
h_zllmc.SetFillStyle(fs)
h_zllmc.SetLineStyle(n_ls)
h_zllmc.SetFillColor(c_zllmc)
h_zllmc.Scale(sf)
h_zllmc.Draw("GOFF")
themax = max( themax, h_zllmc.GetMaximum() )

# MC (ZG->nnJ)
h_znnmc = n_file.Get("h_sig_%s_%s"%(variable,ptrange))
h_znnmc.SetName("h_znnmc")
h_znnmc.SetTitle("")
h_znnmc.SetLineColor(1)
h_znnmc.SetLineWidth(2)
h_znnmc.SetFillStyle(fs)
h_znnmc.SetLineStyle(n_ls)
h_znnmc.SetFillColor(c_znnmc)
h_znnmc.Scale(sf)
h_znnmc.Draw("GOFF")
themax = max( themax, h_znnmc.GetMaximum() )

s_mc = THStack('s_mc','')
s_mc.Add(h_znnmc)
s_mc.Add(h_zllmc)

# Nr Entries MC (ZG->llG)
h_zllent = m_file.Get("h_sig_%s_%s"%(variable,ptrange))
h_zllent.SetName("h_zllent")
nrentries = h_zllent.GetEntries()
h_zllent.SetTitle("")
h_zllent.SetLineStyle(n_ls)
h_zllent.Scale( nrentries / max(h_zllent.Integral(),1.) )
h_zllent.SetLineColor(1)
h_zllent.SetLineWidth(2)
h_zllent.SetFillStyle(fs)
h_zllent.SetLineStyle(n_ls)
h_zllent.SetFillColor(c_zllmc)
h_zllent.Draw("GOFF")
thegmax = h_zllent.GetMaximum()

# GEN Level (ZG->llG)
h_zllgen = m_file.Get("h_gen_%s_%s"%(variable,ptrange))
h_zllgen.SetName("h_zllgen")
h_zllgen.SetTitle("")
h_zllgen.SetLineColor(c_ntot)
h_zllgen.SetLineWidth(4)
h_zllgen.SetFillStyle(0)
h_zllgen.Draw("GOFF")
thegmax = max( thegmax, h_zllgen.GetMaximum() )

# fill legends
leg=TLegend(0.55,0.70,0.88,0.88)
leg.AddEntry(h_data,"Data")
leg.AddEntry(h_zllmc,"Z(ll)G MC","f")
leg.AddEntry(h_znnmc,"Z(nn)Jets MC","f")
leg.SetFillColor(0)
leg.SetBorderSize(0)

legent=TLegend(0.55,0.70,0.88,0.88)
legent.AddEntry(h_zllgen,"Gen Level")
legent.AddEntry(h_zllent,"Reco Level","f")
legent.SetFillColor(0)
legent.SetBorderSize(0)

# and draw


####### Data v MC ##############
h_data.SetMaximum( 1.4*themax )

h_data.Draw("GOFF")
#h_zllmc.Scale(h_data.Integral()/h_zllmc.Integral())
#h_zllmc.Scale(0.041*1.3/117.864)

print("Data:  %.3f"%h_data.Integral())
print("MC:    %.3f"%h_zllmc.Integral())

h_data.GetYaxis().SetTitle("Events / 2 GeV")
#h_data.GetXaxis().SetTitle("dilepton mass (#mu#mu)")
h_data.GetXaxis().SetTitle("dilepton mass (#mu#mu,ee)")

#if(do_log):  h_data.SetMaximum( 50*themax ) 
if(do_log):  h_data.SetMaximum( 11 ) 
c.cd()
h_data.Draw('GOFF')
s_mc.Draw('hist,GOFF,sames')
#h_zllmc.Draw('hist,GOFF,sames')
leg.Draw('sames,GOFF')
h_data.Draw('GOFF,sames')

cpr.prelim_tdr(lumi=2320,hor=0.17)
tex.SetTextSize(0.03)
tex.SetTextAlign(11) #left, bottom
#tex.DrawLatex(0.13,0.75,"pT Range [GeV]: %s"%(ptrange))

c.Update()
time.sleep(1)

c.SaveAs("%s/%s.pdf"%(outdir,out_filename))
c.SaveAs("%s/%s.pdf"%(webdir,out_filename))
print(  "%s/%s.pdf"%(outdir,out_filename))



###### Gen v Raw Reco ##############
h_zllgen.SetMaximum( 1.4*thegmax )

h_zllgen.Draw("GOFF")
#h_zllmc.Scale(h_data.Integral()/h_zllmc.Integral())
#h_zllmc.Scale(0.041*1.3/117.864)

print("Gen:         %.3f"%h_zllgen.Integral())
print("Raw Reco:    %.3f"%h_zllent.Integral())

h_zllgen.GetYaxis().SetTitle("Events / 2 GeV")
#h_zllgen.GetXaxis().SetTitle("dilepton mass (#mu#mu)")
h_zllgen.GetXaxis().SetTitle("dilepton mass (#mu#mu,ee)")

c.cd()
h_zllgen.Draw('GOFF')
h_zllent.Draw('hist,GOFF,sames')
legent.Draw('sames,GOFF')
h_zllgen.Draw('GOFF,sames')

cpr.prelim_tdr(lumi=2240,hor=0.17)
tex.SetTextSize(0.03)
tex.SetTextAlign(11) #left, bottom
#tex.DrawLatex(0.13,0.75,"pT Range [GeV]: %s"%(ptrange))

c.Update()
time.sleep(1)

c.SaveAs("%s/gen%s.pdf"%(outdir,out_filename))
c.SaveAs("%s/gen%s.pdf"%(webdir,out_filename))
print(  "%s/gen%s.pdf"%(outdir,out_filename))

