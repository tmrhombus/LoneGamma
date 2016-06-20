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
parser.add_argument("-ln","--lc_filename",   help="mc filename Z(ll)Jets (.root)")
parser.add_argument("-dn","--data_filename", help="data filename (.root)")
parser.add_argument("-od","--outdir",        help="/path/to/output/directory")
parser.add_argument("-on","--out_filename",  help="output name (no .root)")
parser.add_argument("-var","--variable",     help="variable name")
parser.add_argument("-sel","--selection",    help="selection name")
parser.add_argument("-rng","--ptrange",      help="pT range (ex. 250to400")
parser.add_argument("-log","--do_log",  action='store_true', help="Set Log Scale")
args = parser.parse_args()

version        = args.version
inpdir         = args.inpdir        
mc_filename    = args.mc_filename   
nc_filename    = args.nc_filename   
lc_filename    = args.lc_filename   
data_filename  = args.data_filename   
outdir         = args.outdir  
out_filename   = args.out_filename  
sf = 1
variable       = args.variable
selection      = args.selection
ptrange        = args.ptrange
do_log         = args.do_log
webdir = "/afs/hep.wisc.edu/user/tperry/www/MonoPhoton/znunugPlots/%s"%(version)
www = "http://www.hep.wisc.edu/~tperry/MonoPhoton/znunugPlots/%s/png"%(version)

xn = ""
if(do_log): out_filename+="_log"

m_file = TFile("%s/%s"%(inpdir,mc_filename))
n_file = TFile("%s/%s"%(inpdir,nc_filename))
l_file = TFile("%s/%s"%(inpdir,lc_filename))
d_file = TFile("%s/%s"%(inpdir,data_filename))

outfile=gROOT.FindObject(outdir+"/"+out_filename+'.root')
if outfile : outFile.Close()
outfile = TFile(outdir+"/"+out_filename+'.root','RECREATE','z(ll)g histograms')

#canvas attributes
canx = 1200
cany = 1100 
gStyle.SetOptStat('')
gStyle.SetLineWidth(3)
gStyle.SetPadTickY(1)

c_ntot = ROOT.kBlack
c_zllmc = ROOT.kOrange-4                    #zg->llg
c_zlljmc = ROOT.TColor.GetColor("#FDABCC")  #zj->llj
c_znnmc = ROOT.kRed-10                      #zj->nnj
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
data_obs = d_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
data_obs.SetName("data_obs")
data_obs.SetTitle("")
data_obs.SetLineColor(c_ntot)
data_obs.SetLineWidth(4)
data_obs.SetFillStyle(0)
data_obs.SetLineStyle(n_ls)
data_obs.Draw("GOFF")
#data_obs.Scale( 1. / max(data_obs.Integral(-1,-1),1.) )
data_obs.Write()
themax = data_obs.GetMaximum()

# MC (ZG->llG)
h_zllgmc = m_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
h_zllgmc.SetName("h_zllgmc")
h_zllgmc.SetTitle("")
h_zllgmc.SetLineColor(1)
h_zllgmc.SetLineWidth(2)
h_zllgmc.SetFillStyle(fs)
h_zllgmc.SetLineStyle(n_ls)
h_zllgmc.SetFillColor(c_zllmc)
h_zllgmc.Scale(sf)
h_zllgmc.Draw("GOFF")
h_zllgmc.Write()
themax = max( themax, h_zllgmc.GetMaximum() )

# MC (ZG->nnJ)
h_znnmc = n_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
h_znnmc.SetName("h_znnmc")
h_znnmc.SetTitle("")
h_znnmc.SetLineColor(1)
h_znnmc.SetLineWidth(2)
h_znnmc.SetFillStyle(fs)
h_znnmc.SetLineStyle(n_ls)
h_znnmc.SetFillColor(c_znnmc)
h_znnmc.Scale(sf)
h_znnmc.Draw("GOFF")
h_znnmc.Write()
themax = max( themax, h_znnmc.GetMaximum() )

# MC (ZJ->llJ)
h_zlljmc = l_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
h_zlljmc.SetName("h_zlljmc")
h_zlljmc.SetTitle("")
h_zlljmc.SetLineColor(1)
h_zlljmc.SetLineWidth(2)
h_zlljmc.SetFillStyle(fs)
h_zlljmc.SetLineStyle(n_ls)
h_zlljmc.SetFillColor(c_zlljmc)
h_zlljmc.Scale(sf)
h_zlljmc.Draw("GOFF")
h_zlljmc.Write()
themax = max( themax, h_zlljmc.GetMaximum() )

s_mc = THStack('s_mc','')
s_mc.Add(h_znnmc)
s_mc.Add(h_zlljmc)
s_mc.Add(h_zllgmc)

# Nr Entries MC (ZG->llG)
h_zllent = m_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
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
h_zllgen = m_file.Get("h_gen_%s_%s_%s"%(variable,ptrange,selection))
h_zllgen.SetName("h_zllgen")
h_zllgen.SetTitle("")
h_zllgen.SetLineColor(c_ntot)
h_zllgen.SetLineWidth(4)
h_zllgen.SetFillStyle(0)
h_zllgen.Draw("GOFF")
thegmax = max( thegmax, h_zllgen.GetMaximum() )

# fill legends
leg=TLegend(0.55,0.70,0.88,0.88)
leg.AddEntry(data_obs,"Data")
leg.AddEntry(h_zllgmc,"Z(ll)#gamma MC","f")
leg.AddEntry(h_zlljmc,"Z(ll)Jets MC","f")
#leg.AddEntry(h_znnmc,"Z(nn)Jets MC","f")
leg.SetFillColor(0)
leg.SetBorderSize(0)

legent=TLegend(0.55,0.70,0.88,0.88)
legent.AddEntry(h_zllgen,"Gen Level")
legent.AddEntry(h_zllent,"Reco Level","f")
legent.SetFillColor(0)
legent.SetBorderSize(0)

# and draw


####### Data v MC ##############
data_obs.SetMaximum( 5 )
#data_obs.SetMaximum( 1.4*themax )

data_obs.Draw("GOFF")
#h_zllgmc.Scale(data_obs.Integral()/h_zllgmc.Integral())
#h_zllgmc.Scale(0.041*1.3/117.864)

print("Data:  %.3f"%data_obs.Integral())
print("MC:    %.3f"%h_zllgmc.Integral())
print("Bkg:   %.3f"%h_zlljmc.Integral())

data_obs.GetYaxis().SetTitle("Events / 2 GeV")
#data_obs.GetXaxis().SetTitle("dilepton mass (#mu#mu)")
data_obs.GetXaxis().SetTitle("dilepton mass (#mu#mu,ee)")

#if(do_log):  data_obs.SetMaximum( 50*themax ) 
if(do_log):  data_obs.SetMaximum( 11 ) 
c.cd()
data_obs.Draw('GOFF')
s_mc.Draw('hist,GOFF,sames')
#h_zlljmc.Draw('hist,GOFF,sames')
#h_zllgmc.Draw('hist,GOFF,sames')
leg.Draw('sames,GOFF')
data_obs.Draw('GOFF,sames')

cpr.prelim_tdr(lumi=2200,hor=0.17)
tex.SetTextSize(0.03)
tex.SetTextAlign(11) #left, bottom
#tex.DrawLatex(0.13,0.75,"pT Range [GeV]: %s"%(ptrange))

c.Update()
time.sleep(1)

c.SaveAs("%s/%s.pdf"%(outdir,out_filename))
#c.SaveAs("%s/pdf/%s.pdf"%(webdir,out_filename))
#c.SaveAs("%s/png/%s.png"%(webdir,out_filename))
print(  "%s/%s.png"%(www,out_filename))


outfile.Close()


####### Gen v Raw Reco ##############
#h_zllgen.SetMaximum( 1.4*thegmax )
#
#h_zllgen.Draw("GOFF")
##h_zllgmc.Scale(data_obs.Integral()/h_zllgmc.Integral())
##h_zllgmc.Scale(0.041*1.3/117.864)
#
#print("Gen:         %5.1f"%h_zllgen.Integral())
#print("Raw Reco:    %5.1f"%h_zllent.Integral())
#print("Raw Reco NE: %5.1f"%nrentries)
#
#h_zllgen.GetYaxis().SetTitle("Events / 2 GeV")
##h_zllgen.GetXaxis().SetTitle("dilepton mass (#mu#mu)")
#h_zllgen.GetXaxis().SetTitle("dilepton mass (#mu#mu,ee)")
#
#c.cd()
#h_zllgen.Draw('GOFF')
#h_zllent.Draw('hist,GOFF,sames')
#legent.Draw('sames,GOFF')
#h_zllgen.Draw('GOFF,sames')
#
#cpr.prelim_tdr(lumi=2300,hor=0.17)
#tex.SetTextSize(0.03)
#tex.SetTextAlign(11) #left, bottom
##tex.DrawLatex(0.13,0.75,"pT Range [GeV]: %s"%(ptrange))
#
#c.Update()
#time.sleep(1)
#
#c.SaveAs("%s/gen%s.pdf"%(outdir,out_filename))
#c.SaveAs("%s/pdf/gen%s.pdf"%(webdir,out_filename))
#c.SaveAs("%s/png/gen%s.png"%(webdir,out_filename))
#print(  "%s/gen%s.png"%(www,out_filename))

