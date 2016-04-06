#!/usr/bin/env python
'''
Plots Numerator distributions prefit
Author: T.M.Perry UW
'''
import ROOT
from ROOT import THStack,TH1F,TFile
from ROOT import TLegend,TCanvas,TPad,TLatex,TLine
from ROOT import gROOT,gStyle
import sys
import time
import cmsPrelim as cpr

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v","--version",            help="version name")
parser.add_argument("-id","--inpdir",            help="/full/path/to/input/directory")
parser.add_argument("-mn","--mc_filename",       help="mc filename (.root)")
parser.add_argument("-qn","--qcd_filename",      help="qcd filename (.root)")
parser.add_argument("-dn","--data_filename",     help="data filename (.root)")
parser.add_argument("-en","--data_ele_filename", help="electron data filename (.root)")
parser.add_argument("-od","--outdir",            help="/path/to/output/directory")
parser.add_argument("-on","--out_filename",      help="output name (no .root)")
parser.add_argument("-var","--variable",         help="variable name")
parser.add_argument("-rng","--ptrange",          help="pT range (ex. 250to400")
parser.add_argument("-log","--do_log",  action='store_true', help="Set Log Scale")
args = parser.parse_args()

version            = args.version
inpdir             = args.inpdir        
mc_filename        = args.mc_filename   
qcd_filename       = args.qcd_filename   
data_filename      = args.data_filename   
data_ele_filename  = args.data_ele_filename   
outdir             = args.outdir  
out_filename       = args.out_filename  
variable           = args.variable
ptrange            = args.ptrange
do_log             = args.do_log
webdir = "/afs/hep.wisc.edu/user/tperry/www/MonoPhoton/qcdPlots/%s"%(version)

import errno
import os

try:
 os.mkdir(webdir)
except OSError as exc:
 if exc.errno != errno.EEXIST:
  raise exc
 pass


xn = ""
if(do_log): out_filename+="_log"

m_file = TFile("%s/%s"%(inpdir,mc_filename))
q_file = TFile("%s/%s"%(inpdir,qcd_filename))
d_file = TFile("%s/%s"%(inpdir,data_filename))
e_file = TFile("%s/%s"%(inpdir,data_ele_filename))

#canvas attributes
canx = 1200
cany = 1100 
gStyle.SetOptStat('')
gStyle.SetLineWidth(3)
gStyle.SetPadTickY(1)

#color scheme
c_ntot = ROOT.kBlack
c_etot = 14 #ROOT.kBlack-2
c_nsig = ROOT.kRed
c_nbkg = ROOT.kBlue
c_deno = ROOT.kGreen+1
c_qden = ROOT.kCyan+2
fs = 0
# line style
n_ls = 1
d_ls = 9
q_ls = 7

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

# Numerator Total (Data, Signal Selection)
h_ntot = d_file.Get("h_sig_%s_%s"%(variable,ptrange))
h_ntot.SetName("h_ntot")
h_ntot.SetTitle("")
h_ntot.SetLineColor(c_ntot)
h_ntot.SetLineWidth(4)
h_ntot.SetFillStyle(fs)
h_ntot.SetLineStyle(n_ls)
h_ntot.Draw("GOFF")
h_ntot.Scale( 1. / max(h_ntot.Integral(-1,-1),1.) )
themax = h_ntot.GetMaximum()

# Numerator (electron shape)  Total (Data Ele, Signal Selection)
h_etot = e_file.Get("h_sig_%s_%s"%(variable,ptrange))
h_etot.SetName("h_etot")
h_etot.SetTitle("")
h_etot.SetLineColor(c_etot)
h_etot.SetLineWidth(3)
h_etot.SetFillStyle(fs)
h_etot.SetLineStyle(n_ls)
h_etot.Draw("GOFF")
h_etot.Scale( 1. / max(h_etot.Integral(-1,-1),1.) )
themax = max( themax, h_etot.GetMaximum() )

# Numerator Background (Data, 5 < chIso < 10 sideband)
h_nbkg = d_file.Get("h_bkg_%s_%s"%(variable,ptrange))
h_nbkg.SetName("h_nbkg")
h_nbkg.SetTitle("")
h_nbkg.SetLineColor(c_nbkg)
h_nbkg.SetLineWidth(2)
h_nbkg.SetFillStyle(fs)
h_nbkg.SetLineStyle(n_ls)
h_nbkg.Draw("GOFF")
h_nbkg.Scale( 0.5 / max(h_nbkg.Integral(-1,-1),1.) )
themax = max( themax, h_nbkg.GetMaximum() )

# Numerator Signal (MC, Signal Selection)
h_nsig = m_file.Get("h_sig_%s_%s"%(variable,ptrange))
h_nsig.SetName("h_nsig")
h_nsig.SetTitle("")
h_nsig.SetLineColor(c_nsig)
h_nsig.SetLineWidth(2)
h_nsig.SetFillStyle(fs)
h_nsig.SetLineStyle(n_ls)
h_nsig.Draw("GOFF")
h_nsig.Scale( 0.5 / max(h_nsig.Integral(-1,-1),1.) )
themax = max( themax, h_nsig.GetMaximum() )

# Denominator (Data, Inv. Iso / Pass V.Loose ID/ Fail Loose ID)
h_deno = d_file.Get("h_den_%s_%s"%(variable,ptrange))
h_deno.SetName("h_deno")
h_deno.SetTitle("")
h_deno.SetLineColor(c_deno)
h_deno.SetLineWidth(2)
h_deno.SetFillStyle(fs)
h_deno.SetLineStyle(d_ls)
h_deno.Draw("GOFF")
h_deno.Scale( 0.5 / max(h_deno.Integral(-1,-1),1.) )
themax = max( themax, h_deno.GetMaximum() )

# QCD Denominator (QCD MC, Inv. Iso / Pass V.Loose ID/ Fail Loose ID)
h_qden = q_file.Get("h_den_%s_%s"%(variable,ptrange))
h_qden.SetName("h_qden")
h_qden.SetTitle("")
h_qden.SetLineColor(c_qden)
h_qden.SetLineWidth(2)
h_qden.SetFillStyle(fs)
h_qden.SetLineStyle(q_ls)
h_qden.Draw("GOFF")
h_qden.Scale( 0.5 / max(h_qden.Integral(-1,-1),1.) )
themax = max( themax, h_qden.GetMaximum() )

# fill legend
leg=TLegend(0.5,0.55,0.88,0.88)
leg.AddEntry(h_ntot,"Numerator (photon data)")
leg.AddEntry(h_etot,"Numerator (electron data)")
leg.AddEntry(h_nsig,"Numerator Sig. (GJ MC)")
leg.AddEntry(h_nbkg,"Numerator Bkg. (data)")
leg.AddEntry(h_deno,"Denominator (data)")
leg.AddEntry(h_qden,"Denominator (QCD MC)")
leg.SetFillColor(0)
leg.SetBorderSize(0)

# and draw
h_ntot.SetMaximum( 1.2*themax )

h_ntot.Draw("GOFF")
h_ntot.GetYaxis().SetTitle("Relative Events")
h_ntot.GetXaxis().SetTitle("#sigma(i#eta,i#eta) Full 5x5")
#h_ntot.GetXaxis().SetTitle("%s"%(variable))

#if(do_log):  h_ntot.SetMaximum( 50*themax ) 
if(do_log):  h_ntot.SetMaximum( 11 ) 
c.cd()
h_ntot.Draw('hist,GOFF')
h_etot.Draw('hist,GOFF,sames')
h_nbkg.Draw('hist,GOFF,sames')
h_nsig.Draw('hist,GOFF,sames')
h_qden.Draw('hist,GOFF,sames')
h_deno.Draw('hist,GOFF,sames')
leg.Draw('sames,GOFF')

cpr.prelim_tdr(lumi=2320,hor=0.17)
tex.SetTextSize(0.03)
tex.SetTextAlign(11) #left, bottom
tex.DrawLatex(0.13,0.75,"pT Range [GeV]: %s"%(ptrange))

c.Update()
time.sleep(1)

#c.Print("%s/%s.png"%(outdir,out_filename))
#c.SaveAs("%s/%s.eps"%(outdir,out_filename))
#c.SaveAs("%s/%s.C"%(outdir,out_filename))
c.SaveAs("%s/%s.pdf"%(outdir,out_filename))
c.SaveAs("%s/%s.pdf"%(webdir,out_filename))
print(  "%s/%s.pdf"%(outdir,out_filename))

