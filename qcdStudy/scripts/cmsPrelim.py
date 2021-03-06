import ROOT
from ROOT import TLatex
tex = ROOT.TLatex()
tex.SetNDC(True)
def prelim(lumi=19700,tsize=0.07):
 tex.SetTextSize(tsize)
 tex.SetTextAlign(31) #right, bottom 
 tex.DrawLatex(0.9,0.9,'#font[32]{CMS Preliminary}')
 tex.SetTextSize(0.03)
 tex.SetTextAlign(33)
 lumiTitle = '#int L dt = %.1f fb^{-1}' %(lumi/1000.)
 tex.DrawLatex(0.87,0.87,lumiTitle)

def prelim_alt(lumi=19700,tsize=0.07):
 prlmTitle = 'CMS Preliminary'
 lumiTitle = '#int L dt = %.1f fb^{-1}' %(lumi/1000.)
 engyTitle = '#sqrt{s} = 8 TeV'
 tex.SetTextSize(tsize)
 tex.SetTextAlign(13) #left, top 
 tex.DrawLatex(0.15,0.89,'#font[32]{%s}'%(prlmTitle))
 tex.SetTextSize(0.03)
 tex.DrawLatex(0.13,0.85,lumiTitle)
 tex.DrawLatex(0.17,0.80,engyTitle)

def prelim_noLumi(tsize=0.07):
 tex.SetTextSize(tsize)
 tex.SetTextAlign(13) #left, top 
 tex.DrawLatex(0.11,0.89,'#font[32]{CMS Preliminary}')

def prelim_tdr(lumi=1200,hor=0.17):
 title = 'CMS'
 extra = 'Preliminary'
 lumi  = '%.2f fb^{-1} (13 TeV)' %(lumi/1000.)
 tex.SetTextSize(0.07)
 tex.SetTextAlign(13) #left, top 
 tex.DrawLatex(hor,0.89,'#font[61]{%s}'%(title))
 tex.SetTextSize(0.05)
 tex.DrawLatex(hor,0.83,'#font[52]{%s}'%(extra))

 tex.SetTextAlign(31)
 tex.SetTextSize(0.05)
 tex.DrawLatex(0.9,0.91,'#font[42]{%s}'%(lumi))
 
def public_tdr(lumi=19700,hor=0.17):
 title = 'CMS'
 lumi  = '%.1f fb^{-1} (8 TeV)' %(lumi/1000.)
 tex.SetTextSize(0.07)
 tex.SetTextAlign(13) #left, top 
 tex.DrawLatex(hor,0.89,'#font[61]{%s}'%(title))

 tex.SetTextAlign(31)
 tex.SetTextSize(0.05)
 tex.DrawLatex(0.9,0.91,'#font[42]{%s}'%(lumi))
