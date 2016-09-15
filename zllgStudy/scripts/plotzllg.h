#ifndef plotzllg_h
#define plotzllg_h
 
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include "TH1.h"
#include "TH2.h"
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TList.h>
#include <Riostream.h> 
#include <TGraphAsymmErrors.h>
#include <map>
#include "TRFIOFile.h"
#include "TMath.h"
#include <vector>
#include <TList.h>
#include <set>
#include "TPaveStats.h"
#include "TColor.h"
#include "TFractionFitter.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TStyle.h" 
#include "THStack.h" 


//ROOFIT headers
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooNDKeysPdf.h"
#include "TFile.h"
#include "TLegend.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <algorithm>
#include "RooGlobalFunc.h" 
#include "TLatex.h"

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
 
#ifdef __MAKECINT__                           
#pragma link C++ class map<string,TCanvas*>;
#pragma link C++ class map<string,TPad*>;
#pragma link C++ class map<string,TLegend*>;
#pragma link C++ class map<string,TGraphErrors*>;
#pragma link C++ class map<string,TLine*>;
#pragma link C++ class map<string,TLatex*>;
#pragma link C++ class vector<bool>+;
#pragma link C++ class vector<double>+;
#pragma link C++ class map<TString,TH1D*>+;
#endif
 
using namespace std;        
using namespace ROOT;
using namespace RooFit;
 
 
class plotzllg {
public :

   TString inpath;
   TString outpath;
   TString wwwpath;
   TString Tsubmitbase;
   TString Tversion;
   TString sysname;

   TString outname;
   TString extraname;

   TString name_SinglePhoton ;
   TString name_GJets        ;
   TString name_TTGJets      ;
   TString name_WlnGJets     ;
   TString name_ZllGJets     ;
   TString name_ZllJets      ;

   TFile* file_SinglePhoton ;
   TFile* file_GJets        ;
   TFile* file_TTGJets      ;
   TFile* file_WlnGJets     ;
   TFile* file_ZllGJets     ;
   TFile* file_ZllJets      ;
   TFile* file_out          ;

   ofstream log;

   std::vector<TString> variablenames;
   std::vector<TString> selectionnames;
   std::vector<TString> ptrangenames;
   std::vector<Int_t> rebins;

   bool dolog;

   TH1F* data_obs;
   TH1F* h_GJets   ;
   TH1F* h_TTGJets ;
   TH1F* h_WlnGJets;
   TH1F* h_ZllGJets;
   TH1F* h_ZllJets ;

   TString histname;

   plotzllg();
   virtual ~plotzllg();
   virtual void Loop();
   //My functions
};
#endif
 
#ifdef plotzllg_cxx
plotzllg::plotzllg()
{            
}


plotzllg::~plotzllg()
{
}
#endif 

