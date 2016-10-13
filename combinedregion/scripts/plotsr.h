#ifndef plotsr_h
#define plotsr_h
 
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
 
 
class plotsr {
public :

   TString inpath;
   TString outpath;
   TString wwwpath;
   TString Tsubmitbase;
   TString Tversion;
   TString sysname;

   std::vector<TString> regions ;

   TString outname;
   TString extraname;

   TString name_SinglePhotonData  ;
   TString name_GJets        ;
   TString name_TTGJets      ;
   TString name_TGJets       ;
   TString name_WlnGJets     ;
   TString name_ZllGJets     ;
   TString name_ZllJets      ;
   TString name_GGJets          ;
   TString name_SinglePhotonEle ;
   TString name_SinglePhotonJet ;
   TString name_SinglePhotonHalo ;
   TString name_SinglePhotonSpike ;
   TString name_WWG             ;
   TString name_WZ              ;
   TString name_Wmn             ;
   TString name_Wtn             ;
   TString name_ZZ              ;
   TString name_ZnnGJets        ;

   TFile* file_SinglePhotonData ;
   TFile* file_GJets        ;
   TFile* file_TTGJets      ;
   TFile* file_TGJets       ;
   TFile* file_WlnGJets     ;
   TFile* file_ZllGJets     ;
   TFile* file_ZllJets      ;
   TFile* file_GGJets          ;
   TFile* file_SinglePhotonEle ;
   TFile* file_SinglePhotonJet ;
   TFile* file_SinglePhotonHalo ;
   TFile* file_SinglePhotonSpike ;
   TFile* file_WWG             ;
   TFile* file_WZ              ;
   TFile* file_Wmn             ;
   TFile* file_Wtn             ;
   TFile* file_ZZ              ;
   TFile* file_ZnnGJets        ;
   TFile* file_out          ;

   ofstream log;

   std::vector<TString> variablenames;
   std::vector<TString> selectionnames;
   std::vector<TString> ptrangenames;
   std::vector<Int_t> rebins;
   std::vector<TString> leptons;

   bool dolog;
   bool dobygev;

   TH1F* raw_data_obs;
   TH1F* data_obs;

   TH1F* h_raw_GJets    ;   
   TH1F* h_raw_TTGJets  ;   
   TH1F* h_raw_TGJets   ;   
   TH1F* h_raw_ZllGJets ;   
   TH1F* h_raw_ZllJets  ;   
   TH1F* h_raw_GGJets   ;   
   TH1F* h_raw_EFake    ;   
   TH1F* h_raw_JFake    ;   
   TH1F* h_raw_Halo     ;   
   TH1F* h_raw_Spike    ;   
   TH1F* h_raw_WWG      ;    
   TH1F* h_raw_WZ       ;    
   TH1F* h_raw_WlnGJets ;
   TH1F* h_raw_Wmn      ;    
   TH1F* h_raw_Wtn      ;    
   TH1F* h_raw_ZZ       ;    
   TH1F* h_raw_ZnnGJets ;
   TH1F* h_raw_GWZtt    ;
   TH1F* h_raw_VV       ;

   TH1* h_GJets    ;   
   TH1* h_TTGJets  ;   
   TH1* h_TGJets   ;   
   TH1* h_ZllGJets ;   
   TH1* h_ZllJets  ;   
   TH1* h_GGJets   ;   
   TH1* h_EFake    ;   
   TH1* h_JFake    ;   
   TH1* h_Halo     ;   
   TH1* h_Spike    ;   
   TH1* h_WWG      ;    
   TH1* h_WZ       ;    
   TH1* h_WlnGJets ;
   TH1* h_Wmn      ;    
   TH1* h_Wtn      ;    
   TH1* h_ZZ       ;    
   TH1* h_ZnnGJets ;
   TH1* h_GWZtt    ;
   TH1* h_VV       ;

   TString histname;

   plotsr();
   virtual ~plotsr();
   virtual void Loop();
   //My functions
};
#endif
 
#ifdef plotsr_cxx
plotsr::plotsr()
{            
}


plotsr::~plotsr()
{
}
#endif 

