#ifndef plot_h
#define plot_h
 
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
 
 
class plot {
public :
   // Declaration of leaf types
   TString inpath;
   TString outpath;
   TString wwwpath;
   TString Tsubmitbase;
   TString Tversion;
   TString sysname;

   TFile* datafile;
   TFile* edatafile;
   TFile* mcfile;
   TFile* qcdfile;
   TFile* outfile;

   std::vector<std::vector<Double_t>> qcd_frac;
   std::vector<std::vector<Double_t>> qcd_frac_err;

   std::vector<std::vector<Double_t>> bincenterss;
   std::vector<std::vector<Double_t>> binerrorss;
   std::vector<std::vector<Double_t>> ratioss;
   std::vector<std::vector<Double_t>> ratioerrorss;

   std::vector<Double_t> ms;
   std::vector<Double_t> bs;
   std::vector<Double_t> mes;
   std::vector<Double_t> bes;

   ofstream log;

   std::vector<TString> sysnames;
   std::vector<int> ptbinbounds;
   int nptbins;
   std::vector<TString> ptbinnames;

   Double_t lower[5];

   plot();
   virtual ~plot();
   virtual void Loop();
   //My functions
   void  getFraction(TFile* datafile, 
                     TFile* edatafile,
                     TFile* mcfile,
                     TFile* qcdfile,
                     TFile* outfile,
                     std::vector<double>& fractionQCD,
                     std::vector<double>& fractionQCDErr,
                     std::vector<double>& fractionQCDSam,
                     std::vector<double>& fractionQCDErrSam,
                     TString sysname,
                     TString extraname
                    );

   void  getBinCenters(TFile* datafile, 
                     std::vector<Double_t>& bcenters,
                     std::vector<Double_t>& berrors,
                     TString sysname
                    );

   void  getCorrectedFakeRatio(TFile* datafile,
                               TFile* outfile,
                               std::vector<double> qcdfraction,
                               std::vector<double> qcderr,
                               std::vector<Double_t>& ratios,
                               std::vector<Double_t>& ratioerrors,
                               int sn,
                               TString extraname
                              );

   void getCorrectedFakeRatioErrUpDown(TFile* datafile,  //<----data file
                                 TFile* outfile,   //<---output file
                                 vector<double> corrFactor, //<--input corr. factor = QCD fraction
                                 vector<double> corrErr,    //<--input error
                                 vector<Double_t>& ratios,
                                 vector<Double_t>& ratioerrors,
                                 int sn,
                                 TString extraname
                                );

   void drawAllRates(TString extraname);
   void drawAllRatesRelative(TString extraname);
   void drawAllRatesRelativeHist(TString extraname);

};
#endif
 
#ifdef plot_cxx
plot::plot()
{            
}


plot::~plot()
{
}
#endif 

