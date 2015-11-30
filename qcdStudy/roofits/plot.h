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

 
#ifdef __MAKECINT__                           
#pragma link C++ class map<string,TCanvas*>;
#pragma link C++ class map<string,TPad*>;
#pragma link C++ class map<string,TLegend*>;
#pragma link C++ class map<string,TGraphErrors*>;
#pragma link C++ class map<string,TLine*>;
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
   TFile* f1;
   TFile* f2;
   TFile* f3;

   plot();
   virtual ~plot();
   virtual void Loop();
   //My functions
   void  getFraction(TFile* f1, 
                       TFile* f2,
                        TFile* f3,
                         std::vector<double>& fractionQCD,
                           std::vector<double>& fractionQCDErr,
                            std::vector<double>& fractionQCDSam,
                               std::vector<double>& fractionQCDErrSam
                             );

   void  getCorrectedFakeRatio(TFile* f1,
                                TFile* f3,
                                 std::vector<double> qcdfraction,
                                   std::vector<double> qcderr);

  bool biggerBins;
  bool smallerBins;
  bool normalBins;

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

