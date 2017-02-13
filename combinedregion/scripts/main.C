//#include "plotshapesr.C"
//#include "plottf.C"
#include "plotuncs.C"
//#include "plotsr.C"
#include "TROOT.h"
 
int main(){
 
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <map>");
  gSystem->Load("libDCache.so");
  //plotshapesr m;
  //plottf m;
  plotuncs m;
  //plotsr m;
  m.Loop();
 
}

