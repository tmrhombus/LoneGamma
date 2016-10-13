//#include "plotshapesr.C"
#include "plotsr.C"
#include "TROOT.h"
 
int main(){
 
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <map>");
  gSystem->Load("libDCache.so");
  //plotshapesr m;
  plotsr m;
  m.Loop();
 
}

