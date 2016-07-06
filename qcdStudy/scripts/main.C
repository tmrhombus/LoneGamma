#include "plotfrerrors.C"
//#include "plotfakeratio.C"
#include "TROOT.h"
 
int main(){
 
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <map>");
  gSystem->Load("libDCache.so");
  plotfrerrors m;
  //plotfakeratio m;
  m.Loop();
 
}

