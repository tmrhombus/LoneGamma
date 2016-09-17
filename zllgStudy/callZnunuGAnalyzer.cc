
#include "postAnalyzer_ZnunuG.C"

//
// "ROOT Script" entry point (the same name as the "filename's base").
//
// [bash/csh] root callZnunuGAnalyzer.cxx
// [bash/csh] root callZnunuGAnalyzer.cxx++
// root [0] .x callZnunuGAnalyzer.cxx
// root [0] .x callZnunuGAnalyzer.cxx++
//
void callZnunuGAnalyzer(void)
{

 TChain *theChain = new TChain("ggNtuplizer/EventTree"); ;
 theChain->Reset();
 
 Double_t lumi=2.6;
 Double_t nrEvents=10000;
 Double_t crossSec=40000;
 
 TString path = "../test";

//  Bool_t isMC=kFALSE;
//  Bool_t isZnnG=kFALSE;
//  Bool_t ewkZG=kFALSE;
//  TString inputListName=path+"/filenames_SinglePhoton2016_10.txt";
//  TString outfilename=path+"/ZllG_Data.root";

  Bool_t isMC=kTRUE;
  Bool_t isZnnG=kFALSE;
  Bool_t ewkZG=kTRUE;
  Bool_t isEle=kTRUE;
  Bool_t isJet=kTRUE;
  TString outfilename=path+"/ZllG_WZ.root";
  TString inputListName=path+"/hdfslist_WZ.txt";

// Bool_t isMC=kTRUE;
// Bool_t isZnnG=kTRUE;
// TString outfilename=path+"/ZNNG_ZGMC.root";
// TString inputListName=path+"/hdfslist_ZnnGJets.txt";

 std::cout << "Output File Name: " << outfilename << std::endl;
 
 // get each filename from list and add it to TChain
 std::vector<TString> infilename_dump;
  
 // open file_name_list.txt
 ifstream inputList;
 inputList.open(inputListName);
 if( !inputList.good() ) { 
   std::cerr << "Cannot open the file: \"" << inputListName+"\""<<std::endl;
   abort();
 }
 // we have the file open, start reading lines one by one
 TString infilename = ""; 
 while( !inputList.eof() ) { 
  infilename="";
  inputList >> infilename;
  if( inputList.fail() ) continue; 
  
  std::cout << "Input File Name: "  << infilename <<  std::endl;
 
  theChain->Add( infilename );
 
  infilename_dump.push_back(infilename);
 } //while !inputList.eof()

 postAnalyzer_ZnunuG m;
 m.Init(theChain,isMC);
 m.Loop(outfilename,isMC,lumi,nrEvents,crossSec,isZnnG,ewkZG,isEle,isJet);
}

#if !defined(__CINT__) && !defined(__ACLIC__)
//
// "Standalone Application" entry point ("main").
//
// `root-config --cxx --cflags` -o callZnunuGAnalyzer callZnunuGAnalyzer.cxx `root-config --libs`
// ./callZnunuGAnalyzer
//
int main(int /*argc*/, char ** /*argv*/)
{
  callZnunuGAnalyzer(); // just call the "ROOT Script"
  return 0;
}
#endif /* !defined(__CINT__) && !defined(__ACLIC__) */
