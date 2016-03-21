
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
 
 Double_t lumi=2.24;
 Double_t nrEvents=10000;
 Double_t crossSec=40000;
 
 Bool_t isMC=kTRUE;
 TString path = "./test";
 TString outfilename=path+"/ZNNG_ZGMC.root";
 //TString outfilename=path+"/ZNNG_DataD.root";
 std::cout << "Output File Name: " << outfilename << std::endl;
 
 // get each filename from list and add it to TChain
 TString inputListName=path+"/filenames_ZGMC_10.txt";
 //TString inputListName=path+"/filenames_Data2015D_20.txt";
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
 m.Loop(outfilename,isMC,lumi,nrEvents,crossSec);
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
