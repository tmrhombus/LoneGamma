
#include "postAnalyzer_Signal.C"

//
// "ROOT Script" entry point (the same name as the "filename's base").
//
// [bash/csh] root callSignalAnalyzer.cxx
// [bash/csh] root callSignalAnalyzer.cxx++
// root [0] .x callSignalAnalyzer.cxx
// root [0] .x callSignalAnalyzer.cxx++
//
void callSignalAnalyzer(void)
{

 TChain *theChain = new TChain("ggNtuplizer/EventTree"); ;
 theChain->Reset();
 
 Double_t lumi=2.24;
 Double_t nrEvents=10000;
 Double_t crossSec=40000;
 
 TString path = "../test";

// Bool_t isMC=kFALSE;
// Bool_t isZnnG=kFALSE;
// TString outfilename=path+"/Datatest.root";
// TString inputListName=path+"/SinglePhoton_callpostAnalyzer_Signal-ggtree_data_1117.inputs";

 Bool_t isMC=kTRUE;
 Bool_t isZnnG=kTRUE;
 TString outfilename=path+"/Signal_ZnnGJets.root";
 TString inputListName=path+"/hdfslist_ZnnGJets.txt";

 std::cout << "Input List Name:  " << inputListName << std::endl;
 std::cout << "Output File Name: " << outfilename << std::endl;

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

 postAnalyzer_Signal m;
 m.Init(theChain,isMC);
 m.Loop(outfilename,isMC,lumi,nrEvents,crossSec,isZnnG);
}

#if !defined(__CINT__) && !defined(__ACLIC__)
//
// "Standalone Application" entry point ("main").
//
// `root-config --cxx --cflags` -o callSignalAnalyzer callSignalAnalyzer.cxx `root-config --libs`
// ./callSignalAnalyzer
//
int main(int /*argc*/, char ** /*argv*/)
{
  callSignalAnalyzer(); // just call the "ROOT Script"
  return 0;
}
#endif /* !defined(__CINT__) && !defined(__ACLIC__) */
