
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
 
 Bool_t isMC=kFALSE;
 TString path = "../test";
 //TString outfilename=path+"/Signal_GJMC.root";
 //TString outfilename=path+"/Signal_ZGMC.root";
 TString outfilename=path+"/OneBadSignal_DataD.root";
 std::cout << "Output File Name: " << outfilename << std::endl;
 
 // get each filename from list and add it to TChain
 //TString inputListName=path+"/filenames_GJMC_10.txt";
 //TString inputListName=path+"/filenames_ZGMC_10.txt";
 //TString inputListName=path+"/filenames_Data2015D_Full.txt";
 //TString inputListName=path+"/filenames_Data2015D_20.txt";
 //TString inputListName=path+"/SinglePhoton_badevent.inputs";
 //TString inputListName=path+"/SinglePhoton_callpostAnalyzer_Signal-ggtree_data_1017.inputs";
 TString inputListName=path+"/SinglePhoton_callpostAnalyzer_Signal-ggtree_data_1117.inputs";
 //TString inputListName=path+"/SinglePhoton_callpostAnalyzer_Signal-ggtree_data_37.inputs";
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
 m.Loop(outfilename,isMC,lumi,nrEvents,crossSec);
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
