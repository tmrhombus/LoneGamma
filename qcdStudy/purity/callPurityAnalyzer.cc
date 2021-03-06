
#include "postAnalyzerMC_purity.C"

//
// "ROOT Script" entry point (the same name as the "filename's base").
//
// [bash/csh] root callQCDAnalyzer.cxx
// [bash/csh] root callQCDAnalyzer.cxx++
// root [0] .x callQCDAnalyzer.cxx
// root [0] .x callQCDAnalyzer.cxx++
//
void callPurityAnalyzer(void)
{

 TChain *theChain = new TChain("ggNtuplizer/EventTree"); ;
 theChain->Reset();
 
 Double_t lumi=2.32;
 Double_t nrEvents=10000;
 Double_t crossSec=40000;
 
 Bool_t isMC=kTRUE;
 TString path = "../../test";

// TString outfilename=path+"/Purity_synchfull_GJ.root";
// TString inputListName="/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/gitignore/Madrid/lists/hdfslist_GJets_HT200To400.txt";

// TString outfilename=path+"/Purity_synch_QCD.root";
// TString inputListName=path+"/filenames_synch_QCD.txt";

 TString outfilename=path+"/Purity_GJ.root";
 TString inputListName=path+"/filenames_GJMC_10.txt";

 //TString outfilename=path+"/Purity_synch_GJ.root";
 //TString inputListName=path+"/filenames_synch_GJ.txt";
 //TString inputListName=path+"/GJets_HT100To200_callpostAnalyzerMC_purity-ggtree_mc_280.inputs";

 //TString outfilename=path+"/QCDs_QCDMC.root";
 //TString inputListName=path+"/filenames_MC_QCD_5.txt";

 //TString outfilename=path+"/Purity_GJMC.root";
 //TString outfilename=path+"/Purity_EleData.root";
 std::cout << "Output File Name: " << outfilename << std::endl;
 
 // get each filename from list and add it to TChain
 //TString inputListName=path+"/filenames_GJMC_10.txt";
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

 postAnalyzerMC_purity m;
 m.Init(theChain);
 m.Loop(outfilename,isMC,lumi,nrEvents,crossSec);
}

#if !defined(__CINT__) && !defined(__ACLIC__)
//
// "Standalone Application" entry point ("main").
//
// `root-config --cxx --cflags` -o callQCDAnalyzer callQCDAnalyzer.cxx `root-config --libs`
// ./callQCDAnalyzer
//
int main(int /*argc*/, char ** /*argv*/)
{
  callPurityAnalyzer(); // just call the "ROOT Script"
  return 0;
}
#endif /* !defined(__CINT__) && !defined(__ACLIC__) */
