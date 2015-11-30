
#include "postAnalyzer_QCD.C"

//
// "ROOT Script" entry point (the same name as the "filename's base").
//
// [bash/csh] root callQCDAnalyzer.cxx
// [bash/csh] root callQCDAnalyzer.cxx++
// root [0] .x callQCDAnalyzer.cxx
// root [0] .x callQCDAnalyzer.cxx++
//
void callQCDAnalyzer(void)
{

 TChain *theChain = new TChain("ggNtuplizer/EventTree"); ;
 theChain->Reset();
 
// TString infilename="/hdfs/store/user/gomber/SinglePhoton_2015D_1p2fb1_condor/run_data_2015D_74X-005D4D5A-9871-E511-9942-02163E014303.root";
// theChain->Add(infilename);

 // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
 //GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM
 //hdfs/store/user/jjbuch/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/
 Double_t weight = 20730.; // +- 66 
// //GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM
// //hdfs/store/user/jjbuch/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/
// Double_t weight = 9226.; // +- 36 
// //GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM
// //hdfs/store/user/jjbuch/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/
// Double_t weight = 2300.; // +- 11
// //GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
// //hdfs/store/user/jjbuch/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/
// Double_t weight = 277.4; // +- 1.3
// //GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
// //hdfs/store/user/jjbuch/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/
// Double_t weight = 93.38; // +- 0.46
 
 Bool_t isMC=kFALSE;
 TString path = "./test";
 //TString outfilename=path+"/QCDs_MC.root";
 TString outfilename=path+"/QCDs_Data.root";
 std::cout << "Output File Name: " << outfilename << std::endl;
 
 // get each filename from list and add it to TChain
 //TString inputListName=path+"/filenames_MC_short.txt";
 TString inputListName=path+"/filenames_Data2015D_short.txt";
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

 postAnalyzer_QCD m;
 m.Init(theChain,isMC);
 m.Loop(outfilename,isMC,weight);
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
  callQCDAnalyzer(); // just call the "ROOT Script"
  return 0;
}
#endif /* !defined(__CINT__) && !defined(__ACLIC__) */
