
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
 
 Double_t lumi=12900.;
 Double_t nrEvents=12900.;
 Double_t crossSec=1.;
 
 TString path = "../test";

 //Bool_t isMC=kFALSE;
 //Bool_t isZnnG=kFALSE;
 //Bool_t isEle=kFALSE;
 //Bool_t isHalo=kFALSE;
 //Bool_t isSpike=kTRUE;
 //Bool_t isJet=kFALSE;
 //Bool_t ewkWG=kFALSE;
 //Bool_t ewkZG=kFALSE;

 //TString outfilename=path+"/Datafull.root";
 //TString inputListName=path+"/filenames_data.txt";
 //TString inputListName=path+"/hdfslist_SinglePhoton.txt";
 //TString inputListName=path+"/filenames_SinglePhoton2016_10.txt";

  Bool_t isMC=kFALSE;
  Bool_t isZnnG=kFALSE;
  Bool_t ewkZG=kFALSE;
  Bool_t ewkWG=kFALSE;
  Bool_t isHalo=kFALSE;
  Bool_t isSpike=kFALSE;
  Bool_t isEle=kFALSE;
  Bool_t isJet=kFALSE;
   ///TString outfilename=path+"/ZnnG_Data_204.root";
   ///TString inputListName=path+"/SinglePhotonData_callpostAnalyzer_Signal-ggtree_data_204.inputs";
   TString outfilename=path+"/ZnnG_Data_204_small.root";
   TString inputListName=path+"/SinglePhotonData_callpostAnalyzer_Signal-ggtree_data_204.some";
 // TString outfilename=path+"/ZnnG_Data_2101.root";
 // TString inputListName=path+"/SinglePhotonData_callpostAnalyzer_Signal-ggtree_data_2101.inputs";
//  TString outfilename=path+"/ZnnG_Data_3401.root";
//  TString inputListName=path+"/SinglePhotonData_callpostAnalyzer_Signal-ggtree_data_3401.inputs";

  //TString inputListName=path+"/filenames_SinglePhoton2016_10.txt";



// Bool_t isMC = kTRUE ;
// Bool_t isZnnG = kTRUE ;
// Bool_t isEle = kFALSE;
// Bool_t isHalo = kFALSE;
// Bool_t isSpike = kFALSE;
// Bool_t isJet = kFALSE;
// Bool_t ewkWG = kFALSE;
// Bool_t ewkZG = kTRUE;
//        lumi = 12900. ;
//        nrEvents = 375920 ;
//        crossSec = 0.1903 ;
//  TString outfilename=path+"/ZnnG_ZnnGJets.root";
//  TString inputListName=path+"/hdfslist_ZnnGJets.txt";


//  Bool_t isMC=kTRUE;
//  Bool_t isZnnG=kFALSE;
//  Bool_t ewkZG=kTRUE;
//  Bool_t ewkWG=kFALSE;
//  Bool_t isHalo=kFALSE;
//  Bool_t isSpike=kTRUE;
//  Bool_t isEle=kFALSE;
//  Bool_t isJet=kFALSE;
//  TString outfilename=path+"/ZnnG_GJets.root";
//  TString inputListName=path+"/hdfslist_GJetsHT200to400.txt";

// Bool_t isMC=kTRUE;
// Bool_t isZnnG=kTRUE;
// Bool_t isEle=kFALSE;
// Bool_t isHalo=kFALSE;
// Bool_t isSpike=kFALSE;
// Bool_t isJet=kFALSE;
// TString outfilename=path+"/Signal_ZnnGJets.root";
// TString inputListName=path+"/hdfslist_ZnnGJets.txt";

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
 m.Loop(outfilename,isMC,lumi,nrEvents,crossSec,isZnnG,isEle,isHalo,isSpike,isJet,ewkWG,ewkZG);
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
