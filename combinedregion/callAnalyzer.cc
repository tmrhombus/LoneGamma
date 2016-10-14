
//#include "analyzeSignal.C"
#include "analyzeGenSignal.C"
//#include "analyzeZllG.C"
//#include "analyzeWlnG.C"

void callAnalyzer(void)
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

//  Bool_t isMC=kFALSE;
//  Bool_t isZnnG=kFALSE;
//  Bool_t ewkZG=kFALSE;
//  Bool_t ewkWG=kFALSE;
//  Bool_t isHalo=kFALSE;
//  Bool_t isSpike=kFALSE;
//  Bool_t isEle=kFALSE;
//  Bool_t isJet=kFALSE;
//  TString outfilename=path+"/ZnnG_Datafewer.root";
//  TString inputListName=path+"/SPD.fewer";
////    TString outfilename=path+"/ZnnG_Data_204.root";
////    TString inputListName=path+"/SinglePhotonData_callpostAnalyzer_Base-ggtree_data_204.inputs";
////    TString outfilename=path+"/ZnnG_Data_204_small.root";
////    TString inputListName=path+"/SinglePhotonData_callpostAnalyzer_Signal-ggtree_data_204.some";
////    TString outfilename=path+"/ZnnG_Data_2101.root";
////    TString inputListName=path+"/SinglePhotonData_callpostAnalyzer_Signal-ggtree_data_2101.inputs";
////    TString outfilename=path+"/ZnnG_Data_3401.root";
////    TString inputListName=path+"/SinglePhotonData_callpostAnalyzer_Signal-ggtree_data_3401.inputs";

  //TString inputListName=path+"/filenames_SinglePhoton2016_10.txt";


//  Bool_t isMC=kTRUE;
//  Bool_t isZnnG=kFALSE;
//  Bool_t ewkZG=kTRUE;
//  Bool_t ewkWG=kFALSE;
//  Bool_t isEle=kFALSE;
//  Bool_t isJet=kFALSE;
//  Bool_t isHalo=kFALSE;
//  Bool_t isSpike=kFALSE;
//  TString outfilename=path+"/ZllG_ZllG_v2.root";
//  TString inputListName=path+"/hdfslist_ZllGJets.txt";

//  Bool_t isMC=kTRUE;
//  Bool_t isZnnG=kFALSE;
//  Bool_t ewkZG=kFALSE;
//  Bool_t ewkWG=kTRUE;
//  Bool_t isEle=kFALSE;
//  Bool_t isJet=kFALSE;
//  Bool_t isHalo=kFALSE;
//  Bool_t isSpike=kFALSE;
//  TString outfilename=path+"/WlnG_WG_v2.root";
//  TString inputListName=path+"/hdfslist_WlnGJets.txt";

//  Bool_t isMC=kTRUE;
//  Bool_t isZnnG=kFALSE;
//  Bool_t ewkZG=kFALSE;
//  Bool_t ewkWG=kFALSE;
//  Bool_t isEle=kFALSE;
//  Bool_t isJet=kFALSE;
//  Bool_t isHalo=kFALSE;
//  Bool_t isSpike=kFALSE;
//  TString outfilename=path+"/WlnG_WZ_2.root";
//  TString inputListName=path+"/hdfslist_WZ.txt";

  Bool_t isMC=kTRUE;
  Bool_t isZnnG=kTRUE;
  Bool_t ewkZG=kTRUE;
  Bool_t ewkWG=kFALSE;
  Bool_t isEle=kFALSE;
  Bool_t isJet=kFALSE;
  Bool_t isHalo=kFALSE;
  Bool_t isSpike=kFALSE;
          lumi = 40000. ;
          nrEvents = 375920 ;
          crossSec = 0.1903 ;
  TString outfilename=path+"/Gen_ZnnG_v2.root";
  TString inputListName=path+"/hdfslist_ZnnGJets.txt";


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

//  analyzeSignal m;
//  m.Init(theChain,isMC);

//  analyzeZllG m;
//  m.Init(theChain,isMC);
//  m.InitLep();

//  analyzeWlnG m;
//  m.Init(theChain,isMC);
//  m.InitLep();

  analyzeGenSignal m;
  m.Init(theChain,isMC);
  m.InitGen();

  m.Loop(outfilename,isMC,lumi,nrEvents,crossSec,isZnnG,isEle,isHalo,isSpike,isJet,ewkWG,ewkZG);
}

