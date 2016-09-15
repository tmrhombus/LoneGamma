
#include "postAnalyzer_ZnunuG.C"

void SAMPLENAME_callpostAnalyzer_ZnunuG(void)
{

 TChain *theChain = new TChain("TREENAME"); ;
 theChain->Reset();

 Bool_t isMC = ISMC ;
 Bool_t isZnnG = ISZNNG ;
 Bool_t ewkZG = EWKZG ;
 Double_t lumi = LUMI ;
 Double_t nrEvents = NREVENTS ;
 Double_t crossSec = CROSSSEC ;
 
 TString outfilename = getenv("OUTPUT");
 std::cout << "Output File Name: " << outfilename << std::endl;
 
 // get each filename from list and add it to TChain
 TString inputListName = getenv("INPUT");
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
 m.Init(theChain, isMC);
 m.Loop(outfilename,isMC,lumi,nrEvents,crossSec,isZnnG,ewkZG);
}

#if !defined(__CINT__) && !defined(__ACLIC__)
//
// "Standalone Application" entry point ("main").
//
// `root-config --cxx --cflags` -o callpostAnalyzer_ZnunuG callpostAnalyzer_ZnunuG.cxx `root-config --libs`
// ./callpostAnalyzer_ZnunuG
//
int main(int /*argc*/, char ** /*argv*/)
{
  SAMPLENAME_callpostAnalyzer_ZnunuG(); // just call the "ROOT Script"
  return 0;
}
#endif /* !defined(__CINT__) && !defined(__ACLIC__) */
