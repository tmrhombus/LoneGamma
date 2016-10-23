
#ifndef postAnalyzer_PDFscale_h
#define postAnalyzer_PDFscale_h

#include "postAnalyzer_Gen.h"


class postAnalyzer_PDFscale : public postAnalyzer_Gen {

 public:
   TH1F h_pdfs_uncorret[7],
        h_scale_uncorret[7];

   virtual void     InitPDFscale();
   Bool_t           FillSigHistogramsPDFscale(int ptbin, double weight);
   virtual void     callFillSigHistPDFscale(int selbin, int lastptbin, int inclptbin, int candphotonindex, float event_weight);
   Bool_t           WriteHistogramsPDFscale(int ptbin);

};

void postAnalyzer_PDFscale::InitPDFscale()
{

   std::cout<<"Initializing PDF/scale histograms"<<std::endl;

   ptbins.clear();
   ptbinnames.clear();

   ptbins.push_back(175);
   ptbins.push_back(190);
   ptbins.push_back(250);
   ptbins.push_back(400);
   ptbins.push_back(700);
   ptbins.push_back(1000);

   ptbinnames.push_back("175to190");
   ptbinnames.push_back("190to250");
   ptbinnames.push_back("250to400");
   ptbinnames.push_back("400to700");
   ptbinnames.push_back("700to1000");
   ptbinnames.push_back("175to1000");
   ptbinnames.push_back("allpt");

   for(unsigned int i=0; i<ptbinnames.size(); ++i){
//     // set up names
     TString histname_pdfs_uncorret  = "h_pdfs_uncorret_"+ptbinnames[i];
     TString histname_scale_uncorret = "h_scale_uncorret_"+ptbinnames[i];

     Float_t binspt[] = { 175., 195., 250., 400., 700., 1000. };
     Int_t  binnumpt = sizeof(binspt)/sizeof(Float_t) - 1; // or just = 9

     // reserve histograms
     h_pdfs_uncorret[i].Clear();
     h_pdfs_uncorret[i] = TH1F(histname_pdfs_uncorret,"PDFs",101,0.,101.);
     h_pdfs_uncorret[i].Sumw2();
     //

     h_scale_uncorret[i].Clear();
     h_scale_uncorret[i] = TH1F(histname_scale_uncorret,"#mu_{F},#mu_{R} scales",9,0.,9.);
     h_scale_uncorret[i].Sumw2();
     //

   }  // ptbinnames

}

Bool_t postAnalyzer_PDFscale::FillSigHistogramsPDFscale(int ptbin, double weight){

  for( int i=1; i<102; ++i ){
   h_pdfs_uncorret[ptbin].SetBinContent( i, ( h_pdfs_uncorret[ptbin].GetBinContent(i) + (weight*lheNormalizedWeights->at(110.+i)) ) ) ;
  }
  h_pdfs_uncorret[ptbin].SetBinContent( 0 , h_pdfs_uncorret[ptbin].GetBinContent(0) + weight ) ;

  h_scale_uncorret[ptbin].SetBinContent( 0 , h_scale_uncorret[ptbin].GetBinContent(0) + weight ) ;
  h_scale_uncorret[ptbin].SetBinContent( 1 , h_scale_uncorret[ptbin].GetBinContent(1) + weight*genWeight_QCDscale_muR1_muF1    ) ; 
  h_scale_uncorret[ptbin].SetBinContent( 2 , h_scale_uncorret[ptbin].GetBinContent(2) + weight*genWeight_QCDscale_muR1_muF2    ) ; 
  h_scale_uncorret[ptbin].SetBinContent( 3 , h_scale_uncorret[ptbin].GetBinContent(3) + weight*genWeight_QCDscale_muR1_muF0p5  ) ; 
  h_scale_uncorret[ptbin].SetBinContent( 4 , h_scale_uncorret[ptbin].GetBinContent(4) + weight*genWeight_QCDscale_muR2_muF1    ) ; 
  h_scale_uncorret[ptbin].SetBinContent( 5 , h_scale_uncorret[ptbin].GetBinContent(5) + weight*genWeight_QCDscale_muR2_muF2    ) ; 
  h_scale_uncorret[ptbin].SetBinContent( 6 , h_scale_uncorret[ptbin].GetBinContent(6) + weight*genWeight_QCDscale_muR2_muF0p5  ) ; 
  h_scale_uncorret[ptbin].SetBinContent( 7 , h_scale_uncorret[ptbin].GetBinContent(7) + weight*genWeight_QCDscale_muR0p5_muF1  ) ; 
  h_scale_uncorret[ptbin].SetBinContent( 8 , h_scale_uncorret[ptbin].GetBinContent(8) + weight*genWeight_QCDscale_muR0p5_muF2  ) ; 
  h_scale_uncorret[ptbin].SetBinContent( 9 , h_scale_uncorret[ptbin].GetBinContent(9) + weight*genWeight_QCDscale_muR0p5_muF0p5) ; 

 return kTRUE;
}

Bool_t postAnalyzer_PDFscale::WriteHistogramsPDFscale(int ptbin){

  h_pdfs_uncorret[ptbin].Write();
  h_scale_uncorret[ptbin].Write();

 return kTRUE;

}

//-------------------------callFillSigHistPDFscale
void postAnalyzer_PDFscale::callFillSigHistPDFscale(int selbin, int lastptbin, int inclptbin, int candphotonindex, float event_weight){
 Float_t uncorrectedPhoEt = ((*phoSCRawE)[candphotonindex]/TMath::CosH((*phoSCEta)[candphotonindex]));
 for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
  if(
     ( uncorrectedPhoEt > ptbins[ptb]) &&
     ( uncorrectedPhoEt < ptbins[ptb+1])
    ){
   FillSigHistogramsPDFscale(ptb, event_weight);
  } // end if passes pt cuts then fill
 } // end pt bin loop
 if(  // do an inclusive pT plot from bins
    ( uncorrectedPhoEt > ptbins[0]) &&
    ( uncorrectedPhoEt < ptbins[inclptbin])
   ){
  FillSigHistogramsPDFscale(lastptbin-1, event_weight);
 }
 // and one fully inclusive in pT
 FillSigHistogramsPDFscale(lastptbin, event_weight);
 return;
}


#endif