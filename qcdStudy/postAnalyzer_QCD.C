#define postAnalyzer_QCD_cxx
#include "postAnalyzer_QCD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void postAnalyzer_QCD::Loop(TString outfilename, Bool_t isMC, Bool_t isEle, Double_t lumi, Double_t nrEvents, Double_t crossSec)
{

 TStopwatch sw; 
 sw.Start();

 Int_t nc = 0;
 Int_t dc = 0;

 if (fChain == 0) return;
 //std::cout<<"Numerator/Denominator : run : lumis : event : photonindex "<<std::endl;

 Long64_t nentries = fChain->GetEntriesFast();

 Long64_t nbytes = 0, nb = 0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
 //for (Long64_t jentry=0; jentry<1000;jentry++) {
  if (jentry%20000 == 0)
    {
      std::cout<<"Starting entry "<<jentry<<"/"<<(nentries)<<" at "<<sw.RealTime()<<" RealTime, "<<sw.CpuTime()<<" CpuTime"<<std::endl;
      sw.Continue();
    }
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   nbytes += nb;

  //=1.0 for real data
  event_weight=1.0;
  if(isMC){ event_weight=lumi*crossSec/nrEvents; }
  
   int inclptbin = ptbins.size()-1;
   int lastptbin = ptbinnames.size()-1;
   int lastsysbin = sysbinnames.size();

   // vector of ints, each int corresponds to location in vector of photons of photon passing cut
   // one such vector for each systematic bin
   for(unsigned int sysb=0; sysb<lastsysbin; ++sysb){
    sigPCvint[sysb] = pcPassSel(0,sysb,175,5000,1.4442,isMC,isEle); // passes signal selection (no sieie cut)
    bkgPCvint[sysb] = pcPassSel(1,sysb,175,5000,1.4442,isMC,isEle); // passes background selection (no sieie cut)
    denPCvint[sysb] = pcPassSel(2,sysb,175,5000,1.4442,isMC,isEle); // passes denominator selection (no sieie cut)
   }

   // MET requirements
   bool passMETfilters =  ( metFilters==0 || metFilters==4 ) ;

   // some more filters
   bool passTriggers;
   bool passTriggers175;
   bool passTriggers250;
   bool passMET;

// fill histograms
   for(unsigned int sysb=0; sysb<lastsysbin; ++sysb){
    passMET = pfMET < 30. && passMETfilters ;
    if(sysb==3){ passMET = pfMET < 45. && passMETfilters ;}
    if(sysb==4){ passMET = pfMET < 15. && passMETfilters ;}
    //passMET = pfMET < 75. && passMETfilters ;
    //if(sysb==3){ passMET = pfMET < 90. && passMETfilters ;}
    //if(sysb==4){ passMET = pfMET < 60. && passMETfilters ;}

    // Fill Numerator Signal Histograms ------------------------------------------
    if( sigPCvint[sysb].size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<sigPCvint[sysb].size(); ++k){ // go through each photon passing sel
      int phonr = sigPCvint[sysb].at(k);
      passTriggers = ((HLTPho>>12&1) == 1   && ( (phoFiredSingleTrgs->at(phonr)>>7)!=1 ) ) ;
      passTriggers175 = ((HLTPho>>7&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>8)!=1 ) ) ;
      passTriggers250 = ((HLTPho>>8&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>9)!=1 ) ) ;
      for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
       if(
          (phoEt->at(phonr) > ptbins[ptb]) &&
          (phoEt->at(phonr) < ptbins[ptb+1])
         ){
           if( sysb==0 && ptb==0 ) {
            nc++;
           }
        FillSigHistograms(ptb, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillSigHistograms(ptb, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillSigHistograms(ptb, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillSigHistograms(ptb, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillSigHistograms(ptb, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillSigHistograms(ptb, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillSigHistograms(ptb, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillSigHistograms(ptb, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillSigHistograms(ptb, sysb, 8, phonr, event_weight); }
       }
      }
      if(  // do an inclusive pT plot from bins
         (phoEt->at(phonr) > ptbins[0]) &&
         (phoEt->at(phonr) < ptbins[inclptbin])
        ){
       FillSigHistograms(lastptbin-1, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillSigHistograms(lastptbin-1, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillSigHistograms(lastptbin-1, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillSigHistograms(lastptbin-1, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillSigHistograms(lastptbin-1, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillSigHistograms(lastptbin-1, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillSigHistograms(lastptbin-1, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillSigHistograms(lastptbin-1, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillSigHistograms(lastptbin-1, sysb, 8, phonr, event_weight); }
      }
      // and one fully inclusive in pT
      FillSigHistograms(lastptbin, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillSigHistograms(lastptbin, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillSigHistograms(lastptbin, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillSigHistograms(lastptbin, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillSigHistograms(lastptbin, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillSigHistograms(lastptbin, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillSigHistograms(lastptbin, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillSigHistograms(lastptbin, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillSigHistograms(lastptbin, sysb, 8, phonr, event_weight); }
     }
    }

    // Fill Numerator Background Histograms -----------------------------------------
    if( bkgPCvint[sysb].size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<bkgPCvint[sysb].size(); ++k){ // go through each photon passing sel
      int phonr = bkgPCvint[sysb].at(k);
      passTriggers = ((HLTPho>>12&1) == 1   && ( (phoFiredSingleTrgs->at(phonr)>>7)!=1 ) ) ;
      passTriggers175 = ((HLTPho>>7&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>8)!=1 ) ) ;
      passTriggers250 = ((HLTPho>>8&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>9)!=1 ) ) ;
      for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
       if(
          (phoEt->at(phonr) > ptbins[ptb]) &&
          (phoEt->at(phonr) < ptbins[ptb+1])
         ){
        FillBkgHistograms(ptb, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillBkgHistograms(ptb, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillBkgHistograms(ptb, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillBkgHistograms(ptb, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillBkgHistograms(ptb, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillBkgHistograms(ptb, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillBkgHistograms(ptb, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillBkgHistograms(ptb, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillBkgHistograms(ptb, sysb, 8, phonr, event_weight); }
       }
      }
      if(  // do an inclusive pT plot from bins
         (phoEt->at(phonr) > ptbins[0]) &&
         (phoEt->at(phonr) < ptbins[inclptbin])
        ){
       FillBkgHistograms(lastptbin-1, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillBkgHistograms(lastptbin-1, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillBkgHistograms(lastptbin-1, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillBkgHistograms(lastptbin-1, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillBkgHistograms(lastptbin-1, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillBkgHistograms(lastptbin-1, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillBkgHistograms(lastptbin-1, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillBkgHistograms(lastptbin-1, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillBkgHistograms(lastptbin-1, sysb, 8, phonr, event_weight); }
      }
      // and one fully inclusive in pT
      FillBkgHistograms(lastptbin, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillBkgHistograms(lastptbin, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillBkgHistograms(lastptbin, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillBkgHistograms(lastptbin, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillBkgHistograms(lastptbin, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillBkgHistograms(lastptbin, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillBkgHistograms(lastptbin, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillBkgHistograms(lastptbin, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillBkgHistograms(lastptbin, sysb, 8, phonr, event_weight); }
     }
    }

    // Fill Denominator Histograms
    if( denPCvint[sysb].size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<denPCvint[sysb].size(); ++k){ // go through each photon passing sel
      int phonr = denPCvint[sysb].at(k);
      passTriggers = ((HLTPho>>12&1) == 1   && ( (phoFiredSingleTrgs->at(phonr)>>7)!=1 ) ) ;
      passTriggers175 = ((HLTPho>>7&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>8)!=1 ) ) ;
      passTriggers250 = ((HLTPho>>8&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>9)!=1 ) ) ;
      for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
       if(
          (phoEt->at(phonr) > ptbins[ptb]) &&
          (phoEt->at(phonr) < ptbins[ptb+1])
         ){
           if( sysb==0 && ptb==0 ) {
            dc++;
           }
        FillDenHistograms(ptb, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillDenHistograms(ptb, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillDenHistograms(ptb, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillDenHistograms(ptb, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillDenHistograms(ptb, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillDenHistograms(ptb, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillDenHistograms(ptb, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillDenHistograms(ptb, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillDenHistograms(ptb, sysb, 8, phonr, event_weight); }
       }
      }
      if(  // do an inclusive pT plot from bins
         (phoEt->at(phonr) > ptbins[0]) &&
         (phoEt->at(phonr) < ptbins[inclptbin])
        ){
       FillDenHistograms(lastptbin-1, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillDenHistograms(lastptbin-1, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillDenHistograms(lastptbin-1, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillDenHistograms(lastptbin-1, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillDenHistograms(lastptbin-1, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillDenHistograms(lastptbin-1, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillDenHistograms(lastptbin-1, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillDenHistograms(lastptbin-1, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillDenHistograms(lastptbin-1, sysb, 8, phonr, event_weight); }
      }
      // and one fully inclusive in pT
      FillDenHistograms(lastptbin, sysb, 0, phonr, event_weight);
if(passMET){                                                       FillDenHistograms(lastptbin, sysb, 1, phonr, event_weight); }
if(passTriggers){                                                  FillDenHistograms(lastptbin, sysb, 2, phonr, event_weight); }
if(passTriggers && passMET){                                       FillDenHistograms(lastptbin, sysb, 3, phonr, event_weight); }
if(passTriggers175){                                               FillDenHistograms(lastptbin, sysb, 4, phonr, event_weight); }
if(passTriggers175 && passMET){                                    FillDenHistograms(lastptbin, sysb, 5, phonr, event_weight); }
if(passTriggers250){                                               FillDenHistograms(lastptbin, sysb, 6, phonr, event_weight); }
if(passTriggers250 && passMET){                                    FillDenHistograms(lastptbin, sysb, 7, phonr, event_weight); }
if(passTriggers && passTriggers175 && passTriggers250 && passMET){ FillDenHistograms(lastptbin, sysb, 8, phonr, event_weight); }
     }
    }

   } // for each sysb in sysbins
// end fill histograms

 } //end loop through entries

 // write these histograms to file
  std::cout<<std::endl;
  std::cout<<"Total Passing Numerator: "<<nc<<"  Total Passing Denominator: "<<dc<<std::endl;
  std::cout<<"made it through, about to write"<<std::endl;

 TFile *outfile = new TFile(outfilename,"RECREATE");
 outfile->cd();
 for(unsigned int i=0; i<ptbinnames.size(); ++i){
  for(unsigned int j=0; j<sysbinnames.size(); ++j){
   for(unsigned int k=0; k<cuts.size(); ++k){
    WriteHistograms(i,j,k);
   }
  }
 }
 outfile->Close();
 sw.Stop();
 std::cout<<"Real Time: "<<sw.RealTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"CPU Time: "<<sw.CpuTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"Done"<<std::endl;

} //end Loop()

//   pcPassSel photons passing selections: ( num_sig(0), num_bkg(1), den(2) )
 // systnames   "" "_sbUP" "_sbDown" "_metUP" "_metDown" "_binUP" "_binDown" "_noPiso"
std::vector<int> postAnalyzer_QCD::pcPassSel(int sel, int sys, double phoPtLo, double phoPtHi, double phoEtaMax, Bool_t isMC, Bool_t isEle){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool noncoll;
  bool passKinematics;
  bool passMET;

  bool passHoE;
  bool passPSeed;
  bool passPhoNHMedIso;
  bool passCHMedIso; 
  double chIsoLB, chIsoUB;
  bool passCHBkgIso; 
  double vloosePFCharged;
  double vloosePFPhoton ;
  double vloosePFNeutral;
  bool passLooseIso;
  bool passLoosePIso;
  bool passVLooseIso;

  bool passPhotonID;
  //Loop over photons
  for(int p=0;p<nPho;p++)  // nPho from ntuple (should = phoE->length() )
    {
     // from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#SPRING15_selections_25_ns

     // non collision backgrounds
     //noncoll = kTRUE;
     if(isMC){
      noncoll = (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001; // isMC
     }
     if(!isMC){
      noncoll = fabs((*phoseedTimeFull5x5)[p]) < 3. && (*phomipTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001;
     }

     passKinematics = (
                       ( (*phoEt)[p] > phoPtLo  ) &&
                       ( (*phoEt)[p] <= phoPtHi ) &&
                       ( fabs((*phoSCEta)[p]) < phoEtaMax)
                      );

     // nsig, nbkg, deno
     passHoE = ( (*phoHoverE)[p] < 0.05 );
     passPSeed = false;
     if(isEle){ passPSeed = (*phohasPixelSeed)[p] ==  1; } 
     else     { passPSeed = (*phohasPixelSeed)[p] ==  0; } 
     // nsig, nbkg
     passPhoNHMedIso = (
                        ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < 
                           1.06 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))
                        )  &&
                        ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < 
                           0.28 + (0.0053 * (*phoEt)[p])
                        )
                       );

     // nsig
     passCHMedIso = ( TMath::Max( ( (*phoPFChWorstIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.37 );

     // nbkg
     //if( (*phoEt)[p] <= 190 ){
     // chIsoLB = 9.;
     // chIsoUB = 20.;
     //}
     //if( (*phoEt)[p] > 190 && (*phoEt)[p] <= 250 ){
     // chIsoLB = 6.;
     // chIsoUB = 15.;
     //}
     //if( (*phoEt)[p] > 250 ){
     // chIsoLB = 6.;
     // chIsoUB = 13.;
     //}
     //if(sys==1){chIsoUB = chIsoUB+2.;}
     //if(sys==2){chIsoUB = chIsoLB-2.;}

     chIsoLB = 9.;
     chIsoUB = 24.;
     if(sys==1){chIsoUB = 26.;}
     if(sys==2){chIsoUB = 22.;}

     //passCHBkgIso = (
     //                ( TMath::Max( ( (*phoPFChWorstIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) > chIsoLB )  &&
     //                ( TMath::Max( ( (*phoPFChWorstIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < chIsoUB ) 
     //               );
     passCHBkgIso = (
                     ( TMath::Max( (Float_t)(*phoPFChIso)[p], (Float_t)0.0) > chIsoLB )  &&
                     ( TMath::Max( (Float_t)(*phoPFChIso)[p], (Float_t)0.0) < chIsoUB ) 
                    );

     // deno
     vloosePFCharged= TMath::Min(5.0*(3.32) , 0.20*(*phoEt)[p]);
     vloosePFPhoton = TMath::Min(5.0*(0.81+ (0.0053 * (*phoEt)[p])) , 0.20*(*phoEt)[p]);
     vloosePFNeutral= TMath::Min(5.0*(1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))) , 0.20*(*phoEt)[p]);
     passVLooseIso = ( 
                      ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < vloosePFCharged )  &&  
                      ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < vloosePFNeutral )  &&  
                      ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < vloosePFPhoton )
                     ); 
     passLoosePIso = 
                     ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) <
                       (0.81 + (0.0053 * (*phoEt)[p])) );
     if(sys==7){passLoosePIso= true;}
     passLooseIso = (  // deno must fail this cut
                     ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 3.32 )  &&
                     ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) <
                       (1.92 + (0.014* (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))  &&
                     passLoosePIso
                    );

     if(sel==0){  // numerator signal
      passPhotonID = passHoE && passPSeed && passPhoNHMedIso && passCHMedIso;
     }
     else if(sel==1){ // numerator background
      passPhotonID = passHoE && passPSeed && passPhoNHMedIso && passCHBkgIso;
     }
     else if(sel==2){ // denominator
      passPhotonID = passHoE && passPSeed && !passLooseIso && passVLooseIso;
     }
     if(noncoll && passPhotonID && passKinematics){
      // std::cout<<" Found a photon, pfMET="<<pfMET<<" pT="<<phoEt->at(p)<<" sel: "<<sel<<" sys: "<<sys<<std::endl;
      // std::cout<<"  CHiso:  "<<TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0)<<std::endl;
       tmpCand.push_back(p);
     }
    }
  return tmpCand;
}


// Effective area to be needed in PF Iso for photon ID
// https://indico.cern.ch/event/455258/contribution/0/attachments/1173322/1695132/SP15_253rd.pdf -- slide-5

// worst charged hadron isolation EA
Double_t postAnalyzer_QCD::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.078;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.089;
  return EffectiveArea;
}

Double_t postAnalyzer_QCD::EAneutral(Double_t eta){
  Float_t EffectiveArea = 0.;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0599;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0819;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0696;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0360;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0360;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0462;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0656;

  return EffectiveArea;
}

Double_t postAnalyzer_QCD::EAphoton(Double_t eta){
  Float_t EffectiveArea = 0.;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.1271;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1101;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0756;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.1175;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.1498;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.1857;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.2183;

  return EffectiveArea;
}
