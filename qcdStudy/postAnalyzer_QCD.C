#define postAnalyzer_QCD_cxx
#include "postAnalyzer_QCD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void postAnalyzer_QCD::Loop(TString outfilename, Bool_t isMC, Double_t weight)
{

 TStopwatch sw; 
 sw.Start();

 if (fChain == 0) return;

 //Float_t ptbins[12] = {75.,100.,125.,145.,155.,165.,175.,190.,250.,400.,700.,1000.};
 ptbins.clear();
 //ptbins.push_back(75);
 //ptbins.push_back(100);
 //ptbins.push_back(125);
 //ptbins.push_back(145);
 //ptbins.push_back(155);
 //ptbins.push_back(165);
 ptbins.push_back(175);
 ptbins.push_back(190);
 ptbins.push_back(250);
 //ptbins.push_back(400);
 //ptbins.push_back(700);
 ptbins.push_back(1000);

 Long64_t nentries = fChain->GetEntriesFast();

 Long64_t nbytes = 0, nb = 0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
  if (jentry%20000 == 0)
    {
      std::cout<<"Starting entry "<<jentry<<"/"<<(nentries)<<" at "<<sw.RealTime()<<" RealTime, "<<sw.CpuTime()<<" CpuTime"<<std::endl;
      sw.Continue();
    }
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   nbytes += nb;

  sigPCvint.clear();
  bkgPCvint.clear();
  denPCvint.clear();

  //=1.0 for real data
  event_weight=1.0;
  if(isMC){ event_weight=weight; }
  
  // if event passes MonoPhoton triggers
  if( 
   (HLTPho>>7&1 == 1) ||
   (HLTPho>>8&1 == 1) ||
   (HLTPho>>9&1 == 1) ||
   (HLTPho>>10&1 == 1) ||
   (HLTPho>>11&1 == 1) ||
   (HLTPho>>12&1 == 1) ||
   (HLTPho>>22&1 == 1) )
  {   
   // vector of ints, each int corresponds to location in vector of photons of photon passing cut
   sigPCvint = pcPassSel(0); // passes signal selection (no sieie cut)
   bkgPCvint = pcPassSel(1); // passes background selection (no sieie cut)
   denPCvint = pcPassSel(2); // passes denominator selection (no sieie cut)
  
   // QCD Cut
   passQCD = pfMET < 30; 
  
   // Fill Numerator Signal Histograms
   if( sigPCvint.size()>0 && passQCD ){ // if any photon indexes passed sig selection
    //std::cout<<" passed QCD, pT="<<phoEt->at(sigPCvint[0])<<std::endl;

    for(unsigned int ptb=0; ptb<(ptbins.size()-1); ++ptb){ // break into pT bins
     //std::cout<<"  pt bin="<<ptb<<" for "<<ptbins.at(ptb)<<"-"<<ptbins.at(ptb+1)<<std::endl;
     if(
        (phoEt->at(sigPCvint[0]) > ptbins[ptb]) &&
        (phoEt->at(sigPCvint[0]) < ptbins[ptb+1])
       ){
      //std::cout<<"writing signal"<<std::endl; std::cout<<phoEt->at(sigPCvint[0])<<std::endl;
      FillSigHistograms(ptb, sigPCvint[0], event_weight);
     }
    }
    if(  // also do an inclusive pT plot
       (phoEt->at(sigPCvint[0]) > ptbins[0]) &&
       (phoEt->at(sigPCvint[0]) < ptbins[ptbins.size()-1])
      ){
     FillSigHistograms(ptbins.size()-1, sigPCvint[0], event_weight);
    }
   }

   // Fill Numerator Background Histograms
   if( bkgPCvint.size()>0 && passQCD ){ // if any photon indexes passed bkg selection
    for(unsigned int ptb=0; ptb<(ptbins.size()-1); ++ptb){ // break into pT bins
     if(
        (phoEt->at(bkgPCvint[0]) > ptbins[ptb]) &&
        (phoEt->at(bkgPCvint[0]) < ptbins[ptb+1])
       ){
      FillBkgHistograms(ptb, bkgPCvint[0], event_weight);
     }
    }
    if(  // also do an inclusive pT plot
       (phoEt->at(bkgPCvint[0]) > ptbins[0]) &&
       (phoEt->at(bkgPCvint[0]) < ptbins[ptbins.size()-1])
      ){
     FillBkgHistograms(ptbins.size()-1, bkgPCvint[0], event_weight);
    }
   }

   // Fill Denominator Histograms
   if( denPCvint.size()>0 && passQCD ){ // if any photon indexes passed bkg selection
    for(unsigned int ptb=0; ptb<(ptbins.size()-1); ++ptb){ // break into pT bins
     if(
        (phoEt->at(denPCvint[0]) > ptbins[ptb]) &&
        (phoEt->at(denPCvint[0]) < ptbins[ptb+1])
       ){
      FillDenHistograms(ptb, denPCvint[0], event_weight);
     }
    }
    if(  // also do an inclusive pT plot
       (phoEt->at(denPCvint[0]) > ptbins[0]) &&
       (phoEt->at(denPCvint[0]) < ptbins[ptbins.size()-1])
      ){
     FillDenHistograms(ptbins.size()-1, denPCvint[0], event_weight);
    }
   }

  } //end if passes MonoPhoton triggers
 } //end loop through entries

 // write these histograms to file
 TFile *outfile = new TFile(outfilename,"RECREATE");
 outfile->cd();
 WriteHistograms(ptbins.size());
 outfile->Close();
 sw.Stop();
 std::cout<<"Real Time: "<<sw.RealTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"CPU Time: "<<sw.CpuTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"Done"<<std::endl;

} //end Loop()

//   pcPassSel photons passing selections: ( num_sig(0), num_bkg(1), den(2) )
std::vector<int> postAnalyzer_QCD::pcPassSel(int sel, double phoPtLo, double phoPtHi, double phoEtaMax){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool kinematic;
  bool photonId;
  //Loop over photons
  for(int p=0;p<nPho;p++)  // nPho from ntuple (should = phoE->length() )
    {
     kinematic = (
                  ( (*phoEt)[p] > phoPtLo  ) &&
                  ( (*phoEt)[p] <= phoPtHi ) &&
                  ( fabs((*phoSCEta)[p]) < phoEtaMax)
                 );

     // from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#SPRING15_selections_25_ns
     if(sel==0){  // numerator signal
      photonId = (
                  ((*phoHoverE)[p] < 0.05 ) &&
                  //((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0102 ) &&
                  ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.37 )  &&
                  ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < 
                     1.06 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))
                  )  &&
                  ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < 
                     0.28 + (0.0053 * (*phoEt)[p])
                  )
                 );
     }
     else if(sel==1){ // numerator background
      photonId = (
                  ((*phoHoverE)[p] < 0.05 ) &&
                  //((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0102 ) &&
                  ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) > 5. )  &&
                  ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 10. )  &&
                  ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < 
                     1.06 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))
                  )  &&
                  ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < 
                     0.28 + (0.0053 * (*phoEt)[p])
                  )
                 );
     }
     else if(sel==2){ // denominator

      double  vloosePFCharged= TMath::Min(5.0*(3.32) , 0.20*(*phoEt)[p]);
      double  vloosePFPhoton = TMath::Min(5.0*(0.81+ (0.0053 * (*phoEt)[p])) , 0.20*(*phoEt)[p]);
      double  vloosePFNeutral= TMath::Min(5.0*(1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))) , 0.20*(*phoEt)[p]);
     
      bool passVLooseSel = ( 
                            ((*phoHoverE)[p]                <  0.05   ) &&
                            ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0102 ) &&
                            ((*phohasPixelSeed)[p]              ==  0      ) &&
                            ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < vloosePFCharged )  &&  
                            ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < vloosePFNeutral )  &&  
                            ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < vloosePFPhoton )
                           ); 
      //bool passLooseSel = (
      //                      ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 3.32 )  ||  
      //                      ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) <
      //                        (1.92 + (0.014* (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))  ||  
      //                      ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) <
      //                        (0.81 + (0.0053 * (*phoEt)[p])) )
      //                     ); 
      //photonId = ( !passLooseSel && passVLooseSel );
      bool failLooseSel = (
                            ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) > 3.32 )  ||  
                            ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) >
                              (1.92 + (0.014* (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))  || 
                            ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) >
                              (0.81 + (0.0053 * (*phoEt)[p])) )
                           ); 
      photonId = ( failLooseSel && passVLooseSel );
     }
     if(photonId && kinematic){
       //std::cout<<" Found a photon, pfMET="<<pfMET<<" pT="<<phoEt->at(p)<<" sel: "<<sel<<std::endl;
       tmpCand.push_back(p);
     }
    }
  return tmpCand;
}


// Effective area to be needed in PF Iso for photon ID
// https://indico.cern.ch/event/455258/contribution/0/attachments/1173322/1695132/SP15_253rd.pdf -- slide-5
Double_t postAnalyzer_QCD::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0456;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0500;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0340;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0383;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0339;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0303;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0240;

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
