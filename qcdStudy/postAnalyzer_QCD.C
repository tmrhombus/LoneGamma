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

   // sysbinnames = "" "_sbUP" "_sbDown" "_metUP" "_metDown" "_binUP" "_binDown" "_noPiso"
   // ptbins = 175, 190, 250, 400, 700, 1000
   // ptbinnames = "175to190" "190to250" "250to400" "400to700" "700to1000" "175to1000" "allpt"

   // vector of vector of ints, each int corresponds to location in vector of photons of photon passing cut
   // one such vector for each pt and systematic bin
   // start by filling inclusive (last) pt bin for systematics, then divvy up

   int inclptbin = ptbins.size()-1;
   int lastptbin = ptbinnames.size()-1;
   int lastsysbin = sysbinnames.size();
   for(unsigned int sysb=0; sysb<lastsysbin; ++sysb){
    sigPCvint[lastptbin][sysb] = pcPassSel(0,sysb); // passes signal selection (no sieie cut)
    bkgPCvint[lastptbin][sysb] = pcPassSel(1,sysb); // passes background selection (no sieie cut)
    denPCvint[lastptbin][sysb] = pcPassSel(2,sysb); // passes denominator selection (no sieie cut)
   }

// fill histograms
   for(unsigned int sysb=0; sysb<lastsysbin; ++sysb){
    // Fill Numerator Signal Histograms
    if( sigPCvint[lastptbin][sysb].size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
      if(
         (phoEt->at(sigPCvint[lastptbin][sysb].at(0)) > ptbins[ptb]) &&
         (phoEt->at(sigPCvint[lastptbin][sysb].at(0)) < ptbins[ptb+1])
        ){
       FillSigHistograms(ptb, sysb, sigPCvint[lastptbin][sysb].at(0), event_weight);
      }
     }
     if(  // do an inclusive pT plot from bins
        (phoEt->at(sigPCvint[lastptbin][sysb].at(0)) > ptbins[0]) &&
        (phoEt->at(sigPCvint[lastptbin][sysb].at(0)) < ptbins[inclptbin])
       ){
      FillSigHistograms(lastptbin-1, sysb, sigPCvint[lastptbin][sysb].at(0), event_weight);
     }
     // and one fully inclusive in pT
     FillSigHistograms(lastptbin, sysb, sigPCvint[lastptbin][sysb].at(0), event_weight);
    }

    // Fill Numerator Background Histograms
    if( bkgPCvint[lastptbin][sysb].size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
      if(
         (phoEt->at(bkgPCvint[lastptbin][sysb].at(0)) > ptbins[ptb]) &&
         (phoEt->at(bkgPCvint[lastptbin][sysb].at(0)) < ptbins[ptb+1])
        ){
       FillBkgHistograms(ptb, sysb, bkgPCvint[lastptbin][sysb].at(0), event_weight);
      }
     }
     if(  // do an inclusive pT plot from bins
        (phoEt->at(bkgPCvint[lastptbin][sysb].at(0)) > ptbins[0]) &&
        (phoEt->at(bkgPCvint[lastptbin][sysb].at(0)) < ptbins[inclptbin])
       ){
      FillBkgHistograms(lastptbin-1, sysb, bkgPCvint[lastptbin][sysb].at(0), event_weight);
     }
     // and one fully inclusive in pT
     FillBkgHistograms(lastptbin, sysb, bkgPCvint[lastptbin][sysb].at(0), event_weight);
    }

    // Fill Denominator Histograms
    if( denPCvint[lastptbin][sysb].size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
      if(
         (phoEt->at(denPCvint[lastptbin][sysb].at(0)) > ptbins[ptb]) &&
         (phoEt->at(denPCvint[lastptbin][sysb].at(0)) < ptbins[ptb+1])
        ){
       FillDenHistograms(ptb, sysb, denPCvint[lastptbin][sysb].at(0), event_weight);
      }
     }
     if(  // do an inclusive pT plot from bins
        (phoEt->at(denPCvint[lastptbin][sysb].at(0)) > ptbins[0]) &&
        (phoEt->at(denPCvint[lastptbin][sysb].at(0)) < ptbins[inclptbin])
       ){
      FillDenHistograms(lastptbin-1, sysb, denPCvint[lastptbin][sysb].at(0), event_weight);
     }
     // and one fully inclusive in pT
     FillDenHistograms(lastptbin, sysb, denPCvint[lastptbin][sysb].at(0), event_weight);
    }
   } // for each sysb in sysbins
// end fill histograms

  } //end if passes MonoPhoton triggers
 } //end loop through entries

 // write these histograms to file
  std::cout<<"made it through, about to write"<<std::endl;
 TFile *outfile = new TFile(outfilename,"RECREATE");
 outfile->cd();
 for(unsigned int i=0; i<ptbinnames.size(); ++i){
  for(unsigned int j=0; j<sysbinnames.size(); ++j){
   WriteHistograms(i,j);
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
std::vector<int> postAnalyzer_QCD::pcPassSel(int sel, int sys, double phoPtLo, double phoPtHi, double phoEtaMax){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool passKinematics;
  bool passMET;

  bool passHoEPSeed;
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
     //passCutSieie = ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0102 );  // don't need for us

     passKinematics = (
                       ( (*phoEt)[p] > phoPtLo  ) &&
                       ( (*phoEt)[p] <= phoPtHi ) &&
                       ( fabs((*phoSCEta)[p]) < phoEtaMax)
                      );


     passMET = pfMET < 30.;
     if(sys==3){ passMET = pfMET < 35.;}
     if(sys==4){ passMET = pfMET < 20.;}

     // nsig, nbkg, deno
     passHoEPSeed = (
                     ((*phoHoverE)[p] < 0.05 ) &&
                     ((*phohasPixelSeed)[p] ==  0 )
                    );
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
     passCHMedIso = ( TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.37 );

     // nbkg
     chIsoLB = 5.;
     chIsoUB = 10.;
     if(sys==1){chIsoUB = 12.;}
     if(sys==2){chIsoUB = 8.;}
     passCHBkgIso = (
                     ( TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) > chIsoLB )  &&
                     ( TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < chIsoUB ) 
                    );

     // deno
     vloosePFCharged= TMath::Min(5.0*(3.32) , 0.20*(*phoEt)[p]);
     vloosePFPhoton = TMath::Min(5.0*(0.81+ (0.0053 * (*phoEt)[p])) , 0.20*(*phoEt)[p]);
     vloosePFNeutral= TMath::Min(5.0*(1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))) , 0.20*(*phoEt)[p]);
     passVLooseIso = ( 
                      ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < vloosePFCharged )  &&  
                      ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < vloosePFNeutral )  &&  
                      ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < vloosePFPhoton )
                     ); 
     passLoosePIso = 
                     ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) <
                       (0.81 + (0.0053 * (*phoEt)[p])) );
     if(sys==7){passLoosePIso= true;}
     passLooseIso = (  // deno must fail this cut
                     ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 3.32 )  &&
                     ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) <
                       (1.92 + (0.014* (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))  &&
                     passLoosePIso
                    );

     if(sel==0){  // numerator signal
      passPhotonID = passMET && passHoEPSeed && passPhoNHMedIso && passCHMedIso;
     }
     else if(sel==1){ // numerator background
      passPhotonID = passMET && passHoEPSeed && passPhoNHMedIso && passCHBkgIso;
     }
     else if(sel==2){ // denominator
      passPhotonID = passMET && passHoEPSeed && !passLooseIso && passVLooseIso;
     }
     if(passPhotonID && passKinematics){
       std::cout<<" Found a photon, pfMET="<<pfMET<<" pT="<<phoEt->at(p)<<" sel: "<<sel<<" sys: "<<sys<<std::endl;
       std::cout<<"  CHiso:  "<<TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0)<<std::endl;
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
