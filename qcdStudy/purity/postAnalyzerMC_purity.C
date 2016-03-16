#define postAnalyzerMC_purity_cxx
#include "postAnalyzerMC_purity.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void postAnalyzerMC_purity::Loop(TString outfilename, Bool_t isMC, Double_t lumi, Double_t nrEvents, Double_t crossSec)
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
  
   int incptbin = ptbinnames.size()-1; // since bin labls start at 0 ( = 6 : 175-2000 )
   int selptbin = ptbinnames.size()-2; // ( = 5 : 175-1000  )
   int lstptbin = ptbinnames.size()-3; // ( = 4 : 700-1000)

   // clear collections
   recPCvint.clear(); 
   genPCvint.clear(); 
   mchPCvint.clear(); 
   nmhPCvint.clear(); 

   // Find all reco photons passing selections
   recPCvint = pcPassSel();

   //Find all GEN photons
   for(int i = 0; i < nMC; i++)
   {
    //Find all generated photons with appropriate status flag
    if(mcPID->at(i) == 22 && (((mcStatusFlag->at(i)>>0)&1)==1 || ((mcStatusFlag->at(i)>>1)&1)==1))
    {
     //nGenPho++;
     genPCvint.push_back(i);
    }
   }

   //Loop over all reconstructed photons to find (GEN)matched and nonmatched photons
   for(int i = 0; i < recPCvint.size(); i++)
   {
    int phonr = recPCvint[i];
    //std::cout<<" photon number: "<<phonr<<std::endl;
    double dRmin = 99.99;
    double dRtemp = 99.99;
    // technically a gen photon could be matched to more than one RECO photon in this method..
    for(int k = 0; k < genPCvint.size(); k++) 
    {
     int gennr = genPCvint[k];
     dRtemp = dR(mcEta->at(gennr),mcPhi->at(gennr),phoEta->at(phonr),phoPhi->at(phonr));
     if(dRtemp < dRmin) { dRmin = dRtemp; }
    }
    if(dRmin < 0.1) { mchPCvint.push_back(phonr); }
    if(dRmin > 0.1) { nmhPCvint.push_back(phonr); }
   }

// fill histograms
    // Fill Numerator (matched) Signal Histograms
    if( mchPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<mchPCvint.size(); ++k){ // go through each photon passing sel
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(mchPCvint.at(k)) > ptbins[ptb]) &&
          (phoEt->at(mchPCvint.at(k)) < ptbins[ptb+1])
         ){
            //nc++;
            FillSigHistograms(ptb, mchPCvint.at(k), event_weight);
            FillDenHistograms(ptb, mchPCvint.at(k), event_weight);
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(mchPCvint.at(k)) > ptbins[0]) &&
         (phoEt->at(mchPCvint.at(k)) < ptbins[ptbins.size()-1])
        ){
       FillSigHistograms(selptbin, mchPCvint.at(k), event_weight);
       FillDenHistograms(selptbin, mchPCvint.at(k), event_weight);
      }
      // fully inclusive in pT (depending only on selections for mchPCvint)
      FillSigHistograms(lstptbin, mchPCvint.at(k), event_weight);
      FillDenHistograms(lstptbin, mchPCvint.at(k), event_weight);
    }
   }

    // Fill Numerator (unmatched) Signal Histograms
    if( nmhPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<nmhPCvint.size(); ++k){ // go through each photon passing sel
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(nmhPCvint.at(k)) > ptbins[ptb]) &&
          (phoEt->at(nmhPCvint.at(k)) < ptbins[ptb+1])
         ){
            //nc++;
            FillBkgHistograms(ptb, nmhPCvint.at(k), event_weight);
            FillDenHistograms(ptb, nmhPCvint.at(k), event_weight);
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(nmhPCvint.at(k)) > ptbins[0]) &&
         (phoEt->at(nmhPCvint.at(k)) < ptbins[ptbins.size()-1])
        ){
       FillBkgHistograms(selptbin, nmhPCvint.at(k), event_weight);
       FillDenHistograms(selptbin, nmhPCvint.at(k), event_weight);
      }
      // fully inclusive in pT (depending only on selections for nmhPCvint)
      FillBkgHistograms(lstptbin, nmhPCvint.at(k), event_weight);
      FillDenHistograms(lstptbin, nmhPCvint.at(k), event_weight);
    }
   }

//// end fill histograms

 } //end loop through entries

 // write these histograms to file
//  std::cout<<std::endl;
//  std::cout<<"Total Passing Numerator: "<<nc<<"  Total Passing Denominator: "<<dc<<std::endl;
//  std::cout<<"made it through, about to write"<<std::endl;

//  std::ofstream logfile_numsig (outfilename+"_numsig.log", std::ofstream::out);
//  //logfile_numsig<<"============================="<<std::endl;
//  //logfile_numsig<<" Numerator Signal Selections"<<std::endl;
//  //logfile_numsig<<"============================="<<std::endl;
//  for (unsigned int i=0; i<events_numsig.size(); ++i){ logfile_numsig<<"\t "<<run_numsig[i]<<" : "<<events_numsig[i]<<" : "<<lumi_numsig[i]<<" : 0"<<std::endl; }
//  //for (unsigned int i=0; i<events_numsig.size(); ++i){ logfile_numsig<<"\t "<<run_numsig[i]<<"\t "<<lumi_numsig[i]<<"\t "<<events_numsig[i]<<std::endl; }
//  logfile_numsig.close();
//
//  std::ofstream logfile_numqcd (outfilename+"_numqcd.log", std::ofstream::out);
//  //logfile_numqcd<<"============================="<<std::endl;
//  //logfile_numqcd<<" Numerator QCD Selections"<<std::endl;
//  //logfile_numqcd<<"============================="<<std::endl;
//  for (unsigned int i=0; i<events_numqcd.size(); ++i){ logfile_numqcd<<"\t "<<run_numqcd[i]<<" : "<<events_numqcd[i]<<" : "<<lumi_numqcd[i]<<" : 0"<<std::endl; }
//  //for (unsigned int i=0; i<events_numqcd.size(); ++i){ logfile_numqcd<<"\t "<<run_numqcd[i]<<"\t "<<lumi_numqcd[i]<<"\t "<<events_numqcd[i]<<std::endl; }
//  logfile_numqcd.close();
//
//  std::ofstream logfile_denomi (outfilename+"_denomi.log", std::ofstream::out);
//  //logfile_denomi<<"============================="<<std::endl;
//  //logfile_denomi<<" Denominator Selections"<<std::endl;
//  //logfile_denomi<<"============================="<<std::endl;
//  for (unsigned int i=0; i<events_denomi.size(); ++i){ logfile_denomi<<"\t "<<run_denomi[i]<<" : "<<events_denomi[i]<<" : "<<lumi_denomi[i]<<" : 0"<<std::endl; }
//  //for (unsigned int i=0; i<events_denomi.size(); ++i){ logfile_denomi<<"\t "<<run_denomi[i]<<"\t "<<lumi_denomi[i]<<"\t "<<events_denomi[i]<<std::endl; }
//  logfile_denomi.close();
  

 TFile *outfile = new TFile(outfilename,"RECREATE");
 outfile->cd();
 for(unsigned int i=0; i<ptbinnames.size(); ++i){
   WriteHistograms(i);
 }
 outfile->Close();
 sw.Stop();
 std::cout<<"Real Time: "<<sw.RealTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"CPU Time: "<<sw.CpuTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"Done"<<std::endl;

} //end Loop()

// photon candidates passing signal selection in low MET region - (no sieie or (w)chiso cuts)
std::vector<int> postAnalyzerMC_purity::pcPassSel( double phoPtLo, double phoPtHi, double phoEtaMax ){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool passTriggers;
  bool noncoll;
  bool passKinematics;
  bool passMET;

  bool passHoEPSeed;
  bool passPhoNHMedIso;
  bool passCHMedIso; 

  bool passPhotonID;

  passTriggers = 
   ((HLTPho>>7&1) == 1) ||
   ((HLTPho>>8&1) == 1) ||
   ((HLTPho>>9&1) == 1) ||
   ((HLTPho>>10&1) == 1) ||
   ((HLTPho>>11&1) == 1) ||
   ((HLTPho>>12&1) == 1) ||
   ((HLTPho>>22&1) == 1) ;

  //Loop over photons
  for(int p=0;p<nPho;p++)  // nPho from ntuple (should = phoE->length() )
    {
     // from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#SPRING15_selections_25_ns


     // non collision backgrounds
     noncoll = (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001; // isMC
     //noncoll = fabs((*phoseedTimeFull5x5)[p]) < 3. && (*phomipTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001;

     passKinematics = (
                       ( (*phoEt)[p] > phoPtLo  ) &&
                       ( (*phoEt)[p] <= phoPtHi ) &&
                       ( fabs((*phoSCEta)[p]) < phoEtaMax)
                      );

     // QCD 
     passMET = pfMET < 30.;
     //passMET = kTRUE; //pfMET < 30.;
//     passMET = pfMET < 15.;
//     passMET = pfMET < 45.;

     // H on E, and pixel seed 
     passHoEPSeed = (
                     ((*phoHoverE)[p] < 0.05 ) &&
                     ((*phohasPixelSeed)[p] ==  0 )
                    );
     // photon / neutral hadron isolation
     passPhoNHMedIso = (
                        ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < 
                           1.06 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))
                        )  &&
                        ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < 
                           0.28 + (0.0053 * (*phoEt)[p])
                        )
                       );

//     // (worst) charged hadron isolation
//     passCHMedIso = ( TMath::Max( ( (*phoPFChWorstIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.37 );
     passCHMedIso = kTRUE;

     passPhotonID = passMET && passHoEPSeed && passPhoNHMedIso && passCHMedIso;
     if( passTriggers && noncoll && passPhotonID && passKinematics){
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
Double_t postAnalyzerMC_purity::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.078;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.089;
  return EffectiveArea;
}

Double_t postAnalyzerMC_purity::EAneutral(Double_t eta){
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

Double_t postAnalyzerMC_purity::EAphoton(Double_t eta){
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

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
double postAnalyzerMC_purity::DeltaPhi(double phi1, double phi2)
{
        double pi = TMath::Pi();
        double dphi = fabs(phi1-phi2);
        if(dphi>pi)
                dphi = 2.0*pi - dphi;
        return dphi;
}

double postAnalyzerMC_purity::dR(double eta1, double phi1, double eta2, double phi2)
{
        double deltaeta = abs(eta1 - eta2);
        double deltaphi = DeltaPhi(phi1, phi2);
        double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
        return deltar;
}

