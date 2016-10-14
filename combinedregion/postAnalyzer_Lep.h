
#ifndef postAnalyzer_Lep_h
#define postAnalyzer_Lep_h

#include "postAnalyzer_Base.h"


class postAnalyzer_Lep : public postAnalyzer_Base {

 public:
   TH1F h_ele_et[7],
        h_ele_uncorret[7],
        h_ele_eta[7],
        h_ele_sieieF5x5[7],
        h_ele_pfMET[7],
        h_ele_leptoMET[7],
        h_ele_dilep_mass[7],
        h_ele_diele_mass[7],
        h_ele_dimu_mass[7],
        h_ele_et5[7],
        h_ele_et15[7],
        h_ele_et25[7],
        h_ele_et55[7],
        h_ele_et75[7],
        h_ele_et165[7],
        h_ele_et275[7],
        h_ele_met5[7],
        h_ele_met15[7],
        h_ele_met25[7],
        h_ele_met55[7],
        h_ele_met75[7],
        h_ele_met165[7],
        h_ele_met275[7];


   TH1F h_mu_et[7],
        h_mu_uncorret[7],
        h_mu_eta[7],
        h_mu_sieieF5x5[7],
        h_mu_pfMET[7],
        h_mu_leptoMET[7],
        h_mu_dilep_mass[7],
        h_mu_diele_mass[7],
        h_mu_dimu_mass[7],
        h_mu_et5[7],
        h_mu_et15[7],
        h_mu_et25[7],
        h_mu_et55[7],
        h_mu_et75[7],
        h_mu_et165[7],
        h_mu_et275[7],
        h_mu_met5[7],
        h_mu_met15[7],
        h_mu_met25[7],
        h_mu_met55[7],
        h_mu_met75[7],
        h_mu_met165[7],
        h_mu_met275[7];

   // variables constructed for histograms
   //TLorentzVector fourVec_e1, fourVec_e2; 
   //TLorentzVector fourVec_m1, fourVec_m2; 
   TLorentzVector fourVec_ee, fourVec_mm, fourVec_ll; 
   TLorentzVector fourVec_l1, fourVec_l2; 
   TLorentzVector fourVec_met; 
   TLorentzVector fourVec_leptomet;
   Double_t dilep_mass;

   Double_t leptoMET;
   Double_t leptoMEPhi;

   //TLorentzVector fourVec_e1, fourVec_e2; 
   //TLorentzVector fourVec_m1, fourVec_m2; 
   TLorentzVector fourVec_genee, fourVec_genmm, fourVec_genll; 
   TLorentzVector fourVec_genl1, fourVec_genl2; 
   TLorentzVector fourVec_genMET; 
   TLorentzVector fourVec_genLeptoMET;
   Double_t gen_dilep_mass;

   Double_t genLeptoMET;
   Double_t genLeptoMEPhi;


   virtual void     InitLep();
   Bool_t           FillSigHistogramsLep(int ptbin, int sysbin, int photonIndex, double weight, bool passM);
   virtual void     callFillSigHistLep(int selbin, int lastptbin, int inclptbin, int candphotonindex, float event_weight, bool passM);
   Bool_t           WriteHistogramsLep(int ptbin, int sysbin);

};

void postAnalyzer_Lep::InitLep()
{

   std::cout<<"Initializing leptons"<<std::endl;

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
     TString histname_mu_et          = "h_mu_et_"        +ptbinnames[i];
     TString histname_mu_uncorret    = "h_mu_uncorret_"  +ptbinnames[i];
     TString histname_mu_eta         = "h_mu_eta_"       +ptbinnames[i];
     TString histname_mu_sieieF5x5   = "h_mu_sieieF5x5_" +ptbinnames[i];
     TString histname_mu_pfMET       = "h_mu_pfMET_"     +ptbinnames[i];
     TString histname_mu_leptoMET    = "h_mu_leptoMET_"  +ptbinnames[i];
     TString histname_mu_dilep_mass  = "h_mu_dilep_mass_"+ptbinnames[i];
     TString histname_mu_diele_mass  = "h_mu_diele_mass_"+ptbinnames[i];
     TString histname_mu_dimu_mass   = "h_mu_dimu_mass_" +ptbinnames[i];
     TString histname_mu_et5         = "h_mu_et5_"       +ptbinnames[i];
     TString histname_mu_et15        = "h_mu_et15_"      +ptbinnames[i];
     TString histname_mu_et25        = "h_mu_et25_"      +ptbinnames[i];
     TString histname_mu_et55        = "h_mu_et55_"      +ptbinnames[i];
     TString histname_mu_et75        = "h_mu_et75_"      +ptbinnames[i];
     TString histname_mu_et165       = "h_mu_et165_"     +ptbinnames[i];
     TString histname_mu_et275       = "h_mu_et275_"     +ptbinnames[i];
     TString histname_mu_met5        = "h_mu_met5_"      +ptbinnames[i];
     TString histname_mu_met15       = "h_mu_met15_"     +ptbinnames[i];
     TString histname_mu_met25       = "h_mu_met25_"     +ptbinnames[i];
     TString histname_mu_met55       = "h_mu_met55_"     +ptbinnames[i];
     TString histname_mu_met75       = "h_mu_met75_"     +ptbinnames[i];
     TString histname_mu_met165      = "h_mu_met165_"    +ptbinnames[i];
     TString histname_mu_met275      = "h_mu_met275_"    +ptbinnames[i];

     TString histname_ele_et         = "h_ele_et_"        +ptbinnames[i];
     TString histname_ele_uncorret   = "h_ele_uncorret_"  +ptbinnames[i];
     TString histname_ele_eta        = "h_ele_eta_"       +ptbinnames[i];
     TString histname_ele_sieieF5x5  = "h_ele_sieieF5x5_" +ptbinnames[i];
     TString histname_ele_pfMET      = "h_ele_pfMET_"     +ptbinnames[i];
     TString histname_ele_leptoMET   = "h_ele_leptoMET_"  +ptbinnames[i];
     TString histname_ele_dilep_mass = "h_ele_dilep_mass_"+ptbinnames[i];
     TString histname_ele_diele_mass = "h_ele_diele_mass_"+ptbinnames[i];
     TString histname_ele_dimu_mass  = "h_ele_dimu_mass_" +ptbinnames[i];
     TString histname_ele_et5        = "h_ele_et5_"       +ptbinnames[i];
     TString histname_ele_et15       = "h_ele_et15_"      +ptbinnames[i];
     TString histname_ele_et25       = "h_ele_et25_"      +ptbinnames[i];
     TString histname_ele_et55       = "h_ele_et55_"      +ptbinnames[i];
     TString histname_ele_et75       = "h_ele_et75_"      +ptbinnames[i];
     TString histname_ele_et165      = "h_ele_et165_"     +ptbinnames[i];
     TString histname_ele_et275      = "h_ele_et275_"     +ptbinnames[i];
     TString histname_ele_met5       = "h_ele_met5_"      +ptbinnames[i];
     TString histname_ele_met15      = "h_ele_met15_"     +ptbinnames[i];
     TString histname_ele_met25      = "h_ele_met25_"     +ptbinnames[i];
     TString histname_ele_met55      = "h_ele_met55_"     +ptbinnames[i];
     TString histname_ele_met75      = "h_ele_met75_"     +ptbinnames[i];
     TString histname_ele_met165     = "h_ele_met165_"    +ptbinnames[i];
     TString histname_ele_met275     = "h_ele_met275_"    +ptbinnames[i];

     Float_t binspt[] = { 175., 195., 250., 400., 700., 1000. };
     Int_t  binnumpt = sizeof(binspt)/sizeof(Float_t) - 1; // or just = 9

     // reserve histograms
     // Electron
     h_ele_et[i].Clear();
     h_ele_et[i] = TH1F(histname_ele_et,"Photon Transverse Energy",binnumpt,binspt);
     //h_ele_et[i] = TH1F(histname_ele_et,"Photon Transverse Energy",165,175.,1000.);
     h_ele_et[i].Sumw2();
     //
     h_ele_uncorret[i].Clear();
     h_ele_uncorret[i] = TH1F(histname_ele_uncorret,"Uncorrected Photon Transverse Energy",binnumpt,binspt);
     //h_ele_uncorret[i] = TH1F(histname_ele_uncorret,"Uncorrected Photon Transverse Energy",165,175.,1000.);
     h_ele_uncorret[i].Sumw2();
     //
     h_ele_eta[i].Clear();
     h_ele_eta[i] = TH1F(histname_ele_eta,"Leading Photon Eta",10,-2.,2.);
     h_ele_eta[i].Sumw2();
     //
     h_ele_sieieF5x5[i].Clear();
     h_ele_sieieF5x5[i] = TH1F(histname_ele_sieieF5x5,"Leading Photon SigmaIetaIeta",50,0.,0.025);
     h_ele_sieieF5x5[i].Sumw2();
     //
     h_ele_pfMET[i].Clear();
     h_ele_pfMET[i] = TH1F(histname_ele_pfMET,"ParticleFlow MET",binnumpt,binspt);
     h_ele_pfMET[i].Sumw2();
     //
     h_ele_leptoMET[i].Clear();
     h_ele_leptoMET[i] = TH1F(histname_ele_leptoMET,"PF MET + dilepton",binnumpt,binspt);
     h_ele_leptoMET[i].Sumw2();
     //
     h_ele_dilep_mass[i].Clear();
     h_ele_dilep_mass[i] = TH1F(histname_ele_dilep_mass,"mass dilepton (MET)",15,60.,120.);
     h_ele_dilep_mass[i].Sumw2();
     //
     h_ele_dimu_mass[i].Clear();
     h_ele_dimu_mass[i] = TH1F(histname_ele_dimu_mass,"mass dimuon (MET)",15,60.,120.);
     h_ele_dimu_mass[i].Sumw2();
     //
     h_ele_diele_mass[i].Clear();
     h_ele_diele_mass[i] = TH1F(histname_ele_diele_mass,"mass diele (MET)",15,60.,120.);
     h_ele_diele_mass[i].Sumw2();
     //
     h_ele_et5[i].Clear();
     h_ele_et5[i] = TH1F(histname_ele_et5,"Photon Transverse Energy",165,175.,1000.);
     h_ele_et5[i].Sumw2();
     //
     h_ele_et15[i].Clear();
     h_ele_et15[i] = TH1F(histname_ele_et15,"Photon Transverse Energy",55,175.,1000.);
     h_ele_et15[i].Sumw2();
     //
     h_ele_et25[i].Clear();
     h_ele_et25[i] = TH1F(histname_ele_et25,"Photon Transverse Energy",33,175.,1000.);
     h_ele_et25[i].Sumw2();
     //
     h_ele_et55[i].Clear();
     h_ele_et55[i] = TH1F(histname_ele_et55,"Photon Transverse Energy",15,175.,1000.);
     h_ele_et55[i].Sumw2();
     //
     h_ele_et75[i].Clear();
     h_ele_et75[i] = TH1F(histname_ele_et75,"Photon Transverse Energy",11,175.,1000.);
     h_ele_et75[i].Sumw2();
     //
     h_ele_et165[i].Clear();
     h_ele_et165[i] = TH1F(histname_ele_et165,"Photon Transverse Energy",5,175.,1000.);
     h_ele_et165[i].Sumw2();
     //
     h_ele_et275[i].Clear();
     h_ele_et275[i] = TH1F(histname_ele_et275,"Photon Transverse Energy",3,175.,1000.);
     h_ele_et275[i].Sumw2();
     //
     h_ele_met5[i].Clear();
     h_ele_met5[i] = TH1F(histname_ele_met5,"Missing Transverse Energy",165,175.,1000.);
     h_ele_met5[i].Sumw2();
     //
     h_ele_met15[i].Clear();
     h_ele_met15[i] = TH1F(histname_ele_met15,"Missing Transverse Energy",55,175.,1000.);
     h_ele_met15[i].Sumw2();
     //
     h_ele_met25[i].Clear();
     h_ele_met25[i] = TH1F(histname_ele_met25,"Missing Transverse Energy",33,175.,1000.);
     h_ele_met25[i].Sumw2();
     //
     h_ele_met55[i].Clear();
     h_ele_met55[i] = TH1F(histname_ele_met55,"Missing Transverse Energy",15,175.,1000.);
     h_ele_met55[i].Sumw2();
     //
     h_ele_met75[i].Clear();
     h_ele_met75[i] = TH1F(histname_ele_met75,"Missing Transverse Energy",11,175.,1000.);
     h_ele_met75[i].Sumw2();
     //
     h_ele_met165[i].Clear();
     h_ele_met165[i] = TH1F(histname_ele_met165,"Missing Transverse Energy",5,175.,1000.);
     h_ele_met165[i].Sumw2();
     //
     h_ele_met275[i].Clear();
     h_ele_met275[i] = TH1F(histname_ele_met275,"Missing Transverse Energy",3,175.,1000.);
     h_ele_met275[i].Sumw2();
     //



     // Muon
     h_mu_et[i].Clear();
     h_mu_et[i] = TH1F(histname_mu_et,"Photon Transverse Energy",binnumpt,binspt);
     //h_mu_et[i] = TH1F(histname_mu_et,"Photon Transverse Energy",165,175.,1000.);
     h_mu_et[i].Sumw2();
     //
     h_mu_uncorret[i].Clear();
     h_mu_uncorret[i] = TH1F(histname_mu_uncorret,"Uncorrected Photon Transverse Energy",binnumpt,binspt);
     //h_mu_uncorret[i] = TH1F(histname_mu_uncorret,"Uncorrected Photon Transverse Energy",165,175.,1000.);
     h_mu_uncorret[i].Sumw2();
     //
     h_mu_eta[i].Clear();
     h_mu_eta[i] = TH1F(histname_mu_eta,"Leading Photon Eta",10,-2.,2.);
     h_mu_eta[i].Sumw2();
     //
     h_mu_sieieF5x5[i].Clear();
     h_mu_sieieF5x5[i] = TH1F(histname_mu_sieieF5x5,"Leading Photon SigmaIetaIeta",50,0.,0.025);
     h_mu_sieieF5x5[i].Sumw2();
     //
     h_mu_pfMET[i].Clear();
     h_mu_pfMET[i] = TH1F(histname_mu_pfMET,"ParticleFlow MET",binnumpt,binspt);
     h_mu_pfMET[i].Sumw2();
     //
     h_mu_leptoMET[i].Clear();
     h_mu_leptoMET[i] = TH1F(histname_mu_leptoMET,"PF MET + dilepton",binnumpt,binspt);
     h_mu_leptoMET[i].Sumw2();
     //
     h_mu_dilep_mass[i].Clear();
     h_mu_dilep_mass[i] = TH1F(histname_mu_dilep_mass,"mass dilepton (MET)",15,60.,120.);
     h_mu_dilep_mass[i].Sumw2();
     //
     h_mu_dimu_mass[i].Clear();
     h_mu_dimu_mass[i] = TH1F(histname_mu_dimu_mass,"mass dimuon (MET)",15,60.,120.);
     h_mu_dimu_mass[i].Sumw2();
     //
     h_mu_diele_mass[i].Clear();
     h_mu_diele_mass[i] = TH1F(histname_mu_diele_mass,"mass diele (MET)",15,60.,120.);
     h_mu_diele_mass[i].Sumw2();
     //
     h_mu_et5[i].Clear();
     h_mu_et5[i] = TH1F(histname_mu_et5,"Photon Transverse Energy",165,175.,1000.);
     h_mu_et5[i].Sumw2();
     //
     h_mu_et15[i].Clear();
     h_mu_et15[i] = TH1F(histname_mu_et15,"Photon Transverse Energy",55,175.,1000.);
     h_mu_et15[i].Sumw2();
     //
     h_mu_et25[i].Clear();
     h_mu_et25[i] = TH1F(histname_mu_et25,"Photon Transverse Energy",33,175.,1000.);
     h_mu_et25[i].Sumw2();
     //
     h_mu_et55[i].Clear();
     h_mu_et55[i] = TH1F(histname_mu_et55,"Photon Transverse Energy",15,175.,1000.);
     h_mu_et55[i].Sumw2();
     //
     h_mu_et75[i].Clear();
     h_mu_et75[i] = TH1F(histname_mu_et75,"Photon Transverse Energy",11,175.,1000.);
     h_mu_et75[i].Sumw2();
     //
     h_mu_et165[i].Clear();
     h_mu_et165[i] = TH1F(histname_mu_et165,"Photon Transverse Energy",5,175.,1000.);
     h_mu_et165[i].Sumw2();
     //
     h_mu_et275[i].Clear();
     h_mu_et275[i] = TH1F(histname_mu_et275,"Photon Transverse Energy",3,175.,1000.);
     h_mu_et275[i].Sumw2();
     //
     h_mu_met5[i].Clear();
     h_mu_met5[i] = TH1F(histname_mu_met5,"Missing Transverse Energy",165,175.,1000.);
     h_mu_met5[i].Sumw2();
     //
     h_mu_met15[i].Clear();
     h_mu_met15[i] = TH1F(histname_mu_met15,"Missing Transverse Energy",55,175.,1000.);
     h_mu_met15[i].Sumw2();
     //
     h_mu_met25[i].Clear();
     h_mu_met25[i] = TH1F(histname_mu_met25,"Missing Transverse Energy",33,175.,1000.);
     h_mu_met25[i].Sumw2();
     //
     h_mu_met55[i].Clear();
     h_mu_met55[i] = TH1F(histname_mu_met55,"Missing Transverse Energy",15,175.,1000.);
     h_mu_met55[i].Sumw2();
     //
     h_mu_met75[i].Clear();
     h_mu_met75[i] = TH1F(histname_mu_met75,"Missing Transverse Energy",11,175.,1000.);
     h_mu_met75[i].Sumw2();
     //
     h_mu_met165[i].Clear();
     h_mu_met165[i] = TH1F(histname_mu_met165,"Missing Transverse Energy",5,175.,1000.);
     h_mu_met165[i].Sumw2();
     //
     h_mu_met275[i].Clear();
     h_mu_met275[i] = TH1F(histname_mu_met275,"Missing Transverse Energy",3,175.,1000.);
     h_mu_met275[i].Sumw2();
     //

   }  // ptbinnames

}

Bool_t postAnalyzer_Lep::FillSigHistogramsLep(int ptbin, int selbin, int photonIndex, double weight, bool passM){
 Float_t uncorphoet = ((*phoSCRawE)[photonIndex]/TMath::CosH((*phoSCEta)[photonIndex]));

 if(passM){
  h_mu_uncorret[ptbin]  .Fill(uncorphoet, weight );
  h_mu_et[ptbin]        .Fill( phoEt->at(photonIndex), weight );
  h_mu_eta[ptbin]       .Fill( phoEta->at(photonIndex), weight );
  h_mu_sieieF5x5[ptbin] .Fill( phoSigmaIEtaIEtaFull5x5->at(photonIndex), weight );
  h_mu_pfMET[ptbin]     .Fill( pfMET, weight );
  h_mu_leptoMET[ptbin]  .Fill( leptoMET, weight );
  h_mu_dilep_mass[ptbin].Fill( dilep_mass, weight );
  h_mu_dimu_mass[ptbin] .Fill( fourVec_mm.M(), weight );
  h_mu_diele_mass[ptbin].Fill( fourVec_ee.M(), weight );

  h_mu_et5   [ptbin].Fill(uncorphoet,weight);
  h_mu_et15  [ptbin].Fill(uncorphoet,weight);
  h_mu_et25  [ptbin].Fill(uncorphoet,weight);
  h_mu_et55  [ptbin].Fill(uncorphoet,weight);
  h_mu_et75  [ptbin].Fill(uncorphoet,weight);
  h_mu_et165 [ptbin].Fill(uncorphoet,weight);
  h_mu_et275 [ptbin].Fill(uncorphoet,weight);
  h_mu_met5  [ptbin].Fill(pfMET,weight);
  h_mu_met15 [ptbin].Fill(pfMET,weight);
  h_mu_met25 [ptbin].Fill(pfMET,weight);
  h_mu_met55 [ptbin].Fill(pfMET,weight);
  h_mu_met75 [ptbin].Fill(pfMET,weight);
  h_mu_met165[ptbin].Fill(pfMET,weight);
  h_mu_met275[ptbin].Fill(pfMET,weight);
 }
 else{
  h_ele_uncorret[ptbin]  .Fill(uncorphoet, weight );
  h_ele_et[ptbin]        .Fill( phoEt->at(photonIndex), weight );
  h_ele_eta[ptbin]       .Fill( phoEta->at(photonIndex), weight );
  h_ele_sieieF5x5[ptbin] .Fill( phoSigmaIEtaIEtaFull5x5->at(photonIndex), weight );
  h_ele_pfMET[ptbin]     .Fill( pfMET, weight );
  h_ele_leptoMET[ptbin]  .Fill( leptoMET, weight );
  h_ele_dilep_mass[ptbin].Fill( dilep_mass, weight );
  h_ele_dimu_mass[ptbin] .Fill( fourVec_mm.M(), weight );
  h_ele_diele_mass[ptbin].Fill( fourVec_ee.M(), weight );

  h_ele_et5   [ptbin].Fill(uncorphoet,weight);
  h_ele_et15  [ptbin].Fill(uncorphoet,weight);
  h_ele_et25  [ptbin].Fill(uncorphoet,weight);
  h_ele_et55  [ptbin].Fill(uncorphoet,weight);
  h_ele_et75  [ptbin].Fill(uncorphoet,weight);
  h_ele_et165 [ptbin].Fill(uncorphoet,weight);
  h_ele_et275 [ptbin].Fill(uncorphoet,weight);
  h_ele_met5  [ptbin].Fill(pfMET,weight);
  h_ele_met15 [ptbin].Fill(pfMET,weight);
  h_ele_met25 [ptbin].Fill(pfMET,weight);
  h_ele_met55 [ptbin].Fill(pfMET,weight);
  h_ele_met75 [ptbin].Fill(pfMET,weight);
  h_ele_met165[ptbin].Fill(pfMET,weight);
  h_ele_met275[ptbin].Fill(pfMET,weight);
 }

 return kTRUE;
}

Bool_t postAnalyzer_Lep::WriteHistogramsLep(int ptbin, int sysbin){

  h_ele_uncorret[ptbin]  .Write();
  h_ele_et[ptbin]        .Write();
  h_ele_eta[ptbin]       .Write();
  h_ele_sieieF5x5[ptbin] .Write();
  h_ele_pfMET[ptbin]     .Write();
  h_ele_leptoMET[ptbin]  .Write();
  h_ele_dilep_mass[ptbin].Write();
  h_ele_dimu_mass[ptbin] .Write();
  h_ele_diele_mass[ptbin].Write();
  h_ele_et5[ptbin]       .Write();
  h_ele_et15[ptbin]      .Write();
  h_ele_et25[ptbin]      .Write();
  h_ele_et55[ptbin]      .Write();
  h_ele_et75[ptbin]      .Write();
  h_ele_et165[ptbin]     .Write();
  h_ele_et275[ptbin]     .Write();
  h_ele_met5[ptbin]      .Write();
  h_ele_met15[ptbin]     .Write();
  h_ele_met25[ptbin]     .Write();
  h_ele_met55[ptbin]     .Write();
  h_ele_met75[ptbin]     .Write();
  h_ele_met165[ptbin]    .Write();
  h_ele_met275[ptbin]    .Write();

  h_mu_uncorret[ptbin]  .Write();
  h_mu_et[ptbin]        .Write();
  h_mu_eta[ptbin]       .Write();
  h_mu_sieieF5x5[ptbin] .Write();
  h_mu_pfMET[ptbin]     .Write();
  h_mu_leptoMET[ptbin]  .Write();
  h_mu_dilep_mass[ptbin].Write();
  h_mu_dimu_mass[ptbin] .Write();
  h_mu_diele_mass[ptbin].Write();
  h_mu_et5[ptbin]       .Write();
  h_mu_et15[ptbin]      .Write();
  h_mu_et25[ptbin]      .Write();
  h_mu_et55[ptbin]      .Write();
  h_mu_et75[ptbin]      .Write();
  h_mu_et165[ptbin]     .Write();
  h_mu_et275[ptbin]     .Write();
  h_mu_met5[ptbin]      .Write();
  h_mu_met15[ptbin]     .Write();
  h_mu_met25[ptbin]     .Write();
  h_mu_met55[ptbin]     .Write();
  h_mu_met75[ptbin]     .Write();
  h_mu_met165[ptbin]    .Write();
  h_mu_met275[ptbin]    .Write();

 return kTRUE;

}

//-------------------------callFillSigHistLep
void postAnalyzer_Lep::callFillSigHistLep(int selbin, int lastptbin, int inclptbin, int candphotonindex, float event_weight, bool passM){
 Float_t uncorrectedPhoEt = ((*phoSCRawE)[candphotonindex]/TMath::CosH((*phoSCEta)[candphotonindex]));
 for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
  if(
     ( uncorrectedPhoEt > ptbins[ptb]) &&
     ( uncorrectedPhoEt < ptbins[ptb+1])
    ){
   FillSigHistogramsLep(ptb, selbin, candphotonindex, event_weight, passM);
  } // end if passes pt cuts then fill
 } // end pt bin loop
 if(  // do an inclusive pT plot from bins
    ( uncorrectedPhoEt > ptbins[0]) &&
    ( uncorrectedPhoEt < ptbins[inclptbin])
   ){
  FillSigHistogramsLep(lastptbin-1, selbin, candphotonindex, event_weight, passM);
 }
 // and one fully inclusive in pT
 FillSigHistogramsLep(lastptbin, selbin, candphotonindex, event_weight, passM);
 return;
}


#endif
