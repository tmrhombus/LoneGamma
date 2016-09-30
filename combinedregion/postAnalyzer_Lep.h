
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
        h_ele_dimu_mass[7];

   TH1F h_mu_et[7],
        h_mu_uncorret[7],
        h_mu_eta[7],
        h_mu_sieieF5x5[7],
        h_mu_pfMET[7],
        h_mu_leptoMET[7],
        h_mu_dilep_mass[7],
        h_mu_diele_mass[7],
        h_mu_dimu_mass[7];

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
     TString histname_mu_et  = "h_mu_et_"+ptbinnames[i];
     TString histname_mu_uncorret  = "h_mu_uncorret_"+ptbinnames[i];
     TString histname_mu_eta = "h_mu_eta_"+ptbinnames[i];
     TString histname_mu_sieieF5x5 = "h_mu_sieieF5x5_"+ptbinnames[i];
     TString histname_mu_pfMET = "h_mu_pfMET_"+ptbinnames[i];
     TString histname_mu_leptoMET = "h_mu_leptoMET_"+ptbinnames[i];
     TString histname_mu_dilep_mass = "h_mu_dilep_mass_"+ptbinnames[i];
     TString histname_mu_diele_mass = "h_mu_diele_mass_"+ptbinnames[i];
     TString histname_mu_dimu_mass = "h_mu_dimu_mass_"+ptbinnames[i];

     TString histname_ele_et  = "h_ele_et_"+ptbinnames[i];
     TString histname_ele_uncorret  = "h_ele_uncorret_"+ptbinnames[i];
     TString histname_ele_eta = "h_ele_eta_"+ptbinnames[i];
     TString histname_ele_sieieF5x5 = "h_ele_sieieF5x5_"+ptbinnames[i];
     TString histname_ele_pfMET = "h_ele_pfMET_"+ptbinnames[i];
     TString histname_ele_leptoMET = "h_ele_leptoMET_"+ptbinnames[i];
     TString histname_ele_dilep_mass = "h_ele_dilep_mass_"+ptbinnames[i];
     TString histname_ele_diele_mass = "h_ele_diele_mass_"+ptbinnames[i];
     TString histname_ele_dimu_mass = "h_ele_dimu_mass_"+ptbinnames[i];

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

   }  // ptbinnames

}

Bool_t postAnalyzer_Lep::FillSigHistogramsLep(int ptbin, int selbin, int photonIndex, double weight, bool passM){
 Float_t uncorphoet = ((*phoSCRawE)[photonIndex]/TMath::CosH((*phoSCEta)[photonIndex]));

 if(passM){
  h_mu_uncorret[ptbin].Fill(uncorphoet, weight );
  h_mu_et[ptbin].Fill( phoEt->at(photonIndex), weight );
  h_mu_eta[ptbin].Fill( phoEta->at(photonIndex), weight );
  h_mu_sieieF5x5[ptbin].Fill( phoSigmaIEtaIEtaFull5x5->at(photonIndex), weight );
  h_mu_pfMET[ptbin].Fill( pfMET, weight );
  h_mu_leptoMET[ptbin].Fill( leptoMET, weight );
  h_mu_dilep_mass[ptbin].Fill( dilep_mass, weight );
  h_mu_dimu_mass[ptbin].Fill( fourVec_mm.M(), weight );
  h_mu_diele_mass[ptbin].Fill( fourVec_ee.M(), weight );
 }
 else{
  h_ele_uncorret[ptbin].Fill(uncorphoet, weight );
  h_ele_et[ptbin].Fill( phoEt->at(photonIndex), weight );
  h_ele_eta[ptbin].Fill( phoEta->at(photonIndex), weight );
  h_ele_sieieF5x5[ptbin].Fill( phoSigmaIEtaIEtaFull5x5->at(photonIndex), weight );
  h_ele_pfMET[ptbin].Fill( pfMET, weight );
  h_ele_leptoMET[ptbin].Fill( leptoMET, weight );
  h_ele_dilep_mass[ptbin].Fill( dilep_mass, weight );
  h_ele_dimu_mass[ptbin].Fill( fourVec_mm.M(), weight );
  h_ele_diele_mass[ptbin].Fill( fourVec_ee.M(), weight );
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

  h_mu_uncorret[ptbin]  .Write();
  h_mu_et[ptbin]        .Write();
  h_mu_eta[ptbin]       .Write();
  h_mu_sieieF5x5[ptbin] .Write();
  h_mu_pfMET[ptbin]     .Write();
  h_mu_leptoMET[ptbin]  .Write();
  h_mu_dilep_mass[ptbin].Write();
  h_mu_dimu_mass[ptbin] .Write();
  h_mu_diele_mass[ptbin].Write();

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
