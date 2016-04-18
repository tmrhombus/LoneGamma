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

 int ntot=0;

 int n_mch_all = 0;
 int n_mch_met = 0;
 int n_mch_trg = 0;
 int n_mch_met_trg = 0;
 int n_mch_t175 = 0;
 int n_mch_met_t175 = 0;
 int n_mch_t250 = 0;
 int n_mch_met_t250 = 0;
 int n_mch_allt = 0;
 int n_mch_met_allt = 0;

 int n_nmh_all = 0;
 int n_nmh_met = 0;
 int n_nmh_trg = 0;
 int n_nmh_met_trg = 0;
 int n_nmh_t175 = 0;
 int n_nmh_met_t175 = 0;
 int n_nmh_t250 = 0;
 int n_nmh_met_t250 = 0;
 int n_nmh_allt = 0;
 int n_nmh_met_allt = 0;

 if (fChain == 0) return;
 //std::cout<<"Numerator/Denominator : run : lumis : event : photonindex "<<std::endl;

 printf(" a: Photon ID only (pT > 175) \n");
 printf(" b: ID + MET<30 \n");
 printf(" c: ID + T165_he10 \n");
 printf(" d: ID + MET + T165 \n");
 printf(" e: ID + T175 \n");
 printf(" f: ID + MET + T175 \n");
 printf(" g: ID + T250 \n");
 printf(" h: ID + MET + T250 \n");
 printf(" i: ID + MET + All Triggers \n");

 Long64_t nentries = fChain->GetEntriesFast();

 Long64_t nbytes = 0, nb = 0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
  ntot++;
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
   genPCvint.clear(); 

   // standard (no met or trigger req)
   recPCvint.clear(); 
   mchPCvint.clear(); 
   nmhPCvint.clear(); 
    // denominator
   recDenPCvint.clear(); 
   mchDenPCvint.clear(); 
   nmhDenPCvint.clear(); 

   // Find all reco photons passing selections
   recPCvint = pcPassSel(175.,2000.,1.4442);
   recDenPCvint = pcPassDenSel(175.,2000.,1.4442);

    // Denominator

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

   //Loop over all    rec    photons to find (GEN)matched and nonmatched photons
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


   //Loop over all    den    photons to find (GEN)matched and nonmatched photons
   for(int i = 0; i < recDenPCvint.size(); i++)
   {
    int phonr = recDenPCvint[i];
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
    if(dRmin < 0.1) { mchDenPCvint.push_back(phonr); }
    if(dRmin > 0.1) { nmhDenPCvint.push_back(phonr); }
   }

   //   // lepton rejection
   //   std::vector<int> elelist = electron_passLooseID(candphotonindex, 10.);
   //   std::vector<int> mulist = muon_passLooseID(candphotonindex, 10.);
   //   bool passLepRej = false;
   //   if( elelist.size()==0 && mulist.size()==0 ){ passLepRej = true; }

      // MET requirements
      bool passMETfilters = kTRUE; // ( metFilters==0 ) ; //|| metFilters==4 ) ;

   //   // dPhi( Jets, MET )
   //   std::vector<int>  jetindexvector = selectedJets(candphotonindex);
   //   bool passdPhiJM = passdphiJetMET(&jetindexvector, pfMETPhi);

   //   // dPhi( photon, MET )
   //   bool passdPhiPhoMET = ( DeltaPhi(phoPhi->at(candphotonindex),pfMETPhi)>2.0 ) ;

   //   // spike cleaning
   //   int iphi = 41;
   //   int ieta = 5;
   //   bool passSpike = !(phoIPhi->at(candphotonindex) == iphi && phoIEta->at(candphotonindex) == ieta) ;

   // some more filters
   bool passTriggers;
   bool passTriggers175;
   bool passTriggers250;
   bool passMET;
   passMET = pfMET < 30. && passMETfilters ;

// fill histograms
    // Fill Numerator (matched) Signal Histograms
    if( mchPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<mchPCvint.size(); ++k){ // go through each photon passing sel
      int phonr = mchPCvint.at(k);
      passTriggers = ((HLTPho>>12&1) == 1   && ( (phoFiredSingleTrgs->at(phonr)>>7)!=1 ) ) ;
      passTriggers175 = ((HLTPho>>7&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>8)!=1 ) ) ;
      passTriggers250 = ((HLTPho>>8&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>9)!=1 ) ) ;
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(phonr) > ptbins[ptb]) &&
          (phoEt->at(phonr) < ptbins[ptb+1])
         ){
            //std::cout<<"found event"<<std::endl;
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillSigHistograms(ptb, 0, phonr, event_weight);
            FillDenHistograms(ptb, 0, phonr, event_weight);
            if(passMET){
             FillSigHistograms(ptb, 1, phonr, event_weight);
             FillDenHistograms(ptb, 1, phonr, event_weight);
            }
            if(passTriggers){
             FillSigHistograms(ptb, 2, phonr, event_weight);
             FillDenHistograms(ptb, 2, phonr, event_weight);
            }
            if(passTriggers && passMET){
             FillSigHistograms(ptb, 3, phonr, event_weight);
             FillDenHistograms(ptb, 3, phonr, event_weight);
            }
            if(passTriggers175){
             FillSigHistograms(ptb, 4, phonr, event_weight);
             FillDenHistograms(ptb, 4, phonr, event_weight);
            }
            if(passTriggers175 && passMET){
             FillSigHistograms(ptb, 5, phonr, event_weight);
             FillDenHistograms(ptb, 5, phonr, event_weight);
            }
            if(passTriggers250){
             FillSigHistograms(ptb, 6, phonr, event_weight);
             FillDenHistograms(ptb, 6, phonr, event_weight);
            }
            if(passTriggers250 && passMET){
             FillSigHistograms(ptb, 7, phonr, event_weight);
             FillDenHistograms(ptb, 7, phonr, event_weight);
            }
            if(passTriggers && passTriggers175 && passTriggers250 && passMET){
             FillSigHistograms(ptb, 8, phonr, event_weight);
             FillDenHistograms(ptb, 8, phonr, event_weight);
            }
           }
       }
      // inclusive from pT bins

      if(  
         (phoEt->at(phonr) > 175. ) && //ptbins[0]) &&
         (phoEt->at(phonr) < ptbins[ptbins.size()-1])
        ){
printf("Matched    a  ---  %i:%i:%lli \n", run, lumis, event );
printf("  nVtx                %i  \n", nVtx );
printf("  rho                 %f  \n", rho );
printf("  metFilters          %i  \n", metFilters );
printf("  nPho                %i  \n", nPho );
printf("  phoE                %f  \n", phoE->at(phonr) );
printf("  phoEt               %f  \n", phoEt->at(phonr) );
printf("  phoPhi              %f  \n", phoPhi->at(phonr) );
printf("  sieie               %f  \n", phoSigmaIEtaIEtaFull5x5->at(phonr) );
printf("  sieip               %f  \n", phoSigmaIEtaIPhiFull5x5->at(phonr) );
printf("  sipip               %f  \n", phoSigmaIPhiIPhiFull5x5->at(phonr) );
printf("  phoSCEta            %f  \n", phoSCEta->at(phonr) );
printf("  phoHoverE           %f  \n", phoHoverE->at(phonr) );
printf("  phohasPixelSeed     %i  \n", phohasPixelSeed->at(phonr) );
printf("  phosPFNeuIso        %f  \n", phoPFNeuIso->at(phonr) );
printf("  phoPFPhoIso         %f  \n", phoPFPhoIso->at(phonr) );
printf("  phoPFChIso          %f  \n", phoPFChIso->at(phonr) );
printf("  phoPFChWorstIso     %f  \n", phoPFChWorstIso->at(phonr) );

 n_mch_all++;
          FillSigHistograms(selptbin, 0, phonr, event_weight);
          FillDenHistograms(selptbin, 0, phonr, event_weight);
          if(passMET){
           FillSigHistograms(selptbin, 1, phonr, event_weight);
           FillDenHistograms(selptbin, 1, phonr, event_weight);
printf("Matched    b  ---  %i:%i:%lli \n", run, lumis, event );
 n_mch_met++;
          }
          if(passTriggers){
           FillSigHistograms(selptbin, 2, phonr, event_weight);
           FillDenHistograms(selptbin, 2, phonr, event_weight);
printf("Matched    c  ---  %i:%i:%lli \n", run, lumis, event );
 n_mch_trg++;
          }
          if(passTriggers && passMET){
           FillSigHistograms(selptbin, 3, phonr, event_weight);
           FillDenHistograms(selptbin, 3, phonr, event_weight);
printf("Matched    d  ---  %i:%i:%lli \n", run, lumis, event );
 n_mch_met_trg++;
          }
          if(passTriggers175){
           FillSigHistograms(selptbin, 4, phonr, event_weight);
           FillDenHistograms(selptbin, 4, phonr, event_weight);
printf("Matched    e  ---  %i:%i:%lli \n", run, lumis, event );
          }
          if(passTriggers175 && passMET){
           FillSigHistograms(selptbin, 5, phonr, event_weight);
           FillDenHistograms(selptbin, 5, phonr, event_weight);
printf("Matched    f  ---  %i:%i:%lli \n", run, lumis, event );
          }
          if(passTriggers250){
           FillSigHistograms(selptbin, 6, phonr, event_weight);
           FillDenHistograms(selptbin, 6, phonr, event_weight);
printf("Matched    g  ---  %i:%i:%lli \n", run, lumis, event );
          }
          if(passTriggers250 && passMET){
           FillSigHistograms(selptbin, 7, phonr, event_weight);
           FillDenHistograms(selptbin, 7, phonr, event_weight);
printf("Matched    h  ---  %i:%i:%lli \n", run, lumis, event );
          }
          if(passTriggers && passTriggers175 && passTriggers250 && passMET){
           FillSigHistograms(selptbin, 8, phonr, event_weight);
           FillDenHistograms(selptbin, 8, phonr, event_weight);
printf("Matched    i  ---  %i:%i:%lli \n", run, lumis, event );
          }
         }

      // fully inclusive in pT (depending only on selections for mchPCvint)
      FillSigHistograms(incptbin, 0, phonr, event_weight);
      FillDenHistograms(incptbin, 0, phonr, event_weight);
      if(passMET){
       FillSigHistograms(incptbin, 1, phonr, event_weight);
       FillDenHistograms(incptbin, 1, phonr, event_weight);
      }
      if(passTriggers){
       FillSigHistograms(incptbin, 2, phonr, event_weight);
       FillDenHistograms(incptbin, 2, phonr, event_weight);
      }
      if(passTriggers && passMET){
       FillSigHistograms(incptbin, 3, phonr, event_weight);
       FillDenHistograms(incptbin, 3, phonr, event_weight);
      }
      if(passTriggers175){
       FillSigHistograms(incptbin, 4, phonr, event_weight);
       FillDenHistograms(incptbin, 4, phonr, event_weight);
      }
      if(passTriggers175 && passMET){
       FillSigHistograms(incptbin, 5, phonr, event_weight);
       FillDenHistograms(incptbin, 5, phonr, event_weight);
      }
      if(passTriggers250){
       FillSigHistograms(incptbin, 6, phonr, event_weight);
       FillDenHistograms(incptbin, 6, phonr, event_weight);
      }
      if(passTriggers250 && passMET){
       FillSigHistograms(incptbin, 7, phonr, event_weight);
       FillDenHistograms(incptbin, 7, phonr, event_weight);
      }
      if(passTriggers && passTriggers175 && passTriggers250 && passMET){
       FillSigHistograms(incptbin, 8, phonr, event_weight);
       FillDenHistograms(incptbin, 8, phonr, event_weight);
      }

     }
    }

    // Fill Numerator (unmatched) Bkgnal Histograms
    if( nmhPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<nmhPCvint.size(); ++k){ // go through each photon passing sel
      int phonr = nmhPCvint.at(k);
      passTriggers = ((HLTPho>>12&1) == 1   && ( (phoFiredSingleTrgs->at(phonr)>>7)!=1 ) ) ;
      passTriggers175 = ((HLTPho>>7&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>8)!=1 ) ) ;
      passTriggers250 = ((HLTPho>>8&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>9)!=1 ) ) ;
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(phonr) > ptbins[ptb]) &&
          (phoEt->at(phonr) < ptbins[ptb+1])
         ){
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillBkgHistograms(ptb, 0, phonr, event_weight);
            FillDenHistograms(ptb, 0, phonr, event_weight);
            if(passMET){
             FillBkgHistograms(ptb, 1, phonr, event_weight);
             FillDenHistograms(ptb, 1, phonr, event_weight);
            }
            if(passTriggers){
             FillBkgHistograms(ptb, 2, phonr, event_weight);
             FillDenHistograms(ptb, 2, phonr, event_weight);
            }
            if(passTriggers && passMET){
             FillBkgHistograms(ptb, 3, phonr, event_weight);
             FillDenHistograms(ptb, 3, phonr, event_weight);
            }
            if(passTriggers175){
             FillSigHistograms(ptb, 4, phonr, event_weight);
             FillDenHistograms(ptb, 4, phonr, event_weight);
            }
            if(passTriggers175 && passMET){
             FillSigHistograms(ptb, 5, phonr, event_weight);
             FillDenHistograms(ptb, 5, phonr, event_weight);
            }
            if(passTriggers250){
             FillSigHistograms(ptb, 6, phonr, event_weight);
             FillDenHistograms(ptb, 6, phonr, event_weight);
            }
            if(passTriggers250 && passMET){
             FillSigHistograms(ptb, 7, phonr, event_weight);
             FillDenHistograms(ptb, 7, phonr, event_weight);
            }
            if(passTriggers && passTriggers175 && passTriggers250 && passMET){
             FillSigHistograms(ptb, 8, phonr, event_weight);
             FillDenHistograms(ptb, 8, phonr, event_weight);
            }
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(phonr) > 175. ) && // ptbins[0]) &&
         (phoEt->at(phonr) < ptbins[ptbins.size()-1])
        ){
n_nmh_all++;
printf("Unmatched  a  ---  %i:%i:%lli \n", run, lumis, event );
printf("  nVtx                %i  \n", nVtx );
printf("  rho                 %f  \n", rho );
printf("  metFilters          %i  \n", metFilters );
printf("  nPho                %i  \n", nPho );
printf("  phoE                %f  \n", phoE->at(phonr) );
printf("  phoEt               %f  \n", phoEt->at(phonr) );
printf("  phoPhi              %f  \n", phoPhi->at(phonr) );
printf("  sieie               %f  \n", phoSigmaIEtaIEtaFull5x5->at(phonr) );
printf("  sieip               %f  \n", phoSigmaIEtaIPhiFull5x5->at(phonr) );
printf("  sipip               %f  \n", phoSigmaIPhiIPhiFull5x5->at(phonr) );
printf("  phoSCEta            %f  \n", phoSCEta->at(phonr) );
printf("  phoHoverE           %f  \n", phoHoverE->at(phonr) );
printf("  phohasPixelSeed     %i  \n", phohasPixelSeed->at(phonr) );
printf("  phoPFNeuIso         %f  \n", phoPFNeuIso->at(phonr) );
printf("  phoPFPhoIso         %f  \n", phoPFPhoIso->at(phonr) );
printf("  phoPFChIso          %f  \n", phoPFChIso->at(phonr) );
printf("  phoPFChWorstIso     %f  \n", phoPFChWorstIso->at(phonr) );
          FillBkgHistograms(selptbin, 0, phonr, event_weight);
          FillDenHistograms(selptbin, 0, phonr, event_weight);
          if(passMET){
printf("Unmatched  b  ---  %i:%i:%lli \n", run, lumis, event );
           FillBkgHistograms(selptbin, 1, phonr, event_weight);
           FillDenHistograms(selptbin, 1, phonr, event_weight);
n_nmh_met++;
          }
          if(passTriggers){
printf("Unmatched  c  ---  %i:%i:%lli \n", run, lumis, event );
           FillBkgHistograms(selptbin, 2, phonr, event_weight);
           FillDenHistograms(selptbin, 2, phonr, event_weight);
n_nmh_trg++;
          }
          if(passTriggers && passMET){
printf("Unmatched  d  ---  %i:%i:%lli \n", run, lumis, event );
           FillBkgHistograms(selptbin, 3, phonr, event_weight);
           FillDenHistograms(selptbin, 3, phonr, event_weight);
n_nmh_met_trg++;
          if(passTriggers175){
printf("Unmatched  e  ---  %i:%i:%lli \n", run, lumis, event );
           FillSigHistograms(selptbin, 4, phonr, event_weight);
           FillDenHistograms(selptbin, 4, phonr, event_weight);
          }
          if(passTriggers175 && passMET){
printf("Unmatched  f  ---  %i:%i:%lli \n", run, lumis, event );
           FillSigHistograms(selptbin, 5, phonr, event_weight);
           FillDenHistograms(selptbin, 5, phonr, event_weight);
          }
          if(passTriggers250){
printf("Unmatched  g  ---  %i:%i:%lli \n", run, lumis, event );
           FillSigHistograms(selptbin, 6, phonr, event_weight);
           FillDenHistograms(selptbin, 6, phonr, event_weight);
          }
          if(passTriggers250 && passMET){
printf("Unmatched  h  ---  %i:%i:%lli \n", run, lumis, event );
           FillSigHistograms(selptbin, 7, phonr, event_weight);
           FillDenHistograms(selptbin, 7, phonr, event_weight);
          }
          if(passTriggers && passTriggers175 && passTriggers250 && passMET){
printf("Unmatched  i  ---  %i:%i:%lli \n", run, lumis, event );
           FillSigHistograms(selptbin, 8, phonr, event_weight);
           FillDenHistograms(selptbin, 8, phonr, event_weight);
          }

          }
         }

      // fully inclusive in pT (depending only on selections for nmhPCvint)
      FillBkgHistograms(incptbin, 0, phonr, event_weight);
      FillDenHistograms(incptbin, 0, phonr, event_weight);
      if(passMET){
       FillBkgHistograms(incptbin, 1, phonr, event_weight);
       FillDenHistograms(incptbin, 1, phonr, event_weight);
      }
      if(passTriggers){
       FillBkgHistograms(incptbin, 2, phonr, event_weight);
       FillDenHistograms(incptbin, 2, phonr, event_weight);
      }
      if(passTriggers && passMET){
       FillBkgHistograms(incptbin, 3, phonr, event_weight);
       FillDenHistograms(incptbin, 3, phonr, event_weight);
      }
      if(passTriggers175){
       FillSigHistograms(incptbin, 4, phonr, event_weight);
       FillDenHistograms(incptbin, 4, phonr, event_weight);
      }
      if(passTriggers175 && passMET){
       FillSigHistograms(incptbin, 5, phonr, event_weight);
       FillDenHistograms(incptbin, 5, phonr, event_weight);
      }
      if(passTriggers250){
       FillSigHistograms(incptbin, 6, phonr, event_weight);
       FillDenHistograms(incptbin, 6, phonr, event_weight);
      }
      if(passTriggers250 && passMET){
       FillSigHistograms(incptbin, 7, phonr, event_weight);
       FillDenHistograms(incptbin, 7, phonr, event_weight);
      }
      if(passTriggers && passTriggers175 && passTriggers250 && passMET){
       FillSigHistograms(incptbin, 8, phonr, event_weight);
       FillDenHistograms(incptbin, 8, phonr, event_weight);
      }

     }
    }

    /////////////////////////////// Denominator selections
    // Fill Numerator (matched) Denom Histograms
    if( mchDenPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<mchDenPCvint.size(); ++k){ // go through each photon passing sel
      int phonr = mchDenPCvint.at(k);
      passTriggers = ((HLTPho>>12&1) == 1   && ( (phoFiredSingleTrgs->at(phonr)>>7)!=1 ) ) ;
      passTriggers175 = ((HLTPho>>7&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>8)!=1 ) ) ;
      passTriggers250 = ((HLTPho>>8&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>9)!=1 ) ) ;
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(phonr) > ptbins[ptb]) &&
          (phoEt->at(phonr) < ptbins[ptb+1])
         ){
            //std::cout<<"found event"<<std::endl;
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillDenSigHistograms(ptb, 0, phonr, event_weight);
            FillDenDenHistograms(ptb, 0, phonr, event_weight);
            if(passMET){
             FillDenSigHistograms(ptb, 1, phonr, event_weight);
             FillDenDenHistograms(ptb, 1, phonr, event_weight);
            }
            if(passTriggers){
             FillDenSigHistograms(ptb, 2, phonr, event_weight);
             FillDenDenHistograms(ptb, 2, phonr, event_weight);
            }
            if(passTriggers && passMET){
             FillDenSigHistograms(ptb, 3, phonr, event_weight);
             FillDenDenHistograms(ptb, 3, phonr, event_weight);
            }
            if(passTriggers175){
             FillDenSigHistograms(ptb, 4, phonr, event_weight);
             FillDenDenHistograms(ptb, 4, phonr, event_weight);
            }
            if(passTriggers175 && passMET){
             FillDenSigHistograms(ptb, 5, phonr, event_weight);
             FillDenDenHistograms(ptb, 5, phonr, event_weight);
            }
            if(passTriggers250){
             FillDenSigHistograms(ptb, 6, phonr, event_weight);
             FillDenDenHistograms(ptb, 6, phonr, event_weight);
            }
            if(passTriggers250 && passMET){
             FillDenSigHistograms(ptb, 7, phonr, event_weight);
             FillDenDenHistograms(ptb, 7, phonr, event_weight);
            }
            if(passTriggers && passTriggers175 && passTriggers250 && passMET){
             FillDenSigHistograms(ptb, 8, phonr, event_weight);
             FillDenDenHistograms(ptb, 8, phonr, event_weight);
            }
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(phonr) > 175. ) && //ptbins[0]) &&
         (phoEt->at(phonr) < ptbins[ptbins.size()-1])
        ){
          FillDenSigHistograms(selptbin, 0, phonr, event_weight);
          FillDenDenHistograms(selptbin, 0, phonr, event_weight);
          if(passMET){
           FillDenSigHistograms(selptbin, 1, phonr, event_weight);
           FillDenDenHistograms(selptbin, 1, phonr, event_weight);
          }
          if(passTriggers){
           FillDenSigHistograms(selptbin, 2, phonr, event_weight);
           FillDenDenHistograms(selptbin, 2, phonr, event_weight);
          }
          if(passTriggers && passMET){
           FillDenSigHistograms(selptbin, 3, phonr, event_weight);
           FillDenDenHistograms(selptbin, 3, phonr, event_weight);
          }
            if(passTriggers175){
             FillDenSigHistograms(selptbin, 4, phonr, event_weight);
             FillDenDenHistograms(selptbin, 4, phonr, event_weight);
            }
            if(passTriggers175 && passMET){
             FillDenSigHistograms(selptbin, 5, phonr, event_weight);
             FillDenDenHistograms(selptbin, 5, phonr, event_weight);
            }
            if(passTriggers250){
             FillDenSigHistograms(selptbin, 6, phonr, event_weight);
             FillDenDenHistograms(selptbin, 6, phonr, event_weight);
            }
            if(passTriggers250 && passMET){
             FillDenSigHistograms(selptbin, 7, phonr, event_weight);
             FillDenDenHistograms(selptbin, 7, phonr, event_weight);
            }
            if(passTriggers && passTriggers175 && passTriggers250 && passMET){
             FillDenSigHistograms(selptbin, 8, phonr, event_weight);
             FillDenDenHistograms(selptbin, 8, phonr, event_weight);
            }
         }

      // fully inclusive in pT (depending only on selections for mchDenPCvint)
      FillDenSigHistograms(incptbin, 0, phonr, event_weight);
      FillDenDenHistograms(incptbin, 0, phonr, event_weight);
      if(passMET){
       FillDenSigHistograms(incptbin, 1, phonr, event_weight);
       FillDenDenHistograms(incptbin, 1, phonr, event_weight);
      }
      if(passTriggers){
       FillDenSigHistograms(incptbin, 2, phonr, event_weight);
       FillDenDenHistograms(incptbin, 2, phonr, event_weight);
      }
      if(passTriggers && passMET){
       FillDenSigHistograms(incptbin, 3, phonr, event_weight);
       FillDenDenHistograms(incptbin, 3, phonr, event_weight);
      }
            if(passTriggers175){
             FillDenSigHistograms(incptbin, 4, phonr, event_weight);
             FillDenDenHistograms(incptbin, 4, phonr, event_weight);
            }
            if(passTriggers175 && passMET){
             FillDenSigHistograms(incptbin, 5, phonr, event_weight);
             FillDenDenHistograms(incptbin, 5, phonr, event_weight);
            }
            if(passTriggers250){
             FillDenSigHistograms(incptbin, 6, phonr, event_weight);
             FillDenDenHistograms(incptbin, 6, phonr, event_weight);
            }
            if(passTriggers250 && passMET){
             FillDenSigHistograms(incptbin, 7, phonr, event_weight);
             FillDenDenHistograms(incptbin, 7, phonr, event_weight);
            }
            if(passTriggers && passTriggers175 && passTriggers250 && passMET){
             FillDenSigHistograms(incptbin, 8, phonr, event_weight);
             FillDenDenHistograms(incptbin, 8, phonr, event_weight);
            }

     }
    }


    // FillDen Numerator (unmatched) Den Histograms
    if( nmhDenPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<nmhDenPCvint.size(); ++k){ // go through each photon passing sel
      int phonr = nmhDenPCvint.at(k);
      passTriggers = ((HLTPho>>12&1) == 1   && ( (phoFiredSingleTrgs->at(phonr)>>7)!=1 ) ) ;
      passTriggers175 = ((HLTPho>>7&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>8)!=1 ) ) ;
      passTriggers250 = ((HLTPho>>8&1) == 1 && ( (phoFiredSingleTrgs->at(phonr)>>9)!=1 ) ) ;
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(phonr) > ptbins[ptb]) &&
          (phoEt->at(phonr) < ptbins[ptb+1])
         ){
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillDenBkgHistograms(ptb, 0, phonr, event_weight);
            FillDenDenHistograms(ptb, 0, phonr, event_weight);
            if(passMET){
             FillDenBkgHistograms(ptb, 1, phonr, event_weight);
             FillDenDenHistograms(ptb, 1, phonr, event_weight);
            }
            if(passTriggers){
             FillDenBkgHistograms(ptb, 2, phonr, event_weight);
             FillDenDenHistograms(ptb, 2, phonr, event_weight);
            }
            if(passTriggers && passMET){
             FillDenBkgHistograms(ptb, 3, phonr, event_weight);
             FillDenDenHistograms(ptb, 3, phonr, event_weight);
            }
            if(passTriggers175){
             FillDenBkgHistograms(ptb, 4, phonr, event_weight);
             FillDenDenHistograms(ptb, 4, phonr, event_weight);
            }
            if(passTriggers175 && passMET){
             FillDenBkgHistograms(ptb, 5, phonr, event_weight);
             FillDenDenHistograms(ptb, 5, phonr, event_weight);
            }
            if(passTriggers250){
             FillDenBkgHistograms(ptb, 6, phonr, event_weight);
             FillDenDenHistograms(ptb, 6, phonr, event_weight);
            }
            if(passTriggers250 && passMET){
             FillDenBkgHistograms(ptb, 7, phonr, event_weight);
             FillDenDenHistograms(ptb, 7, phonr, event_weight);
            }
            if(passTriggers && passTriggers175 && passTriggers250 && passMET){
             FillDenBkgHistograms(ptb, 8, phonr, event_weight);
             FillDenDenHistograms(ptb, 8, phonr, event_weight);
            }

           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(phonr) > 175. ) && // ptbins[0]) &&
         (phoEt->at(phonr) < ptbins[ptbins.size()-1])
        ){
          FillDenBkgHistograms(selptbin, 0, phonr, event_weight);
          FillDenDenHistograms(selptbin, 0, phonr, event_weight);
          if(passMET){
           FillDenBkgHistograms(selptbin, 1, phonr, event_weight);
           FillDenDenHistograms(selptbin, 1, phonr, event_weight);
          }
          if(passTriggers){
           FillDenBkgHistograms(selptbin, 2, phonr, event_weight);
           FillDenDenHistograms(selptbin, 2, phonr, event_weight);
          }
          if(passTriggers && passMET){
           FillDenBkgHistograms(selptbin, 3, phonr, event_weight);
           FillDenDenHistograms(selptbin, 3, phonr, event_weight);
          }
            if(passTriggers175){
             FillDenBkgHistograms(selptbin, 4, phonr, event_weight);
             FillDenDenHistograms(selptbin, 4, phonr, event_weight);
            }
            if(passTriggers175 && passMET){
             FillDenBkgHistograms(selptbin, 5, phonr, event_weight);
             FillDenDenHistograms(selptbin, 5, phonr, event_weight);
            }
            if(passTriggers250){
             FillDenBkgHistograms(selptbin, 6, phonr, event_weight);
             FillDenDenHistograms(selptbin, 6, phonr, event_weight);
            }
            if(passTriggers250 && passMET){
             FillDenBkgHistograms(selptbin, 7, phonr, event_weight);
             FillDenDenHistograms(selptbin, 7, phonr, event_weight);
            }
            if(passTriggers && passTriggers175 && passTriggers250 && passMET){
             FillDenBkgHistograms(selptbin, 8, phonr, event_weight);
             FillDenDenHistograms(selptbin, 8, phonr, event_weight);
            }
         }

      // fully inclusive in pT (depending only on selections for nmhDenPCvint)
      FillDenBkgHistograms(incptbin, 0, phonr, event_weight);
      FillDenDenHistograms(incptbin, 0, phonr, event_weight);
      if(passMET){
       FillDenBkgHistograms(incptbin, 1, phonr, event_weight);
       FillDenDenHistograms(incptbin, 1, phonr, event_weight);
      }
      if(passTriggers){
       FillDenBkgHistograms(incptbin, 2, phonr, event_weight);
       FillDenDenHistograms(incptbin, 2, phonr, event_weight);
      }
      if(passTriggers && passMET){
       FillDenBkgHistograms(incptbin, 3, phonr, event_weight);
       FillDenDenHistograms(incptbin, 3, phonr, event_weight);
      }
            if(passTriggers175){
             FillDenBkgHistograms(incptbin, 4, phonr, event_weight);
             FillDenDenHistograms(incptbin, 4, phonr, event_weight);
            }
            if(passTriggers175 && passMET){
             FillDenBkgHistograms(incptbin, 5, phonr, event_weight);
             FillDenDenHistograms(incptbin, 5, phonr, event_weight);
            }
            if(passTriggers250){
             FillDenBkgHistograms(incptbin, 6, phonr, event_weight);
             FillDenDenHistograms(incptbin, 6, phonr, event_weight);
            }
            if(passTriggers250 && passMET){
             FillDenBkgHistograms(incptbin, 7, phonr, event_weight);
             FillDenDenHistograms(incptbin, 7, phonr, event_weight);
            }
            if(passTriggers && passTriggers175 && passTriggers250 && passMET){
             FillDenBkgHistograms(incptbin, 8, phonr, event_weight);
             FillDenDenHistograms(incptbin, 8, phonr, event_weight);
            }

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
  for(unsigned int j=0; j<cuts.size(); ++j){
   WriteHistograms(i,j);
  }
 }
 outfile->Close();
 sw.Stop();
 std::cout<<"Real Time: "<<sw.RealTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"CPU Time: "<<sw.CpuTime()/60.0 <<" minutes"<<std::endl;

 std::cout<<" ntot "<< ntot << std::endl;

 std::cout<<" n_mch_all      "<<  n_mch_all      <<std::endl;    
 std::cout<<" n_mch_met      "<<  n_mch_met      <<std::endl;    
 std::cout<<" n_mch_trg      "<<  n_mch_trg      <<std::endl;    
 std::cout<<" n_mch_met_trg  "<<  n_mch_met_trg  <<std::endl;

 std::cout<<" n_nmh_all      "<<  n_nmh_all      <<std::endl;  
 std::cout<<" n_nmh_met      "<<  n_nmh_met      <<std::endl;  
 std::cout<<" n_nmh_trg      "<<  n_nmh_trg      <<std::endl;  
 std::cout<<" n_nmh_met_trg  "<<  n_nmh_met_trg  <<std::endl;   


 std::cout<<"Done"<<std::endl;

} //end Loop()

//------------------------------------
// photon candidates passing SIGNAL selection in low MET region - (no sieie or (w)chiso cuts)
std::vector<int> postAnalyzerMC_purity::pcPassSel( double phoPtLo, double phoPtHi, double phoEtaMax ){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool noncoll;
  bool passKinematics;

  bool passHoEPSeed;
  bool passPhoNHMedIso;
  bool passCHMedIso; 

  bool passPhotonID;

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

     passPhotonID = passHoEPSeed && passPhoNHMedIso && passCHMedIso;
     if( noncoll && passPhotonID && passKinematics){
      // std::cout<<" Found a photon, pfMET="<<pfMET<<" pT="<<phoEt->at(p)<<" sel: "<<sel<<" sys: "<<sys<<std::endl;
      // std::cout<<"  CHiso:  "<<TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0)<<std::endl;
       tmpCand.push_back(p);
     }
    }
  return tmpCand;
}

//------------------------------------
// photon candidates passing DENOMINATOR selection in low MET region - (no sieie or (w)chiso cuts)
std::vector<int> postAnalyzerMC_purity::pcPassDenSel( double phoPtLo, double phoPtHi, double phoEtaMax ){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool noncoll;
  bool passKinematics;

  bool passHoE;
  bool passPSeed;

  float vloosePFPhoton;
  float vloosePFNeutral;
  bool passVLooseIso;
  bool passLoosePIso;
  bool passLooseIso;

  bool passPhoNHMedIso;
  bool passCHMedIso; 

  bool passPhotonID;

  //Loop over photons
  for(int p=0;p<nPho;p++)  // nPho from ntuple (should = phoE->length() )
    {
     noncoll = (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001; // isMC

     passKinematics = (
                       ( (*phoEt)[p] > phoPtLo  ) &&
                       ( (*phoEt)[p] <= phoPtHi ) &&
                       ( fabs((*phoSCEta)[p]) < phoEtaMax)
                      );


     // H on E, and pixel seed 
     passHoE   = ((*phoHoverE)[p] < 0.05 );
     passPSeed = ((*phohasPixelSeed)[p] ==  0 );

     // deno
     //vloosePFCharged= TMath::Min(5.0*(3.32) , 0.20*(*phoEt)[p]);
     vloosePFPhoton = TMath::Min(5.0*(0.81+ (0.0053 * (*phoEt)[p])) , 0.20*(*phoEt)[p]);
     vloosePFNeutral= TMath::Min(5.0*(1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))) , 0.20*(*phoEt)[p]);
     passVLooseIso = ( 
                      //( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < vloosePFCharged )  &&  
                      ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < vloosePFNeutral )  &&  
                      ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < vloosePFPhoton )
                     );  
     passLoosePIso = 
                     ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) <
                       (0.81 + (0.0053 * (*phoEt)[p])) );
     passLooseIso = (  // deno must fail this cut
                     //( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 3.32 )  &&  
                     ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) <
                       (1.92 + (0.014* (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))  &&  
                     passLoosePIso
                    ); 

     //passPhotonID = passHoE && passPSeed; // && passVLooseIso;
     passPhotonID = passHoE && passPSeed && !passLooseIso && passVLooseIso;
     if( noncoll && passPhotonID && passKinematics){
      // std::cout<<" Found a photon, pfMET="<<pfMET<<" pT="<<phoEt->at(p)<<" sel: "<<sel<<" sys: "<<sys<<std::endl;
      // std::cout<<"  CHiso:  "<<TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0)<<std::endl;
       tmpCand.push_back(p);
     }
    }
  return tmpCand;
}


//-------------------------electron_passLooseID
std::vector<int> postAnalyzerMC_purity::electron_passLooseID(int pho_index, float elePtCut)
{
  //bool veto_passed = true; //pass veto if no good electron found 
  std::vector<int> elelist;

  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false;
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue   
  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
    {
      //Make sure these get reset for every electron  
      pass_SigmaIEtaIEtaFull5x5 = false;
      pass_dEtaIn = false;
      pass_dPhiIn = false;
      pass_HoverE = false;
      pass_iso = false;
      pass_ooEmooP = false;
      pass_d0 = false;
      pass_dz = false;
      pass_missingHits = false;
      pass_convVeto = false;
      //Find EA for corrected relative iso.  
      if(abs(eleSCEta->at(i)) <= 1.0)
        EA = 0.1752;
      else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
        EA = 0.1862;
      else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
        EA = 0.1411;
      else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
        EA = 0.1534;
      else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
        EA = 0.1903;
      else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
        EA = 0.2243;
      else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
        EA = 0.2687;
      EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));

      if(abs(eleSCEta->at(i)) <= 1.479)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
          pass_HoverE = eleHoverE->at(i) < 0.104;
          pass_iso = EAcorrIso < 0.0893;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
          pass_d0 = abs(eleD0->at(i)) < 0.0261;
          pass_dz = abs(eleDz->at(i)) < 0.41;
          pass_missingHits = eleMissHits->at(i) <= 2;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }
      else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
          pass_HoverE = eleHoverE->at(i) < 0.0897;
          pass_iso = EAcorrIso < 0.121;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
          pass_d0 = abs(eleD0->at(i)) < 0.118;
          pass_dz = abs(eleDz->at(i)) < 0.822;
          pass_missingHits = eleMissHits->at(i) <= 1;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }

      //Electron passes Loose Electron ID cuts 
      if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
        {
          //Electron passes pt cut 
          if(elePt->at(i) > elePtCut)
            {
              //Electron does not overlap photon 
              if(dR(eleSCEta->at(i),eleSCPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
                {
                  elelist.push_back(i);
                }
            }
        }
    }
  return elelist;
}


//-------------------------muon_passLooseID
std::vector<int> postAnalyzerMC_purity::muon_passLooseID(int pho_index, float muPtCut)
{
  std::vector<int> mulist;

  bool pass_PFMuon = false;
  bool pass_globalMuon = false;
  bool pass_trackerMuon = false;
  bool pass_iso = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue 
  Float_t zero = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
    {

      pass_PFMuon = muIsPFMuon->at(i);
      pass_globalMuon = muIsGlobalMuon->at(i);
      pass_trackerMuon = muIsTrackerMuon->at(i);
      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.25;
      //if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) )
        {
          //Muon passes pt cut 
          if(muPt->at(i) > muPtCut)
            {
              //Muon does not overlap photon
              if(dR(muEta->at(i),muPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
                {
                 mulist.push_back(i);
                }
            }
        }
    }
  return mulist;
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
//-------------------------passdphiJetMET
bool postAnalyzerMC_purity::passdphiJetMET(std::vector<int> *jets, double mephi)
{

  //pass veto if no jet is found within DeltaPhi(jet,MET) < 0.5                                                                                  
       
  bool passes = false;
  int njetsMax = jets->size();
  //Only look at first four jets                                                                                                                 
       
  if(njetsMax > 4)
    njetsMax = 4;
  int j = 0;
  for(; j < njetsMax; j++)
    {
      //fail veto if a jet is found with DeltaPhi(jet,MET) < 0.5                                                                                 
   
      if(DeltaPhi(jetPhi->at(jets->at(j)), mephi) < 0.5)
        break;
    }
  if(j==njetsMax)
    passes = true;

  return passes;
}


