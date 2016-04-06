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

   // old cuts from james' code
   orecPCvint.clear(); 
   omchPCvint.clear(); 
   onmhPCvint.clear(); 

   // Find all reco photons passing selections
   recPCvint = pcPassSel(175.,2000.,1.4442);
   orecPCvint = pcPassOldSel(175.,2000.,1.4442);
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

   //Loop over all    orec    photons to find (GEN)matched and nonmatched photons
   for(int i = 0; i < orecPCvint.size(); i++)
   {
    int phonr = orecPCvint[i];
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
    if(dRmin < 0.1) { omchPCvint.push_back(phonr); }
    if(dRmin > 0.1) { onmhPCvint.push_back(phonr); }
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

   // some more filters
   bool passTriggers;
   bool passMET;
   passTriggers = ((HLTPho>>12&1) == 1) ;
   passMET = pfMET < 30.;

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
            //std::cout<<"found event"<<std::endl;
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillSigHistograms(ptb, 0, mchPCvint.at(k), event_weight);
            FillDenHistograms(ptb, 0, mchPCvint.at(k), event_weight);
            if(passMET){
             FillSigHistograms(ptb, 1, mchPCvint.at(k), event_weight);
             FillDenHistograms(ptb, 1, mchPCvint.at(k), event_weight);
            }
            if(passTriggers){
             FillSigHistograms(ptb, 2, mchPCvint.at(k), event_weight);
             FillDenHistograms(ptb, 2, mchPCvint.at(k), event_weight);
            }
            if(passTriggers && passMET){
             FillSigHistograms(ptb, 3, mchPCvint.at(k), event_weight);
             FillDenHistograms(ptb, 3, mchPCvint.at(k), event_weight);
            }
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(mchPCvint.at(k)) > 175. ) && //ptbins[0]) &&
         (phoEt->at(mchPCvint.at(k)) < ptbins[ptbins.size()-1])
        ){
          FillSigHistograms(selptbin, 0, mchPCvint.at(k), event_weight);
          FillDenHistograms(selptbin, 0, mchPCvint.at(k), event_weight);
          if(passMET){
           FillSigHistograms(selptbin, 1, mchPCvint.at(k), event_weight);
           FillDenHistograms(selptbin, 1, mchPCvint.at(k), event_weight);
          }
          if(passTriggers){
           FillSigHistograms(selptbin, 2, mchPCvint.at(k), event_weight);
           FillDenHistograms(selptbin, 2, mchPCvint.at(k), event_weight);
          }
          if(passTriggers && passMET){
           FillSigHistograms(selptbin, 3, mchPCvint.at(k), event_weight);
           FillDenHistograms(selptbin, 3, mchPCvint.at(k), event_weight);
          }
         }

      // fully inclusive in pT (depending only on selections for mchPCvint)
      FillSigHistograms(incptbin, 0, mchPCvint.at(k), event_weight);
      FillDenHistograms(incptbin, 0, mchPCvint.at(k), event_weight);
      if(passMET){
       FillSigHistograms(incptbin, 1, mchPCvint.at(k), event_weight);
       FillDenHistograms(incptbin, 1, mchPCvint.at(k), event_weight);
      }
      if(passTriggers){
       FillSigHistograms(incptbin, 2, mchPCvint.at(k), event_weight);
       FillDenHistograms(incptbin, 2, mchPCvint.at(k), event_weight);
      }
      if(passTriggers && passMET){
       FillSigHistograms(incptbin, 3, mchPCvint.at(k), event_weight);
       FillDenHistograms(incptbin, 3, mchPCvint.at(k), event_weight);
      }

     }
    }


    // Fill Numerator (unmatched) Bkgnal Histograms
    if( nmhPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<nmhPCvint.size(); ++k){ // go through each photon passing sel
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(nmhPCvint.at(k)) > ptbins[ptb]) &&
          (phoEt->at(nmhPCvint.at(k)) < ptbins[ptb+1])
         ){
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillBkgHistograms(ptb, 0, nmhPCvint.at(k), event_weight);
            FillDenHistograms(ptb, 0, nmhPCvint.at(k), event_weight);
            if(passMET){
             FillBkgHistograms(ptb, 1, nmhPCvint.at(k), event_weight);
             FillDenHistograms(ptb, 1, nmhPCvint.at(k), event_weight);
            }
            if(passTriggers){
             FillBkgHistograms(ptb, 2, nmhPCvint.at(k), event_weight);
             FillDenHistograms(ptb, 2, nmhPCvint.at(k), event_weight);
            }
            if(passTriggers && passMET){
             FillBkgHistograms(ptb, 3, nmhPCvint.at(k), event_weight);
             FillDenHistograms(ptb, 3, nmhPCvint.at(k), event_weight);
            }
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(nmhPCvint.at(k)) > 175. ) && // ptbins[0]) &&
         (phoEt->at(nmhPCvint.at(k)) < ptbins[ptbins.size()-1])
        ){
          FillBkgHistograms(selptbin, 0, nmhPCvint.at(k), event_weight);
          FillDenHistograms(selptbin, 0, nmhPCvint.at(k), event_weight);
          if(passMET){
           FillBkgHistograms(selptbin, 1, nmhPCvint.at(k), event_weight);
           FillDenHistograms(selptbin, 1, nmhPCvint.at(k), event_weight);
          }
          if(passTriggers){
           FillBkgHistograms(selptbin, 2, nmhPCvint.at(k), event_weight);
           FillDenHistograms(selptbin, 2, nmhPCvint.at(k), event_weight);
          }
          if(passTriggers && passMET){
           FillBkgHistograms(selptbin, 3, nmhPCvint.at(k), event_weight);
           FillDenHistograms(selptbin, 3, nmhPCvint.at(k), event_weight);
          }
         }

      // fully inclusive in pT (depending only on selections for nmhPCvint)
      FillBkgHistograms(incptbin, 0, nmhPCvint.at(k), event_weight);
      FillDenHistograms(incptbin, 0, nmhPCvint.at(k), event_weight);
      if(passMET){
       FillBkgHistograms(incptbin, 1, nmhPCvint.at(k), event_weight);
       FillDenHistograms(incptbin, 1, nmhPCvint.at(k), event_weight);
      }
      if(passTriggers){
       FillBkgHistograms(incptbin, 2, nmhPCvint.at(k), event_weight);
       FillDenHistograms(incptbin, 2, nmhPCvint.at(k), event_weight);
      }
      if(passTriggers && passMET){
       FillBkgHistograms(incptbin, 3, nmhPCvint.at(k), event_weight);
       FillDenHistograms(incptbin, 3, nmhPCvint.at(k), event_weight);
      }

     }
    }

    /////////////////////////////// Denominator selections
    // Fill Numerator (matched) Denom Histograms
    if( mchDenPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<mchDenPCvint.size(); ++k){ // go through each photon passing sel
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(mchDenPCvint.at(k)) > ptbins[ptb]) &&
          (phoEt->at(mchDenPCvint.at(k)) < ptbins[ptb+1])
         ){
            //std::cout<<"found event"<<std::endl;
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillDenSigHistograms(ptb, 0, mchDenPCvint.at(k), event_weight);
            FillDenDenHistograms(ptb, 0, mchDenPCvint.at(k), event_weight);
            if(passMET){
             FillDenSigHistograms(ptb, 1, mchDenPCvint.at(k), event_weight);
             FillDenDenHistograms(ptb, 1, mchDenPCvint.at(k), event_weight);
            }
            if(passTriggers){
             FillDenSigHistograms(ptb, 2, mchDenPCvint.at(k), event_weight);
             FillDenDenHistograms(ptb, 2, mchDenPCvint.at(k), event_weight);
            }
            if(passTriggers && passMET){
             FillDenSigHistograms(ptb, 3, mchDenPCvint.at(k), event_weight);
             FillDenDenHistograms(ptb, 3, mchDenPCvint.at(k), event_weight);
            }
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(mchDenPCvint.at(k)) > 175. ) && //ptbins[0]) &&
         (phoEt->at(mchDenPCvint.at(k)) < ptbins[ptbins.size()-1])
        ){
          FillDenSigHistograms(selptbin, 0, mchDenPCvint.at(k), event_weight);
          FillDenDenHistograms(selptbin, 0, mchDenPCvint.at(k), event_weight);
          if(passMET){
           FillDenSigHistograms(selptbin, 1, mchDenPCvint.at(k), event_weight);
           FillDenDenHistograms(selptbin, 1, mchDenPCvint.at(k), event_weight);
          }
          if(passTriggers){
           FillDenSigHistograms(selptbin, 2, mchDenPCvint.at(k), event_weight);
           FillDenDenHistograms(selptbin, 2, mchDenPCvint.at(k), event_weight);
          }
          if(passTriggers && passMET){
           FillDenSigHistograms(selptbin, 3, mchDenPCvint.at(k), event_weight);
           FillDenDenHistograms(selptbin, 3, mchDenPCvint.at(k), event_weight);
          }
         }

      // fully inclusive in pT (depending only on selections for mchDenPCvint)
      FillDenSigHistograms(incptbin, 0, mchDenPCvint.at(k), event_weight);
      FillDenDenHistograms(incptbin, 0, mchDenPCvint.at(k), event_weight);
      if(passMET){
       FillDenSigHistograms(incptbin, 1, mchDenPCvint.at(k), event_weight);
       FillDenDenHistograms(incptbin, 1, mchDenPCvint.at(k), event_weight);
      }
      if(passTriggers){
       FillDenSigHistograms(incptbin, 2, mchDenPCvint.at(k), event_weight);
       FillDenDenHistograms(incptbin, 2, mchDenPCvint.at(k), event_weight);
      }
      if(passTriggers && passMET){
       FillDenSigHistograms(incptbin, 3, mchDenPCvint.at(k), event_weight);
       FillDenDenHistograms(incptbin, 3, mchDenPCvint.at(k), event_weight);
      }

     }
    }


    // FillDen Numerator (unmatched) Den Histograms
    if( nmhDenPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<nmhDenPCvint.size(); ++k){ // go through each photon passing sel
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(nmhDenPCvint.at(k)) > ptbins[ptb]) &&
          (phoEt->at(nmhDenPCvint.at(k)) < ptbins[ptb+1])
         ){
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillDenBkgHistograms(ptb, 0, nmhDenPCvint.at(k), event_weight);
            FillDenDenHistograms(ptb, 0, nmhDenPCvint.at(k), event_weight);
            if(passMET){
             FillDenBkgHistograms(ptb, 1, nmhDenPCvint.at(k), event_weight);
             FillDenDenHistograms(ptb, 1, nmhDenPCvint.at(k), event_weight);
            }
            if(passTriggers){
             FillDenBkgHistograms(ptb, 2, nmhDenPCvint.at(k), event_weight);
             FillDenDenHistograms(ptb, 2, nmhDenPCvint.at(k), event_weight);
            }
            if(passTriggers && passMET){
             FillDenBkgHistograms(ptb, 3, nmhDenPCvint.at(k), event_weight);
             FillDenDenHistograms(ptb, 3, nmhDenPCvint.at(k), event_weight);
            }
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(nmhDenPCvint.at(k)) > 175. ) && // ptbins[0]) &&
         (phoEt->at(nmhDenPCvint.at(k)) < ptbins[ptbins.size()-1])
        ){
          FillDenBkgHistograms(selptbin, 0, nmhDenPCvint.at(k), event_weight);
          FillDenDenHistograms(selptbin, 0, nmhDenPCvint.at(k), event_weight);
          if(passMET){
           FillDenBkgHistograms(selptbin, 1, nmhDenPCvint.at(k), event_weight);
           FillDenDenHistograms(selptbin, 1, nmhDenPCvint.at(k), event_weight);
          }
          if(passTriggers){
           FillDenBkgHistograms(selptbin, 2, nmhDenPCvint.at(k), event_weight);
           FillDenDenHistograms(selptbin, 2, nmhDenPCvint.at(k), event_weight);
          }
          if(passTriggers && passMET){
           FillDenBkgHistograms(selptbin, 3, nmhDenPCvint.at(k), event_weight);
           FillDenDenHistograms(selptbin, 3, nmhDenPCvint.at(k), event_weight);
          }
         }

      // fully inclusive in pT (depending only on selections for nmhDenPCvint)
      FillDenBkgHistograms(incptbin, 0, nmhDenPCvint.at(k), event_weight);
      FillDenDenHistograms(incptbin, 0, nmhDenPCvint.at(k), event_weight);
      if(passMET){
       FillDenBkgHistograms(incptbin, 1, nmhDenPCvint.at(k), event_weight);
       FillDenDenHistograms(incptbin, 1, nmhDenPCvint.at(k), event_weight);
      }
      if(passTriggers){
       FillDenBkgHistograms(incptbin, 2, nmhDenPCvint.at(k), event_weight);
       FillDenDenHistograms(incptbin, 2, nmhDenPCvint.at(k), event_weight);
      }
      if(passTriggers && passMET){
       FillDenBkgHistograms(incptbin, 3, nmhDenPCvint.at(k), event_weight);
       FillDenDenHistograms(incptbin, 3, nmhDenPCvint.at(k), event_weight);
      }

     }
    }

    /////////////////////////////// old selections
    // Fill Numerator (matched) Signal Histograms
    if( omchPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<omchPCvint.size(); ++k){ // go through each photon passing sel
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(omchPCvint.at(k)) > ptbins[ptb]) &&
          (phoEt->at(omchPCvint.at(k)) < ptbins[ptb+1])
         ){
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillSigHistograms(ptb, 4, omchPCvint.at(k), event_weight);
            FillDenHistograms(ptb, 4, omchPCvint.at(k), event_weight);
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(omchPCvint.at(k)) > 175. ) && // ptbins[0]) &&
         (phoEt->at(omchPCvint.at(k)) < ptbins[ptbins.size()-1])
        ){
          FillSigHistograms(selptbin, 4, omchPCvint.at(k), event_weight);
          FillDenHistograms(selptbin, 4, omchPCvint.at(k), event_weight);
         }

      // fully inclusive in pT (depending only on selections for omchPCvint)
      FillSigHistograms(incptbin, 4, omchPCvint.at(k), event_weight);
      FillDenHistograms(incptbin, 4, omchPCvint.at(k), event_weight);

     }
    }
    // Fill Numerator (unmatched) Bkgnal Histograms
    if( onmhPCvint.size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<onmhPCvint.size(); ++k){ // go through each photon passing sel
      // pt bins
      for(unsigned int ptb=0; ptb<lstptbin; ++ptb){
       if(
          (phoEt->at(onmhPCvint.at(k)) > ptbins[ptb]) &&
          (phoEt->at(onmhPCvint.at(k)) < ptbins[ptb+1])
         ){
            // idnc, idnc_mL30, idnc_trig, idnc_mL30_trig, oldj
            //nc++;
            FillBkgHistograms(ptb, 4, onmhPCvint.at(k), event_weight);
            FillDenHistograms(ptb, 4, onmhPCvint.at(k), event_weight);
           }
       }
      // inclusive from pT bins
      if(  
         (phoEt->at(onmhPCvint.at(k)) > 175. ) && // ptbins[0]) &&
         (phoEt->at(onmhPCvint.at(k)) < ptbins[ptbins.size()-1])
        ){
          FillBkgHistograms(selptbin, 4, onmhPCvint.at(k), event_weight);
          FillDenHistograms(selptbin, 4, onmhPCvint.at(k), event_weight);
         }

      // fully inclusive in pT (depending only on selections for onmhPCvint)
      FillBkgHistograms(incptbin, 4, onmhPCvint.at(k), event_weight);
      FillDenHistograms(incptbin, 4, onmhPCvint.at(k), event_weight);

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


// photon candidates passing signal selection in low MET region - (no sieie or (w)chiso cuts)
std::vector<int> postAnalyzerMC_purity::pcPassOldSel( double phoPtLo, double phoPtHi, double phoEtaMax ){
        bool HEpass = false;
        //bool sigmaietaietapass = false;
        //bool corrChargedIsopass = false;
        bool corrNeutralIsopass = false;
        bool corrPhotonIsopass = false;

        std::vector<int> tmpCand;
        tmpCand.clear();
        //Loop over photons
        for(int p=0;p<nPho;p++)  // nPho from ntuple (should = phoE->length() )
        {

        float photon_pt = phoEt->at(p);
        float SC_eta = phoSCEta->at(p);
        float HE = phoHoverE->at(p);
         //float sigmaietaieta = phoSigmaIEtaIEta->at(p);

         bool passKinematics;
         passKinematics = (
                           ( photon_pt > phoPtLo  ) &&
                           ( photon_pt <= phoPtHi ) &&
                           ( fabs(SC_eta) < phoEtaMax)
                          );
      
         //float EAcharged = 0.0;
         float EAneutral = 0.0;
         float EAphoton = 0.0;
         if(abs(SC_eta) <= 1.0)
         {   
                 //EAcharged = 0.0157;
                 EAneutral = 0.0143;
                 EAphoton = 0.0725;
         }   
         else if(1.0 < abs(SC_eta) && abs(SC_eta) <= 1.479)
         {   
                 //EAcharged = 0.0143;
                 EAneutral = 0.0210;
                 EAphoton = 0.0604;
         }   
         else if(1.479 < abs(SC_eta) && abs(SC_eta) <= 2.0)
         {   
                 //EAcharged = 0.0115;
                 EAneutral = 0.0147;
                 EAphoton = 0.0320;
         }   
         else if(2.0 < abs(SC_eta) && abs(SC_eta) <= 2.2)
         {   
                 //EAcharged = 0.0094;
                 EAneutral = 0.0082;
                 EAphoton = 0.0512;
         }   
         else if(2.2 < abs(SC_eta) && abs(SC_eta) <= 2.3)
         {   
                 //EAcharged = 0.0095;
                 EAneutral = 0.0124;
                 EAphoton = 0.0766;
         }   
         else if(2.3 < abs(SC_eta) && abs(SC_eta) <= 2.4)
         {   
                 //EAcharged = 0.0068;
                 EAneutral = 0.0186;
                 EAphoton = 0.0949;
         }   
         else if(2.4 < abs(SC_eta))
         {   
                 //EAcharged = 0.0053;
                 EAneutral = 0.0320;
                 EAphoton = 0.1160;
         }   

         //Being very explicit about types here, to avoid a conversion ambiguity
         //Float_t corrChargedIso_nonzero = phoPFChIso->at(p) - rho*EAcharged;
         Float_t corrNeutralIso_nonzero = phoPFNeuIso->at(p)  - rho*EAneutral;
         Float_t corrPhotonIso_nonzero = phoPFPhoIso->at(p) - rho*EAphoton;
         Float_t zero = 0.0;
         //Float_t corrChargedIso = TMath::Max(corrChargedIso_nonzero,zero);
         Float_t corrNeutralIso = TMath::Max(corrNeutralIso_nonzero,zero);
         Float_t corrPhotonIso = TMath::Max(corrPhotonIso_nonzero,zero);
      
         //Barrel
         if(SC_eta < 1.4442)
         {   
                 HEpass = HE < 0.05;
                 //sigmaietaietapass = sigmaietaieta < 0.0100;
                 //corrChargedIsopass = corrChargedIso < 1.31;
                 corrNeutralIsopass = corrNeutralIso < (0.60 + exp(0.0044*photon_pt + 0.5809));
                 corrPhotonIsopass = corrPhotonIso < (1.33 + 0.0043*photon_pt);
         }   
         //Endcap
         else if(SC_eta > 1.566)
         {   
                 HEpass = HE < 0.05;
                 //sigmaietaietapass = sigmaietaieta < 0.0267;
                 //corrChargedIsopass = corrChargedIso < 1.25;
                 corrNeutralIsopass = corrNeutralIso < (1.65 + exp(0.0040*photon_pt + 0.9402));
                 corrPhotonIsopass = corrPhotonIso < (1.02 + 0.0041*photon_pt);
         }   
      
         if(HEpass && corrNeutralIsopass && corrPhotonIsopass && passKinematics) {
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

