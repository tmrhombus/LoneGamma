#define postAnalyzer_QCD_cxx
#include "postAnalyzer_QCD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void postAnalyzer_QCD::Loop(TString outfilename, Bool_t isMC, Double_t lumi, Double_t nrEvents, Double_t crossSec)
{

 TStopwatch sw; 
 sw.Start();

 Int_t nc = 0;
 Int_t dc = 0;

 if (fChain == 0) return;
 std::cout<<"Numerator/Denominator : run : lumis : event : photonindex "<<std::endl;

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
  
           if (
              ( run==258446 && event==121132593 && lumis==97  )  
           || ( run==258702 && event==280995575 && lumis==180 ) 
           || ( run==258712 && event==494069084 && lumis==312 ) 
           || ( run==258712 && event==22840508  && lumis==15  )
           || ( run==258712 && event==64168685  && lumis==40  )
           || ( run==258712 && event==88413644  && lumis==55  )
           || ( run==259809 && event==177969654 && lumis==138 ) 
           || ( run==259862 && event==90465340  && lumis==53  )
           || ( run==259862 && event==89817711  && lumis==53  )
             ) {
            std::cout<<"A : "<<run<<" : "<<event<<" : "<<lumis<<std::endl;
       //     std::cout<<"Ashim Event  oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"<<std::endl;
       //     std::cout<<" run           "<<   run            <<std::endl;
       //     std::cout<<" lumis          "<<   lumis           <<std::endl;
       //     std::cout<<" event         "<<   event          <<std::endl;
       //     std::cout<<" (HLTPho>>7&1)  "<< (HLTPho>>7&1) <<std::endl; 
       //     std::cout<<" (HLTPho>>8&1)  "<< (HLTPho>>8&1) <<std::endl; 
       //     std::cout<<" (HLTPho>>9&1)  "<< (HLTPho>>9&1) <<std::endl; 
       //     std::cout<<" (HLTPho>>10&1) "<< (HLTPho>>10&1)<<std::endl; 
       //     std::cout<<" (HLTPho>>11&1) "<< (HLTPho>>11&1)<<std::endl; 
       //     std::cout<<" (HLTPho>>12&1) "<< (HLTPho>>12&1)<<std::endl; 
       //     std::cout<<" (HLTPho>>22&1) "<< (HLTPho>>22&1)<<std::endl; 
       //  for(int p=0;p<nPho;p++)  // nPho from ntuple (should = phoE->length() )
       //    {
       //     std::cout<<"..................................................................."<<std::endl;
       //     std::cout<<" photon index  "<<   p     <<std::endl;

       //     std::cout<<" phoEt          "<<   (*phoEt)[p]          <<std::endl; 
       //     std::cout<<" phoSCEta       "<<   (*phoSCEta)[p]       <<std::endl;
       //     std::cout<<" pfMET          "<<   pfMET          <<std::endl;
       //     std::cout<<" phoPFNeuIso    "<<   (*phoPFNeuIso)[p]    <<std::endl;
       //     std::cout<<" phoPFPhoIso    "<<   (*phoPFPhoIso)[p]    <<std::endl;
       //     std::cout<<" phoPFChIso     "<<   (*phoPFChIso) [p]     <<std::endl;
       //     std::cout<<" rho            "<<   rho            <<std::endl;
       //     std::cout<<" EAneutral      "<<   EAneutral((*phoSCEta)[p]) <<std::endl;
       //     std::cout<<" EAphoton       "<<   EAphoton((*phoSCEta)[p])       <<std::endl<<std::endl;   

       //     std::cout<<"  To Pass Denominator Cut  sssssssssssssssssssssssssssssss"<<std::endl;
       //     std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
       //     std::cout<<" (V)Loose ID Variables"<<std::endl;
       //     std::cout<<"Neutral Iso     NeuIso - rho*EA = "<<
       //      (*phoPFNeuIso)[p] - ( rho * EAneutral((*phoSCEta)[p]) )<<std::endl;
       //     std::cout<<"  1.92 + 0.014pT + 0.000019pT^2 = "<<  
       //      (1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0)))<<std::endl;
       //     std::cout<<"5(1.92 + 0.014pT + 0.000019pT^2)= "<<  
       //      (5 * (1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))<<std::endl;
       //     std::cout<<"                        0.2(pT) = "<<(0.2 * ((*phoEt)[p]))<<std::endl;
       //     std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
       //     std::cout<<"Photon Iso      PhoIso - rho*EA = "<<
       //      (*phoPFPhoIso)[p] - ( rho * EAphoton((*phoSCEta)[p]) )<<std::endl;
       //     std::cout<<"                 0.81 + 0.053pT = "<<  
       //      (0.81 + (0.053 * (*phoEt)[p]) )<<std::endl;
       //     std::cout<<"              5(0.81 + 0.053pT) = "<<  
       //      (5 * (0.81 + (0.053 * (*phoEt)[p]) ))<<std::endl;
       //     std::cout<<"                        0.2(pT) = "<<(0.2 * ((*phoEt)[p]))<<std::endl;
       //     std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
       //     std::cout<<"CHad Iso                        = "<< (*phoPFChIso)[p] <<std::endl;
       //     std::cout<<"                        0.2(pT) = "<<(0.2 * ((*phoEt)[p]))<<std::endl;
       //     std::cout<<"---------------------------------------"<<std::endl;
       //     std::cout<<" Cuts are:"<<std::endl;
       //     std::cout<<"  pT:             175 < "<<(*phoEt)[p]<<" < 190"<<std::endl;
       //     std::cout<<"  eta:            | "<<(*phoSCEta)[p]<<" | < 1.442"<<std::endl;
       //     std::cout<<"  met:            "<<pfMET<<" < 30"<<std::endl;
       //     std::cout<<"  h/e:            "<<(*phoHoverE)[p]<<" < 0.05"<<std::endl;
       //     std::cout<<"  pixel seed:     "<<(*phohasPixelSeed)[p]<<" == 0"<<std::endl;
       //     std::cout<<"  |seed time|     3 > "<<    fabs((*phoseedTimeFull5x5)[p])<<std::endl;
       //     std::cout<<"  MIP             6.3 > "<<    fabs((*phomipTotEnergy)[p])<<std::endl;
       //     std::cout<<"  sieie           0.001 < "<<    fabs((*phoSigmaIEtaIEtaFull5x5)[p])<<std::endl;
       //     std::cout<<"  sipip           0.001 < "<<    fabs((*phoSigmaIPhiIPhiFull5x5)[p])<<std::endl;
       //     std::cout<<"  pass vloose:    must pass all"<<std::endl;

       //      std::cout<<"    PFchIso:        max( 0, "<<
       //      (*phoPFChIso)[p]<< " ) < min( 16.6, "<<
       //      (0.2 * ((*phoEt)[p]))<<" )"<<std::endl;

       //      std::cout<<"    PFphoIso:       max( 0, "<<
       //      (*phoPFPhoIso)[p] - ( rho * EAphoton((*phoSCEta)[p]) )
       //      << " ) < min( "<<
       //      (5 * (0.81 + (0.053 * (*phoEt)[p]) ))
       //      <<", "<<
       //      (0.2 * ((*phoEt)[p]))<<" )"<<std::endl;

       //      std::cout<<"    PFneuIso:       max( 0, "<<
       //      (*phoPFNeuIso)[p] - ( rho * EAneutral((*phoSCEta)[p]) )
       //      << " ) < min( "<<
       //      (5 * (1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))))
       //      <<", "<<
       //      (0.2 * ((*phoEt)[p]))<<" )"<<std::endl;

       //     std::cout<<"  fail loose:    must fail one of"<<std::endl;

       //      std::cout<<"    PFchIso:        max( 0, "<<
       //      (*phoPFChIso)[p]<< " ) < min( 3.32, "<<
       //      (0.2 * ((*phoEt)[p]))<<" )"<<std::endl;

       //      std::cout<<"    PFphoIso:       max( 0, "<<
       //      (*phoPFPhoIso)[p] - ( rho * EAphoton((*phoSCEta)[p]) )
       //      << " ) < min( "<<
       //      (0.81 + (0.053 * (*phoEt)[p]) )
       //      <<", "<<
       //      (0.2 * ((*phoEt)[p]))<<" )"<<std::endl;

       //      std::cout<<"    PFneuIso:       max( 0, "<<
       //      (*phoPFNeuIso)[p] - ( rho * EAneutral((*phoSCEta)[p]) )
       //      << " ) < min( "<<
       //      (1.92 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0)))
       //      <<", "<<
       //      (0.2 * ((*phoEt)[p]))<<" )"<<std::endl<<std::endl;


       //     std::cout<<"  To Pass Numerator Cut  sssssssssssssssssssssssssssssss"<<std::endl;
       //     std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
       //     std::cout<<"Neutral Iso     NeuIso - rho*EA = "<<
       //      (*phoPFNeuIso)[p] - ( rho * EAneutral((*phoSCEta)[p]) )<<std::endl;
       //     std::cout<<"  1.06 + 0.014pT + 0.000019pT^2 = "<<  
       //      (1.06 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0)))<<std::endl;
       //     std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
       //     std::cout<<"Photon Iso      PhoIso - rho*EA = "<<
       //      (*phoPFPhoIso)[p] - ( rho * EAphoton((*phoSCEta)[p]) )<<std::endl;
       //     std::cout<<"                 0.28 + 0.053pT = "<<  
       //      (0.28 + (0.053 * (*phoEt)[p]) )<<std::endl;
       //     std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
       //     std::cout<<"CHad Iso                        = "<< (*phoPFChIso)[p] <<" < 1.37"<<std::endl<<std::endl;

       //     std::cout<<"---------------------------------------"<<std::endl;
       //     std::cout<<" Cuts are:"<<std::endl;
       //     std::cout<<"  pT:             175 < "<<(*phoEt)[p]<<" < 190"<<std::endl;
       //     std::cout<<"  eta:            | "<<(*phoSCEta)[p]<<" | < 1.442"<<std::endl;
       //     std::cout<<"  met:            "<<pfMET<<" < 30"<<std::endl;
       //     std::cout<<"  h/e:            "<<(*phoHoverE)[p]<<" < 0.05"<<std::endl;
       //     std::cout<<"  pixel seed:     "<<(*phohasPixelSeed)[p]<<" == 0"<<std::endl;
       //     std::cout<<"  |seed time|     3 > "<<    fabs((*phoseedTimeFull5x5)[p])<<std::endl;
       //     std::cout<<"  MIP             6.3 > "<<    fabs((*phomipTotEnergy)[p])<<std::endl;
       //     std::cout<<"  sieie           0.001 < "<<    fabs((*phoSigmaIEtaIEtaFull5x5)[p])<<std::endl;
       //     std::cout<<"  sipip           0.001 < "<<    fabs((*phoSigmaIPhiIPhiFull5x5)[p])<<std::endl;
       //     std::cout<<"  pass med:         must pass all"<<std::endl;

       //      std::cout<<"    PFchIso:        max( 0, "<<
       //      (*phoPFChIso)[p]<< " ) < 1.37"<<std::endl;

       //      std::cout<<"    PFphoIso:       max( 0, "<<
       //      (*phoPFPhoIso)[p] - ( rho * EAphoton((*phoSCEta)[p]) )
       //      << " ) < "<<
       //      (0.28 + (0.053 * (*phoEt)[p]) )
       //      <<std::endl;

       //      std::cout<<"    PFneuIso:       max( 0, "<<
       //      (*phoPFNeuIso)[p] - ( rho * EAneutral((*phoSCEta)[p]) )
       //      << " ) < "<<
       //      (1.06 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0)))
       //      <<std::endl<<std::endl;
       //     }
            }

  // if event passes MonoPhoton triggers
  if( 
   ((HLTPho>>7&1) == 1) ||
   ((HLTPho>>8&1) == 1) ||
   ((HLTPho>>9&1) == 1) ||
   ((HLTPho>>10&1) == 1) ||
   ((HLTPho>>11&1) == 1) ||
   ((HLTPho>>12&1) == 1) ||
   ((HLTPho>>22&1) == 1) )
  {   

   // sysbinnames = "" "_sbUP" "_sbDown" "_metUP" "_metDown" "_binUP" "_binDown" "_noPiso"
   // ptbins = 175, 190, 250, 400, 700, 1000
   // ptbinnames = "175to190" "190to250" "250to400" "400to700" "700to1000" "175to1000" "allpt"

   // vector of ints, each int corresponds to location in vector of photons of photon passing cut
   // one such vector for each systematic bin

   int inclptbin = ptbins.size()-1;
   int lastptbin = ptbinnames.size()-1;
   int lastsysbin = sysbinnames.size();

   for(unsigned int sysb=0; sysb<lastsysbin; ++sysb){
    sigPCvint[sysb] = pcPassSel(0,sysb,isMC); // passes signal selection (no sieie cut)
    bkgPCvint[sysb] = pcPassSel(1,sysb,isMC); // passes background selection (no sieie cut)
    denPCvint[sysb] = pcPassSel(2,sysb,isMC); // passes denominator selection (no sieie cut)
   }

// fill histograms
   for(unsigned int sysb=0; sysb<lastsysbin; ++sysb){
    // Fill Numerator Signal Histograms
    if( sigPCvint[sysb].size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<sigPCvint[sysb].size(); ++k){ // go through each photon passing sel
      for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
       if(
          (phoEt->at(sigPCvint[sysb].at(k)) > ptbins[ptb]) &&
          (phoEt->at(sigPCvint[sysb].at(k)) < ptbins[ptb+1])
         ){
           if( sysb==0 && ptb==0 ) {
            nc++;
            std::cout<<"N : "<<run<<" : "<<event<<" : "<<lumis<<" : "<<sigPCvint[sysb].at(k)<<std::endl;
            //std::cout<<"Passed Numerator Cut-------------------------------------------------------------"<<std::endl;
            //std::cout<<" run           "<<   run            <<std::endl;
            //std::cout<<" lumis          "<<   lumis           <<std::endl;
            //std::cout<<" event         "<<   event          <<std::endl;
            //std::cout<<" photon index  "<<   sigPCvint[sysb].at(k)     <<std::endl;
      
            //std::cout<<" phoEt         "<<   (*phoEt)[sigPCvint[sysb].at(k)]          <<std::endl; 
            //std::cout<<" phoSCEta      "<<   (*phoSCEta)[sigPCvint[sysb].at(k)]       <<std::endl;
            //std::cout<<" pfMET         "<<   pfMET          <<std::endl;
            //std::cout<<" phoPFNeuIso   "<<   (*phoPFNeuIso)[sigPCvint[sysb].at(k)]    <<std::endl;
            //std::cout<<" phoPFPhoIso   "<<   (*phoPFPhoIso)[sigPCvint[sysb].at(k)]    <<std::endl;
            //std::cout<<" phoPFChIso    "<<   (*phoPFChIso) [sigPCvint[sysb].at(k)]     <<std::endl;
            //std::cout<<" rho           "<<   rho            <<std::endl;
            //std::cout<<" EAneutral     "<<   EAneutral((*phoSCEta)[sigPCvint[sysb].at(k)]) <<std::endl;
            //std::cout<<" EAphoton      "<<   EAphoton((*phoSCEta)[sigPCvint[sysb].at(k)])       <<std::endl;   
            //std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
            //std::cout<<"Neutral Iso     NeuIso - rho*EA = "<<
            // (*phoPFNeuIso)[sigPCvint[sysb].at(k)] - ( rho * EAneutral((*phoSCEta)[sigPCvint[sysb].at(k)]) )<<std::endl;
            //std::cout<<"  1.06 + 0.014pT + 0.000019pT^2 = "<<  
            // (1.06 + (0.014 * (*phoEt)[sigPCvint[sysb].at(k)]) + (0.000019 * pow((*phoEt)[sigPCvint[sysb].at(k)], 2.0)))<<std::endl;
            //std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
            //std::cout<<"Photon Iso      PhoIso - rho*EA = "<<
            // (*phoPFPhoIso)[sigPCvint[sysb].at(k)] - ( rho * EAphoton((*phoSCEta)[sigPCvint[sysb].at(k)]) )<<std::endl;
            //std::cout<<"                 0.28 + 0.053pT = "<<  
            // (0.28 + (0.053 * (*phoEt)[sigPCvint[sysb].at(k)]) )<<std::endl;
            //std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
            //std::cout<<"CHad Iso                        = "<< (*phoPFChIso)[sigPCvint[sysb].at(k)] <<" < 1.37"<<std::endl<<std::endl;

            //std::cout<<"---------------------------------------"<<std::endl;
            //std::cout<<" Cuts are:"<<std::endl;
            //std::cout<<"  pT:             175 < "<<(*phoEt)[sigPCvint[sysb].at(k)]<<" < 190"<<std::endl;
            //std::cout<<"  eta:            | "<<(*phoSCEta)[sigPCvint[sysb].at(k)]<<" | < 1.442"<<std::endl;
            //std::cout<<"  met:            "<<pfMET<<" < 30"<<std::endl;
            //std::cout<<"  h/e:            "<<(*phoHoverE)[sigPCvint[sysb].at(k)]<<" < 0.05"<<std::endl;
            //std::cout<<"  pixel seed:     "<<(*phohasPixelSeed)[sigPCvint[sysb].at(k)]<<" == 0"<<std::endl;
            //std::cout<<"  |seed time|     3 > "<<    fabs((*phoseedTimeFull5x5)[sigPCvint[sysb].at(k)])<<std::endl;
            //std::cout<<"  MIP             6.3 > "<<    fabs((*phomipTotEnergy)[sigPCvint[sysb].at(k)])<<std::endl;
            //std::cout<<"  sieie           0.001 < "<<    fabs((*phoSigmaIEtaIEtaFull5x5)[sigPCvint[sysb].at(k)])<<std::endl;
            //std::cout<<"  sipip           0.001 < "<<    fabs((*phoSigmaIPhiIPhiFull5x5)[sigPCvint[sysb].at(k)])<<std::endl;
            //std::cout<<"  pass med:         must pass all"<<std::endl;

            // std::cout<<"    PFchIso:        max( 0, "<<
            // (*phoPFChIso)[sigPCvint[sysb].at(k)]<< " ) < 1.37"<<std::endl;

            // std::cout<<"    PFphoIso:       max( 0, "<<
            // (*phoPFPhoIso)[sigPCvint[sysb].at(k)] - ( rho * EAphoton((*phoSCEta)[sigPCvint[sysb].at(k)]) )
            // << " ) < "<<
            // (0.28 + (0.053 * (*phoEt)[sigPCvint[sysb].at(k)]) )
            // <<std::endl;

            // std::cout<<"    PFneuIso:       max( 0, "<<
            // (*phoPFNeuIso)[sigPCvint[sysb].at(k)] - ( rho * EAneutral((*phoSCEta)[sigPCvint[sysb].at(k)]) )
            // << " ) < "<<
            // (1.06 + (0.014 * (*phoEt)[sigPCvint[sysb].at(k)]) + (0.000019 * pow((*phoEt)[sigPCvint[sysb].at(k)], 2.0)))
            // <<std::endl<<std::endl;

           }
        FillSigHistograms(ptb, sysb, sigPCvint[sysb].at(k), event_weight);
       }
      }
      if(  // do an inclusive pT plot from bins
         (phoEt->at(sigPCvint[sysb].at(k)) > ptbins[0]) &&
         (phoEt->at(sigPCvint[sysb].at(k)) < ptbins[inclptbin])
        ){
       FillSigHistograms(lastptbin-1, sysb, sigPCvint[sysb].at(k), event_weight);
      }
      // and one fully inclusive in pT
      FillSigHistograms(lastptbin, sysb, sigPCvint[sysb].at(k), event_weight);
     }
    }

    // Fill Numerator Background Histograms
    if( bkgPCvint[sysb].size()>0 ){ // if any photon indexes passed sig selection
     if( sysb==0 ) {
     }
     for(unsigned int k=0; k<bkgPCvint[sysb].size(); ++k){ // go through each photon passing sel
      for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
       if(
          (phoEt->at(bkgPCvint[sysb].at(k)) > ptbins[ptb]) &&
          (phoEt->at(bkgPCvint[sysb].at(k)) < ptbins[ptb+1])
         ){
        FillBkgHistograms(ptb, sysb, bkgPCvint[sysb].at(k), event_weight);
       }
      }
      if(  // do an inclusive pT plot from bins
         (phoEt->at(bkgPCvint[sysb].at(k)) > ptbins[0]) &&
         (phoEt->at(bkgPCvint[sysb].at(k)) < ptbins[inclptbin])
        ){
       FillBkgHistograms(lastptbin-1, sysb, bkgPCvint[sysb].at(k), event_weight);
      }
      // and one fully inclusive in pT
      FillBkgHistograms(lastptbin, sysb, bkgPCvint[sysb].at(k), event_weight);
     }
    }

    // Fill Denominator Histograms
    if( denPCvint[sysb].size()>0 ){ // if any photon indexes passed sig selection
     for(unsigned int k=0; k<denPCvint[sysb].size(); ++k){ // go through each photon passing sel
      for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
       if(
          (phoEt->at(denPCvint[sysb].at(k)) > ptbins[ptb]) &&
          (phoEt->at(denPCvint[sysb].at(k)) < ptbins[ptb+1])
         ){
           if( sysb==0 && ptb==0 ) {
            dc++;
            std::cout<<"D : "<<run<<" : "<<event<<" : "<<lumis<<" : "<<denPCvint[sysb].at(k)<<std::endl;
            //std::cout<<"Passed Denominator Cut-----------------------------------------------------------"<<std::endl;
            //std::cout<<" run           "<<   run            <<std::endl;
            //std::cout<<" lumis          "<<   lumis           <<std::endl;
            //std::cout<<" event         "<<   event          <<std::endl;
            //std::cout<<" photon index  "<<   denPCvint[sysb].at(k)     <<std::endl;

            //std::cout<<" phoEt         "<<   (*phoEt)[denPCvint[sysb].at(k)]          <<std::endl; 
            //std::cout<<" phoSCEta      "<<   (*phoSCEta)[denPCvint[sysb].at(k)]       <<std::endl;
            //std::cout<<" pfMET         "<<   pfMET          <<std::endl;
            //std::cout<<" phoPFNeuIso   "<<   (*phoPFNeuIso)[denPCvint[sysb].at(k)]    <<std::endl;
            //std::cout<<" phoPFPhoIso   "<<   (*phoPFPhoIso)[denPCvint[sysb].at(k)]    <<std::endl;
            //std::cout<<" phoPFChIso    "<<   (*phoPFChIso) [denPCvint[sysb].at(k)]     <<std::endl;
            //std::cout<<" rho           "<<   rho            <<std::endl;
            //std::cout<<" EAneutral     "<<   EAneutral((*phoSCEta)[denPCvint[sysb].at(k)]) <<std::endl;
            //std::cout<<" EAphoton      "<<   EAphoton((*phoSCEta)[denPCvint[sysb].at(k)])       <<std::endl<<std::endl;   

            //std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
            //std::cout<<" (V)Loose ID Variables"<<std::endl;
            //std::cout<<"Neutral Iso     NeuIso - rho*EA = "<<
            // (*phoPFNeuIso)[denPCvint[sysb].at(k)] - ( rho * EAneutral((*phoSCEta)[denPCvint[sysb].at(k)]) )<<std::endl;
            //std::cout<<"  1.92 + 0.014pT + 0.000019pT^2 = "<<  
            // (1.92 + (0.014 * (*phoEt)[denPCvint[sysb].at(k)]) + (0.000019 * pow((*phoEt)[denPCvint[sysb].at(k)], 2.0)))<<std::endl;
            //std::cout<<"5(1.92 + 0.014pT + 0.000019pT^2)= "<<  
            // (5 * (1.92 + (0.014 * (*phoEt)[denPCvint[sysb].at(k)]) + (0.000019 * pow((*phoEt)[denPCvint[sysb].at(k)], 2.0))))<<std::endl;
            //std::cout<<"                        0.2(pT) = "<<(0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<std::endl;
            //std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
            //std::cout<<"Photon Iso      PhoIso - rho*EA = "<<
            // (*phoPFPhoIso)[denPCvint[sysb].at(k)] - ( rho * EAphoton((*phoSCEta)[denPCvint[sysb].at(k)]) )<<std::endl;
            //std::cout<<"                 0.81 + 0.053pT = "<<  
            // (0.81 + (0.053 * (*phoEt)[denPCvint[sysb].at(k)]) )<<std::endl;
            //std::cout<<"              5(0.81 + 0.053pT) = "<<  
            // (5 * (0.81 + (0.053 * (*phoEt)[denPCvint[sysb].at(k)]) ))<<std::endl;
            //std::cout<<"                        0.2(pT) = "<<(0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<std::endl;
            //std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
            //std::cout<<"CHad Iso                        = "<< (*phoPFChIso)[denPCvint[sysb].at(k)] <<std::endl;
            //std::cout<<"                        0.2(pT) = "<<(0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<std::endl;
            //std::cout<<"---------------------------------------"<<std::endl;
            //std::cout<<" Cuts are:"<<std::endl;
            //std::cout<<"  pT:             175 < "<<(*phoEt)[denPCvint[sysb].at(k)]<<" < 190"<<std::endl;
            //std::cout<<"  eta:            | "<<(*phoSCEta)[denPCvint[sysb].at(k)]<<" | < 1.442"<<std::endl;
            //std::cout<<"  met:            "<<pfMET<<" < 30"<<std::endl;
            //std::cout<<"  h/e:            "<<(*phoHoverE)[denPCvint[sysb].at(k)]<<" < 0.05"<<std::endl;
            //std::cout<<"  pixel seed:     "<<(*phohasPixelSeed)[denPCvint[sysb].at(k)]<<" == 0"<<std::endl;
            //std::cout<<"  |seed time|     3 > "<<    fabs((*phoseedTimeFull5x5)[denPCvint[sysb].at(k)])<<std::endl;
            //std::cout<<"  MIP             6.3 > "<<    fabs((*phomipTotEnergy)[denPCvint[sysb].at(k)])<<std::endl;
            //std::cout<<"  sieie           0.001 < "<<    fabs((*phoSigmaIEtaIEtaFull5x5)[denPCvint[sysb].at(k)])<<std::endl;
            //std::cout<<"  sipip           0.001 < "<<    fabs((*phoSigmaIPhiIPhiFull5x5)[denPCvint[sysb].at(k)])<<std::endl;
            //std::cout<<"  pass vloose:    must pass all"<<std::endl;

            // std::cout<<"    PFchIso:        max( 0, "<<
            // (*phoPFChIso)[denPCvint[sysb].at(k)]<< " ) < min( 16.6, "<<
            // (0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<" )"<<std::endl;

            // std::cout<<"    PFphoIso:       max( 0, "<<
            // (*phoPFPhoIso)[denPCvint[sysb].at(k)] - ( rho * EAphoton((*phoSCEta)[denPCvint[sysb].at(k)]) )
            // << " ) < min( "<<
            // (5 * (0.81 + (0.053 * (*phoEt)[denPCvint[sysb].at(k)]) ))
            // <<", "<<
            // (0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<" )"<<std::endl;

            // std::cout<<"    PFneuIso:       max( 0, "<<
            // (*phoPFNeuIso)[denPCvint[sysb].at(k)] - ( rho * EAneutral((*phoSCEta)[denPCvint[sysb].at(k)]) )
            // << " ) < min( "<<
            // (5 * (1.92 + (0.014 * (*phoEt)[denPCvint[sysb].at(k)]) + (0.000019 * pow((*phoEt)[denPCvint[sysb].at(k)], 2.0))))
            // <<", "<<
            // (0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<" )"<<std::endl;

            //std::cout<<"  fail loose:    must fail one of"<<std::endl;

            // std::cout<<"    PFchIso:        max( 0, "<<
            // (*phoPFChIso)[denPCvint[sysb].at(k)]<< " ) < min( 3.32, "<<
            // (0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<" )"<<std::endl;

            // std::cout<<"    PFphoIso:       max( 0, "<<
            // (*phoPFPhoIso)[denPCvint[sysb].at(k)] - ( rho * EAphoton((*phoSCEta)[denPCvint[sysb].at(k)]) )
            // << " ) < min( "<<
            // (0.81 + (0.053 * (*phoEt)[denPCvint[sysb].at(k)]) )
            // <<", "<<
            // (0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<" )"<<std::endl;

            // std::cout<<"    PFneuIso:       max( 0, "<<
            // (*phoPFNeuIso)[denPCvint[sysb].at(k)] - ( rho * EAneutral((*phoSCEta)[denPCvint[sysb].at(k)]) )
            // << " ) < min( "<<
            // (1.92 + (0.014 * (*phoEt)[denPCvint[sysb].at(k)]) + (0.000019 * pow((*phoEt)[denPCvint[sysb].at(k)], 2.0)))
            // <<", "<<
            // (0.2 * ((*phoEt)[denPCvint[sysb].at(k)]))<<" )"<<std::endl<<std::endl;

           }
        FillDenHistograms(ptb, sysb, denPCvint[sysb].at(k), event_weight);
       }
      }
      if(  // do an inclusive pT plot from bins
         (phoEt->at(denPCvint[sysb].at(k)) > ptbins[0]) &&
         (phoEt->at(denPCvint[sysb].at(k)) < ptbins[inclptbin])
        ){
       FillDenHistograms(lastptbin-1, sysb, denPCvint[sysb].at(k), event_weight);
      }
      // and one fully inclusive in pT
      FillDenHistograms(lastptbin, sysb, denPCvint[sysb].at(k), event_weight);
     }
    }

   } // for each sysb in sysbins
// end fill histograms

  } //end if passes MonoPhoton triggers
 } //end loop through entries

 // write these histograms to file
  std::cout<<std::endl;
  std::cout<<"Total Passing Numerator: "<<nc<<"  Total Passing Denominator: "<<dc<<std::endl;
  std::cout<<"made it through, about to write"<<std::endl;

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
std::vector<int> postAnalyzer_QCD::pcPassSel(int sel, int sys, double phoPtLo, double phoPtHi, double phoEtaMax, Bool_t isMC){
  std::vector<int> tmpCand;
  tmpCand.clear();
  bool noncoll;
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

     // non collision backgrounds
     //noncoll = kTRUE;
     //noncoll = (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001; // isMC
     noncoll = fabs((*phoseedTimeFull5x5)[p]) < 3. && (*phomipTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001;

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
                     ((*phohasPixelSeed)[p] ==  0 ) //&&
                     //((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0102 )
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
Double_t postAnalyzer_QCD::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0; // 0.0456;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0; // 0.0500;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0; // 0.0340;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0; // 0.0383;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0; // 0.0339;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0; // 0.0303;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0; // 0.0240;

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
