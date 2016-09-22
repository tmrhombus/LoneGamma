#define postAnalyzer_Signal_cxx
#include "postAnalyzer_Signal.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void postAnalyzer_Signal::Loop(TString outfilename, Bool_t isMC, Double_t lumi, Double_t nrEvents, Double_t crossSec,
  Bool_t isZnnG, Bool_t isEle, Bool_t isHalo, Bool_t isSpike, Bool_t isJet, Bool_t ewkWG, Bool_t ewkZG)
{

 TStopwatch sw; 
 sw.Start();

 Int_t nc = 0;
 Int_t dc = 0;
 gencount = 0;

 int inclptbin = ptbins.size()-1;
 int lastptbin = ptbinnames.size()-1;
 int lastsysbin = sysbinnames.size();

 int n_initial=0;
 int n_phokin=0;

 int n_passTrig=0;
 int n_passShape=0;
 int n_passSeed=0;
 int n_passSpike=0;
 int n_passNoncoll=0;
 int n_passMIP=0;
 int n_passLepRej=0;
 int n_passMETfilters=0;
 int n_passMET=0;
 int n_passdPhiJM=0;
 int n_passdPhiPhoMET=0;

 int cf_passTrig=0;
 int cf_passShape=0;
 int cf_passSeed=0;
 int cf_passSpike=0;
 int cf_passNoncoll=0;
 int cf_passMIP=0;
 int cf_passLepRej=0;
 int cf_passMETfilters=0;
 int cf_passMET=0;
 int cf_passdPhiJM=0;
 int cf_passdPhiPhoMET=0;

  TFile *f_ewk_corr = new TFile("ewk_corr.root");
  TH1D *ewkZGCorrection = (TH1D*)f_ewk_corr->Get("zg");
  cout<<"ewkZG Correction histogram loaded"<<endl;
  TH1D *ewkWGCorrection = (TH1D*)f_ewk_corr->Get("wg");
  cout<<"ewkWG Correction histogram loaded"<<endl;

 if (fChain == 0) return;

 Long64_t nentries = fChain->GetEntriesFast();

 Long64_t nbytes = 0, nb = 0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
  n_initial++; // 
  if (jentry%10000 == 0)
    {
      std::cout<<"Starting entry "<<jentry<<"/"<<(nentries)<<" at "<<sw.RealTime()<<" RealTime, "<<sw.CpuTime()<<" CpuTime"<<std::endl;
      sw.Continue();
    }

  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   nbytes += nb;


  //// vector of ints, each int corresponds to location in vector of photons of photon passing cut
  //// one such vector for each systematic bin
  //for(unsigned int sysb=0; sysb<lastsysbin; ++sysb){
  // phoCand[sysb] = getPhoCand(175,1.4442);
  //}
  if(isJet){ phoCand = getPhoJetCand(175,1.4442); }
  else{ phoCand = getPhoCand(175,1.4442); }
   if(phoCand.size()>0)
     {
      n_phokin++; //
      int candphotonindex = phoCand.at(0);

      Float_t uncorrectedPhoEt = ((*phoSCRawE)[candphotonindex]/TMath::CosH((*phoSCEta)[candphotonindex]));
      double phopt = phoEt->at(candphotonindex) ;

      //=1.0 for real data
      event_weight=1.0;
      crossSecScl = crossSec;
      if(isZnnG){
       if      ( uncorrectedPhoEt < 190 ) {crossSecScl*=1.39;}  // {crossSecScl=14.4*1.39;}
       else if ( uncorrectedPhoEt < 250 ) {crossSecScl*=1.35;}  // {crossSecScl=29.7*1.35;}
       else if ( uncorrectedPhoEt < 400 ) {crossSecScl*=1.30;}  // {crossSecScl=17.5*1.30;}
       else if ( uncorrectedPhoEt < 700 ) {crossSecScl*=1.23;}  // {crossSecScl=3.7*1.23;}
       else                               {crossSecScl*=1.23;}  // {crossSecScl=0.3*1.23;}
      } 
      if(isMC){ event_weight=0.96*lumi*crossSecScl*(1.013 - 0.0001168*uncorrectedPhoEt)/nrEvents; }
      if(ewkZG){ 
       Double_t EWK_percent_adjustment = ewkZGCorrection->GetBinContent(ewkZGCorrection->GetXaxis()->FindBin(uncorrectedPhoEt));
       event_weight*=(1.0+.01*EWK_percent_adjustment) ; 
      }  
      if(ewkWG){ 
       Double_t EWK_percent_adjustment = ewkWGCorrection->GetBinContent(ewkWGCorrection->GetXaxis()->FindBin(uncorrectedPhoEt));
       event_weight*=(1.0+.01*EWK_percent_adjustment) ; 
      }  
  
      // if event passes MonoPhoton triggers (HLT_Photon165_HE10_v)
      // https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc#L179
      bool passTrig =( 
       ( (HLTPho>>7&1) == 1 ) ||
       ( (HLTPho>>8&1) == 1 ) ||
       ( (HLTPho>>9&1) == 1 ) ||
       ( (HLTPho>>10&1) == 1 ) ||
       ( (HLTPho>>11&1) == 1 ) ||
       ( (HLTPho>>12&1) == 1 )
      );
      //bool passTrig =( (HLTPho>>12&1) == 1);
      if(isMC){ passTrig=true; }

      bool passShape = phoSigmaIEtaIEtaFull5x5->at(candphotonindex)  <  0.0102;
      if(isJet){ passShape = true; }
      if(isHalo){ passShape = phoSigmaIEtaIEtaFull5x5->at(candphotonindex)  <  0.0165; }
      if(isSpike){ passShape = phoSigmaIEtaIEtaFull5x5->at(candphotonindex)  >  0.0102; }

      bool passSeed = phohasPixelSeed->at(candphotonindex) == 0;
      if(isEle){passSeed = phohasPixelSeed->at(candphotonindex) == 1; }

      // spike cleaning
      int iphi = 41; 
      int ieta = 5;
      //bool passSpike = true; // if isMC
      bool passSpike = !(phoIPhi->at(candphotonindex) == iphi && phoIEta->at(candphotonindex) == ieta) ;
      if(isSpike){ passSpike = true ; }

      // isMC
      //bool passNoncoll = (phoSigmaIEtaIEtaFull5x5->at(candphotonindex) > 0.001)
      //                && (phoSigmaIPhiIPhiFull5x5->at(candphotonindex) > 0.001);   

      //bool passNoncoll = (phoSigmaIEtaIEtaFull5x5->at(candphotonindex) > 0.001)
      //                && (phoSigmaIPhiIPhiFull5x5->at(candphotonindex) > 0.001)
      //                && (fabs(phoseedTimeFull5x5->at(candphotonindex)) < 3.);

      //  if(phoSigmaIPhiIPhiFull5x5->size() != phoseedTimeFull5x5->size()){
      //   std::cout<<phoSigmaIEtaIEtaFull5x5->size()<<std::endl;
      //   std::cout<<phoSigmaIPhiIPhiFull5x5->size()<<std::endl;
      //   std::cout<<phoseedTimeFull5x5->size()<<std::endl;
      //   std::cout<<" diff size run:lumis:event "
      //    <<run<<":"<<lumis<<":"<<event<< std::endl; 
      //  }
      bool passSeedTime = false;
      if( phoseedTimeFull5x5->size() > candphotonindex ){ passSeedTime = (fabs(phoseedTimeFull5x5->at(candphotonindex)) < 3.); }
      bool passNoncoll = ( (phoSigmaIEtaIEtaFull5x5->at(candphotonindex) > 0.001)
                        && (phoSigmaIPhiIPhiFull5x5->at(candphotonindex) > 0.001) 
                        && (phoR9->at(candphotonindex) < 1) 
                        && passSeedTime
                        );
                        //&& (fabs(phoseedTimeFull5x5->at(candphotonindex)) < 3.) );


      bool passMIP = phomipTotEnergy->at(candphotonindex) < 4.9;
      if(isHalo){ passMIP = phomipTotEnergy->at(candphotonindex) > 4.9; }

      // lepton rejection
      std::vector<int> elelist = electron_passLooseID(candphotonindex, 10.);
      std::vector<int> mulist = muon_passLooseID(candphotonindex, 10.);
      bool passLepRej = false;
      if( elelist.size()==0 && mulist.size()==0 ){ passLepRej = true; }


      // MET requirements
      bool passMETfilters = ( metFilters==0 ) ;
      if(isMC){ passMETfilters=true; }
      bool passMET = (pfMET > 170.) ;

      // dPhi( Jets, MET )
      std::vector<int>  jetindexvector = selectedJets(candphotonindex);
      bool passdPhiJM = passdphiJetMET(&jetindexvector, pfMETPhi);

      // dPhi( photon, MET )
      bool passdPhiPhoMET = ( DeltaPhi(phoPhi->at(candphotonindex),pfMETPhi)>2.0 ) ;

      if(passTrig       ){ ++n_passTrig       ;}
      if(passShape      ){ ++n_passShape      ;}
      if(passSeed       ){ ++n_passSeed       ;}
      if(passSpike      ){ ++n_passSpike      ;}
      if(passNoncoll    ){ ++n_passNoncoll    ;}
      if(passMIP        ){ ++n_passMIP        ;}
      if(passLepRej     ){ ++n_passLepRej     ;}
      if(passMETfilters ){ ++n_passMETfilters ;}
      if(passMET        ){ ++n_passMET        ;}
      if(passdPhiJM     ){ ++n_passdPhiJM     ;}
      if(passdPhiPhoMET ){ ++n_passdPhiPhoMET ;}

      if( passTrig ){ 
       ++cf_passTrig;
       if( passShape ){ 
        ++cf_passShape;
        if( passSeed ){ 
         ++cf_passSeed;
         if( passSpike ){ 
          ++cf_passSpike;
          if( passNoncoll ){ 
           ++cf_passNoncoll;
           if( passMIP ){ 
            ++cf_passMIP;
            if( passLepRej ){ 
             ++cf_passLepRej;
             if( passMETfilters ){ 
              ++cf_passMETfilters;
              if( passMET ){ 
               ++cf_passMET;
               if( passdPhiJM ){ 
                ++cf_passdPhiJM;
                if( passdPhiPhoMET ){ 
                 ++cf_passdPhiPhoMET;
                }
               }
              }
             }
            }
           }
          }
         }
        }
       }
      }

      // fill histograms
      if (passTrig
      && passShape
      && passSeed
      && passSpike
      && passNoncoll
      && passMIP
      && passLepRej
      && passMETfilters
      && passMET
      && passdPhiJM
      && passdPhiPhoMET){
          printf("%i:%i:%lli \n",run,lumis,event);
          if(!isMC){
           std::cout<<"  uncorret: "<<uncorrectedPhoEt<<"  et: "<<phopt<<"  eta: "<<phoEta->at(candphotonindex)<<std::endl;
           std::cout<<"  phi: "<<phoPhi->at(candphotonindex)<<"  met: "<<pfMET<<std::endl;
          }
          nc++;
          for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
           if(
              (phoEt->at(candphotonindex) > ptbins[ptb]) &&
              (phoEt->at(candphotonindex) < ptbins[ptb+1])
             ){
            FillSigHistograms(ptb, 0, candphotonindex, event_weight);
           } // end if passes pt cuts then fill
          } // end pt bin loop
          if(  // do an inclusive pT plot from bins
             (phoEt->at(candphotonindex) > ptbins[0]) &&
             (phoEt->at(candphotonindex) < ptbins[inclptbin])
            ){
           FillSigHistograms(lastptbin-1, 0, candphotonindex, event_weight);
          }
          // and one fully inclusive in pT
          FillSigHistograms(lastptbin, 0, candphotonindex, event_weight);
         }
      // end fill histograms
     } //end if phoCand[0].size()>0


////////////////////////////////////////////////////////////////////
//  // now do gen part
//    std::vector<int> genPhotonList;
//    std::vector<int> genMuonList;
//    std::vector<int> genEleList;
//    for( int a=0; a<nMC; ++a ){
//     if( mcPID->at(a)==22 && mcPt->at(a)>175. && abs(mcEta->at(a))<1.4442){ genPhotonList.push_back(a); }
//     if( abs(mcPID->at(a))==11 && mcPt->at(a)>10. && abs(mcEta->at(a))<2.5){ genEleList.push_back(a); }
//     if( abs(mcPID->at(a))==13 && mcPt->at(a)>10. && abs(mcEta->at(a))<2.5){ genMuonList.push_back(a); }
//     //printf(" mcPID: %i",mcPID->at(a));
//    }
//
//   int    leadingG;
//   double leadingGpt;
//   leadingG = 0;
//   leadingGpt = 0.;
//   for( int i=0; i<genPhotonList.size(); ++i ){ 
//    if (mcEt->at(genPhotonList.at(i)) > leadingGpt){
//     leadingGpt = mcEt->at(genPhotonList.at(i)) ;
//     leadingG = genPhotonList.at(i) ;
//    } 
//   }
//
//  if( genPhotonList.size()>0 ){
//
//   int candphotonindex = genPhotonList.at(0);
//
//   //printf(" photonindex: %i  pT: %.3f \n",candphotonindex,mcPt->at(candphotonindex));
//
//   fourVec_genl1.SetPtEtaPhiE(0,0,0,0);
//   fourVec_genl2.SetPtEtaPhiE(0,0,0,0);
//
//   // get electrons and muons and put into 4vectors
//   makeGenDilep(candphotonindex, genEleList, genMuonList, &fourVec_genl1, &fourVec_genl2, &fourVec_genee, &fourVec_genmm);
//   // make dilepton object and add to met
//   fourVec_genll = fourVec_genl1 + fourVec_genl2;
//   //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_l1.Pt(), fourVec_l1.Eta(), fourVec_l1.Phi(), fourVec_l1.M() );
//   //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_l2.Pt(), fourVec_l2.Eta(), fourVec_l2.Phi(), fourVec_l2.M() );
//   //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_ll.Pt(), fourVec_ll.Eta(), fourVec_ll.Phi(), fourVec_ll.M() );
//   gen_dilep_mass = fourVec_genll.M();
//
//   fourVec_genMET.SetPtEtaPhiE( genMET, 0., genMETPhi, 0. );
//
//   fourVec_genLeptoMET = fourVec_genMET + fourVec_genll;
//   genLeptoMET = fourVec_genLeptoMET.Pt();
//   genLeptoMEPhi = fourVec_genLeptoMET.Phi();
//
//  bool passGenMET = genLeptoMET > 140.;
//  bool passGendPhiPhoMET = ( DeltaPhi( mcPhi->at(candphotonindex), genLeptoMEPhi ) > 2.0 ) ;
//
//  // fill histograms
//  if ( passGenMET && passGendPhiPhoMET && gen_dilep_mass>60. && gen_dilep_mass<120.)
//     {
//      gencount++;
//      for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
//       if(
//          (mcEt->at(candphotonindex) > ptbins[ptb]) &&
//          (mcEt->at(candphotonindex) < ptbins[ptb+1])
//         ){
//        FillGenHistograms(ptb, 0, candphotonindex, event_weight);
//       } // end if passes pt cuts then fill
//      } // end pt bin loop
//      if(  // do an inclusive pT plot from bins
//         (mcEt->at(candphotonindex) > ptbins[0]) &&
//         (mcEt->at(candphotonindex) < ptbins[inclptbin])
//        ){
//       FillGenHistograms(lastptbin-1, 0, candphotonindex, event_weight);
//      }
//      // and one fully inclusive in pT
//      FillGenHistograms(lastptbin, 0, candphotonindex, event_weight);
//     }
//  // end fill histograms
//  }
/////////////////////////////////////////////////////////////////////


 } //end loop through entries

 // write these histograms to file
  std::cout<<std::endl;
  std::cout<<"Total Passing RECO : "<<nc<<std::endl;
  std::cout<<"Total Passing GEN  : "<<gencount<<std::endl;
  std::cout<<"made it through, about to write"<<std::endl;

  printf("lumi         %f\n", lumi);
  printf("nrEvents     %f\n", nrEvents);
  printf("crossSec     %f\n", crossSec);
  printf("isMC         %i\n", isMC);
  printf("isZnnG       %i\n", isZnnG);
  printf("isEle        %i\n", isEle);
  printf("isHalo       %i\n", isHalo);
  printf("isSpike      %i\n", isSpike);
  printf("isJet        %i\n", isJet);

  printf(" n_initial          %i\n\n",  n_initial);
  printf(" n_phokin           %i\n\n",  n_phokin);

  printf(" n_passTrig         %i\n", n_passTrig       );
  printf(" n_passShape        %i\n", n_passShape      );
  printf(" n_passSeed         %i\n", n_passSeed       );
  printf(" n_passSpike        %i\n", n_passSpike      );
  printf(" n_passNoncoll      %i\n", n_passNoncoll    );
  printf(" n_passMIP          %i\n", n_passMIP        );
  printf(" n_passLepRej       %i\n", n_passLepRej     );
  printf(" n_passMETfilters   %i\n", n_passMETfilters );
  printf(" n_passMET          %i\n", n_passMET        );
  printf(" n_passdPhiJM       %i\n", n_passdPhiJM     );
  printf(" n_passdPhiPhoMET   %i\n\n", n_passdPhiPhoMET );

  printf(" cf_passTrig        %i\n", cf_passTrig       );
  printf(" cf_passShape       %i\n", cf_passShape      );
  printf(" cf_passSeed        %i\n", cf_passSeed       );
  printf(" cf_passSpike       %i\n", cf_passSpike      );
  printf(" cf_passNoncoll     %i\n", cf_passNoncoll    );
  printf(" cf_passMIP         %i\n", cf_passMIP        );
  printf(" cf_passLepRej      %i\n", cf_passLepRej     );
  printf(" cf_passMETfilters  %i\n", cf_passMETfilters );
  printf(" cf_passMET         %i\n", cf_passMET        );
  printf(" cf_passdPhiJM      %i\n", cf_passdPhiJM     );
  printf(" cf_passdPhiPhoMET  %i\n", cf_passdPhiPhoMET );

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


//-------------------------getPhoCand 
std::vector<int> postAnalyzer_Signal::getPhoCand(double phoPtCut, double phoEtaCut){

  std::vector<int> pholist;
  pholist.clear();

  //Loop over photons
  for(int p=0;p<nPho;p++)
    {
      double uncorrectedPhoEt = ((*phoSCRawE)[p]/TMath::CosH((*phoSCEta)[p]));

      bool kinematic = uncorrectedPhoEt > phoPtCut && fabs((*phoSCEta)[p])<phoEtaCut;

      bool photonId = (
                       ((*phoHoverE)[p]                <  0.05   ) &&
                       ( TMath::Max( (*phoPFChIso)[p] - 0.0, 0.0) < 1.37 )  && 
                       ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 1.37 )  &&
                       ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) <
                        (1.06 + (0.014 * uncorrectedPhoEt) + (0.000019 * pow(uncorrectedPhoEt, 2.0))) )  &&
                       ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < 
                        (0.28 + (0.0053 * uncorrectedPhoEt)) ) 
                      );

      if(photonId && kinematic){
        pholist.push_back(p);
      } 
    }  

  return pholist;

}

//-------------------------getPhoJetCand 
std::vector<int> postAnalyzer_Signal::getPhoJetCand(double phoPtCut, double phoEtaCut){

  std::vector<int> pholist;
  pholist.clear();

  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++)
    {
      Float_t uncorrectedPhoEt = ((*phoSCRawE)[p]/TMath::CosH((*phoSCEta)[p]));

      bool kinematic = uncorrectedPhoEt > phoPtCut && fabs((*phoSCEta)[p])<phoEtaCut;

      bool photonId = (
                       ((*phoHoverE)[p]                <  0.05   ) &&
                       ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0102 ) &&
                       ((*phohasPixelSeed)[p]              ==  0      ) &&
                       ( TMath::Max( (*phoPFChIso)[p] - 0.0, 0.0) < 1.37 )  && // wtf TMath - shouldn't do anything since we have worstCHiso ..
                       ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 1.37 )  &&
                       ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) <
                        (1.06 + (0.014 * uncorrectedPhoEt) + (0.000019 * pow(uncorrectedPhoEt, 2.0))) )  &&
                       ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < 
                        (0.28 + (0.0053 * uncorrectedPhoEt)) ) 
                      );


     // denominator selections

     bool passHoE = ( (*phoHoverE)[p] < 0.05 );

     double vloosePFCharged= TMath::Min(5.0*(3.32) , 0.20 *  uncorrectedPhoEt);
     double vloosePFPhoton = TMath::Min(5.0*(0.81+ (0.0053 * uncorrectedPhoEt)) , 0.20*uncorrectedPhoEt);
     double vloosePFNeutral= TMath::Min(5.0*(1.92 + (0.014 * uncorrectedPhoEt) + (0.000019 * pow(uncorrectedPhoEt, 2.0))) , 0.20*uncorrectedPhoEt);
     bool passVLooseIso = ( 
                      ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < vloosePFCharged )  &&  
                      ( TMath::Max( ( (*phoPFChIso)[p] - 0.0 ), 0.0) < vloosePFCharged )  &&  
                      ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < vloosePFNeutral )  &&  
                      ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < vloosePFPhoton )
                     );  
     bool passLoosePIso = 
                     ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) <
                       (0.81 + (0.0053 * uncorrectedPhoEt)) );
     bool passLooseIso = (  // deno must fail this cut
                     ( TMath::Max( ( (*phoPFChIso)[p] - 0.0 ), 0.0) < 3.32 )  &&  
                     ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 3.32 )  &&  
                     ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) <
                       (1.92 + (0.014* uncorrectedPhoEt) + (0.000019 * pow(uncorrectedPhoEt, 2.0))))  &&  
                     passLoosePIso
                    );  

     bool photonID = !passLooseIso && passVLooseIso ;
      if(photonID && kinematic){
        pholist.push_back(p);
      } 
    }  

  return pholist;

}

//-------------------------selectedJets
std::vector<int> postAnalyzer_Signal::selectedJets(int pho_index) {

  bool jetVeto=true;
  std::vector<int> jetindex;
  float value = 0.0;

  for(int i = 0; i < nJet; i++)
    {
      if(0.0 < abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.5)   {value =-0.8;}
      if(2.5 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.75) {value =-0.95;}
      if(2.75 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <3.0) {value =-0.97;}
      if(3.00 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <5.0) {value =-0.99;}

      double deltar = 0.0 ;
      //      std::cout<<"Jet size: "<<nJet<<std::endl;
      //      std::cout<<"Jet no:"<<i<<"coming here pujetid: "<<pfJet_pt[i]<<std::endl;
      //      if(OverlapWithMuon(jetEta->at(i),jetPhi->at(i)))     continue;
      //      std::cout<<"Jet no:"<<i<<"coming here OverlapWithMuon: "<<pfJet_pt[i]<<std::endl;
      //      if(OverlapWithElectron(jetEta->at(i),jetPhi->at(i)))   continue;
      if(pho_index>=0){
        deltar= dR(jetEta->at(i),jetPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index));
        //std::cout<<"Delta R between photon and jet="<<dR<< "jet pt"<<pfJet_pt[i]<<std::endl;
      }
      if(deltar>0.4 && jetPt->at(i) >30.0 && jetPFLooseId->at(i)==1 && jetPUidFullDiscriminant->at(i)>value)
        {
          jetindex.push_back(i);
        }
    }

  //  std::cout<<"Jet size: "<< jetindex.size()<<std::endl;
  //if(jetindex.size()>1)jetVeto = false;
  return jetindex;

}


//-------------------------dR
double postAnalyzer_Signal::dR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

//-------------------------DeltaPhi
double postAnalyzer_Signal::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}

//-------------------------passdphiJetMET
bool postAnalyzer_Signal::passdphiJetMET(std::vector<int> *jets, double mephi)
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



////-------------------------makeGenDilep
//void postAnalyzer_Signal::makeGenDilep(int pho_index, std::vector<int> elelist, std::vector<int> mulist, TLorentzVector *fv_1, TLorentzVector *fv_2, TLorentzVector *fv_ee, TLorentzVector *fv_mm)
//{
//
//  TLorentzVector e1, e2, ee;
//  TLorentzVector m1, m2, mm;
//
//  int    leadingE, leadingM;
//  leadingE = 0;
//  leadingM = 0;
//  double leadingEpt, leadingMpt;
//  leadingEpt = 0.;
//  leadingMpt = 0.;
//
//  for( int i=0; i<mulist.size(); ++i ){ 
//   if (mcEt->at(mulist.at(i)) > leadingMpt){
//    leadingMpt = mcEt->at(mulist.at(i)) ;
//    leadingM = mulist.at(i) ;
//   } 
//  }
//
//  for( int i=0; i<elelist.size(); ++i ){ 
//   if (mcEt->at(elelist.at(i)) > leadingEpt){
//    leadingEpt = mcEt->at(elelist.at(i)) ;
//    leadingE = elelist.at(i) ;
//   } 
//  }
//
//
//  // no pairs
//  if( elelist.size()<2 && mulist.size()<2 ){return;}
//
//  // electrons
//  if( elelist.size()>1 ){
//   for(int i=0; i<elelist.size(); ++i)
//   {
//     if( mcPID->at(leadingE) + mcPID->at(elelist[i])==0 )
//     {
//      //printf(" --we have electrons ");
//      e1.SetPtEtaPhiE( mcPt->at(leadingE), mcEta->at(leadingE), mcPhi->at(leadingE), mcE->at(leadingE) );
//      e2.SetPtEtaPhiE( mcPt->at(elelist[i]), mcEta->at(elelist[i]), mcPhi->at(elelist[i]), mcE->at(elelist[i]) );
//      break;
//     }
//   }
//   ee = e1 + e2;
//   //printf(": dilep mass = %f", ee.M());
//  }
//
//  // muons
//  if( mulist.size()>1 ){
//   for(int i=0; i<mulist.size(); ++i)
//   {
//     if( mcPID->at(leadingM) + mcPID->at(mulist[i])==0 )
//     {
//      //printf(" --we have muons ");
//      m1.SetPtEtaPhiE( mcPt->at(leadingM), mcEta->at(leadingM), mcPhi->at(leadingM), mcE->at(leadingM) );
//      m2.SetPtEtaPhiE( mcPt->at(mulist[i]), mcEta->at(mulist[i]), mcPhi->at(mulist[i]), mcE->at(mulist[i]) );
//      break;
//     }
//   }
//   mm = m1 + m2;
//   //printf(": dilep mass = %f", mm.M());
//  }
//
//  *fv_ee = ee;
//  *fv_mm = mm;
//  // take highest mass dilepton pair
//  if( ee.M()>60. && ee.M()<120. ){ *fv_1 = e1; *fv_2 = e2; }
//  else if( mm.M()>60. && mm.M()<120. ){ *fv_1 = m1; *fv_2 = m2; }
//  //if( ee.M() > mm.M() && ee.M()>60. && ee.M()<120. ){ printf(": using ee\n"); *fv_1 = e1; *fv_2 = e2; }
//  //else if( ee.M() < mm.M() && mm.M()>60. && mm.M()<120. ){ printf(": using mm\n"); *fv_1 = m1; *fv_2 = m2; }
//
//  return;
//  
//}


//-------------------------makeDilep
void postAnalyzer_Signal::makeDilep(int pho_index, TLorentzVector *fv_1, TLorentzVector *fv_2, TLorentzVector *fv_ee, TLorentzVector *fv_mm)
{
  
  std::vector<int> elelist = electron_passLooseID(pho_index, 10.);
  std::vector<int> mulist = muon_passLooseID(pho_index, 10.);

  TLorentzVector e1, e2, ee;
  TLorentzVector m1, m2, mm;

  // no pairs
  if( elelist.size()<2 && mulist.size()<2 ){return;}

  // electrons
  if( elelist.size()>1 ){
   for(int i=1; i<elelist.size(); ++i)
   {
     if( eleCharge->at(0)*eleCharge->at(i)==-1 )
     //if( eleCharge[0]*eleCharge[i]==-1 )
     {
      //printf(" --we have electrons ");
      e1.SetPtEtaPhiE( elePt->at(elelist[0]), eleEta->at(elelist[0]), elePhi->at(elelist[0]), eleEn->at(elelist[0]) );
      e2.SetPtEtaPhiE( elePt->at(elelist[i]), eleEta->at(elelist[i]), elePhi->at(elelist[i]), eleEn->at(elelist[i]) );
      break;
     }
   }
   ee = e1 + e2;
   //printf(": dilep mass = %f", ee.M());
  }
  // muons
  if( mulist.size()>1 ){
   for(int i=1; i<mulist.size(); ++i)
   {
     if( muCharge->at(0)*muCharge->at(i)==-1 )
     {
      //printf(" --we have muons ");
      m1.SetPtEtaPhiE( muPt->at(mulist[0]), muEta->at(mulist[0]), muPhi->at(mulist[0]), muEn->at(mulist[0]) );
      //printf(": charge pass ");
      m2.SetPtEtaPhiE( muPt->at(mulist[i]), muEta->at(mulist[i]), muPhi->at(mulist[i]), muEn->at(mulist[i]) );
      break;
     }
   }
   mm = m1 + m2;
   //printf(": dilep mass = %f", mm.M());
  }

   *fv_ee = ee;
   *fv_mm = mm;
  // take highest mass dilepton pair
  if( mm.M()>60. && mm.M()<120. ){ *fv_1 = m1; *fv_2 = m2; }
  else if( ee.M()>60. && ee.M()<120. ){ *fv_1 = e1; *fv_2 = e2; }
  ///if( ee.M() > mm.M() && ee.M()>60. && ee.M()<120. ){ *fv_1 = e1; *fv_2 = e2; }
  ///else if( ee.M() < mm.M() && mm.M()>60. && mm.M()<120. ){ *fv_1 = m1; *fv_2 = m2; }
  //if( ee.M() > mm.M() && ee.M()>60. && ee.M()<120. ){ printf(": using ee\n"); *fv_1 = e1; *fv_2 = e2; }
  //else if( ee.M() < mm.M() && mm.M()>60. && mm.M()<120. ){ printf(": using mm\n"); *fv_1 = m1; *fv_2 = m2; }

  return;
  
}


//-------------------------electron_passLooseID
std::vector<int> postAnalyzer_Signal::electron_passLooseID(int pho_index, float elePtCut)
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
std::vector<int> postAnalyzer_Signal::muon_passLooseID(int pho_index, float muPtCut)
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
      //Muon passes Loose Muon ID and PF-based combined relative, dBeta-corrected Loose Muon Isolation cuts  

//if( event==767326116 ){
//      printf(" Event: %lli \n",event);
//      printf("  muPt              %f \n",  muPt->at(i)                        );
//      printf("  muIsPFMuon        %f \n",  muIsPFMuon->at(i)                  );
//      printf("  muIsGlobalMuon    %f \n",  muIsGlobalMuon->at(i)              );
//      printf("  muIsTrackerMuon   %f \n",  muIsTrackerMuon->at(i)             );
//      printf("  muPFNeuIso        %f \n",  muPFNeuIso->at(i)                  );
//      printf("  muPFPhoIso        %f \n",  muPFPhoIso->at(i)                  );
//      printf("  muPFPUIso         %f \n",  muPFPUIso->at(i)                   );
//      printf("  muPhoPU           %f \n",  muPhoPU                            );
//      printf("\n\n");
//                   std::cout<<muPt->at(i)<<std::endl;
//                   std::cout<<muEta->at(i)<<std::endl;
//                   std::cout<<"veto Passed!!!!"<<muEta->at(i)<<" "<<muPhi->at(i)<<" "<<phoSCEta->at(pho_index)<<" "<<phoSCPhi->at(pho_index)<<std::endl;
//}

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

////Veto failed if a muon is found that passes Loose Muon ID, Loose Muon Isolation, and muPtcut, and does not overlap the candidate photon within dR of 0.5                                                                                                                                                                    
//bool postAnalyzer::muon_veto_looseID(int pho_index, float muPtCut)
//{
//  bool veto_passed = true; //pass veto if no good muon found                                                                                                                                                                                                                                                           
//  bool pass_PFMuon = false;
//  bool pass_globalMuon = false;
//  bool pass_trackerMuon = false;
//  bool pass_iso = false;
//  //Explicitly stating types to avoid a TMath::Max conversion issue                                                                                                                                                                                                                                                    
//  Float_t zero = 0.0;
//  Float_t muPhoPU = 999.9;
//  Float_t tightIso_combinedRelative = 999.9;
//  for(int i = 0; i < nMu; i++)
//    {   
//      pass_PFMuon = muIsPFMuon->at(i);
//      pass_globalMuon = muIsGlobalMuon->at(i);
//      pass_trackerMuon = muIsTrackerMuon->at(i);
//      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
//      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
//      pass_iso = tightIso_combinedRelative < 0.25;
//      //Muon passes Loose Muon ID and PF-based combined relative, dBeta-corrected Loose Muon Isolation cuts                                                                                                                                                                                                        
//      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon))
//        {   
//          //Muon passes pt cut                                                                                                                                                                                                                                                                                 
//          if(muPt->at(i) > muPtCut)
//            {   
//              //Muon does not overlap photon                                                                                                                                                                                                                                                               
//              if(dR(muEta->at(i),muPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
//                {   
//                  veto_passed = false;
//                  break;
//                }   
//            }   
//        }   
//    }   
//  return veto_passed;
//}




////-------------------------OverlapWithMuon
//bool postAnalyzer_Signal::OverlapWithMuon(double eta, double phi){
//
//  bool overlap = false;
//  //  std::cout<<"No of muon:"<<Muon_n<<std::endl;
//  bool pass_PFMuon = false;
//  bool pass_globalMuon = false;
//  bool pass_trackerMuon = false;
//  bool pass_iso = false;
//
//  Float_t zero1 = 0.0;
//  Float_t muPhoPU = 999.9;
//  Float_t tightIso_combinedRelative = 999.9;
//  for(int i = 0; i < nMu; i++)
//    {
//      pass_PFMuon = muIsPFMuon->at(i);
//      pass_globalMuon = muIsGlobalMuon->at(i);
//      pass_trackerMuon = muIsTrackerMuon->at(i);
//      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
//      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero1,muPhoPU))/(muPt->at(i));
//      pass_iso = tightIso_combinedRelative < 0.25;
//      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
//        {
//          if(muPt->at(i) > 10.)
//            {
//              if(dR(muEta->at(i),muPhi->at(i),eta,phi) < 0.5)
//                {
//                  overlap = true;
//                  break;
//                }
//            }
//        }
//    }
//
//  return overlap;
//
//
//}
//
//
////-------------------------OverlapWithElectron
//bool postAnalyzer_Signal::OverlapWithElectron(double eta, double phi){
//  bool overlap = false;
//
//  bool pass_SigmaIEtaIEtaFull5x5 = false;
//  bool pass_dEtaIn = false;
//  bool pass_dPhiIn = false;
//  bool pass_HoverE = false;
//  bool pass_iso = false;
//  bool pass_ooEmooP = false;
//  bool pass_d0 = false;
//  bool pass_dz = false;
//  bool pass_missingHits = false;
//  bool pass_convVeto = false;
//
//  Float_t EA = 0.0;
//  Float_t zero = 0.0;
//  Float_t EAcorrIso = 999.9;
//  for(int i = 0; i < nEle; i++)
//    {
//      pass_SigmaIEtaIEtaFull5x5 = false;
//      pass_dEtaIn = false;
//      pass_dPhiIn = false;
//      pass_HoverE = false;
//      pass_iso = false;
//      pass_ooEmooP = false;
//      pass_d0 = false;
//      pass_dz = false;
//      pass_missingHits = false;
//      pass_convVeto = false;
//
//      if(abs(eleSCEta->at(i)) <= 1.0)
//        EA = 0.1752;
//      else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
//        EA = 0.1862;
//      else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
//        EA = 0.1411;
//      else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
//        EA = 0.1534;
//      else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
//        EA = 0.1903;
//      else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
//        EA = 0.2243;
//      else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
//        EA = 0.2687;
//      EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));
//
//      if(abs(eleSCEta->at(i)) <= 1.479)
//        {
//          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
//          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
//          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
//          pass_HoverE = eleHoverE->at(i) < 0.104;
//          pass_iso = EAcorrIso < 0.0893;
//          pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
//          pass_d0 = abs(eleD0->at(i)) < 0.0261;
//          pass_dz = abs(eleDz->at(i)) < 0.41;
//          pass_missingHits = eleMissHits->at(i) <= 2;
//          pass_convVeto = eleConvVeto->at(i) == 1;
//        }
//      else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
//        {
//          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
//          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
//          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
//          pass_HoverE = eleHoverE->at(i) < 0.0897;
//          pass_iso = EAcorrIso < 0.121;
//          pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
//          pass_d0 = abs(eleD0->at(i)) < 0.118;
//          pass_dz = abs(eleDz->at(i)) < 0.822;
//          pass_missingHits = eleMissHits->at(i) <= 1;
//          pass_convVeto = eleConvVeto->at(i) == 1;
//        }
//      //Electron passes Loose Electron ID cuts  
//      if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
//        {
//          //Electron passes pt cut 
//          if(elePt->at(i) > 10.)
//            {
//              //Electron does not overlap photon  
//              if(dR(eleSCEta->at(i),eleSCPhi->at(i),eta,phi) < 0.5)
//                {
//                  overlap = true;
//                  break;
//                }
//            }
//        }
//    }
//}


// Effective area to be needed in PF Iso for photon ID
// https://indico.cern.ch/event/455258/contribution/0/attachments/1173322/1695132/SP15_253rd.pdf -- slide-5
// Worst Charged Hadron Isolation EA
Double_t postAnalyzer_Signal::EAchargedworst(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.078;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.089;
  return EffectiveArea;
}

// standard EA
Double_t postAnalyzer_Signal::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0;

  return EffectiveArea;
}

Double_t postAnalyzer_Signal::EAneutral(Double_t eta){
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

Double_t postAnalyzer_Signal::EAphoton(Double_t eta){
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
