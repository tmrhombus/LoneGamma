#define analyzeSignal_cxx
#include "analyzeSignal.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void analyzeSignal::Loop(TString outfilename, Bool_t isMC, Double_t lumi, Double_t nrEvents, Double_t crossSec,
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
      if(isMC){ passSpike = true ; }

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

          callFillSigHist(0, lastptbin, inclptbin, candphotonindex, event_weight);

         }  // end fill histograms
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

