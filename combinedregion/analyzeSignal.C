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

//  TFile *f_ewk_corr = new TFile("ewk_corr.root");
//  TH1D *ewkZGCorrection = (TH1D*)f_ewk_corr->Get("zg");
//  cout<<"ewkZG Correction histogram loaded"<<endl;
//  TH1D *ewkWGCorrection = (TH1D*)f_ewk_corr->Get("wg");
//  cout<<"ewkWG Correction histogram loaded"<<endl;

  f_ewk_corr = new TFile("ewk_corr.root");
  ewkZGCorrection = (TH1D*)f_ewk_corr->Get("zg");
  cout<<"ewkZG Correction histogram loaded"<<endl;
  ewkWGCorrection = (TH1D*)f_ewk_corr->Get("wg");
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
  if(isJet){ phoCand = getPhoJetCand(175., 1.4442); }
  else{ phoCand = getPhoCand(175., 1.4442); }
   if(phoCand.size()>0)
     {
      n_phokin++; //
      int candphotonindex = phoCand.at(0);

      Double_t photonpt = ((*phoSCRawE)[candphotonindex]/TMath::CosH((*phoSCEta)[candphotonindex]));
      //Float_t uncorrectedPhoEt = ((*phoSCRawE)[candphotonindex]/TMath::CosH((*phoSCEta)[candphotonindex]));
      double phopt = phoEt->at(candphotonindex) ;

      event_weight=makeEventWeight(crossSec,lumi,nrEvents,photonpt,isMC,isZnnG,ewkZG,ewkWG,isEle,isJet);
  
      passTrig = askPassTrig(isMC);

      passShape = askPassShape(candphotonindex,isJet,isHalo,isSpike) ;
      //passShape = phoSigmaIEtaIEtaFull5x5->at(candphotonindex)  <  0.0102;
      //if(isJet){ passShape = true; }
      //if(isHalo){ passShape = phoSigmaIEtaIEtaFull5x5->at(candphotonindex)  <  0.0165; }
      //if(isSpike){ passShape = phoSigmaIEtaIEtaFull5x5->at(candphotonindex)  <  0.001; }

      passSeed = askPassSeed(candphotonindex,isEle);
      //passSeed = phohasPixelSeed->at(candphotonindex) == 0;
      //if(isEle){passSeed = phohasPixelSeed->at(candphotonindex) == 1; }

      passSpike = askPassSpike(candphotonindex,isMC,isSpike);
      //// spike cleaning
      //int iphi = 41; 
      //int ieta = 5;
      ////bool passSpike = true; // if isMC
      //passSpike = !(phoIPhi->at(candphotonindex) == iphi && phoIEta->at(candphotonindex) == ieta) ;
      //if(isSpike){ passSpike = true ; }
      //if(isMC){ passSpike = true ; }

      passNoncoll = askPassNonColl(candphotonindex,isSpike);
      //bool passSeedTime = false;
      //if( phoseedTimeFull5x5->size() > candphotonindex ){ passSeedTime = (fabs(phoseedTimeFull5x5->at(candphotonindex)) < 3.); }
      //bool passSpikeShape = ( (phoSigmaIEtaIEtaFull5x5->at(candphotonindex) > 0.001)
      //                     && (phoSigmaIPhiIPhiFull5x5->at(candphotonindex) > 0.001) );
      //if(isSpike){ passSpikeShape = ! ( (phoSigmaIEtaIEtaFull5x5->at(candphotonindex) > 0.001)
      //                     && (phoSigmaIPhiIPhiFull5x5->at(candphotonindex) > 0.001) ); }
      //passNoncoll = ( passSpikeShape 
      //                  && (phoR9->at(candphotonindex) < 1) 
      //                  && passSeedTime
      //                  );
      //                  //&& (fabs(phoseedTimeFull5x5->at(candphotonindex)) < 3.) );


      passMIP = askPassMIP(candphotonindex,isHalo);
      //passMIP = phomipTotEnergy->at(candphotonindex) < 4.9;
      //if(isHalo){ passMIP = phomipTotEnergy->at(candphotonindex) > 4.9; }

      passLepRej = askPassLepRej(candphotonindex);
      //// lepton rejection
      //std::vector<int> elelist = electron_passLooseID(candphotonindex, 10.);
      //std::vector<int> mulist = muon_passLooseID(candphotonindex, 10.);
      //passLepRej = false;
      //if( elelist.size()==0 && mulist.size()==0 ){ passLepRej = true; }


      passMET = askPassMET(pfMET,isMC);
      //// MET SELECTIONS
      //passMETfilters = ( metFilters==0 ) ;
      //if(isMC){ passMETfilters = true ; }
      //passMET = (pfMET > 170.) ;

      // dPhi( Jets, MET )
      passdPhiJM = askPassdPhiJM(candphotonindex,pfMETPhi);
      //std::vector<int>  jetindexvector = selectedJets(candphotonindex);
      //passdPhiJM = passdphiJetMET(&jetindexvector, pfMETPhi);

      // dPhi( photon, MET )
      passdPhiPhoMET = askPassdPhiPhoMET(candphotonindex,pfMETPhi);
      //passdPhiPhoMET = ( DeltaPhi(phoPhi->at(candphotonindex),pfMETPhi)>2.0 ) ;

      if(passTrig       ){ ++n_passTrig       ;}
      if(passShape      ){ ++n_passShape      ;}
      if(passSeed       ){ ++n_passSeed       ;}
      if(passSpike      ){ ++n_passSpike      ;}
      if(passNoncoll    ){ ++n_passNoncoll    ;}
      if(passMIP        ){ ++n_passMIP        ;}
      if(passLepRej     ){ ++n_passLepRej     ;}
      //if(passMETfilters ){ ++n_passMETfilters ;}
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
             //if( passMETfilters ){ 
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
             //}
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
      //&& passMETfilters
      && passMET
      && passdPhiJM
      && passdPhiPhoMET)
      {
       printf("%i:%i:%lli \n",run,lumis,event);
       if(!isMC){
        std::cout<<"  pt: "<<photonpt<<"  eta: "<<phoEta->at(candphotonindex)<<std::endl;
        std::cout<<"  phi: "<<phoPhi->at(candphotonindex)<<"  met: "<<pfMET<<std::endl;
       }
       callFillSigHist(0, lastptbin, inclptbin, candphotonindex, event_weight);

       nc++;
      }
      // end fill histograms
     } //end if phoCand[0].size()>0

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

