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
  for(unsigned int sysb=0; sysb<lastsysbin; ++sysb){
    sysbinname = sysbinnames.at(sysb);

   if(isJet){ 
    phoCand =  getPhoJetCand(175., 1.4442, sysbinname);
   }
   else{
    phoCand = getPhoCand(175., 1.4442, sysbinname); 
   }
   if(phoCand.size()>0)
     {
      n_phokin++; //
      int candphotonindex = phoCand.at(0);

      Float_t phoPt = getPhotonPt(candphotonindex, sysbinname) ;

      event_weight=makeEventWeight(crossSec,lumi,nrEvents,phoPt,isMC,isZnnG,ewkZG,ewkWG,isEle,isJet,sysbinname);
  
      passTrig = askPassTrig(isMC);

      passShape = askPassShape(candphotonindex,isJet,isHalo,isSpike) ;

      passSeed = askPassSeed(candphotonindex,isEle);

      passSpike = askPassSpike(candphotonindex,isMC,isSpike);

      passNoncoll = askPassNonColl(candphotonindex,isSpike);

      passMIP = askPassMIP(candphotonindex,isHalo);

      passLepRej = askPassLepRej(candphotonindex);

      // met stuff
      if(sysbinname==""        ){ theMET=pfMET; }
      if(sysbinname=="_JERUp"  ){ theMET=pfMET_T1JERUp; }
      if(sysbinname=="_JERDown"){ theMET=pfMET_T1JERDo; }
      if(sysbinname=="_JESUp"  ){ theMET=pfMET_T1JESUp; }
      if(sysbinname=="_JESDown"){ theMET=pfMET_T1JESDo; }
      if(sysbinname=="_MESUp"  ){ theMET=pfMET_T1MESUp; }
      if(sysbinname=="_MESDown"){ theMET=pfMET_T1MESDo; }
      if(sysbinname=="_EESUp"  ){ theMET=pfMET_T1EESUp; }
      if(sysbinname=="_EESDown"){ theMET=pfMET_T1EESDo; }
      if(sysbinname=="_PESUp"  ){ theMET=pfMET_T1PESUp; }
      if(sysbinname=="_PESDown"){ theMET=pfMET_T1PESDo; }
      if(sysbinname=="_TESUp"  ){ theMET=pfMET_T1TESUp; }
      if(sysbinname=="_TESDown"){ theMET=pfMET_T1TESDo; }
      if(sysbinname=="_UESUp"  ){ theMET=pfMET_T1UESUp; }
      if(sysbinname=="_UESDown"){ theMET=pfMET_T1UESDo; }
      if(sysbinname=="_EWKUp"  ){ theMET=pfMET; }
      if(sysbinname=="_EWKDown"){ theMET=pfMET; }

      passMET = askPassMET(theMET,isMC);

      // dPhi( Jets, MET )
      passdPhiJM = askPassdPhiJM(candphotonindex,pfMETPhi);

      // dPhi( photon, MET )
      passdPhiPhoMET = askPassdPhiPhoMET(candphotonindex,pfMETPhi);

//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorApplication
//edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
//iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
//JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
//JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);
//for (PFJetCollection::const_iterator jet = jets->begin(); jet != jets->end(); jet++)  {
//  ...............................
//  jecUnc->setJetEta(eta);
//  jecUnc->setJetPt(ptCor); // here you must use the CORRECTED jet pt
//  double unc = jecUnc->getUncertainty(true);
//  double ptCor_shifted = ptCor(1+shift*unc) ; // shift = +1(up), or -1(down)
//}

      if(sysb==0){
       if(passTrig       ){ ++n_passTrig       ;}
       if(passShape      ){ ++n_passShape      ;}
       if(passSeed       ){ ++n_passSeed       ;}
       if(passSpike      ){ ++n_passSpike      ;}
       if(passNoncoll    ){ ++n_passNoncoll    ;}
       if(passMIP        ){ ++n_passMIP        ;}
       if(passLepRej     ){ ++n_passLepRej     ;}
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
      && passMET
      && passdPhiJM
      && passdPhiPhoMET)
      {
       if(sysb==0){ printf("%i:%i:%lli \n",run,lumis,event); 
        if(!isMC){
         std::cout<<"  pt: "<<phoPt<<"  eta: "<<phoEta->at(candphotonindex)<<std::endl;
         std::cout<<"  phi: "<<phoPhi->at(candphotonindex)<<"  met: "<<pfMET<<std::endl;
        }
       nc++;
       }
       callFillSigHist(sysb, lastptbin, inclptbin, candphotonindex, event_weight);
      } // end fill histograms
     } // end if phoCand[sysb][0].size()>0
    } // sysbin

 } //end loop through entries

 // write these histograms to file
  std::cout<<std::endl;
  std::cout<<"Total Passing RECO : "<<nc<<std::endl;
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

