#define analyzePDFscaleWlnG_cxx
#include "analyzePDFscaleWlnG.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void analyzePDFscaleWlnG::Loop(TString outfilename, Bool_t isMC, Double_t lumi, Double_t nrEvents, Double_t crossSec,
  Bool_t isZnnG, Bool_t isEle, Bool_t isHalo, Bool_t isSpike, Bool_t isJet, Bool_t ewkWG, Bool_t ewkZG)
{

 TStopwatch sw; 
 sw.Start();

 Int_t nc = 0;
 Int_t dc = 0;

 nrec = 0;

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
 int n_passZWindow=0;
 int n_passdimass20=0;

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
 int cf_passZWindow=0;
 int cf_passdimass20=0;


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

      // get electrons and muons and put into 4vectors
      bool passM = false;
      bool passE = false;
      bool passL = false;

      tightEles = electron_passTightID(candphotonindex,30.0);
      looseEles = electron_passLooseID(candphotonindex,10.0);
      tightMus  = muon_passTightID(candphotonindex,30.0);
      looseMus  = muon_passLooseID(candphotonindex,10.0);

      passE =( tightEles.size() == 1 && looseEles.size() == 1 && looseMus.size() == 0  ) ;  
      passM =( tightMus.size() == 1  && looseMus.size() == 1  && looseEles.size() == 0 ) ;

      fourVec_e.SetPtEtaPhiE(0,0,0,0);
      fourVec_m.SetPtEtaPhiE(0,0,0,0);
      if(passE){
       //fourVec_e.SetPtEtaPhiE(elePt->at(tightEles[0]),eleEta->at(tightEles[0]),elePhi->at(tightEles[0]),eleEn->at(tightEles[0]));
       fourVec_e.SetPtEtaPhiE( getElectronPt(tightEles[0],sysbinname),
                               eleEta->at(tightEles[0]),
                               elePhi->at(tightEles[0]),
                               getElectronPt(tightEles[0],sysbinname)*TMath::CosH( eleEta->at(tightEles[0]) ) );
      }    

      if(passM){
       //fourVec_m.SetPtEtaPhiE(muPt->at(tightMus[0]),muEta->at(tightMus[0]),muPhi->at(tightMus[0]),muEn->at(tightMus[0]));
       fourVec_m.SetPtEtaPhiE( getMuonPt(tightMus[0],sysbinname),
                               muEta->at(tightMus[0]),
                               muPhi->at(tightMus[0]),
                               getMuonPt(tightMus[0],sysbinname)*TMath::CosH( muEta->at(tightMus[0]) ) );
      }    
      if(passE || passM){passL=true;}

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
      // Lepto MET Creation
      fourVec_met.SetPtEtaPhiE(theMET, 0., pfMETPhi, theMET);
      if(passE){ fourVec_leptomet = fourVec_met + fourVec_e; }
      if(passM){ fourVec_leptomet = fourVec_met + fourVec_m; }

      leptoMET = fourVec_leptomet.Pt();
      leptoMEPhi = fourVec_leptomet.Phi();

      // MET SELECTIONS
      passMET = askPassMET(leptoMET,isMC);

      // dPhi( Jets, MET )
      passdPhiJM = askPassdPhiJM(candphotonindex,leptoMEPhi);

      // dPhi( photon, MET )
      passdPhiPhoMET = askPassdPhiPhoMET(candphotonindex,leptoMEPhi);

      theMET = leptoMET;

      if(sysb==0){
      if(passTrig       ){ ++n_passTrig       ;}
      if(passShape      ){ ++n_passShape      ;}
      if(passSeed       ){ ++n_passSeed       ;}
      if(passSpike      ){ ++n_passSpike      ;}
      if(passNoncoll    ){ ++n_passNoncoll    ;}
      if(passMIP        ){ ++n_passMIP        ;}
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

       //fill histograms
      if (passTrig
      && passShape
      && passSeed
      && passSpike
      && passNoncoll
      && passMIP
      && passMET
      //&& passdPhiJM
      && passdPhiPhoMET
      && passL)
      {
       if(sysb==0){ printf("%i:%i:%lli \n",run,lumis,event); 
        if(!isMC){
         std::cout<<"  pt: "<<phoPt<<"  eta: "<<phoEta->at(candphotonindex)<<std::endl;
         std::cout<<"  phi: "<<phoPhi->at(candphotonindex)<<"  met: "<<pfMET<<std::endl;
        }
       nc++;
       callFillSigHistPDFscale(sysb, lastptbin, inclptbin, candphotonindex, event_weight);
       }
       callFillSigHist(sysb, lastptbin, inclptbin, candphotonindex, event_weight);
       callFillSigHistLep(sysb, lastptbin, inclptbin, candphotonindex, event_weight, passM);
      }
      // end fill histograms
     } //end if phoCand[0].size()>0
   }
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
  printf(" cf_passMETfilters  %i\n", cf_passMETfilters );
  printf(" cf_passMET         %i\n", cf_passMET        );
  printf(" cf_passdPhiJM      %i\n", cf_passdPhiJM     );
  printf(" cf_passdPhiPhoMET  %i\n", cf_passdPhiPhoMET );
  printf(" cf_passZWindow     %i\n", cf_passZWindow    ); 
  printf(" cf_passdimass20    %i\n", cf_passdimass20   ); 

  for(unsigned int i=0; i<ptbinnames.size(); ++i){ 

   double sum_scales = std::accumulate(vec_scales[i].begin(), vec_scales[i].end(), 0.0);
   double mean_scales = sum_scales / vec_scales[i].size();
   double sq_sum_scales = std::inner_product(vec_scales[i].begin(), vec_scales[i].end(), vec_scales[i].begin(), 0.0);
   double stdev_scales = std::sqrt(sq_sum_scales / vec_scales[i].size() - mean_scales * mean_scales);
   
   double sum_pdfs = std::accumulate(vec_pdfs[i].begin(), vec_pdfs[i].end(), 0.0);
   double mean_pdfs = sum_pdfs / vec_pdfs[i].size();
   double sq_sum_pdfs = std::inner_product(vec_pdfs[i].begin(), vec_pdfs[i].end(), vec_pdfs[i].begin(), 0.0);
   double stdev_pdfs = std::sqrt(sq_sum_pdfs / vec_pdfs[i].size() - mean_pdfs * mean_pdfs);
   
   std::cout<<"::::: "<<ptbinnames[i]<<" :::::"<<std::endl;
   std::cout<<std::endl<<"PDFs"<<std::endl;
   std::cout<<"   Mean = "<<mean_pdfs<<" std: "<<stdev_pdfs<<" = "<<stdev_pdfs/mean_pdfs<<std::endl;
   std::cout<<std::endl<<"scales"<<std::endl;
   std::cout<<"   Mean = "<<mean_scales<<" std: "<<stdev_scales<<" = "<<stdev_scales/mean_scales<<std::endl<<std::endl;
   
//   for(int i=0; i<vec_pdfs[i].size(); ++i)
//     std::cout << vec_pdfs[i].at(i) << ' ';
//   
//   for(int i=0; i<vec_scales[i].size(); ++i)
//     std::cout << vec_scales[i].at(i) << ' '<<std::endl<<std::endl;;
 
  }

 TFile *outfile = new TFile(outfilename,"RECREATE");
 outfile->cd();
 for(unsigned int i=0; i<ptbinnames.size(); ++i){
  for(unsigned int j=0; j<sysbinnames.size(); ++j){
   WriteHistograms(i,j);
   WriteHistogramsLep(i,j);
  }
  WriteHistogramsPDFscale(i);
 }
 outfile->Close();
 sw.Stop();
 std::cout<<"Real Time: "<<sw.RealTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"CPU Time: "<<sw.CpuTime()/60.0 <<" minutes"<<std::endl;
 std::cout<<"Done"<<std::endl;

} //end Loop()

