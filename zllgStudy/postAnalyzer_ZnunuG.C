#define postAnalyzer_ZnunuG_cxx
#include "postAnalyzer_ZnunuG.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void postAnalyzer_ZnunuG::Loop(TString outfilename, Bool_t isMC,
  Double_t lumi, Double_t nrEvents, Double_t crossSec,
  Bool_t isZnnG, Bool_t ewkZG)
{

 TStopwatch sw; 
 sw.Start();

 Int_t nc = 0;
 Int_t dc = 0;
 gencount = 0;

 nrec_m110 = 0;
 nrec_m170 = 0;

 Int_t n_total          = 0;
 Int_t n_photons        = 0;
 Int_t n_passSpike      = 0;
 Int_t n_passMET110     = 0;
 Int_t n_passMET170     = 0;
 Int_t n_passMETfilters = 0;
 Int_t n_passTrig       = 0;
 Int_t n_passdPhiJM     = 0;
 Int_t n_passdPhiPhoMET = 0;
 Int_t n_passZWindow    = 0;
 Int_t n_passdimass20   = 0;

 int inclptbin = ptbins.size()-1;
 int lastptbin = ptbinnames.size()-1;
 int lastselbin = selbinnames.size();

  TFile *f_ewk_corr = new TFile("ewk_corr.root");
  TH1D *ewkZGCorrection = (TH1D*)f_ewk_corr->Get("zg");
  cout<<"ewkZG Correction histogram loaded"<<endl;

 if (fChain == 0) return;
 //std::cout<<"Numerator/Denominator : run : lumis : event : photonindex "<<std::endl;

 Long64_t nentries = fChain->GetEntriesFast();

 Long64_t nbytes = 0, nb = 0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
  n_total++;
 //for (Long64_t jentry=0; jentry<1000;jentry++) {
  if (jentry%20000 == 0)
    {
      std::cout<<"Starting entry "<<jentry<<"/"<<(nentries)<<" at "<<sw.RealTime()<<" RealTime, "<<sw.CpuTime()<<" CpuTime"<<std::endl;
      sw.Continue();
    }

  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   nbytes += nb;

  passSpike = false;
  passMET110 = false;
  passMET170 = false;
  passMETfilters = false;
  passTrig = false;
  passdPhiJM = false;
  passdPhiPhoMET = false;
  passZWindow = false;
  passdimass20 = false;

  passGenMET110 = false;
  passGenMET170 = false;
  passGendPhiPhoMET = false;
  passGenZWindow = false;

   // vector of ints, each int corresponds to location in vector of photons of photon passing cut
   // one such vector for each systematic bin

  //for(unsigned int selb=0; selb<lastselbin; ++selb){
  // phoCand[selb] = getPhoCand(175,1.4442);
  //}
  phoCand = getPhoCand(175,1.4442);
   if(phoCand.size()>0)
     {
      n_photons++;
      //printf("\n%llu Event Number", event);
      int candphotonindex = phoCand.at(0);

      double phopt = phoEt->at(candphotonindex) ; 
      Float_t uncorphopt = ((*phoSCRawE)[candphotonindex]/TMath::CosH((*phoSCEta)[candphotonindex]));
      // Event Weight
      //=1.0 for real data
      event_weight=1.0;
      if(isZnnG){
       if      ( uncorphopt < 190 ) {crossSec*=1.39;}  // {crossSec=14.4*1.39;}
       else if ( uncorphopt < 250 ) {crossSec*=1.35;}  // {crossSec=29.7*1.35;}
       else if ( uncorphopt < 400 ) {crossSec*=1.30;}  // {crossSec=17.5*1.30;}
       else if ( uncorphopt < 700 ) {crossSec*=1.23;}  // {crossSec=3.7*1.23;}
       else                         {crossSec*=1.23;}  // {crossSec=0.3*1.23;}
      } 
      if(isMC){ event_weight=lumi*crossSec/nrEvents; }
      //if(isMC){ event_weight=lumi*crossSec*(1.013 - 0.0001168*phopt)/nrEvents; }
      if(ewkZG){ 
       Double_t EWK_percent_adjustment = ewkZGCorrection->GetBinContent(ewkZGCorrection->GetXaxis()->FindBin(uncorphopt));
       event_weight*=(1.0+.01*EWK_percent_adjustment) ; 
      }

      fourVec_l1.SetPtEtaPhiE(0,0,0,0);
      fourVec_l2.SetPtEtaPhiE(0,0,0,0);

      // get electrons and muons and put into 4vectors
      makeDilep(candphotonindex, &fourVec_l1, &fourVec_l2, &fourVec_ee, &fourVec_mm);
      // make dilepton object and add to met
      fourVec_ll = fourVec_l1 + fourVec_l2;
      //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_l1.Pt(), fourVec_l1.Eta(), fourVec_l1.Phi(), fourVec_l1.M() );
      //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_l2.Pt(), fourVec_l2.Eta(), fourVec_l2.Phi(), fourVec_l2.M() );
      //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_ll.Pt(), fourVec_ll.Eta(), fourVec_ll.Phi(), fourVec_ll.M() );
      dilep_mass = fourVec_ll.M();
      //printf(" dilep mass: %0.1f\n", dilep_mass );

      // Dilepton Mass Window
      passZWindow = (dilep_mass>60. && dilep_mass<120.);
      passdimass20 = (dilep_mass>20.);

      // Lepto MET Creation
      fourVec_met.SetPtEtaPhiE(pfMET, 0., pfMETPhi, pfMET);
      fourVec_leptomet = fourVec_met + fourVec_ll;
      leptoMET = fourVec_leptomet.Pt();
      leptoMEPhi = fourVec_leptomet.Phi();

      // MET SELECTIONS
      passMET110 = (leptoMET > 110.) ;
      passMET170 = (leptoMET > 170.) ;
      passMETfilters = ( metFilters==0 ) ;
      if(isMC){ passMETfilters = true ; }

      // Spike Cleaning
      int iphi = 41; 
      int ieta = 5;
      passSpike = !(phoIPhi->at(candphotonindex) == iphi && phoIEta->at(candphotonindex) == ieta) ;
      if(isMC){ passSpike = true ; }
      //passSpike = true; //!(phoIPhi->at(candphotonindex) == iphi && phoIEta->at(candphotonindex) == ieta) ;

      // TRIGGER (HLT_Photon165_HE10_v)
      // https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc#L179
      passTrig = (
       ((HLTPho>>7&1) == 1) ||
       ((HLTPho>>8&1) == 1) ||
       ((HLTPho>>9&1) == 1) ||
       ((HLTPho>>10&1) == 1) ||
       ((HLTPho>>11&1) == 1) ||
       ((HLTPho>>12&1) == 1) ||
       ((HLTPho>>22&1) == 1)
      ) ;
      if(isMC){ passTrig = true ; }
      //passTrig = true; //( (HLTPho>>12&1) == 1);
 
      // dPhi( Jets, MET )
      std::vector<int>  jetindexvector = selectedJets(candphotonindex);
      passdPhiJM = passdphiJetMET(&jetindexvector, leptoMEPhi);

      // dPhi( photon, MET )
      passdPhiPhoMET = ( DeltaPhi(phoPhi->at(candphotonindex),leptoMEPhi)>2.0 ) ;

      // fill histograms
      bool baseline = (passSpike && passMETfilters && passTrig && passdimass20 && passdPhiPhoMET && passZWindow);

      if( passSpike      ) {  n_passSpike      ++ ; }
      if( passMET110     ) {  n_passMET110     ++ ; }
      if( passMET170     ) {  n_passMET170     ++ ; }
      if( passMETfilters ) {  n_passMETfilters ++ ; }
      if( passTrig       ) {  n_passTrig       ++ ; }
      if( passdPhiJM     ) {  n_passdPhiJM     ++ ; }
      if( passdPhiPhoMET ) {  n_passdPhiPhoMET ++ ; }
      if( passZWindow    ) {  n_passZWindow    ++ ; }
      if( passdimass20   ) {  n_passdimass20   ++ ; }

      if( baseline && passMET110 ){
       nrec_m110++; 
       callFillSigHist(0, lastptbin, inclptbin, candphotonindex, event_weight); 
       std::cout<<" nrec_m110<< run:lumis:event "
        <<run<<":"<<lumis<<":"<<event<< std::endl; 
       }
      if( baseline && passMET170 ){
       nrec_m170++; 
       callFillSigHist(1, lastptbin, inclptbin, candphotonindex, event_weight); 
       std::cout<<" nrec_m110<< run:lumis:event "
        <<run<<":"<<lumis<<":"<<event<< std::endl; 
       }
     } //end if phoCand[0].size()>0


//////////////////////////////////////////////////////////////////
// if(isMC){
//  // now do gen part
//    std::vector<int> genPhotonList;
//    std::vector<int> genMuonList;
//    std::vector<int> genEleList;
//
//    for( int a=0; a<nMC; ++a ){
//     if( mcPID->at(a)==22 && (((mcStatusFlag->at(a)>>0)&1)==1 || ((mcStatusFlag->at(a)>>1)&1)==1)
//      && mcPt->at(a)>175. && abs(mcEta->at(a))<1.4442){ genPhotonList.push_back(a); }
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
//   int gencandphotonindex = genPhotonList.at(0);
//
//   //printf(" photonindex: %i  pT: %.3f \n",candphotonindex,mcPt->at(candphotonindex));
//
//   fourVec_genl1.SetPtEtaPhiE(0,0,0,0);
//   fourVec_genl2.SetPtEtaPhiE(0,0,0,0);
//
//   // get electrons and muons and put into 4vectors
//   makeGenDilep(gencandphotonindex, genEleList, genMuonList, &fourVec_genl1, &fourVec_genl2, &fourVec_genee, &fourVec_genmm);
//   // make dilepton object and add to met
//   fourVec_genll = fourVec_genl1 + fourVec_genl2;
//   //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_l1.Pt(), fourVec_l1.Eta(), fourVec_l1.Phi(), fourVec_l1.M() );
//   //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_l2.Pt(), fourVec_l2.Eta(), fourVec_l2.Phi(), fourVec_l2.M() );
//   //printf("  pt: %.4f  eta: %.4f  phi: %.4f  mass: %.4f \n", fourVec_ll.Pt(), fourVec_ll.Eta(), fourVec_ll.Phi(), fourVec_ll.M() );
//   gen_dilep_mass = fourVec_genll.M();
//
//   // Z Mass Window
//   passGenZWindow =( gen_dilep_mass>60. && gen_dilep_mass<120. );
//
//   // Make Leptomet
//   fourVec_genMET.SetPtEtaPhiE( genMET, 0., genMETPhi, 0. );
//   fourVec_genLeptoMET = fourVec_genMET + fourVec_genll;
//   genLeptoMET = fourVec_genLeptoMET.Pt();
//   genLeptoMEPhi = fourVec_genLeptoMET.Phi();
//
//   // Pass MET
//   passGenMET110 = genLeptoMET > 110.;
//   passGenMET170 = genLeptoMET > 170.;
//
//   // dPhi( photon, MET ) 
//   bool passGendPhiPhoMET = ( DeltaPhi( mcPhi->at(gencandphotonindex), genLeptoMEPhi ) > 2.0 ) ;
//
//   ///// dPhi( Jets, MET )
//   /// std::vector<int>  genjetindexvector = selectedGenJets(gencandphotonindex);
//   ///  //for( int i=0; i<genjetindexvector.size(); ++i ){ printf(" %i \n",i);  }
//   /// bool passGendPhiJM = passGendphiJetMET(&genjetindexvector, genLeptoMEPhi);
//  
//  // fill histograms
//  //if ( passGenMET && passGendPhiPhoMET && gen_dilep_mass>60. && gen_dilep_mass<120.)
//  //if ( passGenMET && gen_dilep_mass>60. && gen_dilep_mass<120.)
//  //if ( passGenMET && passGendPhiPhoMET && gen_dilep_mass>60. && gen_dilep_mass<120.)
////  if ( passGenMET && passGendPhiPhoMET && passGendPhiJM && gen_dilep_mass>60. && gen_dilep_mass<120.)
////     {
////         printf(" %lli passed GEN \n", event);
////      gencount++;
////      for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
////       if(
////          (mcEt->at(gencandphotonindex) > ptbins[ptb]) &&
////          (mcEt->at(gencandphotonindex) < ptbins[ptb+1])
////         ){
////        FillGenHistograms(ptb, 0, gencandphotonindex, event_weight);
////       } // end if passes pt cuts then fill
////      } // end pt bin loop
////      if(  // do an inclusive pT plot from bins
////         (mcEt->at(gencandphotonindex) > ptbins[0]) &&
////         (mcEt->at(gencandphotonindex) < ptbins[inclptbin])
////        ){
////       FillGenHistograms(lastptbin-1, 0, gencandphotonindex, event_weight);
////      }
////      // and one fully inclusive in pT
////      FillGenHistograms(lastptbin, 0, gencandphotonindex, event_weight);
////     }
//  // end fill histograms
//
//      // fill histograms
//
////  if( gen_dilep_mass>20. && passGenMET110 && passGenZWindow && passGendPhiPhoMET ){ngen_m110_ywnd_ydphi++; callFillGenHist(0, lastptbin, inclptbin, gencandphotonindex, event_weight); }
////  if( gen_dilep_mass>20. && passGenMET110 && passGenZWindow                      ){ngen_m110_ywnd_ndphi++; callFillGenHist(1, lastptbin, inclptbin, gencandphotonindex, event_weight); }
////  if( gen_dilep_mass>20. && passGenMET110                   && passGendPhiPhoMET ){ngen_m110_nwnd_ydphi++; callFillGenHist(2, lastptbin, inclptbin, gencandphotonindex, event_weight); }
////  if( gen_dilep_mass>20. && passGenMET110                                        ){ngen_m110_nwnd_ndphi++; callFillGenHist(3, lastptbin, inclptbin, gencandphotonindex, event_weight); }
////  if( gen_dilep_mass>20. && passGenMET170 && passGenZWindow && passGendPhiPhoMET ){ngen_m170_ywnd_ydphi++; callFillGenHist(4, lastptbin, inclptbin, gencandphotonindex, event_weight); }
////  if( gen_dilep_mass>20. && passGenMET170 && passGenZWindow                      ){ngen_m170_ywnd_ndphi++; callFillGenHist(5, lastptbin, inclptbin, gencandphotonindex, event_weight); }
////  if( gen_dilep_mass>20. && passGenMET170                   && passGendPhiPhoMET ){ngen_m170_nwnd_ydphi++; callFillGenHist(6, lastptbin, inclptbin, gencandphotonindex, event_weight); }
////  if( gen_dilep_mass>20. && passGenMET170                                        ){ngen_m170_nwnd_ndphi++; callFillGenHist(7, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//
//
//  if( gen_dilep_mass>20. && passGenMET110 && passGenZWindow && passGendPhiPhoMET && passJets){ngen_m110_ywnd_ydphi++; callFillGenHist(0, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//  if( gen_dilep_mass>20. && passGenMET110 && passGenZWindow                      && passJets){ngen_m110_ywnd_ndphi++; callFillGenHist(1, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//  if( gen_dilep_mass>20. && passGenMET110                   && passGendPhiPhoMET && passJets){ngen_m110_nwnd_ydphi++; callFillGenHist(2, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//  if( gen_dilep_mass>20. && passGenMET110                                        && passJets){ngen_m110_nwnd_ndphi++; callFillGenHist(3, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//  if( gen_dilep_mass>20. && passGenMET170 && passGenZWindow && passGendPhiPhoMET && passJets){ngen_m170_ywnd_ydphi++; callFillGenHist(4, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//  if( gen_dilep_mass>20. && passGenMET170 && passGenZWindow                      && passJets){ngen_m170_ywnd_ndphi++; callFillGenHist(5, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//  if( gen_dilep_mass>20. && passGenMET170                   && passGendPhiPhoMET && passJets){ngen_m170_nwnd_ydphi++; callFillGenHist(6, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//  if( gen_dilep_mass>20. && passGenMET170                                        && passJets){ngen_m170_nwnd_ndphi++; callFillGenHist(7, lastptbin, inclptbin, gencandphotonindex, event_weight); }
//
//  } // end gen photon list nonempty
// }
/////////////////////////////////////////////////////////////////


///if( event==1531132 ){
///printf(" %lli \n", event);
/////printf ( "   passSpike           %i \n",      passSpike            );     
/////printf ( "   passMETfilters      %i \n",      passMETfilters       ); 
/////printf ( "   passMET             %i \n",      passMET              );
/////printf ( "   passdPhiPhoMET      %i \n",      passdPhiPhoMET       ); 
/////printf ( "   passdPhiJM          %i \n",      passdPhiJM           );
/////printf ( "   passTrig            %i \n",      passTrig             );
/////printf ( "   dilep_mass>60.      %i \n",      dilep_mass>60.       );  
/////printf ( "   dilep_mass<120.     %i \n",      dilep_mass<120.      );   
/////
/////printf ( "   passGenMET          %i \n",      passGenMET           );
/////printf ( "   passGendPhiPhoMET   %i \n",      passGendPhiPhoMET    );
/////printf ( "   passGendPhiJM       %i \n",      passGendPhiJM        );
/////printf ( "   gen_dilep_mass>60.  %i \n",      gen_dilep_mass>60.   );
/////printf ( "   gen_dilep_mass<120. %i \n",      gen_dilep_mass<120.  );
///
/// printf("N GenJets  = %zu \n",jetGenJetPt->size());
/// printf("N RecoJets = %i \n", nJet);
/// for(int i = 0; i < jetGenJetPt->size(); i++)
/// {
/////  printf(" jetPartonID     %i \n" ,  jetPartonID   ->at(i) ) ; 
/////  printf(" jetGenJetIndex  %i \n" ,  jetGenJetIndex->at(i) ) ; 
/////  printf(" jetGenJetEn     %f \n" ,  jetGenJetEn   ->at(i) ) ; 
///  printf(" jetGenJetPt     %f \n" ,  jetGenJetPt   ->at(i) ) ; 
///  printf(" jetGenJetEta    %f \n" ,  jetGenJetEta  ->at(i) ) ;  
///  printf(" jetGenJetPhi    %f \n\n" ,  jetGenJetPhi  ->at(i) ) ;  
/// }
/// for(int i = 0; i < nJet; i++)
/// {
///  printf("  jetPt                     %f \n", jetPt                  ->at(i) );  
///  printf("  jetEta                    %f \n", jetEta                 ->at(i) );   
///  printf("  jetPhi                    %f \n", jetPhi                 ->at(i) );  
///  printf("  jetPFLooseId              %f \n", jetPFLooseId           ->at(i) );   
///  printf("  jetPUidFullDiscriminant   %f \n\n", jetPUidFullDiscriminant->at(i) );
/// }
///}
         //printf(" %lli passed RECO \n", event);
         //if( jetGenJetPt->size() != nJet){
         // printf("N GenJets  = %zu \n",jetGenJetPt->size());
         // printf("N RecoJets = %i \n", nJet);
         // for(int i = 0; i < jetGenJetPt->size(); i++)
         // {
         //  printf(" jetPartonID     %i \n" ,  jetPartonID   ->at(i) ) ; 
         //  printf(" jetGenJetIndex  %i \n" ,  jetGenJetIndex->at(i) ) ; 
         //  printf(" jetGenJetEn     %f \n" ,  jetGenJetEn   ->at(i) ) ; 
         //  printf(" jetGenJetPt     %f \n" ,  jetGenJetPt   ->at(i) ) ; 
         //  printf(" jetGenJetEta    %f \n" ,  jetGenJetEta  ->at(i) ) ;  
         //  printf(" jetGenJetPhi    %f \n" ,  jetGenJetPhi  ->at(i) ) ;  
         // }
         // for(int i = 0; i < nJet; i++)
         // {
         //  printf("  jetEta                    %f \n", jetEta                 ->at(i) );   
         //  printf("  jetPt                     %f \n", jetPt                  ->at(i) );  
         //  printf("  jetPhi                    %f \n", jetPhi                 ->at(i) );  
         //  printf("  jetPFLooseId              %f \n", jetPFLooseId           ->at(i) );   
         //  printf("  jetPUidFullDiscriminant   %f \n", jetPUidFullDiscriminant->at(i) );
         // }
         //}



 } //end for jentry -  loop through entries

 // write these histograms to file

 printf("nrec_m110 : %i \n", nrec_m110) ;  
 printf("nrec_m170 : %i \n", nrec_m170) ;  

 printf("n_total          : %i \n", n_total           ) ;
 printf("n_photons        : %i \n", n_photons         ) ;
 printf("n_passSpike      : %i \n", n_passSpike       ) ;
 printf("n_passMET110     : %i \n", n_passMET110      ) ;
 printf("n_passMET170     : %i \n", n_passMET170      ) ;
 printf("n_passMETfilters : %i \n", n_passMETfilters  ) ;
 printf("n_passTrig       : %i \n", n_passTrig        ) ;
 printf("n_passdPhiJM     : %i \n", n_passdPhiJM      ) ;
 printf("n_passdPhiPhoMET : %i \n", n_passdPhiPhoMET  ) ;
 printf("n_passZWindow    : %i \n", n_passZWindow     ) ;
 printf("n_passdimass20   : %i \n", n_passdimass20    ) ;

 std::cout<<"made it through, about to write"<<std::endl;

 TFile *outfile = new TFile(outfilename,"RECREATE");
 outfile->cd();
 for(unsigned int i=0; i<ptbinnames.size(); ++i){
  for(unsigned int j=0; j<selbinnames.size(); ++j){
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
std::vector<int> postAnalyzer_ZnunuG::getPhoCand(double phoPtCut, double phoEtaCut){

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
                       ( TMath::Max( ( (*phoPFChIso)[p]       - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.37 )  &&
                       ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 1.37 )  &&
                       ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) <
                        (1.06 + (0.014 * (*phoEt)[p]) + (0.000019 * pow((*phoEt)[p], 2.0))) )  &&
                       ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < 
                        (0.28 + (0.0053 * (*phoEt)[p])) ) 
                      );

                     
      bool noncoll = fabs((*phoseedTimeFull5x5)[p]) < 3. && 
                     (*phomipTotEnergy)[p] < 4.9 && 
                     (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && 
                     (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001;

      if(photonId && kinematic && noncoll){
        pholist.push_back(p);
      } 
    }  

  return pholist;

}


//-------------------------selectedJets
std::vector<int> postAnalyzer_ZnunuG::selectedJets(int pho_index) {

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


////-------------------------selectedGenJets
//std::vector<int> postAnalyzer_ZnunuG::selectedGenJets(int pho_index) {
//
//
//  bool jetVeto=true;
//  std::vector<int> jetindex;
//  float value = 0.0;
//
//  for(int i = 0; i < jetGenJetPt->size(); i++)
//    {
//
//      double deltar = 0.0 ;
//      //      std::cout<<"Jet size: "<<nJet<<std::endl;
//      //      std::cout<<"Jet no:"<<i<<"coming here pujetid: "<<pfJet_pt[i]<<std::endl;
//      //      if(OverlapWithMuon(jetEta->at(i),jetPhi->at(i)))     continue;
//      //      std::cout<<"Jet no:"<<i<<"coming here OverlapWithMuon: "<<pfJet_pt[i]<<std::endl;
//      //      if(OverlapWithElectron(jetEta->at(i),jetPhi->at(i)))   continue;
//      if(pho_index>=0){
//        deltar= dR(jetGenJetEta->at(i),jetGenJetPhi->at(i),mcEta->at(pho_index),mcPhi->at(pho_index));
//        //std::cout<<"Delta R between photon and jet="<<dR<< "jet pt"<<pfJet_pt[i]<<std::endl;
//      }
//      if(deltar>0.4 && jetGenJetPt->at(i) >30.0)
//        {
//          jetindex.push_back(i);
//        }
//    }
//
//  //  std::cout<<"Jet size: "<< jetindex.size()<<std::endl;
//  //if(jetindex.size()>1)jetVeto = false;
//  return jetindex;
//
//}


//-------------------------dR
double postAnalyzer_ZnunuG::dR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

//-------------------------DeltaPhi
double postAnalyzer_ZnunuG::DeltaPhi(double phi1, double phi2)
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
bool postAnalyzer_ZnunuG::passdphiJetMET(std::vector<int> *jets, double mephi)
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


// //-------------------------passGendphiJetMET
// bool postAnalyzer_ZnunuG::passGendphiJetMET(std::vector<int> *jets, double mephi)
// {
// 
//   //pass veto if no jet is found within DeltaPhi(jet,MET) < 0.5                                                                                         
//   bool passes = false;
//   int njetsMax = jets->size();
//   //Only look at first four jets                                                                                                                        
//   if(njetsMax > 4)
//     njetsMax = 4;  
//   int j = 0;
//   for(; j < njetsMax; j++)
//     {   
//       //fail veto if a jet is found with DeltaPhi(jet,MET) < 0.5                                                                                    
//       if(DeltaPhi(jetGenJetPhi->at(jets->at(j)), mephi) < 0.5)
//         break;
//     }    
//   if(j==njetsMax)
//     passes = true;
// 
//   return passes;
// }
// 
// 
// 
// //-------------------------makeGenDilep
// void postAnalyzer_ZnunuG::makeGenDilep(int pho_index, std::vector<int> elelist, std::vector<int> mulist, TLorentzVector *fv_1, TLorentzVector *fv_2, TLorentzVector *fv_ee, TLorentzVector *fv_mm)
// {
// 
//   TLorentzVector e1, e2, ee;
//   TLorentzVector m1, m2, mm;
//   e1.SetPtEtaPhiE( 0,0,0,0 );
//   e2.SetPtEtaPhiE( 0,0,0,0 );
//   m1.SetPtEtaPhiE( 0,0,0,0 );
//   m2.SetPtEtaPhiE( 0,0,0,0 );
// 
//   int    leadingE, leadingM;
//   leadingE = 0;
//   leadingM = 0;
//   double leadingEpt, leadingMpt;
//   leadingEpt = 0.;
//   leadingMpt = 0.;
// 
//   for( int i=0; i<mulist.size(); ++i ){ 
//    if (mcEt->at(mulist.at(i)) > leadingMpt){
//     leadingMpt = mcEt->at(mulist.at(i)) ;
//     leadingM = mulist.at(i) ;
//    } 
//   }
// 
//   for( int i=0; i<elelist.size(); ++i ){ 
//    if (mcEt->at(elelist.at(i)) > leadingEpt){
//     leadingEpt = mcEt->at(elelist.at(i)) ;
//     leadingE = elelist.at(i) ;
//    } 
//   }
// 
// 
//   // no pairs
//   if( elelist.size()<2 && mulist.size()<2 ){return;}
// 
//   // electrons
//   if( elelist.size()>1 ){
//    for(int i=0; i<elelist.size(); ++i)
//    {
//      if( mcPID->at(leadingE) + mcPID->at(elelist[i])==0 )
//      {
//       //printf(" --we have electrons ");
//       e1.SetPtEtaPhiE( mcPt->at(leadingE), mcEta->at(leadingE), mcPhi->at(leadingE), mcE->at(leadingE) );
//       e2.SetPtEtaPhiE( mcPt->at(elelist[i]), mcEta->at(elelist[i]), mcPhi->at(elelist[i]), mcE->at(elelist[i]) );
//       break;
//      }
//    }
//    ee = e1 + e2;
//    //printf(": dilep mass = %f", ee.M());
//   }
// 
//   // muons
//   if( mulist.size()>1 ){
//    for(int i=0; i<mulist.size(); ++i)
//    {
//      if( mcPID->at(leadingM) + mcPID->at(mulist[i])==0 )
//      {
//       //printf(" --we have muons ");
//       m1.SetPtEtaPhiE( mcPt->at(leadingM), mcEta->at(leadingM), mcPhi->at(leadingM), mcE->at(leadingM) );
//       m2.SetPtEtaPhiE( mcPt->at(mulist[i]), mcEta->at(mulist[i]), mcPhi->at(mulist[i]), mcE->at(mulist[i]) );
//       break;
//      }
//    }
//    mm = m1 + m2;
//    //printf(": dilep mass = %f", mm.M());
//   }
// 
//   *fv_ee = ee;
//   *fv_mm = mm;
//   // take highest mass dilepton pair
//   if( mm.M()>ee.M() ){ *fv_1 = m1; *fv_2 = m2; }
//   else               { *fv_1 = e1; *fv_2 = e2; }
//   ///if( mm.M()>60. && mm.M()<120. ){ *fv_1 = m1; *fv_2 = m2; }
//   ///else if( ee.M()>60. && ee.M()<120. ){ *fv_1 = e1; *fv_2 = e2; }
//   //if( ee.M() > mm.M() && ee.M()>60. && ee.M()<120. ){ printf(": using ee\n"); *fv_1 = e1; *fv_2 = e2; }
//   //else if( ee.M() < mm.M() && mm.M()>60. && mm.M()<120. ){ printf(": using mm\n"); *fv_1 = m1; *fv_2 = m2; }
// 
//   return;
//   
// }


//-------------------------makeDilep
void postAnalyzer_ZnunuG::makeDilep(int pho_index, TLorentzVector *fv_1, TLorentzVector *fv_2, TLorentzVector *fv_ee, TLorentzVector *fv_mm)
{
  
  std::vector<int> elelist = electron_passLooseID(pho_index, 10.);
  std::vector<int> mulist = muon_passLooseID(pho_index, 10.);

  TLorentzVector e1, e2, ee;
  TLorentzVector m1, m2, mm;
  e1.SetPtEtaPhiE( 0,0,0,0 );
  e2.SetPtEtaPhiE( 0,0,0,0 );
  m1.SetPtEtaPhiE( 0,0,0,0 );
  m2.SetPtEtaPhiE( 0,0,0,0 );

  // // no pairs
  // if( elelist.size()<2 && mulist.size()<2 ){return;}

  // electrons
  if( elelist.size()>1 ){
   for(int i=1; i<elelist.size(); ++i)
   {
    if( eleCharge->at(0)*eleCharge->at(i)==-1 )
    {
     //printf(" --we have electrons ");
     e1.SetPtEtaPhiE( elePt->at(elelist[0]), eleEta->at(elelist[0]), elePhi->at(elelist[0]), eleEn->at(elelist[0]) );
     e2.SetPtEtaPhiE( elePt->at(elelist[i]), eleEta->at(elelist[i]), elePhi->at(elelist[i]), eleEn->at(elelist[i]) );
     break;
    }
   }
  }
  //printf(": diele mass = %f", ee.M());
  ee = e1 + e2;

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
  }
  //printf(": dimu mass = %f", mm.M());
  mm = m1 + m2;

  *fv_ee = ee;
  *fv_mm = mm;
  // take highest mass dilepton pair
  if( mm.M()>ee.M() ){ *fv_1 = m1; *fv_2 = m2; }
  else               { *fv_1 = e1; *fv_2 = e2; }
  return;
  
}


//-------------------------electron_passLooseID
std::vector<int> postAnalyzer_ZnunuG::electron_passLooseID(int pho_index, float elePtCut)
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
std::vector<int> postAnalyzer_ZnunuG::muon_passLooseID(int pho_index, float muPtCut)
{
  std::vector<int> mulist;

  bool pass_PFMuon      = false;
  bool pass_globalMuon  = false;
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


////-------------------------OverlapWithMuon
//bool postAnalyzer_ZnunuG::OverlapWithMuon(double eta, double phi){
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
//bool postAnalyzer_ZnunuG::OverlapWithElectron(double eta, double phi){
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
Double_t postAnalyzer_ZnunuG::EAcharged(Double_t eta){
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

Double_t postAnalyzer_ZnunuG::EAchargedworst(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.080721;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.063986;

  return EffectiveArea;
}

Double_t postAnalyzer_ZnunuG::EAneutral(Double_t eta){
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

Double_t postAnalyzer_ZnunuG::EAphoton(Double_t eta){
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

//-------------------------callFillSigHist
void postAnalyzer_ZnunuG::callFillSigHist(int selbin, int lastptbin, int inclptbin, int candphotonindex, float event_weight){
 for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
  if(
     (phoEt->at(candphotonindex) > ptbins[ptb]) &&
     (phoEt->at(candphotonindex) < ptbins[ptb+1])
    ){
   FillSigHistograms(ptb, selbin, candphotonindex, event_weight);
  } // end if passes pt cuts then fill
 } // end pt bin loop
 if(  // do an inclusive pT plot from bins
    (phoEt->at(candphotonindex) > ptbins[0]) &&
    (phoEt->at(candphotonindex) < ptbins[inclptbin])
   ){
  FillSigHistograms(lastptbin-1, selbin, candphotonindex, event_weight);
 }
 // and one fully inclusive in pT
 FillSigHistograms(lastptbin, selbin, candphotonindex, event_weight);
 return;
}

// //-------------------------callFillGenHist
// void postAnalyzer_ZnunuG::callFillGenHist(int selbin, int lastptbin, int inclptbin, int gencandphotonindex, float event_weight){
//  for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
//   if(
//      (mcEt->at(gencandphotonindex) > ptbins[ptb]) &&
//      (mcEt->at(gencandphotonindex) < ptbins[ptb+1])
//     ){
//    FillGenHistograms(ptb, selbin, gencandphotonindex, event_weight);
//   } // end if passes pt cuts then fill
//  } // end pt bin loop
//  if(  // do an inclusive pT plot from bins
//     (mcEt->at(gencandphotonindex) > ptbins[0]) &&
//     (mcEt->at(gencandphotonindex) < ptbins[inclptbin])
//    ){
//   FillGenHistograms(lastptbin-1, selbin, gencandphotonindex, event_weight);
//  }
//  // and one fully inclusive in pT
//  FillGenHistograms(lastptbin, selbin, gencandphotonindex, event_weight);
//  return;
// }
