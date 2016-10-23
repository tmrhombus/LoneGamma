
#ifndef analyzePDFscaleZllG_h
#define analyzePDFscaleZllG_h

#include "postAnalyzer_PDFscale.h"

class analyzePDFscaleZllG : public postAnalyzer_PDFscale {

 public:
   virtual void     Loop(TString outfilename, Bool_t isMC,
                         Double_t lumi, Double_t nrEvents,
                         Double_t crossSec, Bool_t isZnnG,
                         Bool_t isEle, Bool_t isHalo,
                         Bool_t isSpike, Bool_t isJet,
                         Bool_t ewkWG, Bool_t ewkZG); 

   virtual void     makeDilep(int pho_index, TLorentzVector *fv_1, TLorentzVector *fv_2,
                                             TLorentzVector *fv_ee, TLorentzVector *fv_mm, bool *passMM);

};

//-------------------------makeDilep
void analyzePDFscaleZllG::makeDilep(int pho_index, TLorentzVector *fv_1, TLorentzVector *fv_2, TLorentzVector *fv_ee, TLorentzVector *fv_mm, bool *passMM)
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
  if( mm.M()>ee.M() ){ *fv_1 = m1; *fv_2 = m2; *passMM = true; }
  else               { *fv_1 = e1; *fv_2 = e2; *passMM = false; }
  return;
  
}

#endif
