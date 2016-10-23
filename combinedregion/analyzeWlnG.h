
#ifndef analyzeWlnG_h
#define analyzeWlnG_h

#include "postAnalyzer_Lep.h"

class analyzeWlnG : public postAnalyzer_Lep {

 public:

   std::vector<int> tightEles;
   std::vector<int> looseEles;
   std::vector<int> tightMus ;
   std::vector<int> looseMus ;

   TLorentzVector fourVec_e, fourVec_m; 

   virtual void     Loop(TString outfilename, Bool_t isMC,
                         Double_t lumi, Double_t nrEvents,
                         Double_t crossSec, Bool_t isZnnG,
                         Bool_t isEle, Bool_t isHalo,
                         Bool_t isSpike, Bool_t isJet,
                         Bool_t ewkWG, Bool_t ewkZG); 

};


#endif
