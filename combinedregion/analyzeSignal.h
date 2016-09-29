
#ifndef analyzeSignal_h
#define analyzeSignal_h

#include "postAnalyzer_Base.h"


class analyzeSignal : public postAnalyzer_Base {

 public:
///   analyzeSignal();
///   virtual ~analyzeSignal();
   virtual void     Loop(TString outfilename, Bool_t isMC,
                         Double_t lumi, Double_t nrEvents,
                         Double_t crossSec, Bool_t isZnnG,
                         Bool_t isEle, Bool_t isHalo,
                         Bool_t isSpike, Bool_t isJet,
                         Bool_t ewkWG, Bool_t ewkZG);
};

#endif
