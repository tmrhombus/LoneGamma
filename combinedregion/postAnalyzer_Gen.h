#ifndef postAnalyzer_Gen_h
#define postAnalyzer_Gen_h

#include "postAnalyzer_Lep.h"

class postAnalyzer_Gen : public postAnalyzer_Lep {

 public:
   TH1F h_gen_et[7],
        h_gen_etres[7],
        h_gen_met[7],
        h_gen_metres[7];

   virtual void     InitGen();
   Bool_t           FillSigHistogramsGen(int ptbin, int sysbin, int photonIndex, int genIndex, double weight);
   virtual void     callFillSigHistGen(int selbin, int lastptbin, int inclptbin, int candphotonindex, int genIndex, float event_weight);
   Bool_t           WriteHistogramsGen(int ptbin, int sysbin);

   vector<float>   *pdf;                             
   Float_t         pthat;                            
   Float_t         processID;                        
   Float_t         genWeight;                        
   Float_t         genHT;                            
   Float_t         pdfWeight;                        
   vector<float>   *pdfSystWeight;                   
   TString         *EventTag;                        
   Int_t           nPUInfo;                          
   vector<int>     *nPU;                             
   vector<int>     *puBX;                            
   vector<float>   *puTrue;                          
   Int_t           nlheWeights;                      
   vector<float>   *lheNormalizedWeights;            
   Float_t         genWeight_QCDscale_muR1_muF1;     
   Float_t         genWeight_QCDscale_muR1_muF2;     
   Float_t         genWeight_QCDscale_muR1_muF0p5;   
   Float_t         genWeight_QCDscale_muR2_muF1;     
   Float_t         genWeight_QCDscale_muR2_muF2;     
   Float_t         genWeight_QCDscale_muR2_muF0p5;   
   Float_t         genWeight_QCDscale_muR0p5_muF1;   
   Float_t         genWeight_QCDscale_muR0p5_muF2;   
   Float_t         genWeight_QCDscale_muR0p5_muF0p5; 
   vector<string>  *lheWeightIDs;                    
   Int_t           nMC;                              
   vector<int>     *mcPID;                           
   vector<float>   *mcVtx;                           
   vector<float>   *mcVty;                           
   vector<float>   *mcVtz;                           
   vector<float>   *mcPt;                            
   vector<float>   *mcMass;                          
   vector<float>   *mcEta;                           
   vector<float>   *mcPhi;                           
   vector<float>   *mcE;                             
   vector<float>   *mcEt;                            
   vector<int>     *mcGMomPID;                       
   vector<int>     *mcMomPID;                        
   vector<float>   *mcMomPt;                         
   vector<float>   *mcMomMass;                       
   vector<float>   *mcMomEta;                        
   vector<float>   *mcMomPhi;                        
   vector<int>     *mcIndex;                         
   vector<unsigned short> *mcStatusFlag;             
   vector<int>     *mcParentage;                     
   vector<int>     *mcStatus;                        
   vector<float>   *mcCalIsoDR03;                    
   vector<float>   *mcTrkIsoDR03;  
   vector<float>   *mcCalIsoDR04;                     
   vector<float>   *mcTrkIsoDR04;    
   Float_t         genMET;                               
   Float_t         genMETPhi;         
   Float_t         pfMETuncorrected;  
   Float_t         pfMETPhiuncorrected;      
   vector<int>     *jetPartonID;                               
   vector<int>     *jetHadFlvr;                                
   vector<int>     *jetGenJetIndex;                            
   vector<float>   *jetGenJetEn;                               
   vector<float>   *jetGenJetPt;                               
   vector<float>   *jetGenJetEta;                              
   vector<float>   *jetGenJetPhi;                              
   vector<int>     *jetGenPartonID;                            
   vector<float>   *jetGenEn;                                  
   vector<float>   *jetGenPt;                                  
   vector<float>   *jetGenEta;                                 
   vector<float>   *jetGenPhi;                                 
   vector<int>     *jetGenPartonMomID;
   vector<int>     *AK8JetPartonID;                            
   vector<int>     *AK8JetHadFlvr;                             
   vector<int>     *AK8JetGenJetIndex;                         
   vector<float>   *AK8JetGenJetEn;                            
   vector<float>   *AK8JetGenJetPt;                            
   vector<float>   *AK8JetGenJetEta;                           
   vector<float>   *AK8JetGenJetPhi;                           
   vector<int>     *AK8JetGenPartonID;                         
   vector<float>   *AK8JetGenEn;                               
   vector<float>   *AK8JetGenPt;                               
   vector<float>   *AK8JetGenEta;                              
   vector<float>   *AK8JetGenPhi;                              
   vector<int>     *AK8JetGenPartonMomID;   
   TBranch        *b_pdf;   //!                                
   TBranch        *b_pthat;   //!                              
   TBranch        *b_processID;   //!                          
   TBranch        *b_genWeight;   //!                          
   TBranch        *b_genHT;   //!                              
   TBranch        *b_pdfWeight;   //!                          
   TBranch        *b_pdfSystWeight;   //!                      
   TBranch        *b_EventTag;   //!                           
   TBranch        *b_nPUInfo;   //!                            
   TBranch        *b_nPU;   //!                                
   TBranch        *b_puBX;   //!                               
   TBranch        *b_puTrue;   //!                             
   TBranch        *b_nlheWeights;   //!                        
   TBranch        *b_lheNormalizedWeights;   //!               
   TBranch        *b_genWeight_QCDscale_muR1_muF1;   //!       
   TBranch        *b_genWeight_QCDscale_muR1_muF2;   //!       
   TBranch        *b_genWeight_QCDscale_muR1_muF0p5;   //!     
   TBranch        *b_genWeight_QCDscale_muR2_muF1;   //!       
   TBranch        *b_genWeight_QCDscale_muR2_muF2;   //!       
   TBranch        *b_genWeight_QCDscale_muR2_muF0p5;   //!     
   TBranch        *b_genWeight_QCDscale_muR0p5_muF1;   //!     
   TBranch        *b_genWeight_QCDscale_muR0p5_muF2;   //!     
   TBranch        *b_genWeight_QCDscale_muR0p5_muF0p5;   //!   
   TBranch        *b_lheWeightIDs;   //!           
   TBranch        *b_nMC;   //!                    
   TBranch        *b_mcPID;   //!                  
   TBranch        *b_mcVtx;   //!                  
   TBranch        *b_mcVty;   //!                  
   TBranch        *b_mcVtz;   //!                  
   TBranch        *b_mcPt;   //!                   
   TBranch        *b_mcMass;   //!                 
   TBranch        *b_mcEta;   //!                  
   TBranch        *b_mcPhi;   //!                  
   TBranch        *b_mcE;   //!                    
   TBranch        *b_mcEt;   //!                   
   TBranch        *b_mcGMomPID;   //!              
   TBranch        *b_mcMomPID;   //!               
   TBranch        *b_mcMomPt;   //!                
   TBranch        *b_mcMomMass;   //!              
   TBranch        *b_mcMomEta;   //!               
   TBranch        *b_mcMomPhi;   //!               
   TBranch        *b_mcIndex;   //!                
   TBranch        *b_mcStatusFlag;   //!           
   TBranch        *b_mcParentage;   //!            
   TBranch        *b_mcStatus;   //!               
   TBranch        *b_mcCalIsoDR03;   //!           
   TBranch        *b_mcTrkIsoDR03;   //!           
   TBranch        *b_mcCalIsoDR04;   //!           
   TBranch        *b_mcTrkIsoDR04;   //! 
   TBranch        *b_genMET;   //!                 
   TBranch        *b_genMETPhi;   //! 
   TBranch        *b_pfMETuncorrected;   //!  
   TBranch        *b_pfMETPhiuncorrected;   //!  
   TBranch        *b_jetPartonID;   //!                
   TBranch        *b_jetHadFlvr;   //!                 
   TBranch        *b_jetGenJetIndex;   //!             
   TBranch        *b_jetGenJetEn;   //!                
   TBranch        *b_jetGenJetPt;   //!                
   TBranch        *b_jetGenJetEta;   //!               
   TBranch        *b_jetGenJetPhi;   //!               
   TBranch        *b_jetGenPartonID;   //!             
   TBranch        *b_jetGenEn;   //!                   
   TBranch        *b_jetGenPt;   //!                   
   TBranch        *b_jetGenEta;   //!                  
   TBranch        *b_jetGenPhi;   //!                  
   TBranch        *b_jetGenPartonMomID;   //!   
   TBranch        *b_AK8JetPartonID;   //!             
   TBranch        *b_AK8JetHadFlvr;   //!              
   TBranch        *b_AK8JetGenJetIndex;   //!          
   TBranch        *b_AK8JetGenJetEn;   //!             
   TBranch        *b_AK8JetGenJetPt;   //!             
   TBranch        *b_AK8JetGenJetEta;   //!            
   TBranch        *b_AK8JetGenJetPhi;   //!            
   TBranch        *b_AK8JetGenPartonID;   //!          
   TBranch        *b_AK8JetGenEn;   //!                
   TBranch        *b_AK8JetGenPt;   //!                
   TBranch        *b_AK8JetGenEta;   //!               
   TBranch        *b_AK8JetGenPhi;   //!               
   TBranch        *b_AK8JetGenPartonMomID;   //!       

};

void postAnalyzer_Gen::InitGen()
{

   std::cout<<"Initializing Gen"<<std::endl;

   for(unsigned int i=0; i<ptbinnames.size(); ++i){
//     // set up names
     TString histname_gen_et  = "h_gen_et_"+ptbinnames[i];
     TString histname_gen_etres  = "h_gen_etres_"+ptbinnames[i];
     TString histname_gen_met  = "h_gen_met_"+ptbinnames[i];
     TString histname_gen_metres  = "h_gen_metres_"+ptbinnames[i];

     Float_t binspt[] = { 175., 195., 250., 400., 700., 1000. };
     Int_t  binnumpt = sizeof(binspt)/sizeof(Float_t) - 1; // or just = 9

     // reserve histograms
     h_gen_et[i].Clear();
     h_gen_et[i] = TH1F(histname_gen_et,"Gen Photon Transverse Energy",binnumpt,binspt);
     //h_gen_et[i] = TH1F(histname_gen_et,"Photon Transverse Energy",165,175.,1000.);
     h_gen_et[i].Sumw2();
     //
     h_gen_etres[i].Clear();
     h_gen_etres[i] = TH1F(histname_gen_etres,"(p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}",201,-0.4,0.4);
     h_gen_etres[i].Sumw2();
     //
     h_gen_met[i].Clear();
     h_gen_met[i] = TH1F(histname_gen_met,"Gen MET",binnumpt,binspt);
     //h_gen_met[i] = TH1F(histname_gen_met,"Photon Transverse Energy",165,175.,1000.);
     h_gen_met[i].Sumw2();
     //
     h_gen_metres[i].Clear();
     h_gen_metres[i] = TH1F(histname_gen_metres,"(MET^{reco}-MET^{gen})/MET^{gen}",200,-1.,15.);
     h_gen_metres[i].Sumw2();
     //

   }  // ptbinnames

     // Set object pointer           
     pdf = 0;                        
     pdfSystWeight = 0;              
     EventTag = 0;                   
     nPU = 0;                        
     puBX = 0;                       
     puTrue = 0;                     
     lheNormalizedWeights = 0;       
     lheWeightIDs = 0;               
     mcPID = 0;                      
     mcVtx = 0;                      
     mcVty = 0;                      
     mcVtz = 0;                      
     mcPt = 0;                       
     mcMass = 0;                     
     mcEta = 0;                      
     mcPhi = 0;                      
     mcE = 0;                        
     mcEt = 0;                       
     mcGMomPID = 0;                  
     mcMomPID = 0;                   
     mcMomPt = 0;                    
     mcMomMass = 0;                  
     mcMomEta = 0;                   
     mcMomPhi = 0;                   
     mcIndex = 0;                    
     mcStatusFlag = 0;               
     mcParentage = 0;                
     mcStatus = 0;                   
     mcCalIsoDR03 = 0;               
     mcTrkIsoDR03 = 0;               
     mcCalIsoDR04 = 0;               
     mcTrkIsoDR04 = 0;               
     jetPartonID = 0;             
     jetHadFlvr = 0;              
     jetGenJetIndex = 0;          
     jetGenJetEn = 0;             
     jetGenJetPt = 0;             
     jetGenJetEta = 0;            
     jetGenJetPhi = 0;            
     jetGenPartonID = 0;          
     jetGenEn = 0;                
     jetGenPt = 0;                
     jetGenEta = 0;               
     jetGenPhi = 0;               
     jetGenPartonMomID = 0;
     AK8JetPartonID = 0;          
     AK8JetHadFlvr = 0;           
     AK8JetGenJetIndex = 0;       
     AK8JetGenJetEn = 0;          
     AK8JetGenJetPt = 0;          
     AK8JetGenJetEta = 0;         
     AK8JetGenJetPhi = 0;         
     AK8JetGenPartonID = 0;       
     AK8JetGenEn = 0;             
     AK8JetGenPt = 0;             
     AK8JetGenEta = 0;            
     AK8JetGenPhi = 0;            
     AK8JetGenPartonMomID = 0;    

     fChain->SetBranchAddress("pdf", &pdf, &b_pdf);                               
     fChain->SetBranchAddress("pthat", &pthat, &b_pthat);                         
     fChain->SetBranchAddress("processID", &processID, &b_processID);             
     fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);             
     fChain->SetBranchAddress("genHT", &genHT, &b_genHT);                         
     fChain->SetBranchAddress("pdfWeight", &pdfWeight, &b_pdfWeight);             
     fChain->SetBranchAddress("pdfSystWeight", &pdfSystWeight, &b_pdfSystWeight); 
     fChain->SetBranchAddress("EventTag", &EventTag, &b_EventTag);                
     fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);                   
     fChain->SetBranchAddress("nPU", &nPU, &b_nPU);                               
     fChain->SetBranchAddress("puBX", &puBX, &b_puBX);                            
     fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);                      
     fChain->SetBranchAddress("nlheWeights", &nlheWeights, &b_nlheWeights);       
     fChain->SetBranchAddress("lheNormalizedWeights", &lheNormalizedWeights, &b_lheNormalizedWeights);                               
     fChain->SetBranchAddress("genWeight_QCDscale_muR1_muF1", &genWeight_QCDscale_muR1_muF1, &b_genWeight_QCDscale_muR1_muF1);       
     fChain->SetBranchAddress("genWeight_QCDscale_muR1_muF2", &genWeight_QCDscale_muR1_muF2, &b_genWeight_QCDscale_muR1_muF2);       
     fChain->SetBranchAddress("genWeight_QCDscale_muR1_muF0p5", &genWeight_QCDscale_muR1_muF0p5, &b_genWeight_QCDscale_muR1_muF0p5); 
     fChain->SetBranchAddress("genWeight_QCDscale_muR2_muF1", &genWeight_QCDscale_muR2_muF1, &b_genWeight_QCDscale_muR2_muF1);       
     fChain->SetBranchAddress("genWeight_QCDscale_muR2_muF2", &genWeight_QCDscale_muR2_muF2, &b_genWeight_QCDscale_muR2_muF2);       
     fChain->SetBranchAddress("genWeight_QCDscale_muR2_muF0p5", &genWeight_QCDscale_muR2_muF0p5, &b_genWeight_QCDscale_muR2_muF0p5); 
     fChain->SetBranchAddress("genWeight_QCDscale_muR0p5_muF1", &genWeight_QCDscale_muR0p5_muF1, &b_genWeight_QCDscale_muR0p5_muF1); 
     fChain->SetBranchAddress("genWeight_QCDscale_muR0p5_muF2", &genWeight_QCDscale_muR0p5_muF2, &b_genWeight_QCDscale_muR0p5_muF2); 
     fChain->SetBranchAddress("genWeight_QCDscale_muR0p5_muF0p5", &genWeight_QCDscale_muR0p5_muF0p5, &b_genWeight_QCDscale_muR0p5_muF0p5);
     fChain->SetBranchAddress("lheWeightIDs", &lheWeightIDs, &b_lheWeightIDs); 
     fChain->SetBranchAddress("nMC", &nMC, &b_nMC);                            
     fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);                      
     fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);                      
     fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);                      
     fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);                      
     fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);                         
     fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);                   
     fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);                      
     fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);                      
     fChain->SetBranchAddress("mcE", &mcE, &b_mcE);                            
     fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);                         
     fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);          
     fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);             
     fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);                
     fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);          
     fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);             
     fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);             
     fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);                
     fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag); 
     fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);    
     fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);             
     fChain->SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03, &b_mcCalIsoDR03); 
     fChain->SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03, &b_mcTrkIsoDR03); 
     fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04); 
     fChain->SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04, &b_mcTrkIsoDR04); 
     fChain->SetBranchAddress("genMET", &genMET, &b_genMET);          
     fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi); 
     fChain->SetBranchAddress("pfMETuncorrected", &pfMETuncorrected, &b_pfMETuncorrected); 
     fChain->SetBranchAddress("pfMETPhiuncorrected", &pfMETPhiuncorrected, &b_pfMETPhiuncorrected);
     fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);                                                                                
     fChain->SetBranchAddress("jetHadFlvr", &jetHadFlvr, &b_jetHadFlvr);                                                                                   
     fChain->SetBranchAddress("jetGenJetIndex", &jetGenJetIndex, &b_jetGenJetIndex);                                                                       
     fChain->SetBranchAddress("jetGenJetEn", &jetGenJetEn, &b_jetGenJetEn);                                                                                
     fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);                                                                                
     fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);                                                                             
     fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);                                                                             
     fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);                                                                       
     fChain->SetBranchAddress("jetGenEn", &jetGenEn, &b_jetGenEn);                                                                                         
     fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);                                                                                         
     fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);                                                                                      
     fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);                                                                                      
     fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
     fChain->SetBranchAddress("AK8JetPartonID", &AK8JetPartonID, &b_AK8JetPartonID);                                                                       
     fChain->SetBranchAddress("AK8JetHadFlvr", &AK8JetHadFlvr, &b_AK8JetHadFlvr);                                                                          
     fChain->SetBranchAddress("AK8JetGenJetIndex", &AK8JetGenJetIndex, &b_AK8JetGenJetIndex);                                                              
     fChain->SetBranchAddress("AK8JetGenJetEn", &AK8JetGenJetEn, &b_AK8JetGenJetEn);                                                                       
     fChain->SetBranchAddress("AK8JetGenJetPt", &AK8JetGenJetPt, &b_AK8JetGenJetPt);                                                                       
     fChain->SetBranchAddress("AK8JetGenJetEta", &AK8JetGenJetEta, &b_AK8JetGenJetEta);                                                                    
     fChain->SetBranchAddress("AK8JetGenJetPhi", &AK8JetGenJetPhi, &b_AK8JetGenJetPhi);                                                                    
     fChain->SetBranchAddress("AK8JetGenPartonID", &AK8JetGenPartonID, &b_AK8JetGenPartonID);                                                              
     fChain->SetBranchAddress("AK8JetGenEn", &AK8JetGenEn, &b_AK8JetGenEn);                                                                                
     fChain->SetBranchAddress("AK8JetGenPt", &AK8JetGenPt, &b_AK8JetGenPt);                                                                                
     fChain->SetBranchAddress("AK8JetGenEta", &AK8JetGenEta, &b_AK8JetGenEta);                                                                             
     fChain->SetBranchAddress("AK8JetGenPhi", &AK8JetGenPhi, &b_AK8JetGenPhi);                                                                             
     fChain->SetBranchAddress("AK8JetGenPartonMomID", &AK8JetGenPartonMomID, &b_AK8JetGenPartonMomID);

}

Bool_t postAnalyzer_Gen::FillSigHistogramsGen(int ptbin, int selbin, int photonIndex, int genIndex,  double weight){
 Float_t uncorphoet = ((*phoSCRawE)[photonIndex]/TMath::CosH((*phoSCEta)[photonIndex]));
 Float_t genphoet = mcPt->at(genIndex);

 Float_t etres  = (uncorphoet - genphoet)/genphoet;
 Float_t metres = (pfMET      - genMET)/genMET;

  h_gen_et[ptbin].Fill( genphoet, weight );
  h_gen_etres[ptbin].Fill( etres, weight );
  h_gen_met[ptbin].Fill( genMET, weight );
  h_gen_metres[ptbin].Fill( metres, weight );

 return kTRUE;
}

Bool_t postAnalyzer_Gen::WriteHistogramsGen(int ptbin, int sysbin){

  h_gen_et[ptbin]    .Write();
  h_gen_etres[ptbin] .Write();
  h_gen_met[ptbin]   .Write();
  h_gen_metres[ptbin].Write();

 return kTRUE;

}

//-------------------------callFillSigHistGen
void postAnalyzer_Gen::callFillSigHistGen(int selbin, int lastptbin, int inclptbin, int candphotonindex, int genIndex, float event_weight){
 Float_t uncorrectedPhoEt = ((*phoSCRawE)[candphotonindex]/TMath::CosH((*phoSCEta)[candphotonindex]));
 for(unsigned int ptb=0; ptb<lastptbin-2; ++ptb){ // break into pT bins
  if(
     ( uncorrectedPhoEt > ptbins[ptb]) &&
     ( uncorrectedPhoEt < ptbins[ptb+1])
    ){
   FillSigHistogramsGen(ptb, selbin, candphotonindex, genIndex, event_weight);
  } // end if passes pt cuts then fill
 } // end pt bin loop
 if(  // do an inclusive pT plot from bins
    ( uncorrectedPhoEt > ptbins[0]) &&
    ( uncorrectedPhoEt < ptbins[inclptbin])
   ){
  FillSigHistogramsGen(lastptbin-1, selbin, candphotonindex, genIndex, event_weight);
 }
 // and one fully inclusive in pT
 FillSigHistogramsGen(lastptbin, selbin, candphotonindex, genIndex, event_weight);
 return;
}


#endif
