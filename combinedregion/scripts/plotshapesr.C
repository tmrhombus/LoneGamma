#define plotshapesr_cxx 
#include "plotshapesr.h"

void plotshapesr::Loop()
{

 char* submitbase;
 submitbase = getenv ("submitbase");
 char* version;
 version = getenv("version");
 Tsubmitbase = TString(submitbase);
 Tversion = TString(version);

 inpath = TString(Tsubmitbase+"/gitignore/"+Tversion+"/analyzed");
 outpath = TString(Tsubmitbase+"/gitignore/"+Tversion+"/plots");

 regions.push_back("GenSignal");
// regions.push_back("Signal");
// regions.push_back("WlnG");
// regions.push_back("ZllG");

// variablenames.push_back("etres"  );
// variablenames.push_back("metres" );

 variablenames.push_back("h_sig_et5"     );     
 variablenames.push_back("h_sig_et15"    );    
 variablenames.push_back("h_sig_et25"    );    
 variablenames.push_back("h_sig_et55"    );    
 variablenames.push_back("h_sig_et75"    );    
 variablenames.push_back("h_sig_et165"   );   
 variablenames.push_back("h_sig_et275"   );   
 variablenames.push_back("h_sig_met5"    );    
 variablenames.push_back("h_sig_met15"   );   
 variablenames.push_back("h_sig_met25"   );   
 variablenames.push_back("h_sig_met55"   );   
 variablenames.push_back("h_sig_met75"   );   
 variablenames.push_back("h_sig_met165"  );  
 variablenames.push_back("h_sig_met275"  );  
 variablenames.push_back("h_gen_et_"     );     
 variablenames.push_back("h_gen_etres_"  );  
 variablenames.push_back("h_gen_met_"    );    
 variablenames.push_back("h_gen_metres_" ); 

// variablenames.push_back("uncorret"  );
// variablenames.push_back("et"        ); 
// variablenames.push_back("eta"       ); 
// variablenames.push_back("sieieF5x5" ); 
// variablenames.push_back("pfMET"     ); 
// variablenames.push_back("leptoMET"  ); 
// variablenames.push_back("dilep_mass"); 
// variablenames.push_back("dimu_mass" ); 
// variablenames.push_back("diele_mass"); 

 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);
 rebins.push_back(1);

// selectionnames.push_back("m110");
// selectionnames.push_back("m170");

 ptrangenames.push_back("allpt");
 //ptrangenames.push_back("175to1000");

 dolog = true;
 //dolog = false;

 dobygev = true;
 //dobygev = false;
 
 doscale = true;

 for(unsigned int r=0; r<regions.size(); ++r){
  TString region = regions[r];

  name_GJets             =  "analyzed"+region+"_GJets.root";
  name_TTGJets           =  "analyzed"+region+"_TTGJets.root";
  name_TGJets            =  "analyzed"+region+"_TGJets.root";
  name_ZllGJets          =  "analyzed"+region+"_ZllGJets.root";
  name_ZllJets           =  "analyzed"+region+"_ZllJets.root";
  name_GGJets            =  "analyzed"+region+"_GGJets.root";
  name_WWG               =  "analyzed"+region+"_WWG.root";
  name_WZ                =  "analyzed"+region+"_WZ.root";
  name_WlnGJets          =  "analyzed"+region+"_WlnGJets.root";
  name_Wmn               =  "analyzed"+region+"_Wmn.root";
  name_Wtn               =  "analyzed"+region+"_Wtn.root";
  name_ZZ                =  "analyzed"+region+"_ZZ.root";
  name_ZnnGJets          =  "analyzed"+region+"_ZnnGJets.root";
  
  file_GJets             = new TFile(inpath+"/"+name_GJets,"READ");
  file_TTGJets           = new TFile(inpath+"/"+name_TTGJets,"READ");
  file_TGJets            = new TFile(inpath+"/"+name_TGJets,"READ");
  file_ZllGJets          = new TFile(inpath+"/"+name_ZllGJets,"READ");
  file_ZllJets           = new TFile(inpath+"/"+name_ZllJets,"READ");
  file_GGJets            = new TFile(inpath+"/"+name_GGJets,"READ");
  file_WWG               = new TFile(inpath+"/"+name_WWG,"READ");
  file_WZ                = new TFile(inpath+"/"+name_WZ,"READ");
  file_WlnGJets          = new TFile(inpath+"/"+name_WlnGJets,"READ");
  file_Wmn               = new TFile(inpath+"/"+name_Wmn,"READ");
  file_Wtn               = new TFile(inpath+"/"+name_Wtn,"READ");
  file_ZZ                = new TFile(inpath+"/"+name_ZZ,"READ");
  file_ZnnGJets          = new TFile(inpath+"/"+name_ZnnGJets,"READ");
 
  extraname = "";
  if(dolog){extraname+="_log";}
  if(dobygev){extraname+="_bygev";}
 
  TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500); 
  gStyle->SetOptStat(0);
  gPad->SetLogy(dolog);
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);
 
  TText* title = new TText(1,1,"") ;
  title->SetTextSize(0.04);
  title->SetTextColor(kBlack);
  title->SetTextAlign(11);
  title->SetTextFont(62);
 
  TText* extra = new TText(1,1,"") ;
  extra->SetTextSize(0.04);
  extra->SetTextColor(kBlack);
  extra->SetTextAlign(11);
  //extra->SetTextAlign(13);
  extra->SetTextFont(52);
 
  TText* lumi = new TText(1,1,"") ;
  lumi->SetTextSize(0.03);
  lumi->SetTextColor(kBlack);
  lumi->SetTextAlign(31);
  lumi->SetTextFont(42);
 
 //  Int_t ci10 = 1010;
 //  TColor *color10 = new TColor(ci10,255./255,127./255, 0./255);
 
   for(unsigned int v=0; v<variablenames.size(); ++v){
   TString variablename = variablenames[v];
   Int_t rebin = rebins[v];
//    for(unsigned int s=0; s<selectionnames.size(); ++s){
//    TString selectionname = selectionnames[s];
     for(unsigned int p=0; p<ptrangenames.size(); ++p){
     TString ptrangename = ptrangenames[p];
 
     canvas->Clear();
 
     TString tailname = variablename+ptrangename; //+"_"+selectionname;
     //TString tailname = variablename+"_"+ptrangename; //+"_"+selectionname;
     outname = "plot"+region+"_"+tailname+extraname; 
     histname = tailname;
     //histname = "h_sig_"+tailname;
     file_out =  new TFile(outpath+"/"+outname+".root","RECREATE");
 
     printf("------------------------\n");
     printf(outname+"\n");
     printf(histname+"\n");
 
     //raw_data_obs = (TH1F*)file_SinglePhotonData->Get(histname)->Clone("raw_data_obs");
 
     h_raw_GJets    = (TH1F*)file_GJets             ->Get(histname)->Clone("h_raw_GJets")    ;   
     h_raw_TTGJets  = (TH1F*)file_TTGJets           ->Get(histname)->Clone("h_raw_TTGJets")  ;   
     h_raw_TGJets   = (TH1F*)file_TGJets            ->Get(histname)->Clone("h_raw_TGJets")   ;   
     h_raw_ZllGJets = (TH1F*)file_ZllGJets          ->Get(histname)->Clone("h_raw_ZllGJets") ;   
     h_raw_ZllJets  = (TH1F*)file_ZllJets           ->Get(histname)->Clone("h_raw_ZllJets")  ;  
     h_raw_GGJets   = (TH1F*)file_GGJets            ->Get(histname)->Clone("h_raw_GGJets")   ;  
     //h_raw_EFake    = (TH1F*)file_SinglePhotonEle   ->Get(histname)->Clone("h_raw_EFake")    ;   
     //h_raw_JFake    = (TH1F*)file_SinglePhotonJet   ->Get(histname)->Clone("h_raw_JFake")    ;   
     //h_raw_Halo     = (TH1F*)file_SinglePhotonHalo  ->Get(histname)->Clone("h_raw_Halo")     ;   
     //h_raw_Spike    = (TH1F*)file_SinglePhotonSpike ->Get(histname)->Clone("h_raw_Spike")    ;   
     h_raw_WWG      = (TH1F*)file_WWG               ->Get(histname)->Clone("h_raw_WWG")      ;     
     h_raw_WZ       = (TH1F*)file_WZ                ->Get(histname)->Clone("h_raw_WZ")       ;      
     h_raw_WlnGJets = (TH1F*)file_WlnGJets          ->Get(histname)->Clone("h_raw_WlnGJets") ;
     h_raw_Wmn      = (TH1F*)file_Wmn               ->Get(histname)->Clone("h_raw_Wmn")      ;     
     h_raw_Wtn      = (TH1F*)file_Wtn               ->Get(histname)->Clone("h_raw_Wtn")      ;     
     h_raw_ZZ       = (TH1F*)file_ZZ                ->Get(histname)->Clone("h_raw_ZZ")       ;      
     h_raw_ZnnGJets = (TH1F*)file_ZnnGJets          ->Get(histname)->Clone("h_raw_ZnnGJets") ;
 
     //raw_data_obs  ->SetMarkerStyle( 20 );  
     //raw_data_obs  ->SetMarkerSize( 1 );  
     //raw_data_obs  ->SetLineWidth( 2 );
 
     h_raw_GJets   ->SetLineColor( 1 ); 
     h_raw_TTGJets ->SetLineColor( 2 ); 
     h_raw_TGJets  ->SetLineColor( 3 ); 
     h_raw_ZllGJets->SetLineColor( 4 ); 
     h_raw_ZllJets ->SetLineColor( 5 ); 
     h_raw_GGJets  ->SetLineColor( 6 );
     h_raw_WWG     ->SetLineColor( 7 );
     h_raw_WZ      ->SetLineColor( 8 );
     h_raw_WlnGJets->SetLineColor( 9 );
     h_raw_Wmn     ->SetLineColor( 10 );
     h_raw_Wtn     ->SetLineColor( 11 );
     h_raw_ZZ      ->SetLineColor( 12 );
     h_raw_ZnnGJets->SetLineColor( 13 );
 
     //raw_data_obs  ->Rebin( rebin ); 
     h_raw_GJets   ->Rebin( rebin ); 
     h_raw_TTGJets ->Rebin( rebin ); 
     h_raw_TGJets  ->Rebin( rebin ); 
     h_raw_ZllGJets->Rebin( rebin ); 
     h_raw_ZllJets ->Rebin( rebin ); 
     h_raw_GGJets  ->Rebin( rebin );
     //h_raw_EFake   ->Rebin( rebin );
     //h_raw_JFake   ->Rebin( rebin );
     //h_raw_Halo    ->Rebin( rebin );
     //h_raw_Spike   ->Rebin( rebin );
     h_raw_WWG     ->Rebin( rebin );
     h_raw_WZ      ->Rebin( rebin );
     h_raw_WlnGJets->Rebin( rebin );
     h_raw_Wmn     ->Rebin( rebin );
     h_raw_Wtn     ->Rebin( rebin );
     h_raw_ZZ      ->Rebin( rebin );
     h_raw_ZnnGJets->Rebin( rebin );
 
     h_raw_GJets   ->SetLineWidth( 2 ); 
     h_raw_TTGJets ->SetLineWidth( 2 ); 
     h_raw_TGJets  ->SetLineWidth( 2 ); 
     h_raw_ZllGJets->SetLineWidth( 2 ); 
     h_raw_ZllJets ->SetLineWidth( 2 ); 
     h_raw_GGJets  ->SetLineWidth( 2 );
     //h_raw_EFake   ->SetLineWidth( 2 );
     //h_raw_JFake   ->SetLineWidth( 2 );
     //h_raw_Halo    ->SetLineWidth( 2 );
     //h_raw_Spike   ->SetLineWidth( 2 );
     h_raw_WWG     ->SetLineWidth( 2 );
     h_raw_WZ      ->SetLineWidth( 2 );
     h_raw_WlnGJets->SetLineWidth( 2 );
     h_raw_Wmn     ->SetLineWidth( 2 );
     h_raw_Wtn     ->SetLineWidth( 2 );
     h_raw_ZZ      ->SetLineWidth( 2 );
     h_raw_ZnnGJets->SetLineWidth( 2 );

     if(doscale){
      h_raw_GJets   ->Scale( 1./max(0.0001,h_raw_GJets   ->Integral(1,-1)) ); 
      h_raw_TTGJets ->Scale( 1./max(0.0001,h_raw_TTGJets ->Integral(1,-1)) ); 
      h_raw_TGJets  ->Scale( 1./max(0.0001,h_raw_TGJets  ->Integral(1,-1)) ); 
      h_raw_ZllGJets->Scale( 1./max(0.0001,h_raw_ZllGJets->Integral(1,-1)) ); 
      h_raw_ZllJets ->Scale( 1./max(0.0001,h_raw_ZllJets ->Integral(1,-1)) ); 
      h_raw_GGJets  ->Scale( 1./max(0.0001,h_raw_GGJets  ->Integral(1,-1)) );
      h_raw_WWG     ->Scale( 1./max(0.0001,h_raw_WWG     ->Integral(1,-1)) );
      h_raw_WZ      ->Scale( 1./max(0.0001,h_raw_WZ      ->Integral(1,-1)) );
      h_raw_WlnGJets->Scale( 1./max(0.0001,h_raw_WlnGJets->Integral(1,-1)) );
      h_raw_Wmn     ->Scale( 1./max(0.0001,h_raw_Wmn     ->Integral(1,-1)) );
      h_raw_Wtn     ->Scale( 1./max(0.0001,h_raw_Wtn     ->Integral(1,-1)) );
      h_raw_ZZ      ->Scale( 1./max(0.0001,h_raw_ZZ      ->Integral(1,-1)) );
      h_raw_ZnnGJets->Scale( 1./max(0.0001,h_raw_ZnnGJets->Integral(1,-1)) );
     }
 
     ///// h_GWZtt
     ///h_raw_GJets   ->SetFillColor( kRed-10 );  
     ///h_raw_Wmn     ->SetFillColor( kRed-10 );
     ///h_raw_Wtn     ->SetFillColor( kRed-10 );
     ///h_raw_ZllGJets->SetFillColor( kRed-10 );
     ///h_raw_TTGJets ->SetFillColor( kRed-10 ); 
     ///h_raw_TGJets  ->SetFillColor( kRed-10 ); 
     ///h_raw_ZllJets ->SetFillColor( kRed-10 );  
     ///// h_VV
     ///h_raw_WZ      ->SetFillColor(TColor::GetColor("#CCCA2A"));
     ///h_raw_ZZ      ->SetFillColor(TColor::GetColor("#CCCA2A"));
     ///h_raw_WWG     ->SetFillColor(TColor::GetColor("#CCCA2A"));

     ///h_raw_Halo    ->SetFillColor( kAzure-9 );
     ///h_raw_Spike   ->SetFillColor(TColor::GetColor("#FF6633"));
     ///h_raw_JFake   ->SetFillColor(TColor::GetColor("#FFFFCC"));
     ///h_raw_EFake   ->SetFillColor( kBlue-8   );
     ///h_raw_WlnGJets->SetFillColor( kRed-6    );  
     ///h_raw_ZnnGJets->SetFillColor( kOrange-4 );
     ///h_raw_GGJets  ->SetFillColor( kBlue-7   );

     ///// h_GWZtt
     ///h_raw_GWZtt = (TH1F*)h_raw_GJets->Clone("h_raw_GWZtt") ;
     ///h_raw_GWZtt -> Add(h_raw_Wmn     ) ;
     ///h_raw_GWZtt -> Add(h_raw_Wtn     ) ;
     ///h_raw_GWZtt -> Add(h_raw_ZllGJets) ;
     ///h_raw_GWZtt -> Add(h_raw_TTGJets ) ;
     ///h_raw_GWZtt -> Add(h_raw_TGJets  ) ;
     ///h_raw_GWZtt -> Add(h_raw_ZllJets ) ;

     ///// h_VV
     ///h_raw_VV = (TH1F*)h_raw_WZ->Clone("h_raw_VV") ;
     ///h_raw_VV -> Add(h_raw_ZZ) ;
     ///h_raw_VV -> Add(h_raw_WWG) ;

     ///h_raw_Halo->Scale(5.5/h_raw_Halo->Integral(1,-1));
     ///h_raw_Spike->Scale(8.5/h_raw_Spike->Integral(1,-1));

     //if(lepton=="mu"){printf("lepton mu\n"  ); h_EFake->Scale( 0.5/h_EFake->Integral() );    }
     //if(lepton=="ele"){printf("lepton ele\n"); h_EFake->Scale( 0.1225/h_EFake->Integral() ); }
     ///h_JFake->Scale( 5.9/h_JFake->Integral(-1,-1) );
     ///h_EFake->Scale( 52.7/h_EFake->Integral(-1,-1) );
 
///     Double_t xbins[6] = {175,190,250,400,700,1000};
///     TH1 *data_obs   = raw_data_obs  ->Rebin(5,"data_obs",xbins);    
///     TH1 *h_GJets    = h_raw_GJets   ->Rebin(5,"h_GJets",   xbins);    
///     TH1 *h_TTGJets  = h_raw_TTGJets ->Rebin(5,"h_TTGJets", xbins);    
///     TH1 *h_TGJets   = h_raw_TGJets  ->Rebin(5,"h_TGJets",  xbins);    
///     TH1 *h_WlnGJets = h_raw_WlnGJets->Rebin(5,"h_WlnGJets",xbins);    
///     TH1 *h_ZllGJets = h_raw_ZllGJets->Rebin(5,"h_ZllGJets",xbins);    
///     TH1 *h_ZllJets  = h_raw_ZllJets ->Rebin(5,"h_ZllJets", xbins);    
///     TH1 *h_GGJets   = h_raw_GGJets  ->Rebin(5,"h_GGJets",  xbins);
///     TH1 *h_EFake    = h_raw_EFake   ->Rebin(5,"h_EFake",   xbins);
///     TH1 *h_JFake    = h_raw_JFake   ->Rebin(5,"h_JFake",   xbins);
///     TH1 *h_Halo     = h_raw_Halo    ->Rebin(5,"h_Halo",    xbins);
///     TH1 *h_Spike    = h_raw_Spike   ->Rebin(5,"h_Spike",   xbins);
///     TH1 *h_WWG      = h_raw_WWG     ->Rebin(5,"h_WWG",     xbins);
///     TH1 *h_WZ       = h_raw_WZ      ->Rebin(5,"h_WZ",      xbins);
///     TH1 *h_Wmn      = h_raw_Wmn     ->Rebin(5,"h_Wmn",     xbins);
///     TH1 *h_Wtn      = h_raw_Wtn     ->Rebin(5,"h_Wtn",     xbins);
///     TH1 *h_ZZ       = h_raw_ZZ      ->Rebin(5,"h_ZZ",      xbins);
///     TH1 *h_ZnnGJets = h_raw_ZnnGJets->Rebin(5,"h_ZnnGJets",xbins);
///     TH1 *h_GWZtt    = h_raw_GWZtt   ->Rebin(5,"h_GWZtt",   xbins);
///     TH1 *h_VV       = h_raw_VV      ->Rebin(5,"h_VV",      xbins);
/// 
///     if(dobygev){
///      if(variablename=="et"      ||
///         variablename=="uncorret"||
///         variablename=="leptoMET"
///        ){
///     //  Double_t xbins[6] = {175,190,250,400,700,1000};
///     //   TH1 *hr_data_obs = data_obs  ->Rebin(5,"hr_data_obs",xbins);    
///     //   TH1 *hr_GJets    = h_GJets   ->Rebin(5,"hr_GJets",   xbins);    
///     //   TH1 *hr_TTGJets  = h_TTGJets ->Rebin(5,"hr_TTGJets", xbins);    
///     //   TH1 *hr_TGJets   = h_TGJets  ->Rebin(5,"hr_TGJets",  xbins);    
///     //   TH1 *hr_WlnGJets = h_WlnGJets->Rebin(5,"hr_WlnGJets",xbins);    
///     //   TH1 *hr_ZllGJets = h_ZllGJets->Rebin(5,"hr_ZllGJets",xbins);    
///     //   TH1 *hr_ZllJets  = h_ZllJets ->Rebin(5,"hr_ZllJets", xbins);    
///     //   TH1 *hr_GGJets   = h_GGJets  ->Rebin(5,"hr_GGJets",  xbins);
///     //   TH1 *hr_EFake    = h_EFake   ->Rebin(5,"hr_EFake",   xbins);
///     //   TH1 *hr_JFake    = h_JFake   ->Rebin(5,"hr_JFake",   xbins);
///     //   TH1 *hr_WWG      = h_WWG     ->Rebin(5,"hr_WWG",     xbins);
///     //   TH1 *hr_WZ       = h_WZ      ->Rebin(5,"hr_WZ",      xbins);
///     //   TH1 *hr_Wmn      = h_Wmn     ->Rebin(5,"hr_Wmn",     xbins);
///     //   TH1 *hr_Wtn      = h_Wtn     ->Rebin(5,"hr_Wtn",     xbins);
///     //   TH1 *hr_ZZ       = h_ZZ      ->Rebin(5,"hr_ZZ",      xbins);
///     //   TH1 *hr_ZnnGJets = h_ZnnGJets->Rebin(5,"hr_ZnnGJets",xbins);
///  //std::cout<<hr_data_obs->GetSize()<<std::endl;
///       for(int i=0; i<5; i++){
///        int bnr = i+1;
///        double divby = xbins[i+1] - xbins[i];
///      
///        double cont_data_obs = data_obs  ->GetBinContent(bnr); 
///        double err_data_obs  = data_obs  ->GetBinError(bnr); 
///        double err_ratio = err_data_obs/cont_data_obs;
///    
///        double cont_GJets    = h_GJets    ->GetBinContent(bnr); 
///        double cont_TTGJets  = h_TTGJets  ->GetBinContent(bnr); 
///        double cont_TGJets   = h_TGJets   ->GetBinContent(bnr); 
///        double cont_WlnGJets = h_WlnGJets ->GetBinContent(bnr); 
///        double cont_ZllGJets = h_ZllGJets ->GetBinContent(bnr); 
///        double cont_ZllJets  = h_ZllJets  ->GetBinContent(bnr); 
///        double cont_GGJets   = h_GGJets   ->GetBinContent(bnr);
///        double cont_EFake    = h_EFake    ->GetBinContent(bnr);
///        double cont_JFake    = h_JFake    ->GetBinContent(bnr);
///        double cont_Halo     = h_Halo     ->GetBinContent(bnr);
///        double cont_Spike    = h_Spike    ->GetBinContent(bnr);
///        double cont_WWG      = h_WWG      ->GetBinContent(bnr);
///        double cont_WZ       = h_WZ       ->GetBinContent(bnr);
///        double cont_Wmn      = h_Wmn      ->GetBinContent(bnr);
///        double cont_Wtn      = h_Wtn      ->GetBinContent(bnr);
///        double cont_ZZ       = h_ZZ       ->GetBinContent(bnr);
///        double cont_ZnnGJets = h_ZnnGJets ->GetBinContent(bnr);
///        double cont_GWZtt    = h_GWZtt    ->GetBinContent(bnr);
///        double cont_VV       = h_VV       ->GetBinContent(bnr);
///      
///        data_obs->SetBinContent(bnr, cont_data_obs / divby );    
///        h_GJets   ->SetBinContent(bnr, cont_GJets    / divby );    
///        h_TTGJets ->SetBinContent(bnr, cont_TTGJets  / divby );    
///        h_TGJets  ->SetBinContent(bnr, cont_TGJets   / divby );    
///        h_WlnGJets->SetBinContent(bnr, cont_WlnGJets / divby );    
///        h_ZllGJets->SetBinContent(bnr, cont_ZllGJets / divby );    
///        h_ZllJets ->SetBinContent(bnr, cont_ZllJets  / divby );    
///        h_GGJets  ->SetBinContent(bnr, cont_GGJets   / divby );
///        h_EFake   ->SetBinContent(bnr, cont_EFake    / divby );
///        h_JFake   ->SetBinContent(bnr, cont_JFake    / divby );
///        h_Halo    ->SetBinContent(bnr, cont_Halo     / divby );
///        h_Spike   ->SetBinContent(bnr, cont_Spike    / divby );
///        h_WWG     ->SetBinContent(bnr, cont_WWG      / divby );
///        h_WZ      ->SetBinContent(bnr, cont_WZ       / divby );
///        h_Wmn     ->SetBinContent(bnr, cont_Wmn      / divby );
///        h_Wtn     ->SetBinContent(bnr, cont_Wtn      / divby );
///        h_ZZ      ->SetBinContent(bnr, cont_ZZ       / divby );
///        h_ZnnGJets->SetBinContent(bnr, cont_ZnnGJets / divby );
///        h_GWZtt   ->SetBinContent(bnr, cont_GWZtt    / divby );
///        h_VV      ->SetBinContent(bnr, cont_VV       / divby );
///      
///        data_obs ->SetBinError(bnr, cont_data_obs*err_ratio / divby );    
///      
///       }
///      }
///     } // dobygev
// 
//     THStack *thestack = new THStack("thestack","");
//     thestack->Add(h_GWZtt   );  
//     thestack->Add(h_VV      );  
//     thestack->Add(h_Halo    );
//     thestack->Add(h_Spike   );
//     thestack->Add(h_JFake   );
//     thestack->Add(h_EFake   );
//     thestack->Add(h_WlnGJets);
//     thestack->Add(h_ZnnGJets);
//     thestack->Add(h_GGJets  );

     TLegend *leg = new TLegend(0.20,0.7,0.88,0.85 );
     //TLegend *leg = new TLegend(0.60,0.5,0.88,0.88 );
     leg-> SetNColumns(2);
     leg->SetBorderSize(0);
     leg->SetFillColor(kWhite);
     leg->AddEntry( h_raw_GJets   , "#gamma+jet","l"); 
     leg->AddEntry( h_raw_TTGJets , "t#bar{t}#gamma","l"); 
     leg->AddEntry( h_raw_TGJets  , "t#gamma","l"); 
     leg->AddEntry( h_raw_ZllGJets, "Z(ll)#gamma","l"); 
     leg->AddEntry( h_raw_ZllJets , "Z(ll)+jets","l"); 
     leg->AddEntry( h_raw_GGJets  , "#gamma#gamma+jets","l"); 
     leg->AddEntry( h_raw_WWG     , "WW#gamma","l"); 
     leg->AddEntry( h_raw_WZ      , "WZ","l"); 
     leg->AddEntry( h_raw_WlnGJets, "W(l#nu)#gamma+jets","l"); 
     leg->AddEntry( h_raw_Wmn     , "W(#mu#nu)","l"); 
     leg->AddEntry( h_raw_Wtn     , "W(#tau#nu)","l"); 
     leg->AddEntry( h_raw_ZZ      , "ZZ","l"); 
     leg->AddEntry( h_raw_ZnnGJets, "Z(#nu#nu)#gamma","l"); 
 
     canvas->cd();
 
     h_raw_ZnnGJets->SetXTitle(TString(h_raw_ZnnGJets->GetTitle())+" [GeV]");
     h_raw_ZnnGJets->SetTitle("");
     h_raw_ZnnGJets->SetYTitle("Events / Bin");
     //hr_h_raw_ZnnGJets->SetMaximum(100);
     if(dobygev){
      h_raw_ZnnGJets->SetYTitle("Events / GeV");
      //hr_h_raw_ZnnGJets->SetMaximum(10);
     }
     if(dolog){
      h_raw_ZnnGJets->SetMaximum(200);
      h_raw_ZnnGJets->SetMinimum(0.0005);
     }
     if(doscale){
      h_raw_ZnnGJets->SetMaximum(1);
     }
     h_raw_ZnnGJets->GetYaxis()->SetTitleOffset(1.4);
 
     h_raw_ZnnGJets->Draw("");
     h_raw_GJets   ->Draw("hist,sames");
     h_raw_TTGJets ->Draw("hist,sames");
     h_raw_TGJets  ->Draw("hist,sames");
     h_raw_ZllGJets->Draw("hist,sames");
     h_raw_ZllJets ->Draw("hist,sames");
     h_raw_GGJets  ->Draw("hist,sames");
     h_raw_WWG     ->Draw("hist,sames");
     h_raw_WZ      ->Draw("hist,sames");
     h_raw_WlnGJets->Draw("hist,sames");
     h_raw_Wmn     ->Draw("hist,sames");
     h_raw_Wtn     ->Draw("hist,sames");
     h_raw_ZZ      ->Draw("hist,sames");
     h_raw_ZnnGJets->Draw("hist,sames");
     h_raw_ZnnGJets->GetXaxis()->Draw();
     leg->Draw("same");
 
     //title->DrawTextNDC(0.17,0.87,"CMS");
     //extra->DrawTextNDC(0.17,0.81,"Preliminary");
     title->DrawTextNDC(0.13,0.91,"CMS");
     extra->DrawTextNDC(0.23,0.91,"Simulation");
     lumi->DrawTextNDC(0.9,0.91,"40 /fb (13 TeV)");
     //lumi->DrawTextNDC(0.9,0.91,"12.9 /fb (13 TeV)");

     gPad->Update();
     gPad->RedrawAxis();
 
     canvas->SaveAs(outpath+"/"+outname+".pdf");
 
     h_raw_GJets    ->Delete();   
     h_raw_TTGJets  ->Delete();   
     h_raw_TGJets   ->Delete();   
     h_raw_ZllGJets ->Delete();   
     h_raw_ZllJets  ->Delete();   
     h_raw_GGJets   ->Delete();   
     h_raw_WWG      ->Delete();    
     h_raw_WZ       ->Delete();    
     h_raw_WlnGJets ->Delete();
     h_raw_Wmn      ->Delete();    
     h_raw_Wtn      ->Delete();    
     h_raw_ZZ       ->Delete();    
     h_raw_ZnnGJets ->Delete();
 
     file_out->Close();
 
     } //  for(unsigned int p=0; p<ptrangenames.size(); ++p)
    //} //  for(unsigned int s=0; s<selectionnames.size(); ++s)
   } //  for(unsigned int v=0; v<variablenames.size(); ++v)
 }
}
