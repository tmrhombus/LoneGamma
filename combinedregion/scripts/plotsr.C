#define plotsr_cxx 
#include "plotsr.h"

void plotsr::Loop()
{

 char* submitbase;
 submitbase = getenv ("submitbase");
 char* version;
 version = getenv("version");
 Tsubmitbase = TString(submitbase);
 Tversion = TString(version);

 inpath = TString(Tsubmitbase+"/gitignore/"+Tversion+"/analyzed");
 outpath = TString(Tsubmitbase+"/gitignore/"+Tversion+"/plots");

//// regions.push_back("GenSignal");
// regions.push_back("Signal");
 regions.push_back("WlnG");
 regions.push_back("ZllG");

// variablenames.push_back("etres"  );  
// variablenames.push_back("metres" ); 

// variablenames.push_back("et5");  
// variablenames.push_back("et15"); 
// variablenames.push_back("et25"); 
// variablenames.push_back("et55"); 
// variablenames.push_back("et75"); 
// variablenames.push_back("et165");
// variablenames.push_back("et275");

 variablenames.push_back("uncorret"  );
// variablenames.push_back("pfMET"     ); 
// variablenames.push_back("leptoMET"  ); 
//// variablenames.push_back("et"        ); 
//// variablenames.push_back("eta"       ); 
//// variablenames.push_back("sieieF5x5" ); 
// variablenames.push_back("dilep_mass"); 
// variablenames.push_back("dimu_mass" ); 
// variablenames.push_back("diele_mass"); 

   selectionnames.push_back("");
   selectionnames.push_back("_JERUp");    // jet energy resolution
   selectionnames.push_back("_JERDown");
   selectionnames.push_back("_JESUp");    // jet energy scale
   selectionnames.push_back("_JESDown");
   selectionnames.push_back("_MESUp");    // muon energy scale
   selectionnames.push_back("_MESDown");
   selectionnames.push_back("_EESUp");    // electron energy scale
   selectionnames.push_back("_EESDown");
   selectionnames.push_back("_PESUp");    // photon energy scale 1.5%
   selectionnames.push_back("_PESDown");
   selectionnames.push_back("_TESUp");    // tau energy scale
   selectionnames.push_back("_TESDown");
   selectionnames.push_back("_UESUp");    // unclustered energy scale
   selectionnames.push_back("_UESDown");
   selectionnames.push_back("_EWKUp");    // electroweak correction
   selectionnames.push_back("_EWKDown");

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

// rebins.push_back(10);
// rebins.push_back(10);
// rebins.push_back(5);
// rebins.push_back(1);
// rebins.push_back(5);
// rebins.push_back(5);
// rebins.push_back(1);
// rebins.push_back(1);
// rebins.push_back(1);


 ptrangenames.push_back("allpt");
 //ptrangenames.push_back("175to1000");

 //dolog = true;
 dolog = false;

 //dobygev = true;
 dobygev = false;
 
 for(unsigned int r=0; r<regions.size(); ++r){
  TString region = regions[r];

  name_SinglePhotonData  =  "analyzed"+region+"_SinglePhotonData.root";
  name_GJets             =  "analyzed"+region+"_GJets.root";
  name_TTGJets           =  "analyzed"+region+"_TTGJets.root";
  name_TGJets            =  "analyzed"+region+"_TGJets.root";
  name_ZllGJets          =  "analyzed"+region+"_ZllGJets.root";
  name_ZllJets           =  "analyzed"+region+"_ZllJets.root";
  name_GGJets            =  "analyzed"+region+"_GGJets.root";
  name_SinglePhotonEle   =  "analyzed"+region+"_SinglePhotonEle.root";
  name_SinglePhotonJet   =  "analyzed"+region+"_SinglePhotonJet.root";
  name_SinglePhotonHalo  =  "analyzed"+region+"_SinglePhotonHalo.root";
  name_SinglePhotonSpike =  "analyzed"+region+"_SinglePhotonSpike.root";
  name_WWG               =  "analyzed"+region+"_WWG.root";
  name_WZ                =  "analyzed"+region+"_WZ.root";
  name_WlnGJets          =  "analyzed"+region+"_WlnGJets.root";
  name_Wmn               =  "analyzed"+region+"_Wmn.root";
  name_Wtn               =  "analyzed"+region+"_Wtn.root";
  name_ZZ                =  "analyzed"+region+"_ZZ.root";
  name_ZnnGJets          =  "analyzed"+region+"_ZnnGJets.root";
  
  file_SinglePhotonData  = new TFile(inpath+"/"+name_SinglePhotonData,"READ");
  file_SinglePhotonEle   = new TFile(inpath+"/"+name_SinglePhotonEle,"READ");
  file_SinglePhotonJet   = new TFile(inpath+"/"+name_SinglePhotonJet,"READ");
  file_SinglePhotonHalo  = new TFile(inpath+"/"+name_SinglePhotonHalo,"READ");
  file_SinglePhotonSpike = new TFile(inpath+"/"+name_SinglePhotonSpike,"READ");
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
    for(unsigned int p=0; p<ptrangenames.size(); ++p){
    TString ptrangename = ptrangenames[p];
    file_out =  new TFile(outpath+"/histos_"+region+"_"+variablename+"_"+ptrangename+".root","RECREATE");
     for(unsigned int s=0; s<selectionnames.size(); ++s){
     TString selectionname = selectionnames[s];
 
     canvas->Clear();
 
     TString tailname = variablename+"_"+ptrangename+selectionname;
     TString tailnamed = variablename+"_"+ptrangename;
     outname = "plot"+region+"_"+tailname+extraname; 
     //histname = "h_mu_"+tailname;
     //TString histnamed = "h_mu_"+tailnamed;
     ////histname = "h_ele_"+tailname;
     ////TString histnamed = "h_ele_"+tailnamed;
     histname = "h_sig_"+tailname;
     TString histnamed = "h_sig_"+tailnamed;

     ofstream log; 
     log.open (outpath+"/"+outname+".log");


//     for(unsigned int sn=0; sn<sysnames.size(); ++sn){
//       log<<boost::format("\nSysname: %s \n") % sysnames[sn];
//      for(unsigned int j=0; j<nptbins; ++j){
//        log<<boost::format(" %s\n") % ptbinnames[j];
//        log<<boost::format("  QCD Fraction: %0.3d +- %0.3d \n") % qcd_frac[sn][j] % qcd_frac_err[sn][j];
//        log<<boost::format("  Bin Center: %0.3d +- %0.3d \n") % bincenterss[sn][j] % binerrorss[sn][j];
//        log<<boost::format("  Ratio: %0.3d +- %0.3d \n") % ratioss[sn][j] % ratioerrorss[sn][j];
//      }
//        log<<boost::format(" y = (%0.3d +- %0.3d)x + (%0.3d +- %0.3d) \n") % ms[sn] % mes[sn] % bs[sn] % bes[sn] ;
//     }
//

     log<<boost::format("------------------------\n");
     log<<boost::format(outname+"\n");
     log<<boost::format(histname+"\n\n");

     printf("------------------------\n");
     printf(outname+"\n");
     printf(histname+"\n");
 
     raw_data_obs   = (TH1F*)file_SinglePhotonData->Get(histnamed)->Clone("raw_data_obs");
 
     h_raw_EFake    = (TH1F*)file_SinglePhotonEle   ->Get(histname)->Clone("h_raw_EFake")    ;   
     h_raw_JFake    = (TH1F*)file_SinglePhotonJet   ->Get(histname)->Clone("h_raw_JFake")    ;   
     h_raw_Halo     = (TH1F*)file_SinglePhotonHalo  ->Get(histname)->Clone("h_raw_Halo")     ;   
     h_raw_Spike    = (TH1F*)file_SinglePhotonSpike ->Get(histname)->Clone("h_raw_Spike")    ;   
     h_raw_GJets    = (TH1F*)file_GJets             ->Get(histname)->Clone("h_raw_GJets")    ;   
     h_raw_TTGJets  = (TH1F*)file_TTGJets           ->Get(histname)->Clone("h_raw_TTGJets")  ;   
     h_raw_TGJets   = (TH1F*)file_TGJets            ->Get(histname)->Clone("h_raw_TGJets")   ;   
     h_raw_ZllGJets = (TH1F*)file_ZllGJets          ->Get(histname)->Clone("h_raw_ZllGJets") ;   
     h_raw_ZllJets  = (TH1F*)file_ZllJets           ->Get(histname)->Clone("h_raw_ZllJets")  ;  
     h_raw_GGJets   = (TH1F*)file_GGJets            ->Get(histname)->Clone("h_raw_GGJets")   ;  
     h_raw_WWG      = (TH1F*)file_WWG               ->Get(histname)->Clone("h_raw_WWG")      ;     
     h_raw_WZ       = (TH1F*)file_WZ                ->Get(histname)->Clone("h_raw_WZ")       ;      
     h_raw_WlnGJets = (TH1F*)file_WlnGJets          ->Get(histname)->Clone("h_raw_WlnGJets") ;
     h_raw_Wmn      = (TH1F*)file_Wmn               ->Get(histname)->Clone("h_raw_Wmn")      ;     
     h_raw_Wtn      = (TH1F*)file_Wtn               ->Get(histname)->Clone("h_raw_Wtn")      ;     
     h_raw_ZZ       = (TH1F*)file_ZZ                ->Get(histname)->Clone("h_raw_ZZ")       ;      
     h_raw_ZnnGJets = (TH1F*)file_ZnnGJets          ->Get(histname)->Clone("h_raw_ZnnGJets") ;
 
//     printf("raw_data_obs    %6.1f \n", raw_data_obs  ->Integral(1,-1)); 
//     printf("h_raw_EFake     %6.1f \n", h_raw_EFake   ->Integral(1,-1)); 
//     printf("h_raw_JFake     %6.1f \n", h_raw_JFake   ->Integral(1,-1)); 
//     printf("h_raw_Halo      %6.1f \n", h_raw_Halo    ->Integral(1,-1)); 
//     printf("h_raw_Spike     %6.1f \n", h_raw_Spike   ->Integral(1,-1)); 
//     printf("h_raw_GJets     %6.1f \n", h_raw_GJets   ->Integral(1,-1)); 
//     printf("h_raw_TTGJets   %6.1f \n", h_raw_TTGJets ->Integral(1,-1)); 
//     printf("h_raw_TGJets    %6.1f \n", h_raw_TGJets  ->Integral(1,-1)); 
//     printf("h_raw_ZllGJets  %6.1f \n", h_raw_ZllGJets->Integral(1,-1)); 
//     printf("h_raw_ZllJets   %6.1f \n", h_raw_ZllJets ->Integral(1,-1)); 
//     printf("h_raw_GGJets    %6.1f \n", h_raw_GGJets  ->Integral(1,-1)); 
//     printf("h_raw_WWG       %6.1f \n", h_raw_WWG     ->Integral(1,-1)); 
//     printf("h_raw_WZ        %6.1f \n", h_raw_WZ      ->Integral(1,-1)); 
//     printf("h_raw_WlnGJets  %6.1f \n", h_raw_WlnGJets->Integral(1,-1)); 
//     printf("h_raw_Wmn       %6.1f \n", h_raw_Wmn     ->Integral(1,-1)); 
//     printf("h_raw_Wtn       %6.1f \n", h_raw_Wtn     ->Integral(1,-1)); 
//     printf("h_raw_ZZ        %6.1f \n", h_raw_ZZ      ->Integral(1,-1)); 
//     printf("h_raw_ZnnGJets  %6.1f \n", h_raw_ZnnGJets->Integral(1,-1)); 
 
 
     raw_data_obs  ->SetMarkerStyle( 20 );  
     raw_data_obs  ->SetMarkerSize( 1 );  
     raw_data_obs  ->SetLineWidth( 2 );
 
     raw_data_obs  ->SetLineColor( kBlack );
     h_raw_GJets   ->SetLineColor( kBlack ); 
     h_raw_TTGJets ->SetLineColor( kBlack ); 
     h_raw_TGJets  ->SetLineColor( kBlack ); 
     h_raw_ZllGJets->SetLineColor( kBlack ); 
     h_raw_ZllJets ->SetLineColor( kBlack ); 
     h_raw_GGJets  ->SetLineColor( kBlack );
     h_raw_EFake   ->SetLineColor( kBlack );
     h_raw_JFake   ->SetLineColor( kBlack );
     h_raw_Halo    ->SetLineColor( kBlack );
     h_raw_Spike   ->SetLineColor( kBlack );
     h_raw_WWG     ->SetLineColor( kBlack );
     h_raw_WZ      ->SetLineColor( kBlack );
     h_raw_WlnGJets->SetLineColor( kBlack );
     h_raw_Wmn     ->SetLineColor( kBlack );
     h_raw_Wtn     ->SetLineColor( kBlack );
     h_raw_ZZ      ->SetLineColor( kBlack );
     h_raw_ZnnGJets->SetLineColor( kBlack );
// 
     raw_data_obs  ->Rebin( rebin ); 
     h_raw_GJets   ->Rebin( rebin ); 
     h_raw_TTGJets ->Rebin( rebin ); 
     h_raw_TGJets  ->Rebin( rebin ); 
     h_raw_ZllGJets->Rebin( rebin ); 
     h_raw_ZllJets ->Rebin( rebin ); 
     h_raw_GGJets  ->Rebin( rebin );
     h_raw_EFake   ->Rebin( rebin );
     h_raw_JFake   ->Rebin( rebin );
     h_raw_Halo    ->Rebin( rebin );
     h_raw_Spike   ->Rebin( rebin );
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
     h_raw_EFake   ->SetLineWidth( 2 );
     h_raw_JFake   ->SetLineWidth( 2 );
     h_raw_Halo    ->SetLineWidth( 2 );
     h_raw_Spike   ->SetLineWidth( 2 );
     h_raw_WWG     ->SetLineWidth( 2 );
     h_raw_WZ      ->SetLineWidth( 2 );
     h_raw_WlnGJets->SetLineWidth( 2 );
     h_raw_Wmn     ->SetLineWidth( 2 );
     h_raw_Wtn     ->SetLineWidth( 2 );
     h_raw_ZZ      ->SetLineWidth( 2 );
     h_raw_ZnnGJets->SetLineWidth( 2 );
 
     // h_GWZtt
     h_raw_GJets   ->SetFillColor( kRed-10 );  
     h_raw_Wmn     ->SetFillColor( kRed-10 );
     h_raw_Wtn     ->SetFillColor( kRed-10 );
     h_raw_TTGJets ->SetFillColor( kRed-10 ); 
     h_raw_TGJets  ->SetFillColor( kRed-10 ); 
     h_raw_ZllJets ->SetFillColor( kRed-10 );  
     h_raw_GGJets  ->SetFillColor( kRed-10 );
     // h_VV
     h_raw_WZ      ->SetFillColor(TColor::GetColor("#CCCA2A"));
     h_raw_ZZ      ->SetFillColor(TColor::GetColor("#CCCA2A"));
     h_raw_WWG     ->SetFillColor(TColor::GetColor("#CCCA2A"));

     h_raw_Halo    ->SetFillColor( kAzure-9 );
     h_raw_Spike   ->SetFillColor(TColor::GetColor("#FF6633"));
     h_raw_JFake   ->SetFillColor(TColor::GetColor("#FFFFCC"));
     h_raw_EFake   ->SetFillColor( kBlue-8   );
     h_raw_WlnGJets->SetFillColor( kRed-6    );  
     h_raw_ZnnGJets->SetFillColor( kOrange-4 );
     h_raw_ZllGJets->SetFillColor( kBlue-7 );

     // h_GWZtt
     h_raw_GWZtt = (TH1F*)h_raw_GJets->Clone("h_raw_GWZtt") ;
     h_raw_GWZtt -> Add(h_raw_Wmn     ) ;
     h_raw_GWZtt -> Add(h_raw_Wtn     ) ;
     h_raw_GWZtt -> Add(h_raw_GGJets  ) ;
     h_raw_GWZtt -> Add(h_raw_TTGJets ) ;
     h_raw_GWZtt -> Add(h_raw_TGJets  ) ;
     h_raw_GWZtt -> Add(h_raw_ZllJets ) ;

     // h_VV
     h_raw_VV = (TH1F*)h_raw_WZ->Clone("h_raw_VV") ;
     h_raw_VV -> Add(h_raw_ZZ) ;
     h_raw_VV -> Add(h_raw_WWG) ;

     h_raw_Halo->Scale(5.5/h_raw_Halo->Integral(1,-1));
     h_raw_Spike->Scale(8.5/h_raw_Spike->Integral(1,-1));

     
     if(
       region=="WlnG" ||
       region=="ZllG"
      ){
       int nbinsx = h_raw_Halo->GetXaxis()->GetNbins();
       for(int ibin=1; ibin<=nbinsx+1; ibin++)
        {     
            h_raw_Halo->SetBinContent(ibin,0.0001);
            h_raw_Halo->SetBinError(ibin,0.0001);
            h_raw_Spike->SetBinContent(ibin,0.0001);
            h_raw_Spike->SetBinError(ibin,0.0001);
        }   
      }

//// for standard binning
//     TH1 *data_obs   = (TH1F*)raw_data_obs  ->Clone("data_obs");
//     TH1 *h_EFake    = (TH1F*)h_raw_EFake   ->Clone("h_EFake");
//     TH1 *h_JFake    = (TH1F*)h_raw_JFake   ->Clone("h_JFake");
//     TH1 *h_Halo     = (TH1F*)h_raw_Halo    ->Clone("h_Halo");
//     TH1 *h_Spike    = (TH1F*)h_raw_Spike   ->Clone("h_Spike");
//     TH1 *h_GJets    = (TH1F*)h_raw_GJets   ->Clone("h_GJets");
//     TH1 *h_TTGJets  = (TH1F*)h_raw_TTGJets ->Clone("h_TTGJets");
//     TH1 *h_TGJets   = (TH1F*)h_raw_TGJets  ->Clone("h_TGJets");
//     TH1 *h_WlnGJets = (TH1F*)h_raw_WlnGJets->Clone("h_WlnGJets");
//     TH1 *h_ZllGJets = (TH1F*)h_raw_ZllGJets->Clone("h_ZllGJets");
//     TH1 *h_ZllJets  = (TH1F*)h_raw_ZllJets ->Clone("h_ZllJets");
//     TH1 *h_GGJets   = (TH1F*)h_raw_GGJets  ->Clone("h_GGJets");
//     TH1 *h_WWG      = (TH1F*)h_raw_WWG     ->Clone("h_WWG");
//     TH1 *h_WZ       = (TH1F*)h_raw_WZ      ->Clone("h_WZ");
//     TH1 *h_Wmn      = (TH1F*)h_raw_Wmn     ->Clone("h_Wmn");
//     TH1 *h_Wtn      = (TH1F*)h_raw_Wtn     ->Clone("h_Wtn");
//     TH1 *h_ZZ       = (TH1F*)h_raw_ZZ      ->Clone("h_ZZ");
//     TH1 *h_ZnnGJets = (TH1F*)h_raw_ZnnGJets->Clone("h_ZnnGJets");
//     TH1 *h_GWZtt    = (TH1F*)h_raw_GWZtt   ->Clone("h_GWZtt");
//     TH1 *h_VV       = (TH1F*)h_raw_VV      ->Clone("h_VV");
 
 
     Double_t xbins[6] = {175,190,250,400,700,1000};
     TH1 *data_obs   = raw_data_obs  ->Rebin(5,"data_obs"+selectionname,  xbins);    
     TH1 *h_GJets    = h_raw_GJets   ->Rebin(5,"h_GJets"+selectionname,   xbins);    
     TH1 *h_TTGJets  = h_raw_TTGJets ->Rebin(5,"h_TTGJets"+selectionname, xbins);    
     TH1 *h_TGJets   = h_raw_TGJets  ->Rebin(5,"h_TGJets"+selectionname,  xbins);    
     TH1 *h_WlnGJets = h_raw_WlnGJets->Rebin(5,"h_WlnGJets"+selectionname,xbins);    
     TH1 *h_ZllGJets = h_raw_ZllGJets->Rebin(5,"h_ZllGJets"+selectionname,xbins);    
     TH1 *h_ZllJets  = h_raw_ZllJets ->Rebin(5,"h_ZllJets"+selectionname, xbins);    
     TH1 *h_GGJets   = h_raw_GGJets  ->Rebin(5,"h_GGJets"+selectionname,  xbins);
     TH1 *h_EFake    = h_raw_EFake   ->Rebin(5,"h_EFake"+selectionname,   xbins);
     TH1 *h_JFake    = h_raw_JFake   ->Rebin(5,"h_JFake"+selectionname,   xbins);
     TH1 *h_Halo     = h_raw_Halo    ->Rebin(5,"h_Halo"+selectionname,    xbins);
     TH1 *h_Spike    = h_raw_Spike   ->Rebin(5,"h_Spike"+selectionname,   xbins);
     TH1 *h_WWG      = h_raw_WWG     ->Rebin(5,"h_WWG"+selectionname,     xbins);
     TH1 *h_WZ       = h_raw_WZ      ->Rebin(5,"h_WZ"+selectionname,      xbins);
     TH1 *h_Wmn      = h_raw_Wmn     ->Rebin(5,"h_Wmn"+selectionname,     xbins);
     TH1 *h_Wtn      = h_raw_Wtn     ->Rebin(5,"h_Wtn"+selectionname,     xbins);
     TH1 *h_ZZ       = h_raw_ZZ      ->Rebin(5,"h_ZZ"+selectionname,      xbins);
     TH1 *h_ZnnGJets = h_raw_ZnnGJets->Rebin(5,"h_ZnnGJets"+selectionname,xbins);
     TH1 *h_GWZtt    = h_raw_GWZtt   ->Rebin(5,"h_GWZtt"+selectionname,   xbins);
     TH1 *h_VV       = h_raw_VV      ->Rebin(5,"h_VV"+selectionname,      xbins);

     log<<boost::format("             total   bin1   bin2   bin3   bin4   bin5   bin6\n");
     log<<boost::format("data_obs    %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %data_obs  ->Integral(1,-1) %data_obs  ->GetBinContent(1) %data_obs  ->GetBinContent(2) %data_obs  ->GetBinContent(3) %data_obs  ->GetBinContent(4) %data_obs  ->GetBinContent(5) %data_obs  ->GetBinContent(6) ;
     log<<boost::format("h_EFake     %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_EFake   ->Integral(1,-1) %h_EFake   ->GetBinContent(1) %h_EFake   ->GetBinContent(2) %h_EFake   ->GetBinContent(3) %h_EFake   ->GetBinContent(4) %h_EFake   ->GetBinContent(5) %h_EFake   ->GetBinContent(6) ;
     log<<boost::format("h_JFake     %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_JFake   ->Integral(1,-1) %h_JFake   ->GetBinContent(1) %h_JFake   ->GetBinContent(2) %h_JFake   ->GetBinContent(3) %h_JFake   ->GetBinContent(4) %h_JFake   ->GetBinContent(5) %h_JFake   ->GetBinContent(6) ;
     log<<boost::format("h_Halo      %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_Halo    ->Integral(1,-1) %h_Halo    ->GetBinContent(1) %h_Halo    ->GetBinContent(2) %h_Halo    ->GetBinContent(3) %h_Halo    ->GetBinContent(4) %h_Halo    ->GetBinContent(5) %h_Halo    ->GetBinContent(6) ;
     log<<boost::format("h_Spike     %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_Spike   ->Integral(1,-1) %h_Spike   ->GetBinContent(1) %h_Spike   ->GetBinContent(2) %h_Spike   ->GetBinContent(3) %h_Spike   ->GetBinContent(4) %h_Spike   ->GetBinContent(5) %h_Spike   ->GetBinContent(6) ;
     log<<boost::format("h_GJets     %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_GJets   ->Integral(1,-1) %h_GJets   ->GetBinContent(1) %h_GJets   ->GetBinContent(2) %h_GJets   ->GetBinContent(3) %h_GJets   ->GetBinContent(4) %h_GJets   ->GetBinContent(5) %h_GJets   ->GetBinContent(6) ;
     log<<boost::format("h_TTGJets   %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_TTGJets ->Integral(1,-1) %h_TTGJets ->GetBinContent(1) %h_TTGJets ->GetBinContent(2) %h_TTGJets ->GetBinContent(3) %h_TTGJets ->GetBinContent(4) %h_TTGJets ->GetBinContent(5) %h_TTGJets ->GetBinContent(6) ;
     log<<boost::format("h_TGJets    %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_TGJets  ->Integral(1,-1) %h_TGJets  ->GetBinContent(1) %h_TGJets  ->GetBinContent(2) %h_TGJets  ->GetBinContent(3) %h_TGJets  ->GetBinContent(4) %h_TGJets  ->GetBinContent(5) %h_TGJets  ->GetBinContent(6) ;
     log<<boost::format("h_ZllGJets  %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_ZllGJets->Integral(1,-1) %h_ZllGJets->GetBinContent(1) %h_ZllGJets->GetBinContent(2) %h_ZllGJets->GetBinContent(3) %h_ZllGJets->GetBinContent(4) %h_ZllGJets->GetBinContent(5) %h_ZllGJets->GetBinContent(6) ;
     log<<boost::format("h_ZllJets   %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_ZllJets ->Integral(1,-1) %h_ZllJets ->GetBinContent(1) %h_ZllJets ->GetBinContent(2) %h_ZllJets ->GetBinContent(3) %h_ZllJets ->GetBinContent(4) %h_ZllJets ->GetBinContent(5) %h_ZllJets ->GetBinContent(6) ;
     log<<boost::format("h_GGJets    %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_GGJets  ->Integral(1,-1) %h_GGJets  ->GetBinContent(1) %h_GGJets  ->GetBinContent(2) %h_GGJets  ->GetBinContent(3) %h_GGJets  ->GetBinContent(4) %h_GGJets  ->GetBinContent(5) %h_GGJets  ->GetBinContent(6) ;
     log<<boost::format("h_WWG       %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_WWG     ->Integral(1,-1) %h_WWG     ->GetBinContent(1) %h_WWG     ->GetBinContent(2) %h_WWG     ->GetBinContent(3) %h_WWG     ->GetBinContent(4) %h_WWG     ->GetBinContent(5) %h_WWG     ->GetBinContent(6) ;
     log<<boost::format("h_WZ        %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_WZ      ->Integral(1,-1) %h_WZ      ->GetBinContent(1) %h_WZ      ->GetBinContent(2) %h_WZ      ->GetBinContent(3) %h_WZ      ->GetBinContent(4) %h_WZ      ->GetBinContent(5) %h_WZ      ->GetBinContent(6) ;
     log<<boost::format("h_WlnGJets  %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_WlnGJets->Integral(1,-1) %h_WlnGJets->GetBinContent(1) %h_WlnGJets->GetBinContent(2) %h_WlnGJets->GetBinContent(3) %h_WlnGJets->GetBinContent(4) %h_WlnGJets->GetBinContent(5) %h_WlnGJets->GetBinContent(6) ;
     log<<boost::format("h_Wmn       %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_Wmn     ->Integral(1,-1) %h_Wmn     ->GetBinContent(1) %h_Wmn     ->GetBinContent(2) %h_Wmn     ->GetBinContent(3) %h_Wmn     ->GetBinContent(4) %h_Wmn     ->GetBinContent(5) %h_Wmn     ->GetBinContent(6) ;
     log<<boost::format("h_Wtn       %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_Wtn     ->Integral(1,-1) %h_Wtn     ->GetBinContent(1) %h_Wtn     ->GetBinContent(2) %h_Wtn     ->GetBinContent(3) %h_Wtn     ->GetBinContent(4) %h_Wtn     ->GetBinContent(5) %h_Wtn     ->GetBinContent(6) ;
     log<<boost::format("h_ZZ        %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_ZZ      ->Integral(1,-1) %h_ZZ      ->GetBinContent(1) %h_ZZ      ->GetBinContent(2) %h_ZZ      ->GetBinContent(3) %h_ZZ      ->GetBinContent(4) %h_ZZ      ->GetBinContent(5) %h_ZZ      ->GetBinContent(6) ;
     log<<boost::format("h_ZnnGJets  %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n") %h_ZnnGJets->Integral(1,-1) %h_ZnnGJets->GetBinContent(1) %h_ZnnGJets->GetBinContent(2) %h_ZnnGJets->GetBinContent(3) %h_ZnnGJets->GetBinContent(4) %h_ZnnGJets->GetBinContent(5) %h_ZnnGJets->GetBinContent(6) ;
 
     if(dobygev){
      if(variablename=="et"      ||
         variablename=="uncorret"||
         variablename=="leptoMET"
        ){
       Double_t xbins[6] = {175,190,250,400,700,1000};
        TH1 *hr_data_obs = data_obs  ->Rebin(5,"hr_data_obs",xbins);    
        TH1 *hr_GJets    = h_GJets   ->Rebin(5,"hr_GJets",   xbins);    
        TH1 *hr_TTGJets  = h_TTGJets ->Rebin(5,"hr_TTGJets", xbins);    
        TH1 *hr_TGJets   = h_TGJets  ->Rebin(5,"hr_TGJets",  xbins);    
        TH1 *hr_WlnGJets = h_WlnGJets->Rebin(5,"hr_WlnGJets",xbins);    
        TH1 *hr_ZllGJets = h_ZllGJets->Rebin(5,"hr_ZllGJets",xbins);    
        TH1 *hr_ZllJets  = h_ZllJets ->Rebin(5,"hr_ZllJets", xbins);    
        TH1 *hr_GGJets   = h_GGJets  ->Rebin(5,"hr_GGJets",  xbins);
        TH1 *hr_EFake    = h_EFake   ->Rebin(5,"hr_EFake",   xbins);
        TH1 *hr_JFake    = h_JFake   ->Rebin(5,"hr_JFake",   xbins);
        TH1 *hr_WWG      = h_WWG     ->Rebin(5,"hr_WWG",     xbins);
        TH1 *hr_WZ       = h_WZ      ->Rebin(5,"hr_WZ",      xbins);
        TH1 *hr_Wmn      = h_Wmn     ->Rebin(5,"hr_Wmn",     xbins);
        TH1 *hr_Wtn      = h_Wtn     ->Rebin(5,"hr_Wtn",     xbins);
        TH1 *hr_ZZ       = h_ZZ      ->Rebin(5,"hr_ZZ",      xbins);
        TH1 *hr_ZnnGJets = h_ZnnGJets->Rebin(5,"hr_ZnnGJets",xbins);
  //std::cout<<hr_data_obs->GetSize()<<std::endl;
       for(int i=0; i<5; i++){
        int bnr = i+1;
        double divby = xbins[i+1] - xbins[i];
      
        double cont_data_obs = data_obs  ->GetBinContent(bnr); 
        double err_data_obs  = data_obs  ->GetBinError(bnr); 
        double err_ratio = err_data_obs/cont_data_obs;
    
        double cont_GJets    = h_GJets    ->GetBinContent(bnr); 
        double cont_TTGJets  = h_TTGJets  ->GetBinContent(bnr); 
        double cont_TGJets   = h_TGJets   ->GetBinContent(bnr); 
        double cont_WlnGJets = h_WlnGJets ->GetBinContent(bnr); 
        double cont_ZllGJets = h_ZllGJets ->GetBinContent(bnr); 
        double cont_ZllJets  = h_ZllJets  ->GetBinContent(bnr); 
        double cont_GGJets   = h_GGJets   ->GetBinContent(bnr);
        double cont_EFake    = h_EFake    ->GetBinContent(bnr);
        double cont_JFake    = h_JFake    ->GetBinContent(bnr);
        double cont_Halo     = h_Halo     ->GetBinContent(bnr);
        double cont_Spike    = h_Spike    ->GetBinContent(bnr);
        double cont_WWG      = h_WWG      ->GetBinContent(bnr);
        double cont_WZ       = h_WZ       ->GetBinContent(bnr);
        double cont_Wmn      = h_Wmn      ->GetBinContent(bnr);
        double cont_Wtn      = h_Wtn      ->GetBinContent(bnr);
        double cont_ZZ       = h_ZZ       ->GetBinContent(bnr);
        double cont_ZnnGJets = h_ZnnGJets ->GetBinContent(bnr);
        double cont_GWZtt    = h_GWZtt    ->GetBinContent(bnr);
        double cont_VV       = h_VV       ->GetBinContent(bnr);
      
        data_obs->SetBinContent(bnr, cont_data_obs / divby );    
        h_GJets   ->SetBinContent(bnr, cont_GJets    / divby );    
        h_TTGJets ->SetBinContent(bnr, cont_TTGJets  / divby );    
        h_TGJets  ->SetBinContent(bnr, cont_TGJets   / divby );    
        h_WlnGJets->SetBinContent(bnr, cont_WlnGJets / divby );    
        h_ZllGJets->SetBinContent(bnr, cont_ZllGJets / divby );    
        h_ZllJets ->SetBinContent(bnr, cont_ZllJets  / divby );    
        h_GGJets  ->SetBinContent(bnr, cont_GGJets   / divby );
        h_EFake   ->SetBinContent(bnr, cont_EFake    / divby );
        h_JFake   ->SetBinContent(bnr, cont_JFake    / divby );
        h_Halo    ->SetBinContent(bnr, cont_Halo     / divby );
        h_Spike   ->SetBinContent(bnr, cont_Spike    / divby );
        h_WWG     ->SetBinContent(bnr, cont_WWG      / divby );
        h_WZ      ->SetBinContent(bnr, cont_WZ       / divby );
        h_Wmn     ->SetBinContent(bnr, cont_Wmn      / divby );
        h_Wtn     ->SetBinContent(bnr, cont_Wtn      / divby );
        h_ZZ      ->SetBinContent(bnr, cont_ZZ       / divby );
        h_ZnnGJets->SetBinContent(bnr, cont_ZnnGJets / divby );
        h_GWZtt   ->SetBinContent(bnr, cont_GWZtt    / divby );
        h_VV      ->SetBinContent(bnr, cont_VV       / divby );
      
        data_obs ->SetBinError(bnr, cont_data_obs*err_ratio / divby );    
      
       }
      }
     } // dobygev
 
     THStack *thestack = new THStack("thestack","");
     thestack->Add(h_GWZtt   );  
     thestack->Add(h_VV      );  
     thestack->Add(h_Halo    );
     thestack->Add(h_Spike   );
     thestack->Add(h_JFake   );
     thestack->Add(h_EFake   );
     thestack->Add(h_WlnGJets);
     thestack->Add(h_ZnnGJets);
     thestack->Add(h_ZllGJets);

     TLegend *leg = new TLegend(0.20,0.7,0.88,0.85 );
     //TLegend *leg = new TLegend(0.60,0.5,0.88,0.88 );
     leg-> SetNColumns(2);
     leg->SetBorderSize(0);
     leg->SetFillColor(kWhite);
     leg->AddEntry( h_GWZtt   ,"#gamma#gamma+jets,#gamma+jet,W(#mu#nu,#tau#nu),t#bar{t}#gamma,t#gamma" , "f");
     leg->AddEntry( h_VV      ,"WZ,ZZ,WW#gamma", "f");
     leg->AddEntry( h_Halo    ,"Beam halo", "f");
     leg->AddEntry( h_Spike   ,"Spikes", "f");
     leg->AddEntry( h_JFake   ,"Jet#rightarrow#gamma misID", "f");  
     leg->AddEntry( h_EFake   ,"e#rightarrow#gamma misID", "f");  
     leg->AddEntry( h_WlnGJets,"W(l#nu)#gamma", "f");  
     leg->AddEntry( h_ZnnGJets,"Z(#nu#nu)#gamma", "f");  
     leg->AddEntry( h_ZllGJets,"Z(ll)#gamma", "f");  
     leg->AddEntry( data_obs  ,"data", "lep");
 
     canvas->cd();
 
       // data_obs carries title

     data_obs->SetXTitle(TString(data_obs->GetTitle())+" [GeV]");
     data_obs->SetTitle("");
     data_obs->SetYTitle("Events / Bin");
     data_obs->SetMaximum(1.3*data_obs->GetMaximum());
     if(dobygev){
      data_obs->SetYTitle("Events / GeV");
      //hr_data_obs->SetMaximum(10);
     }
     if(dolog){
      data_obs->SetMaximum(200);
      data_obs->SetMinimum(0.0005);
     }
     data_obs->GetYaxis()->SetTitleOffset(1.4);

     thestack->Draw("hist,sames");
 
//     h_ZllGJets->Draw("hist");
     data_obs->Draw("");
     thestack->Draw("hist,sames");
     data_obs->Draw("sames");
//     data_obs->GetXaxis()->Draw();
     leg->Draw("same");
 
       // h_ZnnGJets carries title

//     h_ZnnGJets->SetXTitle(TString(h_ZnnGJets->GetTitle())+" [GeV]");
//     h_ZnnGJets->SetTitle("");
//     h_ZnnGJets->SetYTitle("Events / Bin");
//     //hr_h_ZnnGJets->SetMaximum(100);
//     if(dobygev){
//      h_ZnnGJets->SetYTitle("Events / GeV");
//      //hr_h_ZnnGJets->SetMaximum(10);
//     }
//     h_ZnnGJets->SetMaximum(130);
//     h_ZnnGJets->SetMinimum(0);
//     if(dolog){
//      h_ZnnGJets->SetMaximum(10000);
//      h_ZnnGJets->SetMinimum(0.0005);
//     }
//     h_ZnnGJets->GetYaxis()->SetTitleOffset(1.4);
// 
//     h_ZnnGJets->Draw("");
//     thestack->Draw("hist,sames");
//     h_ZnnGJets->GetXaxis()->Draw();
//     leg->Draw("same");
 
     //title->DrawTextNDC(0.17,0.87,"CMS");
     //extra->DrawTextNDC(0.17,0.81,"Preliminary");
     title->DrawTextNDC(0.13,0.91,"CMS");
     extra->DrawTextNDC(0.23,0.91,"Preliminary");
     lumi->DrawTextNDC(0.9,0.91,"12.9 /fb (13 TeV)");

     gPad->Update();
     gPad->RedrawAxis();
 
     data_obs  ->Write();
     h_GWZtt   ->Write();  
     h_VV      ->Write();  
     h_Halo    ->Write();
     h_Spike   ->Write();
     h_JFake   ->Write();
     h_EFake   ->Write();
     h_WlnGJets->Write();
     h_ZnnGJets->Write();
     h_ZllGJets->Write();

     canvas->SaveAs(outpath+"/"+outname+".pdf");
 
     thestack  ->Delete();
     raw_data_obs   ->Delete();
     data_obs       ->Delete();

     h_raw_GJets    ->Delete();   
     h_raw_TTGJets  ->Delete();   
     h_raw_TGJets   ->Delete();   
     h_raw_ZllGJets ->Delete();   
     h_raw_ZllJets  ->Delete();   
     h_raw_GGJets   ->Delete();   
     h_raw_EFake    ->Delete();   
     h_raw_JFake    ->Delete();   
     h_raw_Halo     ->Delete();   
     h_raw_Spike    ->Delete();   
     h_raw_WWG      ->Delete();    
     h_raw_WZ       ->Delete();    
     h_raw_WlnGJets ->Delete();
     h_raw_Wmn      ->Delete();    
     h_raw_Wtn      ->Delete();    
     h_raw_ZZ       ->Delete();    
     h_raw_ZnnGJets ->Delete();
     h_raw_GWZtt    ->Delete();   
     h_raw_VV       ->Delete();   

     h_GJets    ->Delete();   
     h_TTGJets  ->Delete();   
     h_TGJets   ->Delete();   
     h_ZllGJets ->Delete();   
     h_ZllJets  ->Delete();   
     h_GGJets   ->Delete();   
     h_EFake    ->Delete();   
     h_JFake    ->Delete();   
     h_Halo     ->Delete();   
     h_Spike    ->Delete();   
     h_WWG      ->Delete();         
     h_WZ       ->Delete();         
     h_WlnGJets ->Delete();
     h_Wmn      ->Delete();         
     h_Wtn      ->Delete();         
     h_ZZ       ->Delete();         
     h_ZnnGJets ->Delete();
     h_GWZtt    ->Delete();   
     h_VV       ->Delete();
 
     log.close();
 
     } //  for(unsigned int s=0; s<selectionnames.size(); ++s)
    file_out->Close();
    } //  for(unsigned int p=0; p<ptrangenames.size(); ++p)
   } //  for(unsigned int v=0; v<variablenames.size(); ++v)
 }
}
