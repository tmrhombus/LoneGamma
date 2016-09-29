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
 wwwpath = TString("/afs/hep.wisc.edu/user/tperry/www/MonoPhoton/qcdPlots/"+Tversion);

 name_SinglePhotonData  =  "analyzed_SinglePhotonData.root";
 name_GJets             =  "analyzed_GJets.root";
 name_TTGJets           =  "analyzed_TTGJets.root";
 name_TGJets            =  "analyzed_TGJets.root";
 name_ZllGJets          =  "analyzed_ZllGJets.root";
 name_ZllJets           =  "analyzed_ZllJets.root";
 name_GGJets            =  "analyzed_GGJets.root";
 name_SinglePhotonEle   =  "analyzed_SinglePhotonEle.root";
 name_SinglePhotonJet   =  "analyzed_SinglePhotonJet.root";
 name_WWG               =  "analyzed_WWG.root";
 name_WZ                =  "analyzed_WZ.root";
 name_WlnGJets          =  "analyzed_WlnGJets.root";
 name_Wmn               =  "analyzed_Wmn.root";
 name_Wtn               =  "analyzed_Wtn.root";
 name_ZZ                =  "analyzed_ZZ.root";
 name_ZnnGJets          =  "analyzed_ZnnGJets.root";
 
 file_SinglePhotonData  =  new TFile(inpath+"/"+name_SinglePhotonData,"READ");
 file_GJets         =  new TFile(inpath+"/"+name_GJets,"READ");
 file_TTGJets       =  new TFile(inpath+"/"+name_TTGJets,"READ");
 file_TGJets        =  new TFile(inpath+"/"+name_TGJets,"READ");
 file_ZllGJets      =  new TFile(inpath+"/"+name_ZllGJets,"READ");
 file_ZllJets       =  new TFile(inpath+"/"+name_ZllJets,"READ");
 file_GGJets            = new TFile(inpath+"/"+name_GGJets,"READ");
 file_SinglePhotonEle   = new TFile(inpath+"/"+name_SinglePhotonEle,"READ");
 file_SinglePhotonJet   = new TFile(inpath+"/"+name_SinglePhotonJet,"READ");
 file_WWG               = new TFile(inpath+"/"+name_WWG,"READ");
 file_WZ                = new TFile(inpath+"/"+name_WZ,"READ");
 file_WlnGJets          = new TFile(inpath+"/"+name_WlnGJets,"READ");
 file_Wmn               = new TFile(inpath+"/"+name_Wmn,"READ");
 file_Wtn               = new TFile(inpath+"/"+name_Wtn,"READ");
 file_ZZ                = new TFile(inpath+"/"+name_ZZ,"READ");
 file_ZnnGJets          = new TFile(inpath+"/"+name_ZnnGJets,"READ");

 variablenames.push_back("uncorret"  );
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

// rebins.push_back(10);
// rebins.push_back(10);
// rebins.push_back(5);
// rebins.push_back(1);
// rebins.push_back(5);
// rebins.push_back(5);
// rebins.push_back(1);
// rebins.push_back(1);
// rebins.push_back(1);

// selectionnames.push_back("m110");
 selectionnames.push_back("m170");

 ptrangenames.push_back("allpt");
 //ptrangenames.push_back("175to1000");

 dolog = true;
 //dolog = false;

 //dobygev = true;
 dobygev = false;

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
 title->SetTextSize(0.07);
 title->SetTextColor(kBlack);
 title->SetTextAlign(13);
 title->SetTextFont(62);

 TText* extra = new TText(1,1,"") ;
 extra->SetTextSize(0.05);
 extra->SetTextColor(kBlack);
 extra->SetTextAlign(13);
 extra->SetTextFont(52);

 TText* lumi = new TText(1,1,"") ;
 lumi->SetTextSize(0.05);
 lumi->SetTextColor(kBlack);
 lumi->SetTextAlign(31);
 lumi->SetTextFont(42);

//  Int_t ci1  = 1001;
//  Int_t ci2  = 1002;
//  Int_t ci3  = 1003;
//  Int_t ci4  = 1004;
//  Int_t ci5  = 1005;
//  Int_t ci6  = 1006;
//  Int_t ci7  = 1007;
//  Int_t ci8  = 1008;
//  Int_t ci9  = 1009;
//  Int_t ci10 = 1010;
//  TColor *color1  = new TColor( ci1,178./255,223./255,138./255);
//  TColor *color2  = new TColor( ci2, 51./255,160./255, 44./255);
//  TColor *color3  = new TColor( ci3,251./255,154./255,153./255);
//  TColor *color4  = new TColor( ci4,227./255, 26./255, 28./255);
//  TColor *color5  = new TColor( ci5,166./255,206./255,227./255);
//  TColor *color6  = new TColor( ci6, 31./255,120./255,180./255);
//  TColor *color7  = new TColor( ci7,177./255, 89./255, 40./255);
//  TColor *color8  = new TColor( ci8,236./255,208./255,120./255);
//  TColor *color9  = new TColor( ci9,253./255,191./255,111./255);
//  TColor *color10 = new TColor(ci10,255./255,127./255, 0./255);

  for(unsigned int v=0; v<variablenames.size(); ++v){
  TString variablename = variablenames[v];
  Int_t rebin = rebins[v];
   for(unsigned int s=0; s<selectionnames.size(); ++s){
   TString selectionname = selectionnames[s];
    for(unsigned int p=0; p<ptrangenames.size(); ++p){
    TString ptrangename = ptrangenames[p];

    canvas->Clear();

    TString tailname = variablename+"_"+ptrangename; //+"_"+selectionname;
    outname = "plot_"+tailname+extraname; 
    histname = "h_sig_"+tailname;
    //histname = "h_mu_"+tailname;
    file_out =  new TFile(outpath+"/"+outname+".root","RECREATE");

    printf(outname+"\n");
    printf(histname+"\n");

    data_obs = (TH1F*)file_SinglePhotonData->Get(histname)->Clone("data_obs");

    h_GJets    = (TH1F*)file_GJets           ->Get(histname)->Clone("h_GJets")    ;   
    h_TTGJets  = (TH1F*)file_TTGJets         ->Get(histname)->Clone("h_TTGJets")  ;   
    h_TGJets   = (TH1F*)file_TGJets          ->Get(histname)->Clone("h_TGJets")   ;   
    h_ZllGJets = (TH1F*)file_ZllGJets        ->Get(histname)->Clone("h_ZllGJets") ;   
    h_ZllJets  = (TH1F*)file_ZllJets         ->Get(histname)->Clone("h_ZllJets")  ;  
    h_GGJets   = (TH1F*)file_GGJets          ->Get(histname)->Clone("h_GGJets")   ;  
    h_EFake    = (TH1F*)file_SinglePhotonEle ->Get(histname)->Clone("h_EFake")    ;   
    h_JFake    = (TH1F*)file_SinglePhotonJet ->Get(histname)->Clone("h_JFake")    ;   
    h_WWG      = (TH1F*)file_WWG             ->Get(histname)->Clone("h_WWG")      ;     
    h_WZ       = (TH1F*)file_WZ              ->Get(histname)->Clone("h_WZ")       ;      
    h_WlnGJets = (TH1F*)file_WlnGJets        ->Get(histname)->Clone("h_WlnGJets") ;
    h_Wmn      = (TH1F*)file_Wmn             ->Get(histname)->Clone("h_Wmn")      ;     
    h_Wtn      = (TH1F*)file_Wtn             ->Get(histname)->Clone("h_Wtn")      ;     
    h_ZZ       = (TH1F*)file_ZZ              ->Get(histname)->Clone("h_ZZ")       ;      
    h_ZnnGJets = (TH1F*)file_ZnnGJets        ->Get(histname)->Clone("h_ZnnGJets") ;

    data_obs  ->SetMarkerStyle( 20 );  
    data_obs  ->SetMarkerSize( 1 );  
    data_obs  ->SetLineWidth( 2 );

    data_obs  ->SetLineColor( kBlack );
    h_GJets   ->SetLineColor( kBlack ); 
    h_TTGJets ->SetLineColor( kBlack ); 
    h_TGJets  ->SetLineColor( kBlack ); 
    h_ZllGJets->SetLineColor( kBlack ); 
    h_ZllJets ->SetLineColor( kBlack ); 
    h_GGJets  ->SetLineColor( kBlack );
    h_EFake   ->SetLineColor( kBlack );
    h_JFake   ->SetLineColor( kBlack );
    h_WWG     ->SetLineColor( kBlack );
    h_WZ      ->SetLineColor( kBlack );
    h_WlnGJets->SetLineColor( kBlack );
    h_Wmn     ->SetLineColor( kBlack );
    h_Wtn     ->SetLineColor( kBlack );
    h_ZZ      ->SetLineColor( kBlack );
    h_ZnnGJets->SetLineColor( kBlack );

    data_obs  ->Rebin( rebin ); 
    h_GJets   ->Rebin( rebin ); 
    h_TTGJets ->Rebin( rebin ); 
    h_TGJets  ->Rebin( rebin ); 
    h_ZllGJets->Rebin( rebin ); 
    h_ZllJets ->Rebin( rebin ); 
    h_GGJets  ->Rebin( rebin );
    h_EFake   ->Rebin( rebin );
    h_JFake   ->Rebin( rebin );
    h_WWG     ->Rebin( rebin );
    h_WZ      ->Rebin( rebin );
    h_WlnGJets->Rebin( rebin );
    h_Wmn     ->Rebin( rebin );
    h_Wtn     ->Rebin( rebin );
    h_ZZ      ->Rebin( rebin );
    h_ZnnGJets->Rebin( rebin );

    h_GJets   ->SetLineWidth( 2 ); 
    h_TTGJets ->SetLineWidth( 2 ); 
    h_TGJets  ->SetLineWidth( 2 ); 
    h_ZllGJets->SetLineWidth( 2 ); 
    h_ZllJets ->SetLineWidth( 2 ); 
    h_GGJets  ->SetLineWidth( 2 );
    h_EFake   ->SetLineWidth( 2 );
    h_JFake   ->SetLineWidth( 2 );
    h_WWG     ->SetLineWidth( 2 );
    h_WZ      ->SetLineWidth( 2 );
    h_WlnGJets->SetLineWidth( 2 );
    h_Wmn     ->SetLineWidth( 2 );
    h_Wtn     ->SetLineWidth( 2 );
    h_ZZ      ->SetLineWidth( 2 );
    h_ZnnGJets->SetLineWidth( 2 );

  //  leg->AddEntry( h_GJets   ,"#gamma+jets", "f");
  //  leg->AddEntry( h_TTGJets ,"t#bar{t}#gamma+jets", "f");
  //  leg->AddEntry( h_TGJets  ,"t#gamma+jets", "f");
  //  leg->AddEntry( h_WlnGJets,"W(l#nu)#gamma+jets", "f");
  //  leg->AddEntry( h_ZllGJets,"Z(ll)#gamma+jets", "f");
  //  leg->AddEntry( h_ZllJets ,"Z(ll)+jets", "f");
  //  //leg->AddEntry( h_GGJets  , "#gamma#gamma+jets","f" );
  //  leg->AddEntry( h_EFake   , "e#rightarrow fake","f" );
  //  //leg->AddEntry( h_JFake   , "jet#rightarrow fake","f" );
  //  leg->AddEntry( h_WWG     , "WW#gamma","f" );
  //  //leg->AddEntry( h_WZ      , "WZ","f" );
  //  //leg->AddEntry( h_Wmn     , "W(#mu#nu)","f" );
  //  //leg->AddEntry( h_Wtn     , "W(#tau#nu)","f" );
  //  leg->AddEntry( h_ZZ      , "ZZ","f" );
  //  //leg->AddEntry( h_ZnnGJets, "Z(#nu#nu)#gamma+jets","f" );

    h_GJets   ->SetFillColor( kRed-10   );  
    h_TTGJets ->SetFillColor( kAzure-9  ); 
    h_TGJets  ->SetFillColor( kBlue-9   ); 
    h_WlnGJets->SetFillColor( kRed-6    );  
    h_ZllGJets->SetFillColor( kOrange-4 );
    h_ZllJets ->SetFillColor( kBlue-8   );  
    h_GGJets  ->SetFillColor( kRed+1    );
    h_EFake   ->SetFillColor( kOrange-3 );
    h_JFake   ->SetFillColor( kYellow-3 );
    h_WWG     ->SetFillColor( kGreen+1  );
    h_WZ      ->SetFillColor( kGreen-5  );
    h_WlnGJets->SetFillColor( kGreen-9  );
    h_Wmn     ->SetFillColor( kAzure+10 );
    h_Wtn     ->SetFillColor( kBlue+0   );
    h_ZZ      ->SetFillColor( kBlue-9   );
    h_ZnnGJets->SetFillColor( 51        );

    //if(lepton=="mu"){printf("lepton mu\n"  ); h_EFake->Scale( 0.5/h_EFake->Integral() );    }
    //if(lepton=="ele"){printf("lepton ele\n"); h_EFake->Scale( 0.1225/h_EFake->Integral() ); }
    h_JFake->Scale( 5.9/h_JFake->Integral(-1,-1) );
    h_EFake->Scale( 52.7/h_EFake->Integral(-1,-1) );

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

    if(dobygev){
     if(variablename=="et"      ||
        variablename=="uncorret"||
        variablename=="leptoMET"
       ){
    //  Double_t xbins[6] = {175,190,250,400,700,1000};
    //   TH1 *hr_data_obs = data_obs  ->Rebin(5,"hr_data_obs",xbins);    
    //   TH1 *hr_GJets    = h_GJets   ->Rebin(5,"hr_GJets",   xbins);    
    //   TH1 *hr_TTGJets  = h_TTGJets ->Rebin(5,"hr_TTGJets", xbins);    
    //   TH1 *hr_TGJets   = h_TGJets  ->Rebin(5,"hr_TGJets",  xbins);    
    //   TH1 *hr_WlnGJets = h_WlnGJets->Rebin(5,"hr_WlnGJets",xbins);    
    //   TH1 *hr_ZllGJets = h_ZllGJets->Rebin(5,"hr_ZllGJets",xbins);    
    //   TH1 *hr_ZllJets  = h_ZllJets ->Rebin(5,"hr_ZllJets", xbins);    
    //   TH1 *hr_GGJets   = h_GGJets  ->Rebin(5,"hr_GGJets",  xbins);
    //   TH1 *hr_EFake    = h_EFake   ->Rebin(5,"hr_EFake",   xbins);
    //   TH1 *hr_JFake    = h_JFake   ->Rebin(5,"hr_JFake",   xbins);
    //   TH1 *hr_WWG      = h_WWG     ->Rebin(5,"hr_WWG",     xbins);
    //   TH1 *hr_WZ       = h_WZ      ->Rebin(5,"hr_WZ",      xbins);
    //   TH1 *hr_Wmn      = h_Wmn     ->Rebin(5,"hr_Wmn",     xbins);
    //   TH1 *hr_Wtn      = h_Wtn     ->Rebin(5,"hr_Wtn",     xbins);
    //   TH1 *hr_ZZ       = h_ZZ      ->Rebin(5,"hr_ZZ",      xbins);
    //   TH1 *hr_ZnnGJets = h_ZnnGJets->Rebin(5,"hr_ZnnGJets",xbins);
 //std::cout<<hr_data_obs->GetSize()<<std::endl;
      for(int i=0; i<5; i++){
       int bnr = i+1;
       double divby = xbins[i+1] - xbins[i];
     
       double cont_data_obs = hr_data_obs  ->GetBinContent(bnr); 
       double err_data_obs  = hr_data_obs  ->GetBinError(bnr); 
       double err_ratio = err_data_obs/cont_data_obs;
   
       double cont_GJets    = hr_GJets   ->GetBinContent(bnr); 
       double cont_TTGJets  = hr_TTGJets ->GetBinContent(bnr); 
       double cont_TGJets   = hr_TGJets  ->GetBinContent(bnr); 
       double cont_WlnGJets = hr_WlnGJets->GetBinContent(bnr); 
       double cont_ZllGJets = hr_ZllGJets->GetBinContent(bnr); 
       double cont_ZllJets  = hr_ZllJets ->GetBinContent(bnr); 
       double cont_GGJets   = hr_GGJets   ->GetBinContent(bnr);
       double cont_EFake    = hr_EFake    ->GetBinContent(bnr);
       double cont_JFake    = hr_JFake    ->GetBinContent(bnr);
       double cont_WWG      = hr_WWG      ->GetBinContent(bnr);
       double cont_WZ       = hr_WZ       ->GetBinContent(bnr);
       double cont_Wmn      = hr_Wmn      ->GetBinContent(bnr);
       double cont_Wtn      = hr_Wtn      ->GetBinContent(bnr);
       double cont_ZZ       = hr_ZZ       ->GetBinContent(bnr);
       double cont_ZnnGJets = hr_ZnnGJets ->GetBinContent(bnr);
     
       hr_data_obs->SetBinContent(bnr, cont_data_obs / divby );    
       hr_GJets   ->SetBinContent(bnr, cont_GJets    / divby );    
       hr_TTGJets ->SetBinContent(bnr, cont_TTGJets  / divby );    
       hr_TGJets  ->SetBinContent(bnr, cont_TGJets   / divby );    
       hr_WlnGJets->SetBinContent(bnr, cont_WlnGJets / divby );    
       hr_ZllGJets->SetBinContent(bnr, cont_ZllGJets / divby );    
       hr_ZllJets ->SetBinContent(bnr, cont_ZllJets  / divby );    
       hr_GGJets  ->SetBinContent(bnr, cont_GGJets   / divby );
       hr_EFake   ->SetBinContent(bnr, cont_EFake    / divby );
       hr_JFake   ->SetBinContent(bnr, cont_JFake    / divby );
       hr_WWG     ->SetBinContent(bnr, cont_WWG      / divby );
       hr_WZ      ->SetBinContent(bnr, cont_WZ       / divby );
       hr_Wmn     ->SetBinContent(bnr, cont_Wmn      / divby );
       hr_Wtn     ->SetBinContent(bnr, cont_Wtn      / divby );
       hr_ZZ      ->SetBinContent(bnr, cont_ZZ       / divby );
       hr_ZnnGJets->SetBinContent(bnr, cont_ZnnGJets / divby );
     
       hr_data_obs ->SetBinError(bnr, cont_data_obs*err_ratio / divby );    
     
      }
     }
    } // dobygev

    THStack *thestack = new THStack("thestack","");
    thestack->Add( hr_GJets    );
    thestack->Add( hr_TTGJets  );
    //thestack->Add( hr_TGJets  );
    thestack->Add( hr_WlnGJets );
    thestack->Add( hr_ZllGJets );
    thestack->Add( hr_ZllJets  );
    thestack->Add( hr_GGJets   );
    thestack->Add( hr_EFake    );
    thestack->Add( hr_JFake    );
    thestack->Add( hr_WWG      );
    thestack->Add( hr_WZ       );
    thestack->Add( hr_Wmn      );
    thestack->Add( hr_Wtn      );
    thestack->Add( hr_ZZ       );
    thestack->Add( hr_ZnnGJets );

    TLegend *leg = new TLegend(0.60,0.5,0.88,0.88 );
    leg->SetFillColor(kWhite);
    leg->AddEntry( hr_data_obs  ,"data", "lep");
    leg->AddEntry( hr_GJets   ,"#gamma+jets", "f");
    leg->AddEntry( hr_TTGJets ,"t#bar{t}#gamma+jets", "f");
    //leg->AddEntry( hr_TGJets  ,"t#gamma+jets", "f");
    leg->AddEntry( hr_WlnGJets,"W(l#nu)#gamma+jets", "f");
    leg->AddEntry( hr_ZllGJets,"Z(ll)#gamma+jets", "f");
    leg->AddEntry( hr_ZllJets ,"Z(ll)+jets", "f");
    leg->AddEntry( hr_GGJets  , "#gamma#gamma+jets","f" );
    leg->AddEntry( hr_EFake   , "e#rightarrow fake","f" );
    leg->AddEntry( hr_JFake   , "jet#rightarrow fake","f" );
    leg->AddEntry( hr_WWG     , "WW#gamma","f" );
    leg->AddEntry( hr_WZ      , "WZ","f" );
    leg->AddEntry( hr_Wmn     , "W(#mu#nu)","f" );
    leg->AddEntry( hr_Wtn     , "W(#tau#nu)","f" );
    leg->AddEntry( hr_ZZ      , "ZZ","f" );
    leg->AddEntry( hr_ZnnGJets, "Z(#nu#nu)#gamma+jets","f" );

    canvas->cd();

    hr_data_obs->SetXTitle(TString(data_obs->GetTitle())+" [GeV]");
    hr_data_obs->SetTitle("");
    hr_data_obs->SetYTitle("Events / Bin");
    hr_data_obs->SetMaximum(100);
    if(dobygev){
     hr_data_obs->SetYTitle("Events / GeV");
    hr_data_obs->SetMaximum(10);
    }
    hr_data_obs->GetYaxis()->SetTitleOffset(1.4);

    hr_data_obs->Draw("");
    thestack->Draw("hist,sames");
    hr_data_obs->Draw("sames");
    leg->Draw("same");

    title->DrawTextNDC(0.17,0.87,"CMS");
    extra->DrawTextNDC(0.17,0.81,"Preliminary");
    lumi->DrawTextNDC(0.9,0.91,"12.9 /fb (13 TeV)");

    hr_data_obs  ->Write();
    hr_GJets   ->Write();   
    hr_TGJets  ->Write();   
    hr_TTGJets ->Write();   
    hr_WlnGJets->Write();   
    hr_ZllGJets->Write();   
    hr_ZllJets ->Write();  
    hr_GGJets  ->Write();
    hr_EFake   ->Write();
    hr_JFake   ->Write();
    hr_WWG     ->Write();
    hr_WZ      ->Write();
    hr_Wmn     ->Write();
    hr_Wtn     ->Write();
    hr_ZZ      ->Write();
    hr_ZnnGJets->Write();

    canvas->SaveAs(outpath+"/"+outname+".pdf");

    thestack  ->Delete();
    data_obs  ->Delete();
    h_GJets   ->Delete();   
    h_TTGJets ->Delete();   
    h_TGJets  ->Delete();   
    h_WlnGJets->Delete();   
    h_ZllGJets->Delete();   
    h_ZllJets ->Delete();  
    h_GGJets  ->Delete();
    h_EFake   ->Delete();
    h_JFake   ->Delete();
    h_WWG     ->Delete();
    h_WZ      ->Delete();
    h_Wmn     ->Delete();
    h_Wtn     ->Delete();
    h_ZZ      ->Delete();
    h_ZnnGJets->Delete();

    hr_data_obs  ->Delete();
    hr_GJets   ->Delete();   
    hr_TTGJets ->Delete();   
    hr_TGJets  ->Delete();   
    hr_WlnGJets->Delete();   
    hr_ZllGJets->Delete();   
    hr_ZllJets ->Delete();  
    hr_GGJets  ->Delete();
    hr_EFake   ->Delete();
    hr_JFake   ->Delete();
    hr_WWG     ->Delete();
    hr_WZ      ->Delete();
    hr_Wmn     ->Delete();
    hr_Wtn     ->Delete();
    hr_ZZ      ->Delete();
    hr_ZnnGJets->Delete();

    file_out->Close();

    } //  for(unsigned int p=0; p<ptrangenames.size(); ++p)
   } //  for(unsigned int s=0; s<selectionnames.size(); ++s)
  } //  for(unsigned int v=0; v<variablenames.size(); ++v)

}
