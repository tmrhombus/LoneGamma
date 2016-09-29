#define plotwlng_cxx 
#include "plotwlng.h"

void plotwlng::Loop()
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
 //name_WZ                =  "analyzed_WZ.root";
 name_WlnGJets          =  "analyzed_WlnGJets.root";
 name_Wmn               =  "analyzed_Wmn.root";
 name_Wtn               =  "analyzed_Wtn.root";
 name_ZZ                =  "analyzed_ZZ.root";
 name_ZnnGJets          =  "analyzed_ZnnGJets.root";
 
 file_SinglePhotonData  =  new TFile(inpath+"/"+name_SinglePhotonData,"READ");
 file_GJets             =  new TFile(inpath+"/"+name_GJets,"READ");
 file_TTGJets           =  new TFile(inpath+"/"+name_TTGJets,"READ");
 file_TGJets            =  new TFile(inpath+"/"+name_TGJets,"READ");
 file_ZllGJets          =  new TFile(inpath+"/"+name_ZllGJets,"READ");
 file_ZllJets           =  new TFile(inpath+"/"+name_ZllJets,"READ");
 file_GGJets            = new TFile(inpath+"/"+name_GGJets,"READ");
 file_SinglePhotonEle   = new TFile(inpath+"/"+name_SinglePhotonEle,"READ");
 file_SinglePhotonJet   = new TFile(inpath+"/"+name_SinglePhotonJet,"READ");
 file_WWG               = new TFile(inpath+"/"+name_WWG,"READ");
 //file_WZ                = new TFile(inpath+"/"+name_WZ,"READ");
 file_WlnGJets          = new TFile(inpath+"/"+name_WlnGJets,"READ");
 file_Wmn               = new TFile(inpath+"/"+name_Wmn,"READ");
 file_Wtn               = new TFile(inpath+"/"+name_Wtn,"READ");
 file_ZZ                = new TFile(inpath+"/"+name_ZZ,"READ");
 file_ZnnGJets          = new TFile(inpath+"/"+name_ZnnGJets,"READ");

 leptons.push_back("mu");
 leptons.push_back("ele");

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

 ptrangenames.push_back("175to1000");

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

 for(unsigned int l=0; l<leptons.size(); ++l){
  TString lepton = leptons[l];
  for(unsigned int v=0; v<variablenames.size(); ++v){
  TString variablename = variablenames[v];
  Int_t rebin = rebins[v];
   for(unsigned int s=0; s<selectionnames.size(); ++s){
   TString selectionname = selectionnames[s];
    if(lepton=="mu"){ selectionname+="_m"; }
    if(lepton=="ele"){ selectionname+="_e"; }
    for(unsigned int p=0; p<ptrangenames.size(); ++p){
    TString ptrangename = ptrangenames[p];

    canvas->Clear();

    TString tailname = lepton+"_"+variablename+"_"+ptrangename+"_"+selectionname;
    outname = "plot_"+tailname+extraname; 
    histname = "h_"+tailname;
    //histname = "h_mu_"+tailname;
    file_out =  new TFile(outpath+"/"+outname+".root","RECREATE");

    printf(histname+"\n");
    printf(outname+"\n");

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
    //h_WZ       = (TH1F*)file_WZ              ->Get(histname)->Clone("h_WZ")       ;      
    h_WlnGJets = (TH1F*)file_WlnGJets        ->Get(histname)->Clone("h_WlnGJets") ;
    h_Wmn      = (TH1F*)file_Wmn             ->Get(histname)->Clone("h_Wmn")      ;     
    h_Wtn      = (TH1F*)file_Wtn             ->Get(histname)->Clone("h_Wtn")      ;     
    h_ZZ       = (TH1F*)file_ZZ              ->Get(histname)->Clone("h_ZZ")       ;      
    h_ZnnGJets = (TH1F*)file_ZnnGJets        ->Get(histname)->Clone("h_ZnnGJets") ;
    printf(" \n");

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
    //h_WZ      ->SetLineColor( kBlack );
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
    //h_WZ      ->Rebin( rebin );
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
    //h_WZ      ->SetLineWidth( 2 );
    h_WlnGJets->SetLineWidth( 2 );
    h_Wmn     ->SetLineWidth( 2 );
    h_Wtn     ->SetLineWidth( 2 );
    h_ZZ      ->SetLineWidth( 2 );
    h_ZnnGJets->SetLineWidth( 2 );

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
    //h_WZ      ->SetFillColor( kGreen-5  );
    h_WlnGJets->SetFillColor( kGreen-9  );
    h_Wmn     ->SetFillColor( kAzure+10 );
    h_Wtn     ->SetFillColor( kBlue+0   );
    h_ZZ      ->SetFillColor( kBlue-9   );
    h_ZnnGJets->SetFillColor( 51        );

    if(lepton=="mu"){printf("lepton mu\n"  ); h_EFake->Scale( 0.5/h_EFake->Integral() );    }
    if(lepton=="ele"){printf("lepton ele\n"); h_EFake->Scale( 0.1225/h_EFake->Integral() ); }
    //h_JFake   ->Scale( h_JFake->Integral() );

    if(dobygev){
     if(variablename=="et"      ||
        variablename=="uncorret"||
        variablename=="leptoMET"
       ){
      Double_t xbins[6] = {175,190,250,400,700,1000};
      for(int i=0; i<5; i++){
       int bnr = i+1;
       double divby = xbins[i+1] - xbins[i];
     
       double cont_data_obs = data_obs  ->GetBinContent(bnr); 
       double err_data_obs  = data_obs  ->GetBinError(bnr); 
       double err_ratio = err_data_obs/cont_data_obs;
   
       double cont_GJets    = h_GJets   ->GetBinContent(bnr); 
       double cont_TTGJets  = h_TTGJets ->GetBinContent(bnr); 
       double cont_TGJets   = h_TGJets  ->GetBinContent(bnr); 
       double cont_WlnGJets = h_WlnGJets->GetBinContent(bnr); 
       double cont_ZllGJets = h_ZllGJets->GetBinContent(bnr); 
       double cont_ZllJets  = h_ZllJets ->GetBinContent(bnr); 
       double cont_GGJets   = h_GGJets   ->GetBinContent(bnr);
       double cont_EFake    = h_EFake    ->GetBinContent(bnr);
       double cont_JFake    = h_JFake    ->GetBinContent(bnr);
       double cont_WWG      = h_WWG      ->GetBinContent(bnr);
       //double cont_WZ       = h_WZ       ->GetBinContent(bnr);
       double cont_Wmn      = h_Wmn      ->GetBinContent(bnr);
       double cont_Wtn      = h_Wtn      ->GetBinContent(bnr);
       double cont_ZZ       = h_ZZ       ->GetBinContent(bnr);
       double cont_ZnnGJets = h_ZnnGJets ->GetBinContent(bnr);
     
       data_obs  ->SetBinContent(bnr, cont_data_obs / divby );    
       h_GJets   ->SetBinContent(bnr, cont_GJets    / divby );    
       h_TTGJets ->SetBinContent(bnr, cont_TTGJets  / divby );    
       h_TGJets  ->SetBinContent(bnr, cont_TGJets   / divby );    
       h_WlnGJets->SetBinContent(bnr, cont_WlnGJets / divby );    
       h_ZllGJets->SetBinContent(bnr, cont_ZllGJets / divby );    
       h_ZllJets ->SetBinContent(bnr, cont_ZllJets  / divby );    
       h_GGJets  ->SetBinContent(bnr, cont_GGJets   / divby );
       h_EFake   ->SetBinContent(bnr, cont_EFake    / divby );
       h_JFake   ->SetBinContent(bnr, cont_JFake    / divby );
       h_WWG     ->SetBinContent(bnr, cont_WWG      / divby );
       //h_WZ      ->SetBinContent(bnr, cont_WZ       / divby );
       h_Wmn     ->SetBinContent(bnr, cont_Wmn      / divby );
       h_Wtn     ->SetBinContent(bnr, cont_Wtn      / divby );
       h_ZZ      ->SetBinContent(bnr, cont_ZZ       / divby );
       h_ZnnGJets->SetBinContent(bnr, cont_ZnnGJets / divby );
     
       data_obs ->SetBinError(bnr, cont_data_obs*err_ratio / divby );    
     
      }
     }
    } // dobygev

    THStack *thestack = new THStack("thestack","");
    thestack->Add( h_GJets    );
    thestack->Add( h_TTGJets  );
    //thestack->Add( h_TGJets  );
    thestack->Add( h_WlnGJets );
    thestack->Add( h_ZllGJets );
    thestack->Add( h_ZllJets  );
    thestack->Add( h_GGJets   );
    thestack->Add( h_EFake    );
    //thestack->Add( h_JFake    );
    thestack->Add( h_WWG      );
    //thestack->Add( h_WZ       );
    thestack->Add( h_Wmn      );
    thestack->Add( h_Wtn      );
    thestack->Add( h_ZZ       );
    thestack->Add( h_ZnnGJets );

    TLegend *leg = new TLegend(0.60,0.6,0.88,0.88 );
    leg->SetFillColor(kWhite);
    leg->AddEntry( data_obs  ,"data", "lep");
    leg->AddEntry( h_GJets   ,"#gamma+jets", "f");
    leg->AddEntry( h_TTGJets ,"t#bar{t}#gamma+jets", "f");
    //leg->AddEntry( h_TGJets  ,"t#gamma+jets", "f");
    leg->AddEntry( h_WlnGJets,"W(l#nu)#gamma+jets", "f");
    leg->AddEntry( h_ZllGJets,"Z(ll)#gamma+jets", "f");
    leg->AddEntry( h_ZllJets ,"Z(ll)+jets", "f");
    leg->AddEntry( h_GGJets  , "#gamma#gamma+jets","f" );
    leg->AddEntry( h_EFake   , "e#rightarrow fake","f" );
    //leg->AddEntry( h_JFake   , "jet#rightarrow fake","f" );
    leg->AddEntry( h_WWG     , "WW#gamma","f" );
    //leg->AddEntry( h_WZ      , "WZ","f" );
    leg->AddEntry( h_Wmn     , "W(#mu#nu)","f" );
    leg->AddEntry( h_Wtn     , "W(#tau#nu)","f" );
    leg->AddEntry( h_ZZ      , "ZZ","f" );
    leg->AddEntry( h_ZnnGJets, "Z(#nu#nu)#gamma+jets","f" );

    canvas->cd();

    data_obs->SetXTitle(TString(data_obs->GetTitle())+" [GeV]");
    data_obs->SetTitle("");
    data_obs->SetYTitle("Events / Bin");
    data_obs->SetMaximum(1000);
    if(dobygev){
     data_obs->SetYTitle("Events / GeV");
     data_obs->SetMaximum(10);
     data_obs->SetMinimum(0.001);
    }
    data_obs->GetYaxis()->SetTitleOffset(1.4);

    data_obs->Draw("");
    thestack->Draw("hist,sames");
    data_obs->Draw("sames");
    leg->Draw("same");

    title->DrawTextNDC(0.17,0.87,"CMS");
    extra->DrawTextNDC(0.17,0.81,"Preliminary");
    lumi->DrawTextNDC(0.9,0.91,"12.9 /fb (13 TeV)");

    data_obs  ->Write();
    h_GJets   ->Write();   
    h_TGJets  ->Write();   
    h_TTGJets ->Write();   
    h_WlnGJets->Write();   
    h_ZllGJets->Write();   
    h_ZllJets ->Write();  
    h_GGJets  ->Write();
    h_EFake   ->Write();
    h_JFake   ->Write();
    h_WWG     ->Write();
    //h_WZ      ->Write();
    h_Wmn     ->Write();
    h_Wtn     ->Write();
    h_ZZ      ->Write();
    h_ZnnGJets->Write();

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
    //h_WZ      ->Delete();
    h_Wmn     ->Delete();
    h_Wtn     ->Delete();
    h_ZZ      ->Delete();
    h_ZnnGJets->Delete();

    file_out->Close();

    } //  for(unsigned int p=0; p<ptrangenames.size(); ++p)
   } //  for(unsigned int s=0; s<selectionnames.size(); ++s)
  } //  for(unsigned int v=0; v<variablenames.size(); ++v)
 } //  for(unsigned int l=0; l<leptons.size(); ++l)

}
