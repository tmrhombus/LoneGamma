#define plotfrerrors_cxx 
#include "plotfrerrors.h"

void plotfrerrors::Loop()
{

 char* submitbase;
 submitbase = getenv ("submitbase");
 char* version;
 version = getenv("version");
 Tsubmitbase = TString(submitbase);
 Tversion = TString(version);

 outpath = TString(Tsubmitbase+"/gitignore/"+Tversion+"/plots");
 wwwpath = TString("/afs/hep.wisc.edu/user/tperry/www/MonoPhoton/qcdPlots/"+Tversion);

  // Ratio Plot
  //////////////////////////////////////////////////////////////////////////////
  TCanvas* c2 = new TCanvas("c2","c2",900,100,500,500);   

  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  // 2016 (2.6/fb) numbers
 
  // standard
  // sideband Up
  // sideband Down
  // MET Up
  // MET Down
  // binning Up
  // binning Down
  // bkg. no #gamma iso req
  // fit Up
  // fit Down
  // Gluon Sideband
  
  Double_t p_00[4]={0.0557379, 0.0351642, 0.0362076, 0.0219954}; 
  Double_t p_01[4]={0.0579247, 0.0358226, 0.0373663, 0.0265233}; 
  Double_t p_02[4]={0.0535033, 0.0349521, 0.0364459, 0.0219313}; 
  Double_t p_03[4]={0.0541997, 0.0336439, 0.0347376, 0.0227814}; 
  Double_t p_04[4]={0.0694237, 0.030723, 0.0324944, 0.0138439}; 
  Double_t p_05[4]={0.0562571, 0.0354248, 0.0358099, 0.021292}; 
  Double_t p_06[4]={0.0519909, 0.033024, 0.0340734, 0.0204509}; 
  Double_t p_07[4]={0.069507, 0.0436358, 0.044886, 0.0265298}; 
  Double_t p_09[4]={0.0595069, 0.0367678, 0.0379102, 0.0251743}; 
  Double_t p_10[4]={0.0519689, 0.0335605, 0.0345049, 0.0188165}; 
  Double_t p_11[4]={0.0623716, 0.0391031, 0.0368683, 0.0221675}; 

  const Int_t NBINS = 4;
  Double_t edges[NBINS+1] = {175,190,250,400,1000};

  // fill the graphs
  TH1* h_sys00 = new TH1D( "h_sys00", "h_sys00", NBINS, edges );
  TH1* h_sys01 = new TH1D( "h_sys01", "h_sys01", NBINS, edges );
  TH1* h_sys02 = new TH1D( "h_sys02", "h_sys02", NBINS, edges );
  TH1* h_sys03 = new TH1D( "h_sys03", "h_sys03", NBINS, edges );
  TH1* h_sys04 = new TH1D( "h_sys04", "h_sys04", NBINS, edges );
  TH1* h_sys05 = new TH1D( "h_sys05", "h_sys05", NBINS, edges );
  TH1* h_sys06 = new TH1D( "h_sys06", "h_sys06", NBINS, edges );
  TH1* h_sys07 = new TH1D( "h_sys07", "h_sys07", NBINS, edges );
  TH1* h_sys09 = new TH1D( "h_sys09", "h_sys09", NBINS, edges );
  TH1* h_sys10 = new TH1D( "h_sys10", "h_sys10", NBINS, edges );
  TH1* h_sys11 = new TH1D( "h_sys11", "h_sys11", NBINS, edges );

  for(unsigned int p=0; p<4; ++p){ h_sys00->SetBinContent( p+1, (p_00[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys01->SetBinContent( p+1, (p_01[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys02->SetBinContent( p+1, (p_02[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys03->SetBinContent( p+1, (p_03[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys04->SetBinContent( p+1, (p_04[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys05->SetBinContent( p+1, (p_05[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys06->SetBinContent( p+1, (p_06[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys07->SetBinContent( p+1, (p_07[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys09->SetBinContent( p+1, (p_09[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys10->SetBinContent( p+1, (p_10[p]-p_00[p])/p_00[p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys11->SetBinContent( p+1, (p_11[p]-p_00[p])/p_00[p] ); }

  Int_t ci1  = 3001; //TColor::GetFreeColorIndex(); // 1001;
  Int_t ci2  = 3002; //TColor::GetFreeColorIndex(); // 1002;
  Int_t ci3  = 3003; //TColor::GetFreeColorIndex(); // 1003;
  Int_t ci4  = 3004; //TColor::GetFreeColorIndex(); // 1004;
  Int_t ci5  = 3005; //TColor::GetFreeColorIndex(); // 1005;
  Int_t ci6  = 3006; //TColor::GetFreeColorIndex(); // 1006;
  Int_t ci7  = 3007; //TColor::GetFreeColorIndex(); // 1007;
  Int_t ci8  = 3008; //TColor::GetFreeColorIndex(); // 1008;
  Int_t ci9  = 3009; //TColor::GetFreeColorIndex(); // 1009;
  Int_t ci10 = 3010; //TColor::GetFreeColorIndex(); // 1010;
  Int_t ci11 = 3011;
  TColor *color1  = new TColor( ci1,178./255,223./255,138./255);
  TColor *color2  = new TColor( ci2, 51./255,160./255, 44./255);
  TColor *color3  = new TColor( ci3,251./255,154./255,153./255);
  TColor *color4  = new TColor( ci4,227./255, 26./255, 28./255);
  TColor *color5  = new TColor( ci5,166./255,206./255,227./255);
  TColor *color6  = new TColor( ci6, 31./255,120./255,180./255);
  TColor *color7  = new TColor( ci7,177./255, 89./255, 40./255);
  TColor *color8  = new TColor( ci8,236./255,208./255,120./255);
  TColor *color9  = new TColor( ci9,253./255,191./255,111./255);
  TColor *color10 = new TColor(ci10,255./255,127./255,  0./255);
  TColor *color11 = new TColor(ci11,127./255,0./255,  255./255);

  h_sys00 ->SetMarkerColor(kBlack                     );
  h_sys01 ->SetMarkerColor(ci1);
  h_sys02 ->SetMarkerColor(ci2);
  h_sys03 ->SetMarkerColor(ci3);
  h_sys04 ->SetMarkerColor(ci4);
  h_sys05 ->SetMarkerColor(ci5);
  h_sys06 ->SetMarkerColor(ci6);
  h_sys07 ->SetMarkerColor(ci7);
  h_sys09 ->SetMarkerColor(ci9);
  h_sys10->SetMarkerColor(ci10);
  h_sys11->SetMarkerColor(ci11);

  h_sys00 ->SetLineColor(kBlack                     );
  h_sys01 ->SetLineColor(ci1);
  h_sys02 ->SetLineColor(ci2);
  h_sys03 ->SetLineColor(ci3);
  h_sys04 ->SetLineColor(ci4);
  h_sys05 ->SetLineColor(ci5);
  h_sys06 ->SetLineColor(ci6);
  h_sys07 ->SetLineColor(ci7);
  h_sys09 ->SetLineColor(ci9);
  h_sys10->SetLineColor(ci10);
  h_sys11->SetLineColor(ci11);

  h_sys00->SetLineWidth(3);
  h_sys01->SetLineWidth(3);
  h_sys02->SetLineWidth(3);
  h_sys03->SetLineWidth(3);
  h_sys04->SetLineWidth(3);
  h_sys05->SetLineWidth(3);
  h_sys06->SetLineWidth(3);
  h_sys07->SetLineWidth(3);
  h_sys09->SetLineWidth(3);
  h_sys10->SetLineWidth(3);
  h_sys11->SetLineWidth(3);

  TH1F *hs = c2->DrawFrame(175.,-0.5,1000.,0.8,"");
  hs->SetXTitle("photon pT [GeV]");
  hs->SetYTitle("Fake Ratio");
  hs->GetYaxis()->SetTitleOffset(1.4);

   TText* title = new TText(1,1,"") ;
   title->SetTextSize(0.07);
   title->SetTextColor(kBlack);
   title->SetTextAlign(13);
   title->SetTextFont(62);
   title->DrawTextNDC(0.17,0.87,"CMS");

   TText* extra = new TText(1,1,"") ;
   extra->SetTextSize(0.05);
   extra->SetTextColor(kBlack);
   extra->SetTextAlign(13);
   extra->SetTextFont(52);
   extra->DrawTextNDC(0.17,0.81,"Preliminary");

   TText* lumi = new TText(1,1,"") ;
   lumi->SetTextSize(0.05);
   lumi->SetTextColor(kBlack);
   lumi->SetTextAlign(31);
   lumi->SetTextFont(42);
   lumi->DrawTextNDC(0.9,0.91,"2.6 /fb (13 TeV)");

  title->DrawTextNDC(0.17,0.87,"CMS");
  extra->DrawTextNDC(0.17,0.81,"Preliminary");
  lumi->DrawTextNDC(0.9,0.91,"2.6 /fb (13 TeV)");

  h_sys00->Draw("hist,same");
  h_sys01->Draw("hist,same");
  h_sys02->Draw("hist,same");
  h_sys03->Draw("hist,same");
  h_sys04->Draw("hist,same");
  h_sys05->Draw("hist,same");
  h_sys06->Draw("hist,same");
  h_sys07->Draw("hist,same");
  h_sys09->Draw("hist,same");
  h_sys10->Draw("hist,same");
  h_sys11->Draw("hist,same");

  TLegend *leg4 = new TLegend(0.5,0.6,0.88,0.88 );
  leg4->SetFillColor(kWhite);
  leg4->AddEntry( h_sys00,"standard", "L");
  leg4->AddEntry( h_sys07,"bkg. no #gamma iso req", "L");
  leg4->AddEntry( h_sys01,"sideband Up", "L");
  leg4->AddEntry( h_sys02,"sideband Down", "L");
  leg4->AddEntry( h_sys03,"MET Up", "L");
  leg4->AddEntry( h_sys04,"MET Down", "L");
  leg4->AddEntry( h_sys05,"binning Up", "L");
  leg4->AddEntry( h_sys06,"binning Down","L");
//  leg4->AddEntry( h_sys08,"ele template", "L");
  leg4->AddEntry( h_sys09,"fit Up","L");
  leg4->AddEntry( h_sys10,"fit Down","L");
  leg4->AddEntry( h_sys11,"gluon sideband","L");
  leg4->Draw("same");

  c2->Update();

  c2->SaveAs(outpath+"/Graph_HistFakeRatioSystUncs.pdf");
  c2->SaveAs(wwwpath+"/Graph_HistFakeRatioSystUncs.pdf");

  c2->Clear();
}
