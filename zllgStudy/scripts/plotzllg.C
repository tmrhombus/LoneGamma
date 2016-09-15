#define plotzllg_cxx 
#include "plotzllg.h"

void plotzllg::Loop()
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

 name_SinglePhoton  =  "analyzed_SinglePhoton.root";
 name_GJets         =  "analyzed_GJets.root";
 name_TTGJets       =  "analyzed_TTGJets.root";
 //name_WlnGJets      =  "analyzed_WlnGJets.root";
 name_ZllGJets      =  "analyzed_ZllGJets.root";
 name_ZllJets       =  "analyzed_ZllJets.root";

 file_SinglePhoton  =  new TFile(inpath+"/"+name_SinglePhoton,"READ");
 file_GJets         =  new TFile(inpath+"/"+name_GJets,"READ");
 file_TTGJets       =  new TFile(inpath+"/"+name_TTGJets,"READ");
 //file_WlnGJets      =  new TFile(inpath+"/"+name_WlnGJets,"READ");
 file_ZllGJets      =  new TFile(inpath+"/"+name_ZllGJets,"READ");
 file_ZllJets       =  new TFile(inpath+"/"+name_ZllJets,"READ");

 variablenames.push_back("et"        ); 
 variablenames.push_back("uncorret"  );
 variablenames.push_back("eta"       ); 
 variablenames.push_back("sieieF5x5" ); 
 variablenames.push_back("pfMET"     ); 
 variablenames.push_back("leptoMET"  ); 
 variablenames.push_back("dilep_mass"); 
 variablenames.push_back("dimu_mass" ); 
 variablenames.push_back("diele_mass"); 

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

 selectionnames.push_back("m110");
 selectionnames.push_back("m170");

 ptrangenames.push_back("175to1000");

 dolog = false;

 extraname = "";
 if(dolog){extraname+="_log";}

 TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500); 
 gStyle->SetOptStat(0);
 gPad->SetLogy(dolog);
 gPad->SetTickx();
 gPad->SetTicky();
 gStyle->SetLineWidth(3);

 for(unsigned int v=0; v<variablenames.size(); ++v){
 TString variablename = variablenames[v];
 Int_t rebin = rebins[v];
  for(unsigned int s=0; s<selectionnames.size(); ++s){
  TString selectionname = selectionnames[s];
   for(unsigned int p=0; p<ptrangenames.size(); ++p){
   TString ptrangename = ptrangenames[p];

   canvas->Clear();

   TString tailname = variablename+"_"+ptrangename+"_"+selectionname;
   outname = "plot_"+tailname+extraname; 
   histname = "h_sig_"+tailname;
   file_out =  new TFile(outpath+"/"+outname+".root","RECREATE");

   printf(outname+"\n");

   data_obs = (TH1F*)file_SinglePhoton->Get(histname)->Clone("data_obs");
   
   h_GJets    = (TH1F*)file_GJets   ->Get(histname)->Clone("h_GJets")    ;   
   h_TTGJets  = (TH1F*)file_TTGJets ->Get(histname)->Clone("h_TTGJets")  ;   
   //h_WlnGJets = (TH1F*)file_WlnGJets->Get(histname)->Clone("h_WlnGJets") ;   
   h_ZllGJets = (TH1F*)file_ZllGJets->Get(histname)->Clone("h_ZllGJets") ;   
   h_ZllJets  = (TH1F*)file_ZllJets ->Get(histname)->Clone("h_ZllJets")  ;  

   data_obs  ->SetMarkerStyle( 20 );  
   data_obs  ->SetMarkerSize( 1 );  

   h_GJets   ->SetLineColor( kBlack ); 
   h_TTGJets ->SetLineColor( kBlack ); 
   //h_WlnGJets->SetLineColor( kBlack ); 
   h_ZllGJets->SetLineColor( kBlack ); 
   h_ZllJets ->SetLineColor( kBlack ); 

   data_obs  ->Rebin( rebin ); 
   h_GJets   ->Rebin( rebin ); 
   h_TTGJets ->Rebin( rebin ); 
   //h_WlnGJets->Rebin( rebin ); 
   h_ZllGJets->Rebin( rebin ); 
   h_ZllJets ->Rebin( rebin ); 

   h_GJets   ->SetLineWidth( 2 ); 
   h_TTGJets ->SetLineWidth( 2 ); 
   //h_WlnGJets->SetLineWidth( 2 ); 
   h_ZllGJets->SetLineWidth( 2 ); 
   h_ZllJets ->SetLineWidth( 2 ); 

   //hr_SinglePhotonData  ->SetFillColor( kBlack    ); 
   h_GJets   ->SetFillColor( kRed-10   );  
   h_TTGJets ->SetFillColor( kAzure-9  ); 
   //h_WlnGJets->SetFillColor( kRed-6    );  
   h_ZllGJets->SetFillColor( kOrange-4 );
   h_ZllJets ->SetFillColor( kBlue-8   );  

//   if(variablename=="et"      ||
//      variablename=="uncorret"||
//      variablename=="leptoMET"
//     ){
//    Double_t xbins[6] = {175,190,250,400,700,1000};
//    for(int i=0; i<5; i++){
//     int bnr = i+1;
//     double divby = xbins[i+1] - xbins[i];
//  
//     double cont_data_obs = data_obs  ->GetBinContent(bnr); 
//     double cont_GJets    = h_GJets   ->GetBinContent(bnr); 
//     double cont_TTGJets  = h_TTGJets ->GetBinContent(bnr); 
//     //double cont_WlnGJets = h_WlnGJets->GetBinContent(bnr); 
//     double cont_ZllGJets = h_ZllGJets->GetBinContent(bnr); 
//     double cont_ZllJets  = h_ZllJets ->GetBinContent(bnr); 
//  
//     data_obs  ->SetBinContent(bnr, cont_data_obs   / divby );    
//     h_GJets   ->SetBinContent(bnr, cont_GJets    / divby );    
//     h_TTGJets ->SetBinContent(bnr, cont_TTGJets  / divby );    
//     //h_WlnGJets->SetBinContent(bnr, cont_WlnGJets / divby );    
//     h_ZllGJets->SetBinContent(bnr, cont_ZllGJets / divby );    
//     h_ZllJets ->SetBinContent(bnr, cont_ZllJets  / divby );    
//  
//     data_obs ->SetBinError(bnr, sqrt(cont_data_obs / divby) );    
//  
//    }
//   }

   THStack *thestack = new THStack("thestack","");
   thestack->Add( h_GJets    );
   thestack->Add( h_TTGJets  );
   //thestack->Add( h_WlnGJets );
   thestack->Add( h_ZllGJets );
   thestack->Add( h_ZllJets  );

   TLegend *leg = new TLegend(0.60,0.6,0.88,0.88 );
   leg->SetFillColor(kWhite);
   leg->AddEntry( data_obs  ,"data", "lep");
   leg->AddEntry( h_GJets   ,"#gamma+jets", "f");
   leg->AddEntry( h_TTGJets ,"t#bar{t}#gamma+jets", "f");
   //leg->AddEntry( h_WlnGJets,"W(l#nu)#gamma+jets", "f");
   leg->AddEntry( h_ZllGJets,"Z(ll)#gamma+jets", "f");
   leg->AddEntry( h_ZllJets ,"Z(ll)+jets", "f");

   canvas->cd();

   data_obs->Draw("");
   thestack->Draw("hist,sames");
   data_obs->Draw("sames");
   leg->Draw("same");

   h_GJets   ->Write();   
   h_TTGJets ->Write();   
   //h_WlnGJets->Write();   
   h_ZllGJets->Write();   
   h_ZllJets ->Write();  

   canvas->SaveAs(outpath+"/"+outname+".pdf");

   thestack  ->Delete();
   h_GJets   ->Delete();   
   h_TTGJets ->Delete();   
   //h_WlnGJets->Delete();   
   h_ZllGJets->Delete();   
   h_ZllJets ->Delete();  

   file_out->Close();

   } //  for(unsigned int p=0; p<ptrangenames.size(); ++p)
  } //  for(unsigned int s=0; s<selectionnames.size(); ++s)
 } //  for(unsigned int v=0; v<variablenames.size(); ++v)




  //  
  //  #canvas attributes
  //  canx = 1200
  //  cany = 1100 
  //  gStyle.SetOptStat('')
  //  gStyle.SetLineWidth(3)
  //  gStyle.SetPadTickY(1)
  //  
  //  c_ntot = ROOT.kBlack
  //  c_zllmc = ROOT.kOrange-4                    #zg->llg
  //  c_zlljmc = ROOT.TColor.GetColor("#FDABCC")  #zj->llj
  //  c_znnmc = ROOT.kRed-10                      #zj->nnj
  //  fs = 1001
  //  # line style
  //  n_ls = 1
  //  d_ls = 7
  //  
  //  # TLatex
  //  tex = ROOT.TLatex()
  //  tex.SetTextSize(0.07)
  //  tex.SetTextAlign(13)
  //  tex.SetNDC(True)
  //  
  //  rebin = 1
  //  
  //  c = TCanvas('c','Canvas Named c',canx,cany)
  //  c.cd()
  //  c.SetLogy(do_log)
  //  c.SetFrameLineWidth(2)
  //  c.SetCanvasSize(canx,cany)
  //  
  //  # Data
  //  data_obs = d_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
  //  data_obs.SetName("data_obs")
  //  data_obs.SetTitle("")
  //  data_obs.SetLineColor(c_ntot)
  //  data_obs.SetLineWidth(4)
  //  data_obs.SetFillStyle(0)
  //  data_obs.SetLineStyle(n_ls)
  //  data_obs.Draw("GOFF")
  //  #data_obs.Scale( 1. / max(data_obs.Integral(-1,-1),1.) )
  //  data_obs.Write()
  //  themax = data_obs.GetMaximum()
  //  
  //  # MC (ZG->llG)
  //  h_zllgmc = m_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
  //  h_zllgmc.SetName("h_zllgmc")
  //  h_zllgmc.SetTitle("")
  //  h_zllgmc.SetLineColor(1)
  //  h_zllgmc.SetLineWidth(2)
  //  h_zllgmc.SetFillStyle(fs)
  //  h_zllgmc.SetLineStyle(n_ls)
  //  h_zllgmc.SetFillColor(c_zllmc)
  //  h_zllgmc.Scale(sf)
  //  h_zllgmc.Draw("GOFF")
  //  h_zllgmc.Write()
  //  themax = max( themax, h_zllgmc.GetMaximum() )
  //  
  //  # MC (ZG->nnJ)
  //  h_znnmc = n_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
  //  h_znnmc.SetName("h_znnmc")
  //  h_znnmc.SetTitle("")
  //  h_znnmc.SetLineColor(1)
  //  h_znnmc.SetLineWidth(2)
  //  h_znnmc.SetFillStyle(fs)
  //  h_znnmc.SetLineStyle(n_ls)
  //  h_znnmc.SetFillColor(c_znnmc)
  //  h_znnmc.Scale(sf)
  //  h_znnmc.Draw("GOFF")
  //  h_znnmc.Write()
  //  themax = max( themax, h_znnmc.GetMaximum() )
  //  
  //  # MC (ZJ->llJ)
  //  h_zlljmc = l_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
  //  h_zlljmc.SetName("h_zlljmc")
  //  h_zlljmc.SetTitle("")
  //  h_zlljmc.SetLineColor(1)
  //  h_zlljmc.SetLineWidth(2)
  //  h_zlljmc.SetFillStyle(fs)
  //  h_zlljmc.SetLineStyle(n_ls)
  //  h_zlljmc.SetFillColor(c_zlljmc)
  //  h_zlljmc.Scale(sf)
  //  h_zlljmc.Draw("GOFF")
  //  h_zlljmc.Write()
  //  themax = max( themax, h_zlljmc.GetMaximum() )
  //  
  //  s_mc = THStack('s_mc','')
  //  s_mc.Add(h_znnmc)
  //  s_mc.Add(h_zlljmc)
  //  s_mc.Add(h_zllgmc)
  //  
  //  # Nr Entries MC (ZG->llG)
  //  h_zllent = m_file.Get("h_sig_%s_%s_%s"%(variable,ptrange,selection))
  //  h_zllent.SetName("h_zllent")
  //  nrentries = h_zllent.GetEntries()
  //  h_zllent.SetTitle("")
  //  h_zllent.SetLineStyle(n_ls)
  //  h_zllent.Scale( nrentries / max(h_zllent.Integral(),1.) )
  //  h_zllent.SetLineColor(1)
  //  h_zllent.SetLineWidth(2)
  //  h_zllent.SetFillStyle(fs)
  //  h_zllent.SetLineStyle(n_ls)
  //  h_zllent.SetFillColor(c_zllmc)
  //  h_zllent.Draw("GOFF")
  //  thegmax = h_zllent.GetMaximum()
  //  
  //  # GEN Level (ZG->llG)
  //  h_zllgen = m_file.Get("h_gen_%s_%s_%s"%(variable,ptrange,selection))
  //  h_zllgen.SetName("h_zllgen")
  //  h_zllgen.SetTitle("")
  //  h_zllgen.SetLineColor(c_ntot)
  //  h_zllgen.SetLineWidth(4)
  //  h_zllgen.SetFillStyle(0)
  //  h_zllgen.Draw("GOFF")
  //  thegmax = max( thegmax, h_zllgen.GetMaximum() )
  //  
  //  # fill legends
  //  leg=TLegend(0.55,0.70,0.88,0.88)
  //  leg.AddEntry(data_obs,"Data")
  //  leg.AddEntry(h_zllgmc,"Z(ll)#gamma MC","f")
  //  leg.AddEntry(h_zlljmc,"Z(ll)Jets MC","f")
  //  leg.AddEntry(h_znnmc,"t#bar{t}Jets MC","f")
  //  leg.SetFillColor(0)
  //  leg.SetBorderSize(0)
  //  
  //  legent=TLegend(0.55,0.70,0.88,0.88)
  //  legent.AddEntry(h_zllgen,"Gen Level")
  //  legent.AddEntry(h_zllent,"Reco Level","f")
  //  legent.SetFillColor(0)
  //  legent.SetBorderSize(0)
  //  
  //  # and draw
  //  
  //  
  //  ####### Data v MC ##############
  //  #data_obs.SetMaximum( 5 )
  //  data_obs.SetMaximum( 1.4*themax )
  //  
  //  data_obs.Draw("GOFF")
  //  #h_zllgmc.Scale(data_obs.Integral()/h_zllgmc.Integral())
  //  #h_zllgmc.Scale(0.041*1.3/117.864)
  //  
  //  print("Data:  %.3f"%data_obs.Integral())
  //  print("MC:    %.3f"%h_zllgmc.Integral())
  //  print("Bkg:   %.3f"%h_zlljmc.Integral())
  //  
  //  data_obs.GetYaxis().SetTitle("Events / 2 GeV")
  //  #data_obs.GetXaxis().SetTitle("dilepton mass (#mu#mu)")
  //  data_obs.GetXaxis().SetTitle("dilepton mass (#mu#mu,ee)")
  //  
  //  #if(do_log):  data_obs.SetMaximum( 50*themax ) 
  //  if(do_log):  data_obs.SetMaximum( 11 ) 
  //  c.cd()
  //  data_obs.Draw('GOFF')
  //  s_mc.Draw('hist,GOFF,sames')
  //  #h_zlljmc.Draw('hist,GOFF,sames')
  //  #h_zllgmc.Draw('hist,GOFF,sames')
  //  leg.Draw('sames,GOFF')
  //  data_obs.Draw('GOFF,sames')
  //  
  //  cpr.prelim_tdr(lumi=12900,hor=0.17)
  //  tex.SetTextSize(0.03)
  //  tex.SetTextAlign(11) #left, bottom
  //  #tex.DrawLatex(0.13,0.75,"pT Range [GeV]: %s"%(ptrange))
  //  
  //  c.Update()
  //  time.sleep(1)
  //  
  //  c.SaveAs("%s/%s.pdf"%(outdir,out_filename))
  //  #c.SaveAs("%s/pdf/%s.pdf"%(webdir,out_filename))
  //  #c.SaveAs("%s/png/%s.png"%(webdir,out_filename))
  //  print(  "%s/%s.png"%(www,out_filename))
  //  
  //  
  //  outfile.Close()
  //  
  //  
  //  ####### Gen v Raw Reco ##############
  //  #h_zllgen.SetMaximum( 1.4*thegmax )
  //  #
  //  #h_zllgen.Draw("GOFF")
  //  ##h_zllgmc.Scale(data_obs.Integral()/h_zllgmc.Integral())
  //  ##h_zllgmc.Scale(0.041*1.3/117.864)
  //  #
  //  #print("Gen:         %5.1f"%h_zllgen.Integral())
  //  #print("Raw Reco:    %5.1f"%h_zllent.Integral())
  //  #print("Raw Reco NE: %5.1f"%nrentries)
  //  #
  //  #h_zllgen.GetYaxis().SetTitle("Events / 2 GeV")
  //  ##h_zllgen.GetXaxis().SetTitle("dilepton mass (#mu#mu)")
  //  #h_zllgen.GetXaxis().SetTitle("dilepton mass (#mu#mu,ee)")
  //  #
  //  #c.cd()
  //  #h_zllgen.Draw('GOFF')
  //  #h_zllent.Draw('hist,GOFF,sames')
  //  #legent.Draw('sames,GOFF')
  //  #h_zllgen.Draw('GOFF,sames')
  //  #
  //  #cpr.prelim_tdr(lumi=2300,hor=0.17)
  //  #tex.SetTextSize(0.03)
  //  #tex.SetTextAlign(11) #left, bottom
  //  ##tex.DrawLatex(0.13,0.75,"pT Range [GeV]: %s"%(ptrange))
  //  #
  //  #c.Update()
  //  #time.sleep(1)
  //  #
  //  #c.SaveAs("%s/gen%s.pdf"%(outdir,out_filename))
  //  #c.SaveAs("%s/pdf/gen%s.pdf"%(webdir,out_filename))
  //  #c.SaveAs("%s/png/gen%s.png"%(webdir,out_filename))
  //  #print(  "%s/gen%s.png"%(www,out_filename))







/////    // Ratio Plot
/////    //////////////////////////////////////////////////////////////////////////////
/////    TCanvas* c2 = new TCanvas("c2","c2",900,100,500,500);   
/////  
/////  
/////    // 2015 numbers
/////    Double_t size_num_uncor[4] = {20931,34932,12258,1227 };
/////    Double_t size_num_corr[4] = {594.038,707.621,163.113,7.12871 };
/////    Double_t size_den[4] = {10450,16818,4896,352 };
/////    Double_t err_num_corr[4] = {144.675,186.901,110.716,35.0286 };
/////    Double_t err_num_uncor[4] = {6.95299e-310,6.90578e-310,6.90578e-310,6.95299e-310 };
/////    Double_t err_den[4] = {102.225,129.684,69.9714,18.7617}; 
/////  
/////   Double_t tops[4] = {190,250,400,1000};
/////   Double_t bots[4] = {175,190,250,400};
/////   Double_t mid   [4]        = {181.979, 212.927, 292.708, 474.476  };
/////   //Double_t miderr[4]        = {4.31802, 16.5006, 36.2155, 69.6858  };
/////   Double_t *miderrup = new Double_t[4];
/////   Double_t *miderrdn = new Double_t[4];
/////  
/////   miderrup[0] = tops[0]-mid[0];
/////   miderrup[1] = tops[1]-mid[1];
/////   miderrup[2] = tops[2]-mid[2];
/////   miderrup[3] = tops[3]-mid[3];
/////  
/////   miderrdn[0] = mid[0] - bots[0] ;
/////   miderrdn[1] = mid[1] - bots[1] ;
/////   miderrdn[2] = mid[2] - bots[2] ;
/////   miderrdn[3] = mid[3] - bots[3] ;
/////  
/////  //   Double_t x[n]   = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
/////  //   Double_t y[n]   = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
/////  //   Double_t exl[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
/////  //   Double_t eyl[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
/////  //   Double_t exh[n] = {.02,.08,.05,.05,.03,.03,.04,.05,.06,.03};
/////  //   Double_t eyh[n] = {.6,.5,.4,.3,.2,.2,.3,.4,.5,.6};
/////  //   gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
/////  
/////  
/////    gStyle->SetOptStat(0);
/////    gPad->SetLogy();
/////    gPad->SetTickx();
/////    gPad->SetTicky();
/////    gStyle->SetLineWidth(3);
/////  
/////    TGraphAsymmErrors *gr_num_uncor = new TGraphAsymmErrors(4,mid,size_num_uncor,miderrdn,miderrup,err_num_uncor,err_num_uncor);
/////    TGraphAsymmErrors *gr_num_corr  = new TGraphAsymmErrors(4,mid,size_num_corr,miderrdn,miderrup,err_num_corr,err_num_corr);
/////    TGraphAsymmErrors *gr_den       = new TGraphAsymmErrors(4,mid,size_den,miderrdn,miderrup,err_den,err_den);
/////  
/////    // Graphs
/////    gr_num_uncor->SetTitle("");
/////    gr_num_uncor->SetMarkerColor(1);
/////    gr_num_uncor->SetLineColor(1);
/////    gr_num_uncor->SetMarkerStyle(21);
/////  
/////    gr_num_corr->SetMarkerColor(2);
/////    gr_num_corr->SetLineColor(2);
/////    gr_num_corr->SetMarkerStyle(22);
/////  
/////    gr_den->SetMarkerColor(3);
/////    gr_den->SetLineColor(3);
/////    gr_den->SetMarkerStyle(23);
/////  
/////  
/////    TH1F *hr = c2->DrawFrame(175.,0.1,1000.,1000000.,"");
/////    hr->SetXTitle("photon pT [GeV]");
/////    hr->SetYTitle("Events"); 
/////    hr->GetYaxis()->SetTitleOffset(1.4);
/////  
/////     TText* title = new TText(1,1,"") ;
/////     title->SetTextSize(0.07);
/////     title->SetTextColor(kBlack);
/////     title->SetTextAlign(13);
/////     title->SetTextFont(62);
/////     title->DrawTextNDC(0.17,0.87,"CMS");
/////  
/////     TText* extra = new TText(1,1,"") ;
/////     extra->SetTextSize(0.05);
/////     extra->SetTextColor(kBlack);
/////     extra->SetTextAlign(13);
/////     extra->SetTextFont(52);
/////     extra->DrawTextNDC(0.17,0.81,"Preliminary");
/////  
/////     TText* lumi = new TText(1,1,"") ;
/////     lumi->SetTextSize(0.05);
/////     lumi->SetTextColor(kBlack);
/////     lumi->SetTextAlign(31);
/////     lumi->SetTextFont(42);
/////     lumi->DrawTextNDC(0.9,0.91,"2.3 /fb (13 TeV)");
/////  
/////    title->DrawTextNDC(0.17,0.87,"CMS");
/////    extra->DrawTextNDC(0.17,0.81,"Preliminary");
/////    lumi->DrawTextNDC(0.9,0.91,"2.3 /fb (13 TeV)");
/////  
/////    TLegend *leg = new TLegend(0.55,0.6,0.88,0.88 );
/////    leg->SetFillColor(kWhite);
/////    //leg->AddEntry( gr_num_uncor,"Numerator", "L");
/////    leg->AddEntry( gr_num_uncor,"Numerator (uncorrected)", "L");
/////    leg->AddEntry( gr_num_corr,"Numerator (corrected)", "L");
/////    leg->AddEntry( gr_den,"Denominator", "L");
/////    leg->Draw("same");
/////  
/////    gr_num_uncor->Draw("sames,P");
/////    gr_num_corr->Draw("sames,P"); //
/////    gr_den->Draw("sames,P");
/////  
/////    c2->SaveAs(outpath+"/Graph_NuNcD_v2.pdf");
/////  
/////    c2->Clear();
}
