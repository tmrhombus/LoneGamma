#define plotfakeratio_cxx 
#include "plotfakeratio.h"

void plotfakeratio::Loop()
{

 char* submitbase;
 submitbase = getenv ("submitbase");
 char* version;
 version = getenv("version");
 Tsubmitbase = TString(submitbase);
 Tversion = TString(version);

 outpath = TString(Tsubmitbase+"/gitignore/"+Tversion+"/plots");
 wwwpath = TString("/afs/hep.wisc.edu/user/tperry/www/MonoPhoton/qcdPlots/"+Tversion);


std::vector<Double_t> ratios;
std::vector<Double_t> ratioerrors;

  // Ratio Plot
  //////////////////////////////////////////////////////////////////////////////
  TCanvas* c2 = new TCanvas("c2","c2",900,100,500,500);   

  // 2016 numbers  5.6 \fb
 Double_t size_num_corr[4]={1534.77, 1639.05, 473.542, 22.9444};
 Double_t size_den[4]={28398, 45728, 13463, 898};  
 Double_t err_num_corr[4]={224.085, 289.638, 174.063, 52.3737};
 Double_t err_den[4]={168.517, 213.841, 116.03, 29.9666};  
 Double_t mid[4]={182.001, 213.101, 293.495, 470.313}; 
 Double_t miderr[4]={4.31492, 16.5964, 36.7422, 65.3309};

 // 2016 numbers (met < 90)  5.6 \fb
 Double_t v1_size_num_corr[4]={2870.1, 3436.77, 996.909, 48.1623};
 Double_t v1_size_den[4]={51970, 86750, 27820, 2236};
 Double_t v1_err_num_corr[4]={288.076, 380.971, 240.813, 81};
 Double_t v1_err_den[4]={227.969, 294.534, 166.793, 47.2864};
 Double_t v1_mid[4]={182.045, 213.422, 295.302, 474.299};
 Double_t v1_miderr[4]={4.32302, 16.641, 37.4739, 67.8315};

//  // 2016 numbers  2.6 \fb
// Double_t size_num_corr[4]={691.317, 709.649, 212.176, 8.6222} ;
// Double_t size_den[4]={12403, 20181, 5860, 392} ;
// Double_t err_num_corr[4]={150.256, 193.82, 115.304, 34}  ;
// Double_t err_den[4]={111.369, 142.06, 76.5506, 19.799} ;
// Double_t mid[4]={181.968, 213.035, 293.528, 468.739} ;
// Double_t miderr[4]={4.33628, 16.5507, 36.6815, 63.956} ;
//
// // 2016 numbers (met < 90)  2.6 \fb
// Double_t v1_size_num_corr[4]={1405.48, 1470.16, 440.318, 21.17}; 
// Double_t v1_size_den[4]={22534, 37597, 11943, 955}; 
// Double_t v1_err_num_corr[4]={192.255, 253.249, 158.688, 53.066}; 
// Double_t v1_err_den[4]={150.113, 193.899, 109.284, 30.9031}; 
// Double_t v1_mid[4]={182.06, 213.344, 295.18, 473.947}; 
// Double_t v1_miderr[4]={4.33965, 16.6033, 37.4313, 66.8485}; 


//  // 2015 numbers
// Double_t v1_size_num_corr[4] = {594.038, 707.621, 163.113, 7.12871 };
// Double_t v1_size_den     [4] = {10450, 16818, 4896, 352   };
// Double_t v1_err_num_corr [4] = {144.675, 186.901, 110.716, 35.0286  };
// Double_t v1_err_den      [4] = {102.225, 129.684, 69.9714, 18.7617  };
// Double_t v1_mid   [4]        = {181.979, 212.927, 292.708, 474.476  };
// Double_t v1_miderr[4]        = {4.31802, 16.5006, 36.2155, 69.6858  };

  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  //Double_t ratio[3];
  // 2016 Numbers
  Double_t *ratio = new Double_t[4];
  ratio[0] = size_num_corr[0]/size_den[0];
  ratio[1] = size_num_corr[1]/size_den[1];
  ratio[2] = size_num_corr[2]/size_den[2];
  ratio[3] = size_num_corr[3]/size_den[3];

  Double_t *ratioerr = new Double_t[4];
  ratioerr[0]= ratio[0] * pow( pow(err_num_corr[0]/size_num_corr[0] ,2) + pow(err_den[0]/size_den[0] ,2) ,0.5); //0.005;
  ratioerr[1]= ratio[1] * pow( pow(err_num_corr[1]/size_num_corr[1] ,2) + pow(err_den[1]/size_den[1] ,2) ,0.5); //0.005;
  ratioerr[2]= ratio[2] * pow( pow(err_num_corr[2]/size_num_corr[2] ,2) + pow(err_den[2]/size_den[2] ,2) ,0.5); //0.005;
  ratioerr[3]= ratio[3] * pow( pow(err_num_corr[3]/size_num_corr[3] ,2) + pow(err_den[3]/size_den[3] ,2) ,0.5); //0.005;

  // 2015 Numbers
  Double_t *v1_ratio = new Double_t[4];
  v1_ratio[0] = v1_size_num_corr[0]/v1_size_den[0];
  v1_ratio[1] = v1_size_num_corr[1]/v1_size_den[1];
  v1_ratio[2] = v1_size_num_corr[2]/v1_size_den[2];
  v1_ratio[3] = v1_size_num_corr[3]/v1_size_den[3];

  Double_t *v1_ratioerr = new Double_t[4];
  v1_ratioerr[0]= v1_ratio[0] * pow( pow(v1_err_num_corr[0]/v1_size_num_corr[0] ,2) + pow(v1_err_den[0]/v1_size_den[0] ,2) ,0.5); //0.005;
  v1_ratioerr[1]= v1_ratio[1] * pow( pow(v1_err_num_corr[1]/v1_size_num_corr[1] ,2) + pow(v1_err_den[1]/v1_size_den[1] ,2) ,0.5); //0.005;
  v1_ratioerr[2]= v1_ratio[2] * pow( pow(v1_err_num_corr[2]/v1_size_num_corr[2] ,2) + pow(v1_err_den[2]/v1_size_den[2] ,2) ,0.5); //0.005;
  v1_ratioerr[3]= v1_ratio[3] * pow( pow(v1_err_num_corr[3]/v1_size_num_corr[3] ,2) + pow(v1_err_den[3]/v1_size_den[3] ,2) ,0.5); //0.005;

  ratios.push_back(ratio[0]);
  ratios.push_back(ratio[1]);
  ratios.push_back(ratio[2]);
  ratios.push_back(ratio[3]);

  ratioerrors.push_back(ratioerr[0]);
  ratioerrors.push_back(ratioerr[1]);
  ratioerrors.push_back(ratioerr[2]);
  ratioerrors.push_back(ratioerr[3]);

  TGraph *gr_ratio_noe = new TGraph(4,mid,ratio);
  TGraphErrors *gr_ratio = new TGraphErrors(4,mid,ratio,miderr,ratioerr);

  TGraph *v1_gr_ratio_noe = new TGraph(4,v1_mid,v1_ratio);
  TGraphErrors *v1_gr_ratio = new TGraphErrors(4,v1_mid,v1_ratio,v1_miderr,v1_ratioerr);


  v1_gr_ratio->SetMarkerColor(4);
  v1_gr_ratio->SetLineColor(4);
  v1_gr_ratio->SetMarkerStyle(21);

  gr_ratio->SetMarkerColor(2);
  gr_ratio->SetLineColor(2);
  gr_ratio->SetMarkerStyle(22);

  //TH1F *hs = c2->DrawFrame(0.,0.,1000.,0.5,"");
  TH1F *hs = c2->DrawFrame(175.,0.,1000.,0.5,"");
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
   lumi->DrawTextNDC(0.9,0.91,"5.6 /fb (13 TeV)");

  title->DrawTextNDC(0.17,0.87,"CMS");
  extra->DrawTextNDC(0.17,0.81,"Preliminary");
  lumi->DrawTextNDC(0.9,0.91,"5.6 /fb (13 TeV)");

  TLegend *leg3 = new TLegend(0.55,0.7,0.88,0.88 );
  leg3->SetFillColor(kWhite);
  leg3->AddEntry( gr_ratio,"Num./Den.", "L");
  //leg3->AddEntry( gr_ratio,"2016 Num./Den. (5.6 /fb)", "L");
  //leg3->AddEntry( v1_gr_ratio,"2015 Num./Den. (2.3 /fb)", "L");
  //leg3->AddEntry( gr_ratio   ,"Num./Den. ( MET < 30 )", "L");
  //leg3->AddEntry( v1_gr_ratio,"Num./Den. ( MET < 90 )", "L");
  leg3->Draw("same");

////  // 2015
//  v1_gr_ratio->Draw("P");
////  v1_gr_ratio_noe->Draw("P,sames");
////  v1_gr_ratio_noe->Fit("pol1");
////  TF1 *v1_fit_ratio = v1_gr_ratio_noe->GetFunction("pol1");
////  v1_fit_ratio->SetLineColor(13);
////  v1_fit_ratio->SetLineWidth(3);
////  Double_t v1_chi2 = v1_fit_ratio->GetChisquare();
////  Double_t v1_p0   = v1_fit_ratio->GetParameter(0);
////  Double_t v1_p1   = v1_fit_ratio->GetParameter(1);
////  Double_t v1_e0   = v1_fit_ratio->GetParError(0);
////  Double_t v1_e1   = v1_fit_ratio->GetParError(1);

  // 2016
  gr_ratio->Draw("P,sames");
  gr_ratio_noe->Draw("P,sames");
  gr_ratio_noe->Fit("pol0");
  //gr_ratio_noe->Fit("pol1");
  TF1 *fit_ratio = gr_ratio_noe->GetFunction("pol0");
  //TF1 *fit_ratio = gr_ratio_noe->GetFunction("pol1");
  fit_ratio->SetLineColor(1);
  fit_ratio->SetLineWidth(3);
  Double_t chi2 = fit_ratio->GetChisquare();
  Double_t p0 = fit_ratio->GetParameter(0);
  Double_t p1 = 0; //fit_ratio->GetParameter(1);
  Double_t e0 = fit_ratio->GetParError(0);
  Double_t e1 = 0; //fit_ratio->GetParError(1);
  leg3->AddEntry( fit_ratio,"Best Fit Line");

  ms.push_back(p1); mes.push_back(e1);
  bs.push_back(p0); bes.push_back(e0);

  //v1_fit_ratio->Draw("sames");
  fit_ratio->Draw("sames");

  TLatex tex;
  tex.SetTextSize(0.07);
  tex.SetTextColor(kBlack);
  tex.SetTextAlign(11);
  tex.SetTextFont(42);
  tex.DrawLatexNDC(0.45,0.47,"y = m x + b");
  tex.SetTextSize(0.05);
  tex.DrawLatexNDC(0.40,0.40,TString(boost::lexical_cast<string>(boost::format("m = %0.5f +- %0.5f") % p1 % e1)));
  tex.DrawLatexNDC(0.45,0.35,TString(boost::lexical_cast<string>(boost::format("b = %0.3f +- %0.3f") % p0 % e0)));
  ///tex.DrawLatexNDC(0.45,0.30,"#chi^{2}"+TString(boost::lexical_cast<string>(boost::format(" = %0.5f") % chi2)));

  c2->Update();

  //c2->SaveAs(outpath+"/Graph_FakeRatio_yearcomp.pdf");
  //c2->SaveAs(wwwpath+"/Graph_FakeRatio_yearcomp.pdf");
  ///c2->SaveAs(outpath+"/Graph_FakeRatio_metcomp.pdf");
  ///c2->SaveAs(wwwpath+"/Graph_FakeRatio_metcomp.pdf");
  c2->SaveAs(outpath+"/Graph_FakeRatio_fit_pol0.pdf");
  c2->SaveAs(wwwpath+"/Graph_FakeRatio_fit_pol0.pdf");

  c2->Clear();
}
