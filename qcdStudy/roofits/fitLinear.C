#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>


void fitLinear()
{
   // fitting QCD fake ratio

   Int_t n = 4;
   Double_t *x = new Double_t[n];
   Double_t *y = new Double_t[n];
   Double_t *ex = new Double_t[n];

   x[0]=182.5;
   x[1]=220.;
   x[2]=325.;
   x[3]=700.;
 
   ex[0]=7.5  ;
   ex[1]=30.  ;
   ex[2]=75.  ;
   ex[3]=300. ;
 
   y[0] = 0.0901838;
   y[1] = 0.0721998;
   y[2] = 0.0669271;
   y[3] = 0.0753123;

   TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500);
 
   gStyle->SetOptStat(0);
   gPad->SetLogy(kFALSE);
   gPad->SetTickx();
   gPad->SetTicky();
   gStyle->SetLineWidth(3);
 
   TH1F *hr = canvas->DrawFrame(175.,0.,1000.,0.5,"");
   hr->SetXTitle("photon pT [GeV]");
   hr->SetYTitle("Fake Ratio"); 
 
   TGraphErrors *gr_cntrl = new TGraphErrors(n, x, y, 0, 0);
   TGraphErrors *bu_cntrl = new TGraphErrors(n, x, y, ex, 0);
   //hr->cd();
   bu_cntrl->Draw("P,same");
   gr_cntrl->Draw("P*,same");
   //Fit the graph with the predefined "pol3" function
   gr_cntrl->Fit("pol1");
   //Access the fit resuts
   TF1 *fit_cntrl = gr_cntrl->GetFunction("pol1");
   fit_cntrl->SetLineWidth(2);
   Double_t chi2 = fit_cntrl->GetChisquare();
   Double_t p0 = fit_cntrl->GetParameter(0);
   Double_t e0 = fit_cntrl->GetParError(0);
   Double_t p1 = fit_cntrl->GetParameter(1);
   Double_t e1 = fit_cntrl->GetParError(1);
   TString schi2  = TString(boost::lexical_cast<string>( boost::format("%.5f") % chi2 ));
   TString sp0  = TString(boost::lexical_cast<string>( boost::format("%.3f") % p0 ));
   TString se0  = TString(boost::lexical_cast<string>( boost::format("%.3f") % e0 ));
   TString sp1  = TString(boost::lexical_cast<string>( boost::format("%.3f") % p1 ));
   TString se1  = TString(boost::lexical_cast<string>( boost::format("%.3f") % e1 ));

 
   TLegend *leg2 = new TLegend(0.55,0.6,0.88,0.88 );
   leg2->SetFillColor(kWhite);
   leg2->AddEntry( gr_cntrl,"QCD Fake Ratio", "L");
   leg2->AddEntry( fit_cntrl,"best fit result`", "L");
   leg2->Draw("same");

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
   //xframe->addObject(extra);
 
   TText* lumi = new TText(1,1,"") ;
   lumi->SetTextSize(0.05);
   lumi->SetTextColor(kBlack);
   lumi->SetTextAlign(31);
   lumi->SetTextFont(42);
   lumi->DrawTextNDC(0.9,0.91,"2.1 /fb (13 TeV)");
   //xframe->addObject(lumi);

   TText* fitparams = new TText(1,1,"") ;
   fitparams->SetTextSize(0.05);
   fitparams->SetTextColor(kBlack);
   fitparams->SetTextAlign(11);
   fitparams->SetTextFont(42);
   fitparams->DrawTextNDC(0.15,0.6,"Chi2: "+schi2);
   fitparams->DrawTextNDC(0.2,0.5,"y = m x + b");
   fitparams->DrawTextNDC(0.2,0.45,"m = "+sp1+" +- "+se1);
   fitparams->DrawTextNDC(0.2,0.4,"b = "+sp0+" +- "+se0);


   canvas->Update();
   //canvas->SaveAs(outpath+"/Graph_NuNcD.C");
   //canvas->SaveAs(outpath+"/Graph_NuNcD.pdf");
   //canvas->SaveAs(outpath+"/Graph_NuNcD.eps");
   //canvas->SaveAs(wwwpath+"/Graph_NuND.png");
   canvas->SaveAs("test.png");
 
   canvas->Clear();

}

