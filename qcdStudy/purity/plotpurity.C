#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include <vector>
void plotpurity()
{

        std::vector<TString> ptbins;
        ptbins.push_back("90to120");
        ptbins.push_back("120to155");
        ptbins.push_back("155to175");
        ptbins.push_back("175to190");
        ptbins.push_back("190to250");
        ptbins.push_back("250to400");
        ptbins.push_back("400to700");
        ptbins.push_back("700to1000");

        std::vector<TString> isovars;
        isovars.push_back("wchiso");
        isovars.push_back("chiso");

        TString inpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/SezL30met/analyzed";
        TString outpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/SezL30met/plots";
        /////TString inpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/SezAllmet/analyzed";
        /////TString outpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/SezAllmet/plots";
        ////TString inpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/DombAllmet/analyzed";
        ////TString outpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/DombAllmet/plots";
        ///TString inpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/Domb/analyzed";
        ///TString outpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/Domb/plots";
        //TString inpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/Comb/analyzed";
        //TString outpath = "/afs/hep.wisc.edu/cms/tperry/LoneG_slc6_491_CMSSW_7_4_14/src/LoneGamma/qcdStudy/Comb/plots";
        TString inname_GJ  = inpath+"/purity_GJets_Merged.root";
        TString inname_QCD = inpath+"/purity_QCD_Merged.root";
	
	TFile *infile_GJ  = new TFile(inname_GJ);
	TFile *infile_QCD = new TFile(inname_QCD);
	
	Int_t fillcolor = 3;
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetFrameLineWidth(3);
	gStyle->SetLineWidth(2);

        TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500);
        gStyle->SetOptStat(0);
        //gPad->SetLogy();
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

	
        //for(std::vector<TString>::iterator it = ptbins.begin(); it != ptbins.end(); ++it) {
        for(unsigned i=0; i<ptbins.size(); i++) {
         for(unsigned j=0; j<isovars.size(); j++) {

         TString ptrange = ptbins.at(i);
         TString isovar = isovars.at(j);
	 //pt155to190
	 //Signal purity
	 TH1F* h_Nsig_GJ  = (TH1F*)infile_GJ->Get("h_Nsig_"+isovar+"_"+ptrange);
	 TH1F* h_Nsig_QCD = (TH1F*)infile_QCD->Get("h_Nsig_"+isovar+"_"+ptrange);
	 TH1F* h_Nsig  = (TH1F*)h_Nsig_GJ->Clone("h_Nsig");
         h_Nsig->Add(h_Nsig_QCD);
	 //Background purity
	 TH1F* h_Nbkg_GJ  = (TH1F*)infile_GJ->Get("h_Nbkg_"+isovar+"_"+ptrange);
	 TH1F* h_Nbkg_QCD = (TH1F*)infile_QCD->Get("h_Nbkg_"+isovar+"_"+ptrange);
	 TH1F* h_Nbkg  = (TH1F*)h_Nbkg_GJ->Clone("h_Nbkg");
         h_Nbkg->Add(h_Nbkg_QCD);
	 //Denominator
	 TH1F* h_Deno_GJ  = (TH1F*)infile_GJ->Get("h_Deno_"+isovar+"_"+ptrange);
	 TH1F* h_Deno_QCD = (TH1F*)infile_QCD->Get("h_Deno_"+isovar+"_"+ptrange);
	 TH1F* h_Deno  = (TH1F*)h_Deno_GJ->Clone("h_Deno");
         h_Deno->Add(h_Deno_QCD);
	 
	 //Define purity
	 TGraphAsymmErrors *purity_signal     = new TGraphAsymmErrors(h_Nsig,h_Deno);
	 //purity_signal->GetYaxis()->SetTitle("Purity");
	 //purity_signal->GetXaxis()->SetTitle("Worst Charged Hadron Isolation");
	 purity_signal->SetLineWidth(2);
	 purity_signal->SetLineColor(4);
	 purity_signal->SetMarkerColor(4);
	 purity_signal->SetMarkerStyle(20);
	 
	 TGraphAsymmErrors *purity_background = new TGraphAsymmErrors(h_Nbkg,h_Deno);
	 //purity_background->GetYaxis()->SetTitle("Purity");
	 //purity_background->GetXaxis()->SetTitle("Worst Charged Hadron Isolation");
	 purity_background->SetLineWidth(2);
	 purity_background->SetLineColor(2);
	 purity_background->SetMarkerColor(2);
	 purity_background->SetMarkerStyle(24);
 	
         TLegend *leg = new TLegend(0.5,0.75,0.88,0.88);
         leg->SetTextSize(0.04) ;
         leg->SetFillColor(0); leg->SetShadowColor(0);leg->SetBorderSize(0);
         leg->SetTextFont(22.);
         leg->AddEntry(purity_signal,"Signal","lep");
         leg->AddEntry(purity_background,"Background","lep");

         float xmin = 0.;
         float xmax = 25.;
         TLine *lowline = new TLine(xmin,0,xmax,0);
         TLine *hiline = new TLine(xmin,1,xmax,1);
         lowline->SetLineWidth(1); lowline->SetLineColor(1); lowline->SetLineStyle(3);         
         hiline->SetLineWidth(1); hiline->SetLineColor(1); hiline->SetLineStyle(3);         

         TH1F *hframe = canvas->DrawFrame(xmin,-0.05,xmax,1.4,"");
         lowline->Draw();
         hiline->Draw();
         title->DrawTextNDC(0.17,0.87,"CMS");
         extra->DrawTextNDC(0.17,0.81,"Preliminary");
         lumi->DrawTextNDC(0.9,0.91,"2.24 /fb (13 TeV)");
         hframe->SetXTitle(h_Nsig->GetTitle());
         hframe->SetYTitle("Purity");
	 
	 //Draw
	 purity_signal->Draw("P");
	 purity_background->Draw("P");
	 leg->SetHeader("photon p_{T}: "+ptrange);
	 leg->Draw();

         canvas->SaveAs(outpath+"/Purity_"+isovar+"_"+ptrange+".pdf");
         canvas->Clear();
         }
        }
}

