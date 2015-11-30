#define plot_cxx 
#include "plot.h"

void plot::Loop()
{

char* submitbase;
submitbase = getenv ("submitbase");
char* version;
version = getenv("version");
Tsubmitbase = TString(submitbase);
Tversion = TString(version);

inpath = TString(Tsubmitbase+"/"+Tversion+"/analyzed");
outpath = TString(Tsubmitbase+"/"+Tversion+"/plots");

datafile = new TFile(inpath+"/analyzed_SinglePhoton_2015D.root","READ");
mcfile   = new TFile(inpath+"/analyzed_GJets_Merged.root","READ");
outfile  = new TFile(outpath+"/Fitted.root","RECREATE");

//Vector with QCD fractions
std::vector<double> qcdFrac;
 qcdFrac.clear();
std:: vector<double> qcdFracErr;
 qcdFracErr.clear();

//second method without phojet MC
std::vector<double> qcdFracSam;
 qcdFracSam.clear();
std:: vector<double> qcdFracErrSam;
 qcdFracErrSam.clear();

// Do fitting (fill qcdFrac* vectors)
getFraction(
 datafile,
 mcfile,
 outfile,
 qcdFrac,
 qcdFracErr,
 qcdFracSam,
 qcdFracErrSam   
);

////Using RooFit result
//getCorrectedFakeRatio(
// datafile,
// outfile,
// qcdFrac,
// qcdFracErr
//);

}


//------------------------------------
//  get fraction fits from templates
//------------------------------------
void plot::getFraction(
 TFile* datafile,  
 TFile* mcfile,    
 TFile* outfile,   
 std::vector<double>& fractionQCD,      //<<---This gets filled here
 std::vector<double>& fractionQCDErr,   //<<---This gets filled here
 std::vector<double>& fractionQCDSam,   //<<---This gets filled here
 std::vector<double>& fractionQCDErrSam //<<---This gets filled here
 )
{

  TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500);   

  TH1D* hdata;
  TH1D* hqcd;
  TH1D* hphojet;
  TH1D* hqcdfractionSam;
  TH1D* hqcdfraction;

  double lower[6];   
  string shistname;
  TString histname;

  TString datahistname; //data histo
  TString qcdhistname;  //trkflip - > histos will go here
  TString fullhistname; //phojet
  
  //read variables from .txt file  
  char cxvariable[200];
  ifstream infile_xvar;
  infile_xvar.open("xvariables.list", ifstream::in );


//--------Loop over all the histos------------
  while(!infile_xvar.eof())
 {
    infile_xvar >>cxvariable;
    string xvariable(cxvariable);

    if(strncmp(cxvariable,"#",1)==0) // ignore lines starting with #
      {	continue; }

    shistname = xvariable;
    histname  = shistname; 

    //set template histo names
    datahistname = "h_sig_sieieF5x5_"+histname; //data histos with tigh ID and without sigIetaIeta cut
    fullhistname = "h_sig_sieieF5x5_"+histname; //phojet histos with tight ID cut but without sigIetaIeta
    qcdhistname  = "h_bkg_sieieF5x5_"+histname; //data histos with very loose id and  sideband of track iso


    //get Data template and QCD template from data
    datafile->cd();
    hdata = (TH1D*)datafile->Get(datahistname)->Clone();
    //QCD histo
    hqcd  = (TH1D*)datafile->Get(qcdhistname)->Clone();
    double integ_data = hdata->Integral();
    double integ_qcd = hqcd->Integral();
    //get photon+jet signal tempalte
    mcfile->cd();
    hphojet  = (TH1D*)mcfile->Get(fullhistname)->Clone();
    double integ_phojet = hphojet->Integral();
    int nbins = hphojet->GetNbinsX();

  
  //------Lets try to get with method-2 without MC
  cout<<"-----------------------------------------------------"<<endl;                       
  cout<<"                    Running Sam's Method              "<<endl;
  cout<<"-----------------------------------------------------"<<endl;

    cout<<"bin corresponding ot 0.0102 : "<<hdata->FindBin(0.0102)<<endl;
    int myreqbin = (hphojet->FindBin(0.0102))-1;  ////should be 44 for 0.011; 52 if 0.013 not anymore

    double phojet_gt   = 0.0;

    double data_lt     = hdata->Integral(1,myreqbin);
    double data_gt     = hdata->Integral(myreqbin+1, nbins);
    double data_tot    = hdata->Integral();
    
    double uqcd_gt     = hqcd->Integral(myreqbin+1, nbins);
    double uqcd_lt     = hqcd->Integral(1,myreqbin);
    double uqcd_tot    = hqcd->Integral();
    
    double corr_scfac = ( data_gt - phojet_gt*data_tot  )/( uqcd_gt - phojet_gt*uqcd_tot  ); 
    
   ///scaled QCD
   double sqcd_lt = corr_scfac*uqcd_lt;

   double sigperc_lt = (data_lt-sqcd_lt)/data_lt;
   double qcdperc_lt = (sqcd_lt)/data_lt;

   //push results into vector
    if(data_lt >0.0){
     fractionQCDSam.push_back(qcdperc_lt);
     fractionQCDErrSam.push_back(0.01);
    }
    else{ 
     fractionQCDSam.push_back(0.0);
     fractionQCDErrSam.push_back(0.0);
    }

  //avoiding 0 entries in the histograms
  cout<<"-----------------------------------------------------"<<endl;                       
  cout<<"                    Running ROOFIT                   "<<endl;
  cout<<"-----------------------------------------------------"<<endl;

  int ndataentries = hdata->GetEntries();
  float sininmin = 0.000; 
  float sininmax = 0.025;
    
 //Set some value not zero for roofit 
  for(int bincount = 1; bincount <= hqcd->GetNbinsX();bincount++){
   if(hqcd->GetBinContent(bincount) == 0.) hqcd->SetBinContent(bincount,1.e-04);
  }

  for(int bincount = 1; bincount <= hphojet->GetNbinsX();bincount++){
   if(hphojet->GetBinContent(bincount) == 0.) hphojet->SetBinContent(bincount,1.e-04);
  }
                  
  //set variable       
  RooRealVar sinin("sinin","sigieie Title",sininmin,sininmax);
    // why have two ranges here? - one for plotting, one for fitting?
    // does sininmax need to match the range of the input histogram?
  sinin.setRange("sigrange",0.005,0.0102);
                       
  //set histograms pdfs
  RooDataHist faketemplate("faketemplate","fake template",sinin,hqcd);
  RooHistPdf fakepdf("fakepdf","test hist fake pdf",sinin,faketemplate);
                       
  RooDataHist realtemplate("realtemplate","real template",sinin,hphojet);
  RooHistPdf realpdf("realpdf","test hist real pdf",sinin,realtemplate);
                       
  //set data distribution to be fitted to
  RooDataHist data("data","data to be fitted to",sinin,hdata);
                       
  //variables that will contain real and fake estimates
  RooRealVar signum("signum","signum",0,ndataentries);
  RooRealVar fakenum("fakenum","fakenum",0,ndataentries);
                       
  //set extended pdfs  
  RooExtendPdf extpdfsig("Signal","extpdfsig",realpdf,signum,"sigrange");
  RooExtendPdf extpdffake("Background","extpdffake",fakepdf,fakenum,"sigrange");
                       
  //composite model pdf
  RooAddPdf model("model","sig + background",RooArgList(extpdfsig,extpdffake));
                       
  //make fit          
   model.fitTo(data,RooFit::Minos());
  // model.fitTo(data,RooFit::Hesse());


  //plot              
  RooPlot *xframe = sinin.frame();
  xframe->SetTitle("");
  xframe->SetXTitle("#sigma i#eta i#eta");
  xframe->GetXaxis()->CenterTitle(1);
  xframe->SetYTitle("Events");
  xframe->GetYaxis()->CenterTitle(1);
  xframe->SetMaximum(40000.);

  data.plotOn(xframe);
  model.plotOn(xframe);
  model.plotOn(xframe,Components(extpdfsig),LineColor(2),LineStyle(kDashed));
  model.plotOn(xframe,Components(extpdffake),LineColor(8),LineStyle(kDashed));
  xframe->SetMinimum(9.);
  xframe->Draw();     
  xframe->Print();                 
  canvas->cd();
  //gPad
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  //Legend
  TLegend *leg1 = new TLegend(0.655172,0.699153,0.886494,0.883475 );
   leg1->SetFillColor(kWhite);
   leg1->SetTextSize(0.02);
   leg1->SetHeader("");
   leg1->AddEntry(xframe->findObject("h_data"), "data", "P");
   leg1->AddEntry(xframe->findObject("model_Norm[sinin]"), "Fit", "L");
   leg1->AddEntry(xframe->findObject("model_Norm[sinin]_Comp[Signal]"), "signal (MC)", "L");
   leg1->AddEntry(xframe->findObject("model_Norm[sinin]_Comp[Background]"), "fake (data)", "L");
   leg1->Draw("same");

   canvas->SaveAs(outpath+"/Fitted_"+datahistname+".C");
   canvas->SaveAs(outpath+"/Fitted_"+datahistname+".pdf");
   canvas->SaveAs(outpath+"/Fitted_"+datahistname+".eps");
 
  //get estimates and their errors
  float fakevalue = fakenum.getValV();
  float fakeerrorhi = fakenum.getErrorHi();
  float fakeerrorlo = fakenum.getErrorLo();
  float fakeerrormax = max(fabs(fakeerrorhi),fabs(fakeerrorlo));
  //signal fraction                    
  float sigvalue = signum.getValV();
  float sigerrorhi = signum.getErrorHi();
  float sigerrorlo = signum.getErrorLo();
  float sigerrormax = max(fabs(sigerrorhi),fabs(sigerrorlo));
                      
  //find the bin corresponding to 0.011
  int binnr = 44;     
  float numerator = hdata->Integral(0,binnr);

      //get results in a vector
      fractionQCD.push_back(fakevalue/(sigvalue+fakevalue));
      fractionQCDErr.push_back( (fakevalue/(sigvalue+fakevalue))* sqrt( ((fakeerrormax/fakevalue)*(fakeerrormax/fakevalue)) + ((sigerrormax/sigvalue)*(sigerrormax/sigvalue)) ) );


   delete xframe;
  
} //+++++++++ Loop over xvariable.list lines +++++++++++++++++++++++++++



 //-----another histos with several bins added up
    lower[0]=175.0;
    lower[1]=190.0;
    lower[2]=250.0;
    lower[3]=400.0;
    lower[4]=700.0;
    lower[5]=1000.0;

   int nqfbins=5;
    
  //PLOT RooFit fractions
    hqcdfraction = new TH1D("hqcdfraction"," QCD Fraction in Num. (RooFit)",(nqfbins),lower);
   for(int ibin=1; ibin<=5; ibin++)
    {           
        hqcdfraction->SetBinContent(ibin,fractionQCD[ibin-1]);
        hqcdfraction->SetBinError(ibin,fractionQCDErr[ibin-1]);
    } 


  //PLOT Sam fractions   
    hqcdfractionSam = new TH1D("hqcdfractionSam"," QCD Fraction in Num. (Sam) ",(nqfbins),lower);
    for(int ibin=1; ibin<=5; ibin++)
    {           
        hqcdfractionSam->SetBinContent(ibin,fractionQCDSam[ibin-1]);
        hqcdfractionSam->SetBinError(ibin,fractionQCDErrSam[ibin-1]);
     } 
   
 


   outfile->cd(); 
  hqcdfractionSam->Write();
  hqcdfraction->Write();
  return;

 
}




//----------------------------------------------
//   Get the Corrected Fake Ratio
//--------------------------------------------
void plot::getCorrectedFakeRatio(TFile* datafile,  //<----data file
                                  TFile* outfile,   //<---output file
                                    vector<double> corrFactor,  //<--input corr. factor
                                      vector<double> corrErr     //<--input error
                                       ){



  //define histo and root files  
  TH1D  *hdatadeno;
  TH1D  *hfulldeno;
  TH1D  *hdijetdeno;
  TH1D  *hphojetdeno;
  TH1D  *hdiphodeno;
  TH1D  *hzjetdeno;

  TH1D  *hdatanum;
  TH1D  *hfullnum;
  TH1D  *hdijetnum;
  TH1D  *hphojetnum;
  TH1D  *hdiphonum;
  TH1D  *hzjetnum;

  TH1D  *hdataratio;
  TH1D  *hfullratio;
  TH1D  *hdijetratio;
  TH1D  *hphojetratio;
  TH1D  *hdiphoratio;
  TH1D  *hzjetratio;


  string shistname;
  TString histname;
  TString datahistnamedeno;
  TString datahistnamenum;

  //get the variable to be plotted: pTEB or ptEE
  char cxvariable[200];
  ifstream infile_xvar;

  infile_xvar.open("xvariables_pt.list", ifstream::in );
  while(!infile_xvar.eof()){
    infile_xvar >>cxvariable;
    string xvariable(cxvariable);

    if(strncmp(cxvariable,"#",1)==0)
      {	continue; }

    histname  = xvariable;
    shistname = histname; 
    datahistnamedeno    = "data_"+shistname+"_deno";
    datahistnamenum    = "data_"+shistname+"_num";
  
   // c1[xvariable] = new TCanvas(("c"+xvariable).c_str(), "transparent pad",179,30,698,498);
    gStyle->SetOptStat(1);
    gStyle->SetPadGridX(3);
    gStyle->SetPadGridY(3);
    gStyle->SetGridStyle(3);

    //get the numerator and denominator
    datafile->cd();
    //deno
    hdatadeno = (TH1D*)datafile->Get(datahistnamedeno)->Clone();
    //num
    hdatanum = (TH1D*)datafile->Get(datahistnamenum)->Clone();
    
    
    //another histos with several bins added up
    double lower[6];

    lower[0]=175.0;
    lower[1]=190.0;
    lower[2]=250.0;
    lower[3]=400.0;
    lower[4]=700.0;
    lower[5]=1000.0;

    int nbins   = hdatanum->GetNbinsX();


  //Lets first print something
   cout<<""<<endl;
   cout<<"Total Fractions are = "<<corrFactor.size()<<" and should be = 7"<<endl;
   cout<<""<<endl;            


  for(int i=0;i<corrFactor.size(); i++)
  {           
    cout<<lower[i]<<"-"<<lower[i+1]<<"(GeV): "<<corrFactor[i]<<" +/- "<<corrErr[i]<<endl;
   }          

    //also add the entries for denonoweight 
    TH1D *hdatanumnw = (TH1D*)datafile->Get(("data_"+shistname+"noweight_num").c_str())->Clone();
    TH1D *hdatadenonw = (TH1D*)datafile->Get(("data_"+shistname+"noweight_deno").c_str())->Clone();
    
    TH1D *hdataratio;
    hdataratio = new TH1D("hdataratio","ratio : data",(nbins),lower);

    TH1D *hdataratioUncorrected;
    hdataratioUncorrected = new TH1D("hdataratioUncorrected","ratio : data uncorrected",(nbins),lower);

    TH1D *hdatanumCorrected;
    hdatanumCorrected = new TH1D("hdatanumCorrected","Numerator corrected",(nbins),lower);

    //initialize the two histo
    for(int ibin=1; ibin<=5; ibin++)
       {   hdataratio->SetBinContent(ibin,0.);
           hdatanumCorrected->SetBinContent(ibin,0.);
        }


   cout<<" "<<endl;
   cout<<"corrected fake rate"<<endl;
   cout<<" "<<endl;

    //Get the final error and fill histo
    for(int ibin=1; ibin<=5; ibin++)
      {
        //cout<<"Numerator without correction = "<<(hdatanum->GetBinContent(ibin))<<" , fraction = "<<corrFactor[ibin-1]<<endl;
	double iinteg_num = hdatanum->GetBinContent(ibin)*(corrFactor[ibin-1]);
	double iinteg_deno = hdatadeno->GetBinContent(ibin);
	
	double eff;
	if(iinteg_deno!=0)eff = iinteg_num/(iinteg_deno);
	if(iinteg_deno==0)eff = 0;
	
	double binerror = 0;
	double integdeno = 0;
	integdeno   = hdatadenonw->GetBinContent(ibin);//for error calcualtion
	double integnum = hdatanumnw->GetBinContent(ibin);

	/////NOT BINOMIAL BUT BY PROPAGATION OF ERRORS////////
	double firstterm ;
	double secterm;
	double thterm  ;
	
	if(iinteg_deno!=0)
	  {  //when template is there
	    firstterm = ( pow(corrFactor[ibin-1],2)/pow(integdeno,2) )*integnum;
	    secterm = pow((integnum*corrFactor[ibin-1]),2)/pow(integdeno,3);
	    thterm = ( pow(integnum,2)/pow(integdeno,2) )*(pow(corrErr[ibin-1],2));
	    binerror = pow(( firstterm + secterm + thterm ),0.5);
	   }
	
	if(iinteg_deno==0)binerror = 0;
          //fake ratio
          hdataratio->SetBinContent(ibin,eff);
	  hdataratio->SetBinError(ibin,(eff*binerror));
          //Fill corrrected fake ratio	
          hdatanumCorrected->SetBinContent(ibin, iinteg_num);
          hdatanumCorrected->SetBinContent(ibin, (hdatanum->GetBinError(ibin)) ); 
	
          cout<<lower[ibin-1]<<"-"<<lower[ibin]<<" (GeV): Fake Rate = "<<iinteg_num/(hdatadeno->GetBinContent(ibin))<<"+/-"<<binerror<<endl;
        }//get final error and fill histo



   cout<<" "<<endl;
   cout<<"Uncorrected fake rate:"<<endl;
   cout<<" "<<endl;

    //Get the final error and fill histo  
    for(int ibin=1; ibin<=5; ibin++)
      {  
        //JUST SET THE correction factor 1.0 so no correction
         corrFactor[ibin-1]=1.0;
         corrErr[ibin-1]=0.0;
        //cout<<"Numerator without correction = "<<(hdatanum->GetBinContent(ibin))<<" , fraction = "<<corrFactor[ibin-1]<<endl;
        double iinteg_num = hdatanum->GetBinContent(ibin)*(corrFactor[ibin-1]);
        double iinteg_deno = hdatadeno->GetBinContent(ibin);                   
         
        double eff;
        if(iinteg_deno!=0)eff = iinteg_num/(iinteg_deno);                      
        if(iinteg_deno==0)eff = 0;
         
        double binerror = 0;
        double integdeno = 0;
        integdeno   = hdatadenonw->GetBinContent(ibin);//for error calcualtion 
        double integnum = hdatanumnw->GetBinContent(ibin);                     
         
        /////NOT BINOMIAL BUT BY PROPAGATION OF ERRORS////////                 
        double firstterm ;
        double secterm;
        double thterm;
         
        if(iinteg_deno!=0)
            {  //when template is there
              firstterm = ( pow(corrFactor[ibin-1],2)/pow(integdeno,2) )*integnum;
              secterm = pow((integnum*corrFactor[ibin-1]),2)/pow(integdeno,3);   
              thterm = ( pow(integnum,2)/pow(integdeno,2) )*(pow(corrErr[ibin-1],2));
              binerror = pow(( firstterm + secterm + thterm ),0.5);              
            }
         
           if(iinteg_deno==0)binerror = 0;

          //fake ratio
          hdataratioUncorrected->SetBinContent(ibin,eff);                                 
          hdataratioUncorrected->SetBinError(ibin,binerror);                              
          cout<<lower[ibin-1]<<"-"<<lower[ibin]<<" (GeV): Fake Rate = "<<iinteg_num/(hdatadeno->GetBinContent(ibin))<<"+/-"<<binerror<<endl;    
     
        }//get final error and fill histo           


    hdataratio->SetLineColor(2);
    hdataratio->SetMarkerStyle(20);
    hdataratio->Draw("E0");


    hdataratioUncorrected->SetLineColor(2);
    hdataratioUncorrected->SetMarkerStyle(20);
    hdataratioUncorrected->Draw("E0");


   outfile->cd();
   hdataratio->Write();
   hdataratioUncorrected->Write();
   hdatanumCorrected->Write();  
 
  }

    return;

}

