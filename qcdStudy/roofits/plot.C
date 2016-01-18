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
 
 TString extraname="";
 //TString extraname="_clossure";
 
 sysnames.push_back("");
 sysnames.push_back("_sbUp");
 sysnames.push_back("_sbDown");
 sysnames.push_back("_metUp");
 sysnames.push_back("_metDown");
 sysnames.push_back("_binUp");
 sysnames.push_back("_binDown");
 sysnames.push_back("_noPiso");

 ptbinbounds.push_back(175);
 ptbinbounds.push_back(190);
 ptbinbounds.push_back(250);
 ptbinbounds.push_back(400);
 ptbinbounds.push_back(1000);

 lower[0]=175.0;
 lower[1]=190.0;
 lower[2]=250.0;
 lower[3]=400.0;
 lower[4]=1000.0;

 nptbins = ptbinbounds.size() - 1;
 //std::cout<<"nptbins:  "<<nptbins<<std::endl;
 
 //std::vector<TString> ptbinnames;
 for(unsigned int i=0; i<nptbins; ++i){
  ptbinnames.push_back(TString(boost::lexical_cast<string>( boost::format("%ito%i") % ptbinbounds[i] % ptbinbounds[i+1] )));
 }
 
 inpath = TString(Tsubmitbase+"/"+Tversion+"/analyzed");
 outpath = TString(Tsubmitbase+"/"+Tversion+"/plots");
 wwwpath = TString("/afs/hep.wisc.edu/user/tperry/www/MonoPhoton/qcdPlots/"+Tversion);
 
 datafile = new TFile(inpath+"/mrg4bins_DataSP2015D.root","READ");
 mcfile   = new TFile(inpath+"/mrg4bins_GJets.root","READ");
 qcdfile  = new TFile(inpath+"/mrg4bins_QCD.root","READ");
 outfile  = new TFile(outpath+"/Fitted.root","RECREATE");

 //ofstream log;
 log.open (outpath+"/log.txt");
 
 //Vector with QCD fractions
 std::vector<double> qcdFrac;
 std::vector<double> qcdFracErr;
 
 //second method without phojet MC
 std::vector<double> qcdFracSam;
 std::vector<double> qcdFracErrSam;

 std::vector<Double_t> bincenters;
 std::vector<Double_t> binerrors;

 std::vector<Double_t> ratios;
 std::vector<Double_t> ratioerrors;
 
 // Do fitting (fill qcdFrac* vectors)
 for(unsigned int sn=0; sn<sysnames.size(); ++sn){
  qcdFrac.clear();
  qcdFracErr.clear();
  qcdFracSam.clear();
  qcdFracErrSam.clear();
  bincenters.clear();
  binerrors.clear();
  ratios.clear();
  ratioerrors.clear();

  sysname = sysnames[sn];

  std::cout<<endl<<endl<<endl<<"Starting Systematic:   "<<sysname<<std::endl;
  //log<<boost::format("\n\nStarting Systematic:  %s \n") % sysname;

  getFraction(
   datafile,
   mcfile,
   qcdfile,
   outfile,
   qcdFrac,
   qcdFracErr,
   qcdFracSam,
   qcdFracErrSam,
   sysname,
   extraname
  );

  qcd_frac.push_back(qcdFrac);
  qcd_frac_err.push_back(qcdFracErr);

  getBinCenters(datafile, bincenters, binerrors, sysname);

  bincenterss.push_back(bincenters);
  binerrorss.push_back(binerrors);
  
  //Using RooFit result
  getCorrectedFakeRatio(
   datafile,
   outfile,
   qcdFrac,
   qcdFracErr,
   ratios,
   ratioerrors,
   sn
   //sysname
  );

 ratioss.push_back(ratios);
 ratioerrorss.push_back(ratioerrors);

 }

 for(unsigned int sn=0; sn<sysnames.size(); ++sn){
   log<<boost::format("\nSysname: %s \n") % sysnames[sn];
  for(unsigned int j=0; j<nptbins; ++j){
    log<<boost::format(" %s\n") % ptbinnames[j];
    //log<<boost::format("  QCD Fraction: %0.3d +- %0.3d \n") % qcd_frac[sn][j] % qcd_frac_err[sn][j];
    log<<boost::format("  Bin Center: %0.3d +- %0.3d \n") % bincenterss[sn][j] % binerrorss[sn][j];
    log<<boost::format("  Ratio: %0.3d +- %0.3d \n") % ratioss[sn][j] % ratioerrorss[sn][j];
  }
    log<<boost::format(" y = (%0.3d +- %0.3d)x + (%0.3d +- %0.3d) \n") % ms[sn] % mes[sn] % bs[sn] % bes[sn] ;
 }

 drawAllRates();

 log.close();

}


//------------------------------------
//  get fraction fits from templates
//------------------------------------
void plot::getFraction(
 TFile* datafile,  
 TFile* mcfile,    
 TFile* qcdfile,    
 TFile* outfile,   
 std::vector<double>& fractionQCD,      //<<---This gets filled here
 std::vector<double>& fractionQCDErr,   //<<---This gets filled here
 std::vector<double>& fractionQCDSam,   //<<---This gets filled here
 std::vector<double>& fractionQCDErrSam,//<<---This gets filled here
 TString sysname,
 TString extraname
 )
{

  TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500);   
  gStyle->SetOptStat(0);
  //gStyle->SetPadGridX(3);
  //gStyle->SetPadGridY(3);
  //gStyle->SetGridStyle(3);
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  TH1D* hdata;
  TH1D* htemp;
  TH1D* hqcd;
  TH1D* hphojet;
  TH1D* hqcdfractionSam;
  TH1D* hqcdfraction;

  TString datahistname; //data histo
  TString qcdhistname;  //trkflip - > histos will go here
  TString fullhistname; //phojet

  TString binrange;
  for(unsigned int ptb=0; ptb<ptbinnames.size(); ++ptb){
   binrange = ptbinnames[ptb];
   //log<<boost::format(" %s \n") % binrange;
  
    //set template histo names
    datahistname = "h_sig_sieieF5x5_"+binrange+sysname; //data histos with tigh ID and without sigIetaIeta cut
    fullhistname = "h_sig_sieieF5x5_"+binrange+sysname; //phojet histos with tight ID cut but without sigIetaIeta
    qcdhistname  = "h_bkg_sieieF5x5_"+binrange+sysname; //data histos with very loose id and  sideband of track iso


    // for closure test
    //get "Data" template = MC (GJ + QCD) 
    hdata = (TH1D*)mcfile->Get(datahistname)->Clone();
    htemp = (TH1D*)qcdfile->Get(qcdhistname)->Clone();
    hdata->Add(htemp);
    //QCD histo
    hqcd  = (TH1D*)qcdfile->Get(qcdhistname)->Clone();

    // //  used and good
    // //get Data template and QCD template from data
    // datafile->cd();
    // hdata = (TH1D*)datafile->Get(datahistname)->Clone();
    // //QCD histo
    // hqcd  = (TH1D*)datafile->Get(qcdhistname)->Clone();

    double integ_data = hdata->Integral();
    double integ_qcd = hqcd->Integral();
    //get photon+jet signal tempalte
    mcfile->cd();
    hphojet  = (TH1D*)mcfile->Get(fullhistname)->Clone();
    double integ_phojet = hphojet->Integral();
    int nbins = hphojet->GetNbinsX();

    canvas->cd();
    hphojet->Draw();
    canvas->Print(outpath+"/hphojet"+binrange+sysname+".png");
    hdata->Draw();
    canvas->Print(outpath+"/hdata"+binrange+sysname+".png");
    hqcd->Draw();
    canvas->Print(outpath+"/hqcd"+binrange+sysname+".png");


  //------Lets try to get with method-2 without MC
  cout<<"-----------------------------------------------------"<<endl;                       
  cout<<"                    Running Sam's Method              "<<endl;
  cout<<"-----------------------------------------------------"<<endl;
  //log<<boost::format(" Sam's Method \n");

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
  cout<<" corr_scfac = "<<corr_scfac<<endl;
  //log<<boost::format("  corr_scfac:  %.3f \n") % corr_scfac;
  
  ///scaled QCD
  double sqcd_lt = corr_scfac*uqcd_lt;
  cout<<" sqcd_lt = "<<sqcd_lt<<endl;
  //log<<boost::format("  sqcd_lt:  %.3f \n") % sqcd_lt;

  double sigperc_lt = (data_lt-sqcd_lt)/data_lt;
  cout<<" sigperc_lt = "<<sigperc_lt<<endl;
  //log<<boost::format("  sigperc_lt:  %.3f \n") % sigperc_lt;
  double qcdperc_lt = (sqcd_lt)/data_lt;
  cout<<" qcdperc_lt = "<<qcdperc_lt<<endl;
  //log<<boost::format("  qcdperc_lt:  %.3f \n") % qcdperc_lt;

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

  int ndataentries = hdata->Integral(); //hdata->GetEntries();
  float sininmin = 0.000; 
  float sininmax = 0.025;
  float sdataentries = hdata->Integral();
    
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

  // Do some calculations
  cout<<binrange<<endl;
 
  //get estimates and their errors
  float fakevalue = fakenum.getValV();
  float fakeerrorhi = fakenum.getErrorHi();
  float fakeerrorlo = fakenum.getErrorLo();
  float fakeerrormax = max(fabs(fakeerrorhi),fabs(fakeerrorlo));

  cout<<" fakevalue   : "<<fakevalue<<endl;    
  cout<<" fakeerrorhi : "<<fakeerrorhi<<endl; 
  cout<<" fakeerrorlo : "<<fakeerrorlo<<endl; 
  cout<<" fakeerrormax: "<<fakeerrormax<<endl;

  //signal fraction                    
  float sigvalue = signum.getValV();
  float sigerrorhi = signum.getErrorHi();
  float sigerrorlo = signum.getErrorLo();
  float sigerrormax = max(fabs(sigerrorhi),fabs(sigerrorlo));

  cout<<" sigvalue   : "<<sigvalue<<endl;    
  cout<<" sigerrorhi : "<<sigerrorhi<<endl; 
  cout<<" sigerrorlo : "<<sigerrorlo<<endl; 
  cout<<" sigerrormax: "<<sigerrormax<<endl;
                      
  float frQCD = fakevalue/(sigvalue+fakevalue);
  float frQCDerr = ( fakevalue/(sigvalue+fakevalue)) *
                   sqrt(
                    ( (fakeerrormax/fakevalue)*(fakeerrormax/fakevalue) )
                    + ( (sigerrormax/sigvalue)*(sigerrormax/sigvalue) )
                   ) ;

  cout<<" frQCD    : "<<frQCD<<endl;
  cout<<" frQCDerr : "<<frQCDerr<<endl;

  //put results into vector
  fractionQCD.push_back(frQCD);
  fractionQCDErr.push_back(frQCDerr);

  TString sndataentries = TString(boost::lexical_cast<string>( boost::format("%.1f") % ndataentries ));
  TString sdata_lt      = TString(boost::lexical_cast<string>( boost::format("%.1f") % sdata_lt ));
  TString sfakevalue    = TString(boost::lexical_cast<string>( boost::format("%.1f") % fakevalue    )); 
  TString sfakeerrormax = TString(boost::lexical_cast<string>( boost::format("%.1f") % fakeerrormax )); 
  TString ssigvalue     = TString(boost::lexical_cast<string>( boost::format("%.1f") % sigvalue     )); 
  TString ssigerrormax  = TString(boost::lexical_cast<string>( boost::format("%.1f") % sigerrormax  )); 
  TString sfitvalue     = TString(boost::lexical_cast<string>( boost::format("%.1f") % sqrt(sigvalue   *sigvalue   +fakevalue   *fakevalue   )  )); 
  TString sfiterrormax  = TString(boost::lexical_cast<string>( boost::format("%.1f") % sqrt(sigerrormax*sigerrormax+fakeerrormax*fakeerrormax)  )); 

  TString sfrQCD     = TString(boost::lexical_cast<string>( boost::format("%.3f") % frQCD    )); 
  TString sfrQCDerr  = TString(boost::lexical_cast<string>( boost::format("%.3f") % frQCDerr )); 


  //plot              
  RooPlot *xframe = sinin.frame();
  xframe->SetTitle("");
  xframe->SetXTitle("#sigma i#eta i#eta");
  xframe->GetXaxis()->CenterTitle(1);
  xframe->SetYTitle("Events");
  xframe->GetYaxis()->CenterTitle(1);
  xframe->SetMaximum(400000.);

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
  gStyle->SetLineWidth(3);
  //Legend
  TLegend *leg1 = new TLegend(0.55,0.6,0.88,0.88 );
  //TLegend *leg1 = new TLegend(0.655172,0.699153,0.886494,0.883475 );
   leg1->SetFillColor(kWhite);
   //leg1->SetTextSize(0.02);
   //leg1->SetHeader("");
   leg1->AddEntry(xframe->findObject("h_data"), "data", "P");
   //leg1->AddEntry(xframe->findObject("h_data"), "data ("+sdata_lt+")", "P");
   //leg1->AddEntry(xframe->findObject("h_data"), "data ("+sndataentries+")", "P");
   leg1->AddEntry(xframe->findObject("model_Norm[sinin]"), "Fit", "L");
   //leg1->AddEntry(xframe->findObject("model_Norm[sinin]"), 
   // sfitvalue+" +- "+sfiterrormax, "L");
   leg1->AddEntry(xframe->findObject("model_Norm[sinin]_Comp[Signal]"), "signal (GJ MC) ", "L");
   //leg1->AddEntry(xframe->findObject("model_Norm[sinin]_Comp[Signal]"),
   // ssigvalue+" +- "+ssigerrormax, "L");
   //leg1->AddEntry(xframe->findObject("model_Norm[sinin]_Comp[Background]"), "fake (QCD MC)", "L");
   leg1->AddEntry(xframe->findObject("model_Norm[sinin]_Comp[Background]"), "fake (data)", "L");
   //leg1->AddEntry(xframe->findObject("model_Norm[sinin]_Comp[Background]"),
   // sfakevalue+" +- "+sfakeerrormax, "L");
   leg1->Draw("same");

   // Title
   TText* title = new TText(1,1,"") ;
   title->SetTextSize(0.07);
   title->SetTextColor(kBlack);
   title->SetTextAlign(13);
   title->SetTextFont(62);
   title->DrawTextNDC(0.17,0.87,"CMS");
   xframe->addObject(title);

   TText* extra = new TText(1,1,"") ;
   extra->SetTextSize(0.05);
   extra->SetTextColor(kBlack);
   extra->SetTextAlign(13);
   extra->SetTextFont(52);
   extra->DrawTextNDC(0.17,0.81,"Preliminary");
   xframe->addObject(extra);

   TText* lumi = new TText(1,1,"") ;
   lumi->SetTextSize(0.05);
   lumi->SetTextColor(kBlack);
   lumi->SetTextAlign(31);
   lumi->SetTextFont(42);
   lumi->DrawTextNDC(0.9,0.91,"2.1 /fb (13 TeV)");
   xframe->addObject(lumi);

   TText* ptrange = new TText(1,1,"") ;
   ptrange->SetTextSize(0.03);
   ptrange->SetTextColor(kBlack);
   ptrange->SetTextAlign(11);
   ptrange->SetTextFont(42);
   ptrange->DrawTextNDC(0.13,0.73,"pT Range [GeV]: "+binrange);
   xframe->addObject(ptrange);

   TText* sflabel = new TText(1,1,"") ;
   sflabel->SetTextSize(0.07);
   sflabel->SetTextColor(kBlack);
   sflabel->SetTextAlign(11);
   sflabel->SetTextFont(42);
   sflabel->DrawTextNDC(0.450,0.45,"B/(S+B)= ");
   sflabel->DrawTextNDC(0.450,0.38,sfrQCD+" +- "+sfrQCDerr);
   xframe->addObject(sflabel);

   canvas->SaveAs(outpath+"/Fitted_"+datahistname+sysname+extraname+".pdf");
   //canvas->SaveAs(outpath+"/Fitted_"+datahistname+sysname+extraname+".C");
   //canvas->SaveAs(outpath+"/Fitted_"+datahistname+sysname+extraname+".eps");
   //canvas->SaveAs(wwwpath+"/Fitted_"+datahistname+extraname+".png");

  delete xframe;
  
} // Loop over pt bins

  //PLOT RooFit fractions
   hqcdfraction = new TH1D("hqcdfraction"," QCD Fraction in Num. (RooFit)",(nptbins),lower);
   for(int ibin=1; ibin<=nptbins+1; ibin++)
    {           
        hqcdfraction->SetBinContent(ibin,fractionQCD[ibin-1]);
        hqcdfraction->SetBinError(ibin,fractionQCDErr[ibin-1]);
    } 

  //PLOT Sam fractions   
    hqcdfractionSam = new TH1D("hqcdfractionSam"," QCD Fraction in Num. (Sam) ",(nptbins),lower);
    for(int ibin=1; ibin<=nptbins+1; ibin++)
    {           
        hqcdfractionSam->SetBinContent(ibin,fractionQCDSam[ibin-1]);
        hqcdfractionSam->SetBinError(ibin,fractionQCDErrSam[ibin-1]);
     } 

  outfile->cd(); 
  hqcdfractionSam->Write();
  hqcdfraction->Write();
  return;
 
}


void plot::getBinCenters(
 TFile* datafile,  
 std::vector<double>& bcenters, //<<---This gets filled here
 std::vector<double>& berrors,  //<<---This gets filled here
 TString sysname
 )
{

  TH1D* hdata;
  TString datahistname; //data histo
  TString binrange;
  for(unsigned int ptb=0; ptb<ptbinnames.size(); ++ptb){
   binrange = ptbinnames[ptb];
 
   //set template histo names
   datahistname = "h_sig_et_"+binrange+sysname; //data histos with tigh ID and without sigIetaIeta cut

   //// for closure test
   ////get "Data" template = MC (GJ + QCD) 
   //hdata = (TH1D*)mcfile->Get(datahistname)->Clone();
   //htemp = (TH1D*)qcdfile->Get(qcdhistname)->Clone();
   //hdata->Add(htemp);

   //  used and good
   //get Data template and QCD template from data
   datafile->cd();
   hdata = (TH1D*)datafile->Get(datahistname)->Clone();

   bcenters.push_back(hdata->GetMean());
   //berrors.push_back(hdata->GetMeanError());
   berrors.push_back(hdata->GetRMS());
  }
}

//----------------------------------------------
//   Get the Corrected Fake Ratio
//--------------------------------------------
void plot::getCorrectedFakeRatio(TFile* datafile,  //<----data file
                                 TFile* outfile,   //<---output file
                                 vector<double> corrFactor, //<--input corr. factor = QCD fraction
                                 vector<double> corrErr,    //<--input error
                                 vector<Double_t>& ratios,
                                 vector<Double_t>& ratioerrors,
                                 int sn
                                 //TString sysname
                                ){



  sysname = sysnames[sn];

  //define histo and root files  
  TH1D  *hdataden;
 // TH1D  *hfulldeno;
 // TH1D  *hdijetdeno;
 // TH1D  *hphojetdeno;
 // TH1D  *hdiphodeno;
 // TH1D  *hzjetdeno;

  TH1D  *hdatanum;
 // TH1D  *hfullnum;
 // TH1D  *hdijetnum;
 // TH1D  *hphojetnum;
 // TH1D  *hdiphonum;
 // TH1D  *hzjetnum;

 // TH1D  *hdataratio;
 // TH1D  *hfullratio;
 // TH1D  *hdijetratio;
 // TH1D  *hphojetratio;
 // TH1D  *hdiphoratio;
 // TH1D  *hzjetratio;

  string shistname;
  TString binrange;
  TString denhistname;
  TString numhistname;

  //bins 
//  //Double_t lower[4];
//  lower[0]=175.0;
//  lower[1]=190.0;
//  lower[2]=250.0;
//  //lower[3]=1000.0;
//  lower[3]=400.0;
//  lower[4]=1000.0;
//
  //Double_t size_num_uncor[3];
  //Double_t err_num_uncor[3];
  //Double_t size_num_corr[3];
  //Double_t err_num_corr[3];
  //Double_t size_den[3];
  //Double_t err_den[3];
  
  Double_t size_num_uncor[nptbins];
  Double_t err_num_uncor[nptbins];
  Double_t size_num_corr[nptbins];
  Double_t err_num_corr[nptbins];
  Double_t size_den[nptbins];
  Double_t err_den[nptbins];
  
  // temp variables to to in vectors above
  double snu,enu,snc,enc,sd,ed;

  //TString binrange;
  for(unsigned int ptb=0; ptb<ptbinnames.size(); ++ptb){
   binrange = ptbinnames[ptb];

  //set template histo names
  numhistname = "h_sig_sieieF5x5_"+binrange+sysname; //data histos with tight ID cut but without sigIetaIeta
  denhistname = "h_den_sieieF5x5_"+binrange+sysname; //data histos fail loose ID and pass Vloose ID + || inv iso

  //get the numerator and denominator
  datafile->cd();
  //num
  hdatanum = (TH1D*)datafile->Get(numhistname)->Clone();    
  //deno
  hdataden = (TH1D*)datafile->Get(denhistname)->Clone();
   
  // bin corresponding to medium wp sieie cut
  double sieiebin;
  sieiebin = hdatanum->FindBin(0.0102);
  size_num_uncor[ptb] = double(hdatanum->Integral(1,sieiebin));
  size_num_corr[ptb] = double((corrFactor[ptb])) * double(hdatanum->Integral(1,sieiebin));
  size_den[ptb] = double(hdataden->Integral(1,sieiebin));

  double tmp_err_num_corr;
  //hdatanum->IntegralAndError(1, sieiebin, err_num_uncor[ptb]);
  hdatanum->IntegralAndError(1, sieiebin, tmp_err_num_corr);
  //err_num_corr[varcounter] = sqrt( err_num_corr[varcounter]*err_num_corr[varcounter] + 
  err_num_corr[ptb] = sqrt( tmp_err_num_corr*tmp_err_num_corr + 
                              corrErr[ptb]*corrErr[ptb] );
  hdataden->IntegralAndError(1, sieiebin, err_den[ptb]);

  }

  TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500);   

  Double_t *mid = new Double_t[nptbins]; //{182.5, 220., 325., 700. };
  Double_t *miderr = new Double_t[nptbins];

  for(unsigned int k=0; k<nptbins; ++k){
   mid[k]=bincenterss[sn][k];
   miderr[k]=binerrorss[sn][k];
   //std::cout<<bincenterss[sn][k]<<std::endl;
  }

  ////std::vector<double> mid;
  //mid[0]=182.5;
  //mid[1]=220.;
  //mid[2]=325.;
  //mid[3]=700.;
  //miderr[0]=mid[0]-lower[0];
  //miderr[1]=mid[1]-lower[1];
  //miderr[2]=mid[2]-lower[2];
  //miderr[3]=mid[3]-lower[3];

  TGraphErrors *gr_num_uncor = new TGraphErrors(4,mid,size_num_uncor,miderr,err_num_uncor);
  TGraphErrors *gr_num_corr  = new TGraphErrors(4,mid,size_num_corr,miderr,err_num_corr);
  TGraphErrors *gr_den       = new TGraphErrors(4,mid,size_den,miderr,err_den);

  // c1[xvariable] = new TCanvas(("c"+xvariable).c_str(), "transparent pad",179,30,698,498);
  gStyle->SetOptStat(0);
  //gStyle->SetPadGridX(3);
  //gStyle->SetPadGridY(3);
  //gStyle->SetGridStyle(3);
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  gr_num_uncor->SetTitle("");
  gr_num_uncor->SetMarkerColor(1);
  gr_num_uncor->SetLineColor(1);
  gr_num_uncor->SetMarkerStyle(21);

  gr_num_corr->SetMarkerColor(2);
  gr_num_corr->SetLineColor(2);
  gr_num_corr->SetMarkerStyle(22);

  gr_den->SetMarkerColor(3);
  gr_den->SetLineColor(3);
  gr_den->SetMarkerStyle(23);

  TH1F *hr = canvas->DrawFrame(175.,10.,1000.,10000000.,"");
  hr->SetXTitle("photon pT [GeV]");
  hr->SetYTitle("Events"); 

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

  TLegend *leg2 = new TLegend(0.55,0.6,0.88,0.88 );
  leg2->SetFillColor(kWhite);
  //leg2->AddEntry( gr_num_uncor,"Numerator", "L");
  leg2->AddEntry( gr_num_uncor,"Numerator (uncorrected)", "L");
  leg2->AddEntry( gr_num_corr,"Numerator (corrected)", "L");
  leg2->AddEntry( gr_den,"Denominator", "L");
  leg2->Draw("same");

  gr_num_uncor->Draw("P");
  gr_num_corr->Draw("P"); //
  gr_den->Draw("P");

  canvas->SaveAs(outpath+"/Graph_NuNcD"+sysname+".pdf");
  //canvas->SaveAs(outpath+"/Graph_NuNcD.C");
  //canvas->SaveAs(outpath+"/Graph_NuNcD.eps");
  //canvas->SaveAs(wwwpath+"/Graph_NuNcD.png");

  canvas->Clear();

  // Ratio Plot
  //////////////////////////////////////////////////////////////////////////////
  TCanvas* c2 = new TCanvas("c2","c2",900,100,500,500);   

  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  //Double_t ratio[3];
  Double_t *ratio = new Double_t[4];
  ratio[0] = size_num_corr[0]/size_den[0];
  ratio[1] = size_num_corr[1]/size_den[1];
  ratio[2] = size_num_corr[2]/size_den[2];
  ratio[3] = size_num_corr[3]/size_den[3];

  Double_t *ratioerr = new Double_t[4];
  ratioerr[0]= ratio[0] * pow( pow(err_num_corr[0]/size_num_corr[0] + pow(err_den[0]/size_den[0] ,2) ,2) ,0.5); //0.005;
  ratioerr[1]= ratio[1] * pow( pow(err_num_corr[1]/size_num_corr[1] + pow(err_den[1]/size_den[1] ,2) ,2) ,0.5); //0.005;
  ratioerr[2]= ratio[2] * pow( pow(err_num_corr[2]/size_num_corr[2] + pow(err_den[2]/size_den[2] ,2) ,2) ,0.5); //0.005;
  ratioerr[3]= ratio[3] * pow( pow(err_num_corr[3]/size_num_corr[3] + pow(err_den[3]/size_den[3] ,2) ,2) ,0.5); //0.005;

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

  //for(int i=0; i<9; i++){std::cout<<"mid:  "<<mid[i]<<std::endl;}
  //for(int i=0; i<4; i++){std::cout<<"ratio:  "<<ratio[i]<<std::endl;}
  //for(int i=0; i<4; i++){std::cout<<"miderr:  "<<miderr[i]<<std::endl;}
  //for(int i=0; i<4; i++){std::cout<<"ratioerr:  "<<ratioerr[i]<<std::endl;}

  //gr_ratio->Set(4);
  //TF1 *fit_ratio = new TF1("fit_ratio","[0]+[1]*x", 175, 1000);
  //fit_ratio->SetParameter(0,0.05);
  //fit_ratio->SetParameter(1,0.);
  //TF1 *fit_ratio = new TF1("fit_ratio","pol1", 175, 1000);
  //TF1 *fit_ratio = new TF1("fit_ratio","pol0(0)", 175, 1000);
  //TF1 *fit_ratio = new TF1("fit_ratio","pol0", 175, 1000);

  gr_ratio->SetMarkerColor(2);
  gr_ratio->SetLineColor(2);
  gr_ratio->SetMarkerStyle(22);

  //TH1F *hs = c2->DrawFrame(0.,0.,1000.,0.5,"");
  TH1F *hs = c2->DrawFrame(175.,0.,1000.,0.5,"");
  hs->SetXTitle("photon pT [GeV]");
  hs->SetYTitle("Fake Ratio"); 

  title->DrawTextNDC(0.17,0.87,"CMS");
  extra->DrawTextNDC(0.17,0.81,"Preliminary");
  lumi->DrawTextNDC(0.9,0.91,"2.1 /fb (13 TeV)");

  TLegend *leg3 = new TLegend(0.55,0.6,0.88,0.88 );
  leg3->SetFillColor(kWhite);
  leg3->AddEntry( gr_ratio,"Num./Den.", "L");
  leg3->Draw("same");

  gr_ratio->Draw("P");
  gr_ratio_noe->Draw("P,sames");
  gr_ratio_noe->Fit("pol1");
  //gr_ratio->Fit("pol0");
  TF1 *fit_ratio = gr_ratio_noe->GetFunction("pol1");
  fit_ratio->SetLineColor(1);
  fit_ratio->SetLineWidth(3);
  Double_t chi2 = fit_ratio->GetChisquare();
  Double_t p0 = fit_ratio->GetParameter(0);
  Double_t p1 = fit_ratio->GetParameter(1);
  Double_t e0 = fit_ratio->GetParError(0);
  Double_t e1 = fit_ratio->GetParError(1);

  ms.push_back(p1); mes.push_back(e1);
  bs.push_back(p0); bes.push_back(e0);

  fit_ratio->Draw("sames");

  TLatex tex;
  tex.SetTextSize(0.07);
  tex.SetTextColor(kBlack);
  tex.SetTextAlign(11);
  tex.SetTextFont(42);
  tex.DrawLatexNDC(0.45,0.47,"y = m x + b");
  tex.SetTextSize(0.05);
  tex.DrawLatexNDC(0.35,0.40,TString(boost::lexical_cast<string>(boost::format("m = %0.5f +- %0.5f") % p1 % e1))); 
  tex.DrawLatexNDC(0.45,0.35,TString(boost::lexical_cast<string>(boost::format("b = %0.3f +- %0.3f") % p0 % e0))); 
  tex.DrawLatexNDC(0.45,0.30,"#chi^{2}"+TString(boost::lexical_cast<string>(boost::format(" = %0.5f") % chi2))); 

  //TText* eqlabel = new TText(1,1,"") ;
  //eqlabel->SetTextSize(0.07);
  //eqlabel->SetTextColor(kBlack);
  //eqlabel->SetTextAlign(11);
  //eqlabel->SetTextFont(42);
  //eqlabel->DrawTextNDC(0.45,0.47,"y = m x + b");
  //eqlabel->SetTextSize(0.05);
  //eqlabel->DrawTextNDC(0.45,0.40,TString(boost::lexical_cast<string>(boost::format("m = %0.3f +- %0.3f") % p1 % e1))); 
  //eqlabel->DrawTextNDC(0.45,0.35,TString(boost::lexical_cast<string>(boost::format("b = %0.3f +- %0.3f") % p0 % e0))); 
  //eqlabel->DrawTextNDC(0.45,0.30,"#chi^{2}"+TString(boost::lexical_cast<string>(boost::format(" = %0.3f") % chi2))); 
  //log<<boost::format("\n\nStarting Systematic:  %s \n") % sysname;

  //std::cout<<"Npoints: "<<gr_ratio->GetN()<<std::endl;
  //Double_t thex;
  //Double_t they;
  //for(unsigned int b=0; b<gr_ratio->GetN(); ++b){
  // gr_ratio->GetPoint(b,thex,they);
  // std::cout<<"Point "<<b<<": "<<thex<<" "<<they<<std::endl;
  //}
  //std::cout<<"Fit : "<<chi2<<" "<<p0<<" "<<p1<<std::endl;

  c2->Update();
  
  //gr_ratio->Draw("sames,A*");

  c2->SaveAs(outpath+"/Graph_FakeRatio"+sysname+".pdf");
  //c2->SaveAs(outpath+"/Graph_FakeRatio.C");
  //c2->SaveAs(outpath+"/Graph_FakeRatio.eps");
  //c2->SaveAs(wwwpath+"/Graph_FakeRatio.png");

  c2->Clear();

  return;

}



void plot::drawAllRates(){

  TCanvas* c3 = new TCanvas("c3","c3",900,100,500,500);   

  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  //gr_ratio->SetMarkerColor(2);
  //gr_ratio->SetLineColor(2);
  //gr_ratio->SetMarkerStyle(22);

  double xmax;
  xmax = 600.;

  Double_t xlow = 175.;
  Double_t xhi = 600.;
  //TH1F *hs = c2->DrawFrame(0.,0.,1000.,0.5,"");
  //TH1F *hs = c3->DrawFrame(175.,0.,1000.,0.5,"");
  TH1F *hs = c3->DrawFrame(xlow,0.,xhi,0.5,"");
  hs->SetXTitle("photon pT [GeV]");
  hs->SetYTitle("Fake Ratio"); 

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

  title->DrawTextNDC(0.17,0.87,"CMS");
  extra->DrawTextNDC(0.17,0.81,"Preliminary");
  lumi->DrawTextNDC(0.9,0.91,"2.1 /fb (13 TeV)");

  ////////// Draw points first, then lines

  // x
  Double_t *mid0 = new Double_t[nptbins];
  Double_t *mid1 = new Double_t[nptbins];
  Double_t *mid2 = new Double_t[nptbins];
  Double_t *mid3 = new Double_t[nptbins];
  Double_t *mid4 = new Double_t[nptbins];
  Double_t *mid5 = new Double_t[nptbins];
  Double_t *mid6 = new Double_t[nptbins];
  Double_t *mid7 = new Double_t[nptbins];
  for(unsigned int k=0; k<nptbins; ++k){ mid0[k]=bincenterss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid1[k]=bincenterss[1][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid2[k]=bincenterss[2][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid3[k]=bincenterss[3][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid4[k]=bincenterss[4][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid5[k]=bincenterss[5][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid6[k]=bincenterss[6][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid7[k]=bincenterss[7][k]; }

  Double_t *miderr0 = new Double_t[nptbins];
  Double_t *miderr1 = new Double_t[nptbins];
  Double_t *miderr2 = new Double_t[nptbins];
  Double_t *miderr3 = new Double_t[nptbins];
  Double_t *miderr4 = new Double_t[nptbins];
  Double_t *miderr5 = new Double_t[nptbins];
  Double_t *miderr6 = new Double_t[nptbins];
  Double_t *miderr7 = new Double_t[nptbins];
  for(unsigned int k=0; k<nptbins; ++k){ miderr0[k]=binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr1[k]=binerrorss[1][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr2[k]=binerrorss[2][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr3[k]=binerrorss[3][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr4[k]=binerrorss[4][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr5[k]=binerrorss[5][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr6[k]=binerrorss[6][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr7[k]=binerrorss[7][k]; }
  
  // y
  Double_t *points0 = new Double_t[4];
  Double_t *points1 = new Double_t[4];
  Double_t *points2 = new Double_t[4];
  Double_t *points3 = new Double_t[4];
  Double_t *points4 = new Double_t[4];
  Double_t *points5 = new Double_t[4];
  Double_t *points6 = new Double_t[4];
  Double_t *points7 = new Double_t[4];
  for(unsigned int p=0; p<4; ++p){ points0[p]=ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points1[p]=ratioss[1][p]; }
  for(unsigned int p=0; p<4; ++p){ points2[p]=ratioss[2][p]; }
  for(unsigned int p=0; p<4; ++p){ points3[p]=ratioss[3][p]; }
  for(unsigned int p=0; p<4; ++p){ points4[p]=ratioss[4][p]; }
  for(unsigned int p=0; p<4; ++p){ points5[p]=ratioss[5][p]; }
  for(unsigned int p=0; p<4; ++p){ points6[p]=ratioss[6][p]; }
  for(unsigned int p=0; p<4; ++p){ points7[p]=ratioss[7][p]; }

  Double_t *pointerrors0 = new Double_t[4];
  Double_t *pointerrors1 = new Double_t[4];
  Double_t *pointerrors2 = new Double_t[4];
  Double_t *pointerrors3 = new Double_t[4];
  Double_t *pointerrors4 = new Double_t[4];
  Double_t *pointerrors5 = new Double_t[4];
  Double_t *pointerrors6 = new Double_t[4];
  Double_t *pointerrors7 = new Double_t[4];
  for(unsigned int p=0; p<4; ++p){ pointerrors0[p]=ratioerrorss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors1[p]=ratioerrorss[1][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors2[p]=ratioerrorss[2][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors3[p]=ratioerrorss[3][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors4[p]=ratioerrorss[4][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors5[p]=ratioerrorss[5][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors6[p]=ratioerrorss[6][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors7[p]=ratioerrorss[7][p]; }

  // fill the graphs
  TGraphErrors *gr_points0 = new TGraphErrors(4,mid0,points0,miderr0,pointerrors0);
  TGraphErrors *gr_points1 = new TGraphErrors(4,mid1,points1,miderr1,pointerrors1);
  TGraphErrors *gr_points2 = new TGraphErrors(4,mid2,points2,miderr2,pointerrors2);
  TGraphErrors *gr_points3 = new TGraphErrors(4,mid3,points3,miderr3,pointerrors3);
  TGraphErrors *gr_points4 = new TGraphErrors(4,mid4,points4,miderr4,pointerrors4);
  TGraphErrors *gr_points5 = new TGraphErrors(4,mid5,points5,miderr5,pointerrors5);
  TGraphErrors *gr_points6 = new TGraphErrors(4,mid6,points6,miderr6,pointerrors6);
  TGraphErrors *gr_points7 = new TGraphErrors(4,mid7,points7,miderr7,pointerrors7);

  // set colors
  TColor *color0 = gROOT->GetColor(kBlack);    
  TColor *color1 = gROOT->GetColor(kRed+1);
  TColor *color2 = gROOT->GetColor(kYellow-3);
  TColor *color3 = gROOT->GetColor(kGreen+1);
  TColor *color4 = gROOT->GetColor(kGreen-9);
  TColor *color5 = gROOT->GetColor(kAzure+10);
  TColor *color6 = gROOT->GetColor(kBlue-9);
  TColor *color7 = gROOT->GetColor(51);

  //Color_t color0 = gROOT->GetColor(kBlack);
  //Color_t color1 = gROOT->GetColor(kRed+1);
  //Color_t color2 = gROOT->GetColor(kYellow-3);
  //Color_t color3 = gROOT->GetColor(kGreen+1);
  //Color_t color4 = gROOT->GetColor(kGreen-9);
  //Color_t color5 = gROOT->GetColor(kAzure+10);
  //Color_t color6 = gROOT->GetColor(kBlue-9);
  //Color_t color7 = gROOT->GetColor(51);

  gr_points0->SetMarkerColor(kBlack);    //         Color(color0);
  gr_points1->SetMarkerColor(kRed+1);    //         Color(color1);
  gr_points2->SetMarkerColor(kYellow-3); //         Color(color2);
  gr_points3->SetMarkerColor(kGreen+1);  //         Color(color3);
  gr_points4->SetMarkerColor(kGreen-9);  //         Color(color4);
  gr_points5->SetMarkerColor(kAzure+10); //         Color(color5);
  gr_points6->SetMarkerColor(kBlue-9);   //         Color(color6);
  gr_points7->SetMarkerColor(51);        //         Color(color7);

  gr_points0->SetLineWidth(2); //         Color(color0);
  gr_points1->SetLineWidth(2); //         Color(color1);
  gr_points2->SetLineWidth(2); //         Color(color2);
  gr_points3->SetLineWidth(2); //         Color(color3);
  gr_points4->SetLineWidth(2); //         Color(color4);
  gr_points5->SetLineWidth(2); //         Color(color5);
  gr_points6->SetLineWidth(2); //         Color(color6);
  gr_points7->SetLineWidth(2); //         Color(color7);

  gr_points0->SetLineColor(kBlack);    //         Color(color0);
  gr_points1->SetLineColor(kRed+1);    //         Color(color1);
  gr_points2->SetLineColor(kYellow-3); //         Color(color2);
  gr_points3->SetLineColor(kGreen+1);  //         Color(color3);
  gr_points4->SetLineColor(kGreen-9);  //         Color(color4);
  gr_points5->SetLineColor(kAzure+10); //         Color(color5);
  gr_points6->SetLineColor(kBlue-9);   //         Color(color6);
  gr_points7->SetLineColor(51);        //         Color(color7);

//   c_data = 1 
//   c_qcd =   ROOT.EColor.kRed+1
//   c_z =   ROOT.EColor.kOrange-3
//   c_vv =   ROOT.EColor.kYellow-3
//   c_t =  ROOT.EColor.kGreen+1
//   c_tb =  ROOT.EColor.kGreen-5
//   c_ttb = ROOT.EColor.kGreen-9
//   c_wl =  ROOT.EColor.kAzure+10
//   c_wc =  ROOT.EColor.kBlue+0
//   c_wcc = ROOT.EColor.kBlue-9
//   c_wbb = 51

  gr_points1->Draw("P");
  gr_points2->Draw("P,same");
  gr_points3->Draw("P,same");
  gr_points4->Draw("P,same");
  gr_points5->Draw("P,same");
  gr_points6->Draw("P,same");
  gr_points7->Draw("P,same");
  gr_points0->Draw("P,same");

  TLegend *leg4 = new TLegend(0.12,0.4,0.88,0.75 );
  leg4->SetFillColor(kWhite);
  leg4->AddEntry( gr_points0,"standard: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[0] % mes[0] % bs[0] % bes[0])),
   "L");

  leg4->AddEntry( gr_points1,"sideband Up: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[1] % mes[1] % bs[1] % bes[1])),
   "L");

  leg4->AddEntry( gr_points2,"sideband Down: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[2] % mes[2] % bs[2] % bes[2])),
   "L");

  leg4->AddEntry( gr_points3,"MET Up: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[3] % mes[3] % bs[3] % bes[3])),
   "L");

  leg4->AddEntry( gr_points4,"MET Down: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[4] % mes[4] % bs[4] % bes[4])),
   "L");

  leg4->AddEntry( gr_points5,"binning Up: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[5] % mes[5] % bs[5] % bes[5])),
   "L");

  leg4->AddEntry( gr_points6,"binning Down: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[6] % mes[6] % bs[6] % bes[6])),
   "L");

  leg4->AddEntry( gr_points7,"no #gamma iso: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[7] % mes[7] % bs[7] % bes[7])),
   "L");

  leg4->Draw("same");

  TLine *line0 = new TLine(xlow,xlow*ms[0]+bs[0],xhi,xhi*ms[0]+bs[0]);
  TLine *line1 = new TLine(xlow,xlow*ms[1]+bs[1],xhi,xhi*ms[1]+bs[1]);
  TLine *line2 = new TLine(xlow,xlow*ms[2]+bs[2],xhi,xhi*ms[2]+bs[2]);
  TLine *line3 = new TLine(xlow,xlow*ms[3]+bs[3],xhi,xhi*ms[3]+bs[3]);
  TLine *line4 = new TLine(xlow,xlow*ms[4]+bs[4],xhi,xhi*ms[4]+bs[4]);
  TLine *line5 = new TLine(xlow,xlow*ms[5]+bs[5],xhi,xhi*ms[5]+bs[5]);
  TLine *line6 = new TLine(xlow,xlow*ms[6]+bs[6],xhi,xhi*ms[6]+bs[6]);
  TLine *line7 = new TLine(xlow,xlow*ms[7]+bs[7],xhi,xhi*ms[7]+bs[7]);

  line0->SetLineColor(kBlack);   
  line1->SetLineColor(kRed+1);   
  line2->SetLineColor(kYellow-3);
  line3->SetLineColor(kGreen+1); 
  line4->SetLineColor(kGreen-9); 
  line5->SetLineColor(kAzure+10);
  line6->SetLineColor(kBlue-9);  
  line7->SetLineColor(51);       

  line1->Draw();
  line2->Draw();
  line3->Draw();
  line4->Draw();
  line5->Draw();
  line6->Draw();
  line7->Draw();
  line0->Draw();

  c3->Update();
  
  //gr_ratio->Draw("sames,A*");

  c3->SaveAs(outpath+"/Graph_FakeRatios.pdf");
  //c2->SaveAs(outpath+"/Graph_FakeRatio.C");
  //c2->SaveAs(outpath+"/Graph_FakeRatio.eps");
  //c2->SaveAs(wwwpath+"/Graph_FakeRatio.png");

  c3->Clear();

 return;

}
