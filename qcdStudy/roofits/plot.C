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
 sysnames.push_back("_eleTmpl");

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
 
 inpath = TString(Tsubmitbase+"/gitignore/"+Tversion+"/analyzed");
 outpath = TString(Tsubmitbase+"/gitignore/"+Tversion+"/plots");
 wwwpath = TString("/afs/hep.wisc.edu/user/tperry/www/MonoPhoton/qcdPlots/"+Tversion);
 
 datafile  = new TFile(inpath+"/mrg4bins_SinglePhoton.root","READ");
 edatafile = new TFile(inpath+"/mrg4bins_DoubleElectron.root","READ");
 mcfile    = new TFile(inpath+"/mrg4bins_GJets.root","READ");
 //qcdfile   = new TFile(inpath+"/mrg4bins_GJets.root","READ");
 qcdfile   = new TFile(inpath+"/mrg4bins_QCD.root","READ");
 outfile   = new TFile(outpath+"/Fitted.root","RECREATE");

 ofstream log;
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
   edatafile,
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
   sn,
   extraname
  );

 ratioss.push_back(ratios);
 ratioerrorss.push_back(ratioerrors);

 }


// calculate fake ratio error from fit up/down
sysnames.push_back("_errup");
sysnames.push_back("_errdn");
for(unsigned int sn=0; sn<2; ++sn){
 ratios.clear();
 ratioerrors.clear();

 //Using RooFit result
 getCorrectedFakeRatioErrUpDown(
  datafile,
  outfile,
  qcd_frac[0],
  qcd_frac_err[0],
  //qcdFrac,
  //qcdFracErr,
  ratios,
  ratioerrors,
  sn,
  extraname
 );

 // put in by hand frac err up/down
 qcdFrac.clear();
 qcdFracErr.clear();
 for(unsigned int pp=0; pp<nptbins; ++pp){
  double frqcd = qcd_frac[0].at(pp);
  if(sn==0){frqcd = frqcd + qcd_frac_err[0].at(pp);}
  if(sn==1){frqcd = frqcd - qcd_frac_err[0].at(pp);}
  qcdFrac.push_back(frqcd);
  qcdFracErr.push_back(qcd_frac_err[0].at(pp));
 }

 qcd_frac.push_back(qcdFrac);
 qcd_frac_err.push_back(qcdFracErr);
 
 ratioss.push_back(ratios);
 ratioerrorss.push_back(ratioerrors);

 getBinCenters(datafile, bincenters, binerrors, "");

 bincenterss.push_back(bincenters);
 binerrorss.push_back(binerrors);
}


 for(unsigned int sn=0; sn<sysnames.size(); ++sn){
   log<<boost::format("\nSysname: %s \n") % sysnames[sn];
  for(unsigned int j=0; j<nptbins; ++j){
    log<<boost::format(" %s\n") % ptbinnames[j];
    log<<boost::format("  QCD Fraction: %0.3d +- %0.3d \n") % qcd_frac[sn][j] % qcd_frac_err[sn][j];
    log<<boost::format("  Bin Center: %0.3d +- %0.3d \n") % bincenterss[sn][j] % binerrorss[sn][j];
    log<<boost::format("  Ratio: %0.3d +- %0.3d \n") % ratioss[sn][j] % ratioerrorss[sn][j];
  }
    log<<boost::format(" y = (%0.3d +- %0.3d)x + (%0.3d +- %0.3d) \n") % ms[sn] % mes[sn] % bs[sn] % bes[sn] ;
 }

// drawAllRates(extraname);
//// drawAllRatesRelative(extraname);
 drawAllRatesRelativeHist(extraname);

 log.close();

}


//------------------------------------
//  get fraction fits from templates
//------------------------------------
void plot::getFraction(
 TFile* datafile,  
 TFile* edatafile,  
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
  //printf("Starting %15s  %15s");
  std::cout<<"------------------------------------xxxxxxxxxxxxx"<<std::endl;
  std::cout<<"------------------------------------xxxxxxxxxxxxx"<<std::endl;
  std::cout<<"------------------------------------xxxxxxxxxxxxx"<<std::endl;
  std::cout<<"------------------------------------xxxxxxxxxxxxx"<<std::endl;
  std::cout<<"Starting "<<sysname<<" "<<extraname<<std::endl;

  TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500);   
  gStyle->SetOptStat(0);
  //gStyle->SetPadGridX(3);
  //gStyle->SetPadGridY(3);
  //gStyle->SetGridStyle(3);
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  ofstream logg;
  logg.open (outpath+"/logs"+sysname+extraname+".txt");

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
   logg<<boost::format("pT Range:  %s \n\n") % binrange;
  
    //set template histo names
    datahistname = "h_sig_sieieF5x5_"+binrange+sysname; //data histos with tigh ID and without sigIetaIeta cut
    fullhistname = "h_sig_sieieF5x5_"+binrange+sysname; //phojet histos with tight ID cut but without sigIetaIeta
    qcdhistname  = "h_bkg_sieieF5x5_"+binrange+sysname; //data histos with very loose id and  sideband of track iso

    if(sysname.EqualTo("_eleTmpl")){ 
     datahistname = "h_sig_sieieF5x5_"+binrange; //data histos with tigh ID and without sigIetaIeta cut
     fullhistname = "h_sig_sieieF5x5_"+binrange; //phojet histos with tight ID cut but without sigIetaIeta
     qcdhistname  = "h_bkg_sieieF5x5_"+binrange; //data histos with very loose id and  sideband of track iso
    }

    //// for closure test
    ////get "Data" template = MC (GJ + QCD) 
    //hdata = (TH1D*)mcfile->Get(datahistname)->Clone();
    //htemp = (TH1D*)qcdfile->Get(qcdhistname)->Clone();
    //hdata->Add(htemp);
    ////QCD histo
    //hqcd  = (TH1D*)qcdfile->Get(qcdhistname)->Clone();

    //  used and good
    //get Data template and QCD template from data
    datafile->cd();
    if(sysname.EqualTo("_eleTmpl")){ 
     hdata = (TH1D*)edatafile->Get(datahistname)->Clone();
     //hqcd  = (TH1D*)edatafile->Get(qcdhistname)->Clone();
    }
    else{ 
     hdata = (TH1D*)datafile->Get(datahistname)->Clone(); 
     //hqcd  = (TH1D*)datafile->Get(qcdhistname)->Clone();
    }
    //QCD histo
    hqcd  = (TH1D*)datafile->Get(qcdhistname)->Clone();

    double integ_data = hdata->Integral();
    double integ_qcd = hqcd->Integral();
    //get photon+jet signal tempalte
    mcfile->cd();
    hphojet  = (TH1D*)mcfile->Get(fullhistname)->Clone();
    double integ_phojet = hphojet->Integral();
    int nbins = hphojet->GetNbinsX();

    //canvas->cd();
    //hphojet->Draw();
    //canvas->Print(outpath+"/hphojet"+binrange+sysname+".pdf");
    //hdata->Draw();
    //canvas->Print(outpath+"/hdata"+binrange+sysname+".pdf");
    //hqcd->Draw();
    //canvas->Print(outpath+"/hqcd"+binrange+sysname+".pdf");


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
  
  double usig_gt     = hphojet->Integral(myreqbin+1, nbins);
  double usig_lt     = hphojet->Integral(1,myreqbin);
  double usig_tot    = hphojet->Integral();
  
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
  float sinincut = 0.0102;
    
 //Set some value not zero for roofit 
  for(int bincount = 1; bincount <= hqcd->GetNbinsX();bincount++){
   if(hqcd->GetBinContent(bincount) == 0.) hqcd->SetBinContent(bincount,1.e-06);
  }

  for(int bincount = 1; bincount <= hphojet->GetNbinsX();bincount++){
   if(hphojet->GetBinContent(bincount) == 0.) hphojet->SetBinContent(bincount,1.e-06);
  }

  // Variable Being Fit - fit over full range of sinin
  RooRealVar sinin("sinin","sigieie Title",sininmin,sininmax);
  sinin.setRange("fullrange",0.0,0.025);

  //set histograms pdfs - to be fit - Signal,QCD
  RooDataHist faketemplate("faketemplate","fake template",sinin,hqcd);
  RooHistPdf fakepdf("fakepdf","test hist fake pdf",sinin,faketemplate);

  RooDataHist realtemplate("realtemplate","real template",sinin,hphojet);
  RooHistPdf realpdf("realpdf","test hist real pdf",sinin,realtemplate);

   //variables that will contain real and fake estimates - can vary between 0-ndataentries
  RooRealVar signum("signum","signum",0,ndataentries);
  RooRealVar fakenum("fakenum","fakenum",0,ndataentries);

  // set data distribution to be fitted to
  RooDataHist data("data","data to be fitted to",sinin,hdata);
                       
  //set extended pdfs  
   //RooExtendPdf (const char *name, const char *title, const RooAbsPdf &pdf, const RooAbsReal &norm, const char *rangeName=0)
   ///////////////////////////////////////////////////////////////////////////
   //  RooExtendPdf is a wrappper around an existing PDF that adds a 
   //  parameteric extended likelihood term to the PDF, optionally divided by a 
   //  fractional term from a partial normalization of the PDF:
   //
   //  nExpected = N   _or Expected = N / frac 
   //
   //  where N is supplied as a RooAbsReal to RooExtendPdf.
   //  The fractional term is defined as
   //                          _       _ _   _  _
   //            Int(cutRegion[x]) pdf(x,y) dx dy 
   //     frac = ---------------_-------_-_---_--_ 
   //            Int(normRegion[x]) pdf(x,y) dx dy 
   //
   //        _                                                               _
   //  where x is the set of dependents involved in the selection region and y
   //  is the set of remaining dependents.
   //            _
   //  cutRegion[x] is an limited integration range that is contained in
   //  the nominal integration range normRegion[x[]
   //
  RooExtendPdf extpdfsig("Signal","extpdfsig",realpdf,signum,"fullrange");
  RooExtendPdf extpdffake("Background","extpdffake",fakepdf,fakenum,"fullrange");
                       
  //composite model pdf
  RooAddPdf model("model","sig + background",RooArgList(extpdfsig,extpdffake));
                       
  //make fit          
   model.fitTo(data,RooFit::Minos());
  // model.fitTo(data,RooFit::Hesse());


  // Do some calculations
  cout<<binrange<<endl;

  //get estimates and their errors (full range)
  // qcd fraction
  Double_t fakevalue = fakenum.getValV();
  Double_t fakeerrorhi = fakenum.getErrorHi();
  Double_t fakeerrorlo = fakenum.getErrorLo();
  Double_t fakeerrormax = max(fabs(fakeerrorhi),fabs(fakeerrorlo));
  //signal fraction                    
  Double_t sigvalue = signum.getValV();
  Double_t sigerrorhi = signum.getErrorHi();
  Double_t sigerrorlo = signum.getErrorLo();
  Double_t sigerrormax = max(fabs(sigerrorhi),fabs(sigerrorlo));

  // Find fraction of fake/real pdfs in full/partial range (normalized to 1) 
  // Define a range named "signal" in sieie from 0 - 0.0102
  //sinin.setRange("signal",0.0,0.025) ; // TOMTOMTOMTOM
  sinin.setRange("signal",0.0,0.0102) ;
  // create object representing integral over fakepdf
  RooAbsReal* i_full_fakepdf = fakepdf.createIntegral(sinin,NormSet(sinin),Range("fullrange")) ;
  RooAbsReal* i_sgnl_fakepdf = fakepdf.createIntegral(sinin,NormSet(sinin),Range("signal")) ;
  RooAbsReal* i_full_realpdf = realpdf.createIntegral(sinin,NormSet(sinin),Range("fullrange")) ;
  RooAbsReal* i_sgnl_realpdf = realpdf.createIntegral(sinin,NormSet(sinin),Range("signal")) ;

  cout<<boost::format("i_full_fakepdf->getVal()  %0.3d \n") % i_full_fakepdf->getVal()<<endl;
  cout<<boost::format("i_sgnl_fakepdf->getVal()  %0.3d \n") % i_sgnl_fakepdf->getVal()<<endl;
  cout<<boost::format("i_full_realpdf->getVal()  %0.3d \n") % i_full_realpdf->getVal()<<endl;
  cout<<boost::format("i_sgnl_realpdf->getVal()  %0.3d \n") % i_sgnl_realpdf->getVal()<<endl;

  logg<<" Fitted Fraction of pdf in 0<sieie<0.025  |  0<sieie<0.0102 \n";
 // logg<<boost::format("  QCD  %0.3d +- %0.3d | %0.3d +- %0.3d \n") % i_full_fakepdf->getVal() % i_full_fakepdf->getVal() % i_sgnl_fakepdf->getVal() % i_sgnl_fakepdf->getVal() ;
 // logg<<boost::format("  PHO  %0.3d +- %0.3d | %0.3d +- %0.3d \n") % i_full_realpdf->getVal() % i_full_realpdf->getVal() % i_sgnl_realpdf->getVal() % i_sgnl_realpdf->getVal() ;
  logg<<boost::format("  QCD  %0.3d | %0.3d \n") % i_full_fakepdf->getVal() % i_sgnl_fakepdf->getVal() ;
  logg<<boost::format("  PHO  %0.3d | %0.3d \n") % i_full_realpdf->getVal() % i_sgnl_realpdf->getVal() ;

  // fraction of fake/real in 0<sieie<0.025
  Double_t frac_fakeInSig = i_sgnl_fakepdf->getVal();
  Double_t frac_realInSig = i_sgnl_realpdf->getVal();

  // number of fake/real in 0<sieie<0.0102
  Double_t fakeInSig = fakevalue*frac_fakeInSig;
  Double_t realInSig = sigvalue*frac_realInSig;
  // error on fake/real in 0<sieie<0.0102
  Double_t fakeInSig_err = fakeerrormax*frac_fakeInSig;
  Double_t realInSig_err = sigerrormax*frac_realInSig;

  logg<<" Events in 0<sieie<0.025  |  0<sieie<0.0102 \n";
  //logg<<boost::format("  Prefit QCD:  %0.5d | %0.5d \n") %uqcd_tot % uqcd_lt;
  //logg<<boost::format("  Prefit PHO:  %0.5d | %0.5d \n") %usig_tot % usig_lt;
  //logg<<boost::format(" ---------------------------\n");
  logg<<boost::format("  QCD :  %0.5d | %0.5d \n") %fakevalue % fakeInSig ;
  logg<<boost::format("  PHO :  %0.5d | %0.5d \n") %sigvalue  % realInSig ;
  logg<<boost::format(" ---------------------------\n");
  logg<<boost::format("  TOT :  %0.5d | %0.5d \n") % (fakevalue+sigvalue)  % (fakeInSig+realInSig) ;
  logg<<boost::format("  Data:  %0.5d | %0.5d \n\n\n") %data_tot % data_lt;

  float frQCD = fakeInSig/(fakeInSig+realInSig);
  float frQCDerr = ( fakeInSig/(fakeInSig+realInSig) ) *
                   sqrt(
                    ( (fakeInSig_err/fakeInSig)*(fakeInSig_err/fakeInSig) )
                    + ( (realInSig_err/realInSig)*(realInSig_err/realInSig) )
                   ) ;

//  float frQCD = fakevalue/(sigvalue+fakevalue);
//  float frQCDerr = ( fakevalue/(sigvalue+fakevalue)) *
//                   sqrt(
//                    ( (fakeerrormax/fakevalue)*(fakeerrormax/fakevalue) )
//                    + ( (sigerrormax/sigvalue)*(sigerrormax/sigvalue) )
//                   ) ;


  cout<<" frQCD    : "<<frQCD<<endl;
  cout<<" frQCDerr : "<<frQCDerr<<endl;

  //put results into vector
  fractionQCD.push_back(frQCD);
  fractionQCDErr.push_back(frQCDerr);

//  cout<<" a"<<endl;
//
  //TString sndataentries = TString(boost::lexical_cast<string>( boost::format("%.1f") % ndataentries ));
  //TString sdata_lt      = TString(boost::lexical_cast<string>( boost::format("%.1f") % sdata_lt ));

  //TString sfakevalue    = TString(boost::lexical_cast<string>( boost::format("%.1f") % fakevalue    )); 
  //TString sfakeerrormax = TString(boost::lexical_cast<string>( boost::format("%.1f") % fakeerrormax )); 
  //TString ssigvalue     = TString(boost::lexical_cast<string>( boost::format("%.1f") % sigvalue     )); 
  //TString ssigerrormax  = TString(boost::lexical_cast<string>( boost::format("%.1f") % sigerrormax  )); 
  //TString sfitvalue     = TString(boost::lexical_cast<string>( boost::format("%.1f") % sqrt(sigvalue   *sigvalue   +fakevalue   *fakevalue   )  )); 
  //TString sfiterrormax  = TString(boost::lexical_cast<string>( boost::format("%.1f") % sqrt(sigerrormax*sigerrormax+fakeerrormax*fakeerrormax)  )); 

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
  xframe->SetMinimum(.9);
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
   lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");
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

   if(sysname.EqualTo("_eleTmpl")){ 
    datahistname = "h_sig_sieieF5x5_"+binrange+sysname;
   }
   canvas->SaveAs(outpath+"/Fitted_"+datahistname+extraname+".pdf");
   canvas->SaveAs(wwwpath+"/Fitted_"+datahistname+extraname+".pdf");

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
  logg.close();
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
   if(sysname.EqualTo("_eleTmpl")){ 
    datahistname = "h_sig_et_"+binrange; //data histos with tigh ID and without sigIetaIeta cut
   }
   else {datahistname = "h_sig_et_"+binrange+sysname; } //data histos with tigh ID and without sigIetaIeta cut

   //// for closure test
   ////get "Data" template = MC (GJ + QCD) 
   //hdata = (TH1D*)mcfile->Get(datahistname)->Clone();
   //htemp = (TH1D*)qcdfile->Get(qcdhistname)->Clone();
   //hdata->Add(htemp);

   //  used and good
   //get Data template and QCD template from data
   datafile->cd();
   if(sysname.EqualTo("_eleTmpl")){ 
    hdata = (TH1D*)edatafile->Get(datahistname)->Clone();
   }
   else { 
    hdata = (TH1D*)datafile->Get(datahistname)->Clone();
   }

   bcenters.push_back(hdata->GetMean());
   //berrors.push_back(hdata->GetMeanError());
   berrors.push_back(hdata->GetRMS());
   //std::cout<<" aaaa "<<sysname<<std::endl;
   //std::cout<<"aaaaaaaa bincenter: "<<hdata->GetMean()<<std::endl;
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
                                 int sn,
                                 TString extraname
                                ){



  sysname = sysnames[sn];

  ofstream logg;
  logg.open (outpath+"/logs"+sysname+".txt", ios::out | ios::app );


  //define histo and root files  
  TH1D  *hdataden;
  TH1D  *hdatanum;

  string shistname;
  TString binrange;
  TString denhistname;
  TString numhistname;

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
   logg<<boost::format("pT Range:  %s \n\n") % binrange;

  //set template histo names
  if(sysname.EqualTo("_eleTmpl")){ 
   numhistname = "h_sig_sieieF5x5_"+binrange; //data histos with tight ID cut but without sigIetaIeta
   denhistname = "h_den_sieieF5x5_"+binrange; //data histos fail loose ID and pass Vloose ID
  }
  else{
   numhistname = "h_sig_sieieF5x5_"+binrange+sysname; //data histos with tight ID cut but without sigIetaIeta
   denhistname = "h_den_sieieF5x5_"+binrange+sysname; //data histos fail loose ID and pass Vloose ID
  }

  //get the numerator and denominator
  datafile->cd();
  //num
  hdatanum = (TH1D*)datafile->Get(numhistname)->Clone();    
  //deno
  hdataden = (TH1D*)datafile->Get(denhistname)->Clone();
  ////num
  //hdatanum = (TH1D*)edatafile->Get(numhistname)->Clone();    
  ////deno
  //hdataden = (TH1D*)edatafile->Get(denhistname)->Clone();

   
  // bin corresponding to medium wp sieie cut
  double sieiebin;
  //sieiebin = hdatanum->FindBin(0.025);
  sieiebin = hdatanum->FindBin(0.0102);

  std::cout<< double((corrFactor[ptb]))<<std::endl;
  size_num_uncor[ptb] = double(hdatanum->Integral(1,sieiebin));
  size_num_corr[ptb] = double((corrFactor[ptb])) * double(hdatanum->Integral(1,sieiebin));
  size_den[ptb] = double(hdataden->Integral(1,sieiebin));
  logg<<boost::format(" Sizes:   Num(uncorr): %0.1f  Num(corr): %0.1f  Den: %0.1f \n\n")
    % size_num_uncor[ptb] % size_num_corr[ptb] % size_den[ptb] ;

  double tmp_err_num_corr;
  //hdatanum->IntegralAndError(1, sieiebin, err_num_uncor[ptb]);
  hdatanum->IntegralAndError(1, sieiebin, tmp_err_num_corr);
  //err_num_corr[varcounter] = sqrt( err_num_corr[varcounter]*err_num_corr[varcounter] + 
  err_num_corr[ptb] = sqrt( tmp_err_num_corr*tmp_err_num_corr + 
                              corrErr[ptb]*corrErr[ptb] );
  hdataden->IntegralAndError(1, sieiebin, err_den[ptb]);
  //err_den[ptb]=0;

  }

  logg.close();

  TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500);   

  Double_t *mid = new Double_t[nptbins]; //{182.5, 220., 325., 700. };
  Double_t *miderr = new Double_t[nptbins];

  for(unsigned int k=0; k<nptbins; ++k){
   mid[k]=bincenterss[sn][k];
   miderr[k]=binerrorss[sn][k];
   //std::cout<<"xxxxxxxxxxxxxxxxxxxxxxx "<<bincenterss[sn][k]<<std::endl;
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

  TH1F *hr = canvas->DrawFrame(175.,0.1,1000.,1000000.,"");
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
  lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");
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
  canvas->SaveAs(wwwpath+"/Graph_NuNcD"+sysname+".pdf");

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
  lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");

  TLegend *leg3 = new TLegend(0.55,0.6,0.88,0.88 );
  leg3->SetFillColor(kWhite);
  leg3->AddEntry( gr_ratio,"Num./Den.", "L");
  leg3->Draw("same");

  gr_ratio->Draw("P");
  gr_ratio_noe->Draw("P,sames");
  gr_ratio_noe->Fit("pol1");
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

  //std::cout<<"Npoints: "<<gr_ratio->GetN()<<std::endl;
  //Double_t thex;
  //Double_t they;
  //for(unsigned int b=0; b<gr_ratio->GetN(); ++b){
  // gr_ratio->GetPoint(b,thex,they);
  // std::cout<<"Point "<<b<<": "<<thex<<" "<<they<<std::endl;
  //}
  //std::cout<<"Fit : "<<chi2<<" "<<p0<<" "<<p1<<std::endl;

  c2->Update();
  
  c2->SaveAs(outpath+"/Graph_FakeRatio"+sysname+extraname+".pdf");
  c2->SaveAs(wwwpath+"/Graph_FakeRatio"+sysname+extraname+".pdf");

  c2->Clear();

  return;

}



//----------------------------------------------
//   Get the Corrected Fake Ratio
//--------------------------------------------
void plot::getCorrectedFakeRatioErrUpDown(TFile* datafile,  //<----data file
                                 TFile* outfile,   //<---output file
                                 vector<double> corrFactor, //<--input corr. factor = QCD fraction
                                 vector<double> corrErr,    //<--input error
                                 vector<Double_t>& ratios,
                                 vector<Double_t>& ratioerrors,
                                 int sn,
                                 TString extraname
                                ){


  if(sn==0){sysname = "_errup";}
  if(sn==1){sysname = "_errdn";}

  ofstream logg;
  logg.open (outpath+"/logs"+sysname+".txt", ios::out | ios::app );


  //define histo and root files  
  TH1D  *hdataden;
  TH1D  *hdatanum;

  string shistname;
  TString binrange;
  TString denhistname;
  TString numhistname;

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
   logg<<boost::format("pT Range:  %s \n\n") % binrange;

  //set template histo names
  numhistname = "h_sig_sieieF5x5_"+binrange; //data histos with tight ID cut but without sigIetaIeta
  denhistname = "h_den_sieieF5x5_"+binrange; //data histos fail loose ID and pass Vloose ID

  //get the numerator and denominator
  datafile->cd();
  //num
  hdatanum = (TH1D*)datafile->Get(numhistname)->Clone();    
  //deno
  hdataden = (TH1D*)datafile->Get(denhistname)->Clone();

   
  // bin corresponding to medium wp sieie cut
  double sieiebin;
  //sieiebin = hdatanum->FindBin(0.025);
  sieiebin = hdatanum->FindBin(0.0102);

  double cf;
  if(sn==0){cf = corrFactor[ptb] + corrErr[ptb];}
  if(sn==1){cf = corrFactor[ptb] - corrErr[ptb];}
  //cf = corrFactor[ptb];
  std::cout<<"------------------------------------xddddddddxxxx"<<std::endl;
  std::cout<<"------------------------------------xddddddddxxxx"<<std::endl;
  std::cout<<"------------------------------------xddddddddxxxx"<<std::endl;
  std::cout<<"------------------------------------xddddddddxxxx"<<std::endl;
  std::cout<<"Starting "<<sysname<<" "<<extraname<<std::endl;
  std::cout<<corrFactor[ptb]<<std::endl;
  std::cout<<corrErr[ptb]<<std::endl;
  std::cout<<cf<<std::endl;

  size_num_uncor[ptb] = double(hdatanum->Integral(1,sieiebin));
  size_num_corr[ptb] = double(cf) * double(hdatanum->Integral(1,sieiebin));
  size_den[ptb] = double(hdataden->Integral(1,sieiebin));
  logg<<boost::format(" Sizes:   Num(uncorr): %0.1f  Num(corr): %0.1f  Den: %0.1f \n\n")
    % size_num_uncor[ptb] % size_num_corr[ptb] % size_den[ptb] ;

  double tmp_err_num_corr;
  //hdatanum->IntegralAndError(1, sieiebin, err_num_uncor[ptb]);
  hdatanum->IntegralAndError(1, sieiebin, tmp_err_num_corr);
  //err_num_corr[varcounter] = sqrt( err_num_corr[varcounter]*err_num_corr[varcounter] + 
  err_num_corr[ptb] = sqrt( tmp_err_num_corr*tmp_err_num_corr + 
                              corrErr[ptb]*corrErr[ptb] );
  hdataden->IntegralAndError(1, sieiebin, err_den[ptb]);
  //err_den[ptb]=0;

  }

  logg.close();

  TCanvas* canvas = new TCanvas("canvas","canvas",900,100,500,500);   

  Double_t *mid = new Double_t[nptbins]; //{182.5, 220., 325., 700. };
  Double_t *miderr = new Double_t[nptbins];

  for(unsigned int k=0; k<nptbins; ++k){
   mid[k]=bincenterss[sn][k];
   miderr[k]=binerrorss[sn][k];
   //std::cout<<"xxxxxxxxxxxxxxxxxxxxxxx "<<bincenterss[sn][k]<<std::endl;
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

  TH1F *hr = canvas->DrawFrame(175.,0.1,1000.,1000000.,"");
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
  lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");
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
  canvas->SaveAs(wwwpath+"/Graph_NuNcD"+sysname+".pdf");

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
  lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");

  TLegend *leg3 = new TLegend(0.55,0.6,0.88,0.88 );
  leg3->SetFillColor(kWhite);
  leg3->AddEntry( gr_ratio,"Num./Den.", "L");
  leg3->Draw("same");

  gr_ratio->Draw("P");
  gr_ratio_noe->Draw("P,sames");
  gr_ratio_noe->Fit("pol1");
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

  //std::cout<<"Npoints: "<<gr_ratio->GetN()<<std::endl;
  //Double_t thex;
  //Double_t they;
  //for(unsigned int b=0; b<gr_ratio->GetN(); ++b){
  // gr_ratio->GetPoint(b,thex,they);
  // std::cout<<"Point "<<b<<": "<<thex<<" "<<they<<std::endl;
  //}
  //std::cout<<"Fit : "<<chi2<<" "<<p0<<" "<<p1<<std::endl;

  c2->Update();
  
  c2->SaveAs(outpath+"/Graph_FakeRatio"+sysname+extraname+".pdf");
  c2->SaveAs(wwwpath+"/Graph_FakeRatio"+sysname+extraname+".pdf");

  c2->Clear();

  return;

}



void plot::drawAllRates(TString extraname){

  TCanvas* c3 = new TCanvas("c3","c3",900,100,500,500);   

  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  double xmax;
  xmax = 600.;

  Double_t xlow = 175.;
  Double_t xhi = 600.;
  TH1F *hs = c3->DrawFrame(xlow,0.,xhi,0.5,"");
  hs->SetXTitle("photon pT [GeV]");
  hs->SetYTitle("Fake Ratio"); 

  TText* title = new TText(1,1,"") ;
  title->SetTextSize(0.07);
  title->SetTextColor(kBlack);
  title->SetTextAlign(13);
  title->SetTextFont(62);
  //title->DrawTextNDC(0.17,0.87,"CMS");

  TText* extra = new TText(1,1,"") ;
  extra->SetTextSize(0.05);
  extra->SetTextColor(kBlack);
  extra->SetTextAlign(13);
  extra->SetTextFont(52);
  //extra->DrawTextNDC(0.17,0.81,"Preliminary");
  //xframe->addObject(extra);

  TText* lumi = new TText(1,1,"") ;
  lumi->SetTextSize(0.05);
  lumi->SetTextColor(kBlack);
  lumi->SetTextAlign(31);
  lumi->SetTextFont(42);
  //lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");
  //xframe->addObject(lumi);

  title->DrawTextNDC(0.17,0.87,"CMS");
  extra->DrawTextNDC(0.17,0.81,"Preliminary");
  lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");

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
  Double_t *mid8 = new Double_t[nptbins];
  for(unsigned int k=0; k<nptbins; ++k){ mid0[k]=bincenterss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid1[k]=bincenterss[1][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid2[k]=bincenterss[2][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid3[k]=bincenterss[3][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid4[k]=bincenterss[4][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid5[k]=bincenterss[5][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid6[k]=bincenterss[6][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid7[k]=bincenterss[7][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid8[k]=bincenterss[8][k]; }

  Double_t *miderr0 = new Double_t[nptbins];
  Double_t *miderr1 = new Double_t[nptbins];
  Double_t *miderr2 = new Double_t[nptbins];
  Double_t *miderr3 = new Double_t[nptbins];
  Double_t *miderr4 = new Double_t[nptbins];
  Double_t *miderr5 = new Double_t[nptbins];
  Double_t *miderr6 = new Double_t[nptbins];
  Double_t *miderr7 = new Double_t[nptbins];
  Double_t *miderr8 = new Double_t[nptbins];
  for(unsigned int k=0; k<nptbins; ++k){ miderr0[k]=binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr1[k]=binerrorss[1][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr2[k]=binerrorss[2][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr3[k]=binerrorss[3][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr4[k]=binerrorss[4][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr5[k]=binerrorss[5][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr6[k]=binerrorss[6][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr7[k]=binerrorss[7][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr8[k]=binerrorss[8][k]; }
  
  // y
  Double_t *points0 = new Double_t[4];
  Double_t *points1 = new Double_t[4];
  Double_t *points2 = new Double_t[4];
  Double_t *points3 = new Double_t[4];
  Double_t *points4 = new Double_t[4];
  Double_t *points5 = new Double_t[4];
  Double_t *points6 = new Double_t[4];
  Double_t *points7 = new Double_t[4];
  Double_t *points8 = new Double_t[4];
  for(unsigned int p=0; p<4; ++p){ points0[p]=ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points1[p]=ratioss[1][p]; }
  for(unsigned int p=0; p<4; ++p){ points2[p]=ratioss[2][p]; }
  for(unsigned int p=0; p<4; ++p){ points3[p]=ratioss[3][p]; }
  for(unsigned int p=0; p<4; ++p){ points4[p]=ratioss[4][p]; }
  for(unsigned int p=0; p<4; ++p){ points5[p]=ratioss[5][p]; }
  for(unsigned int p=0; p<4; ++p){ points6[p]=ratioss[6][p]; }
  for(unsigned int p=0; p<4; ++p){ points7[p]=ratioss[7][p]; }
  for(unsigned int p=0; p<4; ++p){ points8[p]=ratioss[8][p]; }

  Double_t *pointerrors0 = new Double_t[4];
  Double_t *pointerrors1 = new Double_t[4];
  Double_t *pointerrors2 = new Double_t[4];
  Double_t *pointerrors3 = new Double_t[4];
  Double_t *pointerrors4 = new Double_t[4];
  Double_t *pointerrors5 = new Double_t[4];
  Double_t *pointerrors6 = new Double_t[4];
  Double_t *pointerrors7 = new Double_t[4];
  Double_t *pointerrors8 = new Double_t[4];
  for(unsigned int p=0; p<4; ++p){ pointerrors0[p]=ratioerrorss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors1[p]=ratioerrorss[1][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors2[p]=ratioerrorss[2][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors3[p]=ratioerrorss[3][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors4[p]=ratioerrorss[4][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors5[p]=ratioerrorss[5][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors6[p]=ratioerrorss[6][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors7[p]=ratioerrorss[7][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors8[p]=ratioerrorss[8][p]; }

  // fill the graphs
  TGraphErrors *gr_points0 = new TGraphErrors(4,mid0,points0,miderr0,pointerrors0);
  TGraphErrors *gr_points1 = new TGraphErrors(4,mid1,points1,miderr1,pointerrors1);
  TGraphErrors *gr_points2 = new TGraphErrors(4,mid2,points2,miderr2,pointerrors2);
  TGraphErrors *gr_points3 = new TGraphErrors(4,mid3,points3,miderr3,pointerrors3);
  TGraphErrors *gr_points4 = new TGraphErrors(4,mid4,points4,miderr4,pointerrors4);
  TGraphErrors *gr_points5 = new TGraphErrors(4,mid5,points5,miderr5,pointerrors5);
  TGraphErrors *gr_points6 = new TGraphErrors(4,mid6,points6,miderr6,pointerrors6);
  TGraphErrors *gr_points7 = new TGraphErrors(4,mid7,points7,miderr7,pointerrors7);
  TGraphErrors *gr_points8 = new TGraphErrors(4,mid8,points8,miderr8,pointerrors8);

  // set colors
  //  
  //  _sbUp
  //  _sbDown
  //  _metUp
  //  _metDown
  //  _binUp
  //  _binDown
  //  _noPiso
  //  _eleTmpl
  Int_t ci1  = 1001;
  Int_t ci2  = 1002;
  Int_t ci3  = 1003;
  Int_t ci4  = 1004;
  Int_t ci5  = 1005;
  Int_t ci6  = 1006;
  Int_t ci7  = 1007;
  Int_t ci8  = 1008;
  Int_t ci9  = 1009;
  Int_t ci10 = 1010;
  //Color_t mycolor = new TColor(ci, 0.1,0.2,0.3,"mycolor");
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
  //TColor *mycolor = gROOT->GetColor(kOrange);
  gr_points0->SetLineColor(kBlack);                          
  gr_points1->SetLineColor(ci1);
  gr_points2->SetLineColor(ci2);
  gr_points3->SetLineColor(ci3);
  gr_points4->SetLineColor(ci4);
  gr_points5->SetLineColor(ci5);
  gr_points6->SetLineColor(ci6);
  gr_points7->SetLineColor(ci7);
  gr_points8->SetLineColor(ci8);

  gr_points0->SetMarkerColor(kBlack                     );                          
  gr_points2->SetMarkerColor(ci1);
  gr_points1->SetMarkerColor(ci2);
  gr_points4->SetMarkerColor(ci3);
  gr_points3->SetMarkerColor(ci4);
  gr_points6->SetMarkerColor(ci5);
  gr_points5->SetMarkerColor(ci6);
  gr_points7->SetMarkerColor(ci7);
  gr_points8->SetMarkerColor(ci8);

  gr_points0->SetLineWidth(2); //         Color(color0);
  gr_points1->SetLineWidth(2); //         Color(color1);
  gr_points2->SetLineWidth(2); //         Color(color2);
  gr_points3->SetLineWidth(2); //         Color(color3);
  gr_points4->SetLineWidth(2); //         Color(color4);
  gr_points5->SetLineWidth(2); //         Color(color5);
  gr_points6->SetLineWidth(2); //         Color(color6);
  gr_points7->SetLineWidth(2); //         Color(color7);

  gr_points1->Draw("P");
  gr_points2->Draw("P,same");
  gr_points3->Draw("P,same");
  gr_points4->Draw("P,same");
  gr_points5->Draw("P,same");
  gr_points6->Draw("P,same");
  gr_points7->Draw("P,same");
  gr_points8->Draw("P,same");
  gr_points0->Draw("P,same");

  TLegend *leg4 = new TLegend(0.12,0.4,0.88,0.75 );
  leg4->SetFillColor(kWhite);
  leg4->AddEntry( gr_points0,"standard: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[0] % mes[0] % bs[0] % bes[0])),
   "L");

  leg4->AddEntry( gr_points7,"no #gamma iso: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[7] % mes[7] % bs[7] % bes[7])),
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

  leg4->AddEntry( gr_points8,"ele template: y = "
   +TString(boost::lexical_cast<string>(boost::format(
    "(%0.5f +- %0.5f) x + (%0.3f +- %0.3f)")
    % ms[8] % mes[8] % bs[8] % bes[8])),
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
  TLine *line8 = new TLine(xlow,xlow*ms[8]+bs[8],xhi,xhi*ms[8]+bs[8]);

  line0->SetLineColor(kBlack                     );                          
  line1->SetLineColor(ci1);
  line2->SetLineColor(ci2);
  line3->SetLineColor(ci3);
  line4->SetLineColor(ci4);
  line5->SetLineColor(ci5);
  line6->SetLineColor(ci6);
  line7->SetLineColor(ci7);
  line8->SetLineColor(ci8);

  line1->Draw();
  line2->Draw();
  line3->Draw();
  line4->Draw();
  line5->Draw();
  line6->Draw();
  line7->Draw();
  line8->Draw();
  line0->Draw();

  c3->Update();
  
  //gr_ratio->Draw("sames,A*");

  c3->SaveAs(outpath+"/Graph_FakeRatios"+extraname+".pdf");
  c3->SaveAs(wwwpath+"/Graph_FakeRatios"+extraname+".pdf");

  c3->Clear();

 return;

}

void plot::drawAllRatesRelative(TString extraname){

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
  TH1F *hs = c3->DrawFrame(xlow,-0.3,xhi,0.8,"");
  hs->SetXTitle("photon pT [GeV]");
  hs->SetYTitle("Relative Fake Ratio"); 

  TText* title = new TText(1,1,"") ;
  title->SetTextSize(0.07);
  title->SetTextColor(kBlack);
  title->SetTextAlign(13);
  title->SetTextFont(62);
  //title->DrawTextNDC(0.17,0.87,"CMS");

  TText* extra = new TText(1,1,"") ;
  extra->SetTextSize(0.05);
  extra->SetTextColor(kBlack);
  extra->SetTextAlign(13);
  extra->SetTextFont(52);
  //extra->DrawTextNDC(0.17,0.81,"Preliminary");
  //xframe->addObject(extra);

  TText* lumi = new TText(1,1,"") ;
  lumi->SetTextSize(0.05);
  lumi->SetTextColor(kBlack);
  lumi->SetTextAlign(31);
  lumi->SetTextFont(42);
  //lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");
  //xframe->addObject(lumi);

  title->DrawTextNDC(0.17,0.87,"CMS");
  extra->DrawTextNDC(0.17,0.81,"Preliminary");
  lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");

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
  for(unsigned int k=0; k<nptbins; ++k){ mid1[k]=bincenterss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid2[k]=bincenterss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid3[k]=bincenterss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid4[k]=bincenterss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid5[k]=bincenterss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid6[k]=bincenterss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ mid7[k]=bincenterss[0][k]; }

  Double_t *miderr0 = new Double_t[nptbins];
  Double_t *miderr1 = new Double_t[nptbins];
  Double_t *miderr2 = new Double_t[nptbins];
  Double_t *miderr3 = new Double_t[nptbins];
  Double_t *miderr4 = new Double_t[nptbins];
  Double_t *miderr5 = new Double_t[nptbins];
  Double_t *miderr6 = new Double_t[nptbins];
  Double_t *miderr7 = new Double_t[nptbins];
  for(unsigned int k=0; k<nptbins; ++k){ miderr0[k]=0.;} //binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr1[k]=0.;} //binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr2[k]=0.;} //binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr3[k]=0.;} //binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr4[k]=0.;} //binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr5[k]=0.;} //binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr6[k]=0.;} //binerrorss[0][k]; }
  for(unsigned int k=0; k<nptbins; ++k){ miderr7[k]=0.;} //binerrorss[0][k]; }
  
  // y
  Double_t *points0 = new Double_t[4];
  Double_t *points1 = new Double_t[4];
  Double_t *points2 = new Double_t[4];
  Double_t *points3 = new Double_t[4];
  Double_t *points4 = new Double_t[4];
  Double_t *points5 = new Double_t[4];
  Double_t *points6 = new Double_t[4];
  Double_t *points7 = new Double_t[4];
  for(unsigned int p=0; p<4; ++p){ points0[p]=(ratioss[0][p]-ratioss[0][p])/ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points1[p]=(ratioss[1][p]-ratioss[0][p])/ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points2[p]=(ratioss[2][p]-ratioss[0][p])/ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points3[p]=(ratioss[3][p]-ratioss[0][p])/ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points4[p]=(ratioss[4][p]-ratioss[0][p])/ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points5[p]=(ratioss[5][p]-ratioss[0][p])/ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points6[p]=(ratioss[6][p]-ratioss[0][p])/ratioss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ points7[p]=(ratioss[7][p]-ratioss[0][p])/ratioss[0][p]; }

  Double_t *pointerrors0 = new Double_t[4];
  Double_t *pointerrors1 = new Double_t[4];
  Double_t *pointerrors2 = new Double_t[4];
  Double_t *pointerrors3 = new Double_t[4];
  Double_t *pointerrors4 = new Double_t[4];
  Double_t *pointerrors5 = new Double_t[4];
  Double_t *pointerrors6 = new Double_t[4];
  Double_t *pointerrors7 = new Double_t[4];
  for(unsigned int p=0; p<4; ++p){ pointerrors0[p]=0.;} //ratioerrorss[0][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors1[p]=0.;} //ratioerrorss[1][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors2[p]=0.;} //ratioerrorss[2][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors3[p]=0.;} //ratioerrorss[3][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors4[p]=0.;} //ratioerrorss[4][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors5[p]=0.;} //ratioerrorss[5][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors6[p]=0.;} //ratioerrorss[6][p]; }
  for(unsigned int p=0; p<4; ++p){ pointerrors7[p]=0.;} //ratioerrorss[7][p]; }

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
  gr_points0->SetMarkerColor(kBlack                     );                          
  gr_points2->SetMarkerColor(TColor::GetColor("#33a02c"));
  gr_points1->SetMarkerColor(TColor::GetColor("#b2df8a"));
  gr_points4->SetMarkerColor(TColor::GetColor("#e31a1c"));
  gr_points3->SetMarkerColor(TColor::GetColor("#fb9a99"));
  gr_points6->SetMarkerColor(TColor::GetColor("#1f78b4"));
  gr_points5->SetMarkerColor(TColor::GetColor("#a6cee3"));
  gr_points7->SetMarkerColor(TColor::GetColor("#b15928"));

  gr_points0->SetLineColor(kBlack                     );                          
  gr_points2->SetLineColor(TColor::GetColor("#33a02c"));
  gr_points1->SetLineColor(TColor::GetColor("#b2df8a"));
  gr_points4->SetLineColor(TColor::GetColor("#e31a1c"));
  gr_points3->SetLineColor(TColor::GetColor("#fb9a99"));
  gr_points6->SetLineColor(TColor::GetColor("#1f78b4"));
  gr_points5->SetLineColor(TColor::GetColor("#a6cee3"));
  gr_points7->SetLineColor(TColor::GetColor("#b15928"));

  gr_points0->SetLineWidth(3); //         Color(color0);
  gr_points1->SetLineWidth(3); //         Color(color1);
  gr_points2->SetLineWidth(3); //         Color(color2);
  gr_points3->SetLineWidth(3); //         Color(color3);
  gr_points4->SetLineWidth(3); //         Color(color4);
  gr_points5->SetLineWidth(3); //         Color(color5);
  gr_points6->SetLineWidth(3); //         Color(color6);
  gr_points7->SetLineWidth(3); //         Color(color7);

  gr_points1->Draw("L");
  gr_points2->Draw("L,same");
  gr_points3->Draw("L,same");
  gr_points4->Draw("L,same");
  gr_points5->Draw("L,same");
  gr_points6->Draw("L,same");
  gr_points7->Draw("L,same");
  gr_points0->Draw("L,same");

  TLegend *leg4 = new TLegend(0.5,0.6,0.88,0.88 );
  leg4->SetFillColor(kWhite);
  leg4->AddEntry( gr_points0,"standard", "L");
  leg4->AddEntry( gr_points7,"bkg. no #gamma iso req", "L");
  leg4->AddEntry( gr_points1,"sideband Up", "L");
  leg4->AddEntry( gr_points2,"sideband Down", "L");
  leg4->AddEntry( gr_points3,"MET Up", "L");
  leg4->AddEntry( gr_points4,"MET Down", "L");
  leg4->AddEntry( gr_points5,"binning Up", "L");
  leg4->AddEntry( gr_points6,"binning Down","L");
  leg4->Draw("same");

  //TLine *line0 = new TLine(xlow,xlow*ms[0]+bs[0],xhi,xhi*ms[0]+bs[0]);
  //TLine *line1 = new TLine(xlow,xlow*ms[1]+bs[1],xhi,xhi*ms[1]+bs[1]);
  //TLine *line2 = new TLine(xlow,xlow*ms[2]+bs[2],xhi,xhi*ms[2]+bs[2]);
  //TLine *line3 = new TLine(xlow,xlow*ms[3]+bs[3],xhi,xhi*ms[3]+bs[3]);
  //TLine *line4 = new TLine(xlow,xlow*ms[4]+bs[4],xhi,xhi*ms[4]+bs[4]);
  //TLine *line5 = new TLine(xlow,xlow*ms[5]+bs[5],xhi,xhi*ms[5]+bs[5]);
  //TLine *line6 = new TLine(xlow,xlow*ms[6]+bs[6],xhi,xhi*ms[6]+bs[6]);
  //TLine *line7 = new TLine(xlow,xlow*ms[7]+bs[7],xhi,xhi*ms[7]+bs[7]);

  //line0->SetLineColor(kBlack);   
  //line1->SetLineColor(kRed+1);   
  //line2->SetLineColor(kYellow-3);
  //line3->SetLineColor(kGreen+1); 
  //line4->SetLineColor(kGreen-9); 
  //line5->SetLineColor(kAzure+10);
  //line6->SetLineColor(kBlue-9);  
  //line7->SetLineColor(51);       

  //line1->Draw();
  //line2->Draw();
  //line3->Draw();
  //line4->Draw();
  //line5->Draw();
  //line6->Draw();
  //line7->Draw();
  //line0->Draw();

  c3->Update();
  
  //gr_ratio->Draw("sames,A*");

  c3->SaveAs(outpath+"/Graph_FakeRatiosRelative"+extraname+".pdf");
  c3->SaveAs(wwwpath+"/Graph_FakeRatiosRelative"+extraname+".pdf");

  c3->Clear();

 return;

}

void plot::drawAllRatesRelativeHist(TString extraname){

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
  xmax = 1000.;

  Double_t xlow = 175.;
  Double_t xhi = 1000.;
  //TH1F *hs = c2->DrawFrame(0.,0.,1000.,0.5,"");
  //TH1F *hs = c3->DrawFrame(175.,0.,1000.,0.5,"");
  //TH1F *hs = c3->DrawFrame(xlow,-2.,xhi,2.,"");
  TH1F *hs = c3->DrawFrame(xlow,-0.3,xhi,0.8,"");
  hs->SetXTitle("photon pT [GeV]");
  hs->SetYTitle("Relative Fake Ratio"); 

  TText* title = new TText(1,1,"") ;
  title->SetTextSize(0.07);
  title->SetTextColor(kBlack);
  title->SetTextAlign(13);
  title->SetTextFont(62);
  //title->DrawTextNDC(0.17,0.87,"CMS");

  TText* extra = new TText(1,1,"") ;
  extra->SetTextSize(0.05);
  extra->SetTextColor(kBlack);
  extra->SetTextAlign(13);
  extra->SetTextFont(52);
  //extra->DrawTextNDC(0.17,0.81,"Preliminary");
  //xframe->addObject(extra);

  TText* lumi = new TText(1,1,"") ;
  lumi->SetTextSize(0.05);
  lumi->SetTextColor(kBlack);
  lumi->SetTextAlign(31);
  lumi->SetTextFont(42);
  //lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");
  //xframe->addObject(lumi);

  title->DrawTextNDC(0.17,0.87,"CMS");
  extra->DrawTextNDC(0.17,0.81,"Preliminary");
  lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");

  ////////// Draw points first, then lines

  const Int_t NBINS = 4;
  Double_t edges[NBINS+1] = {175,190,250,400,1000};

  // fill the graphs
  TH1* h_sys0 = new TH1D( "h_sys0", "h_sys0", NBINS, edges );
  TH1* h_sys1 = new TH1D( "h_sys1", "h_sys1", NBINS, edges );
  TH1* h_sys2 = new TH1D( "h_sys2", "h_sys2", NBINS, edges );
  TH1* h_sys3 = new TH1D( "h_sys3", "h_sys3", NBINS, edges );
  TH1* h_sys4 = new TH1D( "h_sys4", "h_sys4", NBINS, edges );
  TH1* h_sys5 = new TH1D( "h_sys5", "h_sys5", NBINS, edges );
  TH1* h_sys6 = new TH1D( "h_sys6", "h_sys6", NBINS, edges );
  TH1* h_sys7 = new TH1D( "h_sys7", "h_sys7", NBINS, edges );
  TH1* h_sys8 = new TH1D( "h_sys8", "h_sys8", NBINS, edges );
  TH1* h_sys9 = new TH1D( "h_sys9", "h_sys9", NBINS, edges );
  TH1* h_sys10 = new TH1D( "h_sys10", "h_sys10", NBINS, edges );

  // y
  Double_t *points0 = new Double_t[4];
  Double_t *points1 = new Double_t[4];
  Double_t *points2 = new Double_t[4];
  Double_t *points3 = new Double_t[4];
  Double_t *points4 = new Double_t[4];
  Double_t *points5 = new Double_t[4];
  Double_t *points6 = new Double_t[4];
  Double_t *points7 = new Double_t[4];
  Double_t *points8 = new Double_t[4];
  Double_t *points9 = new Double_t[4];
  Double_t *points10 = new Double_t[4];
  for(unsigned int p=0; p<4; ++p){ h_sys0->SetBinContent( p+1, (ratioss[0][p]-ratioss[0][p])/ratioss[0][p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys1->SetBinContent( p+1, (ratioss[1][p]-ratioss[0][p])/ratioss[0][p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys2->SetBinContent( p+1, (ratioss[2][p]-ratioss[0][p])/ratioss[0][p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys3->SetBinContent( p+1, (ratioss[3][p]-ratioss[0][p])/ratioss[0][p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys4->SetBinContent( p+1, (ratioss[4][p]-ratioss[0][p])/ratioss[0][p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys5->SetBinContent( p+1, (ratioss[5][p]-ratioss[0][p])/ratioss[0][p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys6->SetBinContent( p+1, (ratioss[6][p]-ratioss[0][p])/ratioss[0][p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys7->SetBinContent( p+1, (ratioss[7][p]-ratioss[0][p])/ratioss[0][p] ); }

  for(unsigned int p=0; p<4; ++p){ h_sys8->SetBinContent( p+1, (ratioss[8][p]-ratioss[0][p])/ratioss[0][p] ); }
  for(unsigned int p=0; p<4; ++p){ h_sys9->SetBinContent( p+1, (ratioss[9][p]-ratioss[0][p])/ratioss[0][p] ); }

//  for(unsigned int p=0; p<4; ++p){ h_sys8->SetBinContent( p+1, (+ratioerrorss[0][p])/ratioss[0][p] ); std::cout<<ratioerrorss[0][p]<<std::endl;}
//  for(unsigned int p=0; p<4; ++p){ h_sys9->SetBinContent( p+1, (-ratioerrorss[0][p])/ratioss[0][p] ); std::cout<<ratioss[0][p]<<std::endl;}
  for(unsigned int p=0; p<4; ++p){ h_sys10->SetBinContent( p+1, (ratioss[10][p]-ratioss[0][p])/ratioss[0][p] ); }


//TColor *mycolor = gROOT->TColor::GetColor("#33a02c");
  //// set colors
  //  
  //  _sbUp
  //  _sbDown
  //  _metUp
  //  _metDown
  //  _binUp
  //  _binDown
  //  _noPiso
  //  _eleTmpl
  Int_t ci1  = 1001;
  Int_t ci2  = 1002;
  Int_t ci3  = 1003;
  Int_t ci4  = 1004;
  Int_t ci5  = 1005;
  Int_t ci6  = 1006;
  Int_t ci7  = 1007;
  Int_t ci8  = 1008;
  Int_t ci9  = 1009;
  Int_t ci10 = 1010;
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

  h_sys0 ->SetMarkerColor(kBlack                     );                          
  h_sys1 ->SetMarkerColor(ci1);
  h_sys2 ->SetMarkerColor(ci2);
  h_sys3 ->SetMarkerColor(ci3);
  h_sys4 ->SetMarkerColor(ci4);
  h_sys5 ->SetMarkerColor(ci5);
  h_sys6 ->SetMarkerColor(ci6);
  h_sys7 ->SetMarkerColor(ci7);
  h_sys8 ->SetMarkerColor(ci8);
  h_sys9 ->SetMarkerColor(ci9);
  h_sys10->SetMarkerColor(ci10);

  h_sys0 ->SetLineColor(kBlack                     );                          
  h_sys1 ->SetLineColor(ci1);
  h_sys2 ->SetLineColor(ci2);
  h_sys3 ->SetLineColor(ci3);
  h_sys4 ->SetLineColor(ci4);
  h_sys5 ->SetLineColor(ci5);
  h_sys6 ->SetLineColor(ci6);
  h_sys7 ->SetLineColor(ci7);
  h_sys8 ->SetLineColor(ci8);
  h_sys9 ->SetLineColor(ci9);
  h_sys10->SetLineColor(ci10);

  h_sys0->SetLineWidth(3); 
  h_sys1->SetLineWidth(3); 
  h_sys2->SetLineWidth(3); 
  h_sys3->SetLineWidth(3); 
  h_sys4->SetLineWidth(3); 
  h_sys5->SetLineWidth(3); 
  h_sys6->SetLineWidth(3); 
  h_sys7->SetLineWidth(3); 
  h_sys8->SetLineWidth(3); 
  h_sys9->SetLineWidth(3); 
  h_sys10->SetLineWidth(3); 

// // sideband binning met gamma
//  TColor *color0 = new TColor( 7430, 0,0,0 ,       "color0" ); 
//  TColor *color1 = new TColor( 7431, 227,26,28  ,  "color1" ); 
//  TColor *color2 = new TColor( 7432, 251,154,153 , "color2" ); 
//  TColor *color3 = new TColor( 7433, 51,160,44 ,   "color3" ); 
//  TColor *color4 = new TColor( 7434, 178,223,138 , "color4" ); 
//  TColor *color5 = new TColor( 7435, 31,120,180 ,  "color5" ); 
//  TColor *color6 = new TColor( 7436, 166,206,227 , "color6" ); 
//  TColor *color7 = new TColor( 7437, 177,89,40 ,   "color7" ); 

  h_sys0->Draw("hist,same");
  h_sys1->Draw("hist,same");
  h_sys2->Draw("hist,same");
  h_sys3->Draw("hist,same");
  h_sys4->Draw("hist,same");
  h_sys5->Draw("hist,same");
  h_sys6->Draw("hist,same");
  h_sys7->Draw("hist,same");
//  h_sys8->Draw("hist,same");
  h_sys9->Draw("hist,same");
  h_sys10->Draw("hist,same");

  TLegend *leg4 = new TLegend(0.5,0.6,0.88,0.88 );
  leg4->SetFillColor(kWhite);
  leg4->AddEntry( h_sys0,"standard", "L");
  leg4->AddEntry( h_sys7,"bkg. no #gamma iso req", "L");
  leg4->AddEntry( h_sys1,"sideband Up", "L");
  leg4->AddEntry( h_sys2,"sideband Down", "L");
  leg4->AddEntry( h_sys3,"MET Up", "L");
  leg4->AddEntry( h_sys4,"MET Down", "L");
  leg4->AddEntry( h_sys5,"binning Up", "L");
  leg4->AddEntry( h_sys6,"binning Down","L");
//  leg4->AddEntry( h_sys8,"ele template", "L");
  leg4->AddEntry( h_sys9,"fit Up","L");
  leg4->AddEntry( h_sys10,"fit Down","L");
  leg4->Draw("same");

  //TLine *line0 = new TLine(xlow,xlow*ms[0]+bs[0],xhi,xhi*ms[0]+bs[0]);
  //TLine *line1 = new TLine(xlow,xlow*ms[1]+bs[1],xhi,xhi*ms[1]+bs[1]);
  //TLine *line2 = new TLine(xlow,xlow*ms[2]+bs[2],xhi,xhi*ms[2]+bs[2]);
  //TLine *line3 = new TLine(xlow,xlow*ms[3]+bs[3],xhi,xhi*ms[3]+bs[3]);
  //TLine *line4 = new TLine(xlow,xlow*ms[4]+bs[4],xhi,xhi*ms[4]+bs[4]);
  //TLine *line5 = new TLine(xlow,xlow*ms[5]+bs[5],xhi,xhi*ms[5]+bs[5]);
  //TLine *line6 = new TLine(xlow,xlow*ms[6]+bs[6],xhi,xhi*ms[6]+bs[6]);
  //TLine *line7 = new TLine(xlow,xlow*ms[7]+bs[7],xhi,xhi*ms[7]+bs[7]);

  //line0->SetLineColor(kBlack);   
  //line1->SetLineColor(kRed+1);   
  //line2->SetLineColor(kYellow-3);
  //line3->SetLineColor(kGreen+1); 
  //line4->SetLineColor(kGreen-9); 
  //line5->SetLineColor(kAzure+10);
  //line6->SetLineColor(kBlue-9);  
  //line7->SetLineColor(51);       

  //line1->Draw();
  //line2->Draw();
  //line3->Draw();
  //line4->Draw();
  //line5->Draw();
  //line6->Draw();
  //line7->Draw();
  //line0->Draw();

  c3->Update();
  
  //gr_ratio->Draw("sames,A*");

  c3->SaveAs(outpath+"/Graph_FakeRatiosRelativeHist"+extraname+".pdf");
  c3->SaveAs(wwwpath+"/Graph_FakeRatiosRelativeHist"+extraname+".pdf");

  c3->Clear();

//// reduced range

  TCanvas* c4 = new TCanvas("c4","c4",900,100,500,500);   

  gStyle->SetOptStat(0);
  //gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  gStyle->SetLineWidth(3);

  Double_t xalow = 175.;
  Double_t xahi = 400.;
  TH1F *hframe = c4->DrawFrame(xalow,-0.5,xahi,0.8,"");
  //TH1F *hframe = c4->DrawFrame(xlow,-0.3,xhi,0.8,"");
  hframe->SetXTitle("photon pT [GeV]");
  hframe->SetYTitle("Relative Fake Ratio"); 

  title->DrawTextNDC(0.17,0.87,"CMS");
  extra->DrawTextNDC(0.17,0.81,"Preliminary");
  lumi->DrawTextNDC(0.9,0.91,"2.32 /fb (13 TeV)");

  h_sys0->Draw("hist,same");
  h_sys1->Draw("hist,same");
  h_sys2->Draw("hist,same");
  h_sys3->Draw("hist,same");
  h_sys4->Draw("hist,same");
  h_sys5->Draw("hist,same");
  h_sys6->Draw("hist,same");
  h_sys7->Draw("hist,same");
  h_sys8->Draw("hist,same");
  h_sys9->Draw("hist,same");
  h_sys10->Draw("hist,same");

  leg4->Draw("same");

  //TLine *line0 = new TLine(xlow,xlow*ms[0]+bs[0],xhi,xhi*ms[0]+bs[0]);
  //TLine *line1 = new TLine(xlow,xlow*ms[1]+bs[1],xhi,xhi*ms[1]+bs[1]);
  //TLine *line2 = new TLine(xlow,xlow*ms[2]+bs[2],xhi,xhi*ms[2]+bs[2]);
  //TLine *line3 = new TLine(xlow,xlow*ms[3]+bs[3],xhi,xhi*ms[3]+bs[3]);
  //TLine *line4 = new TLine(xlow,xlow*ms[4]+bs[4],xhi,xhi*ms[4]+bs[4]);
  //TLine *line5 = new TLine(xlow,xlow*ms[5]+bs[5],xhi,xhi*ms[5]+bs[5]);
  //TLine *line6 = new TLine(xlow,xlow*ms[6]+bs[6],xhi,xhi*ms[6]+bs[6]);
  //TLine *line7 = new TLine(xlow,xlow*ms[7]+bs[7],xhi,xhi*ms[7]+bs[7]);

  //line0->SetLineColor(kBlack);   
  //line1->SetLineColor(kRed+1);   
  //line2->SetLineColor(kYellow-3);
  //line3->SetLineColor(kGreen+1); 
  //line4->SetLineColor(kGreen-9); 
  //line5->SetLineColor(kAzure+10);
  //line6->SetLineColor(kBlue-9);  
  //line7->SetLineColor(51);       

  //line1->Draw();
  //line2->Draw();
  //line3->Draw();
  //line4->Draw();
  //line5->Draw();
  //line6->Draw();
  //line7->Draw();
  //line0->Draw();

  c4->Update();
  
  //gr_ratio->Draw("sames,A*");

  c4->SaveAs(outpath+"/Graph_FakeRatiosRelativeHist_175to400"+extraname+".pdf");
  c4->SaveAs(wwwpath+"/Graph_FakeRatiosRelativeHist_175to400"+extraname+".pdf");

  c4->Clear();

 return;

}
