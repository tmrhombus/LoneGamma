

// matched
double e_mch = h_Nsig_chiso_175to1000_idnc->GetEntries()
double i_mch = h_Nsig_chiso_175to1000_idnc->Integral()

double e_mch_trg = h_Nsig_chiso_175to1000_idnc_trig->GetEntries()
double i_mch_trg = h_Nsig_chiso_175to1000_idnc_trig->Integral()

double e_mch_met = h_Nsig_chiso_175to1000_idnc_mL30->GetEntries()
double i_mch_met = h_Nsig_chiso_175to1000_idnc_mL30->Integral()

double e_mch_met_trg = h_Nsig_chiso_175to1000_idnc_mL30_trig->GetEntries()
double i_mch_met_trg = h_Nsig_chiso_175to1000_idnc_mL30_trig->Integral()


// unmatched
double e_nmh = h_Nbkg_chiso_175to1000_idnc->GetEntries()
double i_nmh = h_Nbkg_chiso_175to1000_idnc->Integral()

double e_nmh_trg = h_Nbkg_chiso_175to1000_idnc_trig->GetEntries()
double i_nmh_trg = h_Nbkg_chiso_175to1000_idnc_trig->Integral()

double e_nmh_met = h_Nbkg_chiso_175to1000_idnc_mL30->GetEntries()
double i_nmh_met = h_Nbkg_chiso_175to1000_idnc_mL30->Integral()

double e_nmh_met_trg = h_Nbkg_chiso_175to1000_idnc_mL30_trig->GetEntries()
double i_nmh_met_trg = h_Nbkg_chiso_175to1000_idnc_mL30_trig->Integral()


// combined hopefully
double e_cmb = h_Deno_chiso_175to1000_idnc->GetEntries()
double i_cmb = h_Deno_chiso_175to1000_idnc->Integral()

double e_cmb_trg = h_Deno_chiso_175to1000_idnc_trig->GetEntries()
double i_cmb_trg = h_Deno_chiso_175to1000_idnc_trig->Integral()

double e_cmb_met = h_Deno_chiso_175to1000_idnc_mL30->GetEntries()
double i_cmb_met = h_Deno_chiso_175to1000_idnc_mL30->Integral()

double e_cmb_met_trg = h_Deno_chiso_175to1000_idnc_mL30_trig->GetEntries()
double i_cmb_met_trg = h_Deno_chiso_175to1000_idnc_mL30_trig->Integral()



printf("ID\n");
printf( " matched    %6.1f  %6.1f  \n" , e_mch, i_mch  ) ;
printf( " unmatched  %6.1f  %6.1f  \n" , e_nmh, i_nmh  ) ;
printf( " combined   %6.1f  %6.1f  \n" , e_cmb, i_cmb  ) ;

printf("ID+Trigger\n");
printf( " matched    %6.1f  %6.1f  \n" , e_mch_trg, i_mch_trg  ) ;
printf( " unmatched  %6.1f  %6.1f  \n" , e_nmh_trg, i_nmh_trg  ) ;
printf( " combined   %6.1f  %6.1f  \n" , e_cmb_trg, i_cmb_trg  ) ;

printf("ID+MET\n");
printf( " matched    %6.1f  %6.1f  \n" , e_mch_met, i_mch_met  ) ;
printf( " unmatched  %6.1f  %6.1f  \n" , e_nmh_met, i_nmh_met  ) ;
printf( " combined   %6.1f  %6.1f  \n" , e_cmb_met, i_cmb_met  ) ;

printf("ID+Trigger+MET\n");
printf( " matched    %6.1f  %6.1f  \n" , e_mch_met_trg, i_mch_met_trg  ) ;
printf( " unmatched  %6.1f  %6.1f  \n" , e_nmh_met_trg, i_nmh_met_trg  ) ;
printf( " combined   %6.1f  %6.1f  \n" , e_cmb_met_trg, i_cmb_met_trg  ) ;


QCD Sample ( Unweighted Entries -- Integral Events )
ID
 matched    16272.0  36954.3  
 unmatched  13824.0  42443.9  
 combined   30096.0  79398.2  
ID+Trigger
 matched    16111.0  36687.6  
 unmatched  3454.0   18863.1  
 combined   19565.0  55550.6  
ID+MET
 matched    7152.0  19442.3  
 unmatched  6023.0  20065.3  
 combined   13175.0 39507.6  
ID+Trigger+MET
 matched    7091.0  19294.1  
 unmatched  1608.0  7072.2  
 combined   8699.0  26366.3  


GJ Sample ( Unweighted Entries -- Integral Events )
ID
 matched    265664.0  97567.7  
 unmatched  6884.0    583.1  
 combined   272548.0  98150.8  
ID+Trigger
 matched    262026.0  96978.3  
 unmatched  2023.0    290.0  
 combined   264049.0  97268.2  
ID+MET
 matched    131816.0  54558.6  
 unmatched  2920.0    283.4  
 combined   134736.0  54842.0  
ID+Trigger+MET
 matched    130374.0  54275.3  
 unmatched  867.0     140.5  
 combined   131241.0  54415.8  

