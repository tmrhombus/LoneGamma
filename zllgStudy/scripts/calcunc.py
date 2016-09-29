import math

 ## 2015 values ##
 # spaces = [
 # "MET < 170, 60 < dilep mass < 120, dPhi(photon,met) > 2",
 # "MET < 170, 60 < dilep mass < 120                      ",
 # "MET < 170, dPhi(photon,met) > 2                       ",
 # "MET < 170                                             ",
 # "MET < 110, 60 < dilep mass < 120, dPhi(photon,met) > 2",
 # "MET < 110, 60 < dilep mass < 120                      ",
 # "MET < 110, dPhi(photon,met) > 2                       ",
 # "MET < 110,                                            ",
 # ]
 # 
 # 
 # datas=  [16.000, 18.000, 16.000, 18.000, 18.000, 21.000, 23.000, 21.000, 23.000]
 # sigs=   [12.696, 13.382, 12.696, 13.382, 13.382, 16.166, 17.437, 16.166, 17.437]
 # bkgs=   [0.347 , 0.398 , 0.347 , 0.398 , 0.398 , 0.541 , 0.649 , 0.541 , 0.64]
 # gens=      [23698.0, 24781.0, 23698.0, 24781.0, 30456.0, 32561.0, 30456.0, 32561.0]
 # rawrecos=  [14604.0, 15393.0, 17572.0, 18509.0, 18596.0, 20058.0, 22320.0, 24045.0]

# using 2015 A*E for no window, yes dphi

spaces = [
 "MET < 170, dPhi(photon,met) > 2  dilep mass > 20 : prefit",
 "MET < 170, dPhi(photon,met) > 2  dilep mass > 20 : fitted",
 ]
 
 
datas=  [ 16.000, 16.000 ]
sigs=   [ 12.310, 15.254 ]
bkgs=   [ 0.448 , 0.467]
gens=      [ 23698.0, 23698.0 ]
rawrecos=  [ 17572.0, 17572.0 ]

for N_data,N_sig,N_bkg,N_gen,N_raw,space in zip(datas,sigs,bkgs,gens,rawrecos,spaces):

 # Expression is
 #
 #  ( N_data - N_bkg ) x  R(Z(nunu)/Z(ll)) / ( Acc x Eff )

 N = N_data - N_bkg
 dN = math.sqrt(N)
 
 R = 2.971
 dR = 0.001

 #  acc * eff
 AE = N_raw/N_gen
 
  # for AxE 
 dN_gen = math.sqrt(N_gen)
 dN_raw = math.sqrt(N_raw)

 dAE = math.sqrt( 
 math.pow((1./N_gen),2)*math.pow((dN_raw),2) +
 math.pow((N_raw/(N_gen*N_gen)),2)*math.pow((dN_gen),2) 
 )
 print("\n\n%s"%space)
 print("Acceptance x Efficiency")
 print(AE)
 print(dAE)
 
 print '{0:.1f} +- {1:.1f}'.format(AE*100, dAE*100)
 

 # Predicted Yield 
 NR = ( N / AE ) * R
 
 dNR = math.sqrt(
 math.pow((R/AE),2)*math.pow(dN,2) +
 math.pow((N*R/(AE*AE)),2)*math.pow(dAE,2) +
 math.pow((N/AE),2)*math.pow(dR,2)
 )
 
 
 print("Total Number of events")
 print(NR)
 print(dNR)
 print '{0:.1f} +- {1:.1f}'.format(NR, dNR)


 print("{0:s} & {1:0.1f} & {2:0.1f} $\pm$ {3:0.1f} & {4:0.1f} $\pm$ {5:0.1f} & {6:0.1f} $\pm$ {7:0.1f} \\\\ ".format(space,N_data-N_bkg,N,dN,AE*100,dAE*100,NR,dNR))

