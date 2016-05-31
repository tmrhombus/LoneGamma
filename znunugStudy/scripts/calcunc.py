import math


#datas = [
#       5.000,
#       6.000,
#       5.000,
#       6.000,
#       6.000,
#       7.000,
#       6.000,
#       7.000,
#]
#
#ns = [
#       5.159,
#       5.442,
#       5.159,
#       5.442,
#       6.537,
#       7.052,
#       6.537,
#       7.052,
#]
#
#gens = [
#       11564.0,
#       12082.0,
#       11564.0,
#       12082.0,
#       14977.0,
#       16020.0,
#       14977.0,
#       16020.0,
#]
#
#rawrecos = [
#       14604.0,
#       15393.0,
#       17572.0,
#       18509.0,
#       18596.0,
#       20058.0,
#       22320.0,
#       24045.0,
#]
#



datas=[  16.000, 18.000, 16.000, 18.000, 18.000, 21.000, 23.000, 21.000, 23.000]
MC:    12.696, 13.382, 12.696, 13.382, 13.382, 16.166, 17.437, 16.166, 17.437]
Bkg:   0.347 , 0.398 , 0.347 , 0.398 , 0.398 , 0.541 , 0.649 , 0.541 , 0.64]
Gen:         23698.0, 24781.0, 23698.0, 24781.0, 30456.0, 32561.0, 30456.0, 32561.0]
Raw Reco:    14604.0, 15393.0, 17572.0, 18509.0, 18596.0, 20058.0, 22320.0, 24045.0]



#  dilepton

datas = [
 16.000,
 18.000,
 16.000,
 18.000,
 21.000,
 23.000,
 21.000,
 23.000,
]

sigs = [
 12.806,
 13.498,
 12.806,
 13.498,
 16.307,
 17.589,
 16.307,
 17.589,
]

gens = [
 23698.0,
 24781.0,
 23698.0,
 24781.0,
 30456.0,
 32561.0,
 30456.0,
 32561.0,
]

rawrecos = [
 14604.0,
 15393.0,
 17572.0,
 18509.0,
 18596.0,
 20058.0,
 22320.0,
 24045.0,
]

spaces = [
"MET < 170, 60 < dilep mass < 120, dPhi(photon,met) > 2",
"MET < 170, 60 < dilep mass < 120                      ",
"MET < 170, dPhi(photon,met) > 2                       ",
"MET < 170                                             ",
"MET < 110, 60 < dilep mass < 120, dPhi(photon,met) > 2",
"MET < 110, 60 < dilep mass < 120                      ",
"MET < 110, dPhi(photon,met) > 2                       ",
"MET < 110,                                            ",
]

#m170_ywnd_ydphi
#m170_ywnd_ndphi
#m170_nwnd_ydphi
#m170_nwnd_ndphi
#m110_ywnd_ydphi
#m110_ywnd_ndphi
#m110_nwnd_ydphi
#m110_nwnd_ndphi

for N,N_gen,N_reco,space,data in zip(ns,gens,rawrecos,spaces,datas):
 
 # define variables
 #N = 17.588 # 13.498  # 12.322 # = N_obs - N_bkg
 dN = math.sqrt(N)
 ddata = math.sqrt(data)
 
 R = 2.971
 dR = 0.001
 
 #N_gen =   32560.0 # 24781.000 # 23698.
 #N_reco =  20057.0 # 15393.000 # 14052.
 
 dN_gen  = math.sqrt(N_gen)
 dN_reco = math.sqrt(N_reco)
 
 
 #  acc * eff
 AE = N_reco/N_gen
 
 dAE = math.sqrt( 
 math.pow((1./N_gen),2)*math.pow((dN_reco),2) +
 math.pow((N_reco/(N_gen*N_gen)),2)*math.pow((dN_gen),2) 
 )
 print("\n\n%s"%space)
 print("Acceptance x Efficiency")
 print(AE)
 print(dAE)
 
 print '{0:.1f} +- {1:.1f}'.format(AE*100, dAE*100)
 
 # total nr
 #NR = ( N / AE ) * R
 NR = ( data / AE ) * R
 
 #dNR = math.sqrt(
 #math.pow((R/AE),2)*math.pow(dN,2) +
 #math.pow((N*R/(AE*AE)),2)*math.pow(dAE,2) +
 #math.pow((N/AE),2)*math.pow(dR,2)
 #)
 
 dNR = math.sqrt(
 math.pow((R/AE),2)*math.pow(ddata,2) +
 math.pow((data*R/(AE*AE)),2)*math.pow(dAE,2) +
 math.pow((data/AE),2)*math.pow(dR,2)
 )
 
 print("Total Number of events")
 print(NR)
 print(dNR)
 print '{0:.1f} +- {1:.1f}'.format(NR, dNR)


 print("{0:s} & {1:0.1f} & {2:0.1f} $\pm$ {3:0.1f} & {4:0.1f} $\pm$ {5:0.1f} & {6:0.1f} $\pm$ {7:0.1f} \\\\ ".format(space,data,N,dN,AE*100,dAE*100,NR,dNR))
