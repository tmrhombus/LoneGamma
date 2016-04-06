import math

# define variables
N = 12.322 # = N_obs - N_bkg
dN = math.sqrt(N)

R = 2.971
dR = 0.001

N_gen =   23698.
N_reco =  14052.

dN_gen  = math.sqrt(N_gen)
dN_reco = math.sqrt(N_reco)


#  acc * eff
AE = N_reco/N_gen

#dAE =( math.sqrt(N_reco) / N_gen) 

dAE = math.sqrt( 
math.pow((1./N_gen),2)*math.pow((dN_reco),2) +
math.pow((N_reco/(N_gen*N_gen)),2)*math.pow((dN_gen),2) 
)
print("Acceptance x Efficiency")
print(AE)
print(dAE)

print '{0:.1f} +- {1:.1f}'.format(AE*100, dAE*100)

# total nr
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
