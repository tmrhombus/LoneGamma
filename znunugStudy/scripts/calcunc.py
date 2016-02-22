import math

# define variables
N = 13.8 # MC weighted
dN = math.sqrt(N)

R = 2.971
dR = 0.001

# random lepton
#N_gen =  27598.
#N_reco = 16299.

N_gen =  27604.
N_reco = 16299.

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
