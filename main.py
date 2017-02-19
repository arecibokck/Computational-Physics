from QHO_PIMC import *

#Initialize
M = 100 #Number of time slices
T = 10 #Imaginary time period
delta = 1.0 #Metropolis step size
nsteps = 300  #Number of Monte Carlo steps

path = Path(M,T,delta)

PIMC(nsteps,path)
