from HO_PIMC import *
from AHO_PIMC import *

#Initialize
M = 100 #Number of time slices
T = 10 #Imaginary time period
delta = 1.0 #Metropolis step size
nsteps = 1000  #Number of Monte Carlo steps

HOpath = HOPath(M,T,delta)
AHOpath = AHOPath(M,T,delta)

HOPIMC(nsteps,HOpath)
