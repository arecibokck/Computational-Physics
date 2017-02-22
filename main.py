from PIMC import *


#Initialize
M = 100 #Number of time slices
T = 10 #Imaginary time period
a = 0.1 #Lattice Spacing
n_steps = 1000  #Number of Monte Carlo steps
x_max = 5
n_bins = 100

path = Path(M,T,a,n_steps,x_max,n_bins)

if __name__=='__main__':
    PIMC(path)
