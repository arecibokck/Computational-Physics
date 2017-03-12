from PIMC import *

#Initialize
M = 100 #Number of time slices
T = 10.0 #Imaginary time period
delta_x = 0.7 #Deflection in position
m = 1.0 #Mass
mu = 1 #Oscillator Frequency?
f = 4.0 #Arbitrary Constant? for Anharmonic Oscillator
char = 'c' #Choose Hot = 'h' or Cold = 'c' Start
thermalize = False #Set True to run Thermalization to pass burn-in phase
lam = 1.0 #Quartic Coupling Constant - lam = 0 -> Harmonic Oscillator; lam > 0 -> Anharmonic Oscillator
n_steps = 1000  #Number of Monte Carlo steps
x_max = 8.0 #Max value corresponding to the last bin
n_bins = 100 #Number of bins to draw the probability density histogram over

path = Path(M,T,delta_x,m,mu,f,char,thermalize,lam,n_steps,x_max,n_bins)

if __name__=='__main__':
    PIMC(path)
