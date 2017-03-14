from PIMC import *
import subprocess
import csv

#Initialize
M = 100 #Number of time slices
T = 10.0 #Imaginary time period
delta = 0.7 #Deflection in position
m = 1.0 #Mass
mu = 1 #Oscillator Frequency?
f = 4.0 #Arbitrary Constant? for Anharmonic Oscillator
char = 'c' #Choose Hot = 'h' or Cold = 'c' Start
thermalize = False #Set True to run Thermalization to pass burn-in phase
lam = 0.0 #Quartic Coupling Constant - lam = 0 -> Harmonic Oscillator; lam > 0 -> Anharmonic Oscillator
n_steps = 1000  #Number of Monte Carlo steps
x_max = 8.0 #Max value corresponding to the last bin
n_bins = 100 #Number of bins to draw the probability density histogram over

if __name__=='__main__':

    filename = ""

    if(lam==0):
        print("Running MCMC for the Harmonic Oscillator with set parameters...")
        filename = "PIMC_data_Harmonic.csv"
    else:
        print("Running MCMC for the Anharmonic Oscillator with set parameters...")
        filename = "PIMC_data_Anharmonic.csv"

    path = Path(M,T,delta,m,mu,f,char,thermalize,lam,n_steps,x_max,n_bins)
    msp, energy, ap = PIMC(path, plots = True)

    print("Dumping all data into a CSV file and running data analysis in R...")
    data_list = [msp,energy.tolist()]
    with open(filename, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(data_list)

    cmd = 'Rscript'
    path2script = './data_analysis.R'
    args = filename

    rcall = [cmd, path2script, args]

    output = subprocess.check_output(rcall)

    print("Finished! - All tasks successfully completed!")
