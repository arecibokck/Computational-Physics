from PIMC import *
from subprocess import check_output
import numpy as np
import matplotlib.pyplot as plt
import csv

#Initialize
M = 100 #Number of time slices
T = 10.0 #Imaginary time period
delta = 0.7 #Deflection in position
m = 1.0 #Mass
mu = 1.0 #Oscillator Frequency?
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

	deltarange = (np.arange(0.1,50.1,0.5)).tolist()
	slopes = []

	for i in deltarange:

		path = Path(M,T,i,m,mu,f,char,thermalize,lam,n_steps,x_max,n_bins)
		msp, energy, ap = PIMC(path, plots = False)

		data_list = [msp]
		with open(filename, "wb") as fi:
			writer = csv.writer(fi)
			writer.writerows(data_list)

		cmd = 'Rscript'
		path2script = './data_analysis.R'
		args = [filename]

		rcall = [cmd, path2script] + args

		output = check_output(rcall)

		if output:

		    slopes.append(output.split()[1])

		else:

			slopes.append("0")

	v=np.array([-1.0/float(i) for i in slopes])
	x = deltarange
	y = v
	y[y==0.0] = float('nan')
	se = [np.std(v)]*len(v)
	plt.plot(x,y, label='Data')
	plt.errorbar(x,y,yerr=se, linestyle="None")
	plt.legend(prop={'size':10})
	plt.xlabel(r"$\mathrm{\Delta}$")
	plt.ylabel("Autocorrelation Time")
	plt.title('AC Time vs Delta')
	plt.grid(True)
	plt.xlim(0.0,max(deltarange))
	plt.show()
	print("Finished! - All tasks successfully completed!")
