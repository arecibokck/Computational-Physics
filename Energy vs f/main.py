from PIMC import *
from subprocess import check_output
import numpy as np
import matplotlib.pyplot as plt
import csv

#Initialize
M = 100 #Number of time slices
T = 10.0 #Imaginary time period
delta = 0.7 #Deflection in position
m = 0.5 #Mass
mu = 1 #Oscillator Frequency?
fs = 1.0 #Arbitrary Constant? (SQUARED) for Anharmonic Oscillator
char = 'c' #Choose Hot = 'h' or Cold = 'c' Start
thermalize = False #Set True to run Thermalization to pass burn-in phase
lam = 1.0 #Quartic Coupling Constant - lam = 0 -> Harmonic Oscillator; lam > 0 -> Anharmonic Oscillator
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

	fsrange = np.arange(-2.0,5.0,1.0)
	means = []
	se = []

	for i in fsrange:

		path = Path(M,T,delta,m,mu,i,char,thermalize,lam,n_steps,x_max,n_bins)
		msp, energy, ap = PIMC(path, plots = False)
		means.append(np.mean(energy))

		data_list = [energy.tolist()]
		with open(filename, "wb") as fi:
			writer = csv.writer(fi)
			writer.writerows(data_list)

		cmd = 'Rscript'
		path2script = './data_analysis.R'
		args = [filename]

		rcall = [cmd, path2script] + args

		output = check_output(rcall)

		se.append(float(output.split()[1]))

	print("Rendering plots...")
	plt.plot(fsrange, means, label='Data')
	plt.errorbar(fsrange, means, yerr=se, linestyle="None")
	plt.hold('on')
	plt.legend(prop={'size':10})
	plt.xlabel(r"$\mathrm{f^{2}}$")
	plt.ylabel("Energy")
	x = (np.arange(-2.0, 5.0, 1.0)).tolist()
	y = [2.677,1.060,1.137,2.289,3.251,3.863,4.366]
	plt.plot(x, y, 'r-', label='Actual')
	#plt.ylim(0.0,max(ls))
	plt.legend(prop={'size':10})
	plt.grid(True)
	plt.hold('off')
	plt.show()
	print("Finished! - All tasks successfully completed!")
