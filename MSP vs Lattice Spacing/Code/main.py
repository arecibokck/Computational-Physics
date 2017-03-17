from PIMC import *
from scipy.optimize import curve_fit
from subprocess import check_output
import numpy as np
import matplotlib.pyplot as plt
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
		ilename = "PIMC_data_Anharmonic.csv"

	Trange = (np.arange(10.0,70.0,1.0)).tolist()
	means = []
	ls = np.divide(Trange,float(M))
	se = []

	for i in Trange:

		path = Path(M,i,delta,m,mu,f,char,thermalize,lam,n_steps,x_max,n_bins)
		msp, energy, ap = PIMC(path, plots = False)
		means.append(np.mean(msp))
		cmd = 'Rscript'
		path2script = './data_analysis.R'
		args = " ".join(str(e) for e in msp)
		rcall = [cmd, path2script, args]
		output = check_output(rcall)
		se.append(float(output.split()[1]))

	print("Dumping all data into a CSV file and rendering plots...")
	data_list = [ls, means, se]
	with open(filename, "wb") as fi:
		writer = csv.writer(fi)
		writer.writerows(data_list)
	plt.scatter(ls,means)
	plt.errorbar(ls,means,yerr=se, linestyle="None")
	plt.hold('on')
	f = lambda x, *p: p[0] * x**p[1]
	popt, pcov = curve_fit(f, ls, means, [0.01,-0.04])
	xfine = np.linspace(0.1, max(ls), 100)  # define values to plot the function for
	plt.plot(xfine, f(xfine, *popt), 'r-', label='Fit')
	plt.legend(prop={'size':10})
	plt.xlabel("Lattice Spacing")
	plt.ylabel(r"$\mathrm{<x^2>}$")
	plt.ylim(0.0,max(ls))
	plt.grid(True)
	plt.hold('off')
	plt.show()
	print("Finished! - All tasks successfully completed!")
