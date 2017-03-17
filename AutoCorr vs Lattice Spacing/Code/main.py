from PIMC import *
from subprocess import check_output
from scipy.optimize import curve_fit
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
	n = 0
	slopes = []

	for i in Trange:

		n+=1

		path = Path(M,i,delta,m,mu,f,char,thermalize,lam,n_steps,x_max,n_bins)
		msp, energy, ap = PIMC(path, plots = False)

		data_list = [msp,energy.tolist()]
		with open(filename, "wb") as fi:
			writer = csv.writer(fi)
			writer.writerows(data_list)

		cmd = 'Rscript'
		path2script = './data_analysis.R'
		args = [filename, str(n)]

		rcall = [cmd, path2script] + args

		output = check_output(rcall)

		if output:

		    slopes.append(output.split()[1])

		else:

			slopes.append("0")

	v=np.array([-1.0/float(i) for i in slopes])

	x = np.divide(Trange,float(M))
	y = v
	y[y==0.0] = float('nan')
	se = [np.std(v)]*len(v)
	f = plt.figure(1,figsize=(10, 5), dpi=80)
	fs = gridspec.GridSpec(1, 2, width_ratios=[1,1])
	fs.update(wspace = 0.3)
	plt.subplot(fs[0])
	plt.scatter(x,y, label='Data')
	plt.errorbar(x,y,yerr=se, linestyle="None")
	plt.hold('on')
	idx = np.isfinite(x) & np.isfinite(y)
	f = lambda z, *p: p[0] * z**p[1]
	popt, pcov = curve_fit(f, [x[i] for i,j in enumerate(idx) if idx[i]], [y[i] for i,j in enumerate(idx) if idx[i]], [0.01,-0.04])
	xfine = np.linspace(0.1, max(x), 100)  # define values to plot the function for
	plt.plot(xfine, f(xfine, *popt), 'r-', label='Fit')
	plt.legend(prop={'size':10})
	plt.xlabel("Lattice Spacing")
	plt.ylabel("Autocorrelation Time")
	plt.title('AC Time vs LS')
	plt.grid(True)

	lx = np.log(x)
	ly = np.log(v)
	se = [np.std(ly)]*len(ly)
	plt.subplot(fs[1])
	plt.xlabel("Log(Lattice Spacing)")
	plt.ylabel("Log(Autocorrelation Time)")
	plt.title("Log(AC Time) vs Log(LS)")
	plt.grid(True)
	plt.legend(prop={'size':10})
	plt.scatter(lx,ly, label='Data')
	plt.errorbar(lx,ly,yerr=se, linestyle="None")
	idx = np.isfinite(x) & np.isfinite(y)
	f = lambda z, *p: p[0] + z*p[1]
	popt, pcov = curve_fit(f, [lx[i] for i,j in enumerate(idx) if idx[i]], [ly[i] for i,j in enumerate(idx) if idx[i]], [-1.0,-0.25])
	xfine = np.linspace(min(lx),0.0, 100)  # define values to plot the function for
	plt.plot(xfine, f(xfine, *popt), 'r-', label='Fit')
	plt.legend(prop={'size':10})
	plt.hold('off')
	plt.show()
	print("Finished! - All tasks successfully completed!")
