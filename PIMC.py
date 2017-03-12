import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv

class Path(object):

    def __init__(self,M,T,delta_x,m,mu,f,char,thermalize,lam,n_steps,x_max,n_bins): #Initialize

        self.M = M
        self.T = T
        self.dT = T/M
        self.delta_x = delta_x
        self.m = m
        self.mu = mu
        self.f = f
        self.thermalize = thermalize
        self.lam = lam
        self.n_steps = n_steps
        self.x_max = x_max
        self.x_min = -x_max
        self.n_bins = n_bins

        if(self.lam==0):
            print("HARMONIC OSCILLATOR")
        else:
            print("ANHARMONIC OSCILLATOR")

        if char=='h': #Initialize a position configuration
            print('Initializing for a hot start...')
            self.X = (np.random.normal(size=(1, self.M)).tolist())[0] #Hot Start
        elif char=='c':
            print('Initializing for a cold start...')
            self.X = np.zeros(self.M, dtype=np.int).tolist() #Cold Start

    def V(self,x):

        return 0.5*(self.mu**2)*(x**2) + self.lam*((x**2-self.f**2)**2) #Potential Energy(PE)
                                                                        #lam = 0 -> Harmonic Oscillator
                                                                        #lam > 0 -> Anharmonic Oscillator

    def DV(self,x):

        return (self.mu**2)*x + 4*self.lam*x*(x**2-self.f**2)

    def K(self,x):

        return 0.5*self.m*(x**2)*(1/self.dT**2) #Kinetic Energy(KE) with m = 1

    def E(self,x):

        return self.V(x) + (0.5*x*self.DV(x)) # Energy from Virial Theorem - V+(0.5*x*dV/dx)


def mcstep(path): #One Markov Chain Monte Carlo Step - Metropolis Algorithm

    k = int(path.M * np.random.uniform()) #Pick an arbitrary index
    x_p = path.X[k] + np.random.uniform(-path.delta_x,path.delta_x) #Pick the position corresponding to random index and..
    k_p = k + 1                                                  #..generate new value by adding a small arbitrary change in value, delta_x.
    if (k_p > path.M-1): k_p = 0
    k_m = k - 1
    if (k_m < 0): k_m = path.M-1
    dVa = path.V(x_p) - path.V(path.X[k]) #Potential Action
    dKa = (path.K(path.X[k_p]-x_p) + path.K(x_p - path.X[k_m])) - (path.K(path.X[k_p] - path.X[k]) + path.K(path.X[k] - path.X[k_m])) #Kinetic Action
    dS = path.dT*(dVa + dKa) # Total Action for Harmonic Oscillator
    if(dS < 0.0 or np.random.normal(0,1) < np.exp(-dS)): #Accept-Reject step
        x_new = path.X[k] = x_p
    else:
        x_new = path.X[k]

    return x_new #New position as determined by the acceptance probability

def gaussian(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def render_plots(path,normed_pdf,ms):

    #Plot of MS values of position
    f = plt.figure(1,figsize=(10, 5), dpi=80)
    fs = gridspec.GridSpec(1, 2, width_ratios=[1,1])
    fs.update(wspace = 0.3)
    plt.subplot(fs[0])
    plt.xlabel('No Of Iterations')
    plt.ylabel(r'$\mathrm{<x^2>}$')
    plt.title(r'$\mathrm{Mean\ Square\ of\ Position}$')
    plt.axis([0,len(ms),0,max(ms)])
    plt.grid(True)
    m = [np.mean(ms)]*len(ms)
    plt.plot(ms)
    plt.hold('on')
    plt.plot(m)

    #Histogram of the probability density
    plt.subplot(fs[1])
    plt.xlabel('Position in x')
    plt.ylabel('Probability')
    plt.title(r'$\mathrm{Probability\ Density\ over\ Position}$')
    plt.hold('off')
    x = np.linspace(path.x_min, path.x_max, path.n_bins) #Range of position values between which path was generated
    y = np.array(normed_pdf)
    plt.scatter(x,y,s=5, label='Sim Data')
    plt.hold('on')
    if (path.lam == 0):
        popt, pcov = curve_fit(gaussian, x, y)
        plt.plot(x, gaussian(x, *popt), 'r-', label='fit')
        plt.legend(prop={'size':10})
    else:
        plt.plot(x, y, 'r-', label='fit')
        plt.legend(prop={'size':8})
    plt.hold('off')
    plt.grid(True)
    plt.ylim((0,(max(y)+0.005)))
    plt.show()

def dump_data(normed_pdf,ms,energy):

    with open('PIMC_data.csv', 'w') as fp:
        a = csv.writer(fp, delimiter=',')
        data = [normed_pdf, ms, energy]
        a.writerows(data)

def PIMC(path):

    if(path.thermalize): #Thermalize
        print("Running Thermalization Steps...")
        for i in range(path.n_steps/2): #Run over a few steps for the first time to pass burn-in phase
            for j in range(path.M): # Run over the full time = dT*M
                mcstep(path)

    #Production Steps
    print("Running Production Steps...")
    ms = []
    energy = []
    ensemble_pdf = [0]*path.n_bins
    for i in range(path.n_steps): # Run over the full range of steps to obtain energy values that can be worked with
        ms.append(np.mean(np.square(path.X))/path.T)
        for j in range(path.M): #Run over the full time = dT*M
            x_new = mcstep(path)
            bins = abs(int((x_new - path.x_min)/ (path.x_max - path.x_min) * path.n_bins))
            if(bins<path.n_bins):ensemble_pdf[bins]+=1
            energy.append(path.E(x_new))
    normed_pdf = [float(i)/sum(ensemble_pdf) for i in ensemble_pdf] #Normalized pdf of the ensemble of paths-sum(ensemble_pdf)=n_steps*M
    print("Rendering Plots...")
    render_plots(path,normed_pdf,ms)
    print("Dumping all data into a CSV file...")
    dump_data(normed_pdf,ms,energy)
    print("Finished! - All tasks successfully completed!")
