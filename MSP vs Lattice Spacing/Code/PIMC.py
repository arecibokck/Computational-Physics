import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class Path(object):

    def __init__(self,M,T,delta,m,mu,f,char,thermalize,lam,n_steps,x_max,n_bins): #Initialize

        self.M = M
        self.T = T
        self.dT = T/M
        self.delta = delta
        self.m = m
        self.mu = mu
        self.f = f
        self.thermalize = thermalize
        self.lam = lam
        self.n_steps = n_steps
        self.x_max = x_max
        self.x_min = -x_max
        self.n_bins = n_bins
        self.n = 0.0

        if char=='h': #Initialize a position configuration
            self.X = (np.random.uniform(size=(1, self.M)).tolist())[0] #Hot Start
        elif char=='c':
            self.X = np.zeros(self.M, dtype=np.int).tolist() #Cold Start

    def V(self,x):

        return 0.5*(self.mu**2)*(x**2) + self.lam*((x**2-self.f**2)**2) # Potential Energy(PE)
                                                                        # lam = 0 -> Harmonic Oscillator
                                                                        # lam > 0 -> Anharmonic Oscillator

    def DV(self,x):

        return (self.mu**2)*x + 4*self.lam*x*(x**2-self.f**2) # First derivative of Potential Energy

    def K(self,x):

        return 0.5*self.m*(x**2)*(1/self.dT**2) # Kinetic Energy(KE) with m = 1

    def E(self,x):

        return self.V(x) + (0.5*x*self.DV(x)) # Energy from Virial Theorem - V+(0.5*x*dV/dx)


def mcstep(path): # One Markov Chain Monte Carlo Step - Metropolis Algorithm

    k = int(path.M * np.random.uniform()) # Pick an arbitrary index
    x_p = path.X[k] + np.random.uniform(-path.delta,path.delta) # Pick the position corresponding to random index and..
    k_p = k + 1                                                  # ..generate new value by adding a small arbitrary change in value, delta_x.
    if (k_p > path.M-1): k_p = 0
    k_m = k - 1
    if (k_m < 0): k_m = path.M-1
    dVa = path.V(x_p) - path.V(path.X[k]) # Potential Action
    dKa = (path.K(path.X[k_p]-x_p) + path.K(x_p - path.X[k_m])) - (path.K(path.X[k_p] - path.X[k]) + path.K(path.X[k] - path.X[k_m])) #Kinetic Action
    dS = path.dT*(dVa + dKa) # Total Action for Harmonic Oscillator
    if(dS < 0.0 or np.random.normal() < np.exp(-dS)): # Accept-Reject step
        x_new = path.X[k] = x_p
        path.n+=1
    else:
        x_new = path.X[k]

    return x_new # New position as determined by the acceptance probability


def render_plots(msp, x, pdf):

    # Plot of mean squared values of position
    f = plt.figure(1,figsize=(10, 5), dpi=80)
    fs = gridspec.GridSpec(1, 2, width_ratios=[1,1])
    fs.update(wspace = 0.3)
    plt.subplot(fs[0])
    plt.xlabel('No Of Iterations')
    plt.ylabel(r'$\mathrm{<x^2>}$')
    plt.title(r'$\mathrm{Mean\ Square\ of\ Position}$')
    plt.axis([0,len(msp),0,max(msp)])
    plt.grid(True)
    m = [np.mean(msp)]*len(msp)
    plt.plot(msp)
    plt.hold('on')
    plt.plot(m, label='Mean')
    plt.legend(prop={'size':10})

    # Plot of Probability density function
    plt.subplot(fs[1])
    plt.xlabel('Position in x')
    plt.ylabel('Probability')
    plt.title(r'$\mathrm{Probability\ Density\ over\ Position}$')
    plt.hold('off')
    plt.plot(x, pdf, 'r-', label='Gaussian \n KDE Fit')
    plt.grid(True)
    plt.legend(prop={'size':10})
    plt.show()


def PIMC(path,plots):

    # Thermalization Steps
    if(path.thermalize):
        for i in range(path.n_steps/2): # Run over a few steps for the first time to pass burn-in phase
            for j in range(path.M): # Run over the full time = dT*M
                mcstep(path)

    # Production Steps
    msp = []
    energy = []
    ensemble_pdf = [0]*path.n_bins
    for i in range(path.n_steps): # Run over the full range of steps to obtain energy values that can be worked with
        msp.append(np.mean(np.square(path.X))/path.T)
        for j in range(path.M): # Run over the full time = dT*M
            x_new = mcstep(path)
            bins = abs(int((x_new - path.x_min)/ (path.x_max - path.x_min) * path.n_bins))
            if(bins<path.n_bins):ensemble_pdf[bins]+=1
            energy.append(path.E(x_new))

    if(plots == True):

        print("Rendering Plots...")
        x = (np.arange(path.x_min,path.x_max,(path.x_max-path.x_min)/path.n_bins)).tolist()
        a = []
        for i in range(len(ensemble_pdf)):
            a.append([x[i]]*ensemble_pdf[i])
        y = np.array([j for i in a for j in i])
        kde = gaussian_kde(y, bw_method=0.2 / y.std(ddof=1)) # Kernel Density Estimation to determine PDF
        kde.covariance_factor = lambda : .25
        kde._compute_covariance()
        pdf = kde.evaluate(x) # PDF of the ensemble of paths
        render_plots(msp, x, pdf)

    return msp, np.divide(energy,path.T), path.n/(path.n_steps*path.M)
