import numpy as np
import matplotlib.pyplot as plt
import csv

class Path(object):

    def __init__(self,M,T,delta):
        #Initialize
        self.M = M #Number of time slices
        self.T = T #Imaginary time period
        self.dT = T/M #Time Step
        self.delta = delta #Metropolis step size

        #Initialize a position configuration
        #self.X = (np.random.uniform(size=(1, self.M)).tolist())[0]
        self.X = np.zeros(self.M, dtype=np.int).tolist()

    def V(self,x):
        return 0.5*(x**2) #Potential Energy(PE) with m = 1, f = 1

    def K(self,x):
        return 0.5*(x**2) #Kinetic Energy(KE) with m = 1

    def E(self,x):
        return self.V(x) + self.K(x) # Energy from Virial Theorem - V+(0.5*x*dV/dx) which gives Total Energy = PE + KE


def mcstep(path): #One Markov Chain Monte Carlo Step - Metropolis Algorithm

    k = int(path.M * np.random.uniform()) #Pick an arbitrary index
    x_p = path.X[k] + (2 * np.random.uniform() - 1) * path.delta #Pick the position corresponding to arbitrary index and..
    k_p = k + 1                                                  #..generate new value by adding delta that increases or reduces the value
    if (k_p > path.M-1): k_p = 0
    k_m = k - 1
    if (k_m < 0): k_m = path.M-1
    v = path.V(x_p) - path.V(path.X[k]) #Potential Action
    k = (path.K(path.X[k_p]-x_p) + path.K(x_p - path.X[k_m])) - (path.K(path.X[k_p] - path.X[k]) + path.K(path.X[k] - path.X[k_m])) #Kinetic Action
    S = v + k # Total Action for Harmonic Oscillator
    if(dS < 0.0 or np.random.uniform(0,1) < np.exp(-dS*path.dT)): #Accept-Reject step
        x_new = path.X[k] = x_p
    else:
        x_new = path.X[k]

    return x_new #New position as determined by the acceptance probability


def render_plots(final_path, rms, energy):

    #the histogram of the data
#    n, bins, patches = plt.hist(final_path, 50, normed=1, facecolor='green', alpha=0.75)

#    plt.xlabel('Expectation Value')
#    plt.ylabel('Probability')
#    plt.title(r'$\mathrm{Histogram\ of\ Position}$')
#    plt.axis([-max(path.X),max(path.X) , 0, 0.5])
#    plt.grid(True)
    plt.plot(rms)
    plt.show()

#def dump_data(final_path,rms,energy):
#    with open('PIMC_data.csv', 'w', newline='') as fp:
#        a = csv.writer(fp, delimiter=',')
#        data = [final_path, rms, energy]
#        a.writerows(data)


def PIMC(nsteps,path):

    #Thermalize
    print("Running Thermalization Steps...")
    for i in range(nsteps/2): #Run over a few steps for the first time to pass burn-in phase
        for j in range(path.M): # Run over the full time = dT*M
            mcstep(path)

    #Production Steps
    print("Running Production Steps...")
    rms = []
    energy = []
    for i in range(nsteps): # Run over the full range of steps to obtain energy values that can be worked with
        rms.append(np.sqrt(np.mean(np.square(path.X))))
        for j in range(path.M): #Run over the full time = dT*M
            x_new = mcstep(path)
            energy.append(path.E(x_new))
    render_plots(path.X,rms,energy)
    #dump_data(rms,energy)
