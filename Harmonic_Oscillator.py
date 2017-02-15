import numpy as np, matplotlib.pyplot as plt

def deltaH(x_p,k):

    k_p = k + 1
    if (k_p > M-1): k_p = 0
    k_m = k - 1
    if (k_m < 0): k_m = M-1
    dV = 0.5*(x_p**2) - 0.5*(x[k]**2)
    dK = ((0.5 * (x[k_p]-x_p))**2 + (0.5 * (x_p - x[k_m]))**2) - ((0.5 * (x[k_p] - x[k]))**2 + (0.5 * (x[k] - x[k_m]))**2)
    dH = dV + dK # Change in Energy for Harmonic Oscillator
    return dH

def mcstep(): #One Monte-Carlo Step - Metropolis Algorithm

    k = int(M * np.random.uniform()) #Pick an arbitrary index
    x_p = x[k] + (2 * np.random.uniform() - 1) * delta #Pick the position corresponding to arbitrary index
    dH = deltaH(x_p,k) #Calculate change in Energy
    if(dH < 0.0 or np.random.uniform(0,1) < np.exp(-dH*dT)): #Accept-Reject step
        x_new = x[k] = x_p
    else:
        x_new = x[k]

    return x_new #New position as determined by the acceptance probability

if __name__=="__main__":

    #Initialize
    M = 100 #Number of time slices
    T = 10 #Imaginary time period
    dT = T/M #Short interval into which the time period is divided
    delta = 1.0 #Metropolis step size
    MC_steps = 1000  #Number of Monte Carlo steps

    #Initialize a position configuration
    x = (np.random.uniform(size=(1, M)).tolist())[0]

    x_new  = 0.0
    E = 0.0

    #Thermalize
    for i in range(MC_steps/2): #Run over a few steps for the first time to pass burn-in phase
        for j in range(M): # Run over the full time = dT*M
            mcstep()

    #Production Steps
    for i in range(MC_steps): # Run over the full range of steps to obtain energy values that can be worked with
        for j in range(M): #Run over the full time = dT*M
            x_new = mcstep()
            e = 0.5*(x_new**2) + (0.5 * x_new * x_new) # Energy from Virial Theorem - V+(0.5*x*dV/dx)
            E += e #Accumulate energy

    E_avg = E/(MC_steps*M) #Average Energy
    print(E_avg)

    #Plot
    #fig = plt.figure(figsize=(10, 10), dpi=80)
    #plt.xlim(x_min,x_max)
    #plt.ylim(0,1)
    #plt.plot(P)
    #plt.show()
