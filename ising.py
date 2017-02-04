import numpy as np, matplotlib.pyplot as plt, matplotlib.animation as animation

def Ising_H(x,y):

    # Sum of products of spin configuration at (x,y) and nearest neighbour spin configurations
    s = L[x,y] * (L[(x+1) % l,y] + L[x, (y+1) % l] + L[(x-1) % l, y] + L[x,(y-1) % l])
    H = -J * s # Hamiltonian for Ising Model free of an external magnetic field
    return H

def mcstep(*args): #One Monte-Carlo Step - Metropolis Algorithm

    x = np.random.randint(l)
    y = np.random.randint(l)
    i = Ising_H(x,y)
    L[x,y] *= -1
    f = Ising_H(x,y)
    deltaH = f - i
    if(np.random.uniform(0,1) > np.exp(-deltaH/T)):
        L[x,y] *= -1

    mesh.set_array(L.ravel())
    return mesh,

def init_spin_config(opt):

    if opt == 'h':
        #Hot Start
        L = np.random.randint(2, size=(l, l)) #lxl Lattice with random spin configuration
        L[L==0] = -1
        return L

    elif opt =='c':
        #Cold Start
        L = np.full((l, l), 1, dtype=int) #lxl Lattice with all +1
        return L

if __name__=="__main__":

    l = 15 #Lattice dimension
    J = 0.3 #Interaction strength
    T = 2.0 #Temperature - change this from a constant to a variable that decreases over time and you have the Simulated Annealing algorithm
    N = 1000 #Number of iterations of MC step
    opt = 'h' #Option to choose between a hot('h') or cold('c') start

    L = init_spin_config(opt) #Initialize a spin configuration

    #Simulation Vizualization
    fig = plt.figure(figsize=(10, 10), dpi=80)
    fig.suptitle("T = %0.1f" % T, fontsize=50)
    X, Y = np.meshgrid(range(l), range(l))
    mesh = plt.pcolormesh(X, Y, L, cmap = plt.cm.RdBu)
    a = animation.FuncAnimation(fig, mcstep, frames = N, interval = 1, blit = True)
    plt.show()
