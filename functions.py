import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt


def rel_pos(N, dim, pos): #Week 1
    """
    Function that returns N ordered lists, each one computing the vectorial distance between particle 0,1,...N-1 and all of them
    For example, the third list coming out of this function will have N entries that will be the position of particle 2 and particle 0,
    particle 2 and particle 1, ..., particle 2 and particle N-1.
    """
    part_positions = np.zeros((N, N , dim)) # (N x N x dim) matrix including the relative positions
    for i in range(N):
        for j in range(N):
            part_positions[i,j,:] = pos[i,:] - pos[j,:]

    return part_positions


def bring_back(L, pos): #Week 2 (modified Week3)
    """
    Function the we apply AFTER getting the next step positions and velocities
    it checks if a particle has exited the box and if it has done, it replicates it in the other side of the box (actualizes
    positions according to periodic bounday consitions, note that we do not touch the velocity).
    
    ----------------------------------- FIRST PART OF THE PERIODIC BOUNDARY CONDITIONS--------------------------------------
    """
    pos = np.remainder(pos,L)
        
    return pos


def closest_rel_pos(N, dim, L, pos): #Week 2 (modified Week3)
    """
    Function that returns the closest distances of a particle to the N others taking into account the virual postition out
    of the domain L.
        
    ------------------------------------- SECOND PART OF THE PERIODIC BOUNDARY CONDITIONS----------------------------------------
    """
    r = rel_pos(N, dim, pos)
    r = np.remainder(r + L/2,L) - L/2

    return(r)


def closest_rel_dis(N, dim, L, pos): # Week 2
    """
    Function that returns the abosulte distance between particle i and particle j for each pair: 
    rel_dis(pos)[i][j]=abs(vec(r_i-r_j)), including
    the virtual ones when necessary. !!! Takes as argument the REAL positions !!!
    """
    r = closest_rel_pos(N, dim, L, pos)
    s = np.zeros((N, N))

    for i in range(dim):
        s[:,:] = s[:,:] + r[:,:,i]**2

    part_distances = np.sqrt(s)

    return part_distances


def potential(N, dim, L, pos): # Week 1 (Improved in Week 2)
    """
    Function that computes the potential energy between two pairs of particles and the total potential
    energy of the system, including the
    contributions of virtual particles. !!! Takes as argument the REAL positions !!!
    """
    d = closest_rel_dis(N, dim, L, pos)
    d_1 = d[d !=0]
    d_2 = np.reshape(d_1, (N, N-1))
    d_new = d_2**(-12) - d_2**(-6)

    P = 4 * np.sum(d_new, axis=1) #Potential on particle -a
    Ptot = 0.5 * np.sum(P) #Total Potential
    
    return P, Ptot


def force(N, dim, L, pos): # Week 4
    """
    Function that computes the force that each particle feels from the N-1 closest other ones
    """
    F = np.zeros((N,dim))
    r = closest_rel_pos(N, dim, L, pos) 
    
    r_2 = np.reshape(r, (N*N, dim))
    r_3 = np.zeros((N*N - N, dim))
    n = 0

    for i in range(N*N):
        if i % (N+1) != 0:
            r_3[n,:] += r_2[i,:]
            n += 1
    

    d = closest_rel_dis(N, dim, L, pos)
    d_1 = d[d !=0]
    d_2 = (2* d_1**(-12) - d_1**(-6)) * d_1**(-2)
    d_3 = np.reshape(d_2, (N*N - N, 1))
    d_4 = np.tile(d_3, (1, dim))

    F1 = np.reshape(24*r_3*d_4, (N, N-1, dim))
    F = np.sum(F1, axis = 1)

    return F


def vel_verlet(N, dim, L, pos, vel, h, loopnum): # Week 3
    """
    Function that utilizes the Velocity - Verlet algorithm to make a loopnum steps time evolution with timestep h.
    It also stores all the positions and velocities in all the different timesteps that have been computed
    (and plots the particles configuration in space every timestep in 2D if the commented block is uncommented).
    Takes the INITIAL positions & velocities as arguments !!!
    """
    lpn = loopnum
    pos_tracker = np.zeros((lpn + 1, N, dim)) # (lpn+1 x N x dim) matrix that tracks the positions of the particles for all the timeslots
    vel_tracker = np.zeros((lpn + 1, N, dim))
    pos_tracker[0,:,:] = pos
    vel_tracker[0,:,:] = vel
    
    for n in range(lpn):
        pos = pos_tracker[n,:,:] + vel_tracker[n,:,:]*h + 0.5 * h**2 * force(N, dim, L, pos_tracker[n,:,:])
        
        pos = bring_back(L, pos) # We need to apply the periodic BCs before plugging the new positions in the velocity computation.
        
        vel = vel_tracker[n,:,:] + 0.5 * h * (force(N, dim, L, pos_tracker[n,:,:]) + force(N, dim, L, pos))

        pos_tracker[n + 1,:,:] = pos
        vel_tracker[n + 1,:,:] = vel
        
        # uncomment this piece of code if you want to see the progress of the script
        """
        print(f"{(n/(lpn-1))*100} % verlet evolution completed")
        """
        
        # piece of code that plots the positions of the particles
        
        """"
        x, y = zip(*pos_tracker[n,:,:])
        plt.scatter(x, y)
        plt.ylim(0,L)
        plt.xlim(0,L)
        plt.title(f"Positions at time {n*h}")
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
        """
    return pos_tracker, vel_tracker


def k_p_energy(N, dim, L, initpos, initvel, h, loopnum): # Week 3 (Altered on week 4 to only return the kinetic energy of the system)
    """
    Fucntion that does a verlet evolution of time and stores the kinetic and potential energy of every timestep as well as the
    positions and the velocities. It takes as input the INITIAL positions and velocities.
    """
    lpn = loopnum
    pos, vel = vel_verlet(N, dim, L, initpos, initvel, h, loopnum) # positions and velocities of the particles
                                                            # for different timeslots, (lpn+1 x N x dim) matrix
    vel_1 = np.reshape(vel, ((lpn + 1) * N, dim))
    
    k = 0.5 * np.sum(vel_1**2, axis=1)
    kinetic_part = np.reshape(k, ((lpn + 1), N)) # A [(lpn+1) x N] matrix that stores the kinetic energy
                                                 # of each individual particle (columns) at each timeslot (lines)
    kinetic_system = np.sum(kinetic_part, axis=1) #The total kinetic energy of the system for each timeslot [(lpn+1) array]
    
    potential_part = np.zeros((lpn + 1, N))
    potential_system = np.zeros((lpn + 1,))

    for n in range(lpn + 1):
        potential_part[n,:], potential_system[n] = potential(N, dim, L, pos[n,:,:])
    
    return kinetic_system, potential_system, pos, vel


def total_energy(N, dim, L, initpos, initvel, h, loopnum): # Week 3 (Altered on week 4 to only return the total energy of the system)
    """
    Function that calls the k_p_energy fucntion and gives the same arguments as it but also the total energy.
    """
    kin_system, pot_system, postrack, veltrack= k_p_energy(N, dim, L, initpos, initvel, h, loopnum)
    total_system = kin_system + pot_system
    
    return total_system, kin_system, pot_system, postrack, veltrack


def fcc_2D(N, L, a): # Week 3 (milestione of Week 4)
    """
    This function returns the positions of a given number of particles N according to an FCC lattice (2 dimensions).
    L is the length of the whole square domain and a the length of the cell.
    vecs is a vector built out of vectors that state the position of every atom in a cell.
    """
    vecs = np.array([[0, 0], [0, a/2], [a/2, 0], [a/2, a/2]]) #4 vectors describing the fcc lattice (2D)
    
    dim = 2
    positions = np.zeros((N,dim))
    x0 = np.linspace(0,L-1,int(L/a))
    y0 = x0
    x,y = np.meshgrid(x0,y0)
    xm = np.reshape(x, int(N/len(vecs)),)
    x_tile=np.tile(xm, (len(vecs),1))
    xm = np.concatenate(x_tile)
    ym = np.reshape(y, int(N/len(vecs)),)
    y_tile = np.tile(ym, (len(vecs),1))
    ym = np.concatenate(y_tile)
    positions[:,0] = xm
    positions[:,1] = ym
    counts = 0
    for vectors in vecs:
        positions[counts*(int(N/len(vecs))):(counts+1)*(int(N/len(vecs)))]+=vectors
        counts += 1
    
    return positions


def fcc_3D(N, L, a): # Week 3 (milestione of Week 4)
    """
    This function returns the positions of a given number of particles N according to an FCC lattice (3 dimensions).
    L is the length of the whole square domain and a the length of the cell.
    vecs is a vector built out of vectors that state the position of every atom in a cell.
    """
    vecs = np.array([[0, 0, 0], [0, a/2, a/2], [a/2, 0, a/2], [a/2, a/2, 0]]) #4 vectors describing the fcc lattice
    
    dim = 3
    positions = np.zeros((N,dim))
    x0 = np.arange(0, L, a)
    y0 = x0
    z0 = x0
    x,y,z = np.meshgrid(x0,y0,z0)
    xm = np.reshape(x, int(N/len(vecs)),)
    x_tile=np.tile(xm, (len(vecs),1))
    xm = np.concatenate(x_tile)
    ym = np.reshape(y, int(N/len(vecs)),)
    y_tile=np.tile(ym, (len(vecs),1))
    ym = np.concatenate(y_tile)
    zm = np.reshape(z, int(N/len(vecs)),)
    z_tile=np.tile(zm, (len(vecs),1))
    zm = np.concatenate(z_tile)
    positions[:,0] = xm
    positions[:,1] = ym
    positions[:,2] = zm
    counts = 0
    for vectors in vecs:
        positions[counts*(int(N/len(vecs))):(counts+1)*(int(N/len(vecs)))] += vectors
        counts += 1
    
    return positions

def initialization_fcc(N, L, a, T, dim): # Week 4
    """ Function that initializes positions according to a 3D fcc lattice and velocities according to a Maxwell - Boltzmann distribution.
    The number of particles N must be given in relation to the length of our domain L and the size, a, of the unit cell by the formula:
    N = 4 * int(L/a)**dim.
    """
    if dim == 3:
        positions = fcc_3D(N, L, a)
    else:
        positions = fcc_2D(N, L, a)
    
    stdev = np.sqrt(T/119.8) # T in adimensional units, T_real = T*119.8K
    velocities = np.random.normal(0,stdev,(N,dim))
    velocities -= np.tile(avg_velocity(N, velocities), (N, 1)) 
    
    return positions, velocities

def avg_velocity(N, vel): #Week 4
    """
    Function that finds the average velocity of the particles (of N number) for each coordinate. Use it at the beginning
    with initial velocities.
    """
    avg = np.sum(vel, axis = 0) / N
    
    return avg

def energy_rescale(N, dim, L, T, initpos, initvel, h, loopnum, precision): #Week 4
    """
    Function that rescales the velocity of the particles such that it matches the desired temperature. Given the initial
    velocities and positions, it runs for loopnum/2 timesteps before computing the mean kinetic energy of the simulation from
    step loopnum/2 to step loopnum. Using this mean kinetic energy it calculates the first lamda value in dimensionless units by
    the formula lamda = sqrt((N-1)*T/(80<Ekin>)). If NOT |lamda - 1| < precision (1), it rescales the velocities as v --> λ*v
    and repeats the whole process until condition (1) is fullfiled.
    """
    n = 0 # Number of rounds
    lamda = 0
    
    pos = initpos
    vel = initvel
    
    kinetic_tracker = np.zeros((30, loopnum+1)) # 30 is arbitrary, the desired value for lamda is reached well before the 10th round
    potential_tracker = np.zeros((30,loopnum+1))
    total_tracker = np.zeros((30,loopnum+1))
    
    while np.abs(lamda - 1) >= precision:
        
        a, b, c, positions, velocities = total_energy(N, dim, L, pos, vel, h, loopnum)
        
        total_tracker[n] += a 
        kinetic_tracker[n] += b
        potential_tracker[n] += c
        
        mean_kin = np.mean(kinetic_tracker[n, int(loopnum/2):loopnum+1])
        mean_kin_squared = np.mean(kinetic_tracker[n, int(loopnum/2):loopnum+1]**2)
        lamda = np.sqrt((N-1)*T*3/(119.8*2*mean_kin))
        pos = positions[loopnum, :, :]
        vel = lamda * velocities[loopnum, :, :]
        
        # uncomment this piece of code if you want to see the progress of the script
        """
        print("rescalation loop nº: ", n+1)
        print("λ = ",lamda)
        """
        
        n += 1
    
    kTfin=mean_kin*2/(3*(N-1))
    a, b, c, positions, velocities = total_energy(N, dim, L, pos, vel, h, loopnum)
    total_tracker[n] += a 
    kinetic_tracker[n] += b
    potential_tracker[n] += c
    
    en_kin = np.reshape(kinetic_tracker[0:n+1], ((loopnum+1)*(n+1),))
    en_pot = np.reshape(potential_tracker[0:n+1], ((loopnum+1)*(n+1),))
    en_total = np.reshape(total_tracker[0:n+1], ((loopnum+1)*(n+1),))
    
    return en_kin, en_pot, en_total, n, pos, vel/lamda, kTfin


def pressuresum(N, dim, L, pos): # Week 4
    """
    Function that given the positions of the particles at a determined timestep, calculates the <...> value that the pressure
    expression depends on.
    """
    
    d = closest_rel_dis(N, dim, L, pos)
    d_1 = d[d != 0]
    d_2 = np.reshape(d_1, (N, N-1))
    d_new = 2*d_2**(-12) - d_2**(-6)

    A = -12 * np.sum(d_new, axis=1)
    Atot = np.sum(A)
    
    return Atot


def pair_correlation_close(N, dim, L, pos, bins): #Week 5
    """
    Function that returns an (histogram) array with the number of particles located with a certain distace d=[r, r+Dr] from
    another one. The distance d fulfills d \in [0, L/2] and the number of bins desired is an input.
    """
    dmax = L/2 # Maximal distance
    bin_size = dmax/bins
    Dr = np.arange(0, dmax, bin_size) # bin array, we choose to have 15 bins of width that depends on dmax
    Nr = np.zeros((bins,)) # The numbers of particles over given bins
    
    d = closest_rel_dis(N, dim, L, pos)# ??? Relative distances or CLOSEST relative distances (taking into acount also the virtual images)
    d_1 = d[d != 0]
    
    for i in range(bins-1):
        distances = d_1[(d_1 > Dr[i]) & (d_1 <= Dr[i+1])] # Finding the distance elements that lie in bin -i
        Nr[i] = len(distances) # Finding the NUMBER of distance elements that lie within this bin
        Nr[i] = Nr[i]/((bin_size*i+(bin_size/2))**2) # We divide by the position of the mid value of the bin (squared)
     
    distances = d_1[d_1 > Dr[bins-1]] # We did not count the elements laying in the last bin
    Nr[bins-1] = len(distances)
    
    Nr = 0.5 * Nr # We double-counted the symmetric elements
    
    return Nr


def func(x, tau): # Week 5
    """
    Exponential function with parameter tau.
    """
    return np.exp(-x/tau)


def autocorrelation(A, h): # Week 5
    """
    Function that takes as argument an array A of the measurements of a physical observable in all timesteps
    (what we call instantaneous measurements) and computes and plots its autocorrelation function.
    It also returns the mean value of all the measurements <A>,
    its correlation time τ and its standard deviation sigmaA. !!! Note that the array A must only include measurements
    of the physical quantity AFTER the rescaling has been already performed and equilibrium has been achieved.
    """
    lpn = len(A) - 1
    chi = np.zeros((lpn,))
    xdata = np.linspace(0,h*(lpn-1),lpn)
    
    for m in range(lpn): # We correlate all times -m (referring to the time of the autocorrelation function)
                         # times -n (timesteps in the array A)
        
        B = np.sum(A[: lpn+1 - m] * A[m:])
        C = A[: lpn+1 - m].sum() * A[m:].sum()
        D = (A[:lpn+1 - m]**2).sum()
        E = (A[:lpn+1 - m].sum())**2
        F = (A[m:]**2).sum()
        G = (A[m:].sum())**2
        
        chi[m] = ((lpn+1 - m)*B - C) / (np.sqrt((lpn+1 - m)*D - E) * np.sqrt((lpn+1 - m)*F - G))
        
    tau, pcov = opt.curve_fit(func, xdata, chi) # chi = ydata
        
    plt.plot(xdata[:100], chi[:100], 'b-', label='measured χ(t)') # We plot only the first 100 elements
                                                                  # assuming tau not greater than 100*h
    plt.plot(xdata[:100], func(xdata[:100], tau), 'r-', label=f'fitted χ(t), τ = {tau[0]}')
    plt.title('Autocorrelation function as a function of time')
    plt.xlabel('t (t*)')
    plt.ylabel('autocorrelation function χ')
    plt.legend()
    plt.show()
    
    Amean = np.mean(A)
    A2mean = np.mean(A**2)
    
    sigmaA = np.sqrt(2 * (tau/h) * (A2mean - Amean**2)/(lpn+1)) # Note that we divide tau by the timestep h because
                                                                # in our definition of tau it has units

    return Amean, tau, sigmaA[0]


def merge(N, dim, L, initpos, initvel, h, loopnum, bins):
    """
    Function that merges the calculations done in pressuresum2, pair_correlation_close2 and total_energy
    (if the functions are deleted please see the commits). It calls verlet 
    one time (for a given timestep and loopnum) and it returns:
    1. The array of the 'instantaneous pressures' computed for every timestep that the system has been through
    2. The array of 'instantaneous correlation functions (histograms) up to L/2 distance' (array of arrays)
    3. The array of 'instantaneous' total, potentail and kinetic energies
    4. The array of 'intantaneous' positions and velocities
    Must be called after rescaling.
    """
    lpn = loopnum
    pos, vel = vel_verlet(N, dim, L, initpos, initvel, h, lpn)
    pb_array = np.zeros((lpn + 1,))
    Nr_tracker = np.zeros((lpn + 1, bins))
    
    vel_1 = np.reshape(vel, ((lpn + 1) * N, dim))
    k = 0.5 * np.sum(vel_1**2, axis=1)
    kinetic_part = np.reshape(k, ((lpn + 1), N)) # A [(lpn+1) x N] matrix that stores the kinetic energy
                                                 # of each individual particle (columns) at each timeslot (lines)
    kinetic_system = np.sum(kinetic_part, axis=1) #The total kinetic energy of the system for each timeslot [(lpn+1) array]
    potential_part = np.zeros((lpn + 1, N))
    potential_system = np.zeros((lpn + 1,))

    
    for n in range(lpn + 1):
        pb_array[n] = pressuresum(N, dim, L, pos[n,:,:])
        Nr_tracker[n] = pair_correlation_close(N, dim, L, pos[n,:,:], bins)
        potential_part[n,:], potential_system[n] = potential(N, dim, L, pos[n,:,:])
        
        # uncomment this piece of code if you want to see the progress of the script
        """
        print(f"{(n/lpn)*100} % instantaneous ensemble values completed")
        """
    
    total_system = kinetic_system + potential_system
    
    return pb_array, Nr_tracker, total_system, kinetic_system, potential_system, pos, vel