"""
This file contains the definitions of the hypha and mycelium classes.
These classes are used to encode the mycelium in the simulation of the process.
"""

import numpy as np
from time import time
import matplotlib.pyplot as plt

# %% Data
"""Angles used to represent the mycelium spatially"""

#angles for the apical branching
theta_e = -10*np.pi/180
theta_0 = 71*np.pi/180

#angle for the lateral branching
theta_p = theta_0
#theta_p = np.pi/2

#angles for the first branching
theta_s = [np.pi/6, 3*np.pi/4, 3*np.pi/2]
#theta_s = [2*i*np.pi/3 for i in range(3)]

# %%
#=============================================================================#


class hypha:
    """
    Definition of the hypha class. 

    Attributes
    ----------
    _len: float
          length of the hypha
    _e: bool
        True if open, False if closed
    _s : bool
        True if alive, False if closed
    _i : int
        index of the mother hypha
    _typeb : bool
            True if coming from an apical branching, False if lateral branching
    """

    def __init__(self, x=0, e=True, s=True, i=0, typeb=False):
        self._len = x
        self._e = e
        self._s = s
        self._i = i
        self._typeb = typeb

    def elongation(self, T=1):
        """
        Updates the length of the hypha after a duration T
        """
        if self._e:
            if self._s:
                self._len += T  # open hyphae lengthen at speed 1

    def length(self):
        """
        Returns the length of the hypha
        """
        return self._len

    #=========================================================================#
    # Update the hypha after branching and create the two new hyphae that result from the branching event
    def apical_branching(self, i):
        """
        Updates the hypha after apical branching
        i : index of the current hypha 
        """
        self._e = False  # Close the hypha
        new1 = hypha(x=0, e=True, s=True, i=i, typeb=True)
        # Create two open hyphae of length 0
        new2 = hypha(x=0, e=True, s=True, i=i, typeb=True)
        return(new1, new2)

    def lateral_branching(self, i):
        """
        Updates the hypha after lateral branching
        i : index of the current hypha
        """
        if self.length() != 0:
            # Draw the branching point (with a margin around the branching points)
            z = np.random.uniform(1e-5, self.length()-1e-5) 

            # create a closed hypha, the root of the branch
            new1 = hypha(x=z, e=False, s=True, i=i, typeb=False)

            # Create a hypha of the same type with length x-i
            new2 = hypha(x=self.length() - z, e=self._e,
                         s=True, i=i, typeb=False)

            # Create an open hypha of length 0
            new3 = hypha(x=0, e=True, s=True, i=i, typeb=False)

            self._s = False  # kill the current hypha
            self._e = False

            return(new1, new2, new3)

    def describe(self):
        """
        Prints a description of the hypha with each of its attributes
        """
        print("len: ", self._len)
        print("e: ", self._e)
        print("state: ", self._s)
        print("mother indx: ", self._i)
        print("type of branching: ", self._typeb)

#=============================================================================#
# %%


class mycelium:
    """
    Definition of the mycelium class

    Attributes
    ----------
    _hyphae: list[hypha]
             list of hyphae
    _b1: float
         apical branching rate
    _b2: float
         lateral branching rate
    _T: float
        lifetime
    """

    def __init__(self, b1, b2):
        self._hyphae = []

        for i in range(3):
            self._hyphae.append(hypha(x=np.random.uniform(0, 0.1)))

        self._b1 = b1
        self._b2 = b2

        self._T = 0

    def run(self, Tmax, non_unif=False):
        """
        Updates the mycelium up to time Tmax
        """

        t = time()
        while self._T < Tmax:
            self.actualisation(non_unif)

        t = time() - t

        return t

    def run_to_n(self, n_max, non_unif=False):
        """
        Updates the mycelium until the number of hyphae is greater than n_max
        """
        from time import time

        t = time()
        while self.n_apex() <= n_max:
            self.actualisation(non_unif)
        t = time() - t

        return t

    def __len__(self):
        """
        Returns the number of hyphae
        """
        return len(self._hyphae)

    def __activelen__(self):
        """
        Returns the number of active (alive) hyphae
        """
        return sum([h._s for h in self._hyphae])

    def __getitem__(self, i):
        """
        Returns the ith hypha
        """
        return self._hyphae[i]
    #=========================================================================#

    def active(self):
        """
        Returns the list of hyphae alive
        """
        return [h for h in self._hyphae if h._s]

    def n_apex(self):
        """
        Returns the number of apex
        """
        return sum([h._e for h in self.active()])

    def total_length_open(self):
        """
        Returns the total length of the open hyphae
        """
        return sum([h.length() for h in self.active() if h._e])

    def total_length(self):
        """
        Returns the total length of the mycelium
        """
        return sum([h.length() for h in self.active()])

    def length_open(self):
        """
        Returns a list of the lengths of each open hypha
        """
        return [h.length() for h in self.active() if h._e]

    def length_closed(self):
        """
        Returns a list of the lengths of each closed hypha
        """
        return [h.length() for h in self.active() if not h._e]

    def length_hyphae(self):
        """
        Returns a list of the lengths of each hypha
        """
        return [h.length() for h in self.active()]
    #=========================================================================#

    def m_apical_branching(self, non_unif=False):
        """
        Draws the waiting time to the next lateral branching event 
        and the index of the hypha that branches
        """
        n = self.n_apex()

        self.__T1 = np.random.exponential(1/(n*self._b1))  # Draw of the waiting time

        # Indices of open hyphae
        id_apex = [i for i in range(len(self)) if self[i]._e and self[i]._s]

        if non_unif:
            p = [h._len/self.total_length_open() for h in self.active() if h._e]

            # non-uniform draw of the branching apex
            self.__id1 = np.random.choice(id_apex, p=p)

        else:
            # uniform draw of the branching apex
            self.__id1 = np.random.choice(id_apex)

    def m_lateral_branching(self):
        """
        Draws the waiting time to the next lateral branching event 
        and the index of the hypha that branches
        """
        def F_inv(u):  # Inverse of the cdf of the waiting time to the next lateral branching event
            L = self.total_length()
            A = self.n_apex()

            delta = (self._b2*L)**2 - 2*self._b2*A * np.log(1-u)

            return (-self._b2*L + np.sqrt(delta)) / (self._b2*A)

        # Draw of the waiting time using inverse cdf method
        self.__T2 = F_inv(np.random.uniform())

        p = [self[i]._s*(self[i].length() + self[i]._e*self.__T2)
             for i in range(len(self))]  # List of updated lengths of hyphae
        p = np.array(p) / sum(p)  # List of probabilities of branching
        if any(p<0):
            print(p)
            
        a = list(range(self.__len__()))  # List of indices of all hyphae

        # Draw of the index of the branching hypha according to p
        self.__id2 = np.random.choice(a, p=p)

    def actualisation(self, non_unif=False):
        """
        Updates the mycelium up to the next branching event
        --------
        non_unif : bool
                True if draw of the branching apex is not uniform
        """
        self.m_apical_branching(non_unif)
        self.m_lateral_branching()

        T = min(self.__T1, self.__T2)  # Waiting time to next branching event
        for i in range(len(self)):
            self[i].elongation(T)  # Update of the lengths of the hyphae

        self._T += T  # Update of the lifetime of the mycelium

        if self.__T1 < self.__T2:  # Apical branching occurs first
            new1, new2 = self[self.__id1].apical_branching(
                i=self.__id1)  # Hyphae created by apical branching
            self._hyphae.append(new1)
            self._hyphae.append(new2)
        else:  # Lateral branching occurs first
            new1, new2, new3 = self[self.__id2].lateral_branching(
                i=self.__id2)  # Hyphae created by lateral branching
            self._hyphae.append(new1)
            self._hyphae.append(new2)
            self._hyphae.append(new3)



#%% Functions to simulate a mycelium

def simulate_mycelium(b1, b2, T_max, print_time=True, non_unif=False):
    """
    Simulates the mycelium up to time T_max

    Args: 
        b1 (float): apical branching rate
        b2 (float): lateral branching rate
        n_run (int): number of simulations
        T_max (float): lifetime at the end of the simulation
        print_time (bool): print time of simulation

    Returns:
        m : simulated mycelium
    """
    if print_time:
        print('start of simulation')

    m = mycelium(b1, b2)
    t = m.run(T_max, non_unif)

    if print_time:
        print("time of sim = " + str(t) + 's')

    return m


def simulate_mycelium_to_n(b1, b2, n_max, print_time=True, non_unif=False):
    """Simulates the mycelium until the number of hyphae is greater than n_max

    Args: 
        b1 (float): apical branching rate
        b2 (float): lateral branching rate
        n_run (int): number of simulations
        n_max (float): number of hyphae at the end of the simulation 
        print_time (bool): print time of simulation

    Returns:
        m : simulated mycelium
    """

    if print_time:
        print('start of simulation')

    m = mycelium(b1, b2)
    t = m.run_to_n(n_max, non_unif)

    if print_time:
        print("time of sim = " + str(t) + 's')

    return m


#%% Preliminary function to plot the mycelium

def coordinate(m, thetas=theta_s, theta0=theta_0,
               thetae=theta_e, thetap=theta_p, print_time=False, ret_fil = False):
    """
    Calculates the coordinates of the hyphae

    Parameters
    ----------
    m : mycelium
        simulated mycelium to plot
    thetas : list, float
        list of the three starting angles. The default is theta_s.
    theta0 : float 
        Angle of the operating branch in apical branching. The default is theta_0.
    thetae : float
        Angle of the exploratory branch in apical branching. The default is theta_e.
    thetap : float
        Angle of lateral branching. The default is theta_p.
    ret_fil : bool
        if True, the function returns the filiation of each hypha G. The default is False.


    Returns
    -------
    S : list
        lists of the coordinates of each hyphe of the mycelium network
    G : list
        indexes of the hyphae in filiation with each hyphe

    """
    t = time()
    L = []  # list to stock ((xa,ya), theta): starting point and starting angle
    indx = []  # list of indexes to link L and S
    S = []  # list to stock ((xs, ys), (xa, ya)) of the selected hyphae
    n = m.__len__()
    if ret_fil : 
        G = []

    # calcul of the hyphae coordinate
    for i in range(3):
        h = m._hyphae[i]
        x = h._len * np.cos(thetas[i])
        y = h._len * np.sin(thetas[i])
        L.append([x, y, thetas[i]])
        S.append([0, 0, x, y])
        indx.append(i)
        
        if ret_fil : 
            G.append([])

    i = 3
    while i in range(3, n):
        h1 = m._hyphae[i]
        Hm = L[h1._i]
        s = 1
        l = len(S)

        if h1._typeb:  # apical branching
        
            x_s, y_s, angle_s = Hm

            x1 = x_s + h1._len * np.cos(s*theta0 + angle_s)
            y1 = y_s + h1._len * np.sin(s*theta0 + angle_s)
            L.append([x1, y1, s*theta0 + angle_s])
            S.append([x_s, y_s, x1, y1])
            indx.append(l)
            if ret_fil : 
                G.append([])

            h2 = m._hyphae[i+1]
            x2 = x_s + h2._len * np.cos(s*thetae + angle_s)
            y2 = y_s + h2._len * np.sin(s*thetae + angle_s)
            L.append([x2, y2, s*thetae + angle_s])
            S.append([x_s, y_s, x2, y2])
            indx.append(l+1)
            if ret_fil: 
                G.append([])
                G[indx[h1._i]] = [l, l+1]
                
            
            i += 2

        else:  # lateral branching
        
            x_s, y_s, angle_s = Hm
            hm = m._hyphae[h1._i]
            h2 = m._hyphae[i+1]
            h3 = m._hyphae[i+2]

            x_s = x_s - hm._len * np.cos(angle_s)
            y_s = y_s - hm._len * np.sin(angle_s)

            x1 = x_s + h1._len * np.cos(angle_s)
            y1 = y_s + h1._len * np.sin(angle_s)
            L.append([x1, y1, angle_s])

            S[indx[h1._i]] = [x_s, y_s, x1, y1]
            indx.append(indx[h1._i])

            x2 = x1 + h2._len * np.cos(angle_s)
            y2 = y1 + h2._len * np.sin(angle_s)
            L.append([x2, y2, angle_s])
            S.append([x1, y1, x2, y2])
            indx.append(l)

            x3 = x1 + h3._len * np.cos(angle_s + s*thetap)

            y3 = y1 + h3._len * np.sin(angle_s + s*thetap)
            L.append([x3, y3, angle_s + s*thetap])
            S.append([x1, y1, x3, y3])
            indx.append(l+1)

            if ret_fil : 
                G.append(G[indx[h1._i]])
                G.append([])
                G[indx[h1._i]] = [l, l+1]
            i += 3
    if print_time:
        print(time() - t)
        
    if ret_fil : 
        return S, G
    
    return S

# %% Function to plot the mycelium in space

def plot_mycelium(m, thetas=theta_s, theta0=theta_0,
                  thetae=theta_e, thetap=theta_p, equal=True):
    """
    Plots the mycelium in space

    Args : 
        m : simulated mycelium
        thetas (list) : list of the three angles to start with
        theta0, thetae (float): angles of the apical branching
        thetap (float): angle of the lateral branching
        equal (bool) : True if the plots has to be orthonormal

    """
    L = []  # list to stock ((x,y), theta)
    n = m.__len__()

    plt.figure()
    for i in range(3):
        h = m._hyphae[i]
        x = h._len * np.cos(thetas[i])
        y = h._len * np.sin(thetas[i])
        L.append([x, y, thetas[i]])
        plt.plot([0, x], [0, y], color="black")
    i = 3
    while i in range(3, n):
        h1 = m._hyphae[i]
        Hm = L[h1._i]
        #s = np.random.choice([-1, 1], p=[1/2, 1/2])
        s = 1

        if h1._typeb:  # apical branching
            x_s, y_s, angle_s = Hm

            x1 = x_s + h1._len * np.cos(s*theta0 + angle_s)
            y1 = y_s + h1._len * np.sin(s*theta0 + angle_s)
            L.append([x1, y1, s*theta0 + angle_s])
            plt.plot([x_s, x1], [y_s, y1])

            h2 = m._hyphae[i+1]
            x2 = x_s + h2._len * np.cos(s*thetae + angle_s)
            y2 = y_s + h2._len * np.sin(s*thetae + angle_s)
            L.append([x2, y2, s*thetae + angle_s])
            plt.plot([x_s, x2], [y_s, y2])

            i += 2

        else:  # lateral branching
            x_s, y_s, angle_s = Hm
            hm = m._hyphae[h1._i]
            h2 = m._hyphae[i+1]
            h3 = m._hyphae[i+2]

            x_s = x_s - hm._len * np.cos(angle_s)
            y_s = y_s - hm._len * np.sin(angle_s)

            x1 = x_s + h1._len * np.cos(angle_s)
            y1 = y_s + h1._len * np.sin(angle_s)
            L.append([x1, y1, angle_s])

            x2 = x1 + h2._len * np.cos(angle_s)
            y2 = y1 + h2._len * np.sin(angle_s)
            L.append([x2, y2, angle_s])
            plt.plot([x1, x2], [y1, y2])

            x3 = x1 + h3._len * np.cos(angle_s + s*thetap)

            y3 = y1 + h3._len * np.sin(angle_s + s*thetap)
            L.append([x3, y3, angle_s + s*thetap])
            plt.plot([x1, x3], [y1, y3])

            i += 3

        plt.title(
            f"Simulation of a Mycelium (up to time T = {m._T}, b1 = {m._b1}, b2 = {m._b2})")
        if equal:
            plt.axis('equal')
        plt.show()


