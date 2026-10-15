"""
This files contains:
  - a simulation to study the evolution of number of apex, crossing and branching during the simulation
  - two simulation to study qualitatively the effect of b1 and b2 on the aspect of the mycelium
"""

import numpy as np
from time import time
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

from mycelium_space import hypha, mycelium, simulate_mycelium, simulate_mycelium_to_n
from mycelium_space import plot_mycelium, coordinate


#%% Data
theta_e = -10*np.pi/180
theta_0 = 71*np.pi/180
theta_dep = [np.pi/6, 3*np.pi/4, 3*np.pi/2]
theta_p = theta_0



#%% Simulation 

def simulate_mycelium_step(b1, b2, Tmax,
                           n_points=10, print_time=True):
    """
    One step of simulation to calculates the evolution of the number 
    of branching, of crossing and of apex during the expansion of the mycelium

    Parameters
    ----------
    b1 : float
        apical branching rate
    b2: float
        lateral branching rate
    Tmax : float
        Maximum lifetime of the mycelium
    n_points : int, optional
        number of discretization points. The default is 10.
    print_time : bool, optional
        if True prints the simulation time. The default is True.

    Returns
    -------
    M : mycelium
        mycelium simulated.

    """
    if print_time:
        print("Start of simulation")
        t = time()
    
    M = np.zeros((3, n_points+1)) #matrix to stock number of branching, of crossing and of apex 
    M[:, 0] = [1, 0, 3]
    m = mycelium(b1, b2)
    
    
    for i in range(1, n_points+1):
        m.run(Tmax=i*Tmax/n_points)
        S = coordinate(m)
        M[:, i] = [m.n_apex() - 2, n_cross(S), m.n_apex()]
    
    if print_time :
        t = time() - t
        print(f'Time of simulation: {t} s')
        
    return M


def simulation(b1, b2, Tmax, n_points, n_myc, to_plot = True, to_print = True): 
    """
    Calculates the evolution of number of branching, crossing and apex during 
    the simulation of a mycelium. Plots and fits the number of apex et branching.

    Parameters
    ----------
    b1 : float
        apical branching rate
    b2: float
        lateral branching rate
    Tmax : float
        Maximum lifetime of the mycelium
    n_points : int, optional
        number of discretization points. The default is 10.
    n_myc : int
        ,number of mycelium to simulate.
    to_plot : bool, optional
        DESCRIPTION. The default is True.
    to_print : bool, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    dict
        DESCRIPTION.

    """
    t = time()
    M = np.zeros((3, n_points+1))
    if to_print:
        print("Start of the simulation")
    
    for i in range(n_myc):
        M += simulate_mycelium_step(b1, b2, Tmax, n_points, print_time = False)
    
    N_branching = M[0, :]/n_myc
    N_cross = M[1, :] /n_myc
    N_apex = M[2, :] /n_myc
    N_int_nodes = N_branching + N_cross
    
    T = np.linspace(0, Tmax, n_points + 1)
    
    t = time() - t
    if to_print:
        print(f"Simulation time: {t} s")
    
    #linear fit of A
    alpha_A, A0 = np.polyfit(T, np.log(N_apex), deg = 1)
    if to_plot : 
        plt.figure()
        plt.plot(T, np.log(N_apex), label = r"$log(A)$")
        plt.plot(T, alpha_A*T+A0, label =r"fit $\alpha_A T + A_0$")
        plt.title(rf"Logaritmic fit of $A$ for {n_myc} mycelium simulated ($b_1 =$ {b1}, $b_2 =$ {b2}, T_max = {Tmax})")
        plt.legend()
        plt.show()
    if to_print :
        print(f"Linear Fit of A : alpha_A = {alpha_A} and A_0 = {A0}")

    #linear fit of N_int_nodes
    alpha_N, N0 = np.polyfit(T, np.log(N_int_nodes), deg = 1)
    if to_plot : 
        plt.figure()
        plt.plot(T, np.log(N_int_nodes), label = r"$log(N)$")
        plt.plot(T, alpha_N*T+N0, label =r"fit $\alpha_N T + N_0$")
        plt.title(rf"Logaritmic fit of $N$ for {n_myc} mycelium simulated ($b_1 =$ {b1}, $b_2 =$ {b2}, T_max = {Tmax})")
        plt.legend()
        plt.show()
    if to_print:
        print(f"Linear Fit of N : alpha_N = {alpha_N} and N0 = {N0}")


    alpha_b = (alpha_A + alpha_N)/2
    alpha_a = (alpha_N - alpha_A)/2
    lam = (b1 + (b1**2+4*b2)**1/2)/2
    
    if to_print:
        print(f"alpha_a = {alpha_a} and alpha_b = {alpha_b}")
    
    return {"N_int_nodes" : N_int_nodes, "N_cross" : N_cross, 
            "N_branching" : N_branching, "N_apex" : N_apex, 
            "lin fit A" : [alpha_A, A0], "lin fit N" : [alpha_N, N0], 
            "alpha_a" : alpha_a, "alpha_b" : alpha_b, "lambda" : lam}


#%% Function to calculate the mycelium's summary statistics


#Radius of the mycelium

def rad(S):
    """
    Calculates the maximal distance between the origin of the mycelium and the 
    apex. 

    Parameters
    ----------
    S : TYPE
        DESCRIPTION.

    """
    dist = 0
    for i in range(len(S)):
        dist = max(dist, S[i][2]**2 + S[i][3]**2)
    return dist

def rad(S):
    return max(S[:][2]**2 + S[:][3]**2)



#Number of crossing of the hypha

def ccw(A, B, C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect


def test(A, B, C, D):
    return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)


def intersect(H1, H2):
    A = [H1[0], H1[1]]
    B = [H1[2], H1[3]]
    C = [H2[0], H2[1]]
    D = [H2[2], H2[3]]
    return test(A, B, C, D) and test(B, A, C, D)


def n_cross(S, print_time=False):
    """
    Counts the number of crossings in the mycelium

    Parameters
    ----------
    S (List): Coordinates of each hyphe [xd, yd, xa, ya].

    Returns
    -------
    n_cross
    """
    t = time()
    n_cross = 0
    l = len(S)
    for i in range(l):
        for j in range(i+1, l):
            n_cross += intersect(S[i], S[j])
    if print_time:
        print(time() - t)
    return n_cross


#%% Qualitative description of the role of b1

n_max = 500
B1 = [0.01, 0.04, 0.07, 0.1, 0.4]
B1 = np.linspace(0.01, 0.5, 10)
n_myc = 5

def role_b1(n_max = 500, B1, n_myc = 5):
  """
    Plots the radius and the number of internal nodes for differents values of b1

    Parameters
    ----------
    n_max (int): maximal number of branching to which we simulate the myc
    B1 (List): list of value of b1
    n_myc (int): Number of mycelium simulated. The default is 5.

    Returns
    -------
    None.

    """
    M = np.zeros((3, len(B1)))
    for i in range(len(B1)):
      t = time()
      print(f"start {i}th it")
      for j in range(n_myc):
          b1 = B1[i]
          m = simulate_mycelium_to_n(b1, 0.006, n_max, print_time=False)
          S = coordinate(m)
          M[:, i] += [rad(S), n_cross(S), m.n_apex() - 2]
      M[:, i] /= n_myc
      print(f"end: {time() - t} s")
    
    N_int_nodes = M[1, :] + M[2, :]
    
    plt.figure()
    plt.plot(B1, M[0, :], label = 'Rayon', linestyle = '--', marker = 'o')
    plt.title("Role of b1 in the radius of the surface covered by the mycelium")
    plt.xlabel('b1')
    plt.ylabel('R(b1)')
    plt.legend()
    plt.show()
    
    
    plt.figure()
    plt.plot(B1, N_int_nodes, label = 'N', linestyle = '--', marker = 'o')
    plt.title("Role of b1 in the number of internal nodes of the mycelium")
    plt.xlabel('b1')
    plt.ylabel('R(b1)')
    plt.legend()
    plt.show()


role_b1(n_max, B1, n_myc)


#%% Qualitative description of the role of b2


n_max = 500
B2 = np.linspace(0.001, 0.03, 10)
n_myc = 5


def role_b2(n_max = n_max, B2 = B2, n_myc = 5):
    """
      Plots the radius and the number of internal nods for differents values of b2

      Parameters
      ----------
      n_max (int): maximal number of branching to which we simulate the myc
      B2 (List): list of value for b2
      n_myc (int): Number of mycelium simulated. The default is 5.

      Returns
      -------
      None.
      
      """
    
     M = np.zeros((3, len(B2)))
     for i in range(len(B2)):
        t = time()
        print(f"start {i}th it")
        for j in range(n_myc):
            b2 = B2[i]
            m = simulate_mycelium_to_n(0.27, b2, n_max, print_time=False)
            S = coordinate(m)
            M[:, i] += [rad(S), n_cross(S), m.n_apex() - 1]
            
        M[:, i] /= n_myc
        print(f"end: {time() - t} s") 
           
           
    N_int_nodes = M[1, :] + M[2, :]
        
    plt.figure()
    plt.plot(B2, M[0, :], label = 'Rayon', linestyle = '--', marker = 'o')
    plt.title("Role of b2 in the radius of the surface covered by the mycelium")
    plt.xlabel('b2')
    plt.ylabel('R(b2)')
    plt.legend()
    plt.show()
        
        
    plt.figure()
    plt.plot(B2,N_int_nodes, label = 'N', linestyle = '--', marker = 'o')
    plt.title("Role of b2 in the number of internal nodes of the mycelium")
    plt.xlabel('b2')
    plt.ylabel('R(b2)')
    plt.legend()
    plt.show()
    
role_b2(n_max, B2, n_myc)
