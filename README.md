# SpatialMycelium
simulation code for the modelling of fungal growth

In mycelium_space.py, find the definitions of two classes : 
  - hypha : represents the segments of the filamentous structures
  - mycelium : mycelium represented by a list of hypha

With these classes we can simulate a mycelium growth up to a lifetime $T$ or to a number of branching event $n$.
You will find also two function acting on object of type mycelium to calculate the coordinate of each hypha in a simulated mycelium and to plot it.


In simulation_mycelium_space.py we run two simulations: 
  - the function "simulation" calculates the evolution of number of branching, crossing and apex during the simulation of a mycelium. It plots and fits the number of apex and branching to linear functions.
    
  - the functions "role_b1" and "role_b2" gives qualitative descriptions of the roles of b1 and b2 on the aspect of the mycelium. It calculates for differents values of the parameters the radius (maximal distance between the origin of the mycelium and its apex) and the number of internal nodes (crossings and branching events).
