System that forms a bicontinuous microemulsion

Contains an A-B diblock, f = 0.5, N = 20, phi = 0.11
In addition, equal amounts of A and B homopolymers with N = 5



Used py-gen-conf utility to initialize the following in the header section,
choosing the kappa*N = 200 and chi*N a little larger than 20

# Insert into py-gen-conf/main.py to generate a similar configuration
Dim = 2
poly_types = 3
rho0 = 3.2
phi = np.array([0.11, 0.89/2., 0.89/2.])

box = np.array([125.0, 125.0, 5.0])

N =     [[10,10],[5],[5]]   # Must be a list of lists
types = [[1 ,2], [1],[2]] # Also must be list of lists

name = "input.data"
