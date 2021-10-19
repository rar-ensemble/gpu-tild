System that should undergo spinodal decomposition

Contains a random 50/50 mixture of N = 10 polymers



Used py-gen-conf utility to initialize the following in the header section,
choosing the kappa*N = 200 and chi*N = 10


# Insert into py-gen-conf/main.py to generate a similar configuration

Dim = 2
poly_types = 2
rho0 = 3.2
phi = np.array([0.5, 0.5])

box = np.array([125.0, 125.0, 5.0])

N =     [[10],[10]]   # Must be a list of lists
types = [[1],[2]] # Also must be list of lists

name = "input.data"
