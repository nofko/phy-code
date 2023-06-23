import numpy as np
import matplotlib.pyplot as plt
import copy

N = 511
L = 6

tau = 0.3

x = np.linspace(-L/2,L/2,N+1)

xs = x*np.sqrt(6/(1-3*tau))

dxs = xs[1]-xs[0]

epsilon2 = 40 
A = -1.1

Au = 3*A/2

# XX = [u_0, ... , u_{N-1}]
# ZZ = full unknowns = [u_0, ... , u_{N-1}, c] 

# SPATIAL OPERATORS hence N points x_i = 1st diagonal block
one_m3 = np.roll(np.eye(N), -3, axis = 1)
one_m2 = np.roll(np.eye(N), -2, axis = 1)
one_m1 = np.roll(np.eye(N), -1, axis = 1)

one_p3 = np.roll(np.eye(N), +3, axis = 1)
one_p2 = np.roll(np.eye(N), +2, axis = 1)
one_p1 = np.roll(np.eye(N), +1, axis = 1)

M1 = -2/120*one_m3 + 18/120*one_m2 - 3/4*one_m1 + 3/4*one_p1 - 18/120*one_p2 + 2/120*one_p3
M3 = 3/24*one_m3 - one_m2 + 13/8*one_m1 - 13/8*one_p1 + one_p2 - 3/24*one_p3
M5 = -1/2*one_m3 + 2*one_m2 - 5/2*one_m1 + 5/2*one_p1 - 2*one_p2 + 1/2*one_p3

M1 /= dxs
M3 /= dxs**3
M5 /= dxs**5


def G(Z):
    """ spatial operator
    """
    GG = np.zeros_like(Z)
    y = Z[:-1]
    # spatial equations
    HH = -Z[-1]*np.matmul(M1,y) + y*np.matmul(M1,y) + np.matmul(M3,y) + epsilon2*np.matmul(M5,y)
    GG[:-1] = HH
    # central boundary condition
    GG[-1]  = Z[int((N+1)/2)]-Au

    return GG


def dG(Z):
    """ computes matrix of dg(Z_n) \cdot zeta
    """
    y = Z[:-1]
    # short_zeta = zeta[:-1]
    # spatial part
    N = len(Z)-1
    block1 = np.zeros((N,N))
    block1 = -Z[-1]*M1  + \
        np.diag(M1@y) + np.diag(y)@M1 + \
        M3 + epsilon2*M5
    #
    full_result = np.zeros((len(Z), len(Z)))
    full_result[:-1, :-1] = block1
    full_result[:-1,-1] = -M1@y
    full_result[-1,int((N+1)/2)] = 1
    #full_result[-1,:] = 1
    return full_result



y0 = Au*np.exp(-xs[:-1]**2*10)

Z0 = np.zeros(N+1)
Z0[:-1] = y0
Z0[-1] = 1

tol = 0.01
i = 0

err = np.linalg.norm(G(Z0))

Z = copy.deepcopy(Z0)
dZ = np.zeros_like(Z)

while err>tol and i<15:

    dZ = np.linalg.solve(dG(Z),G(Z))
    
    Z -= dZ
    
    err = np.linalg.norm(G(Z))
    
    i+=1

    print(err)
    

plt.close('all')
plt.figure(1)
plt.plot(y0,label="Guess")
plt.plot(Z[:-1],label="Solution")
plt.legend()
plt.show(block=False)
