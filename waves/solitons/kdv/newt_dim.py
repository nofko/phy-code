import numpy as np
import matplotlib.pyplot as plt
import copy

N = 1024

L = 0.9

x = np.linspace(-L/2,L/2,N+1)

r = 7
ro = 8.3
w = (ro-r)*2
chi = ro/r
wt = w/2

phi3 = -1
phi5 = 16

bond = 0.85**2/wt**2/chi**4

delta = -phi3/6-bond
gamma = (-phi5/120-bond*phi3/6+delta**2/4)

a_scale = 5/wt/4
x_scale = np.sqrt(2/abs(delta))*ro/wt/chi

A = -0.5
Au = A*a_scale

epsilon2 = -2*chi**2*gamma/delta**2

xs = x*x_scale
dxs = xs[1]-xs[0]

# XX = [u_0, ... , u_{N-1}]
# ZZ = full unknowns = [u_0, ... , u_{N-1}, c] 

# SPATIAL OPERATORS hence N points x_i = 1st diagonal block
one_m3 = np.diag(np.ones(N-3),-3)
one_m2 = np.diag(np.ones(N-2),-2)
one_m1 = np.diag(np.ones(N-1),-1)

one_p3 = np.diag(np.ones(N-3),3)
one_p2 = np.diag(np.ones(N-2),2)
one_p1 = np.diag(np.ones(N-1),1)


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
    HH = -Z[-1]*np.matmul(M1,y) + y*np.matmul(M1,y) - np.matmul(M3,y) + epsilon2*np.matmul(M5,y)
    GG[:-1] = HH
    # central boundary condition
    GG[-1]  = Z[int((N)/2)]-Au

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
        -M3 + epsilon2*M5
    #
    full_result = np.zeros((len(Z), len(Z)))
    full_result[:-1, :-1] = block1
    full_result[:-1,-1] = -M1@y
    full_result[-1,int((N)/2)] = 1
    #full_result[-1,:] = 1
    
    return full_result



y0 = Au*np.exp(-xs[:-1]**2*10)

Z0 = np.zeros(N+1)
Z0[:-1] = y0
Z0[-1] = 1e-5

tol = 0.0001
i = 0

err = np.linalg.norm(G(Z0))

Z = copy.deepcopy(Z0)
dZ = np.zeros_like(Z)

Z[0] = 0
Z[-2] = 0

while err>tol and i<15:

    Z[0] = 0
    Z[-2] = 0

    dZ = np.linalg.solve(dG(Z),G(Z))
    
    Z -= dZ
    
    err = np.linalg.norm(G(Z))
    
    i+=1

    
print("Error: ",err)
print("epsilon^2: ",epsilon2)
#plt.close('all')
plt.figure(1)
plt.plot(y0,label="Guess")
plt.plot(Z[:-1],label="Solution")
#plt.legend()
plt.show(block=False)


plt.figure(2)
#plt.plot(y0,label="Guess")
plt.plot(xs[:-1]/x_scale,Z[:-1]/a_scale,label="Solution")
#plt.legend()
plt.show(block=False)
