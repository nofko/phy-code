import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import fsolve

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)

## ============================================================================
## THEORETICAL ANALYSIS                                                      ##
##                                                                           ##
##             CALCULATION OF DISPERSION OVER 1D PERIODIC BOTTOM             ##
##                                                                           ##
## ============================================================================


### PHYSICS ###

g = 9.81   # M/S^2
sigma = 0.0
rho = 1000

L = 1    # M
om0 = np.sqrt(g/L)

d = 0.2*L

h = L   # M
h1 = 0.05*L

### DISPERSION RELATION ###


def disp_gravity(KK,HH):

    return np.sqrt((g*KK+sigma*KK**3/rho)*np.tanh(KK*HH))


####### FUNCTION TO SOLVE ########

def zero_gravity(KK,HH,OM):

    return disp_gravity(KK,HH)-OM


####### SCALAR SOLUTION #######

def inv_gravity(OM,HH):
    
    sol = fsolve(zero_gravity,0.2,args=(HH,OM))
        
    return sol


####### TRANSMISSION MATRIX #######


def T(k,ki,d,L,x):

    g = k/ki

    t0 = np.array([[(1+g)*np.exp(1j*(ki-k)*x),(1-g)*np.exp(-1j*(ki+k)*x)],
                   [(1-g)*np.exp(1j*(ki+k)*x),(1+g)*np.exp(-1j*(ki-k)*x)]])

    t1 = np.array([[(1+1/g)*np.exp(1j*(k-ki)*(x+d)),(1-1/g)*np.exp(-1j*(ki+k)*(x+d))],
                   [(1-1/g)*np.exp(1j*(ki+k)*(x+d)),(1+1/g)*np.exp(-1j*(k-ki)*(x+d))]])

    return 0.25*t0@t1

### SOLVING ###

N = 100    ## NUMBER OF STEPS

Nf = 100   ## FREQUNCY GRID

M = 100    ## NUMBER OF ENSEMBLES

R = 0.5    ## RANDOMNESS

om = np.linspace(0.1,6,Nf)*om0

def fullT(N,r):

    TTtot = np.zeros((Nf,2,2),dtype=complex)
    
    delta = r*(np.random.rand(N)-0.5)*(L-d)

    x = np.arange(N)*L+delta

    for i in range(Nf):
    
        k = inv_gravity(om[i],h)[0]
        k1 = inv_gravity(om[i],h1)[0]

        TT = np.identity(2)
    
        for j in range(N):

            tt = T(k,k1,d,L,x[j])
        
            TT = TT@tt

        TTtot[i] = TT

    A = 1/np.abs(TTtot[:,0,0])
    gamma = -np.log(A**2)

    return A**2,gamma/N


ST = np.zeros(Nf)
GG = np.zeros(Nf)

for i in range(M):

    a, gamma = fullT(N,R)

    ST += a
    GG += gamma

ST /= M 
GG /= M

A = np.log(ST)

fig, ax = plt.subplots(figsize=[12,8])

ax.plot(om/om0,A)

ax.set_xlabel(r"$\omega/\omega_0$ ",fontsize=35)
ax.set_ylabel("$\ln(T)$",fontsize=35)

ax.tick_params(labelsize=25)


fig, ax = plt.subplots(figsize=[12,8])

ax.plot(om/om0,GG)

ax.set_xlabel(r"$\omega/\omega_0$ ",fontsize=35)
ax.set_ylabel("$\gamma$",fontsize=35)

ax.tick_params(labelsize=25)


fig, ax = plt.subplots(figsize=[12,8])

ax.plot(om/om0,1/GG)

ax.set_xlabel(r"$\omega/\omega_0$ ",fontsize=35)
ax.set_ylabel("$\zeta$",fontsize=35)

ax.tick_params(labelsize=25)


plt.show()
