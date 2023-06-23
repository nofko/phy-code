import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## THEORETICAL ANALYSIS                                                      ##
##                                                                           ##
##                  LOCALIZATION LENGTH OVER A RANDOM 1D BOTTOM              ##
##                                                                           ##
## ============================================================================


############################        PHYSICS        ############################

g = 9.81   # M/S^2
sigma = 0.0
rho = 1000

L = 0.2    # M
om0 = np.sqrt(g/L)

d = 0.08

h = 0.05   # M
h1 = 0.01


############################          SETUP        ############################

N = 18     # NUMBER OF STEPS

Nf = 100    # FREQUNCY GRID

M = 100    # NUMBER OF ENSEMBLES

R = 0.0    # RANDOMNESS

om = np.linspace(0.1,2.1,Nf)*2*np.pi

ST = np.zeros(Nf)
GG = np.zeros(Nf)

############################        FUNCTIONS      ############################


def disp_gravity(KK,HH):
    """
    Dispersion relation
    """

    return np.sqrt((g*KK+sigma*KK**3/rho)*np.tanh(KK*HH))


def zero_gravity(KK,HH,OM):
    """
    Function to solve for K
    """

    return disp_gravity(KK,HH)-OM


def inv_gravity(OM,HH):
    """
    Finding wavenumber for given frequency
    """
    
    sol = fsolve(zero_gravity,0.2,args=(HH,OM))
        
    return sol


def T(k,ki,d,L,x):
    """
    Transfer matrix for each cell
    """

    g = k/ki

    t0 = np.array([[(1+g)*np.exp(1j*(ki-k)*x),(1-g)*np.exp(-1j*(ki+k)*x)],
                   [(1-g)*np.exp(1j*(ki+k)*x),(1+g)*np.exp(-1j*(ki-k)*x)]])

    t1 = np.array([[(1+1/g)*np.exp(1j*(k-ki)*(x+d)),(1-1/g)*np.exp(-1j*(ki+k)*(x+d))],
                   [(1-1/g)*np.exp(1j*(ki+k)*(x+d)),(1+1/g)*np.exp(-1j*(k-ki)*(x+d))]])

    return 0.25*t0@t1

############################         SOLVING      ############################

def fullT(N,r):
    """
    Multiplication of transfer matrices for all N cells, the placement is chosen at random
    RETURN: 
        T[1,1] - first element of total transfer matrix
        gamma/N - Lyapunov exponent
    """
    
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
    gamma = -np.log(A)

    return A**2,gamma/N


for i in range(M):

    a, gamma = fullT(N,R)

    ST += a
    GG += gamma

ST /= M 
GG /= M

A = np.log(ST)


############################         PLOTTING      ############################

fig, ax = plt.subplots(figsize=[12,8])

ax.plot(om/2/np.pi,A)

ax.set_xlabel(r"$\omega/\omega_0$ ",fontsize=35)
ax.set_ylabel("$\ln(T)$",fontsize=35)

ax.tick_params(labelsize=25)


fig, ax = plt.subplots(figsize=[12,8])

ax.plot(om/2/np.pi,GG)

ax.set_xlabel(r"$\omega/\omega_0$ ",fontsize=35)
ax.set_ylabel("$\gamma$",fontsize=35)

ax.tick_params(labelsize=25)


fig, ax = plt.subplots(figsize=[12,8])

ax.plot(om/2/np.pi,1/GG*L)

ax.set_xlabel(r"$f$ [Hz] ",fontsize=35)
ax.set_ylabel("$\zeta$ [m]",fontsize=35)

ax.set_ylim(0,15)

ax.tick_params(labelsize=25)


plt.show()
