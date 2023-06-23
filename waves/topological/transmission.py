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

# L = 1    # M
# om0 = np.sqrt(g/L)

# d = 0.2*L

# h = L   # M
# h1 = 0.05*L


L = 0.2    # M
om0 = np.sqrt(g/L)

d = 0.05

h = 0.05   # M
h1 = 0.01

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

N = 100 

Nf = 200

om = np.linspace(0.1,2.1,Nf)*2*np.pi

TTtot = np.zeros((Nf,2,2),dtype=complex)

x = np.arange(N)*L

for i in range(Nf):
    
    k = inv_gravity(om[i],h)[0]
    k1 = inv_gravity(om[i],h1)[0]

    TT = np.identity(2)
    
    for j in range(N):

        tt = T(k,k1,d,L,x[j])
        
        TT = TT@tt

    TTtot[i] = TT



fig, ax = plt.subplots(figsize=[12,8])

ax.plot(om/2/np.pi,np.log(np.abs(1/TTtot[:,0,0])**2))

ax.set_xlabel(r"$\omega/\omega_0$ ",fontsize=35)
ax.set_ylabel("$\ln(T)$",fontsize=35)

ax.tick_params(labelsize=25)

plt.show()
