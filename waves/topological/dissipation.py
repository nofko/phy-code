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

L = 0.2    # M
om0 = np.sqrt(g/L)

d = 0.05

h = 0.05   # M
h1 = 0.01

b = 0.17

nu = 1e-6

### DISSIPATION ###


def dispersion(k,h):

    return np.sqrt(g*k*np.tanh(k*h))

def group_velocity(k,h):

    return 0.5*(np.tanh(k*h)+k*h/(np.cosh(k*h)**2))*g/dispersion(k,h)

def gamma(k,b,h):

    U = group_velocity(k,h)

    omega = dispersion(k,h)

    return 1/U*(2*nu*k**2+np.sqrt(0.5*nu*omega)*(1/b+k/np.sinh(2*k*h)))


k = np.linspace(1,100,1000)

f = dispersion(k,h)/2/np.pi
f1 = dispersion(k,h1)/2/np.pi

#plt.plot(f,1/gamma(k,b,h))
#plt.plot(f1,1/gamma(k,b,h1))


plt.plot(2*np.pi/k,1/gamma(k,b,h))

plt.show()
