import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from scipy.optimize import fsolve

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ===========================================================================
## Calculation of three-wave resonances for aribtrary dispersion relations  ##
## ===========================================================================



####### DISPERSION RELATIONS #######

def branch_1(k):

    return np.sqrt()


def branch_2(k):

    return



####### PHYSICS #######

g = np.sin(4.5*np.pi/180)*9.81

sigma = 0.055     # SURFACE TENSION
rho = 1000        # FLUID DENSITY

rc = 0.07         # CENTRAL RADIUS IN M
ro = 0.0785       # TORUS OUTER RADIUS
w = 2*(ro-rc)     # TORUS WIDTH
ri =ro-w          # TORUS INNER RADIUS

q = 1             # DISPERSION RELATION WIDTH

#om0 = 42.4  # R=7.9  # CUTOFF FREQUENCY OF SLOSHING BRANCH IN RAD/S
om0 = 44.5   # R=8.1


####### GRID #######

N = 2000

kmax = 150

k1 = np.linspace(-kmax,kmax,N)
k2 = np.linspace(-kmax,kmax,N)

kk1,kk2 = np.meshgrid(k1,k2)

####### PLOT BRANCHES #######

fig, ax = plt.subplots(figsize=[8,6])

ax.set_xlabel('$k$',fontsize=35)
ax.set_ylabel('$\omega$',fontsize=35)

ax.plot(k1,branch_1(k))
ax.plot(k1,branch_2(k))

####### RESONANCE 1+2=2 #######

fig, ax = plt.subplots(figsize=[8,6])

ax.set_xlabel('$k_1$',fontsize=35)
ax.set_ylabel('$k_2$',fontsize=35)
ax.tick_params(labelsize=25)

ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-sloshing(kk1+kk2,om0),[0],colors="b")



plt.show()
