import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from scipy.optimize import fsolve

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ===========================================================================
## Calculation of three-wave resonances in the torus                        ##
## The possible combination for the first three branches are considered     ##
## ===========================================================================



####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))

def varicose_N(k,w,ro,N):

    return N*np.sqrt((g*(k/N)/ro+sigma*(k/N)**3/rho/ro**3)*np.tanh((k/N)*w*ro/2/rc**2))

def sinuous(k,w,ri):

    om0 = np.pi/2*np.sqrt(9.81*np.sin(9*np.pi/180)/2/w)

    return np.sqrt(om0**2+np.abs(sigma/rho*(k/ri)**3))
    
def sloshing(k,om0):
    
    #return np.sqrt((om0**2+(g*k**2/(0.114)+0.9*sigma*k**3/rho/ro**3)))

     return np.sqrt(om0**2+g*k**2/0.076)

def group_slo(k,om0):

    return g*k/sloshing(k,om0)/0.076

def group_var(k,w,ro):

    a = (g/ro+3*sigma*k**2/rho/ro**3)*np.tanh(k*w*ro/2/rc**2)

    b = (g*k/ro+sigma*k**3/rho/ro**3)/np.cosh(k*w*ro/2/rc**2)**2*w/2*ro/rc**2
    
    return 0.5/varicose(k,w,ro)*(a+b)

def phase_var(k,w,ro):

    return varicose(k,w,ro)/k

####### PHYSICS #######

g = np.sin(4.5*np.pi/180)*9.81

sigma = 0.074     # SURFACE TENSION
rho = 1000        # FLUID DENSITY

rc = 0.095         # CENTRAL RADIUS IN M
ro = 0.08       # TORUS OUTER RADIUS
w = 2*(ro-rc)     # TORUS WIDTH
ri =ro-w          # TORUS INNER RADIUS

q = 4             # DISPERSION RELATION WIDTH

#om0 = 42.4  # R=7.9  # CUTOFF FREQUENCY OF SLOSHING BRANCH IN RAD/S
om0 = 44.5   # R=8.1


####### GRID #######

N = 2000

kmax = 150

k1 = np.linspace(-kmax,kmax,N)
k2 = np.linspace(-kmax,kmax,N)

kk1,kk2 = np.meshgrid(k1,k2)


####### RESONANCE VAR+VAR=VAR #######

fig, ax = plt.subplots(figsize=[8,6])

ax.set_title(r"$\omega_V^1+\omega_V^2=\omega_V^3$")

ax.set_xlabel('$k_1$',fontsize=35)
ax.set_ylabel('$k_2$',fontsize=35)
ax.tick_params(labelsize=25)

k3 = 30

# ax.contour(kk1,kk2,-varicose(kk1,w,ro)+varicose(kk2,w,ro)-varicose(-kk1-kk2,w,ro),[0],colors="b")

# ax.contour(kk1,kk2,varicose(kk1,w,ro)-varicose(kk2,w,ro)+varicose(-kk1+kk2,w,ro),[0],colors="b")


#ax.contour(kk1,kk2,varicose(kk1,w,ro)-varicose(kk2,w,ro)+varicose(k3,w,ro)-varicose(kk1+kk2-k3,w,ro),[0],colors="r")



ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-varicose(kk1+kk2,w,ro),[0],colors="b")
ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-varicose(kk1+kk2,w,ro)+q,[0],colors="b")
ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-varicose(kk1+kk2,w,ro)-q,[0],colors="b")

#ax.contour(kk1,kk2,varicose(kk1,w,ro)+sloshing(kk2,om0)-varicose(kk1+kk2,w,ro),[0],colors="b")
#ax.contour(kk1,kk2,varicose(kk1,w,ro)-varicose(kk2,w,ro)+varicose(kk1-kk2,w,ro),[0],colors="b")

#ax.contour(kk1,kk2,varicose(kk1,w,ro)-varicose(kk2,w,ro)-varicose_N(kk1+kk2,w,ro,2),[0],colors="r")

# ax.contour(kk1,kk2,phase_var(kk2,w,ro)-group_slo(kk1+kk2,om0),[0],colors="r")
# ax.contour(kk1,kk2,phase_var(kk1,w,ro)-group_slo(kk1+kk2,om0),[0],colors="r")
# ax.contour(kk1,kk2,kk1+kk2-5.96,[0],colors="g")

# ax.contour(kk1,kk2,varicose(kk1,w,ro)-varicose(kk2,w,ro)-varicose(kk1+kk2,w,ro),[0],colors="k")

# ax.contour(kk1,kk2,-varicose(kk1,w,ro)+varicose(kk2,w,ro)-varicose(kk1+kk2,w,ro),[0],colors="k")


####### RESONANCE VAR+SIN=VAR #######

# fig, ax = plt.subplots(figsize=[8,6])

# ax.set_title(r"$\omega_V^1+\omega_S^2=\omega_V^3$")

# ax.set_xlabel('$k_1$',fontsize=35)
# ax.set_ylabel('$k_2$',fontsize=35)
# ax.tick_params(labelsize=25)

# k3 = 10

# ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ri)+varicose(k3,w,ri)-varicose(kk1+kk2+k3,w,ro),[0],colors="k")

# ax.contour(kk1,kk2,varicose(kk1,w,ro)-sinuous(kk2,w,ri)-varicose(kk1+kk2,w,ro),[0],colors="k")

# ax.contour(kk1,kk2,-varicose(kk1,w,ro)+sinuous(kk2,w,ri)-varicose(kk1+kk2,w,ro),[0],colors="k")


####### RESONANCE VAR+VAR=SLO #######

# fig, ax = plt.subplots(figsize=[8,6])

# ax.set_title(r"$\omega_V^1+\omega_V^2=\omega_\Sigma^3$")

# ax.set_xlabel('$k_1$',fontsize=35)
# ax.set_ylabel('$k_2$',fontsize=35)
# ax.tick_params(labelsize=25)

#ax.contour(kk1,kk2,-varicose(kk1,w,ro)+varicose(kk1-kk2,w,ro)-sloshing(kk1,om0),[0],colors="b")


#ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-sloshing(kk1+kk2,om0),[0],colors="r")

#ax.contour(kk1,kk2,varicose(kk1,w,ro)+sloshing(kk2,om0)-varicose(kk1+kk2,w,ro),[0],colors="r")

# ax.contour(kk1,kk2,varicose(kk1,w,ro)-varicose(kk2,w,ro)-sloshing(kk1+kk2,om0),[0],colors="k")

# ax.contour(kk1,kk2,-varicose(kk1,w,ro)+varicose(kk2,w,ro)-sloshing(kk1+kk2,om0),[0],colors="k")


plt.show()
