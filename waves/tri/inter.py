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

    return g*k/sloshing(k,om0)/0.07

def group_var(k,w,ro):

    a = (g/ro+3*sigma*k**2/rho/ro**3)*np.tanh(k*w*ro/2/rc**2)

    b = (g*k/ro+sigma*k**3/rho/ro**3)/np.cosh(k*w*ro/2/rc**2)**2*w/2*ro/rc**2
    
    return 0.5/varicose(k,w,ro)*(a+b)

def phase_var(k,w,ro):

    return varicose(k,w,ro)/k

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

kmax = 20

k1 = np.linspace(0,kmax,N)
k2 = np.linspace(-kmax,0,N)

kk1,kk2 = np.meshgrid(k1,k2)


####### RESONANCE VAR+VAR=VAR #######

fig, ax = plt.subplots(figsize=[8,6])

ax.set_title(r"$\omega_V^1+\omega_V^2=\omega_V^3$")

ax.set_xlabel('$k_1$',fontsize=35)
ax.set_ylabel('$k_2$',fontsize=35)
ax.tick_params(labelsize=25)

cs = ax.contour(kk1,kk2,varicose(kk1,w,ro)+varicose(kk2,w,ro)-sloshing(kk1+kk2,om0),[0],colors="r")
#ax.contour(kk1,kk2,varicose(kk1,w,ro)+sloshing(kk2,om0)-varicose(kk1+kk2,w,ro),[0],colors="b")


data = cs.allsegs[0][0]

kr1 = data[:,0]
kr2 = data[:,1]


or1 = varicose(kr1,w,ro)
or2 = varicose(kr2,w,ro)


or3 = or1+or2 
kr3 = kr1+kr2



J = -0.25*(or1*or2*(1+np.sign(kr1*kr2))+
           or2*or3*(1+np.sign(kr2*kr3))+
           or3*or1*(1+np.sign(kr3*kr1)))

B1 = J*kr1/or1
B2 = J*kr2/or2
B3 = J*kr3/or3

plt.figure()

plt.plot(or2,np.sqrt(np.abs(B1*B2**2*B3)))
#plt.plot(or2/2/np.pi,B2*B3)
#plt.plot(or1/2/np.pi,B3)

plt.axhline(0)

plt.show()
