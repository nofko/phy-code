import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from scipy.optimize import fsolve

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


####### DISPERSION RELATIONS #######

def varicose(k,w,ro):

    return np.sqrt((g*k/ro+sigma*k**3/rho/ro**3)*np.tanh(k*w*ro/2/rc**2))


def sinuous(k,w,ri):

    om0 = np.pi/2*np.sqrt(9.81*np.sin(9*np.pi/180)/2/w)

    return np.sqrt(om0**2+np.abs(sigma/rho*(k/ri)**3))
    
def sloshing(k,om0):
    
    return np.sqrt(om0**2+g*k**2/rc)

####### PHYSICS #######

g = 9.81*np.sin(4.5*np.pi/180)

sigma = 0.03     # SURFACE TENSION
rho = 1000        # FLUID DENSITY

rc = 0.07         # CENTRAL RADIUS IN M
ro = 0.09       # TORUS OUTER RADIUS
w = 2*(ro-rc)     # TORUS WIDTH
ri =ro-w          # TORUS INNER RADIUS

q = 1             # DISPERSION RELATION WIDTH

#om0 = 42.4  # R=7.9  # CUTOFF FREQUENCY OF SLOSHING BRANCH IN RAD/S
om0 = 44.5   # R=8.1

####### PLOTTING VARICOSE #######


fig, ax = plt.subplots(figsize=[8,6])

k = np.linspace(-40,40,1000)
km = np.linspace(-34,34,1000)

k1 = 10

plt.plot(k,varicose(k,w,ro),c="r")
plt.plot(km+k1,varicose(km,w,ro)+varicose(k1,w,ro),c="b")


####### PLOTTING SLOSHING #######

sigma = 0.055
ro = 0.0785       
w = 2*(ro-rc)

k1 = 17

fig, ax = plt.subplots(figsize=[8,6])
plt.plot(k,varicose(k,w,ro),c="r")
plt.plot(km+k1,varicose(km,w,ro)+varicose(k1,w,ro),c="b")

plt.plot(k,sloshing(k,om0),c="g")


plt.savefig("demo_2.svg")

plt.show()
