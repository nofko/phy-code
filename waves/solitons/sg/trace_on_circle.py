import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipj,ellipk
from scipy.optimize import bisect
from scipy.signal import argrelmin,argrelmax
import sys

import real_trace as sg


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## SINE-GORDON DIRECT SCATTERING                                             ##
##                                                                           ##
##                        SINE-GORDON CIRCLE TRACE                           ##
##                                                                           ##
## ============================================================================




############################          SETUP        ############################


N = 1000
L = 25
T = 1

dx = L/N
x = np.linspace(-L/2,L/2,N)

dt = T/N
t = np.linspace(0,T,N)

xx,tt = np.meshgrid(x,t)

theta = np.linspace(1.5,2.5,1000,dtype=np.complex128)

E = np.exp(1j*theta)*0.25


############################          FUNCS        ############################


def breather(mu,t0,x0):
    """Returns the spatio-temporal form of a breather with a given phase mu of the eigenvalue"""

    
    #mu = np.pi/3
    l0 = 0.5*np.exp(1j*mu)
    u = 4*np.arctan(np.tan(mu)*np.sin((tt-t0)*np.cos(mu))/np.cosh((xx-x0)*np.sin(mu)))

    return u



############################          TRACE        ############################

n = 100

save = 0

#u = wave(np.pi/4,1)

#u = sine(0.01,2,1)

u = breather(np.pi/3+0.05,1,0)
#u += breather(np.pi/3,1,10)

ux = np.diff(u[n])/dx
ut = (u[n+1]-u[n])/dt

print("CALCULATING THE TRACE")


tr_r,tr_i,tr_a = sg.scatter(u[n],ux,ut,E,dx,dt)

############################          PLOT         ############################

    
print("PLOT")



plt.figure()
plt.title("Input")
plt.plot(x,u[n])


plt.figure()
plt.title("Trace - imaginary part")
plt.plot(theta,tr_r)
plt.plot(theta,tr_i)

plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
#plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")
plt.grid()

# plt.figure()
# plt.title("Trace - real part")
# plt.plot(theta,tr_i)
# plt.axhline(1,0,1,ls="--",color="tab:red")
# plt.axhline(-1,0,1,ls="--",color="tab:red")
# #plt.axvline(0,0,1,ls="--",color="tab:red")
# plt.xlabel("E")
# plt.ylabel("Trace")
# plt.grid()


plt.figure()
plt.title("Re-1/2+Im")
plt.plot(theta,tr_a)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
#plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")
plt.grid()


plt.show()
