import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipj,ellipk

import real_trace as sg

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## SINE-GORDON DIRECT SCATTERING                                             ##
##                                                                           ##
##                        SINE-GORDON REAL LINE TRACE                        ##
##                                                                           ##
## ============================================================================




############################          SETUP        ############################


N = 1000
L = 10
T = 1
T0 = 10

dx = L/N
x = np.linspace(-L/2,L/2,N)

dt = T/N
t = np.linspace(T0,T+T0,N)

xx,tt = np.meshgrid(x,t)

E = np.linspace(-0.4,-0.01,N,dtype="complex")



############################          FUNCS        ############################


def kink(K):

    v = (4*K-1)/(1+4*K)

    phi = np.exp((xx-v*tt)/np.sqrt(1-v**2))

    return 4*np.arctan(phi)


def periodic_kink(mm,v):

    lor = np.sqrt(1-v**2)

    Lambda = np.sqrt(mm) * ellipk(mm) * np.sqrt(1-v**2)
    
    Xmin = 0
    Xmax = 2*Lambda
    Xhalf = (Xmax-Xmin)/2.

    X = np.linspace(Xmin, Xmax, N)
    T = np.linspace(0,T0,N)

    XX,TT = np.meshgrid(X,T)
    
    YYper = 2*ellipj((XX-Xhalf-v*TT)/np.sqrt(mm)/lor,mm)[3] + np.pi

    #YYkink = 4*np.arctan(np.exp((XX-Xhalf-v*tt)/lor))


    return YYper,Xmax

############################          TRACE        ############################

n = 0
u,xmax = periodic_kink(0.99999,0.5)
dx = xmax/N

ux = np.diff(u[n])/dx
ut = (u[n+1]-u[n])/dt

trace = sg.scatter(u[n],ux,ut,E,dx,dt)


############################          PLOT         ############################


plt.figure(1)
plt.title("Trace")
plt.plot(E.real,trace)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")

plt.figure(2)
plt.plot(u[0,:] , label ='am function')


plt.show()
