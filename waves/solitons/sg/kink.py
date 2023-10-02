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
##                        SINE-GORDON REAL LINE TRACE                        ##
##                                                                           ##
## ============================================================================




############################          SETUP        ############################


N = 1000


T = 1
T0 = 1

L = 15

dt = T/N
t = np.linspace(T0,T+T0,N)

x = np.linspace(-L/2,L/2,N)
dx = abs(x[1]-x[0])

xx,tt = np.meshgrid(x,t)

E = np.linspace(-2,-0.01,N,dtype="complex")


############################          FUNCS        ############################


def kink(K):

    v = (4*K-1)/(1+4*K)

    phi = np.exp((xx-v*tt)/np.sqrt(1-v**2))

    return 4*np.arctan(phi)


def periodic_kink(mm,v,NN):

    lor = np.sqrt(1-v**2)

    Lambda = NN * np.sqrt(mm) * ellipk(mm) * np.sqrt(1-v**2)
    
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

save = 0

n = 0

m = 0.8
V = 0.3

u = kink(1)

# u, xmax = periodic_kink(m,V,5)

# dx = xmax/N

# x = np.linspace(0,xmax,N)

ux = np.diff(u,axis=1)[n]/dx
ut = np.diff(u,axis=0)[n]/dt

print("CALCULATING THE TRACE")

trace = sg.scatter(u[n],ux,ut,E,dx,dt)


############################          ZEROS        ############################


def trace_finder_plus(E):
    
    Mt = sg.B(u[n],ux,ut,E+0j,dx,dt)

    tr = np.trace(Mt.real)/2

    if np.abs(tr)<=1:
        return tr-1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))-1
    
def trace_finder_minus(E):
    
    Mt = sg.B(u[n],ux,ut,E+0j,dx,dt)

    tr = np.trace(Mt.real)/2
    
    if np.abs(tr)<=1:
        return tr+1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))+1


def find_zeros(E,trace):

    Eigs_plus = []
    Eigs_minus = []
    
    mins = argrelmin(sg.linlog(trace.real))[0]
    maxs = np.concatenate(([int(mins[0]/2)],argrelmax(sg.linlog(trace.real))[0]))


    Emaxs = E[maxs]
    Emins = E[mins]

    Ebr = np.zeros(int(len(Emaxs)+len(Emins)))

    for i in range(len(Emaxs)):

        Ebr[2*i] = Emaxs[i].real

    for i in range(len(Emins)):
    
        Ebr[2*i+1] = Emins[i].real


    for i in range(len(Ebr)-1):

        if (np.sign(trace_finder_plus(Ebr[i]))!=np.sign(trace_finder_plus(Ebr[i+1]))):
    
            Eigs_plus.append(bisect(trace_finder_plus,Ebr[i],Ebr[i+1],xtol=1e-16))

        if (np.sign(trace_finder_minus(Ebr[i]))!=np.sign(trace_finder_minus(Ebr[i+1]))):
        
            Eigs_minus.append(bisect(trace_finder_minus,Ebr[i],Ebr[i+1],xtol=1e-16))


    return Eigs_plus,Eigs_minus



############################          ZEROS        ############################


print("FINDING ZEROS NOW")
Epls, Emns = find_zeros(E.real,trace)

print("+1 eigenvalues: ", Epls)
print("-1 eigenvalues: ", Emns)


############################          SAVE         ############################


print("SAVING DATA")


if save == 1:
    
    name = "m_"+str(m)+"_v_"+str(V)
    wave = "kink-train"

    np.save("/data/sg/"+wave+"/x_"+name+"_N5.npy", x)
    np.save("/data/sg/"+wave+"/eta_"+name+"_N5.npy", u)

    np.save("/data/sg/"+wave+"/E_"+name+"_N5.npy", E)
    np.save("/data/sg/"+wave+"/trace_"+name+"_N5.npy", trace)
    np.save("/data/sg/"+wave+"/Eigs_plus_"+name+"_N5.npy", Epls)
    np.save("/data/sg/"+wave+"/Eigs_minus_"+name+"_N5.npy", Emns)


############################          PLOT         ############################

    
print("PLOT")

plt.figure()
plt.title("Trace")
plt.plot(E.real,sg.linlog(trace.real))
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")
plt.grid()

# plt.scatter(Epls,np.ones_like(Epls))
# plt.scatter(Emns,-np.ones_like(Emns))

plt.figure()
plt.plot(x,u[n])

plt.show()
