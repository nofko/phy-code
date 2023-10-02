import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipj,ellipk
from scipy.optimize import bisect
from scipy.signal import argrelmin,argrelmax


import complex_trace as sg


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)



## ============================================================================
## SINE-GORDON DIRECT SCATTERING                                             ##
##                                                                           ##
##                          SINE-GORDON COMPLEX TRACE                        ##
##                                                                           ##
## ============================================================================




############################          SETUP        ############################


N = 1000
L = 60
T = 1

dx = L/N
x = np.linspace(-L/2,L/2,N)

dt = T/N
t = np.linspace(0,T,N)

xx,tt = np.meshgrid(x,t)

Nr = 400
Ni = 400

Er = np.linspace(-0.3,-0.01,Nr,dtype="complex")

Ei = np.linspace(-0.3,0.3,Ni,dtype="complex")*1j



############################          FUNCS        ############################


def breather(mu,t0,x0):
    """Returns the spatio-temporal form of a breather with a given phase mu of the eigenvalue"""

    
    #mu = np.pi/3
    l0 = 0.5*np.exp(1j*mu)
    u = 4*np.arctan(np.tan(mu)*np.sin((tt-t0)*np.cos(mu))/np.cosh((xx-x0)*np.sin(mu)))

    return u
    
def wave(alpha,A):
    """Returns the spatio-temporal form of a travelling wave with a given phase alpha of the eigenvalue"""
    
    #alpha = np.pi/4
    e1 = A*np.exp(1j*alpha)

    K = np.sqrt(e1*np.conjugate(e1))

    v = ((-4*K-1)/(-4*K+1)).real
    print(v)
    gamma = np.sqrt((1-np.cos(alpha))/2)

    je = ellipj((xx-v*tt)/np.sqrt(v**2-1),gamma**2)
    u = 2*np.arcsin(gamma*je[0])

    return u

def sine(a,k,omega):

    return a*np.sin(k*xx-omega*tt)
    

############################          TRACE        ############################


n = 100

save = 1

#u = wave(np.pi/4,1)

#u = sine(0.01,2,1)

#u = breather(np.pi/3,1)

u = breather(np.pi/3,1,-10)
#u += breather(np.pi/3+0.01,1,10)

ux = np.diff(u[n])/dx
ut = (u[n+1]-u[n])/dt

print("CALCULATING THE TRACE")


plt.figure()
plt.plot(x,u[n])
plt.show()

trace = sg.scatter(u[n],ux,ut,Er,Ei,dx,dt)



############################          SAVE         ############################


print("SAVING DATA")


if save == 1:
    
    name = "breather_2_br"
    wave = "breather"

    np.save("/data/sg/"+wave+"/x_"+name+"_neg.npy", x)
    np.save("/data/sg/"+wave+"/eta_"+name+"_neg.npy", u)

    np.save("/data/sg/"+wave+"/Er_"+name+"_neg.npy", Er)
    np.save("/data/sg/"+wave+"/Ei_"+name+"_neg.npy", Ei)
    np.save("/data/sg/"+wave+"/trace_"+name+"_neg.npy", trace)


############################          PLOT         ############################

plt.figure()
plt.plot(x,u[n])



plt.figure()
plt.pcolormesh(Er.real,Ei.imag,sg.linlogmat(np.transpose(trace.real-2)))
plt.title(r"$\Re(trM)$")
plt.colorbar()
plt.contour(Er.real,Ei.imag,sg.linlogmat(np.transpose(trace.real-2)))
plt.savefig("re.pdf", bbox_inches="tight")

plt.figure()
plt.pcolormesh(Er.real,Ei.imag,sg.linlogmat(np.transpose(trace.imag)))
plt.title(r"$\Im(trM)$")
plt.colorbar()
plt.contour(Er.real,Ei.imag,sg.linlogmat(np.transpose(trace.imag)))
plt.savefig("im.pdf", bbox_inches="tight")



G = (np.abs(trace)-2)**2+trace.imag**2

plt.figure()
plt.pcolormesh(Er.real,Ei.imag,sg.linlogmat(np.transpose(G)))
plt.title(r"$(|\Re(trM)|-2)^2+\Im(trM)^2$")
plt.colorbar()
plt.contour(Er.real,Ei.imag,sg.linlogmat(np.transpose(G)))
plt.savefig("G.pdf", bbox_inches="tight")

plt.show()
