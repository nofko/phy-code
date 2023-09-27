import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipj,ellipk
from scipy.optimize import bisect
from scipy.signal import argrelmin,argrelmax


import real_trace as sg


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
L = 40
T = 1

dx = L/N
x = np.linspace(-L/2,L/2,N)

dt = T/N
t = np.linspace(0,T,N)

xx,tt = np.meshgrid(x,t)

Nr = 1000
Ni = 1000

Emin_x = -0.15
Emax_x = -0.1

Emin_y = 0.2
Emax_y = 0.22

dEx = abs(Emax_x-Emin_x)/Nr
dEy = abs(Emax_y-Emin_y)/Ni

Er = np.linspace(Emin_x-dEx,Emax_x+dEx,Nr+3,dtype="complex")

Ei = np.linspace(Emin_y-dEy,Emax_y+dEy,Ni+3,dtype="complex")


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


def indice_contour(tr1,tr2,tr3,tr4):

    """ 1/(2*j*np.pi) * Integrale de contour de F'/F
    FF est le resultat de F appliquee au meshgrid contenu dans XX et YY
    Le contour est un rectangle sur des points de discretisation
    d'indices ix1, ix2, iy1, iy2
    On calcule F' par differences finies centree
    sens de parcours :
    1: (ix1, iy1) --> (ix2, iy1)
    2: (ix2, iy1) --> (ix2, iy2)
    3: (ix2, iy2) --> (ix1, iy2)
    4: (ix1, iy2) --> (ix1, iy1)
    FF[i,j] = F(xi, yj)
    """
        
    # integration horizontale

    DFF = (np.roll(tr1, 1)-np.roll(tr1, -1))/2/dEx/tr1
    VV = DFF[1 : -1]
    I1 = +dEx*(sum(VV)-(VV[0]+VV[-1])/2.)

    DFF = (np.roll(tr3, 1)-np.roll(tr3, -1))/2/dEx/tr3
    VV = DFF[1 : -1]		 
    I3 = -dEx*(sum(VV)-(VV[0]+VV[-1])/2.) 
	
    # integration verticale
    DFF = (np.roll(tr2, 1)-np.roll(tr2, -1))/2/dEy/tr2
    VV = DFF [1 : -1]
    I2 = +dEy*(sum(VV)-(VV[0]+VV[-1])/2.)

    DFF = (np.roll(tr4, 1)-np.roll(tr4, -1))/2/dEy/tr4
    VV = DFF[1 : -1]		 
    I4 = -dEy*(sum(VV)-(VV[0]+VV[-1])/2.) 
    #
    integrale = I1 + I2 + I3 + I4
    resultat = -1./(2.*1j*np.pi) * integrale
    return resultat


############################          TRACE        ############################


n = 100

save = 0

#u = wave(np.pi/4,1)

#u = sine(0.01,2,1)

u = breather(np.pi/3,1,0)
#u += breather(np.pi/3,1,10)

ux = np.diff(u[n])/dx
ut = (u[n+1]-u[n])/dt

print("CALCULATING THE TRACE")


tr1 = sg.scatter(u[n],ux,ut,Er+Emin_y*1j,dx,dt)+2
tr2 = sg.scatter(u[n],ux,ut,Ei+Emax_x,dx,dt)+2
tr3 = sg.scatter(u[n],ux,ut,Er+Emax_y*1j,dx,dt)+2
tr4 = sg.scatter(u[n],ux,ut,Ei+Emin_x,dx,dt)+2

# x0 = -.125
# y0 = +.21

# FF = lambda x,y : ((x-x0) + 1j * (y-y0))**1

# xx, yy = np.meshgrid(Er,Ei)

# plt.figure()
# plt.pcolormesh(Er.real,Ei.real,FF(xx,yy).imag)
# plt.colorbar()

# #plt.contour(Er.real,Ei.real,FF(xx,yy).imag)


# plt.figure()
# plt.pcolormesh(Er.real,Ei.real,FF(xx,yy).real)
# plt.colorbar()

# #plt.contour(Er.real,(Ei/1j).real,FF(xx,yy).real)


# tr1 = FF(Er.real    , Emin_y )
# tr2 = FF(Emax_x, Ei.real  )
# tr3 = FF(Er.real    , Emax_y )
# tr4 = FF(Emin_x, Ei.real  )

print(indice_contour(tr1,tr2,tr3,tr4))

plt.show()
