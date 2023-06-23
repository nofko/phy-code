import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipe,ellipk
from scipy.optimize import fsolve


Ro = np.linspace(4.5,7,100)
Rc = 4

W = 2*(Ro-Rc)

L = 2*Ro*np.pi

mu = 3.3

Bo = 0.85**2*Rc**4/W**2/Ro**4

deltaBo = (1/23.6 - Bo)

#Velocity

D = W**2*Ro**2*deltaBo/2/Rc**4
B = mu/W


Nmodes = 2

PERIOD = 2*np.pi/Nmodes

v=np.zeros(len(Ro))

A = np.linspace(0.4,0.4,100)

#A = 0.08*(W-3.04)**2+0.05

th = np.linspace(0,2*np.pi,1000)

amps=np.zeros(len(Ro))

Ms = np.zeros(len(Ro))
Hs = np.zeros(len(Ro))

def H_zeros(M,D,B,A):

    EE = ellipe(M)
    KK = ellipk(M)

    PERIOD = 2*np.pi/Nmodes
    H = 1/((PERIOD/4/KK)**2/3/abs(D)/M*B)

    return H-A

def H(M,D,B):

    EE = ellipe(M)
    KK = ellipk(M)

    PERIOD = 2*np.pi/Nmodes
    H = 1/((PERIOD/4/KK)**2/3/abs(D)/M*B)

    return H


for i in range(len(Ro)):

    M = fsolve(H_zeros,1-1e-4,args=(D[i],B[i],A[i]),xtol=1e-12)[0]
    
    print(M)

    Ms[i] = M
    Hs[i] = H(M,D[i],B[i])
    
    EE = ellipe(M)
    KK = ellipk(M)
    
    v[i] = 1 + np.sign(D[i]) * 2*B[i]*H(M,D[i],B[i])/3/M * (1-M/2-3/2*EE/KK)

    amps[i]=A[i]*np.sign(D[i])

    
plt.scatter(amps/W,v)


def profil(m,D,B):

    h = H(m,D,B)
    
    eta = np.sign(D)*h*ellipj(np.sqrt(B*h/12/abs(D)/m)*th,m)[1]**2

    return eta

