import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipe,ellipk,ellipj
from mpl_toolkits.mplot3d import Axes3D  

# Physics

pxcm = 19/1503


Ro = np.linspace(4.5,7,10)
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

Mtab_fine   = np.linspace(-14,-1, 51,dtype="float64")
Mtab_fine = 10**Mtab_fine
eps = 1e-3
Mtab_coarse = np.linspace(.1+eps, .9-eps, 51,dtype="float64")

Mtab = np.concatenate((Mtab_fine, Mtab_coarse, 1.-Mtab_fine[::-1]))

th = np.linspace(0,2*np.pi,1000,dtype="float64")

profils = np.zeros((len(Mtab),1000),dtype="float64")

tabRES = np.zeros((len(Mtab), 5,len(W)),dtype="float64")

for i in range(len(Mtab)):
    M=Mtab[i]

    EE = ellipe(M)
    KK = ellipk(M)

    PERIOD = 2*np.pi/Nmodes
    H = 1/((PERIOD/4/KK)**2/3/abs(D)/M*B)
  
    LAMBDA = 2 * np.sqrt(3*abs(D)*M/(B*H))
    checkPERIOD = 4*np.sqrt(3*abs(D)*M/(B*H))*KK/2/np.pi
  
    v = + np.sign(D) * 2*B*H/3/M * (1-M/2-3/2*EE/KK)

    
    tabRES[i,0] = D
    tabRES[i,1] = B
    tabRES[i,2] = M
    tabRES[i,3] = H*np.sign(D)
    tabRES[i,4] = v



plt.figure(5)
plt.semilogx(1-tabRES[:,2], tabRES[:,3], '-o')
plt.xlabel('$\mu = 1-m$')
plt.ylabel('$h$')

plt.grid()

###############################################################
plt.figure(6)
plt.plot(tabRES[:,3], tabRES[:,4]+1, '-')
plt.xlabel('A [cm]')
plt.ylabel('$\Omega/\Omega_0$')

plt.grid()

#plt.figure()
#plt.polar(th,Ro+profils[-2])
#plt.polar(th,Ro+profils[-10])
#plt.polar(th,Ro+profils[-5])
#plt.polar(th,Ro+profils[])

n = 100
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

hs = tabRES[:n,3,:]

X, Y = np.meshgrid(W, hs[:n,0])

Z = 1+tabRES[:n,4,:]

#ax.set_ylim([0,1])

ax.plot_surface(X, hs, Z)
