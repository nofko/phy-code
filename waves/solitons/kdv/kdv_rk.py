import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation



def rhs(tt,VV,uu):
    
    return -1j*kx*0.5*np.fft.fft(uu**2)*np.exp(-1j*kx**3*delta**2*tt)

    

#----- Constructing the grid -----
L   = 2.
nx  = 512
x   = np.linspace(0.,L, nx+1)
x   = x[:nx]  

kx1 = np.linspace(0,int(nx/2-1),int(nx/2))
kx2 = np.linspace(1,int(nx/2),  int(nx/2))
kx2 = -1*kx2[::-1]
kx  = (2.* np.pi/L)*np.concatenate((kx1,kx2))

#----- Parameters -----
delta  = 0.022
delta2 = delta**2

T = 3
nt = 100000
t = np.linspace(0,T,nt)

dt = t[1]-t[0]

u = 0.8*np.exp(-(x-L/2)**2/0.01)

uhat = np.fft.fft(u)
vhat = uhat

sol = np.zeros((nt,nx))

for i in range(nt-1):

    k1 = dt*rhs(t[i],vhat,u)
    k2 = dt*rhs(t[i]+0.5*dt, vhat+0.5*k1,u)
    k3 = dt*rhs(t[i]+0.5*dt, vhat+0.5*k2,u)
    k4 = dt*rhs(t[i]+dt,     vhat+k3,u)
    
    vhat = vhat + (k1+2*k2+2*k3+k4)/6
    
    uhat = vhat*np.exp(1j*kx**3*delta**2*t[i+1])

    u = np.real(np.fft.ifft(uhat))

    sol[i] = u


plt.figure()
# for i in sol:
#     plt.clf()
#     plt.plot(i)
#     plt.pause(0.001)
plt.imshow(sol,aspect=nx/nt)
