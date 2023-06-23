import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import get_window


N = 20000
T = 2000

fs =  N/T
print(fs)
t = np.linspace(0,T,N)

dt = t[1]-t[0]

u = np.array([0.001,0.0])

s = np.zeros((N,2))

omega0 = 1

A = 0.1
C = 1

omega = 1.8

def rhs(t,u):

    x1_dot = u[1]
    v1_dot = -omega0**2*(1+(-A*np.cos(omega*t)+u[1])/(1-u[0]+A*np.sin(omega*t)))*u[0]#-u[1]**3*1-u[1]*0.1

    return np.array([x1_dot,v1_dot])


for i in range(N):

    s[i] = u

    a = dt*rhs(t[i],u)
    b = dt*rhs(t[i]+dt/2,u+a/2)
    c = dt*rhs(t[i]+dt/2,u+b/2)
    d = dt*rhs(t[i]+dt,u+c)

    u = u+(a+2*(b+c)+d)/6



plt.figure(1)
plt.plot(t,s[:,0])
#plt.plot(t,0.1*np.sin(t))


# plt.figure(2)
#n = 10000


#st = np.sin(t)*np.cos(0.8*t)
#plt.plot(s[n:-1,0],np.diff(s[n:,0]*fs),lw=0.2)

#f = np.linspace(0,N/T,N)*np.pi/2

# w = get_window("flattop",N)

# plt.figure()
#plt.loglog(f,np.abs(np.fft.fft(s[:,0])))
 
plt.show()
