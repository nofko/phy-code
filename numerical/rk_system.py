import matplotlib.pyplot as plt
import numpy as np

N = 1000
T = 10

t = np.linspace(0,T,N)

dt = t[1]-t[0]

u = np.array([0.5,-0.5,0,0])

s = np.zeros((N,4))

k1 = 1
k2 = 2

def rhs(t,u):

    x1_dot = u[2]
    x2_dot = u[3]
    v1_dot = -(k1+k2)*u[0]+k2*u[1]
    v2_dot = -(k1+k2)*u[1]+k2*u[0]
    
    return np.array([x1_dot,x2_dot,v1_dot,v2_dot])


for i in range(N):

    s[i] = u

    a = dt*rhs(t[i],u)
    b = dt*rhs(t[i]+dt/2,u+a/2)
    c = dt*rhs(t[i]+dt/2,u+b/2)
    d = dt*rhs(t[i]+dt,u+c)

    u = u+(a+2*(b+c)+d)/6



plt.figure(1)
plt.plot(t,s[:,0])

omega1 = np.sqrt((k1+2*k2))
omega2 = np.sqrt((k2))

plt.plot(t, 0.5 * np.cos(omega1 * t))


plt.show()
