import matplotlib.pyplot as plt
import numpy as np


N = 1000

x = np.linspace(0,50,N)

dx = x[1]-x[0]

#s_exact = -0.75*np.exp(-2*x)+0.5*x+1.75

def rhs(x,u):

    return -u

u = 1

s = np.zeros(N)

for i in range(N):

    s[i] = u

    a = dx*rhs(x[i],u)
    b = dx*rhs(x[i]+dx/2,u+a/2)
    c = dx*rhs(x[i]+dx/2,u+b/2)
    d = dx*rhs(x[i]+dx,u+c)

    u = u+(a+2*(b+c)+d)/6



plt.figure(1)
plt.plot(s)
#plt.plot(s_exact)

#plt.figure(2)
#plt.plot(s-s_exact)

plt.show()
