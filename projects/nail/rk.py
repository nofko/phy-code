import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                          RUNGE-KUTTA FALLING ROD                          ##
##                                                                           ##
## ============================================================================


N = 100000
T = 500

t = np.linspace(0,T,N)

dt = t[1]-t[0]

u = np.array([0.1,0.0])

s = np.zeros((N,2))

omega = 9
A = 9

def rhs(t,u):

    x1_dot = u[1]

    #v1_dot = (A*np.cos(omega*t)-1)*np.sin(u[0])
    v1_dot = (A*np.cos(omega*t)*np.sin(u[0])+np.sin(u[0])-u[1]**2*np.cos(u[0])*np.sin(u[0]))/(1/3+np.sin(u[0])**2)
    
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

plt.axhline(np.pi,linestyle="--",color="k")

plt.show()
