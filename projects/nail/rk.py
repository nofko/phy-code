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


N = 10000
T = 2

t = np.linspace(0,T,N)

dt = t[1]-t[0]

th0 = 1*np.pi/180
u = np.array([th0,0.0])

s = np.zeros((N,2))

g = 9.81
L = 0.1
M = 0.01
mu = 0.15
I = 4*M*L**2/3
om0 = np.sqrt(M*g*L/I)
OM = 9*2*np.pi
B = 120

ft = np.zeros((N))
nt = np.zeros((N))

slide = 0

def FN(th,om,alpha,ch):
    
    NN = M*g-M*L*(np.sin(th)*alpha+om**2*np.cos(th))

    F = M*L*(np.cos(th)*alpha-om**2**np.sin(th))

    if ch==0:
        
        return F, NN
    else:

        return  mu*NN, NN

def rhs(t,u):

    x1_dot = u[1]


    v1_dot = (om0**2+B*np.sin(OM*t))*np.sin(u[0])

    # else:
    #     print(1)
    #     v1_dot = M*L*(g-L*u[1]**2*np.cos(u[0]))*(np.sin(u[0])-mu1*np.cos(u[0]))/(I+M*L**2*np.sin(u[0])*(np.sin(u[0])-mu1*np.cos(u[0])));
        
    
    return np.array([x1_dot,v1_dot])


def rhs_slide(t,u):

    x1_dot = u[1]

    v1_dot = M*L*(g-L*u[1]**2*np.cos(u[0])+B*np.sin(OM*t))*(np.sin(u[0])-mu*np.cos(u[0]))/(I+M*L**2*np.sin(u[0])*(np.sin(u[0])-mu*np.cos(u[0])));
    
    return np.array([x1_dot,v1_dot])


for i in range(N):

    s[i] = u
    

    if slide == 0:

        
        
        a = dt*rhs(t[i],u)
        b = dt*rhs(t[i]+dt/2,u+a/2)
        c = dt*rhs(t[i]+dt/2,u+b/2)
        d = dt*rhs(t[i]+dt,u+c)

        u = u+(a+2*(b+c)+d)/6

        ft[i], nt[i] = FN(u[0],u[1],rhs(t[i],u)[1],slide)
        
        if ft[i]>=nt[i]*mu:
        
            slide += 1
        
    else:

        
        a = dt*rhs_slide(t[i],u)
        b = dt*rhs_slide(t[i]+dt/2,u+a/2)
        c = dt*rhs_slide(t[i]+dt/2,u+b/2)
        d = dt*rhs_slide(t[i]+dt,u+c)

        u = u+(a+2*(b+c)+d)/6

        ft[i], nt[i] = FN(u[0],u[1],rhs_slide(t[i],u)[1],slide)
        
    

th = s[:,0]

plt.figure(1)
plt.plot(t,s[:,0])

plt.axhline(np.pi/2,linestyle="--",color="k")

plt.figure(2)


# F = M*L*om0**2*np.sin(th)*(3*np.cos(th)-2*np.cos(th0))
# N = M*g-M*L*om0**2*(1+2*np.cos(th)*np.cos(th0)-3*np.cos(th)**2)

# plt.plot(t,F*1000)
# plt.plot(t,N*mu*1000)
# plt.plot(t,N*1000)


plt.plot(t,ft*1000)
plt.plot(t,nt*1000)
plt.plot(t,nt*1000*mu)

print(ft)
plt.show()
