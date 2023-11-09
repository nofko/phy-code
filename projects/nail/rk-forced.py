import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                          RUNGE-KUTTA FALLING ROD                          ##
##                                                                           ##
## ============================================================================




############################         SETUP         ############################


NT = 10000
T = 16

t = np.linspace(0,T,NT)

dt = t[1]-t[0]
print(dt)

th0 = 89*np.pi/180
u = np.array([th0, 0, 0, 0])

s = np.zeros((NT,4))

g = 9.81
L = 0.2
M = 0.12

muk = 0.5
mus = 1.4

muk = 0.15
mus = 0.26

J = 5
omega = 2*np.pi*40
fi0 = np.pi

ft = np.zeros((NT))
nt = np.zeros((NT))

slide = 0


############################         FUNCS         ############################


def FN(th,om,alpha):
    global slide
    
    # NN = M*g*(1+9*np.sin(th)**2-6*np.sin(th0)*np.sin(th))/4

    # F = M*g*(9*np.sin(th)*np.cos(th)-6*np.sin(th0)*np.cos(th))/4

    NN = M*g+M*L/2*(alpha*np.cos(th)-om**2*np.sin(th))

    F = -M*L/2*(alpha*np.sin(th)+om**2*np.cos(th))
    

    if abs(slide)>0:
        
        return NN*muk, NN

    else:

        return F, NN

def rhs(t,u):

    global slide

    th = u[0]
    x1_dot = u[1]

    x2_dot = u[3]

    
    if slide == 0:
        
        v1_dot = -3*g/L*np.cos(u[0])/2+2*J/2/M/L*np.cos(u[0])*np.sin(omega*t)
        v2_dot = 0
    
    else:

        st = np.sin(th)
        ct = np.cos(th)

        v1_dot = -((g-L/2*u[1]**2*st)*(ct-np.sign(-u[3])*muk*st)+J/M/L* 2*ct*np.sin(omega*t+fi0))/(L/6+L/2*ct*(ct-np.sign(-u[3])*muk*st))

        v2_dot = -u[1]**2*L/2*ct - v1_dot*L/2*st + (np.sign(-u[3])*muk*(L/2*v1_dot-6*J*ct/L/M*np.sin(omega*t+fi0))/(3*ct-3*np.sign(-u[3])*muk*st))
        
    f, n = FN(u[0],u[1],v1_dot)

    
    if abs(f)>=mus*n and slide == 0:
        
        slide = np.sign(f)*1


    # if u[0]>np.pi or u[0]<0:

    #     x1_dot = 0
    #     x2_dot = 0
        
    #     v2_dot = 0
    #     v2_dot = 0

    return np.array([x1_dot,v1_dot,x2_dot,-v2_dot])



############################         SOLVE         ############################



for i in range(NT):

    s[i] = u
        
    a = dt*rhs(t[i],u)
    b = dt*rhs(t[i]+dt/2,u+a/2)
    c = dt*rhs(t[i]+dt/2,u+b/2)
    d = dt*rhs(t[i]+dt,u+c)
        
    u = u+(a+2*(b+c)+d)/6

    if abs(slide)>0 and abs(u[3])<0.01:
            
        slide = 0

    ft[i], nt[i] = FN(u[0],u[1],rhs(t[i],u)[1])

    

############################         PLOT          ############################


th = s[:,0]

plt.figure(1)
plt.plot(t,s[:,0]*180/np.pi)

plt.axhline(0,linestyle="--",color="k")
plt.axhline(180,linestyle="--",color="k")

#plt.figure(2)

#print(s[:,0]*180/3.14)
#print(s[:,2])

#plt.plot(s[:,0]*180/3.14,s[:,2]*1000)

#plt.plot(s[:,0]*180/3.14,s[:,3]*1000)

#plt.axhline(0,linestyle="--",color="k")

# F = M*L*om0**2*np.sin(th)*(3*np.cos(th)-2*np.cos(th0))
# N = M*g-M*L*om0**2*(1+2*np.cos(th)*np.cos(th0)-3*np.cos(th)**2)

# plt.plot(t,F*1000)
# plt.plot(t,N*mu*1000)
# plt.plot(t,N*1000)

plt.figure()
plt.plot(t,s[:,2]*1000)
plt.plot(t,s[:,3]*100)

plt.axhline(0,linestyle="--",color="k")

# plt.plot(t,ft*1000)
# plt.plot(t,nt*1000)
#plt.plot(t,nt*1000*mu)
plt.show()
