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


NT = 1000

T = 0.22

t = np.linspace(0,T,NT)

dt = t[1]-t[0]


th0 = 70*np.pi/180

u = np.array([th0, 0, 0, 0])

s = np.zeros((NT,4))

g = 9.81
L = 0.025
M = 0.00018

# muk = 0.15
# mus = 0.26

muk = 0.8
mus = 1.6

J = 0.0000014
omega = 2*np.pi*2

ft = np.zeros((NT))
nt = np.zeros((NT))

slide = 0

print(J/M/L)
print(np.sqrt(3*J/2/M/L)/2/np.pi)
print(np.sqrt(g/L)/2/np.pi)

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
    

def rhs(t,u,fr):

    th = u[0]
    
    x1_dot = u[1]
    
    if slide == 0:

        x2_dot = 0
        
        v1_dot = -3*g/L*np.cos(u[0])/2 + 3*J**2/M/L*np.sin(2*u[0])*np.cos(omega*t)**2
        v2_dot = 0
    
    else:

        st = np.sin(th)
        ct = np.cos(th)
        mu = - muk * fr
        
        x2_dot = u[3]

        v1_dot = ( -(g-L/2*u[1]**2*st) * (ct-mu*st) + 6*J/M/L*np.sin(2*th)*np.cos(omega*t)**2) / (L/6+L/2*ct*(ct-mu*st))
        
        v2_dot = -u[1]**2*L/2*ct - v1_dot*L/2*st + mu*(L/2*v1_dot-6*J*np.sin(2*th)/L/M*np.cos(omega*t)**2)/(3*ct-3*mu*st)

        print(th,mu,ct,st,L/6+L/2*ct*(ct-mu*st))
        
        # v1_dot = -(2*g/L-u[1]**2*np.sin(u[0]))/(np.cos(u[0])+1/(3*np.cos(u[0])-3*np.sign(-u[3])*muk*np.sin(u[0])))

        # v2_dot = -u[1]**2*L/2*np.cos(u[0])-v1_dot*L/2*(np.sin(u[0])-np.sign(-u[3])*muk/(3*np.cos(u[0])-3*np.sign(-u[3])*muk*np.sin(u[0])))

    

    return np.array([x1_dot,v1_dot,x2_dot,-v2_dot])



############################         SOLVE         ############################



for i in range(NT):

    s[i] = u

    fr = np.sign(u[3])
    
    f, n = FN(u[0],u[1],rhs(t[i],u,fr)[1])

    if abs(f)>=mus*n and slide == 0:

        slide = 1
        print(t[i])

    a = dt*rhs(t[i],u,fr)
    b = dt*rhs(t[i]+dt/2,u+a/2,fr)
    c = dt*rhs(t[i]+dt/2,u+b/2,fr)
    d = dt*rhs(t[i]+dt,u+c,fr)
        
    u = u+(a+2*(b+c)+d)/6

    # if abs(slide)>0 and abs(u[3])<0.001:
            
    #     slide = 0

    #     u[3] = 0

    ft[i], nt[i] = f,n

    

############################         PLOT          ############################


th = s[:,0]

plt.figure(1)
plt.plot(t,s[:,0]*180/np.pi)
#plt.plot(t,s[:,1]*180/np.pi/100)

plt.axhline(0,linestyle="--",color="k")
plt.axhline(180,linestyle="--",color="k")
plt.ylim([-200,200])

#plt.figure(2)

#print(s[:,0]*180/3.14)
#print(s[:,2])

#plt.plot(s[:,0]*180/3.14,s[:,2]*1000)

#plt.plot(s[:,0]*180/3.14,s[:,3]*1000)

#plt.axhline(0,linestyle="--",color="k")

# F = M*L*om0**2*np.sin(th)*(3*np.cos(th)-2*np.cos(th0))
# N = M*g-M*L*om0**2*(1+2*np.cos(th)*np.cos(th0)-3*np.cos(th)**2)

# plt.plot(t,ft*1000)
# #plt.plot(t,muk*nt*1000)
# plt.plot(t,nt*1000)

plt.figure()

plt.plot(t,s[:,1])
plt.plot(t,s[:,3]*100)

plt.axhline(0,linestyle="--",color="k")
plt.ylim([-200,200])

# plt.plot(t,ft*1000)
# plt.plot(t,nt*1000)
#plt.plot(t,nt*1000*mu)
plt.show()
