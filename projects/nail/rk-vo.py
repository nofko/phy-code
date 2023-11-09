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




############################         SETUP         ############################



N = 10000
T = 0.5

t = np.linspace(0,T,N)

dt = t[1]-t[0]

th0 = 74*np.pi/180
u = np.array([th0, 0, 0, 0])

s = np.zeros((N,4))

g = 9.81
L = 0.22
M = 0.055


muk = 0.4
mus = 1.9


ft = np.zeros((N))
nt = np.zeros((N))

slide = 0


############################         FUNCS         ############################


def FN(th,om,alpha):
    global slide
    
    NN = M*g*(1+9*np.sin(th)**2-6*np.sin(th0)*np.sin(th))/4

    F = M*g*(9*np.sin(th)*np.cos(th)-6*np.sin(th0)*np.cos(th))/4

    NN = M*g+M*L/2*(alpha*np.cos(th)-om**2*np.sin(th))

    F = -M*L/2*(alpha*np.sin(th)+om**2*np.cos(th))
    
    if abs(slide)>0:
        
        return slide*NN*muk, NN

    else:

        return F, NN

def rhs(t,u):

    global slide
    
    x1_dot = u[1]

    x2_dot = u[3]
    
        
    if slide == 0:
        
        v1_dot = -3*g/L*np.cos(u[0])/2
        v2_dot = 0
    
    else:

        v1_dot = -(2*g/L-u[1]**2*np.sin(u[0]))/(np.cos(u[0])+1/(3*np.cos(u[0])-3*np.sign(-u[3])*muk*np.sin(u[0])))

        v2_dot = -u[1]**2*L/2*np.cos(u[0])-v1_dot*L/2*(np.sin(u[0])-np.sign(-u[3])*muk/(3*np.cos(u[0])-3*np.sign(-u[3])*muk*np.sin(u[0])))

        # if abs(u[3])<1e-3:

        #     slide = 0
        

    f, n = FN(u[0],u[1],v1_dot)
    
    if abs(f)>=mus*n:
        
        slide = np.sign(f)*1
        
    return np.array([x1_dot,v1_dot,x2_dot,-v2_dot])



############################         SOLVE         ############################



for i in range(N):

    s[i] = u
        
    a = dt*rhs(t[i],u)
    b = dt*rhs(t[i]+dt/2,u+a/2)
    c = dt*rhs(t[i]+dt/2,u+b/2)
    d = dt*rhs(t[i]+dt,u+c)
        
    u = u+(a+2*(b+c)+d)/6

    ft[i], nt[i] = FN(u[0],u[1],rhs(t[i],u)[1])

    

############################         PLOT          ############################


th = s[:,0]

plt.figure(1)
plt.plot(t,s[:,0]*180/np.pi)

plt.axhline(0,linestyle="--",color="k")

plt.figure(2)

print(s[:,0]*180/3.14)
print(s[:,2])

plt.plot(s[:,0]*180/3.14,s[:,2]*1000)

#plt.plot(s[:,0]*180/3.14,s[:,3]*1000)

plt.axhline(0,linestyle="--",color="k")

# F = M*L*om0**2*np.sin(th)*(3*np.cos(th)-2*np.cos(th0))
# N = M*g-M*L*om0**2*(1+2*np.cos(th)*np.cos(th0)-3*np.cos(th)**2)

# plt.plot(t,F*1000)
# plt.plot(t,N*mu*1000)
# plt.plot(t,N*1000)

plt.figure()
plt.plot(t,s[:,2]*1000)
plt.plot(t,s[:,3]*1000)

plt.figure()
plt.plot(t,ft)
plt.plot(t,nt*mus)

plt.axhline(0,linestyle="--",color="k")

# plt.plot(t,ft*1000)
# plt.plot(t,nt*1000)
#plt.plot(t,nt*1000*mu)
plt.show()
