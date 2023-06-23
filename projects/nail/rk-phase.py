import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ==========================================================================##
## EXPERIMENTAL DATA ANALYSIS                                                ##
##                                                                           ##
##                           RUNGE-KUTTA INTEGRATION                         ##
##                                                                           ##
## ==========================================================================##



############################          SETUP        ############################


N = 10000                  # NUMBER OF POINTS
T = 200                   # DURATION

fs =  int(N/T)

t = np.linspace(0,T,N)    
dt = t[1]-t[0]


############################          PARAMS       ############################

ic = np.array([-0.1,0.0]) # INITIAL CONDITIONS


############################        FUNCTIONS      ############################


def rhs(t,u,A,omega):

    x1_dot = u[1]

    v1_dot = (A*np.sin(omega*t)*np.sin(u[0])+np.sin(u[0])-u[1]**2*np.cos(u[0])*np.sin(u[0]))/(1/3+np.sin(u[0])**2)

    return np.array([x1_dot,v1_dot])


def solve(u0,omega,A):

    s = np.zeros((N,2))
    u = u0

    for i in range(N):

        s[i] = u

        a = dt*rhs(t[i],u,A,omega)
        b = dt*rhs(t[i]+dt/2,u+a/2,A,omega)
        c = dt*rhs(t[i]+dt/2,u+b/2,A,omega)
        d = dt*rhs(t[i]+dt,u+c,A,omega)

        u = u+(a+2*(b+c)+d)/6

    return s


def bifurcation(omega,A):

    ddt = 1/(omega/2/np.pi)
    ddn = ddt*fs

    th = np.zeros((len(A),2))
    
    for i in range(len(A)):

        sx =  solve(ic,omega,A[i])[:,1]

        pm = sx[-fs:-1]

        if (np.max(np.abs(pm))>np.pi/2):
        
            th[i] = [100,100] 

        else:

            th[i] = [0,0]  
        
        
    return th


def build(omega,AS):

    D = np.zeros((len(omega),len(AS),2))
    
    for i in range(len(omega)):

        print(int(100*i/len(omega)))
        
        ths = bifurcation(omega[i],AS)

        D[i,:,:] = ths

    return D
    

    

############################         PLOTTING      ############################


AS = np.linspace(0.1,5,50)
oms = np.linspace(0.1,10,50)*2*np.pi


plt.figure()

th = build(oms,AS)

np.save("rkp.npy",th)

plt.imshow(np.abs(th[:,::-1,0].transpose()))


plt.show()
