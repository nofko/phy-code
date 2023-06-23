import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipj


plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,'font.family' : 'lmodern','svg.fonttype':'none'}
plt.rcParams.update(params)


## ============================================================================
##                                                                           ##
##                                                                           ##
##                        SINE-GORDON REAL LINE TRACE                        ##
##                                                                           ##
## ============================================================================

N = 1000
L = 100
T = 1
T0 = 10

dx = L/N
x = np.linspace(-L/2,L/2,N)

dt = T/N
t = np.linspace(T0,T+T0,N)

xx,tt = np.meshgrid(x,t)

E = np.linspace(-3,-0.01,N,dtype="complex")


def kink(K):

    v = (4*K-1)/(1+4*K)

    phi = np.exp((xx-v*tt)/np.sqrt(1-v**2))

    return 4*np.arctan(phi)


def kink_train(e1,e2,x0):

    e = 0.5*(np.sqrt(e1/e2)+np.sqrt(e2/e1))
    K = np.sqrt(e1*e2)
    v = (4*K-1)/(1+4*K)
    print(e)
    gamma = 2/(e+1)

    je = ellipj((xx-v*tt+x0)/np.sqrt(gamma)/np.sqrt(1-v**2),gamma) 

    return np.pi+2*je[3]


def two_kink(K1,K2):

    v1=(4*K1-1)/(1+4*K1)
    v2=(4*K2-1)/(1+4*K2)

    print(v1)
    print(v2)
    
    phi1 = np.exp((xx-v1*tt)/np.sqrt(1-v1**2))
    phi2 = np.exp((xx-v2*tt+10)/np.sqrt(1-v2**2))

    a1 = np.sqrt((1-v1)/(1+v1))
    a2 = np.sqrt((1-v2)/(1+v2))
    k = (a1-a2)/(a1+a2)

    u = 4*np.arctan(k*(phi1-1/phi2)/(1+phi1/phi2))

    return u
    

def linlog(y):

    N = len(y)
    out = np.zeros(N)

    for i in range(N):
        if np.abs(y[i])<=1:
            out[i] = y[i]
            
        else:
            out[i] = np.sign(y[i])*(1+np.log(np.abs(y[i])))

    return out


def B(eta,etax,etat,e):

    Mi = np.identity(2,dtype="complex")
    Mi1 = np.identity(2,dtype="complex")
    m = np.identity(2,dtype="complex")
    

    w = etax-etat[:-1]
    
    for i in range(N-1):
    
        
        #e11 = -1j*w[i]/4
        #e12 = 1j*np.exp(-1j*eta[i])/16/np.sqrt(abs(e)) + 1j*np.sqrt(abs(e))
        #e21 = -1j*np.exp(1j*eta[i])/16/np.sqrt(abs(e)) - 1j*np.sqrt(abs(e))
        #e22 = 1j*w[i]/4

        #e11 = 1j*np.sqrt(e)/2-1j*np.cos(eta[i])/8/np.sqrt(e)
        #e22 = -1j*np.sqrt(e)/2+1j*np.cos(eta[i])/8/np.sqrt(e)

        #e12 = 1j/4*w[i]-np.sin(eta[i])/8/np.sqrt(e)
        #e21 = 1j/4*w[i]+np.sin(eta[i])/8/np.sqrt(e)

        e11 = np.cos(dx)+1j*np.sin(dx)*(np.sqrt(e)/2-np.cos(eta[i])/8/np.sqrt(e))
        e22 = np.cos(dx)-1j*np.sin(dx)*(np.sqrt(e)/2-np.cos(eta[i])/8/np.sqrt(e))

        e12 = 1j*np.sin(dx)*((w[i])/4+1j*np.sin(eta[i])/8/np.sqrt(e))
        e21 = 1j*np.sin(dx)*((w[i])/4-1j*np.sin(eta[i])/8/np.sqrt(e))
        
        m = np.array([[e11,e12],[e21,e22]])

        #m = linalg.expm(dx*A)


        Mi1 = Mi
        Mi = np.matmul(m,Mi)
    
    return Mi


def scatter(u,ux,ut):

    M = np.zeros((N,2,2),dtype="complex")

    for j in range(N-1):

        print(j/N)

        Mi = B(u,ux,ut,E[j])
    
        M[j] = Mi

    trM = np.trace(M,axis1=1,axis2=2)
        
    trace = linlog(trM.real/2)

    return trace


n = 0

u = kink_train(0.5,0.7,0)
u = kink(1)

ux = np.diff(u[n])/dx
ut = (u[n+1]-u[n])/dt

trace = scatter(u[n],ux,ut)

plt.figure(1)
plt.title("Trace")
plt.plot(E.real,trace)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")

plt.show()
