import numpy as np
import sys

## ============================================================================
## SINE-GORDON DIRECT SCATTERING                                             ##
##                                                                           ##
##                        SINE-GORDON REAL LINE TRACE                        ##
##                                                                           ##
## ============================================================================




def linlog(y):

    N = len(y)
    out = np.zeros(N)

    for i in range(N):
        if np.abs(y[i])<=1:
            out[i] = y[i]
            
        else:
            out[i] = np.sign(y[i])*(1+np.log(np.abs(y[i])))

    return out


def B(eta,etax,etat,e,dx,dt):


    N = len(eta)
    
    Mi = np.identity(2,dtype="complex")
    Mi1 = np.identity(2,dtype="complex")
    m = np.identity(2,dtype="complex")


    w = etax-etat[:-1]
    
    for i in range(N-1):
    
        
        e11 = np.cos(dx)+1j*np.sin(dx)*(np.sqrt(e)/2-np.cos(eta[i])/8/np.sqrt(e))
        e22 = np.cos(dx)-1j*np.sin(dx)*(np.sqrt(e)/2-np.cos(eta[i])/8/np.sqrt(e))

        e12 = 1j*np.sin(dx)*((w[i])/4+1j*np.sin(eta[i])/8/np.sqrt(e))
        e21 = 1j*np.sin(dx)*((w[i])/4-1j*np.sin(eta[i])/8/np.sqrt(e))
        
        m = np.array([[e11,e12],[e21,e22]])

        Mi1 = Mi
        Mi = np.matmul(m,Mi)

    
    return Mi


def scatter(u,ux,ut,E,dx,dt):

    N = len(u)
    Ne = len(E)
    
    M = np.zeros((Ne,2,2),dtype="complex")

    for j in range(Ne):

        sys.stdout.flush()
        sys.stdout.write("Progress: "+str(int(100*j/(Ne-1)))+" %\r")

        Mi = B(u,ux,ut,E[j],dx,dt)
    
        M[j] = Mi

    trM = np.trace(M,axis1=1,axis2=2)
            
    return trM
