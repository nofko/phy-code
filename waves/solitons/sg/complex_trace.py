import numpy as np
import sys

## ============================================================================
## SINE-GORDON DIRECT SCATTERING                                             ##
##                                                                           ##
##                          SINE-GORDON COMPLEX TRACE                        ##
##                                                                           ##
## ============================================================================



def linlogscalar(y):
    
    """Scaling of a scalar value - linear between +-1 and log outside"""
    
    if np.abs(y)<=1:
        return(y)
            
    else:
        return np.sign(y)*(1+np.log(np.abs(y)))


def B(eta,etax,etat,e,dx,dt):
    
    """Calculation of monodromy matrix for given values of the signal, derivative in x and t and for a given energy"""
    
    Mi = np.identity(2,dtype="complex")
    Mi1 = np.identity(2,dtype="complex")
    m = np.identity(2,dtype="complex")

    N = len(eta)

    w = etax-etat[:-1]
    #er =1/(err.real**2+err.imag**2)*(err.real-1j*err.imag)
    
    for i in range(N-1):
    
        

        e11 = np.cos(dx)+1j*np.sin(dx)*(np.sqrt(e)/2-np.cos(eta[i])/8/np.sqrt(e))
        e22 = np.cos(dx)-1j*np.sin(dx)*(np.sqrt(e)/2-np.cos(eta[i])/8/np.sqrt(e))

        e12 = 1j*np.sin(dx)*((w[i])/4+1j*np.sin(eta[i])/8/np.sqrt(e))
        e21 = 1j*np.sin(dx)*((w[i])/4-1j*np.sin(eta[i])/8/np.sqrt(e))
        
        m = np.array([[e11,e12],[e21,e22]])

        #m = linalg.expm(dx*A)

        Mi1 = Mi
        Mi = np.matmul(m,Mi)
    
    return Mi


def scatter(u,ux,ut,Er,Ei,dx,dt):
    
    """Scattering of a given signal over a mesh of energies - returns the trace of the monodromy matrix
    for a given range of energies on the mesh

    """

    N = len(u)
    Nr = len(Er)
    Ni = len(Ei)
    
    M = np.zeros((N,2,2),dtype="complex")

    tr = np.zeros((Nr,Ni),dtype="complex")

    for j in range(Ni):

        sys.stdout.flush()
        sys.stdout.write("Progress: "+str(int(100*j/(Ni)))+" %\r")

        for k in range(Nr):

            Mi = B(u[:-1],ux,ut,Er[k]+Ei[j],dx,dt)
    
            tr[k,j] = np.trace(Mi)

    return tr


def linlogmat(y):
    """Scaling of a matrix - linear between +-1 and log outside"""
    
    Nx,Ny = y.shape
    out = np.zeros((Nx,Ny))

    for i in range(Nx):
        
        for j in range(Ny):
            
            if np.abs(y[i,j])<=1:
                out[i,j] = y[i,j]
            
            else:
                out[i,j] = np.sign(y[i,j])*(1+np.log(np.abs(y[i,j])))

    return out
