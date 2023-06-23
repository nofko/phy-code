from scipy.special import ellipj

N = 1000
L = 60
T = 1

dx = L/N
x = np.linspace(-L/2,L/2,N)

dt = T/N
t = np.linspace(0,T,N)

xx,tt = np.meshgrid(x,t)

Nr = 400
Ni = 400

Er = np.linspace(0.01,1.2,Nr,dtype="complex")

Ei = np.linspace(0.01,1.2,Ni,dtype="complex")*1j


def breather(mu,t0):
    """Returns the spatio-temporal form of a breather with a given phase mu of the eigenvalue"""

    
    #mu = np.pi/3
    l0 = 0.5*np.exp(1j*mu)
    u = 4*np.arctan(np.tan(mu)*np.sin((tt-t0)*np.cos(mu))/np.cosh(xx*np.sin(mu)))

    return u
    
def wave(alpha,A):
    """Returns the spatio-temporal form of a travelling wave with a given phase alpha of the eigenvalue"""
    
    #alpha = np.pi/4
    e1 = A*np.exp(1j*alpha)

    K = np.sqrt(e1*np.conjugate(e1))

    v = ((-4*K-1)/(-4*K+1)).real
    print(v)
    gamma = np.sqrt((1-np.cos(alpha))/2)

    je = ellipj((xx-v*tt)/np.sqrt(v**2-1),gamma**2)
    u = 2*np.arcsin(gamma*je[0])

    return u


def linlogscalar(y):
    """Scaling of a scalar value - linear between +-1 and log outside"""
    
    if np.abs(y)<=1:
        return(y)
            
    else:
        return np.sign(y)*(1+np.log(np.abs(y)))


def B(eta,etax,etat,e):
    """Calculation of monodromy matrix for given values of the signal, derivative in x and t and for a given energy"""
    
    Mi = np.identity(2,dtype="complex")
    Mi1 = np.identity(2,dtype="complex")
    m = np.identity(2,dtype="complex")
    

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


def scatter(u,ux,ut):
    """Scattering of a given signal over a mesh of energies - returns the trace of the monodromy matrix for a given 
    range of energies on the mesh"""

    M = np.zeros((N,2,2),dtype="complex")

    tr = np.zeros((Nr,Ni),dtype="complex")

    for j in range(Ni):

        print(j/Ni)

        for k in range(Nr):

            Mi = B(u[:-1],ux,ut,Er[k]+Ei[j])
    
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

n = 100

u = wave(np.pi/4,1)

ux = np.diff(u[n])/dx
ut = (u[n+1]-u[n])/dt

trace = scatter(u[n],ux,ut)

plt.figure()
plt.pcolormesh(Er.real,Ei.imag,linlogmat(np.transpose(trace.real)))
plt.title(r"$\Re(trM)$")

plt.figure()
plt.pcolormesh(Er.real,Ei.imag,linlogmat(np.transpose(trace.imag)))
plt.title(r"$\Im(trM)$")
