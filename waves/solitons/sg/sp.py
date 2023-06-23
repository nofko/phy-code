
N = 1000
L = 3

dx = L/N
dy = L/N

x = np.linspace(-L,L,N)
y = np.linspace(-L,L,N)

xx,yy = np.meshgrid(x,y)

l = 0.5

u = 4*np.arctanh(l*np.cos(np.sqrt(l**2+1)*(xx-1.5))/(np.sqrt(1+l**2)*np.cosh(l*yy)))

ux = np.diff(u[500])/dx
uy = (u[501]-u[500])/dy

Nr = 100
Ni = 100

Er = np.linspace(-0.03,0.03,Nr,dtype="complex")

Ei = np.linspace(-0.03,0.03,Ni,dtype="complex")*1j

def linlog(y):

    N = len(y)
    out = np.zeros(N)

    for i in range(N):
        if np.abs(y[i])<=1:
            out[i] = y[i]
            
        else:
            out[i] = np.sign(y[i])*(1+np.log(np.abs(y[i])))

    return out


M = np.zeros((N,2,2),dtype="complex")

tr = np.zeros((Nr,Ni),dtype="complex")

def B(eta,etax,etat,e):

    Mi = np.identity(2,dtype="complex")
    Mi1 = np.identity(2,dtype="complex")
    m = np.identity(2,dtype="complex")
    

    w = etax*1j+etat[:-1]
    #er =1/(err.real**2+err.imag**2)*(err.real-1j*err.imag)
    
    for i in range(N-1):
    
        

        e11 = np.cos(dx)+1j*np.sin(dx)*(np.sqrt(e)/2*1j-np.cosh(eta[i])/8/np.sqrt(e)*1j)
        e22 = np.cos(dx)-1j*np.sin(dx)*(np.sqrt(e)/2*1j-np.cosh(eta[i])/8/np.sqrt(e)*1j)

        e12 = 1j*np.sin(dx)*((w[i])/4-1j*np.sinh(eta[i])/8/np.sqrt(e))
        e21 = 1j*np.sin(dx)*((w[i])/4+1j*np.sinh(eta[i])/8/np.sqrt(e))
        
        m = np.array([[e11,e12],[e21,e22]])

        #m = linalg.expm(dx*A)


        Mi1 = Mi
        Mi = np.matmul(m,Mi)
    
    return Mi


for j in range(Ni):

    print(j/Ni)

    for k in range(Nr):

        Mi = B(u[500][:-1],ux,uy,Er[k]+Ei[j])
    
        tr[k,j] = np.trace(Mi)

#trM = np.trace(M,axis1=1,axis2=2)

#trace = linlog(trM.real/2)


def linlogmat(y):

    Nx,Ny = y.shape
    out = np.zeros((Nx,Ny))

    for i in range(Nx):
        
        for j in range(Ny):
            
            if np.abs(y[i,j])<=1:
                out[i,j] = y[i,j]
            
            else:
                out[i,j] = np.sign(y[i,j])*(1+np.log(np.abs(y[i,j])))

    return out
