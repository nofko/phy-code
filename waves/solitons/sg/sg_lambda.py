
N = 500
L = 100
T = 1

dx = L/N
x = np.linspace(-L/2,L/2,N)

dt = T/N
t = np.linspace(0,1,N)

xx,tt = np.meshgrid(x,t)

K = 1
v = 8*np.sqrt(K)/(1+16*np.sqrt(K))
v=(4*K-1)/(1+4*K)
phi1 = np.exp((xx-v*tt+10)/np.sqrt(1-v**2))

K = 1.5
v = 8*np.sqrt(K)/(1+16*np.sqrt(K))
v=(4*K-1)/(1+4*K)
phi2 = np.exp(-(xx-v*tt-10)/np.sqrt(1-v**2))

u = 4*np.arctan(phi1)#+4*np.arctan(phi2)


#uxt = 4*phi[0]/np.sqrt(1-v**2)/(1+phi[0]**2)
#utt = -4*v*phi[0]/np.sqrt(1-v**2)/(1+phi[0]**2)

mu = np.pi/4
l0 = -0.5*np.exp(-1j*mu)

#u = 4*np.arctan(np.tan(mu)*np.sin(tt*np.cos(mu)+0.1)/np.cosh(xx*np.sin(mu)))

ux = np.diff(u[0])/dx
ut = (u[1]-u[0])/dt

lam = np.linspace(-2,2,N)*1j+0.2


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

def B(eta,etax,etat,l):

    Mi = np.identity(2,dtype="complex")
    Mi1 = np.identity(2,dtype="complex")
    m = np.identity(2,dtype="complex")
    
    w = etax-etat[:-1]
    
    for i in range(N-1):
    
        
        #e11 = -1j*w[i]/4
        #e12 = 1j*np.exp(-1j*eta[i])/16/np.sqrt(abs(e)) + 1j*np.sqrt(abs(e))
        #e21 = -1j*np.exp(1j*eta[i])/16/np.sqrt(abs(e)) - 1j*np.sqrt(abs(e))
        #e22 = 1j*w[i]/4

        e11 = 1j*l/2-1j*np.cos(eta[i])/8/l
        e22 = -1j*l/2+1j*np.cos(eta[i])/8/l

        e12 = 1j/4*w[i]-np.sin(eta[i])/8/l
        e21 = 1j/4*w[i]+np.sin(eta[i])/8/l
        
        A = np.array([[e11,e12],[e21,e22]])

        m = linalg.expm(dx*A)


        Mi1 = Mi
        Mi = np.matmul(m,Mi)
    
    return Mi


for j in range(N-1):

    print(j/N)

    Mi = B(u[0],ux,ut,lam[j])
    
    M[j] = Mi

trM = np.trace(M,axis1=1,axis2=2)

trace = linlog(trM.real/2)


plt.figure()
plt.title("Trace")
plt.plot(lam.imag,trace)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")
