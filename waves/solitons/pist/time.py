pxcm = 19/1525
Ro =8.17

W =2.49


Rc = 7

L = 2*Ro*np.pi

## Setup

N = 2705

x = np.linspace(0,L,N)
dx = L/N

## Physics

mu = 3

Bo = 0.85**2*Rc**4/W**2/Ro**4
deltaBo = 1/23.6 - Bo

lam = mu*Rc**4/3/W**3/Ro**4/abs(deltaBo)

u = -np.load("data/signal_21_10_4.npy")*pxcm


## Potential and energies



## Calculating M



def linlog(y):

    N = len(y)
    out = np.zeros(N)

    for i in range(N):
        if np.abs(y[i])<=1:
            out[i] = y[i]

        else:
            out[i] = np.sign(y[i])*(1+np.log(np.abs(y[i])))

    return out


## The M matrix for a given energy E

def B(q):

    Mi = np.identity(2)
    #Mi1 = np.identity(2)
    
    s = 0

    for m in range(N):

        if q[m]<0:

            cos = np.cosh(dx*np.sqrt(-q[m]))
            sh = np.sinh(dx*np.sqrt(-q[m]))
            
            exp12 = sh/np.sqrt(-q[m])
            exp21 = sh*np.sqrt(-q[m])
            
        else:

            cos = np.cos(dx*np.sqrt(q[m]))
            sin = np.sin(dx*np.sqrt(q[m]))

            exp12 = sin/np.sqrt(q[m])
            exp21 = -sin*np.sqrt(q[m])
            
        exp = np.array([[cos,exp12],[exp21,cos]])

        #Mi1 = Mi
        Mi = np.dot(exp,Mi)
        
        #s += np.abs(np.sign(Mi[0,0])-np.sign(Mi1[0,0]))/2

    return Mi,s

def M12_finder(E):
    Mt, sst = B(lam*u+E)
    
    return Mt[0,1]

def trace_finder_plus(E,n):
    Mt, sst = B(lam*u[n]+E)

    tr = np.trace(Mt)/2
    
    if np.abs(tr)<=1:
        return tr-1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))-1
    
def trace_finder_minus(E,n):
    Mt, sst = B(lam*u[n]+E)

    tr = np.trace(Mt)/2
    
    if np.abs(tr)<=1:
        return tr+1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))+1




E_all=[]

for i in range(len(u)):
    
    print(i/len(u))

    M = np.zeros((N,2,2))
    S = np.zeros(N)

    E_min = -lam*(np.max(u[i])-np.mean(u[i]))

    E_max =  E_min + 4*np.abs(E_min) # It can also be set to the Nyquist frequency np.pi/dx

    E = np.linspace(E_min,E_max,N)
    
    for j in range(N):
        
        
        q = lam*u[i][:N]+E[j]

        Mi, s = B(q)
    
        S[j] = s
        M[j] = Mi

    trM = np.trace(M,axis1=1,axis2=2)

    trace = linlog(trM/2)
    M11 = linlog(M[:,0,0])
    M12 = linlog(M[:,0,1])

    Eigs_plus = []
    Eigs_minus = []

    mins = argrelmin(trace)[0]
    maxs = np.concatenate(([int(mins[0]/2)],argrelmax(trace)[0]))


    Emaxs = E[maxs]
    Emins = E[mins]

    Ebr = np.zeros(int(len(Emaxs)+len(Emins)))

    for k in range(len(Emaxs)):

        Ebr[2*k] = Emaxs[k]

    for k in range(len(Emins)):
    
        Ebr[2*k+1] = Emins[k]




    for k in range(len(Ebr)-1):

        if (np.sign(trace_finder_plus(Ebr[k],i))!=np.sign(trace_finder_plus(Ebr[k+1],i))):

            Eigs_plus.append(bisect(trace_finder_plus,Ebr[k],Ebr[k+1],args=(i,),xtol=1e-13))
            
        if (np.sign(trace_finder_minus(Ebr[k],i))!=np.sign(trace_finder_minus(Ebr[k+1],i))):
            
            Eigs_minus.append(bisect(trace_finder_minus,Ebr[k],Ebr[k+1],args=(i,),xtol=1e-13))


    Eref =  0

    Eigs_plus = np.array(Eigs_plus)-Eref
    Eigs_minus = np.array(Eigs_minus)-Eref
    E = E-Eref

    Eigs = np.concatenate((Eigs_plus,Eigs_minus))
    Eigs = np.sort(Eigs)

    E_all.append(Eigs[:8])
