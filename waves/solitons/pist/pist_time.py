mat = scipy.io.loadmat("data/canal/24-1-2022/CumSum_five_solitons_A11cm_08_03_02_05_h3cm_L4m_Deph_Pos150_Mes2m_10s_plage_P.mat")
sig = mat["tableau"]

N0 = len(sig)


u0 = sig[:,1]*1/1.7
#u0 = u0-np.mean(u0)

sos = signal.butter(4,5,'low',fs=2000,output="sos")
u = signal.sosfilt(sos,u0)
#sos2 = signal.butter(4,0.3,'hp',fs=2000,output="sos")
#u = signal.sosfilt(sos2,u)

u = u[4000:]

## Setup

N = len(u)
#N = 10000
T = 10
t = sig[:,0]

#t = np.linspace(0,T,N)
dx = t[1]-t[0]
#dx = T/N

## Physics
h = 3
c0 = np.sqrt(9.81*h*100)
lam = 3/2/h**3*c0**2


##Soliton-theory
v = c0*(1+A/2/h)
delta = np.sqrt(4*h**3/A/3)
#u=A/np.cosh(((t-5)*v)/delta)**2



## Potential and energies

E_min = -lam*(np.max(u)-np.mean(u))

E_max =  E_min + 2*np.abs(E_min) # It can also be set to the Nyquist frequency np.pi/dx

E = np.linspace(E_min,E_max,1000)


## Useful plotting function

def linlog(y):

    N = len(y)
    out = np.zeros(N)

    for i in range(N):
        if np.abs(y[i])<=1:
            out[i] = y[i]

        else:
            out[i] = np.sign(y[i])*(1+np.log(np.abs(y[i])))

    return out


## Calculating M

M = np.zeros((len(E),2,2))
S = np.zeros(N)


## The M matrix for a given energy E

def B(q):

    Mi = np.identity(2)    
    
    for i in range(len(q)):
    
        if q[i]<0:
            sq =np.sqrt(-q[i])
            
            cos = np.cosh(dx*sq)
            sh = np.sinh(dx*sq)
            
            exp12 = sh/sq
            exp21 = sh*sq
            
        elif q[i]>0:
            sq =np.sqrt(q[i])
            
            cos = np.cos(dx*sq)
            sin = np.sin(dx*sq)

            exp12 = sin/sq
            exp21 = -sin*sq

        else:

            cos = 1
            exp12 = 1
            exp21 = 1
            
        exp = np.array([[cos,exp12],[exp21,cos]])

        Mi = np.dot(exp,Mi)
        
    return Mi


## Calculating the matrix M for the range of values between Emin and Emax

for j in range(len(E)):
    print(j/len(E))
    q = lam*u+E[j]
    
    M[j] = B(q)

trM = np.trace(M,axis1=1,axis2=2)

trace = linlog(trM/2)
M11 = linlog(M[:,0,0])
M12 = linlog(M[:,0,1])



plt.figure(1)
plt.title("Trace")
plt.plot(E,trace)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")
